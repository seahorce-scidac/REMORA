#!/bin/sh
#
# Assert how much an integrated quantity drifts over a run.
#
# usage: check_conservation.sh <datalog> <column> <bound> <below|above> [<column> <bound> <below|above> ...]
#
# Reads a remora.data_log table -- a header row of column names, then one row per
# remora.sum_interval -- and compares the last row against the first as a relative drift,
# |last - first| / |first|.
#
# The "above" mode is the one that keeps these tests honest. A conservation test that only
# ever asserts "drift is small" passes just as well when the case has nothing to conserve:
# no coarse-fine interface, no gradient across it, a correction that has quietly become a
# no-op. Pairing a run that must conserve with a control that must NOT is what makes the
# first assertion mean something.
#
# The sums come from REMORA::sum_integrated_quantities, which weights by Hz/(pm*pn) and masks
# out coarse cells covered by finer ones, so they are totals over the whole hierarchy rather
# than per level. Needs remora.v >= 1; the routine returns immediately otherwise.

set -eu

if [ "$#" -lt 4 ]; then
    echo "usage: $0 <datalog> <column> <bound> <below|above> [...]" >&2
    exit 2
fi

LOG="$1"; shift

if [ ! -f "$LOG" ]; then
    echo "FAIL: data log not found: $LOG" >&2
    exit 1
fi

nrows=$(awk 'NR > 1 && NF > 0' "$LOG" | wc -l)
if [ "$nrows" -lt 2 ]; then
    echo "FAIL: $LOG has $nrows data rows, need at least 2." >&2
    echo "      Check remora.v >= 1, remora.sum_interval > 0 and remora.data_log." >&2
    exit 1
fi

status=0
while [ "$#" -ge 3 ]; do
    COLUMN="$1"; BOUND="$2"; MODE="$3"; shift 3

    RESULT=$(awk -v col="$COLUMN" -v bound="$BOUND" -v mode="$MODE" '
        NR == 1 {
            for (i = 1; i <= NF; i++) { if ($i == col) { c = i } }
            if (!c) { print "NOCOL"; exit }
            next
        }
        NF == 0 { next }
        !seen  { first = $c; seen = 1 }
                 { last = $c }
        END {
            if (!c) { exit }
            if (first == 0) { print "ZEROFIRST"; exit }
            drift = (last - first) / first
            if (drift < 0) { drift = -drift }
            ok = (mode == "above") ? (drift > bound) : (drift < bound)
            printf "%s %.6e\n", (ok ? "PASS" : "FAIL"), drift
        }' "$LOG")

    case "$RESULT" in
        NOCOL)
            echo "FAIL: no column '$COLUMN' in $LOG" >&2
            status=1
            continue
            ;;
        ZEROFIRST)
            echo "FAIL: '$COLUMN' starts at zero in $LOG, so relative drift is undefined" >&2
            status=1
            continue
            ;;
    esac

    VERDICT=${RESULT% *}
    DRIFT=${RESULT#* }

    if [ "$VERDICT" = "PASS" ]; then
        echo "PASS: $COLUMN drift $DRIFT is $MODE $BOUND"
    else
        echo "FAIL: $COLUMN drift $DRIFT is not $MODE $BOUND" >&2
        status=1
    fi
done

exit $status
