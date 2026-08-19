#!/bin/sh
#
# Compare variables' min and max in a plotfile against expected values.
#
# usage: check_extrema.sh <fextrema_exe> <pltfile> <tol> <var> <min> <max> [<var> <min> <max> ...]
#
# Used by add_test_extrema in CTestList.cmake. Unlike an fcompare gold file, this asserts what
# the value should BE -- from a closed-form reference, or from a constant that an initial
# condition is required to reproduce -- so it also catches a baseline that was wrong when it
# was blessed. Takes any number of variables so one model run can carry several assertions.
#
# amrex_fextrema prints one "<name> <min> <max>" line per variable at 11 significant digits.
# Comparison is absolute; pick <tol> to suit the magnitudes being checked. Every mismatch is
# reported before exiting, so one run tells you about all of them.

set -eu

if [ "$#" -lt 6 ]; then
    echo "usage: $0 <fextrema_exe> <pltfile> <tol> <var> <min> <max> [<var> <min> <max> ...]" >&2
    exit 2
fi

FEXTREMA="$1"; shift
PLTFILE="$1"; shift
TOL="$1"; shift

if [ ! -x "$FEXTREMA" ]; then
    echo "FAIL: fextrema not found or not executable: $FEXTREMA" >&2
    exit 1
fi
if [ ! -d "$PLTFILE" ]; then
    echo "FAIL: plotfile not found: $PLTFILE" >&2
    exit 1
fi

# One fextrema call for all variables; the plotfile is read once.
EXTREMA=$("$FEXTREMA" "$PLTFILE")

status=0
while [ "$#" -ge 3 ]; do
    var="$1"; want_min="$2"; want_max="$3"; shift 3
    line=$(printf '%s\n' "$EXTREMA" | awk -v v="$var" '$1 == v { print; exit }')
    if [ -z "$line" ]; then
        echo "FAIL: variable $var not present in $PLTFILE"
        status=1
        continue
    fi
    if printf '%s\n' "$line" | awk -v v="$var" -v mn="$want_min" -v mx="$want_max" -v tol="$TOL" '
        function abs(x) { return x < 0 ? -x : x }
        {
            dmin = abs($2 - mn); dmax = abs($3 - mx)
            if (dmin > tol || dmax > tol) {
                printf "FAIL %s: got [%.11g, %.11g], want [%.11g, %.11g], tol %g (dmin %.3g, dmax %.3g)\n",
                       v, $2, $3, mn, mx, tol, dmin, dmax
                exit 1
            }
            printf "OK %s: [%.11g, %.11g] within %g of [%.11g, %.11g]\n", v, $2, $3, tol, mn, mx
        }'
    then
        :
    else
        status=1
    fi
done

if [ "$#" -ne 0 ]; then
    echo "FAIL: trailing arguments; variables must come in <var> <min> <max> triples" >&2
    exit 2
fi

exit "$status"
