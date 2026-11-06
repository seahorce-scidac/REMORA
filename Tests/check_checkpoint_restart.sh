#!/bin/sh
# Checkpoint restart and mismatch checks. Run in an isolated test directory.
# Usage: check_checkpoint_restart.sh <input> <fcompare> <checkpoint> <plotfile>
#            <missing BiologyHeader: accept|reject> <launcher...> <exe>
# MPI launcher words precede the executable.
set -eu
input=$1
compare=$2
chk=$3
plt=$4
missing_header=$5
shift 5
mkdir -p runA runB mismatch

# A restart from the checkpoint must reproduce the continuous run exactly.
(
    cd runA
    "$@" "$input" > ../run.log 2>&1
)
(
    cd runB
    "$@" "$input" remora.restart="../runA/$chk" remora.check_int=-1 > ../restart.log 2>&1
)
"$compare" -a -r 0 --abs_tol 0 "runA/$plt" "runB/$plt"

# Save the metadata so rerunning this test leaves the original checkpoint intact.
cp "runA/$chk/BiologyHeader" saved_biology_header
cp "runA/$chk/Header" saved_header
restore_header() {
    cp saved_biology_header "runA/$chk/BiologyHeader"
    cp saved_header "runA/$chk/Header"
}
trap restore_header EXIT

# Each mutation must be rejected with its message, except a missing
# BiologyHeader, which models that accept legacy checkpoints restart from.
for mutation in identity ordering missing; do
    expect=reject
    case "$mutation" in
        identity)
            sed '2s/.*/none/' saved_biology_header > "runA/$chk/BiologyHeader"
            message='Checkpoint biology identity differs'
            ;;
        ordering)
            # Same count, different order: the Header layout check.
            sed '/^cons_names /s/ temp salt / salt temp /' saved_header > "runA/$chk/Header"
            if cmp -s saved_header "runA/$chk/Header"; then
                echo "Header has no cons_names line to reorder" >&2
                exit 1
            fi
            message='Checkpoint cell-centered components differ'
            ;;
        missing)
            mv "runA/$chk/BiologyHeader" missing_biology_header
            message='Selected biology requires checkpoint identity metadata'
            expect=$missing_header
            ;;
    esac
    if (cd mismatch && "$@" "$input" remora.restart="../runA/$chk" remora.check_int=-1 > "../$mutation.log" 2>&1); then
        if [ "$expect" = reject ]; then
            echo "Restart incorrectly accepted $mutation metadata" >&2
            exit 1
        fi
    else
        if [ "$expect" = accept ]; then
            cat "$mutation.log" >&2
            echo "Restart incorrectly rejected $mutation metadata" >&2
            exit 1
        fi
        if ! grep -q "$message" "$mutation.log"; then
            cat "$mutation.log" >&2
            exit 1
        fi
    fi
    restore_header
done
