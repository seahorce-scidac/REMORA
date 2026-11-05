#!/bin/sh
# Run in an isolated test directory. MPI launcher words precede the executable.
# Usage: check_biology_restart.sh <input> <fcompare> <gold plotfile> <launcher...> <exe>
set -eu
input=$1
compare=$2
gold=$3
shift 3
mkdir -p runA runB alias mismatch
(
    cd runA
    "$@" "$input" remora.check_mask_consistency=true > ../run.log 2>&1
)
# The answer itself, at the tolerance add_test_r uses; the checks below only
# prove the run is self-consistent.
"$compare" -a -r 1e-11 --abs_tol 1.0e-11 "$gold" runA/plt00003
(
    cd runB
    "$@" "$input" remora.restart=../runA/chk00002 remora.check_int=-1 > ../restart.log 2>&1
)
"$compare" -a -r 0 --abs_tol 0 runA/plt00003 runB/plt00003
(
    cd alias
    "$@" "$input" 'remora.plot_vars_3d=salt temp tracer x_velocity y_velocity z_velocity ecb' > ../alias.log 2>&1
)
"$compare" -a -r 0 --abs_tol 0 runA/plt00003 alias/plt00003
# Small CPU tiles split each box horizontally; the column kernel must give
# the same answer whatever the tiling.
mkdir -p tiled
(
    cd tiled
    "$@" "$input" 'remora.omp_tile_size=2 3 1024000' > ../tiled.log 2>&1
)
"$compare" -a -r 0 --abs_tol 0 runA/plt00003 tiled/plt00003
# Save the metadata so rerunning this test leaves the original checkpoint intact.
cp runA/chk00002/BiologyHeader saved_biology_header
cp runA/chk00002/Header saved_header
restore_header() {
    cp saved_biology_header runA/chk00002/BiologyHeader
    cp saved_header runA/chk00002/Header
}
trap restore_header EXIT
for mutation in identity ordering missing; do
    case "$mutation" in
        identity)
            sed '2s/.*/none/' saved_biology_header > runA/chk00002/BiologyHeader
            message='Checkpoint biology identity differs'
            ;;
        ordering)
            # Same count, different order: the generic Header layout check.
            sed '/^cons_names /s/ temp salt / salt temp /' saved_header > runA/chk00002/Header
            if cmp -s saved_header runA/chk00002/Header; then
                echo "Header has no cons_names line to reorder" >&2
                exit 1
            fi
            message='Checkpoint cell-centered components differ'
            ;;
        missing)
            mv runA/chk00002/BiologyHeader missing_biology_header
            message='Selected biology requires checkpoint identity metadata'
            ;;
    esac
    if (cd mismatch && "$@" "$input" remora.restart=../runA/chk00002 > "../$mutation.log" 2>&1); then
        echo "Restart incorrectly accepted $mutation metadata" >&2
        exit 1
    fi
    if ! grep -q "$message" "$mutation.log"; then
        cat "$mutation.log" >&2
        exit 1
    fi
    restore_header
 done
