#!/bin/sh
#
# Assert how many cells a refinement level covers.
#
# usage: check_level_cells.sh <pltfile> <level> <expected_cells>
#
# Used by add_test_cells in CTestList.cmake. The other two grid assertions answer "does this
# level exist" (check_max_level.sh) and "how far does it reach" (check_level_extent.sh).
# Neither can see a change that removes cells from the middle of a refined region while
# leaving its bounding box alone, which is exactly what a mask stencil that excludes the
# water cells along a coast does. This counts them.
#
# The count is read from <pltfile>/Level_<level>/Cell_H, which records the level's BoxArray in
# index space as lines of "((ilo,jlo,klo) (ihi,jhi,khi) (..))". It is summed over boxes, so it
# is invariant under how AMReX chops that region for load balance and therefore under the
# number of MPI ranks -- unlike the box count or the BoxArray itself, which are not.

set -eu

if [ "$#" -ne 3 ]; then
    echo "usage: $0 <pltfile> <level> <expected_cells>" >&2
    exit 2
fi

PLTFILE="$1"
LEVEL="$2"
EXPECTED="$3"

CELL_H="$PLTFILE/Level_$LEVEL/Cell_H"

if [ ! -f "$CELL_H" ]; then
    # An absent level is a real answer -- zero cells -- not a broken plotfile, so long as the
    # plotfile itself is there. Distinguish the two.
    if [ ! -f "$PLTFILE/Header" ]; then
        echo "check_level_cells: no such plotfile: $PLTFILE" >&2
        exit 2
    fi
    ACTUAL=0
else
    ACTUAL=$(awk '
        # Box lines look like ((39,0,0) (62,44,15) (0,0,0)); anything else in the header is
        # skipped. gsub to spaces lets the default field splitting pull out the nine integers.
        /^\(\([0-9-]+,[0-9-]+,[0-9-]+\) \([0-9-]+,[0-9-]+,[0-9-]+\)/ {
            line = $0
            gsub(/[(),]/, " ", line)
            n = split(line, f, " ")
            if (n < 6) { next }
            total += (f[4]-f[1]+1) * (f[5]-f[2]+1) * (f[6]-f[3]+1)
            nbox++
        }
        END {
            if (!nbox) {
                print "check_level_cells: no box lines found in the header" > "/dev/stderr"
                exit 1
            }
            printf "%d\n", total
        }
    ' "$CELL_H")
fi

if [ "$ACTUAL" -ne "$EXPECTED" ]; then
    echo "check_level_cells: FAIL: level $LEVEL covers $ACTUAL cells, expected $EXPECTED" >&2
    echo "  (from $CELL_H)" >&2
    exit 1
fi

echo "check_level_cells: level $LEVEL covers $ACTUAL cells as expected"
exit 0
