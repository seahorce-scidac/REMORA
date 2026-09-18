#!/bin/sh
#
# Assert how far a refinement level reaches in one direction.
#
# usage: check_level_extent.sh <pltfile> <level> <dim> <expected_lo> <expected_hi> [tol]
#
#   <level>  refinement level to measure, 1 or higher
#   <dim>    0 for x, 1 for y, 2 for z
#
# Used by add_test_extent in CTestList.cmake. check_max_level.sh says which levels exist;
# this says where they are, which is what distinguishes "the region the user asked for was
# refined" from "some part of it was". The case it was written for is a static refinement box
# drawn across a coastline: REMORA used to clip such a box back to its water part, and the
# only visible difference is the extent of the level it built.
#
# A plotfile Header lists each level as "<lev> <nboxes> <x>", a time, then three lines of
# "<lo> <hi>" physical coordinates per box (x, y, then z), and closes with "Level_<lev>/Cell".
# Level <L>'s block is therefore what lies between the "Level_<L-1>/Cell" line and the
# "Level_<L>/Cell" line.
#
# What is compared is the min of the los and the max of the his over that level's boxes, to
# within <tol> (default 1e-9, absolute). That is an extent, not a proof of coverage: boxes
# could in principle span the range and leave a hole between them. For a static box, which is
# tiled without holes, the extent is the assertion that matters.

set -eu

if [ "$#" -lt 5 ] || [ "$#" -gt 6 ]; then
    echo "usage: $0 <pltfile> <level> <dim> <expected_lo> <expected_hi> [tol]" >&2
    exit 2
fi

PLTFILE="$1"
LEVEL="$2"
DIM="$3"
EXP_LO="$4"
EXP_HI="$5"
TOL="${6:-1e-9}"

HEADER="$PLTFILE/Header"

if [ ! -f "$HEADER" ]; then
    echo "check_level_extent: no such plotfile header: $HEADER" >&2
    exit 1
fi

if [ "$LEVEL" -lt 1 ]; then
    echo "check_level_extent: level must be 1 or higher (got $LEVEL)" >&2
    exit 1
fi

PREV=$((LEVEL - 1))

awk -v lev="$LEVEL" -v prev="$PREV" -v dim="$DIM" \
    -v exp_lo="$EXP_LO" -v exp_hi="$EXP_HI" -v tol="$TOL" '
    $0 == "Level_" prev "/Cell" { inblock = 1; nline = 0; next }
    $0 == "Level_" lev  "/Cell" { inblock = 0 }
    inblock {
        nline++
        # 1: "<lev> <nboxes> <x>", 2: time, then 3 coordinate lines per box.
        if (nline == 1) { nboxes = $2; next }
        if (nline == 2) { next }
        idx = nline - 3            # 0-based index among the coordinate lines
        if (idx % 3 != dim) { next }
        if (seen++ == 0) { lo = $1; hi = $2 }
        if ($1 < lo) { lo = $1 }
        if ($2 > hi) { hi = $2 }
    }
    END {
        if (!seen) {
            printf "check_level_extent: FAIL: level %d has no boxes in the header\n", lev > "/dev/stderr"
            exit 1
        }
        dlo = lo - exp_lo; if (dlo < 0) { dlo = -dlo }
        dhi = hi - exp_hi; if (dhi < 0) { dhi = -dhi }
        if (dlo > tol || dhi > tol) {
            printf "check_level_extent: FAIL: level %d dim %d spans [%.10g, %.10g], expected [%.10g, %.10g]\n", \
                   lev, dim, lo, hi, exp_lo, exp_hi > "/dev/stderr"
            printf "  (%d boxes at that level)\n", nboxes > "/dev/stderr"
            exit 1
        }
        printf "check_level_extent: level %d dim %d spans [%.10g, %.10g] as expected\n", lev, dim, lo, hi
    }
' "$HEADER"
