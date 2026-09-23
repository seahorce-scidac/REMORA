#!/bin/sh
#
# Assert the highest refinement level a run actually created.
#
# usage: check_max_level.sh <logfile> <expected_max_level>
#
# Used by add_test_nlevels in CTestList.cmake. Nothing else in the suite looks at grid
# structure: fcompare reads the plotfile's values, and a tagging bug leaves those perfectly
# self-consistent. A criterion that has quietly stopped tagging, or one that has started
# tagging the coastline, changes only which levels exist and where -- so that is what this
# checks. An expected level of 0 asserts that no criterion tagged anything at all.
#
# REMORA prints a "GRIDS AT LEVEL n ARE" line whenever it builds or rebuilds a level's grids
# -- from MakeNewLevelFromScratch, MakeNewLevelFromCoarse and RemakeLevel -- so a regrid that
# leaves the BoxArray unchanged prints nothing and the line count is not the regrid count.
# What the highest n in the log does mean is the deepest level the run ever built, which is
# what is asserted; taking the maximum rather than the last means a level that appeared and
# later vanished still counts.
#
# A run always builds level 0, so at least one line must match. Requiring that matters: this
# script is used to assert an expected level of 0, and treating an unparsable log as 0 would
# make such an assertion pass for the wrong reason if the print were reworded or if rank
# output under MPI landed mid-line.

set -eu

if [ "$#" -ne 2 ]; then
    echo "usage: $0 <logfile> <expected_max_level>" >&2
    exit 2
fi

LOGFILE="$1"
EXPECTED="$2"

if [ ! -f "$LOGFILE" ]; then
    echo "check_max_level: no such log file: $LOGFILE" >&2
    exit 1
fi

LEVELS=$(sed -n 's/^GRIDS AT LEVEL \([0-9][0-9]*\) ARE.*/\1/p' "$LOGFILE" | sort -n)

if [ -z "$LEVELS" ]; then
    echo "check_max_level: FAIL: no 'GRIDS AT LEVEL n ARE' line in $LOGFILE" >&2
    echo "  Every run builds level 0, so this is a parse failure, not a level count." >&2
    exit 1
fi

ACTUAL=$(echo "$LEVELS" | tail -1)

if [ "$ACTUAL" -ne "$EXPECTED" ]; then
    echo "check_max_level: FAIL: expected max level $EXPECTED, run built $ACTUAL" >&2
    echo "  (from $LOGFILE)" >&2
    exit 1
fi

echo "check_max_level: max level $ACTUAL as expected"
exit 0
