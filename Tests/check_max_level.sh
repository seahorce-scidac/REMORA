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
# REMORA prints one "GRIDS AT LEVEL n ARE" line per level on every regrid. The highest n
# anywhere in the log is the deepest level the run ever built, which is what is asserted;
# taking the maximum rather than the last means a level that appeared and later vanished
# still counts.

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

# No match at all is level 0, not an error: a single-level run prints the level 0 line, and
# grep exiting 1 under `set -e` would otherwise abort before the comparison.
ACTUAL=$(sed -n 's/^GRIDS AT LEVEL \([0-9][0-9]*\) ARE.*/\1/p' "$LOGFILE" | sort -n | tail -1)
ACTUAL=${ACTUAL:-0}

if [ "$ACTUAL" -ne "$EXPECTED" ]; then
    echo "check_max_level: FAIL: expected max level $EXPECTED, run built $ACTUAL" >&2
    echo "  (from $LOGFILE)" >&2
    exit 1
fi

echo "check_max_level: max level $ACTUAL as expected"
exit 0
