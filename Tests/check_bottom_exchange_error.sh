#!/bin/sh
# An invalid boundary contract must fail with the intended diagnostic.
set -eu
expected="$1"
shift
if "$@" > bottom_error.log 2>&1; then
    cat bottom_error.log
    echo "FAIL: invalid bottom exchange was accepted"
    exit 1
fi
if grep -F "$expected" bottom_error.log; then
    exit 0
fi
cat bottom_error.log
echo "FAIL: expected bottom exchange diagnostic was absent"
exit 1
