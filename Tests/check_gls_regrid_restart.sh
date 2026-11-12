#!/bin/sh
# Compare every checkpoint array, including both GLS history slots and scratch.
set -eu
input=$1
exe=$2
subcycle=$3
shift 3
rm -rf continuous restarted
mkdir -p continuous restarted
cd continuous
"$@" "$exe" "$input" 'remora.n_cell=20 24 12' remora.fixed_dt=1 \
  remora.max_step=3 remora.plot_int=-1 remora.check_int=2 amr.regrid_int=-1 \
  amrex.fpe_trap_invalid=0 remora.do_substep="$subcycle" > run.log 2>&1
cd ../restarted
"$@" "$exe" "$input" 'remora.n_cell=20 24 12' remora.fixed_dt=1 \
  remora.max_step=3 remora.plot_int=-1 remora.check_int=2 amr.regrid_int=-1 \
  amrex.fpe_trap_invalid=0 remora.do_substep="$subcycle" test.skip_transfer=1 \
  remora.restart=../continuous/chk00002 > run.log 2>&1
for level in 0 1; do
  diff -r ../continuous/chk00003/Level_"$level" chk00003/Level_"$level"
done
