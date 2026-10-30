#!/bin/bash
# Phase 3 (land layouts) test plan on Perlmutter, end to end: builds, the masked grid, the
# CPU and GPU jobs, the answer checks, and a summary table.
#
#   source Build/saul-env.sh
#   Tests/perf/perlmutter_phase3.sh -A m1234 [-d $SCRATCH/REMORA-data] [-w dir] [-f 128]
#                                   [-s "A B C D"] [-b 64] [-j 16] [-q debug] [-S]
#   Tests/perf/perlmutter_phase3.sh -w dir -R        # print the summary of a finished run
#
#   -A account   NERSC allocation (_g is added for GPU jobs)
#   -d dir       REMORA-data checkout (must contain IdealMiniGrid)  [$SCRATCH/REMORA-data]
#   -w dir       work directory; reuse it to resume                 [$SCRATCH/remora_p3]
#   -f factor    grid factor: (10f) x (16f) x 20 cells, 36 % land   [128]
#   -s stages    which tests to submit                             ["A B C D"]
#                A  dropped vs full answers, bitwise (CPU 128x1 and 16x8, GPU)
#                B  CPU time, full vs dropped, boxes 32/64/128 x mixes 128x1/32x4/16x8
#                C  CPU land weight 1/0.5/0 with dropped boxes (box size -b)
#                D  GPU time, boxes 256/128/64, full vs dropped
#   -b n         box size for stage C                               [64]
#   -j n         make -j                                            [16]
#   -q qos       Slurm QOS (each job is under 30 minutes)           [debug]
#   -S           build and write the job scripts, do not submit
#   -R           only print the summary of the work directory
#
# The builds are CMake builds with the Cray wrappers (cc/CC), PnetCDF, TinyProfiler and the
# remora_subset_compare tool; they go in <work>/build_cpu and <work>/build_gpu. Every run is
# in <work>/runs/<name>/ with its log; answer checks are in <work>/runs/cmp_*.log. Runs
# that finished are skipped when a job is resubmitted, so a timed-out job can just be rerun.
# Stage E of PERLMUTTER_PHASE3_TEST_PLAN.md uses Tests/perf/perlmutter_phase1.sh; the
# command is printed at the end.
set -u

ACCOUNT=""; DATA=${SCRATCH:-$HOME}/REMORA-data; WORK=${SCRATCH:-$HOME}/remora_p3; FACTOR=128
STAGES="A B C D"; CBOX=64; JOBS=16; QOS=debug; SUBMIT=1; SUMMARY_ONLY=0
while getopts "A:d:w:f:s:b:j:q:SR" opt; do
    case $opt in
        A) ACCOUNT=$OPTARG;; d) DATA=$OPTARG;; w) WORK=$OPTARG;; f) FACTOR=$OPTARG;;
        s) STAGES=$OPTARG;; b) CBOX=$OPTARG;; j) JOBS=$OPTARG;; q) QOS=$OPTARG;;
        S) SUBMIT=0;; R) SUMMARY_ONLY=1;;
        *) echo "unknown option"; exit 1;;
    esac
done
mkdir -p "$WORK" && WORK=$(cd "$WORK" && pwd)
RUNS=$WORK/runs

# ---------------------------------------------------------------- summary
summary () {
    echo "== $WORK"
    printf '%-34s %10s  %s\n' run "Evolve s" "level-0 boxes kept"
    for d in "$RUNS"/*/; do
        local n log ev kept
        n=$(basename "$d"); log=$d/log
        [ -f "$log" ] || continue
        ev=$(awk '/Incl\. Min/ {incl=1} incl && v == "" && /^REMORA::Evolve\(\)/ {v=$(NF-1)}
                  END {print (v == "" ? "--" : v)}' "$log")
        kept=$(grep -m1 -o "keeps [0-9]* of [0-9]* boxes" "$log" || true)
        grep -q "finalized" "$log" || kept="$kept (did not finish)"
        printf '%-34s %10s  %s\n' "$n" "$ev" "${kept:--}"
    done
    echo
    echo "answer checks (dropped vs full, every kept cell, zero tolerance):"
    local any=0
    for c in "$RUNS"/cmp_*.log; do
        [ -f "$c" ] || continue; any=1
        printf '  %-40s %s\n' "$(basename "$c" .log)" "$(grep -m1 -oE 'SUBSET (AGREE|DIFFER)' "$c" || echo 'no verdict, see log')"
    done
    [ $any = 1 ] || echo "  none yet"
}
if [ $SUMMARY_ONLY = 1 ]; then summary; exit 0; fi

[ -n "$ACCOUNT" ] || { echo "-A account is required"; exit 1; }
[ -d "$DATA/IdealMiniGrid" ] || { echo "-d: no IdealMiniGrid in $DATA"; exit 1; }
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(git -C "$HERE" rev-parse --show-toplevel)
echo "testing $(git -C "$REPO" log -1 --format='%h %ad %s' --date=short)"
if [ -f "$REPO/Build/saul-env.sh" ] && command -v module > /dev/null; then
    source "$REPO/Build/saul-env.sh"
fi
# REMORA's CMake finds PnetCDF with pkg-config; the Cray module sets PNETCDF_DIR.
if [ -n "${PNETCDF_DIR:-}" ]; then
    export PKG_CONFIG_PATH=$PNETCDF_DIR/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}
fi

# ---------------------------------------------------------------- builds
COMMON="-DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DCMAKE_BUILD_TYPE=Release
        -DREMORA_ENABLE_MPI=ON -DREMORA_ENABLE_PNETCDF=ON -DREMORA_ENABLE_TESTS=ON
        -DAMReX_TINY_PROFILE=ON"
build () {   # build <cpu|gpu>
    local b=$WORK/build_$1 extra targets
    if [ $1 = gpu ]; then
        module load cudatoolkit 2>/dev/null || true
        extra="-DREMORA_ENABLE_CUDA=ON -DAMReX_CUDA_ARCH=8.0 -DCMAKE_CUDA_HOST_COMPILER=CC"
        targets="remora_exec"
    else
        extra="-DREMORA_ENABLE_OPENMP=ON"
        targets="remora_exec remora_subset_compare"
    fi
    if [ ! -f "$b/CMakeCache.txt" ]; then
        echo "configuring $1 build in $b"
        cmake -S "$REPO" -B "$b" $COMMON $extra > "$WORK/configure_$1.log" 2>&1 \
            || { echo "configure ($1) failed, see $WORK/configure_$1.log"; tail -5 "$WORK/configure_$1.log"; exit 1; }
    fi
    echo "building $1 (log: $WORK/build_$1.log)"
    cmake --build "$b" -j "$JOBS" --target $targets > "$WORK/build_$1.log" 2>&1 \
        || { echo "build ($1) failed, see $WORK/build_$1.log"; tail -15 "$WORK/build_$1.log"; exit 1; }
}
need_gpu=0; for s in $STAGES; do case $s in A|D) need_gpu=1;; esac; done
build cpu
[ $need_gpu = 1 ] && build gpu
X=$WORK/build_cpu/Exec/remora_exec
XG=$WORK/build_gpu/Exec/remora_exec
CMP=$WORK/build_cpu/Tests/remora_subset_compare
for f in "$X" "$CMP"; do [ -x "$f" ] || { echo "missing $f"; exit 1; }; done
[ $need_gpu = 1 ] && { [ -x "$XG" ] || { echo "missing $XG"; exit 1; }; }

# ---------------------------------------------------------------- grid
GRID=$WORK/grid_f$FACTOR
if [ ! -s "$GRID/grid.nc" ] || [ ! -s "$GRID/ini.nc" ]; then
    module load python 2>/dev/null || true
    echo "generating the factor-$FACTOR masked grid in $GRID"
    python3 "$REPO/Tests/perf/make_masked_grid.py" "$DATA/IdealMiniGrid" "$GRID" "$FACTOR" \
        || { echo "grid generation failed (needs numpy and netCDF4)"; exit 1; }
fi
NX=$((10*FACTOR)); NY=$((16*FACTOR))
IN=$REPO/Tests/perf/inputs/masked_grid_perf.in
mkdir -p "$RUNS"

# ---------------------------------------------------------------- job scripts
cat > "$WORK/gpu_wrap.sh" <<'EOF'
#!/bin/bash
export CUDA_VISIBLE_DEVICES=$((3-SLURM_LOCALID))
exec "$@"
EOF
chmod +x "$WORK/gpu_wrap.sh"

# Body shared by every job: run <name> <launcher...> -- <exe> <args...>, skipped when its
# log already ends cleanly; check <name> <subset run> <full run> for the answer checks.
common_body () {
cat <<EOF
set -u
export MPICH_OFI_NIC_POLICY=GPU
RUNS=$RUNS; GRID=$GRID; IN=$IN; CMP=$CMP
SIZE="remora.n_cell=$NX $NY 20"
run () {
    local name=\$1; shift
    local launch=(); while [ "\$1" != "--" ]; do launch+=("\$1"); shift; done; shift
    local exe=\$1; shift   # the inputs file must be the first argument of the executable
    local d=\$RUNS/\$name
    if [ -f \$d/log ] && grep -q finalized \$d/log; then echo "=== \$name: done before, skipped"; return; fi
    mkdir -p \$d && cd \$d
    ln -sf \$GRID/grid.nc grid.nc; ln -sf \$GRID/ini.nc ini.nc
    echo "=== \$name  \$(date +%T)"
    "\${launch[@]}" \$exe \$IN \$SIZE "\$@" > log 2>&1 || echo "   RUN FAILED: \$d/log"
    awk '/Incl\\. Min/ {incl=1} incl && v == "" && /^REMORA::Evolve\\(\\)/ {v=\$(NF-1)}
         END {print "   Evolve " (v == "" ? "n/a" : v " s")}' log
    grep -m1 -o "keeps [0-9]* of [0-9]* boxes" log | sed 's/^/   /' || true
}
check () {   # check <label> <dropped run> <full run> [--chk]
    local lab=\$1 a=\$RUNS/\$2 b=\$RUNS/\$3 mode=\${4:-}
    local pa pb
    if [ "\$mode" = --chk ]; then pa=\$(ls -d \$a/chk[0-9]* 2>/dev/null | tail -1); pb=\$(ls -d \$b/chk[0-9]* 2>/dev/null | tail -1)
    else pa=\$(ls -d \$a/plt[0-9]* 2>/dev/null | tail -1); pb=\$(ls -d \$b/plt[0-9]* 2>/dev/null | tail -1); fi
    [ -n "\$pa" ] && [ -n "\$pb" ] || { echo "   \$lab: missing output" | tee \$RUNS/cmp_\$lab.log; return; }
    srun -n 4 -c 2 \$CMP \$mode \$pa \$pb > \$RUNS/cmp_\$lab.log 2>&1
    echo "   \$lab: \$(grep -m1 -oE 'SUBSET (AGREE|DIFFER)' \$RUNS/cmp_\$lab.log || echo 'no verdict')"
}
EOF
}

# Launchers: CPU ranks x threads on one node (2 hardware threads per core), and 4 GPUs.
cpu_launch () {   # cpu_launch <ranks> <threads>
    echo "env OMP_NUM_THREADS=$2 OMP_PROC_BIND=spread OMP_PLACES=cores srun -n $1 -c $((256/$1)) --cpu-bind=cores"
}
GPU_LAUNCH="srun -n 4 --cpus-per-task=32 --cpu-bind=cores $WORK/gpu_wrap.sh"
OUT="remora.plot_int=20 remora.check_int=20"
TILE="remora.omp_tile_size=16 16 1024000"   # pinned so tile defaults cannot confuse the comparison

header () {   # header <cpu|gpu> <name>
    local acct=$ACCOUNT extra=""
    [ $1 = gpu ] && { acct=${ACCOUNT}_g; extra="#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=none"; }
    cat <<EOF
#!/bin/bash
#SBATCH -A $acct
#SBATCH -C $1
#SBATCH -q $QOS
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -J p3_$2
#SBATCH -o $WORK/job_$2.out
$extra
EOF
    [ $1 = gpu ] && echo "module load cudatoolkit 2>/dev/null || true"
    common_body
}

# write_job runs at the end of a pipeline, i.e. in a subshell, so the list of jobs written by
# this invocation is kept in a file.
JOBLIST=$WORK/.jobs_this_run; : > "$JOBLIST"
write_job () {   # write_job <cpu|gpu> <name> ; body on stdin
    local f=$WORK/job_$2.sbatch
    { header $1 $2; cat; } > "$f"
    echo "$f" >> "$JOBLIST"
}

for s in $STAGES; do case $s in
A)  # answers; also one dropped run with signalling NaNs and the invalid trap on
    { for mix in "128 1" "16 8"; do
        set -- $mix; L=$(cpu_launch $1 $2); t=${1}x$2
        echo "run A_cpu_${t}_full $L -- $X \"$TILE\" $OUT amr.max_grid_size_x=128 amr.max_grid_size_y=128"
        echo "run A_cpu_${t}_drop $L -- $X \"$TILE\" $OUT amr.max_grid_size_x=128 amr.max_grid_size_y=128 remora.drop_dry_boxes=1"
        echo "check A_cpu_${t}_plt A_cpu_${t}_drop A_cpu_${t}_full"
        echo "check A_cpu_${t}_chk A_cpu_${t}_drop A_cpu_${t}_full --chk"
      done
      L=$(cpu_launch 128 1)
      echo "run A_cpu_128x1_drop_snan $L -- $X \"$TILE\" amr.max_grid_size_x=128 amr.max_grid_size_y=128 remora.drop_dry_boxes=1 amrex.init_snan=1 amrex.fpe_trap_invalid=1"
    } | write_job cpu A_cpu
    { echo "run A_gpu_full $GPU_LAUNCH -- $XG $OUT amr.max_grid_size_x=128 amr.max_grid_size_y=128 amrex.use_gpu_aware_mpi=1"
      echo "run A_gpu_drop $GPU_LAUNCH -- $XG $OUT amr.max_grid_size_x=128 amr.max_grid_size_y=128 amrex.use_gpu_aware_mpi=1 remora.drop_dry_boxes=1"
      echo "check A_gpu_plt A_gpu_drop A_gpu_full"
      echo "check A_gpu_chk A_gpu_drop A_gpu_full --chk"
    } | write_job gpu A_gpu
    ;;
B)  # CPU time: one job per rank/thread mix, 6 runs each
    for mix in "128 1" "32 4" "16 8"; do
        set -- $mix; L=$(cpu_launch $1 $2); t=${1}x$2
        { for b in 32 64 128; do
            echo "run B_${t}_b${b}_full $L -- $X \"$TILE\" amr.max_grid_size_x=$b amr.max_grid_size_y=$b"
            echo "run B_${t}_b${b}_drop $L -- $X \"$TILE\" amr.max_grid_size_x=$b amr.max_grid_size_y=$b remora.drop_dry_boxes=1"
          done; } | write_job cpu B_$t
    done
    ;;
C)  # land weighting with dropped boxes
    { for mix in "128 1" "16 8"; do
        set -- $mix; L=$(cpu_launch $1 $2); t=${1}x$2
        for w in 1 0.5 0; do
            echo "run C_${t}_b${CBOX}_w$w $L -- $X \"$TILE\" amr.max_grid_size_x=$CBOX amr.max_grid_size_y=$CBOX remora.drop_dry_boxes=1 remora.load_balance_land_weight=$w"
        done
      done; } | write_job cpu C
    ;;
D)  # GPU: large boxes full vs smaller boxes with and without dropping
    { echo "run D_gpu_b256_full $GPU_LAUNCH -- $XG amr.max_grid_size_x=256 amr.max_grid_size_y=256 amrex.use_gpu_aware_mpi=1"
      for b in 128 64; do
        echo "run D_gpu_b${b}_full $GPU_LAUNCH -- $XG amr.max_grid_size_x=$b amr.max_grid_size_y=$b amrex.use_gpu_aware_mpi=1"
        echo "run D_gpu_b${b}_drop $GPU_LAUNCH -- $XG amr.max_grid_size_x=$b amr.max_grid_size_y=$b amrex.use_gpu_aware_mpi=1 remora.drop_dry_boxes=1"
      done; } | write_job gpu D
    ;;
*) echo "unknown stage $s"; exit 1;;
esac; done

mapfile -t JOBFILES < "$JOBLIST"
echo "job scripts:"; printf '  %s\n' "${JOBFILES[@]}"
if [ $SUBMIT = 1 ]; then
    for f in "${JOBFILES[@]}"; do echo "  $(sbatch "$f")  $(basename "$f")"; done
    echo
    echo "When the jobs are done:  $0 -w $WORK -R"
else
    echo "not submitted (-S); submit with: for f in $WORK/job_*.sbatch; do sbatch \$f; done"
fi
echo
echo "Stage E (Upwelling, Phase 1-3 at 16x16 tiles, all mixes) runs with the existing script:"
echo "  Tests/perf/perlmutter_phase1.sh -A $ACCOUNT -b c8bc0b8 \\"
echo "    -r \"base:c8bc0b8 p1:perf-phase1-omp p2:perf-phase2 p3:perf-phase3\" \\"
echo "    -s \"1024 2048 40\" -t 10 -c \"128x1 64x2 32x4 16x8\" -T \"16x16\""
