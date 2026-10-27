#!/bin/bash
# Build two REMORA revisions on Perlmutter and time them against each other on CPU and GPU
# nodes, checking that their answers agree bitwise.
#
#   Tests/perf/perlmutter_phase1.sh -A m1234 [options]
#
#   -A account     NERSC allocation (the _g suffix is added for GPU jobs)
#   -a arch        cpu | gpu | both                         [both]
#   -b ref         baseline git ref                         [merge-base with upstream/development]
#   -n ref         new git ref                              [HEAD]
#   -w dir         work directory for sources, binaries and runs
#                                                           [$SCRATCH/remora_perf/<date>]
#   -g dir         directory holding the GS_toy*.nc Gulf Stream files; the nested case is
#                  skipped when not given
#   -s "NX NY NZ"  Upwelling grid                           ["1024 2048 40"]
#   -t steps       steps per run                            [20]
#   -q qos         Slurm QOS                                [regular]
#   -j n           make -j                                  [16]
#   -S             build only, do not submit
#
# Each job runs the baseline and the new binary back to back on the same node(s), one
# configuration at a time, so the pair always shares hardware. CPU jobs sweep ranks x threads
# on one 128-core node; the baseline is also run with the new tile default so the tiling
# and the threading effects separate. GPU jobs use one node with four GPUs. After the runs
# the job compares the last plotfiles with amrex fcompare. Summarise with
#   python3 Tests/perf/summarize.py <workdir>/runs
set -euo pipefail

ACCOUNT=""; ARCH=both; BASE=""; NEW=HEAD; WORK=""; GULF=""; SIZE="1024 2048 40"
NSTEP=20; QOS=regular; JOBS=16; SUBMIT=1
while getopts "A:a:b:n:w:g:s:t:q:j:S" opt; do
    case $opt in
        A) ACCOUNT=$OPTARG;; a) ARCH=$OPTARG;; b) BASE=$OPTARG;; n) NEW=$OPTARG;;
        w) WORK=$OPTARG;; g) GULF=$OPTARG;; s) SIZE=$OPTARG;; t) NSTEP=$OPTARG;;
        q) QOS=$OPTARG;; j) JOBS=$OPTARG;; S) SUBMIT=0;;
        *) echo "unknown option"; exit 1;;
    esac
done
[ -n "$ACCOUNT" ] || { echo "-A account is required"; exit 1; }

HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(git -C "$HERE" rev-parse --show-toplevel)
[ -n "$BASE" ] || BASE=$(git -C "$REPO" merge-base HEAD upstream/development 2>/dev/null \
                        || git -C "$REPO" merge-base HEAD origin/development)
[ -n "$WORK" ] || WORK=${SCRATCH:-$HOME}/remora_perf/$(date +%Y%m%d_%H%M)
mkdir -p "$WORK/bin" "$WORK/runs"
echo "baseline $BASE  new $NEW  work $WORK"

# ---------------------------------------------------------------- environment and builds
if [ -f "$REPO/Build/saul-env.sh" ] && command -v module > /dev/null; then
    source "$REPO/Build/saul-env.sh"
fi
build_one () {   # build_one <tag> <ref> <cpu|gpu>
    local tag=$1 ref=$2 arch=$3 src=$WORK/src_$1 flags
    if [ ! -d "$src" ]; then
        git -C "$REPO" worktree add --detach "$src" "$ref"
        git -C "$src" submodule update --init Submodules/AMReX
    fi
    if [ $arch = gpu ]; then
        module load cudatoolkit 2>/dev/null || true
        flags="USE_OMP=FALSE USE_CUDA=TRUE"
    else
        flags="USE_OMP=TRUE USE_CUDA=FALSE"
    fi
    (cd "$src/Exec" && make -j "$JOBS" COMP=gnu USE_MPI=TRUE USE_PNETCDF=TRUE TINY_PROFILE=TRUE \
        USE_PARTICLES=FALSE $flags > "$WORK/build_${tag}_${arch}.log" 2>&1) \
        || { echo "build of $tag ($arch) failed, see $WORK/build_${tag}_${arch}.log"; exit 1; }
    local exe
    exe=$(ls -t "$src"/Exec/REMORA3d.gnu.*.ex | head -1)
    cp "$exe" "$WORK/bin/remora_${tag}_${arch}"
    echo "built $WORK/bin/remora_${tag}_${arch} from $(git -C "$src" rev-parse --short HEAD)"
}
build_fcompare () {
    [ -x "$WORK/bin/fcompare" ] && return
    local src=$WORK/src_new
    (cd "$src/Submodules/AMReX/Tools/Plotfile" && make -j "$JOBS" programs=fcompare \
        > "$WORK/build_fcompare.log" 2>&1)
    cp "$src"/Submodules/AMReX/Tools/Plotfile/fcompare*.ex "$WORK/bin/fcompare"
}
for arch in cpu gpu; do
    if [ $ARCH = both ] || [ $ARCH = $arch ]; then
        build_one base "$BASE" $arch
        build_one new  "$NEW"  $arch
    fi
done
build_fcompare

# ---------------------------------------------------------------- inputs
read -r NX NY NZ <<< "$SIZE"
# Upwelling's domain is 41 km wide: keep dtfast = dt/30 below dx / sqrt(g*150 m)
DT=$(python3 -c "import math; dx=41000.0/$NX; print(round(min(50.0, 0.9*30*dx/math.sqrt(9.81*150)), 2))")
fill () {   # fill <template> <out> MGX MGY
    sed -e "s/@NSTEP@/$NSTEP/g" -e "s/@NX@/$NX/g" -e "s/@NY@/$NY/g" -e "s/@NZ@/$NZ/g" \
        -e "s/@MGX@/$3/g" -e "s/@MGY@/$4/g" -e "s/@DT@/$DT/g" "$1" > "$2"
}
fill "$HERE/inputs/upwelling_perf.in" "$WORK/runs/upwelling_cpu.in" 128 128
fill "$HERE/inputs/upwelling_perf.in" "$WORK/runs/upwelling_gpu.in" 256 256
if [ -n "$GULF" ]; then
    sed -e "s/@NSTEP@/$((NSTEP*5))/g" "$HERE/inputs/gulf_perf.in" > "$WORK/runs/gulf.in"
fi

# ---------------------------------------------------------------- job scripts
# run_pair <arch> <case> <inputs> <config tag> <srun prefix> [extra inputs args]
# is defined inside the job scripts: it runs base, base with the new tile default (cpu
# only) and new, then compares the last plotfiles of base and new.
common_body () {
cat <<EOF
FC=$WORK/bin/fcompare
run_pair () {
    local arch=\$1 case=\$2 inp=\$3 cfg=\$4 srun_prefix=\$5; shift 5
    local variants="base new"; [ \$arch = cpu ] && variants="base base_j16 new"
    for v in \$variants; do
        local bin=$WORK/bin/remora_\${v%%_*}_\$arch
        local extra=""; [ \$v = base_j16 ] && extra="remora.omp_tile_size=1024000 16 1024000"
        local out=$WORK/runs/\${arch}_\${case}_\${cfg}_\$v; mkdir -p \$out
        echo "=== \$arch \$case \$cfg \$v  \$(date +%T)"
        \$srun_prefix \$bin \$inp remora.plot_file=\$out/plt \$extra "\$@" > \$out/log 2>&1 \\
            || echo "RUN FAILED: \$out"
        grep -m1 "TinyProfiler total time" \$out/log | sed "s/^/   /"
        awk '/^REMORA::Evolve\\(\\)/ {v=\$(NF-1)} END {print "   Evolve " v " s"}' \$out/log
    done
    local a=\$(ls -d $WORK/runs/\${arch}_\${case}_\${cfg}_base/plt* 2>/dev/null | tail -1)
    local b=\$(ls -d $WORK/runs/\${arch}_\${case}_\${cfg}_new/plt*  2>/dev/null | tail -1)
    if [ -n "\$a" ] && [ -n "\$b" ]; then
        \$FC \$a \$b > $WORK/runs/\${arch}_\${case}_\${cfg}_fcompare.log 2>&1 \\
            && echo "   base vs new: \$(grep -q AGREE $WORK/runs/\${arch}_\${case}_\${cfg}_fcompare.log && echo IDENTICAL || echo DIFFER)" \\
            || echo "   base vs new: DIFFER (see fcompare log)"
    fi
}
EOF
}

cpu_job=$WORK/cpu.sbatch
cat > "$cpu_job" <<EOF
#!/bin/bash
#SBATCH -A $ACCOUNT
#SBATCH -C cpu
#SBATCH -q $QOS
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH -J remora_perf_cpu
#SBATCH -o $WORK/runs/cpu_%j.out
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
cd $WORK/runs
$(common_body)
# ranks x threads on one 128-core node (two threads per core are left to the OS)
for cfg in 128x1 64x2 32x4 16x8 8x16; do
    R=\${cfg%x*}; T=\${cfg#*x}
    export OMP_NUM_THREADS=\$T
    run_pair cpu upwelling $WORK/runs/upwelling_cpu.in \$cfg "srun -n \$R -c \$((256/R)) --cpu-bind=cores"
done
EOF
if [ -n "$GULF" ]; then
cat >> "$cpu_job" <<EOF
cd $GULF
for cfg in 8x1 8x4; do
    R=\${cfg%x*}; T=\${cfg#*x}
    export OMP_NUM_THREADS=\$T
    run_pair cpu gulf $WORK/runs/gulf.in \$cfg "srun -n \$R -c \$((256/R)) --cpu-bind=cores"
done
EOF
fi

# one rank per GPU; CUDA_VISIBLE_DEVICES inverse to the local task id for GPU-aware MPI
cat > "$WORK/gpu_wrap.sh" <<'EOF'
#!/bin/bash
export CUDA_VISIBLE_DEVICES=$((3-SLURM_LOCALID))
exec "$@"
EOF
chmod +x "$WORK/gpu_wrap.sh"

gpu_job=$WORK/gpu.sbatch
cat > "$gpu_job" <<EOF
#!/bin/bash
#SBATCH -A ${ACCOUNT}_g
#SBATCH -C gpu
#SBATCH -q $QOS
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=none
#SBATCH -t 00:40:00
#SBATCH -J remora_perf_gpu
#SBATCH -o $WORK/runs/gpu_%j.out
export MPICH_OFI_NIC_POLICY=GPU
cd $WORK/runs
$(common_body)
SRUN="srun -n 4 --cpus-per-task=32 --cpu-bind=cores $WORK/gpu_wrap.sh"
run_pair gpu upwelling $WORK/runs/upwelling_gpu.in 4gpu "\$SRUN" amrex.use_gpu_aware_mpi=1
EOF
if [ -n "$GULF" ]; then
cat >> "$gpu_job" <<EOF
cd $GULF
run_pair gpu gulf $WORK/runs/gulf.in 4gpu "\$SRUN" amrex.use_gpu_aware_mpi=1
EOF
fi

echo "job scripts: $cpu_job $gpu_job"
if [ $SUBMIT = 1 ]; then
    { [ $ARCH = both ] || [ $ARCH = cpu ]; } && sbatch "$cpu_job"
    { [ $ARCH = both ] || [ $ARCH = gpu ]; } && sbatch "$gpu_job"
    echo "when the jobs finish: python3 $HERE/summarize.py $WORK/runs"
fi
