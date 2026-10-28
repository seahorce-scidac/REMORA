#!/bin/bash
# git-bisect test for "does this commit still agree bitwise with the baseline" on Perlmutter.
#
# Copy this script outside the repository first (git bisect checks out old commits):
#   cp Tests/perf/perlmutter_bisect.sh $SCRATCH/bisect.sh
#   git bisect start <bad> <good>
#   git bisect run $SCRATCH/bisect.sh -A m1234 -a cpu -b <good> [-w dir] [-s "NX NY NZ"] [-t 3]
#
#   -A account   NERSC allocation (_g is added for GPU)
#   -a arch      cpu | gpu                                        [cpu]
#   -b ref       baseline revision; built and run once, its plotfile is the reference
#   -w dir       work directory                                   [$SCRATCH/remora_bisect]
#   -s "NX NY NZ" Upwelling grid                                  ["1024 2048 40"]
#   -t steps     steps per run                                    [3]
#   -q qos       Slurm QOS                                        [debug]
#   -j n         make -j                                          [16]
#
# For each commit git bisect checks out, this builds it with GNU make (incremental, same
# tmp_build_dir), submits a short job with sbatch --wait (CPU: 128 ranks x 1 thread on one
# node; GPU: 4 ranks on one node), and compares the last plotfile with the baseline's using
# fcompare. Exit 0 = identical (good), 1 = differs (bad), 125 = could not build/run (skip).
set -u

ACCOUNT=""; ARCH=cpu; BASE=""; WORK=""; SIZE="1024 2048 40"; NSTEP=3; QOS=debug; JOBS=16
while getopts "A:a:b:w:s:t:q:j:" opt; do
    case $opt in
        A) ACCOUNT=$OPTARG;; a) ARCH=$OPTARG;; b) BASE=$OPTARG;; w) WORK=$OPTARG;;
        s) SIZE=$OPTARG;; t) NSTEP=$OPTARG;; q) QOS=$OPTARG;; j) JOBS=$OPTARG;;
        *) echo "unknown option"; exit 125;;
    esac
done
[ -n "$ACCOUNT" ] && [ -n "$BASE" ] || { echo "-A and -b are required"; exit 125; }
[ -n "$WORK" ] || WORK=${SCRATCH:-$HOME}/remora_bisect
mkdir -p "$WORK"
REPO=$(git rev-parse --show-toplevel) || exit 125
SHA=$(git rev-parse --short HEAD)
PERF=$REPO/Tests/perf
[ -f "$PERF/inputs/upwelling_perf.in" ] || PERF=$(dirname "$(readlink -f "$0")")   # old commit: use the copy's dir

if [ -f "$REPO/Build/saul-env.sh" ] && command -v module > /dev/null; then
    source "$REPO/Build/saul-env.sh"
fi
if [ $ARCH = gpu ]; then
    module load cudatoolkit 2>/dev/null || true
    MAKEFLAGS_ARCH="USE_OMP=FALSE USE_CUDA=TRUE"
    SBATCH_ACCT=${ACCOUNT}_g; SBATCH_C=gpu
    SBATCH_EXTRA="#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=none"
    MGX=256; MGY=256
else
    MAKEFLAGS_ARCH="USE_OMP=TRUE USE_CUDA=FALSE"
    SBATCH_ACCT=$ACCOUNT; SBATCH_C=cpu
    SBATCH_EXTRA=""
    MGX=128; MGY=128
fi

build () {   # build <src dir> <out binary>
    (cd "$1/Exec" && make -j "$JOBS" COMP=gnu USE_MPI=TRUE USE_PNETCDF=TRUE TINY_PROFILE=TRUE \
        USE_PARTICLES=FALSE $MAKEFLAGS_ARCH > "$2.build.log" 2>&1) || return 1
    cp "$(ls -t "$1"/Exec/REMORA3d.gnu.*.ex | head -1)" "$2"
}

# ---------------------------------------------------------------- inputs (once)
read -r NX NY NZ <<< "$SIZE"
INP=$WORK/upwelling_$ARCH.in
if [ ! -f "$INP" ]; then
    DT=$(python3 -c "import math; dx=41000.0/$NX; print(round(min(50.0, 30*0.35*dx/math.sqrt(9.81*150)), 2))")
    sed -e "s/@NSTEP@/$NSTEP/g" -e "s/@NX@/$NX/g" -e "s/@NY@/$NY/g" -e "s/@NZ@/$NZ/g" \
        -e "s/@MGX@/$MGX/g" -e "s/@MGY@/$MGY/g" -e "s/@DT@/$DT/g" "$PERF/inputs/upwelling_perf.in" > "$INP"
fi

# ---------------------------------------------------------------- one run as a job
run_job () {   # run_job <binary> <out dir>
    local bin=$1 out=$2 job=$2.sbatch
    mkdir -p "$out"
    if [ $ARCH = gpu ]; then
        cat > "$WORK/gpu_wrap.sh" <<'EOF'
#!/bin/bash
export CUDA_VISIBLE_DEVICES=$((3-SLURM_LOCALID))
exec "$@"
EOF
        chmod +x "$WORK/gpu_wrap.sh"
        SRUN="srun -n 4 --cpus-per-task=32 --cpu-bind=cores $WORK/gpu_wrap.sh"
        EXTRA_ARGS="amrex.use_gpu_aware_mpi=1"
    else
        SRUN="srun -n 128 -c 2 --cpu-bind=cores"
        EXTRA_ARGS=""
    fi
    cat > "$job" <<EOF
#!/bin/bash
#SBATCH -A $SBATCH_ACCT
#SBATCH -C $SBATCH_C
#SBATCH -q $QOS
#SBATCH -N 1
#SBATCH -t 00:25:00
#SBATCH -J remora_bisect
#SBATCH -o $out/slurm.out
$SBATCH_EXTRA
export OMP_NUM_THREADS=1
export MPICH_OFI_NIC_POLICY=GPU
cd $out
$SRUN $bin $INP remora.plot_file=$out/plt $EXTRA_ARGS > $out/log 2>&1
EOF
    sbatch --wait "$job" > /dev/null 2>&1 || return 1
    ls -d "$out"/plt[0-9]* > /dev/null 2>&1
}

# ---------------------------------------------------------------- baseline (once)
REF=$WORK/ref_$ARCH
if [ ! -d "$REF/run" ] || ! ls -d "$REF"/run/plt[0-9]* > /dev/null 2>&1; then
    echo "bisect: building and running the baseline $BASE"
    [ -d "$REF/src" ] || { git worktree add --detach "$REF/src" "$BASE" > /dev/null 2>&1 &&
                           git -C "$REF/src" submodule update --init Submodules/AMReX > /dev/null 2>&1; } || exit 125
    build "$REF/src" "$REF/remora_exec" || { echo "bisect: baseline build failed"; exit 125; }
    (cd "$REF/src/Submodules/AMReX/Tools/Plotfile" && make -j "$JOBS" programs=fcompare > /dev/null 2>&1) || exit 125
    cp "$REF"/src/Submodules/AMReX/Tools/Plotfile/fcompare*.ex "$WORK/fcompare"
    run_job "$REF/remora_exec" "$REF/run" || { echo "bisect: baseline run failed"; exit 125; }
fi
REFPLT=$(ls -d "$REF"/run/plt[0-9]* | tail -1)

# ---------------------------------------------------------------- this commit
OUT=$WORK/$SHA
echo "bisect: $SHA building"
build "$REPO" "$OUT.exe" || { echo "bisect: $SHA build failed (skip)"; exit 125; }
echo "bisect: $SHA running"
run_job "$OUT.exe" "$OUT" || { echo "bisect: $SHA run failed (skip)"; exit 125; }
PLT=$(ls -d "$OUT"/plt[0-9]* | tail -1)
if "$WORK/fcompare" "$REFPLT" "$PLT" > "$OUT/fcompare.log" 2>&1 && grep -q AGREE "$OUT/fcompare.log"; then
    echo "bisect: $SHA IDENTICAL to baseline (good)"; exit 0
else
    echo "bisect: $SHA DIFFERS from baseline (bad)"; grep -E "^ [a-z]" "$OUT/fcompare.log" | head -8; exit 1
fi
