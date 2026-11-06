#!/bin/bash
# Build several REMORA revisions on Perlmutter and time them against each other on CPU and
# GPU nodes, checking that their answers agree bitwise.
#
#   Tests/perf/perlmutter_phase1.sh -A m1234 [options]
#
#   -A account     NERSC allocation (the _g suffix is added for GPU jobs)
#   -a arch        cpu | gpu | both                         [both]
#   -r "L:ref ..." revisions to compare as label:ref pairs; the first is the baseline
#                  that the others are compared with, e.g. "base:c8bc0b8 p1:perf-phase1-omp
#                  p2:perf-phase2"                          ["base:<merge-base> new:HEAD"]
#   -b ref         baseline ref (shorthand for the first pair, label "base")
#   -n ref         new ref (shorthand for a second pair, label "new")
#   -w dir         work directory for job scripts, binaries and runs
#   -B dir         build directory for the source trees and the built binaries, shared
#                  between work directories so a revision built once is reused  [= -w]
#                                                           [$SCRATCH/remora_perf/<date>]
#   -g dir         directory holding the GS_toy*.nc Gulf Stream files; the nested case is
#                  skipped when not given
#   -s "NX NY NZ"  Upwelling grid                           ["512 1024 40"]
#   -t steps       steps per run                            [10]
#   -c "RxT ..."   CPU ranks x threads configurations        ["128x1 64x2 32x4"]
#   -T "spec ..."  extra CPU tile sizes for every revision      []
#   -F mode        -ffp-contract for every CPU build (off|on|fast) [off]
#   -x "k=v ..."   extra inputs added to every run, e.g. a mixing scheme  []
#   -N             also profile each GPU revision with Nsight Systems (one extra run each)
#                  (-x "$(cat Tests/perf/inputs/gls.in)" turns GLS on)
#   -q qos         Slurm QOS                                [debug]
#   -j n           make -j                                  [16]
#   -S             build only, do not submit
#
# Each job runs every revision back to back on the same node(s), one configuration at a
# time, so they always share hardware. CPU jobs sweep ranks x threads on one 128-core node;
# the baseline is also run with 16-row tiles (`<base>_j16`), so the tiling and
# the threading effects separate; use -T "16x16" to run every revision at the current
# CPU default (16x16 columns). GPU jobs use one node with four GPUs. After the runs the
# job compares every revision's last plotfile with the baseline's using amrex fcompare.
# Summarise with
#   python3 Tests/perf/summarize.py <workdir>/runs
set -euo pipefail

ACCOUNT=""; ARCH=both; REVS=""; BASE=""; NEW=""; WORK=""; GULF=""; SIZE="512 1024 40"
NSTEP=10; QOS=debug; JOBS=16; SUBMIT=1; CPU_CFGS="128x1 64x2 32x4"; TILES=""; FPCONTRACT=off
XINPUTS=""; NSYS=0; BUILD=""
while getopts "A:a:r:b:n:w:g:s:t:q:j:Sc:T:F:x:NB:" opt; do
    case $opt in
        A) ACCOUNT=$OPTARG;; a) ARCH=$OPTARG;; r) REVS=$OPTARG;; b) BASE=$OPTARG;; n) NEW=$OPTARG;;
        w) WORK=$OPTARG;; g) GULF=$OPTARG;; s) SIZE=$OPTARG;; t) NSTEP=$OPTARG;;
        q) QOS=$OPTARG;; j) JOBS=$OPTARG;; S) SUBMIT=0;; c) CPU_CFGS=$OPTARG;; T) TILES=$OPTARG;;
        F) FPCONTRACT=$OPTARG;; x) XINPUTS=$OPTARG;; N) NSYS=1;; B) BUILD=$OPTARG;;
        *) echo "unknown option"; exit 1;;
    esac
done
case $FPCONTRACT in off|on|fast) ;; *) echo "-F takes off, on or fast"; exit 1;; esac
# -T "8x8 16x16 j8 none": extra CPU runs of every revision with these tile sizes
# (AxB = A by B columns, jN = N rows, iN = N columns, none = untiled), each checked
# bitwise against the baseline. The baseline is always also run with 16-row tiles.
tile_string () {   # tile_string <spec>  -> the remora.omp_tile_size value, spaces as _
    case $1 in
        none)  echo "1024000_1024000_1024000";;
        j[0-9]*) echo "1024000_${1#j}_1024000";;
        i[0-9]*) echo "${1#i}_1024000_1024000";;
        [0-9]*x[0-9]*) echo "${1%x*}_${1#*x}_1024000";;
        *) echo "bad tile spec $1" >&2; exit 1;;
    esac
}
[ -n "$ACCOUNT" ] || { echo "-A account is required"; exit 1; }

HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(git -C "$HERE" rev-parse --show-toplevel)
resolve_ref () {   # resolve_ref <ref> -> commit SHA
    # A commit SHA or a tag is taken as given. A branch name is fetched fresh from upstream,
    # then origin, so a cluster checkout follows what was pushed; a local branch of the same
    # name that points elsewhere is reported, since a stale one would otherwise be built in
    # silence. A name on no remote (a local-only branch, origin/xyz) is resolved locally.
    local ref=$1 sha remote local_sha
    if [[ $ref =~ ^[0-9a-f]{7,40}$ ]] || git -C "$REPO" rev-parse --verify -q "refs/tags/$ref" > /dev/null 2>&1; then
        git -C "$REPO" rev-parse --verify -q "$ref^{commit}" && return
    fi
    for remote in upstream origin; do
        git -C "$REPO" remote get-url $remote > /dev/null 2>&1 || continue
        git -C "$REPO" fetch -q $remote "$ref" 2>/dev/null || continue
        sha=$(git -C "$REPO" rev-parse --verify -q 'FETCH_HEAD^{commit}') || continue
        local_sha=$(git -C "$REPO" rev-parse --verify -q "refs/heads/$ref" 2>/dev/null)
        if [ -n "$local_sha" ] && [ "$local_sha" != "$sha" ]; then
            echo "note: $ref is taken from $remote ($(git -C "$REPO" rev-parse --short "$sha")); the local branch is at $(git -C "$REPO" rev-parse --short "$local_sha")" >&2
        fi
        echo "$sha"; return
    done
    sha=$(git -C "$REPO" rev-parse --verify -q "$ref^{commit}") && { echo "$sha"; return; }
    echo "cannot resolve revision '$ref': not a commit, tag, local branch or a branch on upstream/origin" >&2
    return 1
}
if [ -z "$REVS" ]; then
    if [ -z "$BASE" ]; then
        # The merge base with development, from a fresh fetch: a stale remote-tracking
        # branch would silently pick a months-old baseline.
        for remote in upstream origin; do
            if git -C "$REPO" remote get-url $remote > /dev/null 2>&1; then
                git -C "$REPO" fetch -q $remote development 2>/dev/null \
                    && BASE=$(git -C "$REPO" merge-base HEAD FETCH_HEAD) && break
            fi
        done
        [ -n "$BASE" ] || { echo "could not determine the baseline; pass -b <ref>"; exit 1; }
    fi
    REVS="base:$BASE new:${NEW:-HEAD}"
fi
LABELS=(); REFS=()
for pair in $REVS; do LABELS+=("${pair%%:*}"); REFS+=("${pair#*:}"); done
BASE_LABEL=${LABELS[0]}
# The runs of a configuration, as tag:label:tile triples (tile empty = the binary's default)
GPU_VARIANTS=""; CPU_VARIANTS="$BASE_LABEL:$BASE_LABEL: ${BASE_LABEL}_j16:$BASE_LABEL:1024000_16_1024000"
for l in "${LABELS[@]}"; do
    GPU_VARIANTS="$GPU_VARIANTS $l:$l:"
    [ "$l" = "$BASE_LABEL" ] || CPU_VARIANTS="$CPU_VARIANTS $l:$l:"
done
for spec in $TILES; do
    for l in "${LABELS[@]}"; do
        [ "$l" = "$BASE_LABEL" ] && [ "$spec" = j16 ] && continue   # already in the list
        CPU_VARIANTS="$CPU_VARIANTS ${l}_$spec:$l:$(tile_string $spec)"
    done
done
if [ -n "$GULF" ]; then
    GULF=$(cd "$GULF" 2>/dev/null && pwd) || { echo "-g: no such directory"; exit 1; }
    ls "$GULF"/GS_toyGrd_lev0_v1.3.nc "$GULF"/GS_toyIni_lev0_v1.3.nc > /dev/null \
        || { echo "-g: $GULF has no GS_toy*.nc files"; exit 1; }
fi
[ -n "$WORK" ] || WORK=${SCRATCH:-$HOME}/remora_perf/$(date +%Y%m%d_%H%M)
mkdir -p "$WORK/bin" "$WORK/runs" && WORK=$(cd "$WORK" && pwd)
# Source trees and built binaries live in the build directory; each work directory keeps a
# copy of the binaries it runs, so a later build cannot change a job that is still queued.
# Two invocations sharing a build directory must not build at the same time.
[ -n "$BUILD" ] || BUILD=$WORK
mkdir -p "$BUILD/bin" && BUILD=$(cd "$BUILD" && pwd)
for i in "${!LABELS[@]}"; do
    sha=$(resolve_ref "${REFS[$i]}") || exit 1
    printf '%-8s %s\n' "${LABELS[$i]}" "$(git -C "$REPO" log -1 --format='%h %ad %s' --date=short "$sha")"
done
echo "revisions: $REVS  (baseline $BASE_LABEL)  work $WORK  builds ${BUILD:-$WORK}"

# ---------------------------------------------------------------- environment and builds
if [ -f "$REPO/Build/saul-env.sh" ] && command -v module > /dev/null; then
    source "$REPO/Build/saul-env.sh"
fi
build_one () {   # build_one <label> <ref> <cpu|gpu>
    # One source tree per revision and architecture: `make clean` empties the whole
    # tmp_build_dir, so a shared tree would lose the other architecture's objects.
    local tag=$1 ref=$2 arch=$3 src=$BUILD/src_${1}_$3 flags sha
    sha=$(resolve_ref "$ref") || exit 1
    if [ ! -d "$src" ]; then
        git -C "$REPO" worktree add --detach "$src" "$sha"
    elif [ "$(git -C "$src" rev-parse HEAD)" != "$sha" ]; then
        echo "$src is at $(git -C "$src" rev-parse --short HEAD); moving it to $(git -C "$REPO" rev-parse --short "$sha")"
        git -C "$src" checkout -q --detach "$sha"
    fi
    git -C "$src" submodule update --init Submodules/AMReX
    if [ $arch = gpu ]; then
        module load cudatoolkit 2>/dev/null || true
        flags="USE_OMP=FALSE USE_CUDA=TRUE"
    else
        # Every CPU build, the baseline included, gets the same contraction setting, so the
        # bitwise comparison is not confused by a baseline that predates the default of
        # -ffp-contract=off (Docs: Performance tuning). -F fast allows FMA fusion everywhere.
        flags="USE_OMP=TRUE USE_CUDA=FALSE XTRA_CXXFLAGS=-ffp-contract=$FPCONTRACT"
    fi
    local cmd="make COMP=gnu USE_MPI=TRUE USE_PNETCDF=TRUE TINY_PROFILE=TRUE USE_PARTICLES=FALSE $flags"
    # GNU make rebuilds an object when its sources change, not when its flags or the
    # makefiles do, so a tree reused for another revision or another command line would
    # keep stale objects. Start that configuration over when either differs from the
    # record of the last build here.
    local record="$sha $cmd" recfile=$BUILD/bin/remora_${tag}_${arch}.build
    if [ ! -f "$recfile" ] || [ "$(cat "$recfile")" != "$record" ]; then
        (cd "$src/Exec" && $cmd clean > /dev/null 2>&1) || true
    fi
    (cd "$src/Exec" && $cmd -j "$JOBS" > "$BUILD/build_${tag}_${arch}.log" 2>&1) \
        || { echo "build of $tag ($arch) failed, see $BUILD/build_${tag}_${arch}.log"; exit 1; }
    echo "$record" > "$recfile"
    local exe
    exe=$(ls -t "$src"/Exec/REMORA3d.gnu.*.ex | head -1)
    cp "$exe" "$BUILD/bin/remora_${tag}_${arch}"
    [ "$BUILD" = "$WORK" ] || cp "$BUILD/bin/remora_${tag}_${arch}" "$WORK/bin/remora_${tag}_${arch}"
    echo "built $WORK/bin/remora_${tag}_${arch} from $(git -C "$src" rev-parse --short HEAD)"
}
build_fcompare () {
    if [ ! -x "$BUILD/bin/fcompare" ]; then
        local src   # any source tree built in this run has the AMReX submodule
        src=$(ls -d "$BUILD"/src_${LABELS[-1]}_* 2>/dev/null | head -1)
        [ -n "$src" ] || { echo "no source tree to build fcompare from"; exit 1; }
        (cd "$src/Submodules/AMReX/Tools/Plotfile" && make -j "$JOBS" programs=fcompare \
            > "$BUILD/build_fcompare.log" 2>&1) \
            || { echo "fcompare build failed, see $BUILD/build_fcompare.log"; exit 1; }
        cp "$src"/Submodules/AMReX/Tools/Plotfile/fcompare*.ex "$BUILD/bin/fcompare"
    fi
    [ "$BUILD" = "$WORK" ] || cp "$BUILD/bin/fcompare" "$WORK/bin/fcompare"
}
for arch in cpu gpu; do
    if [ $ARCH = both ] || [ $ARCH = $arch ]; then
        for i in "${!LABELS[@]}"; do build_one "${LABELS[$i]}" "${REFS[$i]}" $arch; done
    fi
done
build_fcompare

# ---------------------------------------------------------------- inputs
read -r NX NY NZ <<< "$SIZE"
# Upwelling's domain is 41 km wide and 150 m deep. The barotropic step dtfast = dt/30 is
# limited by the fast wave c = sqrt(g*150 m) across the cell diagonal: keep its Courant
# number dtfast*c*sqrt(2)/dx near 0.5 (ROMS blows up above about 0.8), so dtfast = 0.35 dx/c.
DT=$(python3 -c "import math; dx=41000.0/$NX; print(round(min(50.0, 30*0.35*dx/math.sqrt(9.81*150)), 2))")
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
# run_set <arch> <case> <inputs> <config tag> <srun prefix> [extra inputs args] is defined
# inside the job scripts: it runs every revision (and, on CPU, the baseline with the current
# tile default), then compares each revision's last plotfile with the baseline's.
common_body () {
cat <<EOF
FC=$WORK/bin/fcompare
BASE_LABEL=$BASE_LABEL
XINPUTS="$XINPUTS"
CPU_VARIANTS="$CPU_VARIANTS"
GPU_VARIANTS="$GPU_VARIANTS"
run_set () {
    local arch=\$1 case=\$2 inp=\$3 cfg=\$4 srun_prefix=\$5; shift 5
    local triples="\$GPU_VARIANTS"; [ \$arch = cpu ] && triples="\$CPU_VARIANTS"
    local variants=""
    for t in \$triples; do
        local v=\${t%%:*} rest=\${t#*:}; local label=\${rest%%:*} tile=\${rest#*:}
        variants="\$variants \$v"
        local bin=$WORK/bin/remora_\${label}_\$arch
        local extra=""; [ -n "\$tile" ] && extra="remora.omp_tile_size=\${tile//_/ }"
        local out=$WORK/runs/\${arch}_\${case}_\${cfg}_\$v; mkdir -p \$out
        echo "=== \$arch \$case \$cfg \$v  \$(date +%T)"
        \$srun_prefix \$bin \$inp remora.plot_file=\$out/plt \$XINPUTS \$extra "\$@" > \$out/log 2>&1 \\
            || echo "RUN FAILED: \$out"
        grep -m1 "TinyProfiler total time" \$out/log | sed "s/^/   /"
        awk '/Incl\\. Min/ {incl=1} incl && v == "" && /^REMORA::Evolve\\(\\)/ {v=\$(NF-1)}
             /TinyProfiler total time/ {t=\$NF}
             END {if (v != "") print "   Evolve " v " s"; else print "   Evolve n/a (no inclusive table); TinyProfiler total " t " s"}' \$out/log
    done
    local a=\$(ls -d $WORK/runs/\${arch}_\${case}_\${cfg}_\$BASE_LABEL/plt* 2>/dev/null | tail -1)
    for v in \$variants; do
        [ \$v = \$BASE_LABEL ] && continue
        local b=\$(ls -d $WORK/runs/\${arch}_\${case}_\${cfg}_\$v/plt* 2>/dev/null | tail -1)
        [ -n "\$a" ] && [ -n "\$b" ] || continue
        local log=$WORK/runs/\${arch}_\${case}_\${cfg}_fcompare_\$v.log
        if \$FC \$a \$b > \$log 2>&1 && grep -q AGREE \$log; then
            echo "   \$BASE_LABEL vs \$v: IDENTICAL"
        else
            echo "   \$BASE_LABEL vs \$v: DIFFER (see \$log)"
        fi
    done
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
#SBATCH -t 00:30:00
#SBATCH -J remora_perf_cpu
#SBATCH -o $WORK/runs/cpu_%j.out
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
cd $WORK/runs
$(common_body)
# ranks x threads on one 128-core node; the few-rank configurations take hours for a
# revision whose threads idle, so they are opt-in through -c
for cfg in $CPU_CFGS; do
    R=\${cfg%x*}; T=\${cfg#*x}
    export OMP_NUM_THREADS=\$T
    run_set cpu upwelling $WORK/runs/upwelling_cpu.in \$cfg "srun -n \$R -c \$((256/R)) --cpu-bind=cores"
done
EOF
if [ -n "$GULF" ]; then
cat >> "$cpu_job" <<EOF
cd $GULF
for cfg in 8x1 8x4; do
    R=\${cfg%x*}; T=\${cfg#*x}
    export OMP_NUM_THREADS=\$T
    run_set cpu gulf $WORK/runs/gulf.in \$cfg "srun -n \$R -c \$((256/R)) --cpu-bind=cores"
done
EOF
fi

# one rank per GPU; CUDA_VISIBLE_DEVICES inverse to the local task id for GPU-aware MPI
# With NSYS_OUT set, rank 0 runs under Nsight Systems: one rank keeps the report small, and the
# ranks do the same work. TinyProfiler regions appear as NVTX ranges in CUDA builds.
cat > "$WORK/gpu_wrap.sh" <<'EOF'
#!/bin/bash
export CUDA_VISIBLE_DEVICES=$((3-SLURM_LOCALID))
if [ -n "${NSYS_OUT:-}" ] && [ "$SLURM_PROCID" = 0 ]; then
    exec nsys profile --trace=cuda,nvtx,mpi,osrt --mpi-impl=mpich --force-overwrite true -o "$NSYS_OUT" "$@"
fi
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
#SBATCH -t 00:30:00
#SBATCH -J remora_perf_gpu
#SBATCH -o $WORK/runs/gpu_%j.out
export MPICH_OFI_NIC_POLICY=GPU
cd $WORK/runs
$(common_body)
SRUN="srun -n 4 --cpus-per-task=32 --cpu-bind=cores $WORK/gpu_wrap.sh"
run_set gpu upwelling $WORK/runs/upwelling_gpu.in 4gpu "\$SRUN" amrex.use_gpu_aware_mpi=1
EOF
if [ -n "$GULF" ]; then
cat >> "$gpu_job" <<EOF
cd $GULF
run_set gpu gulf $WORK/runs/gulf.in 4gpu "\$SRUN" amrex.use_gpu_aware_mpi=1
EOF
fi

# -N: one more Upwelling run of each revision under Nsight Systems, after the timed runs and
# kept out of the timings (the run directories end in _nsys, which summarize.py skips).
if [ $NSYS = 1 ]; then
cat > "$WORK/nsys_runs.sh" <<'NSYS_EOF'
#!/bin/bash
# nsys_runs.sh <work dir> <inputs> "<revision labels>" "<extra inputs>"
# Writes runs/gpu_upwelling_4gpu_<label>_nsys/profile.nsys-rep (open it in the Nsight Systems
# GUI for the timeline), CSV summaries of GPU kernels, CUDA API calls and NVTX ranges
# (stats_*.csv), and prints the ten costliest kernels of each revision.
W=$1 inp=$2 labels=$3 xin=$4
SRUN="srun -n 4 --cpus-per-task=32 --cpu-bind=cores $W/gpu_wrap.sh"
module load cudatoolkit 2>/dev/null || true
command -v nsys > /dev/null || { echo "nsys not found; load the cudatoolkit module"; exit 1; }
for v in $labels; do
    out=$W/runs/gpu_upwelling_4gpu_${v}_nsys; mkdir -p "$out"
    echo "=== nsys $v  $(date +%T)"
    NSYS_OUT=$out/profile $SRUN "$W/bin/remora_${v}_gpu" "$inp" remora.plot_int=-1 \
        remora.check_int=-1 $xin amrex.use_gpu_aware_mpi=1 > "$out/log" 2>&1 \
        || echo "   NSYS RUN FAILED: $out/log"
    [ -f "$out/profile.nsys-rep" ] || { echo "   no report written"; continue; }
    nsys stats --report cuda_gpu_kern_sum,cuda_api_sum,nvtx_sum --format csv \
        --force-export true -o "$out/stats" "$out/profile.nsys-rep" > "$out/stats.log" 2>&1
    k=$(ls "$out"/stats*cuda_gpu_kern_sum*.csv 2>/dev/null | head -1)
    [ -n "$k" ] || { echo "   no kernel summary, see $out/stats.log"; continue; }
    echo "   top GPU kernels: % of GPU kernel time, calls, name"
    python3 - "$k" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
# the column is "Time (%)" in recent nsys and "Time(%)" in older ones
pct = next((c for c in (rows[0].keys() if rows else []) if c.replace(" ", "") == "Time(%)"), None)
for r in rows[:10]:
    p = f"{float(r[pct]):6.1f}" if pct else "     ?"
    print(f"   {p} %  {int(r['Instances']):8d}  {r['Name'][:90]}")
total = sum(int(r['Instances']) for r in rows)
print(f"   {len(rows)} distinct kernels, {total} launches")
PY
done
NSYS_EOF
chmod +x "$WORK/nsys_runs.sh"
labels=""; for t in $GPU_VARIANTS; do labels="$labels ${t%%:*}"; done
printf 'bash %s/nsys_runs.sh %s %s/runs/upwelling_gpu.in "%s" "%s"\n' \
    "$WORK" "$WORK" "$WORK" "$labels" "$XINPUTS" >> "$gpu_job"
fi

echo "job scripts: $cpu_job $gpu_job"
if [ $SUBMIT = 1 ]; then
    { [ $ARCH = both ] || [ $ARCH = cpu ]; } && sbatch "$cpu_job"
    { [ $ARCH = both ] || [ $ARCH = gpu ]; } && sbatch "$gpu_job"
    echo "when the jobs finish: python3 $HERE/summarize.py $WORK/runs"
fi
