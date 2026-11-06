# Performance comparison on Perlmutter

`perlmutter_phase1.sh` builds two revisions of REMORA with GNU make (MPI + OpenMP for CPU
nodes, MPI + CUDA for GPU nodes, PnetCDF and TinyProfiler on), generates the timing inputs
and submits one CPU job and one GPU job. Each job runs the baseline and the new binary back
to back on the same node for every configuration and compares their last plotfiles with
`fcompare`, so the pair always shares hardware and the answers are checked bitwise.

```
source Build/saul-env.sh                                    # modules
Tests/perf/perlmutter_phase1.sh -A m1234 -g $SCRATCH/REMORA-data/GulfRefinementTest
...
python3 Tests/perf/summarize.py $SCRATCH/remora_perf/<date>/runs
```

Cases:

- `upwelling`: `inputs/upwelling_perf.in`, the Upwelling problem at `-s "NX NY NZ"`
  (default 512 x 1024 x 40, 21 M cells, 10 steps) with 128x128 boxes on CPU and 256x256 on GPU and a
  time step that keeps the barotropic Courant number near 0.5. CPU jobs sweep 128x1, 64x2 and 32x4 ranks x threads on one
  node by default (`-c "128x1 64x2 32x4 16x8"` to change it); the GPU job uses four GPUs.
- `gulf` (with `-g`): `inputs/gulf_perf.in`, the nested Gulf Stream toy case of
  `Exec/GulfRefinementTest` with level 1 active, two-way coupling and NetCDF grid and
  initial data; it needs the `GS_toy*.nc` files from the REMORA-data repository.

Neither case has land. For the land-layout options (`remora.drop_dry_boxes`,
`remora.load_balance`), `make_masked_grid.py <REMORA-data>/IdealMiniGrid <dir> <factor>`
writes a staircase-coastline grid and initial file (36 % land, reaching the domain edges)
that `inputs/masked_grid_perf.in` runs from `<dir>`; it is not yet wired into the script.
Run it with and without the options and compare the outputs with `remora_subset_compare`
(built with the tests), since the BoxArrays differ.

The CPU job also runs the baseline with `remora.omp_tile_size=1024000 16 1024000`
(`<baseline>_j16`, the 16-row strips that were an earlier default; the current CPU default is
16x16 columns), so the effect of the tile shape and of the threading can be told
apart, and `-T "8x8 16x16 j8 none"` adds runs of every revision at those tile sizes
(`AxB` columns, `jN` rows, `iN` columns, `none` untiled) for a tile-size scan. The summary
prints `IDENTICAL` or the largest relative difference `fcompare` reported.
`-x "k=v ..."` adds inputs to every run; `-x "$(cat Tests/perf/inputs/gls.in)"` turns on GLS
vertical mixing (the settings of `Exec/Upwelling/inputs_gls`), the costliest kernel, and
use a separate `-w` directory so those runs do not mix with the default ones.
`-B dir` keeps the source trees and built binaries in a build directory shared between
work directories. A revision is rebuilt only when it resolves to a different commit or the
make command changes (a branch name follows its remote tip; a commit SHA never moves), so
one build directory with a fresh `-w` per comparison reuses old builds and keeps each
summary table to the runs of that comparison. Each work directory gets its own copy of the
binaries it runs, so rebuilding for a later comparison cannot change a queued job; two
invocations sharing a build directory must not build at the same time. Without `-B`, the
work directory holds the builds too.
`-N` adds one Nsight Systems run of each revision to the GPU job, after the timed runs, with
rank 0 under `nsys profile`, tracing CUDA calls and NVTX ranges (TinyProfiler regions show
up as NVTX ranges in CUDA builds). `NSYS_TRACE=cuda,nvtx,mpi,osrt` in the environment also
traces MPI (as Cray MPICH) and OS-runtime calls, at a large cost in run time; each profiled
run is stopped after `NSYS_TIMEOUT` seconds (default 600) so one slow profile cannot use up
the job. Each writes `runs/gpu_upwelling_4gpu_<label>_nsys/profile.nsys-rep` for the
Nsight Systems timeline and `stats_*.csv` summaries of GPU kernels, CUDA API calls and NVTX
ranges, and the job output lists the ten costliest kernels and the launch count. These runs
are slowed by the profiler and are left out of the `summarize.py` table.
Builds use `TEST=FALSE USE_ASSERTION=FALSE` whatever `Exec/GNUmakefile` says: its `TEST = TRUE` turns
on AMReX's signalling-NaN fill of every new array (`amrex.init_snan`), which a profile showed as
the costliest GPU kernel, so timings made before this change include that fill.
Every CPU build, the baseline included, is compiled with `-ffp-contract=off` (`-F fast` to
allow FMA fusion in all of them instead), so a baseline that predates REMORA's default does
not differ at roundoff for that reason alone. A source tree that is reused for another
revision or another build command is cleaned first, since GNU make does not rebuild
objects when only flags or makefiles change.
Any number of revisions can be compared: `-r "base:c8bc0b8 p1:perf-phase1-omp p2:perf-phase2"`
builds each (the first is the baseline the others are compared with) and
runs them back to back in every configuration; `-b`/`-n` are the two-revision shorthand,
defaulting to the merge base with `upstream/development` and `HEAD`. `-S` builds without
submitting; the job scripts are left in the work directory.
