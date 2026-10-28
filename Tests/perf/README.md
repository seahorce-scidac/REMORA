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
  (default 1024 x 2048 x 40, 84 M cells) with 128x128 boxes on CPU and 256x256 on GPU and a
  time step that keeps the barotropic Courant number near 0.5. CPU jobs sweep 128x1, 64x2, 32x4, 16x8 and 8x16
  ranks x threads on one node; the GPU job uses four GPUs.
- `gulf` (with `-g`): `inputs/gulf_perf.in`, the nested Gulf Stream toy case of
  `Exec/GulfRefinementTest` with level 1 active, two-way coupling and NetCDF grid and
  initial data; it needs the `GS_toy*.nc` files from the REMORA-data repository.

The CPU job also runs the baseline with `remora.omp_tile_size=1024000 16 1024000`
(`base_j16`), so the effect of the tile default and of the threading can be told apart.
Baseline defaults to the merge base with `upstream/development`; `-b` and `-n` choose other
refs. `-S` builds without submitting; the job scripts are left in the work directory.
