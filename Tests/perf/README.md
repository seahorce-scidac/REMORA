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
(`<baseline>_j16`), so the effect of the tile default and of the threading can be told
apart, and `-T "8x8 16x16 j8 none"` adds runs of every revision at those tile sizes
(`AxB` columns, `jN` rows, `iN` columns, `none` untiled) for a tile-size scan. The summary
prints `IDENTICAL` or the largest relative difference `fcompare` reported.
Any number of revisions can be compared: `-r "base:c8bc0b8 p1:perf-phase1-omp p2:perf-phase2"`
builds each (the first is the baseline the others are compared with) and
runs them back to back in every configuration; `-b`/`-n` are the two-revision shorthand,
defaulting to the merge base with `upstream/development` and `HEAD`. `-S` builds without
submitting; the job scripts are left in the work directory.
