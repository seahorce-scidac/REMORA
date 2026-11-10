# Repository Guidelines

REMORA is an AMReX-based ocean model (C++20 in CMake builds) that ports ROMS. Read the code
around a change and match it. The rules below are conventions that are easy to miss.

## Project Structure

- `Source/`: model code, by subsystem (`IO/`, `Initialization/`, `TimeIntegration/`,
  `BoundaryConditions/`, `Particles/`, `Biology/`, `Utils/`). The `REMORA` class is in
  `Source/REMORA.H`; runtime options are read in `Source/REMORA_DataStruct.H`.
- `Exec/<Problem>/`: problem setups and runnable input decks. `Exec/Generic/inputs_generic` is
  generated (see Git and PRs).
- `Tests/`: `CTestList.cmake` registers the tests, `test_files/<Name>/` holds their inputs,
  `REMORA_Gold_Files/<Name>/` their gold plotfiles, and `Unit/` the unit tests.
- `Docs/sphinx_doc/`: user documentation (Sphinx). Don't hand-edit `_build/` or other
  generated output.
- `Submodules/AMReX/`: vendored AMReX. A fresh worktree needs
  `git submodule update --init --recursive`.

## Build and Test

- `cmake -S . -B Build -DCMAKE_BUILD_TYPE=Release -DREMORA_ENABLE_MPI=ON
  -DREMORA_ENABLE_ALL_WARNINGS=ON -DREMORA_ENABLE_TESTS=ON -DREMORA_ENABLE_FCOMPARE=ON`:
  a test build. Without `REMORA_ENABLE_FCOMPARE` every comparison test fails as "not found".
  CI (`.github/workflows/ci.yml`) uses the same flags with `Debug`.
- `cmake --build Build -j 4`, then `ctest --test-dir Build -j 2 --output-on-failure`. Each
  regression test runs on 2 MPI ranks.
- `-DREMORA_ENABLE_PNETCDF=ON` is the only way to define `REMORA_USE_NETCDF`. There is no
  `REMORA_ENABLE_NETCDF`.
- `make -C Exec -f GNUmakefile`: the GNU make build. Pass options on the command line
  (`USE_MPI=TRUE PRECISION=FLOAT ...`) instead of editing `Exec/GNUmakefile`.
- `make -C Docs/sphinx_doc html`: build the docs.

## Coding Style

CI does not run the repository `.clang-format`, and much of the code predates it, so don't
reformat existing code. CI checks only tabs, trailing whitespace and codespell. Write new code
like its neighbours:

- Four-space indentation, no tabs. Lines run to about 100 columns.
- AMReX spacing: a space before the parameter list in declarations and definitions
  (`void advance_3d (int lev, ...)`), none in calls (`stretch_transform(lev)`). Definitions put
  the return type on its own line, and open the function brace on the next line.
- Doxygen block above each function, `@param[in   ]` / `@param[inout]` style, saying what it
  does and what corresponds in ROMS (`/** Corresponds to step3d_uv.F ... */`).
- Kernels are `ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) { ... });` on
  `Array4`s, inside an `MFIter` loop with `TilingIfNotGPU()` under
  `#pragma omp parallel if (Gpu::notInLaunchRegion())` (guarded by `#ifdef AMREX_USE_OMP`).
  Write `Real(0.25)` for literals.
- Initializers fill ghost cells too (`setVal` with `nGrowVect()`, not 0 ghosts).
- Look up state components by name (`component_index`), not by assumed position or a private
  enum.
- Use `amrex::` types, existing helpers and AMReX routines (FillPatch, interpolaters,
  `average_down`, ParallelFor reductions) instead of writing new versions. Before adding a
  routine, search for one that already does it, including in AMReX. Pull duplicated code into
  one helper rather than copying it again.
- Read time-dependent NetCDF input (forcing, boundaries, climatology, rivers) through
  `NCTimeSeries` (`Source/IO/REMORA_NCTimeSeries*.H`), which also combines records split
  across files; don't open the files yourself to decide how to read them. Read a static field
  with `BuildFABsFromNetCDFFile` and copy it into the MultiFab, as
  `REMORA::init_bathymetry_from_netcdf` does.
- Fill ghost cells of a per-level field (a vector over levels on the level layouts, like
  `vec_h`) with `REMORA::FillPatch` and a BC set such as `foextrap_periodic_bc()`. It covers
  periodic, physical-edge and coarse-to-fine ghosts; `FillBoundary` alone fills only periodic
  and internal ones. Don't write your own extrapolation loops, and don't use
  `PhysBCFunctNoOp` on a level that touches a physical boundary.
- Files are `REMORA_<Name>.cpp/.H`. Reuse existing names for existing quantities (e.g.
  `z_velocity`, as in plotfiles, not a new `w`), and don't use names that clash with ROMS
  notation (`xi`/`eta` are ROMS grid directions).
- Informational `amrex::Print()` output goes behind a verbosity check. A startup message must
  make clear whether it reports a limitation or a feature.

GPU and build-configuration traps that CPU builds miss:

- A `REMORA` member that launches kernels must be declared `public` in `REMORA.H`. nvcc rejects
  extended device lambdas in private or protected members.
- Don't pass `REMORA_Constants.H` constants (`zero`, `one`) to `amrex::max`/`min`/`Clamp` in
  device code: nvcc reports them undefined. Use `Real(0.0)`.
- Code in `Source/Initialization/REMORA_init_from_netcdf.cpp` exists only in PnetCDF builds,
  because the whole file is `#ifdef REMORA_USE_NETCDF`. Put generic helpers elsewhere, and
  compile touched files both with and without PnetCDF.
- `PRECISION=FLOAT` makes `Real` float, but `ParticleReal` stays double, so particle-to-mesh
  code needs explicit casts. `std::pow(Real, int)` also returns double, so `amrex::min`/`max`
  can't deduce its type; write `x*x` for squares.

## Runtime Options

- New options go in the `remora.` namespace. AMReX owns `amr.`. When moving an option out of
  `amr.`, keep the `amr.` name as a fallback.
- Read options with `pp.query`/`pp.queryAdd` beside the related ones in `REMORA_DataStruct.H`,
  or in the component's own parameter reader.
- A choice of behaviour becomes a runtime option whose default matches ROMS (or the physical
  constant), so runs can still be compared with ROMS.
- New initializers and forcing follow the existing pattern: one constant/analytic/netcdf
  selector, with variables set by name from the inputs file.
- Don't add options, switches or hard stops nobody asked for. But a component that doesn't
  support a shared option (such as a Fennel-only biology switch) rejects it rather than
  silently ignoring it.

## Design

- The goal is to match ROMS. Before designing a numerical change, read how ROMS does it, and
  ERF for AMR and subcycling.
- When porting, copy the physics faithfully, including its branching structure and quirks;
  fix only real bugs such as race conditions, and list suspected bugs apart from intentional
  differences. Don't transliterate the Fortran itself: use zero-based layer indices with no
  `k-1` offsets, declare `const` locals at first use, and drop dead variables.
- REMORA does not always use ROMS units. `srflx` is W/m² (ROMS: kinematic, degC m/s), while
  the other surface fluxes are kinematic. When porting a formula, check the units of every
  REMORA field it reads, and drop ROMS's `rho0*Cp` conversions where REMORA stores W/m².
- Structure a new instance of an existing kind of component (biology model, forcing type)
  like its sibling, e.g. Fennel for biology: scratch arrays, parameter reader, tests and docs.
- Keep checkpoint metadata generic: no commit hashes or option-set strings in identity
  fields, so a version bump doesn't reject older checkpoints.
- Prefer the simpler design. Ask the user before making a design decision; otherwise go as
  far as you can without intervention.
- Stay on the task. Don't chase edge cases or side issues; note them and open issues instead.
- Size performance work for real user problems, which are much larger and longer than the
  regression tests.
- Experimental features get a warning in both the docs and stdout.

## Comments and Documentation

- Write tight and clear the first time; when concision and clarity conflict, choose clarity.
  This applies to code, comments, test inputs, CTestList comments, docs, commits and PRs.
- Say what the code does and the one non-obvious reason. Don't record how the author found it
  out, development phases or test tiers, issue numbers, or how results shifted under other
  options. Describe the current state, not the history.
- Match the comment density of the surrounding code. Comments about a function's internals go
  inside the function, not in the block above it.
- No jargon or private shorthand ("lane", "index-safe", "hot"). Use plain, defined terms and
  full sentences.
- User documentation is Sphinx (`Docs/sphinx_doc/*.rst`), not Doxygen. When a change adds a
  test, an option or a behaviour change, update the docs (including `testing.rst` /
  `RegressionTests.rst` for tests) in the same PR. Extend an existing table or section rather
  than adding a parallel one.

## Testing and Verification

- Register tests in `Tests/CTestList.cmake` with the existing helpers (`add_test_r`,
  `add_test_equiv`, `add_test_abort`, ...). Run `ctest -R <Name>` first, then the full suite.
- A new model or process needs a gold-file test of a full REMORA run. Restart, alias and
  layout comparisons only show a run agrees with itself, and fixtures fed straight from ROMS
  skip REMORA's input conventions, so neither catches unit mismatches.
- Don't commit standalone reference fixtures, generator scripts or planning documents. Tests
  go in the existing structure (`CTestList.cmake`, gold files, `Tests/Unit/`).
- A kernel with per-tile state needs a tiled-versus-untiled bitwise test (`add_test_layout`)
  on a case that varies in both x and y, such as Seamount.
- Before pushing or opening a PR, a local Release build and ctest is enough; a local Debug
  build takes too long. After pushing, check CI: its Debug build traps FPEs on signalling
  NaNs that Release misses (reproduce one in Release with
  `amrex.init_snan=1 amrex.fpe_trap_invalid=1`), and its CUDA/HIP/SYCL jobs are the first
  run of device code on GPU.
- CI has no NetCDF; NetCDF cases run only in the nightly suite. Prefer tests that need no
  NetCDF. Also build the no-PnetCDF configuration when the work was done in a PnetCDF build.
- Before reporting a failure or bug, reproduce it on unchanged `development`. CI passes there.
- For a refactor, show the output is unchanged: build a reference binary from the base commit
  and compare plotfiles with `fcompare` (success prints `PLOTFILE AGREE`). Add a unit test if
  none covers the code.
- Watch for degenerate checks. Passive tracers are identically zero in several cases, open
  boundary conditions can overwrite edge faces, and some paths have no test at all: the CFL
  timestep (every input sets `fixed_dt`), NetCDF history output, NetCDF forcing and
  climatology, and (unless CI builds them) particles. Confirm a field is non-trivial with
  `fextrema`, and say when a change was verified by hand.
- State what each number measures and in what units, and check record or time alignment
  before trusting a comparison.

## Git and Pull Requests

- Branch from and target `development`. One topic per PR; put unrelated fixes in their own
  PRs. Use one commit per logical change, so each justification is easy to follow.
- Check GitHub before starting: the work may already be merged, open, or covered by an issue.
- Don't push, comment on GitHub, or open a PR without the user's approval. Force-push only
  after tests pass.
- Never commit `Exec/Generic/inputs_generic` on a feature branch; a bot regenerates it after
  each merge. Preview a branch's effect with
  `.github/workflows/style/check_generic_inputs.sh --preview <remote>/development`.
- Commit subjects and PR titles are short and imperative ("Read only each rank's part of a
  lat-lon forcing record"); PRs are squash-merged, so the title becomes the commit summary.
- PR body: `## Summary` (the problem, then what the PR changes, as bullets) and `## Testing`
  (build configuration, cases, results with numbers, gold-file changes, and what was not run).
  Close the issues the PR resolves, and open one issue per deferred finding.

## Reporting

- Cite file and line for claims about the code, and give sources for names and formulas.
- Explain any shorthand in replies.
