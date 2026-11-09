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
  `Array4`s taken from `MFIter`. Write `Real(0.25)` for literals.
- Use `amrex::` types, existing helpers and AMReX routines (FillPatch, interpolaters,
  `average_down`, ParallelFor reductions) instead of writing new versions. Before adding a
  routine, search for one that already does it, including in AMReX. Pull duplicated code into
  one helper rather than copying it again.
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
  code needs explicit casts.

## Runtime Options

- New options go in the `remora.` namespace. AMReX owns `amr.`. When moving an option out of
  `amr.`, keep the `amr.` name as a fallback.
- Read options with `pp.query`/`pp.queryAdd` beside the related ones in `REMORA_DataStruct.H`.
- A choice of behaviour becomes a runtime option whose default matches ROMS (or the physical
  constant), so runs can still be compared with ROMS.
- New initializers and forcing follow the existing pattern: one constant/analytic/netcdf
  selector, with variables set by name from the inputs file.
- Don't add options, switches or hard stops nobody asked for.
- Document every new option in `Docs/sphinx_doc/` in the same PR.

## Design

- The goal is to match ROMS. Before designing a numerical change, read how ROMS does it, and
  ERF for AMR and subcycling. When porting, copy the reference faithfully, including its
  branching structure and quirks; fix only real bugs such as race conditions, and list
  suspected bugs apart from intentional differences.
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

- REMORA is tested mainly by running the model. A test is a small input deck in
  `Tests/test_files/<Name>/` (a few steps on a small grid, seconds on 2 MPI ranks) plus an
  assertion from `Tests/CTestList.cmake`. Most model code works on `MultiFab`s owned by a
  `REMORA` object built from an inputs file, so a short run tests it where it is used, at
  little cost. `testing.rst` describes each assertion.
- Choose the assertion that says what must be true, not just that nothing changed. A gold
  file (`add_test_r`) shows only that a run still matches the output its gold was made from:
  it cannot show that output was right, or that a feature still takes effect. Where they fit,
  prefer `add_test_extrema` against a value known from outside REMORA, `add_test_equiv` or
  `add_test_r_selfcompare` for two routes that must agree, `add_test_r_differ` for a feature
  meant to change the answer, `add_test_conservation` with an `above` control, and
  `add_test_abort` for inputs that must be rejected. Add a new gold file only when none of
  these fits.
- Write a unit test (`Tests/Unit`) only for a self-contained function that can be called
  without a `REMORA` object or a run, such as header-only arithmetic, calendar or
  interpolation code. Take the expected values from outside the implementation (ROMS, a
  closed form, an identity). A unit test is a plain `main()` that returns nonzero on failure,
  with no framework. Don't restructure model code just to make it unit-testable.
- Make variants of an existing case instead of copying it: use `compose_test_input`, the
  `INPUT_NAME` argument, or command-line options. `remora.max_step = 0` with
  `remora.plot_int = 1` sets up the run and writes a plotfile in about a second, which is
  enough to check an initial condition.
- Check that a new test fails when the feature breaks: break the change, rebuild, and confirm
  the test fails.
- Register tests in `Tests/CTestList.cmake` with the existing helpers (`add_test_r`,
  `add_test_equiv`, `add_test_abort`, ...). Run `ctest -R <Name>` first, then the full suite.
- Before pushing or opening a PR, a local Release or test build and ctest is enough; don't
  run a local Debug build, which takes too long. CI's Debug build traps FPEs on signalling
  NaNs that Release misses, so check its result after pushing. To reproduce a Debug trap in a
  Release run, add `amrex.init_snan=1 amrex.fpe_trap_invalid=1`.
- CI has no NetCDF; NetCDF cases run only in the nightly suite. Prefer tests that need no
  NetCDF. Also build the no-PnetCDF configuration when the work was done in a PnetCDF build.
- Before reporting a failure or bug, reproduce it on unchanged `development`. CI passes there.
- For a refactor, show the output is unchanged: build a reference binary from the base commit
  and compare plotfiles with `fcompare` (success prints `PLOTFILE AGREE`). If no test covers the
  code, add a regression test, or a unit test if the code is self-contained.
- Watch for degenerate checks. Passive tracers are identically zero in several cases, open
  boundary conditions can overwrite edge faces, and some paths have no test at all: the CFL
  timestep (every input sets `fixed_dt`), NetCDF history output, NetCDF forcing and
  climatology, and (unless CI builds them) particles. Confirm a field is non-trivial with
  `fextrema`, and say when a change was verified by hand.
- State what each number measures and in what units, and check record or time alignment
  before trusting a comparison.
- Device code must run on GPU. After pushing, check the CUDA/HIP/SYCL CI jobs.

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
