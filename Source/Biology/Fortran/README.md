# Fennel Fortran Bridge — Build and Run Guide

Operator guide for the Fennel biology parity campaign. The diagnostic record
schema and tag table are frozen separately in [`tag_map.md`](tag_map.md); this
file covers how to build the thing and drive it.

## What this is

Two implementations of the same kernel, in one executable, selected at run time:

| | Path | Selected by |
|---|---|---|
| **Path A** | ROMS `biology_tile`, the oracle | `remora.use_biology_cpp_answer = 0` |
| **Path B** | Native `REMORA::advance_biology` | `remora.use_biology_cpp_answer = 1` (default) |

Path A runs the tracked copy of the ROMS Fennel kernel
(`REMORA_fennel.h`, instrumented only with tag calls at block seams) through an
`iso_c_binding` wrapper. Comparing the two paths never requires a rebuild.

Files:

| File | Role |
|---|---|
| `REMORA_fennel.h` | Tracked copy of ROMS `fennel.h` + tag calls. Diff against `../fennel.h` should show additions only. |
| `REMORA_fennel_mod.h` | Tracked copy of ROMS `fennel_mod.h`, verbatim |
| `REMORA_fennel_roms.F` | ROMS module stubs, the two tracked copies, and the isohelper — one translation unit |
| `REMORA_FennelBridge.cpp` | Pack/unpack between REMORA state and ROMS-shaped buffers |
| `REMORA_Fennel_Fortran_Interface.H` | `extern "C"` declarations |
| `remora_bio_cppdefs.h` | **The compiled ROMS option set.** Edit this to validate a different configuration. |
| `set_bounds.h`, `tile.h` | Minimal replacements for the ROMS includes of the same name |

## Build

Two independent flags, both off by default:

| GNUmake | CMake | Effect |
|---|---|---|
| `USE_FENNEL_FORT=TRUE` | `-DREMORA_ENABLE_FENNEL_FORT=ON` | Builds the Path A oracle. Defines `REMORA_USE_FENNEL_FORT`, and on the GNUmake side flips `BL_NO_FORT=FALSE` so AMReX compiles Fortran at all. |
| `USE_BIOLOGY_DIAG=TRUE` | `-DREMORA_ENABLE_BIOLOGY_DIAG=ON` | Compiles in the tag emitters on both paths. |

They are deliberately separate. The native kernel carries tag instrumentation
whether or not the oracle is built, and the C++ emitters sit inside the
`AMREX_GPU_DEVICE` column lambda — leaving them in costs device `printf`
machinery and closure captures in every production build, so they are gone
unless asked for. With the flag on, `remora.biology_debug` still selects
verbosity at run time; no rebuild is needed to change paths or levels during a
campaign.

Asking for something that was not compiled in aborts with a message naming the
flag to rebuild with. Neither option silently degrades: `use_biology_cpp_answer
= 0` without the bridge does not fall back to Path B, and `biology_debug > 0`
without the diagnostics does not produce an empty log.

```sh
# GNUmake, full validation build
cd Exec
make -j8 USE_FENNEL_FORT=TRUE USE_BIOLOGY_DIAG=TRUE
```

```sh
# CMake equivalent
cmake -S . -B build -DREMORA_ENABLE_FENNEL_FORT=ON -DREMORA_ENABLE_BIOLOGY_DIAG=ON
cmake --build build -j8
```

`REMORA_ENABLE_FENNEL_FORT` must be on the **initial** configure command line:
`project()` consults it to decide whether to enable the Fortran language,
before the `option()` declaration is reached. Turning it on in an existing
build directory that was configured without it produces a `FATAL_ERROR` telling
you to reconfigure; a fresh build directory always works.

The GNUmake executable name reflects `Exec/GNUmakefile` settings, e.g.
`REMORA3d.gnu.TEST.MPI.PNC.ex` for `COMP=gnu TEST=TRUE USE_MPI=TRUE
USE_PNETCDF=TRUE`.

### Parity build

For first-parity work, disable optimization and FMA contraction so neither
side reassociates:

```sh
cd Exec
make -B -j8 USE_FENNEL_FORT=TRUE USE_BIOLOGY_DIAG=TRUE \
     GENERIC_GNU_FLAGS="-O0 -ffp-contract=off"
```

`GENERIC_GNU_FLAGS` reaches `CXXFLAGS`, `FFLAGS` and `F90FLAGS`, so both paths
get it.

In practice the residuals are the same at O0 and O3, so the optimized build is
fine for routine work; the O0 build exists to remove compiler reassociation as
an explanation when a divergence is under investigation.

### `make -B` when toggling flags

GNUmake keys off timestamps and does not notice a flag-only change. Two
distinct failure modes:

- Changing `GENERIC_GNU_FLAGS` silently mixes O0 and O3 objects.
- Toggling `USE_FENNEL_FORT` changes `BL_NO_FORT`, which changes how AMReX's
  **C++** is compiled. Going from a no-Fortran build to a bridge build without
  `-B` fails at link with `undefined reference to bl_proffortfuncstart_cpp_`
  and friends: the stale `AMReX_BLProfiler.o` was compiled with `-DBL_NO_FORT`
  and omits the symbols the Fortran objects need.

Use `-B` whenever the flag set changes. It is unnecessary when you have only
edited sources.

## Run

Set up a scratch directory — never run validation in `Exec/BioToy` itself, or
you will overwrite prior evidence:

```sh
mkdir -p /tmp/fennel_sc && cd /tmp/fennel_sc
B=$REMORA_HOME/Exec/BioToy
cp $B/inputs .
ln -sf $B/bio_toy_ini_fennel_classic64.nc .
ln -sf $B/bio_toy_grd_classic64.nc .
ln -sf $B/bio_toy_frc_classic64.nc .
```

Note the three NetCDF inputs are excluded by `.gitignore:18` (`Exec/*/*.nc`)
and no `.nc` file is tracked anywhere in the repo, so a fresh clone cannot run
BioToy without obtaining them separately.

Then run each path into its own log:

```sh
EXE=$REMORA_HOME/Exec/REMORA3d.gnu.TEST.MPI.PNC.ex

# Path A (oracle)
$EXE inputs remora.use_biology_cpp_answer=0 remora.biology_debug=1 \
     remora.max_step=1 remora.plot_int=-1 > log.fortran

# Path B (native)
$EXE inputs remora.use_biology_cpp_answer=1 remora.biology_debug=1 \
     remora.max_step=1 remora.plot_int=-1 > log.cpp
```

Serial is the default acceptance scope. MPI and GPU parity are separate lanes
and are not opened by this campaign; the bridge aborts under `AMREX_USE_GPU`
with a message saying so.

### Runtime options

| Option | Meaning |
|---|---|
| `remora.use_biology_cpp_answer` | `0` = Path A oracle, `1` = Path B native (default). Needs the bridge compiled in. |
| `remora.biology_debug` | `0` off, `1` target column only, `2` every column. Any value > 0 needs the diagnostics compiled in. |
| `remora.biology_debug_i`, `_j` | Target column for level 1, REMORA 0-based, default `2 2` |

For BioToy (4x4x30, one biology call per step) level 1 emits 3630 records per
call and level 2 emits 58080.

## Compare

At `biology_debug=1` a plain diff is valid:

```sh
grep '^FENNEL-' log.fortran | sed 's/^FENNEL-FORT /FENNEL/' > a.txt
grep '^FENNEL-' log.cpp     | sed 's/^FENNEL-CPP  /FENNEL/' > b.txt
diff a.txt b.txt | head -40
```

Both prefixes are 12 characters, so the two `sed` patterns differ in their
space count and are not interchangeable. `diff` should report exactly two
lines per mismatching record.

**At `biology_debug=2` a positional diff is meaningless.** The two paths emit
the same records in a different order, because Path A tags at ROMS block seams
(outside the `i` loop, so all columns per tag) while Path B tags inside the
per-column `ParallelFor` (all tags per column). Match on the
`(tag, iter, i, j, k, var)` key, or sort both logs first.

Because the physics is sequential, a divergence at an early tag propagates
through every later one. Stop at the first divergent tag; do not read past it,
and do not patch logic without line-level evidence at that site.

## Validating a different option set

ROMS resolves biology options at compile time; REMORA resolves them at run
time. The bridge can therefore only speak one configuration, pinned in
`remora_bio_cppdefs.h` — currently matching `Exec/BioToy/inputs`:

```
CARBON, DENITRIFICATION, BIO_SEDIMENT, BULK_FLUXES, MASKING
(PO4, OXYGEN, ODU off)
```

To validate something else, edit that file and rebuild. Only this directory
recompiles; Path B is untouched.

`fennel_check_config` compares the compiled CPP set against REMORA's runtime
flags on every call and stops if they disagree, so a mismatch can never
produce a plausible-looking but meaningless comparison:

```
 FENNEL BRIDGE: runtime option set does not match the compiled ROMS CPP set.
 Comparing these two paths would be invalid.

   option           runtime   compiled
   carbon                 0          1
   tracer count           7         11

 Edit Source/Biology/Fortran/remora_bio_cppdefs.h to match the case, then rebuild.
```

It exits non-zero. Verify with
`... remora.use_biology_cpp_answer=0 remora.fennel.carbon=0`, which should
abort with `EXIT=1`.

The guard checks the **option set only**. It cannot detect that the cppdefs
select a numerical variant the native port does not implement — e.g.
`RW14_CO2_SC` against the Wanninkhof 1992 coefficients in the C++. Those are
not threaded through the interface.

## Traps

**The preprocessor eats padding inside Fortran character literals.** AMReX
compiles `.F` through `cpp -E -traditional-cpp -P`, which collapses runs of two
or more spaces even inside quotes. `'G00_PRE   '` reaches the compiler as
`'G00_PRE '`, and `Aw` right-justifies anything shorter than `w`. This silently
destroyed the tag column alignment and made every log line appear to differ.
Both emitters now assign into fixed-length locals, and format literals use `nX`
descriptors rather than space runs. Do not "simplify" either back.

**`STOP '<string>'` exits zero in gfortran.** A character stop code is written
to stderr with exit status 0; only an integer code or `ERROR STOP` gives
non-zero. Any new abort path here must use `ERROR STOP`, or an exit-status gate
will score a failed run as a pass.

**Do not modify `../fennel.h` or `../fennel_mod.h`.** They are the pristine
upstream baseline. Verify the tracked copy is additions-only:

```sh
diff ../fennel.h REMORA_fennel.h | grep -c '^<'   # expect 0
```

## Scope

What the current evidence covers: serial, one tile, BioToy's option set, and
tag-level rather than plotfile comparison. BioToy's 4x4 grid is doubly periodic
with uniform IC and forcing, so all 16 columns are **degenerate replicas** —
`biology_debug=2` widens the record count but adds no independent evidence.
Genuine multi-column coverage needs a case with horizontal structure.

See `FENNEL_BRIDGE_HANDOFF.md` at the repo root for current parity results and
the open lane list.
