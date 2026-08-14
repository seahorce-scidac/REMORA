# Fennel Parity Tag Map (frozen)

Tier-1 canonical diagnostic surface for the REMORA Fennel bridge-vs-native
parity campaign. This file is the contract: **tag names, tag order, variable
names, and field order are frozen**. Any divergence in tag naming between
Path A (Fortran bridge) and Path B (native C++) is forbidden. Changing the
schema invalidates prior ledger rows and requires a new campaign generation.

## Vocabulary

- `path_a_fortran_bridge` — `Source/Biology/Fortran/REMORA_fennel.h`, the
  tracked instrumented copy of ROMS `fennel.h`. Oracle.
- `path_b_native_cpp` — `Source/Biology/REMORA_Biology.cpp::advance_biology`.
- Selected at run time by `remora.use_biology_cpp_answer` (no rebuild).

## Record schema

One record per (tag, iteration, column, level, tracer). Emitted by
`fennel_tag` (Fortran) and `remora_fennel_tag` (C++). The two must produce
byte-identical output apart from the leading path token — verified by the
plain `diff` in the sidecar protocol below reporting exactly two lines per
mismatching record. See the preprocessing trap below for what breaks this.

```
FENNEL-<PATH> <TAG> iter=%3d i=%5d j=%5d k=%4d <VAR> <VALUE>
```

| Field   | Width / format        | Notes                                        |
|---------|-----------------------|----------------------------------------------|
| `PATH`  | `FORT` or `CPP`       | Only differing token on an otherwise equal line |
| `TAG`   | 10 chars, left-just.  | From the table below                          |
| `iter`  | `%3d`                 | ROMS `Iter`, 1-based, both paths              |
| `i`,`j` | `%5d`                 | REMORA 0-based index on both paths            |
| `k`     | `%4d`                 | REMORA 0-based index on both paths (`k=0` = bottom) |
| `VAR`   | 4 chars, left-just.   | From the variable table below                 |
| `VALUE` | `ES26.17E2` / `%26.17E` | 18 significant digits, 2-digit exponent     |

Index convention is **REMORA's on both paths**. The Fortran side subtracts one
from its `k` before printing so the two logs diff without an index map. This
is the single normalization applied to Path A output; values are untouched.

### Trap: cpp eats padding inside Fortran character literals

AMReX compiles `.F` sources through `cpp -E -traditional-cpp -P`, which
collapses runs of two or more spaces **even inside Fortran character
literals**. The 10-character tag literals in `REMORA_fennel.h` do not survive
it intact: `'G00_PRE   '` reaches the compiler as `'G00_PRE '`, while
`'G03_ZMETAB'` is unaffected because it has no padding.

Edit descriptor `Aw` right-justifies a string shorter than `w`, so writing the
dummy argument straight into `a10` put the tag at a column that varied from tag
to tag and never matched the C++ `%-10s`. Every line then differed under a
plain `diff`, and the whole log looked divergent.

`fennel_tag_line` therefore assigns the tag and variable name into fixed-length
locals (`character(len=10)`, `character(len=4)`) before writing. Character
assignment blank-pads on the right, which restores left justification without
depending on any literal surviving preprocessing.

Do not "simplify" that back to writing the dummies directly, and do not try to
fix it by padding the literals harder — the preprocessor will eat the padding
again.

## Tags

Emitted in this order. `G01`–`G09` are inside `ITER_LOOP` and repeat per
iteration; `G00`, `G10`, `G11` are emitted once per column per call.

| #   | Tag          | Emitted after                                              | `fennel.h` anchor |
|-----|--------------|------------------------------------------------------------|-------------------|
| G00 | `G00_PRE`    | Bio/Bio_old extraction, temp/salt clamp, `PARsur`           | 615               |
| G01 | `G01_LIGHT`  | fused light-attenuation / Eppley / NO3+NH4 uptake / chlorophyll / TIC / nitrification loop | 871 |
| G02 | `G02_GRAZE`  | zooplankton grazing, assimilation, egestion, phyto mortality | 917              |
| G03 | `G03_ZMETAB` | zooplankton basal metabolism, mortality, excretion           | 970              |
| G04 | `G04_COAG`   | coagulation SDeN + Phyt to LDeN                              | 994              |
| G05 | `G05_REMINN` | N-detritus remineralization                                  | 1083             |
| G06 | `G06_O2FLX`  | surface O2 gas exchange (`OXYGEN` only)                      | 1142             |
| G07 | `G07_REMINC` | C-detritus remineralization and diagnostic alkalinity (`CARBON` only) | 1181    |
| G08 | `G08_CO2FLX` | surface CO2 gas exchange (`CARBON` only)                     | 1273             |
| G09 | `G09_SINK`   | `SINK_LOOP` (all sinking constituents)                       | 1521             |
| G10 | `G10_POST`   | `ITER_LOOP` exit, before the global update                   | 1522             |
| G11 | `G11_UPDATE` | `t(nnew)` update; prints the updated tracer, not `Bio`       | 1562             |

`G06` and `G08`/`G07` are absent from both logs when the corresponding option
is off, so `diff` stays aligned.

`G01` fuses what would naturally be four blocks. ROMS computes them in a
single `DO k=N,1,-1` loop with a serial PAR descent, so they cannot be split
without restructuring the Fortran. Narrowing inside `G01` is a Tier-1.5/Tier-2
retreat activity, not a Tier-1 tag.

## Variables

Printed in this order at every tag. Inactive tracers are skipped on **both**
paths, so field order stays contractual across option sets.

| Order | `VAR`  | ROMS index | REMORA component | Active when |
|-------|--------|------------|------------------|-------------|
| 1     | `NO3`  | `iNO3_`    | `bio_comp.no3`   | always      |
| 2     | `NH4`  | `iNH4_`    | `bio_comp.nh4`   | always      |
| 3     | `CHLO` | `iChlo`    | `bio_comp.chlo`  | always      |
| 4     | `PHYT` | `iPhyt`    | `bio_comp.phyt`  | always      |
| 5     | `ZOOP` | `iZoop`    | `bio_comp.zoop`  | always      |
| 6     | `LDEN` | `iLDeN`    | `bio_comp.lden`  | always      |
| 7     | `SDEN` | `iSDeN`    | `bio_comp.sden`  | always      |
| 8     | `PO4`  | `iPO4_`    | `bio_comp.po4`   | `po4`       |
| 9     | `LDEC` | `iLDeC`    | `bio_comp.ldec`  | `carbon`    |
| 10    | `SDEC` | `iSDeC`    | `bio_comp.sdec`  | `carbon`    |
| 11    | `TIC`  | `iTIC_`    | `bio_comp.tic`   | `carbon`    |
| 12    | `TALK` | `iTAlk`    | `bio_comp.talk`  | `carbon`    |
| 13    | `OXYG` | `iOxyg`    | `bio_comp.oxyg`  | `oxygen`    |
| 14    | `ODU`  | `iODU_`    | `bio_comp.odu`   | `odu`       |

ROMS tracer index equals REMORA component + 1 throughout
(`itemp`=1/`Temp_comp`=0, `isalt`=2/`Salt_comp`=1, `iNO3_`=3/`Tracer_comp`=2),
verified against `fennel_mod.h:523-561` and `REMORA_Biology.H:172-191`.

## Debug levels

Two-level gating.

`USE_BIOLOGY_DIAG=TRUE` / `-DREMORA_ENABLE_BIOLOGY_DIAG=ON` decides at **compile
time** whether the emitters exist at all. Off by default: the Path B emitters
sit inside the `AMREX_GPU_DEVICE` column lambda, where a `printf` costs device
print machinery and closure captures even when never reached. With the flag
off, `remora.biology_debug > 0` aborts rather than producing an empty log.

`remora.biology_debug` then selects verbosity at **run time**, threaded through
the `bind(C)` signature as `integer(c_int), value`, so switching paths or
levels within a campaign never requires a rebuild.

This is a deliberate, narrowed deviation from the source corpus, which retired
compile-time diagnostic gating (Rule 30 addendum) on the grounds that runtime
control must reach the Fortran. It still does — the compile flag only decides
whether the code is present, not how it is controlled. WSM6's instrumentation
was CPU-side and did not pay device-print costs. Record the deviation in the
decision ledger rather than treating this paragraph as the authority.

| Level | Scope                                                           |
|-------|-----------------------------------------------------------------|
| 0     | Off (production)                                                  |
| 1     | Tier 1: all tags, all levels, all tracers, target column only     |
| 2     | Tier 1 for every column in the tile                               |

Target column: `remora.biology_debug_i` / `remora.biology_debug_j`
(REMORA 0-based, default `2 2`). Level 2 on anything larger than the BioToy
4x4 grid produces unusable volume; that is intentional.

## Sidecar protocol

```
# Path A (oracle)
remora.use_biology_cpp_answer = 0
remora.biology_debug          = 1
mpirun -n 1 ./REMORA3d.gnu.ex inputs > log.fortran

# Path B (native)
remora.use_biology_cpp_answer = 1
remora.biology_debug          = 1
mpirun -n 1 ./REMORA3d.gnu.ex inputs > log.cpp

grep '^FENNEL-' log.fortran | sed 's/^FENNEL-FORT /FENNEL/' > a.txt
grep '^FENNEL-' log.cpp     | sed 's/^FENNEL-CPP  /FENNEL/' > b.txt
diff a.txt b.txt | head -40
```

Both path prefixes are 12 characters (`FENNEL-FORT ` and `FENNEL-CPP  `), so
the two `sed` patterns above are not interchangeable — note the differing
space counts. With the emitters correct, `diff` reports exactly two lines per
mismatching record and nothing else.

Because the physics is sequential, a divergence at `G01_LIGHT` propagates
through every later tag. Stop at the first divergent tag; do not read past it.
Do not patch logic without first-line evidence at the scoped retreat site.

Serial (`-n 1`) is the default acceptance scope. MPI and GPU parity are
separate lanes and are not opened by this campaign.

### `biology_debug = 2` requires key-based comparison, not `diff`

The plain `diff` above is only valid at `biology_debug = 1`, where a single
column is printed. At level 2 the two paths emit the **same record set in a
different order**, because the tag calls sit in different loop nests:

- Path A calls `fennel_tag` at each ROMS block seam, which is outside the `i`
  loop, so it emits *all columns* for one tag before moving to the next tag.
- Path B calls `tag_state` inside the per-column `ParallelFor`, so it emits
  *all tags* for one column before moving to the next column.

This is inherent to where the seams are and is not a defect. Compare level-2
logs by matching on the `(tag, iter, i, j, k, var)` key, or sort both logs
first. A positional `diff` will report almost every line as changed.

## Known scoped gaps and bridge normalizations

Recorded here so they are not rediscovered as "bugs" mid-campaign.

1. **`BULK_FLUXES`** — BioToy runs `remora.bulk_fluxes = true`, so ROMS proper
   takes the `Uwind`/`Vwind` branch for gas exchange. The bridge is compiled
   with `BULK_FLUXES` and validates the runtime option set before calling
   ROMS, so BioToy compares the native port against the same branch.
2. **`srflx` units** — `fennel.h:614` computes
   `PARsur = PARfrac*srflx*rho0*Cp` (ROMS `srflx` is kinematic, degC m/s);
   `REMORA_Biology.cpp:564` computes `PARsur = PARfrac*srflx` (assumes W/m2).
   The bridge packs `vec_srflx/(rho0*Cp)` into the ROMS buffer while still
   passing the real `rho0` and `Cp` for gas-exchange code, so both paths see
   the same PAR surface forcing.
3. **`pH`** — not a divergence source. `fennel.h` fixes `DoNewton = 0`, and
   the bisection brackets are constants (`pH_hi=10`, `pH_lo=5`), matched by
   `REMORA_Biology.cpp:143-144`. `pH` is write-only output on this path, so
   the bridge passes per-call scratch rather than persistent state.
4. **`z_r`** — a `biology_tile` dummy argument with zero uses in the routine
   body. The bridge passes a correctly shaped buffer for interface fidelity.
