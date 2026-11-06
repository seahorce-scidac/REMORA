# Plan: NetCDF-4 read/write in REMORA via PnetCDF

The plan is to go entirely through PnetCDF. PnetCDF 1.14.0 has a NetCDF-4
driver that forwards its calls to netCDF-C's parallel API (`nc_open_par`,
`nc_create_par`). Once PnetCDF is built with that driver, REMORA's existing
`ncmpi_*` code can open NetCDF-4 files, and parallel reads and writes both
work. Reading needs almost no code changes. Writing needs the history writer's
write path rebuilt, because it relies on calls the driver doesn't support.

## What was checked

Sources: `~/Downloads/pnetcdf-1.14.0/src/drivers/nc4io`,
`doc/README.NetCDF4.md`, and REMORA's `ncmpi_*` usage.

- **The local install has the driver off.** `~/PnetCDF` reports
  `enable_netcdf4="0"` and `pnetcdf.h` has `PNETCDF_DRIVER_NETCDF4 0`. The
  `/usr/local` netCDF-C 4.9.3 is serial (`--has-parallel4 no`).
  `/usr/bin/h5pcc` exists, so parallel HDF5 may already be installed.
- **Opening needs no code change.** `ncmpi_open` checks the file's format and
  picks the right driver itself.
- **What the driver supports:**
  - Supported: `begin/end_indep_data` (it sets each variable's parallel
    access mode), `def_var_fill`, attributes, inquiry calls, and
    `get/put_vara/vars[_all]`.
  - Not supported: `varn`, all nonblocking calls (`iput`/`wait_all`),
    `vard`, and chunking or compression (PnetCDF has no API for them).
- **Driver rules REMORA's writer must follow:**
  - Every rank in one collective call has to use the same kind of call.
  - A new file can have only one unlimited dimension.
  - Reads see only the root group, and NC_STRING variables can't be read.
- **REMORA's current usage:**
  - Every input read (init, clim/forcing time series, boundary, river,
    lat-lon regrid) has rank 0 open the file and read in independent mode,
    then broadcast. All of that is supported.
  - The history writer (`Source/IO/REMORA_NCPlotFile.cpp`) writes each
    variable with one `ncmpi_put_varn_*_all` call, which returns
    `NC_ENOTSUPPORT` on NetCDF-4.
  - `NCVar::iput` and `NCFile::wait_all` are wrappers with no callers.

## Phase 0: Toolchain

1. Build or locate parallel HDF5 (1.10.4 or newer), then netCDF-C with
   `--enable-parallel4`. Don't build netCDF-C with PnetCDF support; the
   PnetCDF README says that combination isn't supported.
2. Rebuild PnetCDF 1.14.x with `--enable-netcdf4 --with-netcdf4=<path>`.
3. Linking:
   - CMake: `CMake/FindPNetCDF.cmake` uses pkg-config. Shared builds should
     link as is. For static builds, check whether `Libs.private` picks up
     `-lnetcdf -lhdf5_hl -lhdf5`; if not, add them.
   - GNU Make: check whether `Exec/Make.REMORA` needs the same libraries
     added.
4. Perlmutter: check whether `cray-parallel-netcdf` was built with the
   NetCDF-4 driver. If it wasn't, document a self-built PnetCDF that sits on
   top of `cray-hdf5-parallel` and `cray-netcdf-hdf5parallel`.

### Phase 0 status (local machine)

Done. Ubuntu's `libhdf5-mpi-dev` (HDF5 1.10.10) and `libnetcdf-mpi-dev`
(netCDF-C 4.9.0, `NC_HAS_PARALLEL4 1`, `NC_HAS_PNETCDF 0`) are enough; nothing
had to be built besides PnetCDF. PnetCDF 1.14.0 is installed with the driver at
`~/PnetCDF-nc4` (source in `~/software/pnetcdf-nc4-src`), configured with

    ./configure --prefix=$HOME/PnetCDF-nc4 \
        --with-netcdf4=/usr/lib/x86_64-linux-gnu/netcdf/mpi \
        CPPFLAGS="-I/usr/lib/x86_64-linux-gnu/netcdf/mpi/include -I/usr/include/hdf5/openmpi" \
        LDFLAGS="-L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi" \
        LIBS=-lnetcdf_mpi --disable-fortran

Two traps found:
- Without `LIBS=-lnetcdf_mpi`, `libpnetcdf.so` leaves `nc_open_par` and the
  other `nc_*` symbols undefined and records no dependency (`Libs.private` is
  empty), so an application link could pick up the serial `/usr/local`
  netCDF-C.
- The shell profile puts `~/PnetCDF/lib` on `LD_LIBRARY_PATH`, which beats the
  executable's RUNPATH, so a REMORA built against `~/PnetCDF-nc4` still loads
  the driverless library unless `LD_LIBRARY_PATH` is changed.

Perlmutter is still open.

## Phase 1: Reading (small PR)

### Phase 1 findings

- **Compressed NetCDF-4 variables cannot be read in independent mode.**
  netCDF-C (`libsrc4/nc4var.c`, 4.9.2 and 4.9.3) refuses independent access to
  any variable with a zlib, shuffle or fletcher32 filter, reads included.
  PnetCDF's `begin_indep_data` then fails at the first such variable (REMORA
  ignored the return code), and the I/O rank's independent read aborts with
  "Operation not allowed in collective data mode". netCDF-C also turns on
  HDF5's collective metadata reads for every parallel open.
- So the readers now open each file on the I/O rank alone
  (`NCFile::open_local`, `MPI_COMM_SELF`) and read it with the collective calls,
  which are trivially collective on one rank. This replaces the plan's item 4
  ("leave the rank-0-reads-and-broadcasts pattern alone"): the pattern stays,
  but the other ranks no longer open the file.
- Numeric type conversion needed no work: every read goes through
  `ncmpi_get_var*_double/float/int`, which convert from any numeric type.
  NC_STRING attributes cannot be read through PnetCDF at all, so they now abort
  with the attribute, variable and file named.


1. **Clear error when the driver is missing.** In `NCFile::open`, call
   `ncmpi_inq_file_format` first. If the file is `NC_FORMAT_NETCDF4` or
   `NC_FORMAT_NETCDF4_CLASSIC` and `PNETCDF_DRIVER_NETCDF4 == 0`, abort with
   "rebuild PnetCDF with `--enable-netcdf4`" instead of failing with
   `NC_ENOTBUILT`.
2. **Open read-only.** Read paths currently open files with `NC_NOCLOBBER`,
   which is a create flag. Switch them to `NC_NOWRITE` so the netCDF-C open
   doesn't depend on how it treats stray flags.
3. **Data types.** Add readers for the types NetCDF-4 files tend to contain
   (`NC_UBYTE`, `NC_INT64`, and so on) wherever `get_var` switches on the
   variable's type. Give a clear error for NC_STRING and for variables inside
   groups.
4. **Leave the rank-0-reads-and-broadcasts pattern alone.** It works
   unchanged, and it avoids HDF5's strict collective rules. Truly collective
   reads are a separate, format-independent optimisation.
5. **Tests:**
   - Use `nccopy -k nc4` and `-k nc7` to make NetCDF-4 and NetCDF-4-classic
     copies of the existing NetCDF inputs (grid/init, clim, forcing,
     boundary, river, OCCAMS lat-lon).
   - Check the plotfiles match the CDF-5 runs bit for bit with `fcompare`, on
     1, 2 and 4 ranks.
   - Add a ctest gated on a CMake check for the driver.
   - The NetCDF ctests don't cover NCTimeSeries (forcing/clim) reads, so
     check those by hand with the IdealMiniGrid/BioToy recipe.

## Phase 1b: Parallel reads (done)

Phase 1 left every read on the I/O rank, followed by a broadcast of the whole
field to every rank. Now:

- `BuildLocalFABsFromNetCDFFile` (`Source/IO/REMORA_NCFile.H`) opens the file
  on every rank. Each rank reads, in one collective `get_vara_all` per
  variable, the hyperslab covering the region it needs, into FABs in the same
  index space as `BuildFABsFromNetCDFFile`. So the callers' copy on
  intersection is unchanged. One hyperslab per rank (the bounding box of its
  FABs) avoids `varn`, which the NetCDF-4 driver lacks. A rank with no boxes
  reads nothing but still takes part in the collective call.
- `nc_local_region` gives that region from the MultiFabs a caller fills: the
  cell-centered bounding box of the FABs the rank owns. A rank whose FABs lie
  far apart reads everything between them. AMReX's default space-filling-curve
  distribution keeps a rank's boxes together. If measurements show otherwise,
  the fix is one read per FAB, with every rank making the same number of calls.
- Switched: `NCTimeSeries` (forcing and climatology records) and every
  init/grid reader in `REMORA_ReadFromInitNetcdf.cpp`, including the hires
  readers. The hires full-domain MultiFabs are one box per level, so on that
  path the rank that owns the box reads the whole field and the others read
  nothing. `unwrap_to_first_cell` relies on that and asserts it.
- Still on one rank: boundary files (`NCTimeSeriesBoundary`, small slabs),
  rivers, lat-lon forcing (its grid is not the model's decomposition), and
  metadata queries.

Verification: bit-identical plotfiles against the base-commit binary on CDF-5
and compressed NetCDF-4 inputs, at 1, 2 and 4 ranks and with 9 boxes on 3 and
4 ranks. Decks used: IdealMiniGrid `inputs`, `inputs_frc`, `inputs_clim_nudg`
(with a perturbed climatology), `inputs_synth_hires`, and BioToy. Negative
controls (empty or shrunken regions) fail. The `NCFormats` unit test now
compares the two readers cell for cell on a 6-box MultiFab, for a record
variable, a static one and a face-centered one; a one-cell shrink fails it.

The same runs also pass in a Debug build with assertions, `amrex.init_snan=1`
and FPE traps.

Local measurement (`Tests/tools/make_read_scaling_data.py`, 768 x 768 x 16
cells, start-up only, 4 ranks): peak RSS per rank 980-1047 MB with the
branch against 1152-1441 MB at the base commit. Start-up took about 2.0 s
against 2.5 s, but those timings overlapped a build. At 1 rank memory is the
same, as it should be. The real test is at scale on Perlmutter, below.

### Perlmutter test

1. **Modules.** Use `PrgEnv-gnu`, `cray-mpich`, `cray-hdf5-parallel` and
   `cray-netcdf-hdf5parallel`, and do not load `cray-parallel-netcdf`.
   Check `nc-config --has-parallel4` (must be `yes`) and
   `nc-config --has-pnetcdf`. If the second is `yes`, PnetCDF's README says the
   combination is unsupported; try it anyway, and if it fails, build netCDF-C
   with `--enable-parallel4 --disable-pnetcdf` against `cray-hdf5-parallel`.
2. **PnetCDF with the driver**, into `$PNC=$CFS/<project>/$USER/pnetcdf-nc4`
   or `$HOME`. Configure with `CC=cc CXX=CC --with-netcdf4=$NETCDF_DIR
   --disable-fortran --prefix=$PNC`. Check that the configure summary says
   `NetCDF-4 support - enabled`, that `$PNC/bin/pnetcdf-config --netcdf4`
   prints `enabled`, and that `ldd $PNC/lib/libpnetcdf.so` lists `libnetcdf`.
   If it doesn't, rebuild with `LIBS=-lnetcdf`.
3. **REMORA** with CMake, so the unit test is built: put `$PNC/lib/pkgconfig`
   first on `PKG_CONFIG_PATH` and `$PNC/lib` first on `LD_LIBRARY_PATH`, then
   configure with `-DREMORA_ENABLE_MPI=ON -DREMORA_ENABLE_PNETCDF=ON
   -DREMORA_ENABLE_TESTS=ON -DREMORA_ENABLE_FCOMPARE=ON`. Check with
   `ldd Build/Exec/remora_exec | grep pnetcdf`. Build the base commit
   (`a54d74b7`) the same way, in a second directory, as the reference.
4. **Correctness**, in an interactive CPU allocation: `ctest -R NCFormats`.
   Then compare the IdealMiniGrid and BioToy decks, the base binary on
   `_classic64` inputs against the branch on CDF-5 and `nccopy -k nc4 -d 4 -s`
   copies, with `amrex_fcompare`, at 4, 64 and 256 ranks. Use a perturbed
   climatology for `inputs_clim_nudg`.
5. **Scale**: generate 1536 x 1536 x 32 (and 3072 x 3072 x 32) inputs with
   `make_read_scaling_data.py`, in CDF-5 and NetCDF-4 (`--format NETCDF4`).
   Then run `inputs_frc` with `remora.max_step=0` at 128, 512 and 2048
   ranks (1, 4 and 16 nodes), base versus branch. Record the wall time and
   each rank's peak RSS (`srun ... /usr/bin/time -f "%M %e"`). Put the input
   files on `$SCRATCH`, striped (`lfs setstripe -c 16`).

Found along the way, all present at the base commit, not fixed:
- `inputs_clim_nudg`'s climatology equals its initial state (and is constant
  in time for temp, u, ubar), so that deck cannot detect a wrong climatology
  read. Make a perturbed copy to test with.
- `amr.max_grid_size=4` (which also splits the vertical) segfaults, and
  `inputs_clim_nudg` with 9 column boxes on 3 or 4 ranks goes NaN by step 100.
- `inputs_synth_hires` with small boxes at the hires level segfaults.

## Phase 2: Writing (main work)

1. **Add a `remora.nc_format` option** taking `cdf5` (the default) or
   `netcdf4`. With `netcdf4` the file is created with
   `NC_CLOBBER|NC_NETCDF4`, plus `NC_CLASSIC_MODEL` so files stay readable by
   classic tools. Document it in the Sphinx docs (`Plotfiles.rst`).
2. **Replace the `varn` path in `REMORA_NCPlotFile.cpp`** (around line 106 and
   the flush at line 1791). Each rank writes a different number of boxes, so:
   1. Find the largest box count across ranks with `MPI_Allreduce`.
   2. Every rank makes that many `put_vara_*_all` calls, filling its spare
      calls with zero-length requests.
   3. Every rank uses the same call kind, to satisfy the driver's rule.

   Keep `varn` for CDF-5, where it's faster. Hide the choice inside `NCVar`,
   for example a `put_boxes_all()` method, so the writer keeps one code path.
3. **Audit writes in independent mode.**
   - Coordinate writes on rank 0 (around line 937) are fine for fixed-size
     variables.
   - Any independent write that grows the unlimited `ocean_time` dimension is
     a problem, because netCDF-4/HDF5 can only grow a dataset collectively.
     `ocean_time` itself already uses `put_all`, but every record variable
     needs checking.
4. **Appending to existing files** (the `NC_WRITE` open around line 256, and
   `ncmpi_redef`) works through the driver. Test that a restart-and-append
   run gives the same file as one continuous run.
5. **Delete `iput`/`wait_all`.** They have no callers and would fail on
   NetCDF-4.
6. **Tests:**
   - The history writer has no ctest. Compare `cdf5` and `netcdf4` output by
     hand with `ncmpidiff` (it reads both formats once the driver is built) on
     1, 2 and 4 ranks, including multilevel output.
   - Time both on a medium case. HDF5's collective metadata cost and the extra
     zero-length calls will make it slower; check by how much.

## Phase 3 (optional): compression

PnetCDF can't set chunking or deflate, so Phase 2 output is uncompressed,
contiguous HDF5. The gain is compatibility with tools that expect NetCDF-4,
not smaller files. Parallel compressed output needs a small direct netCDF-C
writer: `nc_create_par`, then `nc_def_var_chunking` and `nc_def_var_deflate`,
then collective `nc_put_vara`. HDF5 supports parallel filtered writes from
1.10.2 if they're collective. That gives up "PnetCDF only", so decide on it
after Phase 2 is measured.

## Risks

- **Linking:** static builds on Perlmutter are where the extra NetCDF and HDF5
  libraries most likely go missing.
- **Speed of repeated reads:** opening an HDF5 file is slower than opening
  CDF-5, and NCTimeSeries reopens files every time it reads. That matters if
  forcing is read often.
- **Unreadable files:** NetCDF-4 inputs that use groups or strings (some newer
  datasets) still won't read through PnetCDF. The fallback is to convert them
  once with `nccopy -k cdf5`, which is also the workaround available today.

## Suggested order

Do Phase 0 and Phase 1 first as one small PR. Phase 2 is where the real effort
and the speed questions are.
