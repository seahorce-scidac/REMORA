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

## Phase 1: Reading (small PR)

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
