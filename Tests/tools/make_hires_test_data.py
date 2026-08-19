#!/usr/bin/env python3
"""Generate small high-resolution NetCDF input files for the hires_grid_level /
hires_init_level developer lanes.

Writes two files in **64-bit offset classic** (CDF-2) format, matching the ``*_classic64.nc``
convention used throughout REMORA's inputs. PnetCDF cannot read NETCDF4/HDF5, so the format is
not negotiable -- the writer below pins it explicitly and refuses to guess.

    <dir>/hires_grd.nc   h, pm, pn                    (SN_WE)
    <dir>/hires_ini.nc   temp, salt, u, v, zeta       (Time_BT_SN_WE / Time_SN_WE)

Only the variables the hires readers actually consume are written. With both hires levels set,
plus ``remora.use_coriolis = false`` and ``remora.mask_type = none``, the level-0 files are never
dereferenced, so no ``x_rho``, ``f`` or ``mask_*`` are needed anywhere.

The grow-cell rule (Docs/sphinx_doc/Inputs.rst) fixes the dimensions::

    xi_rho  = n_cell_x * r + 2r
    eta_rho = n_cell_y * r + 2r

and REMORA now aborts if a file is smaller than that, so the sizes here are load-bearing.

Fields are deliberately asymmetric: ``pm != pn`` and ``u`` unlike ``v``, so a transposition or a
pm/pn swap in the reader shows up as a wrong answer rather than an identical one.

Usage
-----
    python3 Tests/tools/make_hires_test_data.py --dir Exec/GulfRefinementTest \\
            --n-cell 10 10 4 --ref-ratio 3

This script is a developer tool: it is not wired into CMake and is not run by CI. The generated
files are the artifacts of record -- regenerating them invalidates any baseline taken against
them.
"""

import argparse
import hashlib
import os
import sys

import numpy as np


def open_cdf2(path, dims):
    """Open ``path`` for writing in 64-bit offset classic (CDF-2) and define ``dims``.

    Returns (handle, kind) where kind is 'scipy' or 'netCDF4'; the two have different enough
    APIs for variable creation that the caller has to know which it got.
    """
    try:
        from scipy.io import netcdf_file
        fh = netcdf_file(path, "w", version=2)   # version=2 == 64-bit offset
        for name, size in dims.items():
            fh.createDimension(name, size)
        return fh, "scipy"
    except ImportError:
        pass

    try:
        import netCDF4
    except ImportError:
        sys.exit("need either scipy or netCDF4 to write NetCDF; install one "
                 "(scipy needs only numpy and writes CDF-2 natively)")
    fh = netCDF4.Dataset(path, "w", format="NETCDF3_64BIT_OFFSET")
    for name, size in dims.items():
        fh.createDimension(name, size)
    return fh, "netCDF4"


def put(fh, kind, name, dims, data, **attrs):
    if kind == "scipy":
        var = fh.createVariable(name, "d", dims)
        var[:] = data
    else:
        var = fh.createVariable(name, "f8", dims)
        var[:] = data
    for key, value in attrs.items():
        setattr(var, key, value)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=".", help="directory to write the two files into")
    ap.add_argument("--n-cell", nargs=3, type=int, default=[10, 10, 4],
                    metavar=("NX", "NY", "NZ"), help="level-0 remora.n_cell")
    ap.add_argument("--ref-ratio", type=int, default=3,
                    help="cumulative refinement ratio at the hires level")
    ap.add_argument("--dx", type=float, default=512.0,
                    help="level-0 cell size in metres (dyadic by default, so pm is exact)")
    ap.add_argument("--with-biology", action="store_true",
                    help="also write the seven Fennel tracers as constants, matching the "
                         "analytic BioToy profile so expected values stay exact")
    args = ap.parse_args()

    nx, ny, nz = args.n_cell
    r = args.ref_ratio
    if r < 1:
        sys.exit("--ref-ratio must be >= 1")

    # Grow-cell rule: the file covers the refined domain plus r rings on each side.
    xi_rho = nx * r + 2 * r
    eta_rho = ny * r + 2 * r
    xi_u = xi_rho - 1
    eta_v = eta_rho - 1
    dx_f = args.dx / r          # fine-level cell size

    os.makedirs(args.dir, exist_ok=True)

    # Cell-centre coordinates of the refined grid, including the grow rings, measured from the
    # lower-left corner of the level-0 domain.
    xi = (np.arange(xi_rho) - r + 0.5) * dx_f
    eta = (np.arange(eta_rho) - r + 0.5) * dx_f
    X, Y = np.meshgrid(xi, eta)                       # (eta, xi), the NetCDF axis order
    lx, ly = nx * args.dx, ny * args.dx

    # --- grid file ------------------------------------------------------------------------
    # A smooth bowl, non-separable in x and y, strictly positive and well away from zero.
    h = 120.0 + 380.0 * (1.0 - np.exp(-(((X - 0.5 * lx) / (0.30 * lx)) ** 2
                                        + ((Y - 0.5 * ly) / (0.45 * ly)) ** 2)))
    # pm and pn differ from each other by construction: a swap must not go unnoticed. The +-6%
    # modulation also exercises the curvilinear terms derived from the metrics.
    pm = (1.0 / dx_f) * (1.0 + 0.06 * np.sin(2.0 * np.pi * X / lx))
    pn = (1.0 / dx_f) * (1.0 - 0.06 * np.cos(2.0 * np.pi * Y / ly))

    grd = os.path.join(args.dir, "hires_grd.nc")
    fh, kind = open_cdf2(grd, {"eta_rho": eta_rho, "xi_rho": xi_rho})
    put(fh, kind, "h", ("eta_rho", "xi_rho"), h, units="meter",
        long_name="bathymetry at RHO-points")
    put(fh, kind, "pm", ("eta_rho", "xi_rho"), pm, units="meter-1",
        long_name="curvilinear coordinate metric in XI")
    put(fh, kind, "pn", ("eta_rho", "xi_rho"), pn, units="meter-1",
        long_name="curvilinear coordinate metric in ETA")
    fh.close()

    # --- initial file ---------------------------------------------------------------------
    # Stratified in k with a horizontal front, so the pressure gradient sees something real.
    k = np.arange(nz)
    zfrac = (k + 0.5) / nz                            # 0 at the bottom, 1 at the surface
    temp = (10.0 + 8.0 * zfrac[:, None, None]
            + 1.5 * np.tanh((X[None, :, :] - 0.5 * lx) / (0.15 * lx)))
    salt = 35.0 + 0.25 * np.sin(2.0 * np.pi * Y[None, :, :] / ly) * (1.0 - zfrac[:, None, None])
    # u and v are NOT mirror images of each other, so the staggered box types are validated.
    u = 0.05 * np.tanh((Y[None, :, :xi_u] - 0.5 * ly) / (0.20 * ly)) * zfrac[:, None, None]
    v = -0.03 * np.sin(2.0 * np.pi * X[None, :eta_v, :] / lx) * (0.5 + 0.5 * zfrac[:, None, None])
    zeta = 0.02 * np.sin(2.0 * np.pi * (X + Y) / (lx + ly))

    ini = os.path.join(args.dir, "hires_ini.nc")
    dims = {"ocean_time": 1, "s_rho": nz, "eta_rho": eta_rho, "xi_rho": xi_rho,
            "xi_u": xi_u, "eta_v": eta_v}
    fh, kind = open_cdf2(ini, dims)
    put(fh, kind, "ocean_time", ("ocean_time",), np.zeros(1),
        units="seconds since 2000-01-01 00:00:00")
    put(fh, kind, "temp", ("ocean_time", "s_rho", "eta_rho", "xi_rho"), temp[None, ...],
        units="Celsius", long_name="potential temperature")
    put(fh, kind, "salt", ("ocean_time", "s_rho", "eta_rho", "xi_rho"), salt[None, ...],
        units="PSU", long_name="salinity")
    put(fh, kind, "u", ("ocean_time", "s_rho", "eta_rho", "xi_u"), u[None, ...],
        units="meter second-1", long_name="u-momentum component")
    put(fh, kind, "v", ("ocean_time", "s_rho", "eta_v", "xi_rho"), v[None, ...],
        units="meter second-1", long_name="v-momentum component")
    put(fh, kind, "zeta", ("ocean_time", "eta_rho", "xi_rho"), zeta[None, ...],
        units="meter", long_name="free-surface")

    if args.with_biology:
        # The analytic BioToy profile (Source/Prob/REMORA_InitAnalyticBiology_BioToy.H) sets six
        # of these as constants. Writing the same constants here keeps the expected level-0
        # values exact after the average-down, so a NetCDF-sourced biology lane can assert on
        # values rather than on a snapshot. NO3 there is a function of temperature; a constant
        # stands in for it, which is all the plumbing test needs.
        for name, value in (("NO3", 2.0), ("NH4", 0.1), ("chlorophyll", 0.02),
                            ("phytoplankton", 0.08), ("zooplankton", 0.06),
                            ("LdetritusN", 0.02), ("SdetritusN", 0.04)):
            put(fh, kind, name, ("ocean_time", "s_rho", "eta_rho", "xi_rho"),
                np.full((1, nz, eta_rho, xi_rho), value),
                long_name="Fennel tracer (constant, for exact assertions)")
    fh.close()

    print("wrote 64-bit offset classic (CDF-2) files for n_cell = %d %d %d, ref_ratio = %d"
          % (nx, ny, nz, r))
    print("  required dims: xi_rho = %d*%d + 2*%d = %d, eta_rho = %d"
          % (nx, r, r, xi_rho, eta_rho))
    for path in (grd, ini):
        print("  %-40s %8.1f KB  sha256=%s"
              % (path, os.path.getsize(path) / 1024.0, sha256(path)))


if __name__ == "__main__":
    main()
