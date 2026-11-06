#!/usr/bin/env python3
"""Generate level-0 NetCDF inputs of any size, for timing REMORA's NetCDF readers.

Writes, for a domain of NX x NY x NZ cells, ROMS-style files with every variable the
level-0 readers take:

    <dir>/scaling_grd.nc   pm, pn, h, f, angle, mask_rho, x/y at rho, u, v and psi points
    <dir>/scaling_ini.nc   temp, salt, u, v (Time_BT_SN_WE); ubar, vbar, zeta (Time_SN_WE)
    <dir>/scaling_frc.nc   Uwind, Vwind on two wind_time records

Fields vary in both horizontal directions and in depth, so a reader that places data on
the wrong cells or ranks changes the answer. Sizes follow ROMS: xi_rho = NX + 2,
eta_rho = NY + 2, one ring of ghost cells.

Run them with Exec/IdealMiniGrid/inputs_frc, overriding the grid, the walls and the files on
the command line. For example, 512 x 512 x 32 cells:

    python3 Tests/tools/make_read_scaling_data.py --dir scaling --n-cell 512 512 32
    python3 Tests/tools/make_read_scaling_data.py --dir scaling --n-cell 512 512 32 \\
            --format NETCDF4 --suffix _nc4z

    remora_exec Exec/IdealMiniGrid/inputs_frc remora.n_cell=512 512 32 \\
        remora.prob_lo=0 0 -5000 remora.prob_hi=512000 512000 0 \\
        remora.bc.xlo.type=slipwall remora.bc.xhi.type=slipwall \\
        remora.bc.ylo.type=slipwall remora.bc.yhi.type=slipwall \\
        remora.fixed_dt=60 remora.ndtfast=30 \\
        remora.nc_grid_file_0=scaling/scaling_grd.nc remora.nc_init_file_0=scaling/scaling_ini.nc \\
        remora.nc_frc_file=scaling/scaling_frc.nc

prob_hi is n_cell times --dx (1000 m by default) and prob_lo's depth is --depth. Keep the
barotropic step fixed_dt/ndtfast below about dx / sqrt(9.8 * depth) (4.5 s for the defaults):
with a larger one the run blows up without stopping, and its plotfiles say nothing. With
remora.max_step=0 a run reads every file and stops, which isolates the cost of reading.

Formats: NETCDF3_64BIT_DATA (CDF-5, the default), NETCDF4 (compressed with zlib level 4
and shuffle), or NETCDF4_CLASSIC. Needs numpy and netCDF4 (python).

This script is a developer tool: it is not wired into CMake and is not run by CI.
"""
import argparse
import os
import sys

import numpy as np
import netCDF4


def create(path, fmt, dims):
    ds = netCDF4.Dataset(path, "w", format=fmt)
    for name, size in dims.items():
        ds.createDimension(name, size)
    return ds


def var(ds, name, dims, data, fmt, units=None):
    kw = {"zlib": True, "complevel": 4, "shuffle": True} if fmt == "NETCDF4" else {}
    v = ds.createVariable(name, "f8", dims, **kw)
    if units is not None:
        v.units = units
    v[:] = data
    return v


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=".", help="directory to write the files into")
    ap.add_argument("--n-cell", nargs=3, type=int, required=True, metavar=("NX", "NY", "NZ"),
                    help="remora.n_cell")
    ap.add_argument("--dx", type=float, default=1000.0, help="cell size in metres")
    ap.add_argument("--depth", type=float, default=5000.0, help="maximum depth in metres")
    ap.add_argument("--format", default="NETCDF3_64BIT_DATA",
                    choices=["NETCDF3_64BIT_DATA", "NETCDF4", "NETCDF4_CLASSIC"])
    ap.add_argument("--suffix", default="", help="appended to each file name")
    args = ap.parse_args()

    nx, ny, nz = args.n_cell
    dx, fmt = args.dx, args.format
    if min(nx, ny, nz) < 2:
        sys.exit("--n-cell must be at least 2 in each direction")
    os.makedirs(args.dir, exist_ok=True)

    # Index grids at each staggering; i, j count from the first ghost cell as in ROMS.
    jr, ir = np.indices((ny + 2, nx + 2), dtype=float)
    ju, iu = np.indices((ny + 2, nx + 1), dtype=float)
    jv, iv = np.indices((ny + 1, nx + 2), dtype=float)
    jp, ip = np.indices((ny + 1, nx + 1), dtype=float)
    Lx, Ly = nx * dx, ny * dx
    xr, yr = (ir - 0.5) * dx, (jr - 0.5) * dx

    grid_dims = {"xi_rho": nx + 2, "eta_rho": ny + 2, "xi_u": nx + 1, "eta_u": ny + 2,
                 "xi_v": nx + 2, "eta_v": ny + 1, "xi_psi": nx + 1, "eta_psi": ny + 1}
    rho, u2, v2 = ("eta_rho", "xi_rho"), ("eta_u", "xi_u"), ("eta_v", "xi_v")
    psi = ("eta_psi", "xi_psi")

    path = os.path.join(args.dir, "scaling_grd" + args.suffix + ".nc")
    with create(path, fmt, grid_dims) as ds:
        var(ds, "pm", rho, np.full(xr.shape, 1.0 / dx), fmt, "meter-1")
        var(ds, "pn", rho, np.full(xr.shape, 1.0 / dx), fmt, "meter-1")
        h = args.depth * (0.6 + 0.25 * np.sin(2 * np.pi * xr / Lx) * np.cos(np.pi * yr / Ly))
        var(ds, "h", rho, h, fmt, "meter")
        var(ds, "f", rho, 1.0e-4 + 1.0e-9 * yr / dx, fmt, "second-1")
        var(ds, "angle", rho, np.zeros(xr.shape), fmt, "radians")
        var(ds, "mask_rho", rho, np.ones(xr.shape), fmt)
        var(ds, "x_rho", rho, xr, fmt, "meter")
        var(ds, "y_rho", rho, yr, fmt, "meter")
        var(ds, "x_u", u2, iu * dx, fmt, "meter")
        var(ds, "y_u", u2, (ju - 0.5) * dx, fmt, "meter")
        var(ds, "x_v", v2, (iv - 0.5) * dx, fmt, "meter")
        var(ds, "y_v", v2, jv * dx, fmt, "meter")
        var(ds, "x_psi", psi, ip * dx, fmt, "meter")
        var(ds, "y_psi", psi, jp * dx, fmt, "meter")
    print(path)

    k = np.arange(nz, dtype=float).reshape(1, nz, 1, 1) / nz
    ini_dims = dict(grid_dims, s_rho=nz, ocean_time=None)
    path = os.path.join(args.dir, "scaling_ini" + args.suffix + ".nc")
    with create(path, fmt, ini_dims) as ds:
        t = ds.createVariable("ocean_time", "f8", ("ocean_time",))
        t.units = "days since 2020-01-01 00:00:00"
        t[:] = [0.0]
        r3 = ("ocean_time", "s_rho") + rho
        var(ds, "temp", r3, 10.0 + 8.0 * k + 0.5 * np.sin(2 * np.pi * xr / Lx)[None, None]
                            * np.cos(2 * np.pi * yr / Ly)[None, None], fmt, "Celsius")
        var(ds, "salt", r3, 34.0 + 0.5 * k + 0.1 * np.cos(2 * np.pi * xr / Lx)[None, None],
            fmt, "PSU")
        var(ds, "u", ("ocean_time", "s_rho") + u2,
            0.01 * np.sin(np.pi * ju / (ny + 1))[None, None] * (1 - k), fmt, "meter second-1")
        var(ds, "v", ("ocean_time", "s_rho") + v2,
            -0.01 * np.sin(np.pi * iv / (nx + 1))[None, None] * (1 - k), fmt, "meter second-1")
        var(ds, "ubar", ("ocean_time",) + u2, np.zeros((1,) + iu.shape), fmt, "meter second-1")
        var(ds, "vbar", ("ocean_time",) + v2, np.zeros((1,) + iv.shape), fmt, "meter second-1")
        var(ds, "zeta", ("ocean_time",) + rho,
            0.01 * np.sin(2 * np.pi * xr / Lx)[None] * np.sin(np.pi * yr / Ly)[None], fmt, "meter")
    print(path)

    frc_dims = dict(grid_dims, wind_time=None)
    path = os.path.join(args.dir, "scaling_frc" + args.suffix + ".nc")
    with create(path, fmt, frc_dims) as ds:
        t = ds.createVariable("wind_time", "f8", ("wind_time",))
        t.units = "days since 2020-01-01"
        t[:] = [0.0, 1.0]
        w = np.stack([5.0 * np.cos(np.pi * yr / Ly), 6.0 * np.cos(np.pi * yr / Ly)])
        var(ds, "Uwind", ("wind_time",) + rho, w, fmt, "meter second-1")
        var(ds, "Vwind", ("wind_time",) + rho, 0.2 * w, fmt, "meter second-1")
    print(path)


if __name__ == "__main__":
    main()
