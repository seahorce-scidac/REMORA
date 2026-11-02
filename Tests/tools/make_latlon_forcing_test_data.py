#!/usr/bin/env python3
"""Make and check synthetic lat-lon forcing for REMORA's forcing regrid (#680).

Two subcommands:

    make   write a forcing file on a coarse rectilinear lat-lon grid covering a model grid
           file, and optionally a copy of the grid file with a nonzero, varying ``angle``
    check  compare the forcing REMORA wrote to its NetCDF history file with the fields the
           make step wrote, interpolated and rotated exactly

The fields are chosen so that REMORA's interpolation and rotation reproduce them to roundoff,
which makes ``check`` an exact test rather than a tolerance on a resolution error:

    Tair  = 10 + 0.3 L + 0.5 P + 0.02 L P                  (bilinear)
    Uwind = 1 + 0.1 L + 0.2 P,      Vwind = 2 - 0.05 L + 0.1 P        (east, north)
    sustr = 1e-4 (0.1 + 0.01 L - 0.005 P),  svstr = 1e-4 (-0.05 + 0.002 L + 0.008 P)

with P the latitude and L the longitude taken within 180 degrees of the model grid's mean.
Winds and stresses are given east/north and REMORA rotates them by the grid ``angle``, so the
expected values are rotated the same way; the stress is then averaged to the u and v faces.

Usage, from a directory holding the IdealMiniGrid files of the remora-data repository
(see Exec/IdealMiniGrid/inputs_frc_latlon):

    python3 Tests/tools/make_latlon_forcing_test_data.py make \\
            --grid idmini_grd_v1_classic64.nc --out idmini_frc_latlon.nc \\
            --rotated-grid idmini_grd_rotated.nc
    remora_exec Exec/IdealMiniGrid/inputs_frc_latlon
    python3 Tests/tools/make_latlon_forcing_test_data.py check plt_latlon_his_d01.nc

Files are written as 64-bit offset classic (CDF-2), which PnetCDF reads. Needs numpy and
netCDF4. This is a developer tool: it is not wired into CMake and is not run by CI.
"""

import argparse
import sys

import numpy as np

try:
    import netCDF4
except ImportError:
    sys.exit("make_latlon_forcing_test_data.py needs the netCDF4 Python package")


def continuous_lon(lon, center):
    """Longitude taken within 180 degrees of center"""
    return lon - 360.0 * np.round((lon - center) / 360.0)


def fields(L, P):
    """The east/north fields the file holds, as functions of continuous lon L and lat P"""
    return {
        "Tair":  10 + 0.3 * L + 0.5 * P + 0.02 * L * P,
        "Uwind": 1 + 0.1 * L + 0.2 * P,
        "Vwind": 2 - 0.05 * L + 0.1 * P,
        "sustr": 1e-4 * (0.1 + 0.01 * L - 0.005 * P),
        "svstr": 1e-4 * (-0.05 + 0.002 * L + 0.008 * P),
    }


def make(args):
    grid = netCDF4.Dataset(args.grid)
    # Extents from the points that hold coordinates: a fill value would otherwise be taken
    # for one
    lon_rho = np.ma.masked_invalid(grid["lon_rho"][:]).compressed()
    lat_rho = np.ma.masked_invalid(grid["lat_rho"][:]).compressed()
    if lon_rho.size == 0 or lat_rho.size == 0:
        sys.exit(f"{args.grid} has no valid lon_rho/lat_rho")
    center = float(np.mean(lon_rho))
    L_grid = continuous_lon(lon_rho, center)

    # Source axes covering the grid, its boundary ring included, with a margin
    d, m = args.spacing, args.margin
    lon = np.arange(L_grid.min() - m, L_grid.max() + m + 0.5 * d, d)
    lat = np.arange(lat_rho.min() - m, lat_rho.max() + m + 0.5 * d, d)
    LON, LAT = np.meshgrid(lon, lat)
    f = fields(LON, LAT)

    lon_store = lon % 360.0 if args.lon360 else lon
    if args.flip_lat:
        lat = lat[::-1]
        f = {k: v[::-1, :] for k, v in f.items()}

    out = netCDF4.Dataset(args.out, "w", format="NETCDF3_64BIT_OFFSET")
    out.createDimension("wind_time", 2)
    out.createDimension("lat", len(lat))
    out.createDimension("lon", len(lon_store))
    t = out.createVariable("wind_time", "f8", ("wind_time",))
    t[:] = [0.0, 3650.0]
    t.units = "days since 2020-01-01"
    v = out.createVariable("lon", "f8", ("lon",)); v[:] = lon_store; v.units = "degrees_east"
    v = out.createVariable("lat", "f8", ("lat",)); v[:] = lat; v.units = "degrees_north"
    for name, data in f.items():
        v = out.createVariable(name, "f8", ("wind_time", "lat", "lon"))
        v[:] = np.broadcast_to(data, (2,) + data.shape)
        v.coordinates = "lon lat"
    out.close()
    print(f"wrote {args.out}: {len(lon_store)} x {len(lat)} lat-lon grid, "
          f"lon {lon_store[0]:g} to {lon_store[-1]:g}, lat {lat[0]:g} to {lat[-1]:g}")

    if args.rotated_grid:
        ny, nx = grid["lon_rho"].shape
        J, I = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
        angle = 0.1 + 0.01 * I + 0.02 * J
        rot = netCDF4.Dataset(args.rotated_grid, "w", format="NETCDF3_64BIT_OFFSET")
        rot.setncatts({a: grid.getncattr(a) for a in grid.ncattrs()})
        for name, dim in grid.dimensions.items():
            rot.createDimension(name, None if dim.isunlimited() else len(dim))
        for name, var in grid.variables.items():
            # A faithful copy: raw values, char arrays as characters, scalars as scalars
            var.set_auto_maskandscale(False)
            var.set_auto_chartostring(False)
            fill = var.getncattr("_FillValue") if "_FillValue" in var.ncattrs() else None
            o = rot.createVariable(name, var.dtype, var.dimensions, fill_value=fill)
            o.set_auto_maskandscale(False)
            o.set_auto_chartostring(False)
            for a in var.ncattrs():
                if a != "_FillValue":
                    o.setncattr(a, var.getncattr(a))
            if name == "angle":
                o[:] = angle
            elif var.ndim == 0:
                o.assignValue(var.getValue())
            else:
                o[:] = var[:]
        if "angle" not in grid.variables:
            o = rot.createVariable("angle", "f8", ("eta_rho", "xi_rho"))
            o[:] = angle
        rot.close()
        print(f"wrote {args.rotated_grid}: angle {angle.min():g} to {angle.max():g} rad")


def check(args):
    h = netCDF4.Dataset(args.history)
    lon = np.asarray(h["lon_rho"][:])
    lat = np.asarray(h["lat_rho"][:])
    a = np.asarray(h["angle"][:]) if "angle" in h.variables else np.zeros_like(lon)
    L = continuous_lon(lon, float(np.mean(lon)))
    f = fields(L, lat)
    c, s = np.cos(a), np.sin(a)
    rec = args.record
    worst = 0.0
    lines = []

    def report(name, err, scale):
        nonlocal worst
        rel = err / scale
        worst = max(worst, rel)
        lines.append(f"{name:6s} max |REMORA - exact| = {err:.3e}  (relative {rel:.1e})")

    # Interior only: the ring is the boundary fill's, as for a file on the model grid, and
    # record 0 is written before any forcing is read.
    if "Tair" in h.variables:
        got = np.asarray(h["Tair"][rec])
        report("Tair", np.abs(got - f["Tair"])[1:-1, 1:-1].max(), np.abs(f["Tair"]).max())
    if "Uwind" in h.variables:
        xi, eta = f["Uwind"] * c + f["Vwind"] * s, -f["Uwind"] * s + f["Vwind"] * c
        for name, exp in (("Uwind", xi), ("Vwind", eta)):
            got = np.asarray(h[name][rec])
            report(name, np.abs(got - exp)[1:-1, 1:-1].max(), np.abs(exp).max())
    if args.stress:
        tx, ty = f["sustr"] * c + f["svstr"] * s, -f["sustr"] * s + f["svstr"] * c
        su_exp = 0.5 * (tx[:, :-1] + tx[:, 1:])
        sv_exp = 0.5 * (ty[:-1, :] + ty[1:, :])
        for name, exp in (("sustr", su_exp), ("svstr", sv_exp)):
            got = np.asarray(h[name][rec])
            report(name, np.abs(got - exp)[1:-1, 1:-1].max(), np.abs(exp).max())

    if not lines:
        sys.exit(f"{args.history} holds none of Tair, Uwind/Vwind (set remora.output_forcing = true)")
    print("\n".join(lines))
    # Roundoff of the precision REMORA ran in, which the history file is written in
    tol = args.tol if args.tol is not None else 1000 * float(np.finfo(h["lon_rho"].dtype).eps)
    ok = worst <= tol
    print(("PASS" if ok else "FAIL") + f": worst relative error {worst:.1e} (tolerance {tol:.1e})")
    sys.exit(0 if ok else 1)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="command", required=True)

    pm = sub.add_parser("make", help="write the forcing file (and a rotated grid)")
    pm.add_argument("--grid", required=True, help="model grid file with lon_rho and lat_rho")
    pm.add_argument("--out", required=True, help="forcing file to write")
    pm.add_argument("--spacing", type=float, default=2.5, help="source grid spacing, degrees")
    pm.add_argument("--margin", type=float, default=1.5, help="extent beyond the model grid, degrees")
    pm.add_argument("--lon360", action="store_true", help="store longitude in [0, 360)")
    pm.add_argument("--flip-lat", action="store_true", help="store latitude decreasing, as ERA5 does")
    pm.add_argument("--rotated-grid", help="also write a copy of the grid with angle = 0.1 + 0.01 i + 0.02 j")
    pm.set_defaults(func=make)

    pc = sub.add_parser("check", help="compare a history file with the exact fields")
    pc.add_argument("history", help="REMORA NetCDF history file (remora.output_forcing = true)")
    pc.add_argument("--record", type=int, default=-1, help="history record to check (default last)")
    pc.add_argument("--stress", action="store_true", help="also check sustr and svstr")
    pc.add_argument("--tol", type=float, default=None,
                    help="largest relative error to accept (default 1000 machine epsilons of the "
                         "history file's precision: 2e-13 in double, 1e-4 in single)")
    pc.set_defaults(func=check)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
