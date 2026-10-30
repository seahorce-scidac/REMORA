"""Make a masked NetCDF grid and initial file for land-layout timing (remora.drop_dry_boxes,
remora.load_balance).

    python3 Tests/perf/make_masked_grid.py <REMORA-data>/IdealMiniGrid <out dir> <factor>

Upsamples the masked IdealMiniGrid grid (10 x 16 rho points plus the one-cell ring, a diagonal
coastline with 36 % land) by an integer factor in x and y: every rho field is repeated
factor x factor, pm and pn are multiplied by the factor, so the domain stays 320 km square
and the coastline becomes a staircase that reaches the domain edges. The initial file holds
zero velocity, a 0.5 m free-surface bump and a 3-degree east-west temperature front on the
water, so the run has something to do. Writes grid.nc and ini.nc, which
Tests/perf/inputs/masked_grid_perf.in reads from the run directory (factor 16 there; set
remora.n_cell = 10*factor 16*factor 20 for another).
"""
import os
import sys

import numpy as np
import netCDF4

src, out, f = sys.argv[1], sys.argv[2], int(sys.argv[3])
os.makedirs(out, exist_ok=True)
g = netCDF4.Dataset(src + "/idmini_grd_masked_v1_classic64.nc")
ini = netCDF4.Dataset(src + "/idmini_ini_v1_classic64.nc")


def up_rho(a):
    """(..., eta_rho, xi_rho) with a one-cell ring -> upsampled interior, same ring."""
    inner = np.repeat(np.repeat(a[..., 1:-1, 1:-1], f, axis=-2), f, axis=-1)
    ny, nx = inner.shape[-2:]
    res = np.empty(a.shape[:-2] + (ny + 2, nx + 2), dtype=a.dtype)
    res[..., 1:-1, 1:-1] = inner
    res[..., 0, 1:-1] = np.repeat(a[..., 0, 1:-1], f, axis=-1)
    res[..., -1, 1:-1] = np.repeat(a[..., -1, 1:-1], f, axis=-1)
    res[..., 1:-1, 0] = np.repeat(a[..., 1:-1, 0], f, axis=-1)
    res[..., 1:-1, -1] = np.repeat(a[..., 1:-1, -1], f, axis=-1)
    for jj in (0, -1):
        for ii in (0, -1):
            res[..., jj, ii] = a[..., jj, ii]
    return res


NY = (g.dimensions["eta_rho"].size - 2) * f + 2
NX = (g.dimensions["xi_rho"].size - 2) * f + 2
nz = g.dimensions["s_rho"].size
dx = 1.0 / (np.asarray(g.variables["pm"][:]).mean() * f)
dy = 1.0 / (np.asarray(g.variables["pn"][:]).mean() * f)


def make(path, with_ini):
    d = netCDF4.Dataset(path, "w", format="NETCDF3_64BIT_OFFSET")
    for name, n in (("eta_rho", NY), ("xi_rho", NX), ("eta_u", NY), ("xi_u", NX - 1),
                    ("eta_v", NY - 1), ("xi_v", NX), ("eta_psi", NY - 1), ("xi_psi", NX - 1),
                    ("s_rho", nz), ("s_w", nz + 1)):
        d.createDimension(name, n)
    rho = ("eta_rho", "xi_rho")
    for name in ("pm", "pn", "f", "h", "mask_rho", "angle", "dmde", "dndx"):
        if name in g.variables:
            a = np.asarray(g.variables[name][:])
            d.createVariable(name, "f8", rho)[:] = up_rho(a * f if name in ("pm", "pn") else a)
    jj, ii = np.meshgrid(np.arange(NY), np.arange(NX), indexing="ij")
    xr, yr = (ii - 0.5) * dx, (jj - 0.5) * dy
    d.createVariable("x_rho", "f8", rho)[:] = xr
    d.createVariable("y_rho", "f8", rho)[:] = yr
    d.createVariable("x_u", "f8", ("eta_u", "xi_u"))[:] = 0.5 * (xr[:, 1:] + xr[:, :-1])
    d.createVariable("y_u", "f8", ("eta_u", "xi_u"))[:] = 0.5 * (yr[:, 1:] + yr[:, :-1])
    d.createVariable("x_v", "f8", ("eta_v", "xi_v"))[:] = 0.5 * (xr[1:, :] + xr[:-1, :])
    d.createVariable("y_v", "f8", ("eta_v", "xi_v"))[:] = 0.5 * (yr[1:, :] + yr[:-1, :])
    quad = lambda a: 0.25 * (a[1:, 1:] + a[:-1, 1:] + a[1:, :-1] + a[:-1, :-1])
    d.createVariable("x_psi", "f8", ("eta_psi", "xi_psi"))[:] = quad(xr)
    d.createVariable("y_psi", "f8", ("eta_psi", "xi_psi"))[:] = quad(yr)
    m = np.asarray(d.variables["mask_rho"][:])
    d.createVariable("mask_u", "f8", ("eta_u", "xi_u"))[:] = m[:, 1:] * m[:, :-1]
    d.createVariable("mask_v", "f8", ("eta_v", "xi_v"))[:] = m[1:, :] * m[:-1, :]
    d.createVariable("mask_psi", "f8", ("eta_psi", "xi_psi"))[:] = quad(m) == 1.0
    for name in ("s_rho", "Cs_r"):
        d.createVariable(name, "f8", ("s_rho",))[:] = g.variables[name][:]
    for name in ("s_w", "Cs_w"):
        d.createVariable(name, "f8", ("s_w",))[:] = g.variables[name][:]
    for name in ("theta_s", "theta_b", "Tcline", "hc", "Vtransform", "xl", "el"):
        if name in g.variables:
            d.createVariable(name, "f8", ())[...] = g.variables[name][...]
    if with_ini:
        d.createDimension("ocean_time", None)
        t = d.createVariable("ocean_time", "f8", ("ocean_time",))
        for k in ini.variables["ocean_time"].ncattrs():
            if k != "_FillValue":
                t.setncattr(k, ini.variables["ocean_time"].getncattr(k))
        t[0:1] = ini.variables["ocean_time"][:]
        temp = up_rho(np.asarray(ini.variables["temp"][:]))
        temp[0] = 20.0 + 3.0 * (ii / NX)[None, :, :]
        d.createVariable("temp", "f8", ("ocean_time", "s_rho") + rho)[:] = temp
        d.createVariable("salt", "f8", ("ocean_time", "s_rho") + rho)[:] = up_rho(
            np.asarray(ini.variables["salt"][:]))
        bump = 0.5 * np.exp(-(((ii - 0.7 * NX) / (0.1 * NX)) ** 2
                              + ((jj - 0.4 * NY) / (0.1 * NY)) ** 2))
        d.createVariable("zeta", "f8", ("ocean_time",) + rho)[:] = (bump * m)[None]
        d.createVariable("u", "f8", ("ocean_time", "s_rho", "eta_u", "xi_u"))[:] = 0.0
        d.createVariable("v", "f8", ("ocean_time", "s_rho", "eta_v", "xi_v"))[:] = 0.0
        d.createVariable("ubar", "f8", ("ocean_time", "eta_u", "xi_u"))[:] = 0.0
        d.createVariable("vbar", "f8", ("ocean_time", "eta_v", "xi_v"))[:] = 0.0
    d.close()


make(out + "/grid.nc", False)
make(out + "/ini.nc", True)
m = np.asarray(netCDF4.Dataset(out + "/grid.nc")["mask_rho"][:])
print(f"{out}/grid.nc: {NX - 2} x {NY - 2} rho points, wet fraction {m[1:-1, 1:-1].mean():.3f}")
