#!/usr/bin/env python3
"""Generate small CDF-2 tidal and subtidal fixtures; requires numpy and netCDF4."""

from pathlib import Path
import sys

import netCDF4
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "Tools"))
from make_boundary_tides import write_boundary_tides


def main():
    directory = Path(__file__).resolve().parents[1] / "test_files" / "TidalBoundary"
    directory.mkdir(exist_ok=True)
    nx, ny = 8, 10
    epoch = "2000-01-02 00:00:00"
    periods = np.asarray([240.0, 600.0])
    sides = ("west", "south", "east", "north")
    scales = {"zeta": 0.1, "ubar": 0.01, "vbar": 0.02}
    fields = {}
    for b, side in enumerate(sides):
        length = ny + 2 if side in ("west", "east") else nx + 2
        p = np.arange(length)
        normal = "ubar" if side in ("west", "east") else "vbar"
        for quantity in ("zeta", normal):
            amplitude = scales[quantity] * (1 + 0.01 * p) * np.asarray([[1], [0.3]])
            phase = np.asarray([np.where(p % 2 == 0, 359.0, 1.0),
                                np.full(length, 120.0 + 10.0 * b)])
            fields[f"{quantity}_{side}"] = (amplitude, phase)
    write_boundary_tides(directory / "tides.nc", nx, ny, ["test240", "test600"],
                         periods, epoch, fields)
    with netCDF4.Dataset(directory / "tides.nc") as source:
        with netCDF4.Dataset(directory / "bad_phase.nc", "w",
                             format="NETCDF3_64BIT_OFFSET") as target:
            for name, dim in source.dimensions.items():
                target.createDimension(name, len(dim))
            for name, var in source.variables.items():
                copy = target.createVariable(name, var.dtype, var.dimensions)
                copy.setncatts({key: var.getncattr(key) for key in var.ncattrs()})
                copy[...] = var[...]
            target["tide_zeta_phase_west"].units = "radians"

    for filename, times, include_tides in (
        ("subtidal.nc", np.asarray([0.0, 86400.0]), False),
        ("reference.nc", np.arange(0.0, 301.0, 30.0), True),
    ):
        with netCDF4.Dataset(directory / filename, "w", format="NETCDF3_64BIT_OFFSET") as nc:
            for dim, length in {"xi_rho": nx+2, "eta_rho": ny+2,
                                "xi_u": nx+1, "eta_v": ny+1,
                                "ocean_time": len(times)}.items():
                nc.createDimension(dim, length)
            time = nc.createVariable("ocean_time", "f8", ("ocean_time",))
            time.units = f"days since {epoch}"
            time[:] = times / 86400.0
            for quantity in scales:
                for b, side in enumerate(sides):
                    dim = (("eta_v" if quantity == "vbar" else "eta_rho")
                           if side in ("west", "east") else
                           ("xi_u" if quantity == "ubar" else "xi_rho"))
                    p = np.arange(len(nc.dimensions[dim]))
                    values = (0.02 + 0.002*b + 0.0001*p
                              + 0.003*times[:, None]/86400.0)
                    if include_tides and f"{quantity}_{side}" in fields:
                        amplitude, phase = fields[f"{quantity}_{side}"]
                        values = values + np.sum(
                            amplitude[None, :, :] * np.cos(
                                2*np.pi*times[:, None, None]/periods[None, :, None]
                                - np.deg2rad(phase)[None, :, :]), axis=1)
                    var = nc.createVariable(f"{quantity}_{side}", "f8", ("ocean_time", dim))
                    var[:] = values


if __name__ == "__main__":
    main()
