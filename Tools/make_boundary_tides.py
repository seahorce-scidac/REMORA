#!/usr/bin/env python3
"""Write REMORA boundary tides from prepared, model-aligned harmonic arrays.

Import write_boundary_tides to supply real boundary harmonics. Running this file writes
a synthetic M2/S2 example, not an atlas-derived forcing file. Requires numpy and netCDF4.
"""

import argparse

import netCDF4
import numpy as np


def write_boundary_tides(path, nx, ny, names, periods, epoch, fields):
    """Write CDF-2 amplitude/phase pairs indexed by (constituent, boundary point).

    fields maps e.g. 'zeta_west' or 'ubar_east' to (amplitude, phase_degrees).
    Elevation is in metres, velocity in metres per second, and periods in seconds.
    epoch is a UTC date string defining zero elapsed time for the supplied phase lags.
    Boundary arrays include both endpoint entries: nx+2 or ny+2 values.
    """
    with netCDF4.Dataset(path, "w", format="NETCDF3_64BIT_OFFSET") as nc:
        nc.createDimension("tide", len(names))
        nc.createDimension("name_strlen", max(len(name) for name in names))
        nc.createDimension("xi_rho", nx + 2)
        nc.createDimension("eta_rho", ny + 2)
        labels = nc.createVariable("tide_name", "S1", ("tide", "name_strlen"))
        labels[:] = np.asarray([list(name.ljust(labels.shape[1])) for name in names], "S1")
        period = nc.createVariable("tide_period", "f8", ("tide",))
        period.units = "s"
        period[:] = periods
        origin = nc.createVariable("tide_epoch", "f8")
        origin.units = f"days since {epoch}"
        origin.calendar = "proleptic_gregorian"
        origin.assignValue(0.0)
        for field, (amplitude, phase) in fields.items():
            quantity, side = field.split("_")
            if side not in ("west", "east", "south", "north"):
                raise ValueError(f"{field}: unknown boundary side")
            dimension = "eta_rho" if side in ("west", "east") else "xi_rho"
            if quantity not in ("zeta", "ubar" if dimension == "eta_rho" else "vbar"):
                raise ValueError(f"{field}: supply elevation or the normal velocity component")
            for suffix, values, units in (
                ("amp", amplitude, "m" if quantity == "zeta" else "m s-1"),
                ("phase", phase, "degrees"),
            ):
                var = nc.createVariable(f"tide_{quantity}_{suffix}_{side}", "f8",
                                        ("tide", dimension))
                var.units = units
                var[:] = values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, required=True, help="interior level-0 x cells")
    parser.add_argument("--ny", type=int, required=True, help="interior level-0 y cells")
    parser.add_argument("--epoch", default="2020-01-01 00:00:00")
    parser.add_argument("--output", default="boundary_tides.nc")
    args = parser.parse_args()
    fields = {}
    for side in ("west", "east", "south", "north"):
        length = args.ny + 2 if side in ("west", "east") else args.nx + 2
        shape = (2, length)
        fields[f"zeta_{side}"] = (np.broadcast_to([[0.1], [0.03]], shape),
                                    np.broadcast_to([[0.0], [45.0]], shape))
        velocity = "ubar" if side in ("west", "east") else "vbar"
        fields[f"{velocity}_{side}"] = (np.broadcast_to([[0.01], [0.003]], shape),
                                       np.broadcast_to([[30.0], [90.0]], shape))
    write_boundary_tides(args.output, args.nx, args.ny, ["M2", "S2"],
                         [44714.164, 43200.0], args.epoch, fields)
    print(f"Wrote synthetic example {args.output}; replace its harmonics for a real simulation.")


if __name__ == "__main__":
    main()
