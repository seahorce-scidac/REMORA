#!/usr/bin/env python3
"""Independent reference for the high-resolution average-down.

Computes what ``remora.sum_interval``'s ``VOLUME`` diagnostic must report at t = 0 for a run with
``hires_grid_level`` / ``hires_init_level`` set, straight from the NetCDF input files -- without
running REMORA. Comparing the two checks the average-down arithmetic against an implementation
that shares no code with it.

Why VOLUME is the right quantity: at t = 0, ``sum_k Hz(i,j,k) == zeta(i,j) + h(i,j)`` exactly for
any vertical stretching, because ``Cs_w(0) == sc_w(0) == -1`` gives ``z_w(i,j,0) == -h``
(Source/Utils/REMORA_DepthStretchTransform.H). So

    VOLUME(t=0) == sum_ij (zeta + h) / (pm * pn)

is a direct readout of the averaged-down bathymetry, free surface and grid metrics.

What this pins down, all of which a "the run completed" test would miss:
  * ``h`` and ``zeta`` are averaged down as a plain arithmetic mean over r x r blocks
  * ``pm`` and ``pn`` are averaged down and then divided by the refinement ratio -- coarsening
    multiplies dx by r, so pm = 1/dx must shrink by r
  * the block alignment: which fine cells belong to which coarse cell, and that the grow rings
    are excluded from the interior sum

Usage
-----
    python3 Tests/reference/hires_volume_reference.py --grid hires_grd.nc --init hires_ini.nc \\
            --ref-ratio 3

then compare against a run with ``remora.v = 1 remora.sum_interval = 1
remora.sum_precision = 15`` (the default 6 digits is too coarse: the hires-vs-non-hires signal
can be ~1e-5 relative).

Verified against Exec/GulfRefinementTest/inputs_synth_hires on the files
``make_hires_test_data.py --n-cell 10 10 4 --ref-ratio 3`` produces: expected and reported
VOLUME agree to 4e-16 relative.

Note this is an exact identity only because it compares like with like. Summing
``h/(pm*pn)`` over the *fine* grid instead is NOT the same number when pm varies in space --
averaging a product is not the product of averages -- so a small mismatch there is quadrature,
not a bug.
"""

import argparse

import netCDF4
import numpy as np


def coarsen(field, r, interior):
    """Arithmetic mean over r x r blocks of the interior of ``field``."""
    a = field[interior]
    n0, n1 = a.shape[0] // r, a.shape[1] // r
    if n0 * r != a.shape[0] or n1 * r != a.shape[1]:
        raise SystemExit("interior %s is not divisible by the refinement ratio %d"
                         % (a.shape, r))
    return a.reshape(n0, r, n1, r).mean(axis=(1, 3))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grid", required=True, help="high-resolution grid file (h, pm, pn)")
    ap.add_argument("--init", help="high-resolution initial file (zeta); zeta = 0 if omitted")
    ap.add_argument("--ref-ratio", type=int, required=True,
                    help="cumulative refinement ratio at the hires level")
    args = ap.parse_args()

    r = args.ref_ratio
    grid = netCDF4.Dataset(args.grid)
    h = np.array(grid.variables["h"][:])
    pm = np.array(grid.variables["pm"][:])
    pn = np.array(grid.variables["pn"][:])

    # The file carries the refined domain plus r grow rings on every side; the sum is over the
    # interior only, since that is what maps onto the level-0 domain.
    interior = (slice(r, -r), slice(r, -r))

    if args.init:
        init = netCDF4.Dataset(args.init)
        zeta = np.array(init.variables["zeta"][:])
        zeta = zeta[0] if zeta.ndim == 3 else zeta
    else:
        zeta = np.zeros_like(h)

    h_c = coarsen(h, r, interior)
    zeta_c = coarsen(zeta, r, interior)
    pm_c = coarsen(pm, r, interior) / r
    pn_c = coarsen(pn, r, interior) / r

    volume = ((zeta_c + h_c) / (pm_c * pn_c)).sum()
    area = (1.0 / (pm_c * pn_c)).sum()

    print("hires interior   : %d x %d cells" % (h[interior].shape[1], h[interior].shape[0]))
    print("level-0 equivalent: %d x %d cells" % (h_c.shape[1], h_c.shape[0]))
    print("expected AREA     = %.15g" % area)
    print("expected VOLUME   = %.15g" % volume)
    print()
    print("compare with:  remora.v=1 remora.sum_interval=1 remora.sum_precision=15")


if __name__ == "__main__":
    main()
