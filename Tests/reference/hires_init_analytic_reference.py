#!/usr/bin/env python3
"""Independent reference for analytic high-resolution initialization (Seamount).

With remora.hires_init_level = L, level-0 temperature must be the r x r block mean of the
analytic profile at level-L depths. This rebuilds that field in NumPy, sharing no code with
REMORA (Seamount bathymetry, Vtransform 2 with calc_stretch_coeffs' stretching,
T = T0 + 7.5 exp(z/1000)), and compares it cell by cell with level 0 of a plotfile. It also
prints the gap to the plain level-0 evaluation, which is what a run ignoring the flag gives.

Usage
-----
    python3 Tests/reference/hires_init_analytic_reference.py plt00000 --ref-ratio 3

--ref-ratio is the cumulative x/y ratio from level 0 to hires_init_level. The bathymetry is
evaluated analytically there, as the code does unless hires_grid_level is higher; then pass
--grid-ratio, the cumulative ratio to hires_grid_level, and it is block-averaged down first.
"""

import argparse
import os
import re

import numpy as np


def read_plotfile_level0(plt, var):
    """Return (field[z, y, x], prob_lo, prob_hi, ncell) for one variable on level 0."""
    with open(os.path.join(plt, "Header")) as f:
        lines = f.read().split("\n")
    nvar = int(lines[1])
    names = lines[2:2 + nvar]
    icomp = names.index(var)
    i = 2 + nvar
    i += 3                              # dim, time, finest level
    prob_lo = np.array(lines[i].split(), dtype=float)
    prob_hi = np.array(lines[i + 1].split(), dtype=float)
    i += 3                              # prob_lo, prob_hi, ref ratios (empty for one level)
    nums = [int(v) for v in re.findall(r"-?\d+", lines[i])]
    lo, hi = np.array(nums[0:3]), np.array(nums[3:6])
    ncell = hi - lo + 1

    field = np.full(ncell[::-1], np.nan)
    with open(os.path.join(plt, "Level_0", "Cell_H")) as f:
        cell_h = f.read()
    boxes = [tuple(int(v) for v in re.findall(r"-?\d+", b))
             for b in re.findall(r"\(\([^)]*\) \([^)]*\) \([^)]*\)\)", cell_h)]
    fods = re.findall(r"FabOnDisk: (\S+) (\d+)", cell_h)
    for box, (fname, offset) in zip(boxes, fods):
        blo, bhi = np.array(box[0:3]), np.array(box[3:6])
        with open(os.path.join(plt, "Level_0", fname), "rb") as f:
            f.seek(int(offset))
            header = f.readline().decode()
            assert "(8 7 6 5 4 3 2 1)" in header, "expected little-endian doubles: " + header
            ncomp = int(header.split()[-1])
            shape = bhi - blo + 1
            data = np.fromfile(f, dtype="<f8", count=int(np.prod(shape)) * ncomp)
        data = data.reshape(ncomp, shape[2], shape[1], shape[0])
        field[blo[2]:bhi[2] + 1, blo[1]:bhi[1] + 1, blo[0]:bhi[0] + 1] = data[icomp]
    assert not np.isnan(field).any(), "level 0 not fully covered by the boxes read"
    return field, prob_lo, prob_hi, ncell


def stretching(N, theta_s, theta_b):
    """Cs at rho points, as calc_stretch_coeffs computes it."""
    k = np.arange(N)
    s_r = (k - N + 0.5) / N
    if theta_s > 0:
        csur = (1.0 - np.cosh(theta_s * s_r)) / (np.cosh(theta_s) - 1.0)
    else:
        csur = -s_r**2
    if theta_b > 0:
        return s_r, (np.exp(theta_b * csur) - 1.0) / (1.0 - np.exp(-theta_b))
    return s_r, csur


def seamount_h(nx, ny, prob_lo, prob_hi):
    """Seamount bathymetry on an nx x ny grid, as [y, x]."""
    dx = (prob_hi[0] - prob_lo[0]) / nx
    dy = (prob_hi[1] - prob_lo[1]) / ny
    x = prob_lo[0] + (np.arange(nx) + 0.5) * dx
    y = prob_lo[1] + (np.arange(ny) + 0.5) * dy
    X, Y = np.meshgrid(x, y)
    return 5000.0 - 4500.0 * np.exp(-(((X - 160000.0) / 40000.0)**2 +
                                      ((Y - 160000.0) / 40000.0)**2))


def seamount_temp(nx, ny, N, prob_lo, prob_hi, args, h_ratio=1):
    """Seamount temperature on an nx x ny x N grid, as [z, y, x]. The bathymetry is evaluated
    h_ratio times finer and block-averaged onto the grid."""
    h = seamount_h(nx * h_ratio, ny * h_ratio, prob_lo, prob_hi)
    h = h.reshape(ny, h_ratio, nx, h_ratio).mean(axis=(1, 3))
    hc = -min(prob_hi[2], -args.tcline)
    s_r, cs_r = stretching(N, args.theta_s, args.theta_b)
    zeta = 0.0
    z_r = zeta + (zeta + h[None]) * (hc * s_r[:, None, None] + cs_r[:, None, None] * h[None]) \
        / (hc + h[None])
    return args.T0 + 7.5 * np.exp(z_r / 1000.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfile")
    ap.add_argument("--ref-ratio", type=int, required=True)
    ap.add_argument("--grid-ratio", type=int, default=None)
    ap.add_argument("--theta-s", type=float, default=6.5)
    ap.add_argument("--theta-b", type=float, default=2.0)
    ap.add_argument("--tcline", type=float, default=100.0)
    ap.add_argument("--T0", type=float, default=10.0)
    args = ap.parse_args()

    temp, prob_lo, prob_hi, ncell = read_plotfile_level0(args.plotfile, "temp")
    nx, ny, N = ncell
    r = args.ref_ratio

    g = args.grid_ratio if args.grid_ratio else r
    if g % r:
        raise SystemExit("--grid-ratio must be a multiple of --ref-ratio")
    fine = seamount_temp(nx * r, ny * r, N, prob_lo, prob_hi, args, h_ratio=g // r)
    expected = fine.reshape(N, ny, r, nx, r).mean(axis=(2, 4))
    plain = seamount_temp(nx, ny, N, prob_lo, prob_hi, args)

    err = np.abs(temp - expected).max()
    gap = np.abs(temp - plain).max()
    print("max |plotfile - hires reference|     = %.3e" % err)
    print("max |plotfile - level-0 evaluation|  = %.3e" % gap)
    print("reference extrema: min %.12f max %.12f" % (expected.min(), expected.max()))


if __name__ == "__main__":
    main()
