#!/usr/bin/env python3
"""Check the physical z coordinates a viewer reconstructs from a REMORA AMReX plotfile.

usage: check_mapped_grid.py PLOTFILE [--periodic x y] [--cf-tol TOL] [--tol 1e-9]

A viewer (amrexplorer, for one) places node (i,j,k) of level L at

    prob_lo + (i*dx_L, j*dy_L, k*dz_L) + nu(i,j,k)

using that level's cell-size line in the Header and the level's Nu_nd set (amrexvec_nu_x/y/z),
and ignores the Header's refinement-ratio line. This script rebuilds z the same way, level by
level, and checks:

  layout     every Nu_nd box is the matching Cell box grown by one on the high side (nodal), and
             the nodal k extent is the Header's domain k extent plus one;
  bottom/top the lowest node layer is the psi-point average of -h and the highest is that of
             zeta, where the plotfile carries h and zeta (remora.plot_vars_2d = "h zeta");
  monotone   z increases with k in every column;
  coarse-fine at every fine node that coincides with a coarse node, |z_fine - z_coarse|, split
             into the patch perimeter and its interior. Reported always; asserted only with
             --cf-tol.

On a refined level the bottom/top checks are asserted at interior nodes only -- a perimeter node
averages fine ghost columns that are not in the plotfile -- and the perimeter is reported.
Columns next to a masked cell (mask_rho = 0) are skipped for the top check, since the plotfile's
zeta is a fill value there. Exit status: 0 pass, 1 fail, 2 unusable input.
"""
import os
import re
import sys

import numpy as np

BOX_RE = re.compile(r"\(\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\) \((\d+),(\d+),(\d+)\)\)")


def parse_box(s):
    g = BOX_RE.search(s).groups()
    return (tuple(int(x) for x in g[0:3]), tuple(int(x) for x in g[3:6]), tuple(int(x) for x in g[6:9]))


def read_header(plt):
    """Header -> dict: names, finest, prob_lo, dx[lev], domain[lev], sets {setname: (names, paths)}."""
    with open(os.path.join(plt, "Header")) as f:
        lines = [ln.rstrip("\n") for ln in f]
    p = 0
    p += 1                                            # version
    ncomp = int(lines[p]); p += 1
    names = lines[p:p + ncomp]; p += ncomp
    dim = int(lines[p]); p += 1
    if dim != 3:
        sys.exit("only 3-d plotfiles are supported")
    p += 1                                            # time
    finest = int(lines[p]); p += 1
    prob_lo = [float(x) for x in lines[p].split()]; p += 1
    p += 1                                            # prob_hi
    p += 1                                            # refinement ratios: read and ignored, as a viewer does
    domains = [parse_box(s) for s in re.findall(r"\(\([^)]*\) \([^)]*\) \([^)]*\)\)", lines[p])]; p += 1
    p += 1                                            # level steps
    dx = []
    for _ in range(finest + 1):
        dx.append([float(x) for x in lines[p].split()]); p += 1
    p += 2                                            # coord sys, "0"
    cell_paths = []
    for lev in range(finest + 1):
        nbox = int(lines[p].split()[1]); p += 1
        p += 1                                        # level steps
        p += 3 * nbox                                 # physical box bounds
        cell_paths.append(lines[p]); p += 1
    sets = {"Cell": (names, cell_paths)}
    if p < len(lines) and lines[p].strip():
        p += 1                                        # num_extra_mfs: an undercount, so read to the end
    while p < len(lines) and lines[p].strip():
        n = int(lines[p]); p += 1
        snames = lines[p:p + n]; p += n
        paths = lines[p:p + finest + 1]; p += finest + 1
        key = os.path.basename(paths[0])              # Nu_nd, UFace, rho2d, ...
        sets[key] = (snames, paths)
    return dict(names=names, finest=finest, prob_lo=prob_lo, dx=dx, domains=domains, sets=sets)


def read_mf_header(plt, path):
    """<path>_H -> (boxes [(lo,hi,type)], fabs [(file, offset)], ncomp, nghost)."""
    with open(os.path.join(plt, path + "_H")) as f:
        lines = f.read().splitlines()
    ncomp = int(lines[2])
    nghost = int(lines[3].split()[0].strip("()"))
    nbox = int(re.match(r"\((\d+) ", lines[4]).group(1))
    boxes = [parse_box(lines[5 + k]) for k in range(nbox)]
    idx = 5 + nbox + 1
    assert int(lines[idx]) == nbox, lines[idx]
    fabs = []
    for k in range(nbox):
        parts = lines[idx + 1 + k].split()
        fabs.append((parts[1], int(parts[2])))
    return boxes, fabs, ncomp, nghost


def read_set(plt, path, icomp):
    """One component of one level's MultiFab, assembled on its bounding box; NaN where no box covers.
    Returns (array, lo_of_bounding_box, boxes)."""
    boxes, fabs, ncomp, nghost = read_mf_header(plt, path)
    level_dir = os.path.dirname(os.path.join(plt, path))
    glo = tuple(min(b[0][a] for b in boxes) for a in range(3))
    ghi = tuple(max(b[1][a] for b in boxes) for a in range(3))
    arr = np.full(tuple(ghi[a] - glo[a] + 1 for a in range(3)), np.nan)
    for (lo, hi, _), (fname, off) in zip(boxes, fabs):
        with open(os.path.join(level_dir, fname), "rb") as f:
            f.seek(off)
            hdr = b""
            while not hdr.endswith(b"\n"):
                hdr += f.read(1)
            flo, fhi, _ = parse_box(hdr.decode())     # the stored box, ghost-grown if nghost > 0
            n = tuple(fhi[a] - flo[a] + 1 for a in range(3))
            npts = n[0] * n[1] * n[2]
            f.seek(off + len(hdr) + icomp * npts * 8)
            data = np.frombuffer(f.read(npts * 8), dtype="<f8").reshape(n[2], n[1], n[0]).transpose(2, 1, 0)
        sl_src = tuple(slice(lo[a] - flo[a], hi[a] - flo[a] + 1) for a in range(3))
        sl_dst = tuple(slice(lo[a] - glo[a], hi[a] - glo[a] + 1) for a in range(3))
        arr[sl_dst] = data[sl_src]
    return arr, glo, boxes


def psi_average(cell, periodic, dom_lo, dom_hi):
    """Node (i,j) <- mean of cells (i-1..i, j-1..j); clamp at non-periodic edges, wrap on periodic
    ones. `cell` is indexed from dom_lo. NaN propagates (an uncovered cell poisons its nodes)."""
    nx, ny = cell.shape
    ii = np.arange(nx + 1); jj = np.arange(ny + 1)
    def pick(idx, n, per):
        a = idx - 1; b = idx
        if per:
            a = a % n; b = b % n
        else:
            a = np.clip(a, 0, n - 1); b = np.clip(b, 0, n - 1)
        return a, b
    ia, ib = pick(ii, nx, periodic[0]); ja, jb = pick(jj, ny, periodic[1])
    return 0.25 * (cell[np.ix_(ia, ja)] + cell[np.ix_(ib, ja)] + cell[np.ix_(ia, jb)] + cell[np.ix_(ib, jb)])


def main():
    args = [a for a in sys.argv[1:]]
    if not args:
        print(__doc__); sys.exit(2)
    plt = args[0]
    periodic = [False, False]
    cf_tol = None; tol = 1e-9
    i = 1
    while i < len(args):
        if args[i] == "--periodic":
            i += 1
            while i < len(args) and args[i] in ("x", "y"):
                periodic["xy".index(args[i])] = True; i += 1
            continue
        if args[i] == "--cf-tol": cf_tol = float(args[i + 1]); i += 2; continue
        if args[i] == "--tol": tol = float(args[i + 1]); i += 2; continue
        sys.exit(f"unknown argument {args[i]}")

    H = read_header(plt)
    finest = H["finest"]
    if "Nu_nd" not in H["sets"]:
        print("no Nu_nd set in the Header: nothing to check (remora.plot_nodal_data = 0?)"); sys.exit(2)
    nd_names, nd_paths = H["sets"]["Nu_nd"]
    if nd_names != ["amrexvec_nu_x", "amrexvec_nu_y", "amrexvec_nu_z"]:
        print(f"unexpected nodal component names {nd_names}"); sys.exit(2)
    rho2d = H["sets"].get("rho2d")
    have_hz = rho2d is not None and "h" in rho2d[0] and "zeta" in rho2d[0]
    have_mask = rho2d is not None and "mask_rho" in rho2d[0]
    if not have_hz:
        print("note: no h/zeta in the plotfile (set remora.plot_vars_2d = \"h zeta\"); bottom/top checks skipped")

    failed = False
    z_levels = []        # (z array on nodal bounding box, glo, dom_lo, dx)
    cell_cover = []      # (bool array of covered cells on the Cell bounding box, glo)
    print(f"{plt}: {finest + 1} level(s); periodic x={periodic[0]} y={periodic[1]}")
    for lev in range(finest + 1):
        dxl = H["dx"][lev]
        dom_lo, dom_hi, _ = H["domains"][lev]
        cell_boxes, _, _, _ = read_mf_header(plt, H["sets"]["Cell"][1][lev])
        nd_boxes, _, _, _ = read_mf_header(plt, nd_paths[lev])
        # ---- layout
        ok = len(cell_boxes) == len(nd_boxes)
        if ok:
            for (clo, chi, _), (nlo, nhi, nty) in zip(cell_boxes, nd_boxes):
                if nty != (1, 1, 1) or nlo != clo or tuple(h + 1 for h in chi) != nhi:
                    ok = False; break
        nd_k = max(b[1][2] for b in nd_boxes) - min(b[0][2] for b in nd_boxes) + 1
        dom_k = dom_hi[2] - dom_lo[2] + 1
        layout_msg = (f"Nu_nd boxes {len(nd_boxes)} vs Cell {len(cell_boxes)}; nodal k extent {nd_k}, "
                      f"Header domain k extent {dom_k} (+1 expected)")
        if not ok or nd_k != dom_k + 1:
            failed = True
            print(f"  level {lev}: layout FAIL -- {layout_msg}")
        else:
            print(f"  level {lev}: layout ok -- {layout_msg}")
        # ---- reconstruct z exactly as the viewer does
        nu_z, glo, _ = read_set(plt, nd_paths[lev], 2)
        k_idx = glo[2] + np.arange(nu_z.shape[2])
        z = H["prob_lo"][2] + (k_idx - dom_lo[2])[None, None, :] * dxl[2] + nu_z
        z_levels.append((z, glo, dom_lo, dxl))
        cells, cglo, _ = read_set(plt, H["sets"]["Cell"][1][lev], 0)
        cell_cover.append((~np.isnan(cells[:, :, 0]), cglo))
        # ---- monotone in k, over covered columns
        col_ok = ~np.isnan(z).any(axis=2)
        dz_k = np.diff(z, axis=2)
        bad = (dz_k[col_ok] <= 0).sum()
        if bad:
            failed = True
            print(f"  level {lev}: monotone FAIL -- {bad} non-increasing node intervals")
        else:
            print(f"  level {lev}: monotone ok -- min layer thickness {np.nanmin(dz_k[col_ok]):.4g}")
        # ---- bottom and top against h and zeta
        if have_hz:
            names2d, paths2d = rho2d
            h2d, hlo, _ = read_set(plt, paths2d[lev], names2d.index("h"))
            zeta2d, _, _ = read_set(plt, paths2d[lev], names2d.index("zeta"))
            h2d = h2d[:, :, 0]; zeta2d = zeta2d[:, :, 0]
            assert hlo[0] == glo[0] and hlo[1] == glo[1], (hlo, glo)
            scale = np.nanmax(np.abs(h2d))
            # A node on the edge of the bounding box averages columns the plotfile does not carry:
            # ghost h/zeta at a physical boundary on level 0, interpolated ghosts on a patch
            # perimeter. Wrap on a periodic axis of level 0 (every column is in the file); elsewhere
            # clamp for the report and leave the edge out of the assertion.
            per = periodic if lev == 0 else [False, False]
            hb = psi_average(-h2d, per, None, None)
            zt = psi_average(zeta2d, per, None, None)
            interior = ~np.isnan(hb)
            if not per[0]: interior[0, :] = False; interior[-1, :] = False
            if not per[1]: interior[:, 0] = False; interior[:, -1] = False
            d_bot = np.abs(z[:, :, 0] - hb); d_top = np.abs(z[:, :, -1] - zt)
            if have_mask:
                m2d, _, _ = read_set(plt, paths2d[lev], names2d.index("mask_rho"))
                wet = psi_average(m2d[:, :, 0], per, None, None) > 0.999
            else:
                wet = np.ones_like(interior)
            bot_int = np.nanmax(d_bot[interior]) if interior.any() else 0.0
            top_int = np.nanmax(d_top[interior & wet]) if (interior & wet).any() else 0.0
            per_mask = ~interior & ~np.isnan(hb)
            bot_per = np.nanmax(d_bot[per_mask]) if per_mask.any() else 0.0
            top_per = np.nanmax(d_top[per_mask & wet]) if (per_mask & wet).any() else 0.0
            verdict = "ok" if (bot_int <= tol * scale and top_int <= tol * scale) else "FAIL"
            if verdict == "FAIL": failed = True
            worst = d_bot if bot_int >= top_int else d_top
            iw, jw = np.unravel_index(np.nanargmax(np.where(interior & wet, worst, -1.0)), worst.shape)
            edge = "perimeter" if lev > 0 else "domain edge"
            print(f"  level {lev}: bottom/top {verdict} -- interior max |z0+h| {bot_int:.3e}, |zN-zeta| {top_int:.3e}"
                  f" at node ({iw + glo[0]},{jw + glo[1]}); {edge} {bot_per:.3e}, {top_per:.3e}  (tol {tol * scale:.1e})")
    # ---- coarse-fine agreement
    for lev in range(1, finest + 1):
        zf, flo, fdom, fdx = z_levels[lev]; zc, clo, cdom, cdx = z_levels[lev - 1]
        r = [int(round(cdx[a] / fdx[a])) for a in range(3)]
        fi = flo[0] + np.arange(zf.shape[0]); fj = flo[1] + np.arange(zf.shape[1]); fk = flo[2] + np.arange(zf.shape[2])
        si = np.where(fi % r[0] == 0)[0]; sj = np.where(fj % r[1] == 0)[0]; sk = np.where(fk % r[2] == 0)[0]
        ci = fi[si] // r[0] - clo[0]; cj = fj[sj] // r[1] - clo[1]; ck = fk[sk] // r[2] - clo[2]
        keep_i = (ci >= 0) & (ci < zc.shape[0]); keep_j = (cj >= 0) & (cj < zc.shape[1]); keep_k = (ck >= 0) & (ck < zc.shape[2])
        si, ci = si[keep_i], ci[keep_i]; sj, cj = sj[keep_j], cj[keep_j]; sk, ck = sk[keep_k], ck[keep_k]
        d = np.abs(zf[np.ix_(si, sj, sk)] - zc[np.ix_(ci, cj, ck)])
        # perimeter: a fine node with an uncovered neighbouring fine cell
        cov, cglo = cell_cover[lev]
        padded = np.zeros((cov.shape[0] + 2, cov.shape[1] + 2), dtype=bool); padded[1:-1, 1:-1] = cov
        # node (i,j) touches cells (i-1..i, j-1..j); in padded coordinates cell c is at c+1
        ni = fi[si] - cglo[0]; nj = fj[sj] - cglo[1]
        touch = (padded[np.ix_(ni, nj)] & padded[np.ix_(ni + 1, nj)] & padded[np.ix_(ni, nj + 1)] & padded[np.ix_(ni + 1, nj + 1)])
        interior = touch[:, :, None] & ~np.isnan(d)
        perim = ~touch[:, :, None] & ~np.isnan(d)
        d_int = np.nanmax(d[interior]) if interior.any() else 0.0
        d_per = np.nanmax(d[perim]) if perim.any() else 0.0
        worst_k = int(fk[sk][np.unravel_index(np.nanargmax(np.where(np.isnan(d), -1, d)), d.shape)[2]])
        verdict = ""
        if cf_tol is not None:
            verdict = "ok" if d_int <= cf_tol else "FAIL"
            if verdict == "FAIL": failed = True
        print(f"  levels {lev - 1}/{lev}: coarse-fine {verdict} ratio {r} -- max |dz| interior {d_int:.3e} "
              f"({int(interior.sum())} nodes), perimeter {d_per:.3e} ({int(perim.sum())} nodes), worst fine k {worst_k}")
    print("RESULT:", "FAIL" if failed else "PASS")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
