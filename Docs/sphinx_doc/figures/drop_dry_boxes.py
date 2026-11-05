"""Draw drop_dry_boxes.svg: which level-0 boxes remora.drop_dry_boxes leaves out on an example
grid, with the ghost footprint of two boxes and the refinement footprint of one indicator.

The rules follow drop_dry_level0_boxes and level0_refine_footprint in
Source/Initialization/REMORA_land_layout.cpp. Run with python3; it writes the SVG next to itself.
"""
import math
import os

NX, NY, S, OX, OY = 96, 64, 10, 40, 20          # cells, px per cell, plot origin
NGROW = 2
GHOST = NGROW + 1                                 # ghost rings the mask carries
B, NEB, NPROPER, BF1, RR = 8, 1, 1, 8, 2         # box size and AMR settings shown
IND = ((52, 26), (60, 36))                        # refinement indicator, level-0 cells
LEGEND_W = 300

C = dict(bg="#ffffff", fg="#15202a", muted="#5a6875", water="#cfe2ec", land="#d3c4a3",
         kept="#1f5f7a", accent="#7b3fb3", accent_soft="#7b3fb3", drop="#8a3b2e")


def is_land(i, j):
    coast = 63 + 6 * math.sin(j / 7.5) + 3.5 * math.sin(j / 3.1 + 1.0)
    if i >= coast:
        return True
    if ((i - 22) / 7.5) ** 2 + ((j - 13) / 5.5) ** 2 < 1:
        return True
    return ((i - 30) / 4) ** 2 + ((j - 49) / 3) ** 2 < 1


def land_at(i, j):   # the ring beyond the domain copies the edge cell
    return is_land(min(NX - 1, max(0, i)), min(NY - 1, max(0, j)))


def X(i):
    return OX + i * S


def Y(j):
    return OY + (NY - j) * S


def rect(lo, hi, style):
    return (f'<rect x="{X(lo[0])}" y="{Y(hi[1] + 1)}" width="{(hi[0] - lo[0] + 1) * S}" '
            f'height="{(hi[1] - lo[1] + 1) * S}" {style}/>')


def footprint():
    bf = max(1, BF1 // RR)
    stencil = (NGROW + 2 + RR - 1) // RR + 1
    lo, hi = list(IND[0]), list(IND[1])
    layers = [(lo[:], hi[:])]
    lo = [v - NEB for v in lo]; hi = [v + NEB for v in hi]
    layers.append((lo[:], hi[:]))
    lo = [(v // bf) * bf for v in lo]; hi = [(v // bf) * bf + bf - 1 for v in hi]
    layers.append((lo[:], hi[:]))
    g = (NPROPER + 1) * bf
    lo = [v - g for v in lo]; hi = [v + g for v in hi]
    layers.append((lo[:], hi[:]))
    lo = [max(0, lo[0] - stencil), max(0, lo[1] - stencil)]
    hi = [min(NX - 1, hi[0] + stencil), min(NY - 1, hi[1] + stencil)]
    layers.append((lo, hi))
    return layers, bf, stencil


def main():
    layers, bf, stencil = footprint()
    prot = layers[-1]
    boxes = []
    for bj in range(0, NY, B):
        for bi in range(0, NX, B):
            lo, hi = (bi, bj), (min(NX, bi + B) - 1, min(NY, bj + B) - 1)
            rlo = (max(-1, lo[0] - GHOST), max(-1, lo[1] - GHOST))
            rhi = (min(NX, hi[0] + GHOST), min(NY, hi[1] + GHOST))
            wet = any(not land_at(i, j) for j in range(rlo[1], rhi[1] + 1)
                      for i in range(rlo[0], rhi[0] + 1))
            own = any(not is_land(i, j) for j in range(lo[1], hi[1] + 1)
                      for i in range(lo[0], hi[0] + 1))
            touch = (lo[0] <= prot[1][0] and prot[0][0] <= hi[0] and
                     lo[1] <= prot[1][1] and prot[0][1] <= hi[1])
            kind = "kept" if wet else ("prot" if touch else "drop")
            boxes.append(dict(lo=lo, hi=hi, rlo=rlo, rhi=rhi, kind=kind, own=own))

    W = OX + NX * S + 30 + LEGEND_W
    H = OY + NY * S + 46
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
         f'font-family="Helvetica, Arial, sans-serif">',
         '<defs><pattern id="hatch" width="8" height="8" patternUnits="userSpaceOnUse" '
         f'patternTransform="rotate(45)"><line x1="0" y1="0" x2="0" y2="8" stroke="{C["drop"]}" '
         'stroke-width="1.4" stroke-opacity="0.55"/></pattern></defs>',
         f'<rect width="{W}" height="{H}" fill="{C["bg"]}"/>',
         f'<rect x="{X(0)}" y="{Y(NY)}" width="{NX * S}" height="{NY * S}" fill="{C["water"]}"/>']
    for j in range(NY):
        i = 0
        while i < NX:
            if not is_land(i, j):
                i += 1
                continue
            k = i
            while k < NX and is_land(k, j):
                k += 1
            o.append(f'<rect x="{X(i)}" y="{Y(j + 1)}" width="{(k - i) * S}" height="{S}" '
                     f'fill="{C["land"]}"/>')
            i = k
    style = {"kept": f'fill="none" stroke="{C["kept"]}" stroke-width="1.6"',
             "prot": f'fill="{C["accent_soft"]}" fill-opacity="0.16" stroke="{C["accent"]}" stroke-width="1.6"',
             "drop": f'fill="url(#hatch)" stroke="{C["drop"]}" stroke-width="1.2" stroke-dasharray="3 3"'}
    for b in boxes:
        o.append(rect(b["lo"], b["hi"], style[b["kind"]]))
    dashes = ["", "2 3", "6 3", "10 4", ""]
    for k, (lo, hi) in enumerate(layers):
        if k == 0:
            o.append(rect(lo, hi, f'fill="none" stroke="{C["fg"]}" stroke-width="2.4"'))
        else:
            w = 2.6 if k == len(layers) - 1 else 1.4
            da = f' stroke-dasharray="{dashes[k]}"' if dashes[k] else ""
            o.append(rect(lo, hi, f'fill="none" stroke="{C["accent"]}" stroke-width="{w}"{da}'))
    # ghost footprints of two example boxes: one land box kept by water in its ghosts, one dropped
    # examples far enough apart that their ghost footprints do not overlap
    ex_ghost = max((b for b in boxes if b["kind"] == "kept" and not b["own"]),
                   key=lambda b: b["lo"][1], default=None)
    ex_prot = min((b for b in boxes if b["kind"] == "prot"),
                  key=lambda b: abs(b["lo"][1] - 32), default=None)
    ex_drop = min((b for b in boxes if b["kind"] == "drop"),
                  key=lambda b: (abs(b["lo"][1] - 32), -b["lo"][0]), default=None)
    for b, lab in ((ex_ghost, "A"), (ex_prot, "B"), (ex_drop, "C")):
        if b is None:
            continue
        o.append(rect(b["rlo"], b["rhi"], f'fill="none" stroke="{C["fg"]}" stroke-width="2" '
                                           'stroke-dasharray="6 4"'))
        cx = (X(b["lo"][0]) + X(b["hi"][0] + 1)) / 2
        cy = (Y(b["hi"][1] + 1) + Y(b["lo"][1])) / 2
        o.append(f'<circle cx="{cx}" cy="{cy}" r="11" fill="{C["bg"]}" stroke="{C["fg"]}"/>'
                 f'<text x="{cx}" y="{cy + 5}" text-anchor="middle" font-size="14" font-weight="700" '
                 f'fill="{C["fg"]}">{lab}</text>')
    for i in range(0, NX + 1, 16):
        o.append(f'<text x="{X(i)}" y="{Y(0) + 18}" text-anchor="middle" font-size="12" '
                 f'fill="{C["muted"]}">{i}</text>')
    for j in range(0, NY + 1, 16):
        o.append(f'<text x="{X(0) - 6}" y="{Y(j) + 4}" text-anchor="end" font-size="12" '
                 f'fill="{C["muted"]}">{j}</text>')

    # legend
    lx, ly = OX + NX * S + 30, OY + 4
    lines = [("swatch", C["water"], "water"), ("swatch", C["land"], "land"),
             ("box", "kept", "kept: water in its ghost footprint"),
             ("box", "prot", "kept: touches the refinement footprint"),
             ("box", "drop", "dropped"), ("gap",), ("ind", None, "refinement indicator"),
             ("layer", 1, f"+ tag buffer ({NEB})"), ("layer", 2, f"rounded to blocks of {bf}"),
             ("layer", 3, f"+ (n_proper+1) blocks ({(NPROPER + 1) * bf})"),
             ("layer", 4, f"+ fill stencil ({stencil}): protects"), ("gap",),
             ("ghost", None, "ghost footprints of A, B, C"), ("gap",),
             ("note", None, f"boxes {B}x{B}, n_error_buf {NEB},"),
             ("note", None, f"n_proper {NPROPER}, blocking factor {BF1}"),
             ("note", None, f"(level 1), ref_ratio {RR}")]
    y = ly
    for item in lines:
        if item[0] == "gap":
            y += 10
            continue
        kind, a, label = item
        if kind == "swatch":
            o.append(f'<rect x="{lx}" y="{y}" width="16" height="16" fill="{a}"/>')
        elif kind == "box":
            o.append(f'<rect x="{lx}" y="{y}" width="16" height="16" {style[a]}/>')
        elif kind == "ind":
            o.append(f'<rect x="{lx}" y="{y + 1}" width="16" height="14" fill="none" stroke="{C["fg"]}" '
                     'stroke-width="2.4"/>')
        elif kind == "layer":
            w = 2.6 if a == 4 else 1.4
            da = f' stroke-dasharray="{dashes[a]}"' if dashes[a] else ""
            o.append(f'<line x1="{lx}" y1="{y + 8}" x2="{lx + 22}" y2="{y + 8}" stroke="{C["accent"]}" '
                     f'stroke-width="{w}"{da}/>')
        elif kind == "ghost":
            o.append(f'<line x1="{lx}" y1="{y + 8}" x2="{lx + 22}" y2="{y + 8}" stroke="{C["fg"]}" '
                     'stroke-width="2" stroke-dasharray="6 4"/>')
        tx = lx if kind == "note" else lx + 30
        col = C["muted"] if kind == "note" else C["fg"]
        o.append(f'<text x="{tx}" y="{y + 13}" font-size="13" fill="{col}">{label}</text>')
        y += 24
    o.append("</svg>")
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "drop_dry_boxes.svg")
    with open(out, "w") as f:
        f.write("\n".join(o) + "\n")
    n_drop = sum(b["kind"] == "drop" for b in boxes)
    print(f"wrote {out}: {len(boxes) - n_drop} of {len(boxes)} boxes kept")


if __name__ == "__main__":
    main()
