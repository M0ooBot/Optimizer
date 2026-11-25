"""Render zebra (reflection-stripe) images of wing.step from CAD normals, no edge lines.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/fin40/blunt/render_zebra.py

Orthographic camera and fixed stripe lighting.  Each face of the reimported
STEP is sampled on a dense (u, v) grid; every sample's colour comes from its own
CAD normal (reflection of parallel light bands), and samples are splatted into a
depth-buffered image, so no mesh facets or edge lines are drawn.  The flat
trailing-edge face is shaded blue-grey.
Writes zebra_side.png, zebra_tip_oblique.png, zebra_join_mid.png, zebra_join_le.png and
zebra_te_closure.png.
"""
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from OCP.gp import gp_Pnt, gp_Vec  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import occ_util as ou  # noqa: E402


def camera(az, el):
    a, e = np.radians(az), np.radians(el)
    eye = np.array([np.sin(a) * np.cos(e), np.cos(a) * np.cos(e), np.sin(e)])
    d = -eye
    right = np.cross(d, [0, 0, 1.0])
    right /= np.linalg.norm(right)
    up = np.cross(right, d)
    return d, right, up


def sample(F, ubox, vbox, nu, nv):
    us = np.linspace(*ubox, nu)
    vs = np.linspace(*vbox, nv)
    P = np.full((nu, nv, 3), np.nan)
    N = np.full((nu, nv, 3), np.nan)
    p, su, sv = gp_Pnt(), gp_Vec(), gp_Vec()
    for i, u in enumerate(us):
        for j, v in enumerate(vs):
            if not F.inside(u, v):
                continue
            F.surf.D1(u, v, p, su, sv)
            n = np.cross([su.X(), su.Y(), su.Z()], [sv.X(), sv.Y(), sv.Z()])
            nn = np.linalg.norm(n)
            if nn < 1e-14:
                continue
            P[i, j] = (p.X(), p.Y(), p.Z())
            N[i, j] = F.sign * n / nn
    return P.reshape(-1, 3), N.reshape(-1, 3)


def stripes(N, d, axis, K):
    a = axis / np.linalg.norm(axis)
    b = np.cross(a, d)
    b /= np.linalg.norm(b)
    c = np.cross(a, b)
    r = d - 2 * (N @ d)[:, None] * N
    return (np.sin(K * np.arctan2(r @ c, r @ b)) > 0).astype(float)


def render(F, boxes, view, axis_deg, K, center, half, fname, title, px=900, flat=("te",)):
    d, right, up = camera(*view)
    th = np.radians(axis_deg)
    axis = np.cos(th) * right + np.sin(th) * up
    cx, cy = center @ right, center @ up
    img = np.ones((px, px, 3)) * np.array([0.93, 0.94, 0.96])
    zbuf = np.full((px, px), np.inf)
    scale = px / (2 * half)
    for name, (ubox, vbox, nu, nv) in boxes.items():
        P, N = sample(F[name], ubox, vbox, nu, nv)
        ok = ~np.isnan(P[:, 0])
        P, N = P[ok], N[ok]
        vis = (N @ d) < 0
        P, N = P[vis], N[vis]
        X = ((P @ right) - cx) * scale + px / 2
        Y = px / 2 - ((P @ up) - cy) * scale
        D = P @ d
        ix, iy = np.floor(X).astype(int), np.floor(Y).astype(int)
        inb = (ix >= 0) & (ix < px) & (iy >= 0) & (iy < px)
        ix, iy, D, N = ix[inb], iy[inb], D[inb], N[inb]
        shade = 0.75 + 0.25 * np.clip(-(N @ d), 0, 1)
        if name in flat:
            col = np.stack([0.5 * shade, 0.55 * shade, 0.68 * shade], 1)
        else:
            z = stripes(N, d, axis, K)
            g = (0.1 + 0.85 * z) * shade
            col = np.stack([g, g, g], 1)
        order = np.argsort(-D)          # far first, near last wins
        ix, iy, D, col = ix[order], iy[order], D[order], col[order]
        # depth test: keep nearest per pixel
        lin = iy * px + ix
        cur = zbuf.ravel()[lin]
        closer = D <= cur + 1e-9
        lin, D, col = lin[closer], D[closer], col[closer]
        # later (nearer) samples overwrite earlier ones for the same pixel
        flat_img = img.reshape(-1, 3)
        flat_img[lin] = col
        zb = zbuf.ravel()
        zb[lin] = np.minimum(zb[lin], D)
    # fill pixels that received no sample (sub-pixel gaps between samples) from the
    # nearest sampled pixel within 2 px; wider gaps stay background
    from scipy import ndimage
    empty = ~np.isfinite(zbuf)
    dist, (iy_n, ix_n) = ndimage.distance_transform_edt(empty, return_indices=True)
    fill = empty & (dist <= 2.0)
    img[fill] = img[iy_n[fill], ix_n[fill]]
    fig, ax = plt.subplots(figsize=(7, 7), dpi=px / 7)
    ax.imshow(img, extent=(cx - half, cx + half, cy - half, cy + half), interpolation="nearest")
    ax.set_title(title, fontsize=8)
    ax.set_xlabel("screen right (mm)", fontsize=7)
    ax.set_ylabel("screen up (mm)", fontsize=7)
    ax.tick_params(labelsize=6)
    fig.tight_layout()
    fig.savefig(HERE / fname)
    plt.close(fig)


def main():
    F = ou.classify(ou.read_step(HERE / "wing.step"))
    tu = sorted(F["te"].uv(p)[0] for p in ((150.2, 0, 232.0), (149.9, 0, 241.5)))
    # 1. whole fin, oblique from +Y
    render(F, {"parent_up": ((0, 1), (0, 240), 900, 1400), "cap_up": ((0, np.pi), (0, 1), 1200, 900)},
           (-25, 15), 80, 24, np.array([80.0, 0, 140.0]), 150.0, "zebra_side.png",
           "whole fin, view az -25 el 15 (+Y side); stripe axis 80 deg, 24 per turn")
    # 2. tip, oblique from +Y, forward-high
    render(F, {"parent_up": ((0, 1), (190, 240), 1600, 1200), "cap_up": ((0, np.pi), (0, 1), 2000, 1600)},
           (35, 25), 60, 24, np.array([108.0, 0, 255.0]), 50.0, "zebra_tip_oblique.png",
           "tip, view az 35 el 25 (+Y side); stripe axis 60 deg, 24 per turn")
    # 3. attachment join at mid-chord, close
    render(F, {"parent_up": ((0.37, 0.53), (236, 240), 1800, 700), "cap_up": ((1.22, 1.58), (0, 0.03), 1800, 700)},
           (20, -20), 100, 60, np.array([100.0, 0, 240.0]), 3.0, "zebra_join_mid.png",
           "attachment join (z = 240) near mid-chord, view az 20 el -20; stripe axis 100 deg, 60 per turn")
    # 4. attachment join near the leading edge, close
    render(F, {"parent_up": ((0, 0.4), (232, 240), 1600, 900), "cap_up": ((0, 0.9), (0, 0.2), 1800, 1400)},
           (-40, 10), 75, 30, np.array([68.0, 0, 241.0]), 6.0, "zebra_join_le.png",
           "attachment join near the leading edge (z = 240), view az -40 el 10; stripe axis 75 deg, 30 per turn")
    # 5. trailing-edge face closure, from aft and +Y
    render(F, {"parent_up": ((0.9, 1.0), (234.0, 240.0), 1200, 900), "cap_up": ((np.pi - 0.35, np.pi), (0.0, 0.08), 1200, 1200),
               "te": (tuple(tu), (-0.2, 0.2), 900, 300)},
           (70, 15), 30, 24, np.array([149.95, 0.0, 240.6]), 1.8, "zebra_te_closure.png",
           "trailing-edge face closure (flat face blue-grey); view az 70 el 15; stripe axis 30 deg, 24 per turn")


if __name__ == "__main__":
    main()
