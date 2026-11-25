"""Zebra renderings of the reimported wing.step (verification step 5).

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/flat/blunt/render_zebra.py

Each pixel's stripe value is computed from the exact CAD normal of the STEP
face (B-spline poles and knots), sampled densely in each face's parameter
space and splatted with a depth buffer. No mesh normals, no edge lines.
Stripes: reflected view ray r, stripe value = floor(asin(r . a) / period + 1/4) mod 2.
For the left half the stripe axis a is mirrored (a_y -> -a_y), so the two tips
should look like mirror images with a centred camera.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import verify_cap as V  # noqa: E402
from OCP.TopAbs import TopAbs_FACE
from OCP.TopoDS import TopoDS
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.GeomAbs import GeomAbs_Plane

OUTDIR = os.path.join(HERE, 'images')


def load():
    shape = V.read(V.STEP)
    items = []
    for f in V.explore(shape, TopAbs_FACE, TopoDS.Face_s):
        ad = BRepAdaptor_Surface(f)
        if ad.GetType() == GeomAbs_Plane:
            items.append(('plane', f))
        else:
            items.append(('bs', V.Surf(ad.BSpline(), f)))
    return items


def samples(items, nu, nv, region=None):
    pts, nrm = [], []
    for kind, S in items:
        if kind == 'plane':
            # closure faces: tiny, flat (normal +X); sample their bounding rectangle inside the boundary
            continue
        u0, u1 = S.urange; v0, v1 = S.vrange
        is_cap = abs(v0 + 1) < 1e-12 and abs(v1 - 1) < 1e-12 and S.ku == 5
        if is_cap:
            u = np.unique(np.concatenate([np.linspace(u0, u1, nu), 1 - np.geomspace(1e-7, 0.05, nu // 8), np.geomspace(1e-6, 0.05, nu // 8)]))
            e = np.linspace(0, 1, nv)
            v = np.unique(np.concatenate([np.linspace(v0, v1, nv), np.sqrt(1 - e), -np.sqrt(1 - e), 1 - np.geomspace(1e-7, 0.05, nv // 4), -1 + np.geomspace(1e-7, 0.05, nv // 4)]))
        elif S.ku == 5:  # parent skin: s x Y
            u = np.unique(np.concatenate([np.linspace(u0, u1, nu), 1 - np.geomspace(1e-7, 0.05, nu // 8), np.geomspace(1e-6, 0.05, nu // 8)]))
            vv = np.linspace(v0, v1, 4 * nv)
            if region is not None:
                vv = vv[(np.abs(vv) >= region[0])]
                vv = np.unique(np.concatenate([vv, np.sign(v1) * np.linspace(region[0], 350, nv), -np.linspace(region[0], 350, nv)]))
            v = vv
        else:  # TE face: Y x t
            u = np.linspace(u0, u1, 4 * nv); v = np.linspace(v0, v1, 9)
            if region is not None:
                u = np.unique(np.concatenate([np.linspace(region[0], 350, nv), -np.linspace(region[0], 350, nv)]))
        d = S.grid(u, v)
        n = np.cross(d['Pu'], d['Pv']); ln = np.linalg.norm(n, axis=-1, keepdims=True)
        n = n / np.where(ln > 0, ln, 1)
        # orient outward: the solid is roughly the slab |Z| < 6, |Y| < 400 around the chord line
        P = d['P']
        c = np.stack([np.clip(P[..., 0], 2, 118), np.clip(P[..., 1], -340, 340), np.zeros_like(P[..., 0])], -1)
        if (np.nansum((P - c) * n) < 0):
            n = -n
        pts.append(P.reshape(-1, 3)); nrm.append(n.reshape(-1, 3))
    return np.concatenate(pts), np.concatenate(nrm)


def render(P, N, view, up, center, width_mm, px=(1600, 900), period_deg=4.0, stripe_deg=0.0, mode='zebra', light=None):
    view = np.asarray(view, float); view /= np.linalg.norm(view)        # direction toward the eye
    up = np.asarray(up, float); right = np.cross(up, view); right /= np.linalg.norm(right); upv = np.cross(view, right)
    W, H = px; scale = W / width_mm
    q = P - center
    x = (q @ right) * scale + W / 2; y = H / 2 - (q @ upv) * scale; depth = q @ view
    ok = np.isfinite(x) & (x >= 0) & (x < W) & (y >= 0) & (y < H) & np.all(np.isfinite(N), -1)
    facing = (N @ view) > 0
    ok &= facing
    x, y, depth, Pn, Nn = x[ok], y[ok], depth[ok], P[ok], N[ok]
    eye = -view
    if mode == 'zebra':
        vn = Nn @ eye
        r = eye[None] - 2 * vn[:, None] * Nn
        ang = np.radians(stripe_deg)
        b1 = right; b2 = upv
        a = np.cos(ang) * b1 + np.sin(ang) * b2
        A = np.broadcast_to(a, r.shape).copy()
        A[Pn[:, 1] < 0, 1] *= -1     # mirrored stripe lighting for the left half
        # mirror the view as well for the left half so a centred camera sees mirror-image stripes
        val = (np.floor(np.arcsin(np.clip((r * A).sum(-1), -1, 1)) / np.radians(period_deg) + 0.25) % 2)
        col = 0.12 + 0.8 * val
    else:
        L = np.asarray(light, float); L /= np.linalg.norm(L)
        Lm = np.broadcast_to(L, Nn.shape).copy(); Lm[Pn[:, 1] < 0, 1] *= -1
        col = 0.15 + 0.8 * np.clip((Nn * Lm).sum(-1), 0, 1)
    ix = x.astype(int); iy = y.astype(int); pix = iy * W + ix
    order = np.lexsort((-depth, pix))       # per pixel, nearest (largest depth toward eye) last
    pix, col = pix[order], col[order]
    last = np.r_[pix[1:] != pix[:-1], True]
    img = np.full(W * H, np.nan)
    img[pix[last]] = col[last]
    img = img.reshape(H, W)
    # fill isolated holes from the 4-neighbourhood
    for _ in range(2):
        hole = np.isnan(img)
        nb = np.stack([np.roll(img, s, a) for s in (1, -1) for a in (0, 1)])
        cnt = np.sum(~np.isnan(nb), 0)
        fill = hole & (cnt >= 3)
        img[fill] = np.nanmean(nb[:, fill], 0)
    return img


def save(img, path, title):
    H, W = img.shape
    fig = plt.figure(figsize=(W / 100, H / 100 + 0.35), dpi=100)
    ax = fig.add_axes([0, 0, 1, H / (H + 35)])
    cm = matplotlib.colormaps['gray'].copy(); cm.set_bad((0.55, 0.62, 0.7))
    ax.imshow(img, cmap=cm, vmin=0, vmax=1, interpolation='nearest'); ax.axis('off')
    fig.text(0.01, 1 - 0.3 / (H / 100 + 0.35), title, fontsize=10, va='center')
    fig.savefig(path); plt.close(fig)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    items = load()
    out = []
    # 1. whole wing, top view, centred camera
    P, N = samples(items, 1400, 500)
    img = render(P, N, (0, 0, 1), (-1, 0, 0), np.array([60, 0, 0]), 840, px=(1800, 300), stripe_deg=90)
    p = os.path.join(OUTDIR, 'zebra_top_full.png'); save(img, p, 'wing.step top view, centred camera, mirrored stripes (stripe axis across the span)'); out.append(p)
    # 2. tips, top view, close (both tips, mirrored), two stripe angles
    P, N = samples(items, 3000, 900, region=(300,))
    for k, (sd, el) in enumerate([(90, 90), (0, 90), (90, 55)]):
        th = np.radians(el)
        view = np.array([-np.cos(th), 0, np.sin(th)])
        imgs = []
        for sgn in (-1, 1):
            imgs.append(render(P, N, view * np.array([1, sgn, 1]), (-1, 0, 0) if el == 90 else (0, 0, 1), np.array([60, sgn * 375, 0]), 140, px=(900, 700), stripe_deg=sd))
        img = np.concatenate([imgs[0], np.full((700, 6), np.nan), imgs[1]], 1)
        p = os.path.join(OUTDIR, f'zebra_tips_{k}.png'); save(img, p, f'left tip | right tip, camera elevation {el} deg, stripe angle {sd} deg, mirrored lighting'); out.append(p)
    # 3. oblique views of the right tip: from forward-outboard-above and aft-outboard-above
    for k, view in enumerate([(-0.6, 0.55, 0.58), (0.6, 0.55, 0.58), (0.0, 1.0, 0.25)]):
        img = render(P, N, view, (0, 0, 1), np.array([60, 375, 0]), 150, px=(1100, 700), stripe_deg=60)
        p = os.path.join(OUTDIR, f'zebra_right_oblique_{k}.png'); save(img, p, f'right tip, view toward eye {view}, stripes 60 deg'); out.append(p)
    # trailing-edge corner close-ups: see render_te_corner.py
    print('\n'.join(out))


if __name__ == '__main__':
    main()
