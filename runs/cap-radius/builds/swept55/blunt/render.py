"""Zebra renders of the reimported wing.step (review in context, skill verification step 5).

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/swept55/blunt/render.py
Writes figures/*.png.  Every face of the STEP is a B-spline surface or a plane; the B-spline
faces are evaluated from their STEP poles (NumPy), so the stripes come from the exported
geometry.  Orthographic camera, no edge lines, no overlays.  Stripes: reflected view ray r,
phase asin(r . e) with stripe width 4 deg; e is mirrored in Y (mirrored stripe lighting).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bspline_fit import basis_matrix  # noqa: E402
import verify as V  # noqa: E402

from OCP.BRepAdaptor import BRepAdaptor_Surface  # noqa: E402
from OCP.GeomAbs import GeomAbs_BSplineSurface  # noqa: E402

FIG = os.path.join(HERE, "figures")


def bspline_data(face):
    a = BRepAdaptor_Surface(face)
    s = a.BSpline()
    P = np.array([[(s.Pole(i, j).X(), s.Pole(i, j).Y(), s.Pole(i, j).Z()) for j in range(1, s.NbVPoles() + 1)]
                  for i in range(1, s.NbUPoles() + 1)])

    def kv(n, val, mul):
        return np.concatenate([[val(i)] * mul(i) for i in range(1, n + 1)])
    ku = kv(s.NbUKnots(), s.UKnot, s.UMultiplicity)
    kvv = kv(s.NbVKnots(), s.VKnot, s.VMultiplicity)
    flip = face.Orientation().name == "TopAbs_REVERSED"
    return P, ku, kvv, s.UDegree(), s.VDegree(), flip, (a.FirstUParameter(), a.LastUParameter(),
                                                         a.FirstVParameter(), a.LastVParameter())


def dbasis(t, k, deg):
    """First-derivative basis matrix."""
    Nm = basis_matrix(t, k, deg - 1)
    n = len(k) - deg - 1
    D = np.zeros((len(t), n))
    for i in range(n):
        a = k[i + deg] - k[i]
        b = k[i + deg + 1] - k[i + 1]
        if a > 0:
            D[:, i] += deg / a * Nm[:, i]
        if b > 0:
            D[:, i] -= deg / b * Nm[:, i + 1]
    return D


def sample_face(face, nu, nv, urange=None, vrange=None):
    P, ku, kv, du, dv, flip, (u0, u1, v0, v1) = bspline_data(face)
    s = np.linspace(0, 1, nu); t = np.linspace(0, 1, nv)
    ua, ub = urange or (u0, u1); va, vb = vrange or (v0, v1)
    U = ua + (ub - ua) * s; Vv = va + (vb - va) * t
    # normalise to the knot range used by basis_matrix ([k0, k_end] -> works for any range)
    Bu, Bv = basis_matrix(U, ku, du), basis_matrix(Vv, kv, dv)
    Du, Dv = dbasis(U, ku, du), dbasis(Vv, kv, dv)
    def tp(A, B):
        tmp = np.tensordot(A, P, axes=(1, 0))            # (i, b, c)
        return np.einsum("ibc,jb->ijc", tmp, B, optimize=True)
    X = tp(Bu, Bv)
    Xu = tp(Du, Bv)
    Xv = tp(Bu, Dv)
    N = np.cross(Xu, Xv)
    N /= np.maximum(np.linalg.norm(N, axis=-1, keepdims=True), 1e-300)
    if flip:
        N = -N
    ok = np.isfinite(N).all(-1) & (np.linalg.norm(np.cross(Xu, Xv), axis=-1) > 1e-14)
    return X[ok], N[ok]


def view_frame(d, right=None):
    d = np.asarray(d, float); d /= np.linalg.norm(d)
    if right is None:
        up = np.array([0, 0, 1.0]) if abs(d[2]) < 0.9 else np.array([-1.0, 0, 0])
        x = np.cross(up, -d)
    else:
        x = np.asarray(right, float) - np.dot(right, d) * d
    x /= np.linalg.norm(x)
    y = np.cross(-d, x)
    return d, x, y


def adaptive_samples(faces, d, center, width_mm, px, right=None, cap=8_000_000):
    """Sample each (face, urange) densely enough for ~0.35 px spacing in this view."""
    d_, ex, ey = view_frame(d, right)
    scale = 2 * px / width_mm
    out = []
    for f, ur in faces:
        X, _ = sample_face(f, 80, 80, urange=ur)
        g = X.reshape(80, 80, 3) if len(X) == 6400 else None
        if g is None:
            nu, nv = 1500, 600
        else:
            q = np.stack([(g - center) @ ex, (g - center) @ ey], -1) * scale
            du = np.linalg.norm(np.diff(q, axis=0), axis=-1).max() * 80
            dv = np.linalg.norm(np.diff(q, axis=1), axis=-1).max() * 80
            nu, nv = int(np.clip(du / 0.3, 60, 12000)), int(np.clip(dv / 0.3, 60, 12000))
            while nu * nv > cap:
                nu, nv = int(nu * 0.9), int(nv * 0.9)
        out.append(sample_face(f, nu, nv, urange=ur))
    return out


def render(samples, d, e_axis, center, width_mm, px=1400, aspect=0.76, stripe_deg=4.0, name="view", right=None):
    d, ex, ey = view_frame(d, right)
    H = int(px * aspect)
    ss = 2                                      # supersampling
    W2, H2 = px * ss, H * ss
    zbuf = np.full((H2, W2), np.inf)
    col = np.full((H2, W2), 0.82)
    scale = W2 / width_mm
    for X, N in samples:
        rel = X - center
        sx = (rel @ ex) * scale + W2 / 2
        sy = H2 / 2 - (rel @ ey) * scale
        depth = rel @ d
        front = (N @ d) < 0
        inside = front & (sx >= 0) & (sx < W2) & (sy >= 0) & (sy < H2)
        if not inside.any():
            continue
        ix, iy = sx[inside].astype(int), sy[inside].astype(int)
        dep = depth[inside]
        n = N[inside]
        e = np.tile(np.asarray(e_axis, float), (len(n), 1))
        e[:, 1] *= np.where(X[inside][:, 1] >= 0, 1, -1)          # mirrored lighting
        r = d - 2 * (n @ d)[:, None] * n
        ph = np.degrees(np.arcsin(np.clip((r * e).sum(-1) / np.linalg.norm(e, axis=1), -1, 1))) / stripe_deg
        val = np.where(np.floor(ph) % 2 == 0, 0.08, 0.97)
        order = np.argsort(-dep)                 # far first, near last overwrite
        lin = iy[order] * W2 + ix[order]
        zb = zbuf.ravel(); cb = col.ravel()
        # keep the nearest sample per pixel
        uniq, first = np.unique(lin[::-1], return_index=True)
        keep = order[::-1][first]
        lin_k = iy[keep] * W2 + ix[keep]
        closer = dep[keep] < zb[lin_k]
        zb[lin_k[closer]] = dep[keep][closer]
        cb[lin_k[closer]] = val[keep][closer]
    img = col.reshape(H2, W2).reshape(H, ss, px, ss).mean(axis=(1, 3))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    plt.imsave(os.path.join(FIG, name + ".png"), img, cmap="gray", vmin=0, vmax=1)


def main():
    import cap_model as cm
    shape = V.read(V.STEP)
    _, faces = V.topology(shape)
    cls = V.classify(faces)
    bs = [f for f in faces if BRepAdaptor_Surface(f).GetType() == GeomAbs_BSplineSurface]
    allf = [(f, None) for f in bs]
    views = [
        # name, faces, view dir, stripe axis, centre, width, right-vector
        ("wing_top_zebra", allf, (0, 0, -1), (1, 0.35, 0), (345.0, 0, 0), 980, (0, 1, 0)),
        ("wing_top_zebra_spanwise", allf, (0, 0, -1), (0.25, 1, 0), (345.0, 0, 0), 980, (0, 1, 0)),
        ("tip_top_zebra", allf, (0, 0, -1), (1, 0.35, 0), (625.0, 418.0, 0), 160, (0, 1, 0)),
        ("tip_top_zebra_spanwise", allf, (0, 0, -1), (0.25, 1, 0), (625.0, 418.0, 0), 160, (0, 1, 0)),
        ("tip_bottom_zebra", allf, (0, 0, 1), (1, 0.35, 0), (625.0, 418.0, 0), 160, (0, 1, 0)),
        ("tip_oblique_aft_zebra", allf, (-0.55, -0.45, -0.70), (1, 0.6, 0.3), (640.0, 425.0, 0), 140, None),
    ]
    for name, ff, d, e, c, w, right in views:
        if os.path.exists(os.path.join(FIG, name + ".png")) and "--all" not in sys.argv:
            continue
        S = adaptive_samples(ff, d, np.array(c), w, 1400, right)
        render(S, d, e, np.array(c), w, name=name, right=right)
        print(name, flush=True)
    # TE corner (half-round) close-up
    up, lo = cls["cap"]["R_upper"], cls["cap"]["R_lower"]
    a = BRepAdaptor_Surface(up)
    u0, u1 = a.FirstUParameter(), a.LastUParameter()
    ur = (u1 - 0.25 * (u1 - u0), u1)
    te_par = [f for f in cls["parent"] if f in bs and V.bbox(f)[4] > 1]
    ff = [(up, ur), (lo, ur)] + [(f, (0.9, 1.0)) for f in te_par]
    c0 = np.array([cm.P0_XY[0] + 0.3, cm.P0_XY[1] + 0.35, 0.0])
    # the planar half-round face (flat: uniform stripe), sampled in polar coordinates
    rr, th = np.meshgrid(np.linspace(0, cm.TE_HALF, 400), np.linspace(-np.pi / 2, np.pi / 2, 800), indexing="ij")
    ea = np.r_[cm.E_A, 0.0]; ez = np.array([0, 0, 1.0])
    P0 = np.r_[cm.P0_XY, 0.0]
    disc = (P0 + rr[..., None] * (np.cos(th)[..., None] * ea + np.sin(th)[..., None] * ez)).reshape(-1, 3)
    disc_n = np.tile(-np.r_[cm.E_B, 0.0], (len(disc), 1))
    for name, d, e, w, right in [("te_corner_oblique_zebra", (-0.45, 0.55, -0.70), (0.3, 1, 0.4), 2.0, None),
                                 ("te_corner_from_aft_zebra", (-0.66, 0.75, -0.25), (0, 0.3, 1), 2.0, None),
                                 ("te_corner_top_zebra", (0, 0, -1), (1, 0.35, 0), 2.0, (0, 1, 0))]:
        S = adaptive_samples(ff, d, c0, w, 1400, right) + [(disc, disc_n)]
        render(S, d, e, c0, w, name=name, right=right)
        print(name, flush=True)


if __name__ == "__main__":
    main()
