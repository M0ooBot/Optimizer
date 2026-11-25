"""Review figures from the reimported wing.step (skill: Verify step 5).

Run from the project root:
  .venv/bin/python runs/cap-radius/builds/swept55/sharp/render.py
Every face is evaluated from the B-spline data read back from wing.step
(poles and knots of the reimported faces, checked against OCC's own D1 at
random points), sampled far finer than the pixel grid, and drawn with a
per-pixel depth buffer: no mesh, no mesh normals, no edge lines. Zebra
stripes use the exact CAD normal at every sample, orthographic views, the
stripe axis mirrored for the left half. Writes PNGs into figures/.
"""
import pathlib
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.interpolate import BSpline  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import verify as V  # noqa: E402

warnings.simplefilter("ignore")
OUT = HERE / "figures"


class Evaluator:
    """Vectorised evaluation of a reimported B-spline face from its own poles and knots."""

    def __init__(self, face):
        bs = face.ad.BSpline()
        self.rev = face.rev
        nu, nv = bs.NbUPoles(), bs.NbVPoles()
        self.P = np.array([[[bs.Pole(i, j).X(), bs.Pole(i, j).Y(), bs.Pole(i, j).Z()] for j in range(1, nv + 1)]
                           for i in range(1, nu + 1)])
        self.du, self.dv = bs.UDegree(), bs.VDegree()
        self.tu = np.repeat([bs.UKnot(i) for i in range(1, bs.NbUKnots() + 1)], [bs.UMultiplicity(i) for i in range(1, bs.NbUKnots() + 1)])
        self.tv = np.repeat([bs.VKnot(i) for i in range(1, bs.NbVKnots() + 1)], [bs.VMultiplicity(i) for i in range(1, bs.NbVKnots() + 1)])
        assert not bs.IsURational() and not bs.IsVRational()
        k, t = self.du, self.tu
        self.Pu = k * (self.P[1:] - self.P[:-1]) / (t[k + 1:-1] - t[1:-k - 1])[:, None, None]
        k, t = self.dv, self.tv
        self.Pv = k * (self.P[:, 1:] - self.P[:, :-1]) / (t[k + 1:-1] - t[1:-k - 1])[None, :, None]

    @staticmethod
    def _bm(x, t, k):
        x = np.clip(x, t[0], t[-1])
        return BSpline.design_matrix(x, t, k).toarray()

    def grid(self, us, vs):
        Bu, Bv = self._bm(us, self.tu, self.du), self._bm(vs, self.tv, self.dv)
        Bu1, Bv1 = self._bm(us, self.tu[1:-1], self.du - 1), self._bm(vs, self.tv[1:-1], self.dv - 1)
        ev = lambda A, Q, B: np.stack([A @ Q[..., c] @ B.T for c in range(3)], -1)
        p, su, sv = ev(Bu, self.P, Bv), ev(Bu1, self.Pu, Bv), ev(Bu, self.Pv, Bv1)
        n = np.cross(su, sv)
        nn = np.linalg.norm(n, axis=-1, keepdims=True)
        n = n / np.where(nn == 0, 1, nn)
        return p, (-n if self.rev else n), nn[..., 0]


def check_evaluator(face, ev):
    rng = np.random.default_rng(1)
    worst_p, worst_n = 0.0, 0.0
    for _ in range(40):
        u = rng.uniform(face.u0, face.u1)
        v = rng.uniform(face.v0, face.v1)
        J = face.jet(np.array([u]), np.array([v]))
        p, n, _ = ev.grid(np.array([u]), np.array([v]))
        worst_p = max(worst_p, float(np.linalg.norm(p[0, 0] - J["p"][0])))
        worst_n = max(worst_n, float(np.linalg.norm(n[0, 0] - J["n"][0])))
    return worst_p, worst_n


def camera(elev, azim):
    e, a = np.radians(elev), np.radians(azim)
    back = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])   # towards the viewer
    right = np.array([-np.sin(a), np.cos(a), 0.0])
    up = np.cross(back, right)
    return right, up, -back


def render(samples, elev, azim, deg, freq, frame, px=1400):
    """samples: list of (p, n). frame: (cx, cy, half_width) in screen mm. Returns an RGB image."""
    right, up, view = camera(elev, azim)
    cx, cy, hw = frame
    W = px
    H = int(px * 0.72)
    hh = hw * H / W
    img = np.ones((H, W)) * 1.0
    depth = np.full(H * W, np.inf)
    val = np.full(H * W, np.nan)
    for p, n in samples:
        p = p.reshape(-1, 3); n = n.reshape(-1, 3)
        front = (n @ view) < 0
        p, n = p[front], n[front]
        xs = (p @ right - cx + hw) / (2 * hw) * W
        ys = (cy + hh - p @ up) / (2 * hh) * H
        ok = (xs >= 0) & (xs < W) & (ys >= 0) & (ys < H)
        xs, ys, p, n = xs[ok], ys[ok], p[ok], n[ok]
        idx = ys.astype(int) * W + xs.astype(int)
        d = p @ view
        r = view - 2 * (n @ view)[:, None] * n
        a = V.stripe_axis(deg, p[:, 1])
        phi = freq * (r * a).sum(-1)
        z = (np.floor(2 * phi) % 2)
        lam = 0.75 + 0.25 * np.abs(n @ view)
        order = np.lexsort((d, idx))
        idx, d, z, lam = idx[order], d[order], z[order], lam[order]
        first = np.ones(len(idx), bool)
        first[1:] = idx[1:] != idx[:-1]
        idx, d, z, lam = idx[first], d[first], z[first], lam[first]
        better = d < depth[idx]
        depth[idx[better]] = d[better]
        val[idx[better]] = (0.1 + 0.85 * z[better]) * lam[better]
    img = val.reshape(H, W)
    # fill pixels that no sample landed on from the nearest sampled pixel (at most 2 px away);
    # anything further is background (pale blue, so white stripes stay visible)
    from scipy.ndimage import distance_transform_edt
    empty = np.isnan(img)
    dist, (iy, ix) = distance_transform_edt(empty, return_indices=True)
    filled = img[iy, ix]
    img = np.where(empty & (dist <= 2.0), filled, img)
    rgb = np.repeat(np.nan_to_num(img, nan=0.0)[..., None], 3, -1)
    bg = np.isnan(img)
    rgb[bg] = (0.78, 0.84, 0.92)
    return rgb


def main():
    OUT.mkdir(exist_ok=True)
    F = V.classify(V.read(V.STEP))
    evs = {("cap", s): Evaluator(f) for s, f in F["caps"].items()}
    evs.update({("skin", k): Evaluator(f) for k, f in F["skins"].items()})
    checks = {f"{a}_{b}": check_evaluator(F["caps"][b] if a == "cap" else F["skins"][b], e) for (a, b), e in evs.items()}
    print("evaluator vs OCC D1 (max point mm, max normal):", {k: [f"{x:.1e}" for x in v] for k, v in checks.items()})

    def cap(s, nu=2400, nw=2400):
        us = np.concatenate([np.linspace(0, 0.03, nu // 4), np.linspace(0.03, 0.97, nu), np.linspace(0.97, 1, nu // 4)])
        ws = np.linspace(0, 1, nw)
        p, n, _ = evs[("cap", s)].grid(us, ws)
        return p, n

    def skin(k, v0, nu=2400, nv=600):
        us = np.concatenate([np.linspace(0, 0.03, nu // 4), np.linspace(0.03, 1, nu)])
        p, n, _ = evs[("skin", k)].grid(us, np.linspace(v0, 390, nv))
        return p, n

    capR, capL = cap("R"), cap("L")
    near = [skin(k, 280) for k in ("Rup", "Rlo")]
    full = [skin(k, 0, nu=2000, nv=1500) for k in ("Rup", "Lup", "Rlo", "Llo")]
    figs = []
    jobs = [
        ("wing_top_zebra_0", full + [capR, capL], 90, -90, 0, 14, (0, 0, 0)),
        ("wing_top_zebra_60", full + [capR, capL], 90, -90, 60, 14, (0, 0, 0)),
        ("wing_front_oblique_zebra", full + [capR, capL], 35, 180, 90, 14, (0, 0, 0)),
        ("tip_top_zebra", near + [capR], 90, -90, 0, 24, (630, 400, 85)),
        ("tip_top_45_zebra", near + [capR], 90, -90, 45, 24, (630, 400, 85)),
        ("tip_front_oblique_zebra", near + [capR], 35, -160, 90, 24, (0, 0, 0)),
        ("tip_rear_oblique_zebra", near + [capR], 30, 30, 0, 24, (0, 0, 0)),
        ("tip_outboard_low_zebra", near + [capR], 10, 80, 135, 10, (0, 0, 0)),
        ("tip_aft_lip_top_zebra", near + [capR], 90, -90, 30, 24, (670, 425, 25)),
        ("tip_bottom_zebra", near + [capR], -90, -90, 0, 24, (630, -400, 85)),
    ]
    for name, samples, elev, azim, deg, freq, frame in jobs:
        if frame[2] == 0:   # frame automatically: the right tip, or the whole wing
            right, up, _ = camera(elev, azim)
            pts = capR[0].reshape(-1, 3)[::97] if not name.startswith("wing") else \
                np.concatenate([s[0].reshape(-1, 3)[::997] for s in samples])
            xs, ys = pts @ right, pts @ up
            frame = (0.5 * (xs.min() + xs.max()), 0.5 * (ys.min() + ys.max()),
                     max(0.54 * (xs.max() - xs.min()), 0.54 * (ys.max() - ys.min()) / 0.72))
        img = render(samples, elev, azim, deg, freq, frame)
        fig, ax = plt.subplots(figsize=(14, 10.08), dpi=100)
        ax.imshow(img, interpolation="nearest")
        ax.set_axis_off()
        ax.set_title(f"{name}: elev {elev}, azim {azim}, stripe axis {deg} deg, {freq} stripes/unit (CAD normals, depth-buffered)")
        fn = OUT / f"{name}.png"
        fig.savefig(fn, bbox_inches="tight"); plt.close(fig); figs.append(fn)
    # printable-radius map on the right cap (plan view of the upper half)
    m = np.load(HERE / "map_right.npz")
    u, w, ratio = m["u"], m["w"], m["ratio"]
    sel = w <= 0.5
    P, _, _ = evs[("cap", "R")].grid(u, w[sel])
    fig, ax = plt.subplots(figsize=(12, 7), dpi=110)
    cs = ax.pcolormesh(P[..., 0], P[..., 1], np.log10(np.minimum(ratio[:, sel], 1e3)), shading="gouraud", cmap="viridis", vmin=0, vmax=1)
    ax.contour(P[..., 0], P[..., 1], ratio[:, sel], levels=[1.05, 1.5, 3.0], colors=["r", "w", "0.3"], linewidths=0.8)
    ax.set_aspect("equal")
    ax.set_xlabel("X (mm)"); ax.set_ylabel("Y (mm)")
    fig.colorbar(cs, label="log10(smallest principal radius / rule minimum at that chord position)")
    ax.set_title("right cap, upper half: radius over the printable rule (contours 1.05 red, 1.5 white, 3 grey)")
    fn = OUT / "radius_ratio_map.png"
    fig.savefig(fn, bbox_inches="tight"); plt.close(fig); figs.append(fn)
    print("\n".join(str(x) for x in figs))


if __name__ == "__main__":
    main()
