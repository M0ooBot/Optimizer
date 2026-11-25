"""Verify wing.step on the reimported STEP (skill: Verify the finished cap, steps 1-4).

Run from the project root:
  .venv/bin/python runs/cap-radius/builds/swept55/sharp/verify.py
Writes readback.json (and the figures via render.py separately).

Everything below is measured on the STEP read back from disk, using the CAD
surfaces' own derivatives (OCC BRepAdaptor_Surface.D2), not on cap_model.py.
Chord position of a cap point = its u parameter (the cap's u is the parent's
chord parameter q; chord fraction U = 3u^2 - 2u^3), judged against the
parent's attachment thickness t(U) from parent.py.
"""
import hashlib
import json
import pathlib
import sys
import warnings

import numpy as np
from scipy.optimize import minimize

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cap_model as cm  # noqa: E402

from OCP.BRepAdaptor import BRepAdaptor_Surface  # noqa: E402
from OCP.BRepCheck import BRepCheck_Analyzer  # noqa: E402
from OCP.BRepGProp import BRepGProp  # noqa: E402
from OCP.BRepBndLib import BRepBndLib  # noqa: E402
from OCP.Bnd import Bnd_Box  # noqa: E402
from OCP.GProp import GProp_GProps  # noqa: E402
from OCP.GeomAbs import GeomAbs_BSplineSurface  # noqa: E402
from OCP.IFSelect import IFSelect_RetDone  # noqa: E402
from OCP.STEPControl import STEPControl_Reader  # noqa: E402
from OCP.ShapeAnalysis import ShapeAnalysis_FreeBounds  # noqa: E402
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_REVERSED, TopAbs_SHELL, TopAbs_SOLID  # noqa: E402
from OCP.TopExp import TopExp_Explorer  # noqa: E402
from OCP.TopoDS import TopoDS  # noqa: E402
from OCP.gp import gp_Pnt, gp_Vec  # noqa: E402

warnings.simplefilter("ignore")
STEP = HERE / "wing.step"
PARENT_STEP = cm.PARENT_DIR / "parent.step"
P = cm.parent


def U(q):
    return 3 * q**2 - 2 * q**3


def t_att(q):
    return cm.thickness_at_attachment(q)


def rule_need(t):
    """Smallest radius the rule allows at attachment thickness t (thin: the half-round t/2)."""
    return np.where(t >= 0.8, 0.4, t / 2)


# ------------------------------------------------------------------ CAD access
def read(path):
    r = STEPControl_Reader()
    assert r.ReadFile(str(path)) == IFSelect_RetDone
    r.TransferRoots()
    return r.OneShape()


def count(shape, kind):
    ex, n = TopExp_Explorer(shape, kind), 0
    while ex.More():
        n += 1
        ex.Next()
    return n


class Face:
    def __init__(self, face):
        self.face = face
        self.ad = BRepAdaptor_Surface(face)
        self.rev = face.Orientation() == TopAbs_REVERSED
        self.u0, self.u1 = self.ad.FirstUParameter(), self.ad.LastUParameter()
        self.v0, self.v1 = self.ad.FirstVParameter(), self.ad.LastVParameter()
        self.kind = self.ad.GetType()

    def d2(self, u, v):
        p, a, b, aa, bb, ab = gp_Pnt(), gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec()
        self.ad.D2(float(u), float(v), p, a, b, aa, bb, ab)
        f = lambda x: np.array([x.X(), x.Y(), x.Z()])
        return f(p), f(a), f(b), f(aa), f(ab), f(bb)

    def jet(self, us, vs):
        us, vs = np.broadcast_arrays(np.asarray(us, float), np.asarray(vs, float))
        out = np.array([np.concatenate(self.d2(u, v)) for u, v in zip(us.ravel(), vs.ravel())]).reshape(us.shape + (6, 3))
        J = {"p": out[..., 0, :], "su": out[..., 1, :], "sv": out[..., 2, :], "suu": out[..., 3, :],
             "suv": out[..., 4, :], "svv": out[..., 5, :]}
        n = np.cross(J["su"], J["sv"])
        n /= np.linalg.norm(n, axis=-1, keepdims=True)
        J["n"] = -n if self.rev else n          # outward for a correctly oriented solid
        return J

    def point(self, u, v):
        return self.d2(u, v)[0]


def fundamental(J):
    su, sv, n = J["su"], J["sv"], J["n"]
    E = (su * su).sum(-1); F = (su * sv).sum(-1); G = (sv * sv).sum(-1)
    L = (J["suu"] * n).sum(-1); M = (J["suv"] * n).sum(-1); N = (J["svv"] * n).sum(-1)
    return E, F, G, L, M, N


def principal(J):
    E, F, G, L, M, N = fundamental(J)
    den = E * G - F * F
    H = (E * N - 2 * F * M + G * L) / (2 * den)
    K = (L * N - M * M) / den
    d = np.sqrt(np.maximum(H * H - K, 0))
    return H + d, H - d, K


def II_of(J, a, b):
    """Second fundamental form (outward normal) of 3D tangent vectors a, b."""
    E, F, G, L, M, N = fundamental(J)
    su, sv = J["su"], J["sv"]
    den = E * G - F * F

    def coords(t):
        ts, tv = (t * su).sum(-1), (t * sv).sum(-1)
        return (G * ts - F * tv) / den, (E * tv - F * ts) / den
    a1, a2 = coords(a)
    b1, b2 = coords(b)
    return L * a1 * b1 + M * (a1 * b2 + a2 * b1) + N * a2 * b2


def unit(x):
    return x / np.linalg.norm(x, axis=-1, keepdims=True)


# ------------------------------------------------------------------ checks
def classify(shape):
    faces = []
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        faces.append(Face(TopoDS.Face_s(ex.Current())))
        ex.Next()
    out = {"skins": {}, "caps": {}}
    for f in faces:
        assert f.kind == GeomAbs_BSplineSurface
        pm = f.point(0.5 * (f.u0 + f.u1), 0.5 * (f.v0 + f.v1))
        if abs(f.v1 - 390.0) < 1e-9 and f.v0 == 0.0:
            side = "R" if f.point(0.5, 390.0)[1] > 0 else "L"
            out["skins"][side + ("up" if pm[2] > 0 else "lo")] = f
        else:
            out["caps"]["R" if pm[1] > 0 else "L"] = f
    assert len(out["skins"]) == 4 and len(out["caps"]) == 2, out
    return out


def topology(shape):
    an = BRepCheck_Analyzer(shape)
    fb = ShapeAnalysis_FreeBounds(shape, 1e-7, False, False)
    free_wires = count(fb.GetClosedWires(), TopAbs_EDGE) + count(fb.GetOpenWires(), TopAbs_EDGE)
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props, 1e-9)
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape, box, False, False)
    b = box.Get()
    return {"solids": count(shape, TopAbs_SOLID), "shells": count(shape, TopAbs_SHELL), "faces": count(shape, TopAbs_FACE),
            "brepcheck_valid": bool(an.IsValid()), "free_edges": int(free_wires),
            "volume_mm3": props.Mass(), "bounding_box_mm": [float(x) for x in b]}


def parent_preserved(F):
    ref = classify_parent()
    worst = 0.0
    us = np.linspace(0, 1, 41)
    vs = np.linspace(0, 390, 14)
    for k, f in F["skins"].items():
        g = ref[k]
        for u in us:
            for v in vs:
                worst = max(worst, float(np.linalg.norm(f.point(u, v) - g.point(u, v))))
    # and against the exact parent definition
    dev_py = 0.0
    for k, f in F["skins"].items():
        side = +1 if k.endswith("up") else -1
        sgn = +1 if k.startswith("R") else -1
        for u in np.linspace(0, 1, 101):
            for v in np.linspace(0, 390, 14):
                dev_py = max(dev_py, float(np.linalg.norm(f.point(u, v) - P.skin(u, sgn * v, side))))
    return {"retained_skin_vs_parent_step_max_mm": worst, "retained_skin_vs_parent_py_max_mm": dev_py,
            "samples_per_face": len(us) * len(vs)}


def classify_parent():
    shape = read(PARENT_STEP)
    out = {}
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        f = Face(TopoDS.Face_s(ex.Current()))
        if f.kind == GeomAbs_BSplineSurface:
            pm = f.point(0.5, 195.0)
            out[("R" if f.point(0.5, 390.0)[1] > 0 else "L") + ("up" if pm[2] > 0 else "lo")] = f
        ex.Next()
    return out


def outline_and_dims(F):
    res = {}
    for side, f in F["caps"].items():
        sg = 1 if side == "R" else -1
        us = np.linspace(0, 1, 2001)
        lip = np.array([f.point(u, 0.5) for u in us])
        ref = P.tip_outline(np.linspace(0, 1, 200001), side=sg)
        from scipy.spatial import cKDTree
        d, _ = cKDTree(ref).query(lip)
        same_p = np.linalg.norm(lip - P.tip_outline(us, side=sg), axis=1)
        J = f.jet(us[1:-1], np.full(len(us) - 2, 0.5))
        tang = unit(J["su"])
        nz = np.abs(J["n"][:, 2]).max()
        n_dot_t = np.abs((J["n"] * tang).sum(-1)).max()
        # the lip must be the plan extreme: no cap point beyond the outline (|Y| above the outline's span there)
        grid_u = np.linspace(0, 1, 301)
        grid_w = np.linspace(0, 1, 301)
        pts = np.array([[f.point(u, w) for w in grid_w] for u in grid_u])
        # plan containment: every cap point lies inboard of the outline at its own chord fraction
        r = np.abs(pts[..., 1]).ravel()
        xi = (pts[..., 0].ravel() - cm.SWEEP * r) / (cm.ROOT_CHORD - cm.TAPER * r)
        xi = np.clip(xi, 0, 1)
        pp = np.array([float(np.real(np.roots([-2, 3, 0, -x])[np.argmin(np.abs(np.roots([-2, 3, 0, -x]) - 0.5) + 10 * (np.abs(np.imag(np.roots([-2, 3, 0, -x]))) > 1e-9) + 10 * ((np.real(np.roots([-2, 3, 0, -x])) < -1e-9) | (np.real(np.roots([-2, 3, 0, -x])) > 1 + 1e-9)))])) for x in xi])
        excess = r - (390 + 240 * pp * (1 - pp))
        res[side + "_plan_excess_beyond_outline_max_mm"] = float(excess.max())
        res[side] = {"lip_to_outline_max_mm": float(d.max()), "lip_vs_outline_same_parameter_max_mm": float(same_p.max()),
                     "lip_normal_z_max": float(nz), "lip_normal_dot_lip_tangent_max": float(n_dot_t),
                     "max_abs_Y_mm": float(np.abs(pts[..., 1]).max()), "lip_samples": len(us)}
    # attachment station: the shared edges lie at |Y| = 390
    att = 0.0
    for side, f in F["caps"].items():
        for w in (0.0, 1.0):
            for u in np.linspace(0, 1, 501):
                att = max(att, abs(abs(f.point(u, w)[1]) - 390.0))
    res["attachment_edges_max_abs_Y_minus_390_mm"] = att
    return res


# ------------------------------------------------------------------ joins
def joins(F):
    """Separation, tangent-plane angle and curvature (common frame) along all four wing-to-cap joins."""
    us = np.concatenate([np.geomspace(1e-4, 0.02, 40), np.linspace(0.02, 0.98, 481), 1 - np.geomspace(0.02, 1e-4, 40)])
    out = {}
    for side in ("R", "L"):
        cap = F["caps"][side]
        for skin, w in (("up", 0.0), ("lo", 1.0)):
            par = F["skins"][side + skin]
            Jp = par.jet(us, np.full_like(us, 390.0))
            Jc = cap.jet(us, np.full_like(us, w))
            sep = np.linalg.norm(Jp["p"] - Jc["p"], axis=-1)
            ang = np.degrees(np.arccos(np.clip((Jp["n"] * Jc["n"]).sum(-1), -1, 1)))
            e1 = unit(Jp["su"])                       # along the join
            e2 = unit(np.cross(Jp["n"], e1))          # across it, in the common tangent plane
            Kp = np.stack([II_of(Jp, e1, e1), II_of(Jp, e1, e2), II_of(Jp, e2, e2)], -1)
            Kc = np.stack([II_of(Jc, e1, e1), II_of(Jc, e1, e2), II_of(Jc, e2, e2)], -1)
            dk = np.abs(Kp - Kc).max(-1)
            scale = np.maximum(np.abs(Kp).max(-1), 1e-6)
            i = int(np.argmax(dk))
            core = (us > 0.01) & (us < 0.99)
            out[side + skin] = {
                "samples": len(us), "u_range": [float(us[0]), float(us[-1])],
                "endpoint_distance_mm": [float(np.linalg.norm(Jp["p"][0] - P.skin(0, 390, 1))),
                                         float(np.linalg.norm(Jp["p"][-1] - P.skin(1, 390, 1)))],
                "separation_max_mm": float(sep.max()),
                "tangent_angle_max_deg": float(ang.max()), "tangent_angle_max_at_u": float(us[np.argmax(ang)]),
                "curvature_diff_max_per_mm": float(dk.max()), "at_u": float(us[i]),
                "parent_curvature_there_per_mm": float(np.abs(Kp[i]).max()),
                "curvature_diff_rel_max_core": float((dk / scale)[core].max()),
                "curvature_diff_rel_max_all": float((dk / scale).max()),
            }
    return out


# ------------------------------------------------------------------ sections
FLAT = 1e-5   # 1/mm: curvature below this (radius over 100 m) counts as straight when counting reversals
def trace_plane(F, side_faces, A, n_p, guess_u):
    """Points and signed section curvature along the plane through A with normal n_p:
    upper skin |Y| 370 -> 390, cap w 0 -> 1, lower skin 390 -> 370 (right tip)."""
    segs = [(side_faces["up"], np.linspace(370, 390, 81)[:-1], "v"),
            (side_faces["cap"], np.concatenate([np.linspace(0, 0.4, 400), np.linspace(0.4, 0.6, 1601), np.linspace(0.6, 1, 400)]), "w"),
            (side_faces["lo"], np.linspace(390, 370, 81)[1:], "v")]
    rows = []
    u = guess_u
    for face, ss, kind in segs:
        for s in ss:
            for _ in range(40):
                p, su, sv, *_ = face.d2(u, s)
                f = float(n_p @ (p - A))
                du = -f / float(n_p @ su)
                u = min(max(u + du, face.u0), face.u1)
                if abs(du) < 1e-14:
                    break
            J = face.jet(np.array([u]), np.array([s]))
            rows.append((kind, face is side_faces["cap"], u, s, J))
    pts = np.array([r[4]["p"][0] for r in rows])
    res = []
    for kind, is_cap, u, s, J in rows:
        n = J["n"][0]
        t = unit(np.cross(n, n_p))
        Nc = unit(np.cross(n_p, t))
        if Nc @ n < 0:
            Nc = -Nc
        kn = float(II_of({k: v[0] for k, v in J.items()}, t, t))
        k_sec = kn / float(Nc @ n)          # Meusnier: curvature of the planar section, signed (outward n)
        k1, k2, _ = principal({k: v[0] for k, v in J.items()})
        res.append({"cap": is_cap, "u": u, "s": s, "p": J["p"][0], "convexity": -k_sec, "kn": -kn,
                    "kmax": max(abs(k1), abs(k2))})
    return res, pts


def sections(F):
    side_faces = {"up": F["skins"]["Rup"], "lo": F["skins"]["Rlo"], "cap": F["caps"]["R"]}
    cap = F["caps"]["R"]
    uk = sorted(set(np.round([cap.ad.BSpline().UKnot(i) for i in range(1, cap.ad.BSpline().NbUKnots() + 1)], 12)))
    uk = np.array(uk)
    at_knots = uk[(uk > 0) & (uk < 1)]
    between = 0.5 * (uk[:-1] + uk[1:])
    thin = np.array([0.003, 0.006, 0.0125, 0.9, 0.93, 0.96, 0.98, 0.99, 0.995])
    qs = np.unique(np.concatenate([at_knots, between, thin]))
    out = []
    for q in qs:
        xi = U(q)
        A = P.skin(q, 390.0, +1)
        g = np.array([cm.SWEEP - cm.TAPER * xi, 1.0, 0.0])     # guide direction (follows the sweep)
        n_p = unit(np.array([g[1], -g[0], 0.0]))
        res, pts = trace_plane(F, side_faces, A, n_p, q)
        conv = np.array([r["convexity"] for r in res])
        kn = np.array([r["kn"] for r in res])
        cap_mask = np.array([r["cap"] for r in res])
        rad = 1 / np.abs(conv)
        rad_n = 1 / np.maximum(np.abs(kn), 1e-12)
        sgn = np.sign(np.where(np.abs(conv) < FLAT, 0, conv))   # |k| < FLAT counts as straight
        nz = sgn != 0
        s_nz = sgn[nz]
        rev_idx = np.where(s_nz[1:] * s_nz[:-1] < 0)[0]
        rev_pts = pts[nz][rev_idx]
        i = int(np.argmin(np.where(cap_mask, rad, np.inf)))
        lip = int(np.argmin(np.where(cap_mask, np.abs(pts[:, 2]), np.inf)))
        tt = float(t_att(q))
        need = float(rule_need(tt))
        out.append({
            "u": float(q), "chord_fraction": float(xi), "attachment_thickness_mm": tt, "rule_min_radius_mm": need,
            "min_section_radius_mm": float(rad[i]), "at_xyz": pts[i].round(4).tolist(),
            "min_section_radius_at_lip": bool(abs(i - lip) <= 2),
            "lip_section_radius_mm": float(rad[lip]),
            "min_normal_radius_along_trace_mm": float(rad_n[cap_mask].min()),
            "lip_normal_radius_along_trace_mm": float(rad_n[lip]),
            "min_principal_radius_on_trace_mm": float(1 / max(r["kmax"] for r, c in zip(res, cap_mask) if c)),
            "reversals": int(len(rev_idx)), "reversal_xyz": rev_pts.round(4).tolist()[:6],
            "max_concave_curvature_per_mm": float(max(0.0, -conv.min())),
            "max_concave_at_xyz": pts[int(np.argmin(conv))].round(4).tolist() if conv.min() < 0 else None,
            "max_concave_on_cap": bool(cap_mask[int(np.argmin(conv))]) if conv.min() < 0 else None,
            "section_meets_rule": bool(rad[i] >= need * (1 - 1e-6)),
            "surface_meets_rule_on_trace": bool(rad_n[cap_mask].min() >= need * (1 - 1e-6)),
        })
    return out


def lip_normal_sections(F):
    """Cuts perpendicular to the lip through lip points: the section that sees the lip's true radius.
    Traced on the cap face from the lip outwards (upper and lower) until the cut leaves the face."""
    f = F["caps"]["R"]
    out = []
    for p in np.concatenate([np.geomspace(0.003, 0.05, 6), np.linspace(0.08, 0.92, 22), 1 - np.geomspace(0.05, 0.003, 6)]):
        A = f.point(p, 0.5)
        Jl = f.jet(np.array([p]), np.array([0.5]))
        n_p = unit(Jl["su"][0])                      # plane normal = lip tangent
        rows = []
        for direction in (-1, +1):
            u = p
            for k in range(1, 4000):
                w = 0.5 + direction * k * 2.5e-4
                if not 0 <= w <= 1:
                    break
                ok = False
                for _ in range(40):
                    q, su, *_ = f.d2(u, w)
                    du = -float(n_p @ (q - A)) / float(n_p @ su)
                    u = u + du
                    if not 0 < u < 1:
                        break
                    if abs(du) < 1e-14:
                        ok = True
                        break
                if not ok:
                    break
                rows.append((w, u))
        rows.append((0.5, p))
        rows.sort()
        ws, us = np.array(rows).T
        J = f.jet(us, ws)
        n = J["n"]
        t = unit(np.cross(n, n_p))
        Nc = unit(np.cross(n_p[None, :], t))
        Nc *= np.sign((Nc * n).sum(-1))[:, None]
        kn = II_of(J, t, t)
        conv = -kn / (Nc * n).sum(-1)
        rad = 1 / np.abs(conv)
        il = int(np.argmin(np.abs(ws - 0.5)))
        i = int(np.argmin(rad))
        s = np.where(np.abs(conv) < FLAT, 0, np.sign(conv))
        s = s[s != 0]
        span = np.linalg.norm(J["p"][-1] - J["p"][0])
        tt = float(t_att(p))
        need = float(rule_need(tt))
        out.append({"lip_u": float(p), "attachment_thickness_mm": tt, "rule_min_radius_mm": need,
                    "lip_radius_mm": float(rad[il]), "min_radius_mm": float(rad[i]), "min_at_lip": bool(abs(i - il) <= 2),
                    "min_at_xyz": J["p"][i].round(4).tolist(), "reversals": int(np.sum(s[1:] * s[:-1] < 0)),
                    "max_concave_curvature_per_mm": float(max(0.0, -conv.min())),
                    "traced_w_range": [float(ws[0]), float(ws[-1])], "chord_of_trace_mm": float(span),
                    "meets_rule": bool(rad[i] >= need * (1 - 1e-6))})
    return out


# ------------------------------------------------------------------ principal map
def cap_map(f, nu=500, nw=500):
    us = np.unique(np.concatenate([np.geomspace(1e-4, 0.02, 60), np.linspace(0.02, 0.98, nu), 1 - np.geomspace(0.02, 1e-4, 60)]))
    ws = np.linspace(0, 1, nw + 1)
    Uu, Ww = np.meshgrid(us, ws, indexing="ij")
    J = f.jet(Uu, Ww)
    k1, k2, K = principal(J)
    kmax = np.maximum(np.abs(k1), np.abs(k2))
    rad = 1 / kmax
    tt = t_att(Uu)
    need = rule_need(tt)
    ratio = rad / need

    def ratio_at(x):
        u, w = float(np.clip(x[0], 1e-4, 1 - 1e-4)), float(np.clip(x[1], 0, 1))
        Jx = f.jet(np.array([u]), np.array([w]))
        a, b, _ = principal(Jx)
        r = 1 / max(abs(a[0]), abs(b[0]))
        return r / float(rule_need(t_att(u))), r
    # refine every local minimum of the ratio (local maximum of curvature relative to the rule)
    from scipy.ndimage import minimum_filter
    loc = (ratio == minimum_filter(ratio, size=5, mode="nearest"))
    cand = np.argwhere(loc)
    cand = cand[np.argsort(ratio[loc])][:40]
    refined = []
    for i, j in cand:
        r0 = minimize(lambda x: ratio_at(x)[0], [Uu[i, j], Ww[i, j]], method="Nelder-Mead",
                      options={"xatol": 1e-7, "fatol": 1e-9, "maxiter": 400})
        rr, rabs = ratio_at(r0.x)
        u, w = float(np.clip(r0.x[0], 1e-4, 1 - 1e-4)), float(np.clip(r0.x[1], 0, 1))
        refined.append({"u": u, "w": w, "xyz": f.point(u, w).round(4).tolist(), "radius_mm": rabs, "ratio_to_rule": rr,
                        "attachment_thickness_mm": float(t_att(u))})
    refined.sort(key=lambda d: d["ratio_to_rule"])
    thick = tt >= 0.8
    it = np.unravel_index(np.argmin(np.where(thick, rad, np.inf)), rad.shape)
    ia = np.unravel_index(np.argmin(rad), rad.shape)
    # lip radius along the thin trailing-edge zone vs t/2
    thin_te = (us > 0.5) & (t_att(us) < 0.8)
    lipj = len(ws) // 2
    lip_ratio = (rad[thin_te, lipj] / (t_att(us[thin_te]) / 2))
    minw_ratio = (rad[thin_te].min(1) / (t_att(us[thin_te]) / 2))
    pole_le = float(np.linalg.norm(f.point(us[0], 0.5) - f.point(0, 0.5)))
    pole_te = float(np.linalg.norm(f.point(us[-1], 0.5) - f.point(1, 0.5)))
    return {
        "grid_u_w": [len(us), len(ws)], "u_range": [float(us[0]), float(us[-1])],
        "grid_reach_to_collapsed_points_mm": {"leading_edge": pole_le, "trailing_edge": pole_te},
        "smallest_ratio_radius_to_rule": refined[0] if refined else None,
        "next_smallest_ratios": refined[1:6],
        "smallest_radius_where_t_ge_0p8": {"radius_mm": float(rad[it]), "u": float(Uu[it]), "w": float(Ww[it]),
                                            "xyz": f.point(Uu[it], Ww[it]).round(4).tolist(),
                                            "attachment_thickness_mm": float(tt[it])},
        "smallest_radius_anywhere_on_grid": {"radius_mm": float(rad[ia]), "u": float(Uu[ia]), "w": float(Ww[ia]),
                                              "attachment_thickness_mm": float(tt[ia])},
        "grid_points_below_rule": int((ratio < 1 - 1e-6).sum()),
        "thin_te_lip_radius_over_half_t": [float(lip_ratio.min()), float(lip_ratio.max())],
        "thin_te_min_wrap_radius_over_half_t": [float(minw_ratio.min()), float(minw_ratio.max())],
        "saddle_fraction": float((K < 0).mean()),
    }, (us, ws, rad, ratio)


# ------------------------------------------------------------------ zebra
CAMERAS = {"top": np.array([340.0, 0.0, 1400.0]), "front_above": np.array([-500.0, 0.0, 700.0]),
           "rear_above": np.array([1400.0, 0.0, 700.0]), "right_oblique": np.array([900.0, 1100.0, 500.0])}
STRIPE_ANGLES = [0, 45, 90, 135]
FREQ = 12.0


def stripe_axis(deg, y):
    a = np.array([np.cos(np.radians(deg)), np.sin(np.radians(deg)), 0.35])
    a = np.broadcast_to(a, y.shape + (3,)).copy()
    a[..., 1] *= np.where(y < 0, -1, 1)         # mirrored lighting for the mirrored wing
    return unit(a)


def zebra_fields(J, cam, deg):
    """Stripe phase and its surface gradient from CAD normals and Weingarten derivatives."""
    p, n, su, sv = J["p"], J["n"], J["su"], J["sv"]
    E, F_, G, L, M, N = fundamental(J)
    den = E * G - F_ * F_
    nu = (((F_ * M - G * L) / den)[..., None] * su + ((F_ * L - E * M) / den)[..., None] * sv)
    nv = (((F_ * N - G * M) / den)[..., None] * su + ((F_ * M - E * N) / den)[..., None] * sv)
    d = p - cam
    dist = np.linalg.norm(d, axis=-1, keepdims=True)
    v = d / dist
    vu = (su - (v * su).sum(-1, keepdims=True) * v) / dist
    vv = (sv - (v * sv).sum(-1, keepdims=True) * v) / dist
    vn = (v * n).sum(-1, keepdims=True)
    r = v - 2 * vn * n
    ru = vu - 2 * (((vu * n).sum(-1, keepdims=True) + (v * nu).sum(-1, keepdims=True)) * n + vn * nu)
    rv = vv - 2 * (((vv * n).sum(-1, keepdims=True) + (v * nv).sum(-1, keepdims=True)) * n + vn * nv)
    a = stripe_axis(deg, p[..., 1])
    phi = FREQ * (r * a).sum(-1)
    pu, pv = FREQ * (ru * a).sum(-1), FREQ * (rv * a).sum(-1)
    grad = (((G * pu - F_ * pv) / den)[..., None] * su + ((E * pv - F_ * pu) / den)[..., None] * sv)
    return phi, grad


def zebra(F):
    us = np.concatenate([np.geomspace(1e-3, 0.02, 30), np.linspace(0.02, 0.98, 961), 1 - np.geomspace(0.02, 1e-3, 30)])
    res = {"cameras_mm": {k: v.tolist() for k, v in CAMERAS.items()}, "stripe_angles_deg": STRIPE_ANGLES,
           "stripes_per_unit_reflection": FREQ, "u_range": [float(us[0]), float(us[-1])], "joins": {}}
    worst = {"phase": (0, None), "grad_rel": (0, None), "dir_deg_at_crossing": (0, None), "grad_abs": (0, None)}
    for side in ("R", "L"):
        cap = F["caps"][side]
        for skin, w in (("up", 0.0), ("lo", 1.0)):
            Jp = F["skins"][side + skin].jet(us, np.full_like(us, 390.0))
            Jc = cap.jet(us, np.full_like(us, w))
            jr = {"phase_max": 0.0, "grad_rel_max": 0.0, "grad_abs_max": 0.0, "crossings": 0, "dir_change_max_deg": 0.0}
            for cname, cam in CAMERAS.items():
                for deg in STRIPE_ANGLES:
                    fp, gp = zebra_fields(Jp, cam, deg)
                    fc, gc = zebra_fields(Jc, cam, deg)
                    dphi = np.abs(fp - fc)
                    gnorm = np.linalg.norm(gp, axis=-1)
                    dg = np.linalg.norm(gp - gc, axis=-1)
                    rel = dg / np.maximum(gnorm, 1e-9)
                    sig = gnorm > 0.05 * np.median(gnorm)
                    tag = f"{side}{skin} {cname} {deg}deg"
                    for key, val, idx in (("phase", dphi.max(), int(np.argmax(dphi))),
                                          ("grad_rel", rel[sig].max(), int(np.flatnonzero(sig)[np.argmax(rel[sig])])),
                                          ("grad_abs", dg.max(), int(np.argmax(dg)))):
                        if val > worst[key][0]:
                            worst[key] = (float(val), {"where": tag, "u": float(us[idx]), "xyz": Jp["p"][idx].round(3).tolist()})
                    jr["phase_max"] = max(jr["phase_max"], float(dphi.max()))
                    jr["grad_rel_max"] = max(jr["grad_rel_max"], float(rel[sig].max()))
                    jr["grad_abs_max"] = max(jr["grad_abs_max"], float(dg.max()))
                    # light/dark crossings on the join: phase through a half-integer
                    h = 2 * fp
                    idx = np.where(np.floor(h[1:]) != np.floor(h[:-1]))[0]
                    for k in idx:
                        sp = unit(np.cross(Jp["n"][k], gp[k]))      # stripe (iso-phase) direction
                        sc = unit(np.cross(Jc["n"][k], gc[k]))
                        dd = float(np.degrees(np.arccos(np.clip(abs(sp @ sc), -1, 1))))
                        jr["crossings"] += 1
                        jr["dir_change_max_deg"] = max(jr["dir_change_max_deg"], dd)
                        if dd > worst["dir_deg_at_crossing"][0]:
                            worst["dir_deg_at_crossing"] = (dd, {"where": tag, "u": float(us[k]), "xyz": Jp["p"][k].round(3).tolist()})
            res["joins"][side + skin] = jr
    res["worst"] = {k: {"value": v[0], **(v[1] or {})} for k, v in worst.items()}
    return res


def main():
    data = STEP.read_bytes()
    shape = read(STEP)
    F = classify(shape)
    out = {"file": "wing.step", "step_sha256": hashlib.sha256(data).hexdigest(), "step_bytes": len(data)}
    out["topology"] = topology(shape)
    out["parent_skin"] = parent_preserved(F)
    out["locked"] = outline_and_dims(F)
    out["joins"] = joins(F)
    out["zebra"] = zebra(F)
    out["sections_right_tip"] = sections(F)
    out["lip_normal_sections_right_tip"] = lip_normal_sections(F)
    out["principal_map"] = {}
    for side in ("R", "L"):
        m, raw = cap_map(F["caps"][side])
        out["principal_map"][side] = m
        if side == "R":
            np.savez_compressed(HERE / "map_right.npz", u=raw[0], w=raw[1], radius=raw[2].astype(np.float32),
                                ratio=raw[3].astype(np.float32))
    secs = out["sections_right_tip"]
    out["sections_summary"] = {
        "count": len(secs),
        "all_sections_meet_rule_by_section_radius": all(s["section_meets_rule"] for s in secs),
        "all_sections_meet_rule_by_surface_radius": all(s["surface_meets_rule_on_trace"] for s in secs),
        "worst_section_radius_ratio": min(s["min_section_radius_mm"] / s["rule_min_radius_mm"] for s in secs),
        "worst_surface_radius_ratio": min(s["min_normal_radius_along_trace_mm"] / s["rule_min_radius_mm"] for s in secs),
        "sections_with_reversals": [s["u"] for s in secs if s["reversals"]],
        "sections_min_not_at_lip": [s["u"] for s in secs if not s["min_section_radius_at_lip"]],
        "max_concave_curvature_per_mm": max(s["max_concave_curvature_per_mm"] for s in secs),
    }
    ln = out["lip_normal_sections_right_tip"]
    out["lip_normal_summary"] = {
        "count": len(ln), "all_meet_rule": all(s["meets_rule"] for s in ln),
        "worst_ratio": min(s["min_radius_mm"] / s["rule_min_radius_mm"] for s in ln),
        "min_not_at_lip": [s["lip_u"] for s in ln if not s["min_at_lip"]],
        "with_reversals": [s["lip_u"] for s in ln if s["reversals"]],
        "max_concave_curvature_per_mm": max(s["max_concave_curvature_per_mm"] for s in ln),
    }
    (HERE / "readback.json").write_text(json.dumps(out, indent=1, default=float))
    summary = {k: out[k] for k in ("step_sha256", "topology", "parent_skin", "locked", "sections_summary", "lip_normal_summary")}
    summary["joins"] = {k: {kk: v[kk] for kk in ("separation_max_mm", "tangent_angle_max_deg", "curvature_diff_max_per_mm", "curvature_diff_rel_max_core")} for k, v in out["joins"].items()}
    summary["zebra_worst"] = out["zebra"]["worst"]
    summary["map"] = {s: {k: m[k] for k in ("smallest_ratio_radius_to_rule", "smallest_radius_where_t_ge_0p8", "grid_points_below_rule", "thin_te_lip_radius_over_half_t", "thin_te_min_wrap_radius_over_half_t", "grid_reach_to_collapsed_points_mm")} for s, m in out["principal_map"].items()}
    print(json.dumps(summary, indent=1, default=float))


if __name__ == "__main__":
    main()
