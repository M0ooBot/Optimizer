"""Verify wing.step (reimported) for the twisted/blunt rounded-tip build.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/twisted/blunt/verify_cap.py

Writes readback.json next to this script.  Every measurement is taken on the
surfaces read back from wing.step (B-spline poles/knots evaluated exactly with
SciPy; spot-checked against OCCT's own evaluator).
"""
import hashlib
import json
import os
import sys
import time
import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import minimize, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cap_model as cm  # noqa: E402

from OCP.STEPControl import STEPControl_Reader
from OCP.IFSelect import IFSelect_RetDone
from OCP.TopExp import TopExp_Explorer, TopExp
from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID, TopAbs_SHELL, TopAbs_EDGE
from OCP.TopoDS import TopoDS
from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
from OCP.BRep import BRep_Tool
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.GeomAdaptor import GeomAdaptor_Surface
from OCP.GeomAbs import GeomAbs_BSplineSurface, GeomAbs_Plane
from OCP.GeomLProp import GeomLProp_SLProps
from OCP.Geom import Geom_BSplineSurface
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.gp import gp_Pln, gp_Pnt, gp_Dir

STEP = os.path.join(HERE, "wing.step")
OUT = os.path.join(HERE, "readback.json")
PARENT_STEP = os.path.join(cm.PARENT_DIR, "parent.step")


# ---------------------------------------------------------------- B-spline surfaces
class Surf:
    """Exact evaluator for a (non-rational) Geom_BSplineSurface read back from STEP."""

    def __init__(self, gs):
        self.gs = gs
        nu, nv = gs.NbUPoles(), gs.NbVPoles()
        P = np.empty((nu, nv, 3))
        for i in range(nu):
            for j in range(nv):
                q = gs.Pole(i + 1, j + 1)
                P[i, j] = (q.X(), q.Y(), q.Z())
        self.P = P
        self.rational = gs.IsURational() or gs.IsVRational()
        self.du, self.dv = gs.UDegree(), gs.VDegree()

        def flat(n, kget, mget):
            out = []
            for k in range(1, n + 1):
                out += [kget(k)] * mget(k)
            return np.array(out)
        self.tu = flat(gs.NbUKnots(), gs.UKnot, gs.UMultiplicity)
        self.tv = flat(gs.NbVKnots(), gs.VKnot, gs.VMultiplicity)
        self.u0, self.u1 = self.tu[self.du], self.tu[-self.du - 1]
        self.v0, self.v1 = self.tv[self.dv], self.tv[-self.dv - 1]
        self._bu = BSpline(self.tu, np.eye(nu), self.du)
        self._bv = BSpline(self.tv, np.eye(nv), self.dv)

    def mats(self, b, x, lo, hi):
        x = np.clip(np.asarray(x, float), lo, hi)
        return [b(x, nu=d) for d in range(3)]

    def grid(self, u, v):
        """Derivatives on a tensor grid: returns dict of (Nu, Nv, 3) arrays."""
        Bu = self.mats(self._bu, u, self.u0, self.u1)
        Bv = self.mats(self._bv, v, self.v0, self.v1)
        nu, nv = self.P.shape[:2]
        T = [np.tensordot(B, self.P, axes=(1, 0)) for B in Bu]          # (Nu, nv, 3)
        f = lambda a, b: np.transpose(np.tensordot(T[a], Bv[b], axes=(1, 1)), (0, 2, 1))
        return dict(S=f(0, 0), Su=f(1, 0), Sv=f(0, 1), Suu=f(2, 0), Suv=f(1, 1), Svv=f(0, 2))

    def pts(self, u, v):
        """Derivatives at scattered points."""
        u = np.atleast_1d(np.asarray(u, float)); v = np.atleast_1d(np.asarray(v, float))
        out = {k: [] for k in ("S", "Su", "Sv", "Suu", "Suv", "Svv")}
        nu, nv = self.P.shape[:2]
        Pf = self.P.reshape(nu, -1)
        for c0 in range(0, len(u), 20000):
            uu, vv = u[c0:c0 + 20000], v[c0:c0 + 20000]
            Bu = self.mats(self._bu, uu, self.u0, self.u1)
            Bv = self.mats(self._bv, vv, self.v0, self.v1)
            T = [(B @ Pf).reshape(len(uu), nv, 3) for B in Bu]
            f = lambda a, b: np.einsum("njk,nj->nk", T[a], Bv[b])
            for k, (a, b) in zip(out, ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2))):
                out[k].append(f(a, b))
        return {k: np.concatenate(vl) for k, vl in out.items()}


def fundamental(d):
    Su, Sv = d["Su"], d["Sv"]
    n = np.cross(Su, Sv)
    nn = np.linalg.norm(n, axis=-1, keepdims=True)
    N = n / np.maximum(nn, 1e-300)
    E = np.sum(Su * Su, -1); F = np.sum(Su * Sv, -1); G = np.sum(Sv * Sv, -1)
    L = np.sum(d["Suu"] * N, -1); M = np.sum(d["Suv"] * N, -1); Nn = np.sum(d["Svv"] * N, -1)
    return E, F, G, L, M, Nn, N, nn[..., 0]


def principal(d):
    E, F, G, L, M, Nn, N, jac = fundamental(d)
    det = E * G - F * F
    H = (E * Nn - 2 * F * M + G * L) / (2 * det)
    K = (L * Nn - M * M) / det
    disc = np.sqrt(np.maximum(H * H - K, 0.0))
    kmax = np.abs(H) + disc
    return kmax, H, K, N, jac


def II_in_frame(d, e1, e2):
    """Second fundamental form components in a common tangent frame (e1, e2): k11, k12, k22."""
    E, F, G, L, M, Nn, N, _ = fundamental(d)
    Su, Sv = d["Su"], d["Sv"]
    # express e in (Su, Sv) by least squares
    def coords(e):
        A = np.stack([Su, Sv], -1)                         # (n, 3, 2)
        ATA = np.einsum("nki,nkj->nij", A, A)
        ATe = np.einsum("nki,nk->ni", A, e)
        det = ATA[:, 0, 0] * ATA[:, 1, 1] - ATA[:, 0, 1] ** 2
        scale = (ATA[:, 0, 0] + ATA[:, 1, 1]) ** 2
        bad = det <= 1e-14 * np.maximum(scale, 1e-300)
        ATA = np.where(bad[:, None, None], np.eye(2)[None], ATA)
        x = np.linalg.solve(ATA, ATe[..., None])[..., 0]
        x[bad] = np.nan                                     # collapsed (degenerate) parameter point
        return x
    a, b = coords(e1), coords(e2)
    def II(x, y):
        return L * x[:, 0] * y[:, 0] + M * (x[:, 0] * y[:, 1] + x[:, 1] * y[:, 0]) + Nn * x[:, 1] * y[:, 1]
    return II(a, a), II(a, b), II(b, b), N


# ---------------------------------------------------------------- shape helpers
def read_step(path):
    r = STEPControl_Reader()
    if r.ReadFile(path) != IFSelect_RetDone:
        raise RuntimeError("cannot read " + path)
    r.TransferRoots()
    return r.OneShape()


def faces_of(shape):
    out = []
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        out.append(TopoDS.Face_s(ex.Current()))
        ex.Next()
    return out


def count(shape, kind):
    n = 0
    ex = TopExp_Explorer(shape, kind)
    while ex.More():
        n += 1
        ex.Next()
    return n


def bspline_of(face):
    s = BRep_Tool.Surface_s(face)
    ad = GeomAdaptor_Surface(s)
    if ad.GetType() != GeomAbs_BSplineSurface:
        return None, ad.GetType()
    return ad.BSpline(), ad.GetType()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def loc(P):
    return [round(float(x), 5) for x in P]


# ---------------------------------------------------------------- checks
def topology(shape):
    m = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(shape, TopAbs_EDGE, TopAbs_FACE, m)
    bad = 0; degen = 0
    for i in range(1, m.Extent() + 1):
        e = TopoDS.Edge_s(m.FindKey(i))
        if BRep_Tool.Degenerated_s(e):
            degen += 1
            continue
        if m.FindFromIndex(i).Extent() != 2:
            bad += 1
    pr = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, pr, 1e-9)
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape, box, False, False)
    xmin, ymin, zmin, xmax, ymax, zmax = box.Get()
    return dict(solids=count(shape, TopAbs_SOLID), shells=count(shape, TopAbs_SHELL), faces=count(shape, TopAbs_FACE),
                edges=m.Extent(), degenerated_edges=degen, edges_not_shared_by_two_faces=bad,
                brepcheck_valid=bool(BRepCheck_Analyzer(shape).IsValid()),
                volume_mm3=float(pr.Mass()), bounding_box_mm=[xmin, xmax, ymin, ymax, zmin, zmax])


def classify(shape):
    pshape = read_step(PARENT_STEP)
    pf = []
    for f in faces_of(pshape):
        b, t = bspline_of(f)
        if b is not None:
            pf.append(Surf(b))
    faces = faces_of(shape)
    info = []
    for f in faces:
        b, t = bspline_of(f)
        if b is None:
            info.append(dict(face=f, kind="non-bspline", gtype=str(t)))
            continue
        S = Surf(b)
        match = None
        for k, P in enumerate(pf):
            if P.P.shape == S.P.shape:
                d = float(np.max(np.abs(P.P - S.P)))
                if d < 1e-6:
                    match = (k, d)
        info.append(dict(face=f, surf=S, parent=match))
    return info, pf


def match_on_iso(S, isov, P, ngs=8001):
    """For points P (n,3), the parameter u on the iso-curve v = isov of S closest to each point."""
    us = np.linspace(S.u0, S.u1, ngs)
    C = S.grid(us, [isov])["S"][:, 0]
    from scipy.spatial import cKDTree
    _, idx = cKDTree(C).query(P)
    u = us[idx]
    for _ in range(30):
        d = S.pts(u, np.full_like(u, isov))
        r = d["S"] - P
        g1 = np.sum(r * d["Su"], -1)
        g2 = np.sum(d["Su"] * d["Su"], -1) + np.sum(r * d["Suu"], -1)
        u = np.clip(u - g1 / np.where(np.abs(g2) > 1e-300, g2, 1e-300), S.u0, S.u1)
    return u


def join_check(Sa, ua, va, Sb, ub, vb, tangent_from="u"):
    """Compare face a (reference) and b at matched parameters. Frame: e1 along the join, e2 = N x e1."""
    da = Sa.pts(ua, va); db = Sb.pts(ub, vb)
    Ta = da["Su"] if tangent_from == "u" else da["Sv"]
    e1 = Ta / np.linalg.norm(Ta, axis=-1, keepdims=True)
    _, _, _, _, _, _, Na, _ = fundamental(da)
    e2 = np.cross(Na, e1)
    ka11, ka12, ka22, Na = II_in_frame(da, e1, e2)
    kb11, kb12, kb22, Nb = II_in_frame(db, e1, e2)
    sgn = np.sign(np.sum(Na * Nb, -1))
    kb11, kb12, kb22, Nb = kb11 * sgn, kb12 * sgn, kb22 * sgn, Nb * sgn[:, None]
    sep = np.linalg.norm(da["S"] - db["S"], axis=-1)
    ang = np.degrees(np.arccos(np.clip(np.sum(Na * Nb, -1), -1, 1)))
    return dict(P=da["S"], sep=sep, ang=ang, e1=e1, e2=e2, Na=Na, Nb=Nb,
                ka=(ka11, ka12, ka22), kb=(kb11, kb12, kb22))


def summarize_join(j, name, exclude=None):
    sel = np.ones(len(j["sep"]), bool) if exclude is None else ~exclude
    dk22 = np.abs(j["ka"][2] - j["kb"][2]); dk11 = np.abs(j["ka"][0] - j["kb"][0]); dk12 = np.abs(j["ka"][1] - j["kb"][1])
    ref = np.maximum(np.maximum(np.abs(j["ka"][2]), np.abs(j["ka"][0])), 1e-3)
    def mx(a):
        a2 = np.where(sel & np.isfinite(a), a, -np.inf)
        i = int(np.argmax(a2))
        return float(a2[i]), loc(j["P"][i])
    out = {"samples": int(sel.sum()), "excluded_samples": int((~sel).sum()),
           "samples_with_undefined_derivatives": int(np.sum(~np.isfinite(dk22)))}
    for nm, a in (("separation_mm", j["sep"]), ("tangent_plane_angle_deg", j["ang"]),
                  ("cross_join_normal_curvature_jump_per_mm", dk22), ("along_join_normal_curvature_jump_per_mm", dk11),
                  ("twist_jump_per_mm", dk12), ("cross_join_curvature_jump_relative", dk22 / ref)):
        v, l = mx(a)
        out[nm] = {"max": v, "at": l}
    return out


def zebra_join(j, name, c=2 * np.pi / 0.25, excl=None):
    """Numerical zebra across a join: phase mismatch and stripe direction change at actual crossings."""
    views = []
    for th in np.radians([0, 35, 70, 110, 145, 180]):
        for ph in np.radians(np.arange(0, 360, 45)):
            views.append([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
    views = np.unique(np.round(np.array(views), 12), axis=0)
    sel = np.ones(len(j["sep"]), bool) if excl is None else ~excl
    e1, e2 = j["e1"], j["e2"]
    worst_phase = (0.0, None); worst_dir = (0.0, None); worst_spacing = (0.0, None); ncross = 0
    for v in views:
        v = v / np.linalg.norm(v)
        a = np.cross(v, [0.0, 0.0, 1.0]) if abs(v[2]) < 0.9 else np.cross(v, [1.0, 0.0, 0.0])
        a /= np.linalg.norm(a); bvec = np.cross(v, a)
        for ang in np.radians(np.arange(0, 180, 30)):
            d = np.cos(ang) * a + np.sin(ang) * bvec
            vals = []
            for N, (k11, k12, k22) in ((j["Na"], j["ka"]), (j["Nb"], j["kb"])):
                vn = N @ v
                visible = vn < -0.05          # surface facing the camera (view direction v)
                r = v[None, :] - 2 * vn[:, None] * N
                psi = c * (r @ d)
                g = []
                for (ki1, ki2) in ((k11, k12), (k12, k22)):
                    dN = -(ki1[:, None] * e1 + ki2[:, None] * e2)
                    dr = -2 * ((dN @ v)[:, None] * N + vn[:, None] * dN)
                    g.append(c * (dr @ d))
                vals.append((psi, np.stack(g, -1), visible))
            (pa_, ga, va_), (pb_, gb, vb_) = vals
            ok = sel & va_ & vb_ & np.isfinite(ga).all(-1) & np.isfinite(gb).all(-1)
            if not np.any(ok):
                continue
            dph = np.abs(np.angle(np.exp(1j * (pa_ - pb_)))) / (2 * np.pi)
            i = int(np.argmax(np.where(ok, dph, -1)))
            if dph[i] > worst_phase[0]:
                worst_phase = (float(dph[i]), dict(at=loc(j["P"][i]), view=loc(v), stripe_dir=loc(d)))
            # actual light/dark boundaries crossing the join: psi passes a multiple of pi between samples
            k = np.floor(pa_ / np.pi)
            cross = np.where(ok[:-1] & ok[1:] & (k[1:] != k[:-1]))[0]
            if len(cross):
                ncross += len(cross)
                gA, gB = ga[cross], gb[cross]
                nA, nB = np.linalg.norm(gA, axis=-1), np.linalg.norm(gB, axis=-1)
                okc = (nA > 1e-9) & (nB > 1e-9)
                if np.any(okc):
                    cosang = np.abs(np.sum(gA * gB, -1)) / np.maximum(nA * nB, 1e-300)
                    dang = np.where(okc, np.degrees(np.arccos(np.minimum(1.0, cosang))), -1)
                    k1 = int(np.argmax(dang))
                    if dang[k1] > worst_dir[0]:
                        worst_dir = (float(dang[k1]), dict(at=loc(j["P"][cross[k1]]), view=loc(v), stripe_dir=loc(d)))
                    sp = np.where(okc, np.abs(nB / np.maximum(nA, 1e-300) - 1), -1)
                    k2 = int(np.argmax(sp))
                    if sp[k2] > worst_spacing[0]:
                        worst_spacing = (float(sp[k2]), dict(at=loc(j["P"][cross[k2]]), view=loc(v), stripe_dir=loc(d)))
    return {"stripe_period_in_r_dot_d": 0.25, "views": int(len(views)), "stripe_directions_per_view": 6,
            "max_phase_mismatch_fraction_of_period": worst_phase[0], "phase_worst": worst_phase[1],
            "crossings_examined": int(ncross),
            "max_stripe_direction_change_at_crossings_deg": worst_dir[0], "direction_worst": worst_dir[1],
            "max_stripe_spacing_change_at_crossings_relative": worst_spacing[0], "spacing_worst": worst_spacing[1]}


def param_grid(t, deg, nper):
    u = np.unique(t[deg:-deg])
    out = [np.linspace(u[i], u[i + 1], nper + 1)[:-1] for i in range(len(u) - 1)] + [[u[-1]]]
    return np.concatenate(out)


def chord_info(P):
    Pm = P * np.array([1.0, np.sign(P[..., 1:2].mean()) if np.ndim(P) else 1.0, 1.0]) if False else P
    Ya = np.abs(P[..., 1])
    Q = P.copy(); Q[..., 1] = Ya
    chi = cm.chord_fraction(Q, Ya)
    return chi, cm.t_att(chi), cm.rule_radius(chi), Ya - cm.CUT


def curvature_map(S, name, lo_excl, hi_excl, nper=4, refine=12, trim=None):
    u = param_grid(S.tu, S.du, nper); v = param_grid(S.tv, S.dv, nper)
    if trim is None:
        u = u[(u > lo_excl) & (u <= hi_excl)]
    d = S.grid(u, v)
    kmax, H, K, N, jac = principal(d)
    P = d["S"]
    R = 1.0 / np.maximum(kmax, 1e-300)
    chi, t, req, eta = chord_info(P)
    ratio = R / req
    if trim is not None:
        UU, VV = np.meshgrid(u, v, indexing="ij")
        ins = inside_trim(trim, UU, VV)
        R = np.where(ins, R, np.inf); ratio = np.where(ins, ratio, np.inf)
    corner = (eta < cm.ETA_F) & (chi > 0.99)
    res = {"grid": [int(len(u)), int(len(v))], "u_range_evaluated": [float(u[0]), float(u[-1])]}
    def worst(mask, key):
        a = np.where(mask, key, np.inf)
        i, j = np.unravel_index(np.argmin(a), a.shape)
        return i, j
    # local refinement of the worst grid points (by ratio) in and outside the corner zone
    cands = []
    for mask in (~corner, corner):
        if not np.any(mask):
            continue
        a = np.where(mask, ratio, np.inf)
        flat = np.argsort(a, axis=None)[:refine]
        for f in flat:
            i, j = np.unravel_index(f, a.shape)
            if np.isfinite(a[i, j]):
                cands.append((i, j))
    du = np.diff(u).max() if len(u) > 1 else 1e-3
    dv = np.diff(v).max()
    refined = []
    for i, j in cands:
        def fobj(x):
            dd = S.pts([x[0]], [x[1]])
            km = principal(dd)[0][0]
            c_, t_, r_, e_ = chord_info(dd["S"])
            return (1.0 / km) / r_[0]
        x0 = [u[i], v[j]]
        bnds = [(max(u[max(i - 2, 0)], S.u0), min(u[min(i + 2, len(u) - 1)], S.u1)),
                (max(v[max(j - 2, 0)], S.v0), min(v[min(j + 2, len(v) - 1)], S.v1))]
        r = minimize(fobj, x0, method="Nelder-Mead", options=dict(xatol=1e-9, fatol=1e-10, maxiter=400))
        x = np.clip(r.x, [b[0] for b in bnds], [b[1] for b in bnds])
        if trim is not None and not inside_trim(trim, np.array([x[0]]), np.array([x[1]]))[0]:
            continue
        dd = S.pts([x[0]], [x[1]])
        km = principal(dd)[0][0]
        c_, t_, r_, e_ = chord_info(dd["S"])
        refined.append(dict(u=float(x[0]), v=float(x[1]), radius_mm=float(1 / km), ratio=float((1 / km) / r_[0]),
                            P=dd["S"][0], chi=float(c_[0]), t_att=float(t_[0]), rule=float(r_[0]), eta=float(e_[0])))
    def rep(mask, label):
        if not np.any(mask):
            return None
        i, j = worst(mask, ratio)
        best = dict(radius_mm=float(R[i, j]), ratio_to_rule=float(ratio[i, j]), at=loc(P[i, j]),
                    chord_fraction=float(chi[i, j]), attachment_thickness_mm=float(t[i, j]), rule_radius_mm=float(req[i, j]),
                    span_beyond_attachment_mm=float(eta[i, j]), uv=[float(u[i]), float(v[j])])
        for rf in refined:
            inzone = (rf["eta"] < cm.ETA_F and rf["chi"] > 0.99) == (label == "corner")
            if inzone and rf["ratio"] < best["ratio_to_rule"]:
                best = dict(radius_mm=rf["radius_mm"], ratio_to_rule=rf["ratio"], at=loc(rf["P"]), chord_fraction=rf["chi"],
                            attachment_thickness_mm=rf["t_att"], rule_radius_mm=rf["rule"],
                            span_beyond_attachment_mm=rf["eta"], uv=[rf["u"], rf["v"]], refined=True)
        i2, j2 = worst(mask, R)
        best["smallest_radius_in_zone_mm"] = float(min(R[i2, j2], min([rf["radius_mm"] for rf in refined
                                                                         if ((rf["eta"] < cm.ETA_F and rf["chi"] > 0.99) == (label == "corner"))] + [np.inf])))
        best["smallest_radius_at"] = loc(P[i2, j2])
        best["smallest_radius_chord_fraction"] = float(chi[i2, j2])
        best["smallest_radius_attachment_thickness_mm"] = float(t[i2, j2])
        best["points_below_rule"] = int(np.sum(mask & (ratio < 1.0)))
        best["points"] = int(mask.sum())
        return best
    res["outside_corner_zone"] = rep(~corner, "main")
    res["te_corner_zone"] = rep(corner, "corner")
    thick = (~corner) & (t >= cm.RULE_T)
    thin = (~corner) & (t < cm.RULE_T)
    res["thick_zone_min_radius_mm"] = float(R[thick].min()) if np.any(thick) else None
    res["thin_zone_min_ratio_to_half_thickness"] = float(ratio[thin].min()) if np.any(thin) else None
    return res, dict(P=P, R=R, ratio=ratio, chi=chi, corner=corner, u=u, v=v)


def section_curvature(S, u, v, n_pl, orient):
    d = S.pts(u, v)
    E, F, G, L, M, Nn, N, _ = fundamental(d)
    N = N * orient
    L, M, Nn = L * orient, M * orient, Nn * orient
    t = np.cross(N, n_pl[None, :]); t /= np.linalg.norm(t, axis=-1, keepdims=True)
    A = np.stack([d["Su"], d["Sv"]], -1)
    ATA = np.einsum("nki,nkj->nij", A, A)
    a = np.linalg.solve(ATA, np.einsum("nki,nk->ni", A, t)[..., None])[..., 0]
    kn = L * a[:, 0] ** 2 + 2 * M * a[:, 0] * a[:, 1] + Nn * a[:, 1] ** 2
    nc = np.cross(n_pl[None, :], t)
    cth = np.sum(nc * N, -1)
    bend = -kn / np.where(np.abs(cth) > 1e-9, np.abs(cth), np.nan)     # > 0 convex
    return d["S"], bend


def span_iso(S, target):
    """Which parameter runs along the span and its value where |Y| = target."""
    c = S.grid([S.u0, S.u1], [S.v0, S.v1])["S"]
    dYu = abs(abs(c[1, 0, 1]) - abs(c[0, 0, 1])); dYv = abs(abs(c[0, 1, 1]) - abs(c[0, 0, 1]))
    if dYv > dYu:
        val = S.v0 if abs(abs(c[0, 0, 1]) - target) < abs(abs(c[0, 1, 1]) - target) else S.v1
        return "v", val
    val = S.u0 if abs(abs(c[0, 0, 1]) - target) < abs(abs(c[1, 0, 1]) - target) else S.u1
    return "u", val


def project(S, P, guess_grid=(300, 60)):
    from scipy.spatial import cKDTree
    ug = np.linspace(S.u0, S.u1, guess_grid[0]); vg = np.linspace(S.v0, S.v1, guess_grid[1])
    G = S.grid(ug, vg)["S"]
    tree = cKDTree(G.reshape(-1, 3))
    dist, idx = tree.query(P)
    i, j = np.unravel_index(idx, G.shape[:2])
    u, v = ug[i].copy(), vg[j].copy()
    for _ in range(40):
        d = S.pts(u, v)
        r = d["S"] - P
        Ju = d["Su"]; Jv = d["Sv"]
        a11 = np.sum(Ju * Ju, -1) + np.sum(r * d["Suu"], -1)
        a12 = np.sum(Ju * Jv, -1) + np.sum(r * d["Suv"], -1)
        a22 = np.sum(Jv * Jv, -1) + np.sum(r * d["Svv"], -1)
        b1 = np.sum(r * Ju, -1); b2 = np.sum(r * Jv, -1)
        det = a11 * a22 - a12 * a12
        det = np.where(np.abs(det) > 1e-300, det, 1e-300)
        du = (a22 * b1 - a12 * b2) / det; dv = (a11 * b2 - a12 * b1) / det
        u = np.clip(u - du, S.u0, S.u1); v = np.clip(v - dv, S.v0, S.v1)
    d = S.pts(u, v)
    return u, v, np.linalg.norm(d["S"] - P, axis=-1)


def run_section(shape, cands, P0, n, keep, center, nu_dir, npts=6000):
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    sec = BRepAlgoAPI_Section(shape, gp_Pln(gp_Pnt(*map(float, P0)), gp_Dir(*map(float, n))), False)
    sec.Approximation(True)
    sec.Build()
    pts = []
    ex = TopExp_Explorer(sec.Shape(), TopAbs_EDGE)
    while ex.More():
        e = TopoDS.Edge_s(ex.Current())
        ca = BRepAdaptor_Curve(e)
        u0, u1 = ca.FirstParameter(), ca.LastParameter()
        for tt in np.linspace(u0, u1, 400):
            q = ca.Value(tt)
            pts.append((q.X(), q.Y(), q.Z()))
        ex.Next()
    pts = np.array(pts)
    pts = pts[keep(pts)]
    if len(pts) == 0:
        return None
    # assign to faces and project
    best = np.full(len(pts), np.inf); fid = np.zeros(len(pts), int); U = np.zeros(len(pts)); V = np.zeros(len(pts))
    for k, (S, orient, trim) in enumerate(cands):
        u, v, dd = project(S, pts)
        if trim is not None:
            dd = np.where(inside_trim(trim, u, v), dd, np.inf)
        sel = dd < best
        best[sel] = dd[sel]; fid[sel] = k; U[sel] = u[sel]; V[sel] = v[sel]
    # densify: order by angle around the centre, then resample each face's run on the plane
    mvec = np.cross(n, nu_dir); mvec /= np.linalg.norm(mvec)
    rel = pts - center
    th = np.arctan2(rel @ mvec, -(rel @ nu_dir))
    o = np.argsort(th)
    pts, fid, U, V, th = pts[o], fid[o], U[o], V[o], th[o]
    # refine every point exactly onto the plane by Newton along the surface gradient of F
    for k, (S, orient, trim) in enumerate(cands):
        s = fid == k
        if not np.any(s):
            continue
        u, v = U[s], V[s]
        for _ in range(20):
            d = S.pts(u, v)
            F = (d["S"] - P0) @ n
            gu = d["Su"] @ n; gv = d["Sv"] @ n
            g2 = gu * gu + gv * gv
            u = np.clip(u - F * gu / np.maximum(g2, 1e-300), S.u0, S.u1)
            v = np.clip(v - F * gv / np.maximum(g2, 1e-300), S.v0, S.v1)
        U[s], V[s] = u, v
    P = np.zeros_like(pts); bend = np.zeros(len(pts))
    for k, (S, orient, trim) in enumerate(cands):
        s = fid == k
        if np.any(s):
            P[s], bend[s] = section_curvature(S, U[s], V[s], n, orient)
    return dict(P=P, bend=bend, face=fid, theta=th, proj_err=best)


def band_polyline(face):
    """Band-edge pcurve of a trimmed flat face as (s, Y) samples sorted by s (read from the STEP topology)."""
    from OCP.BRepAdaptor import BRepAdaptor_Curve2d
    best = None
    ex = TopExp_Explorer(face, TopAbs_EDGE)
    while ex.More():
        e = TopoDS.Edge_s(ex.Current())
        c2 = BRepAdaptor_Curve2d(e, face)
        t0, t1 = c2.FirstParameter(), c2.LastParameter()
        pts = np.array([[c2.Value(t).X(), c2.Value(t).Y()] for t in np.linspace(t0, t1, 20001)])
        if np.ptp(pts[:, 1]) > 1.0:          # the attachment edge is the straight iso-line Y = 440
            best = pts
        ex.Next()
    o = np.argsort(best[:, 0])
    return best[o]


def inside_trim(poly, s, Y):
    return Y <= np.interp(s, poly[:, 0], poly[:, 1]) + 1e-9


def outward_orient(face, S):
    from OCP.TopAbs import TopAbs_REVERSED
    return -1.0 if face.Orientation() == TopAbs_REVERSED else 1.0


def main():
    t0 = time.time()
    out = {"model": "twisted/blunt: parent + rounded caps at both tips (right cap mirrored to the left)",
           "units": "mm", "step": {"file": os.path.relpath(STEP, cm.ROOT), "sha256": sha256(STEP),
                                   "bytes": os.path.getsize(STEP)}}
    shape = read_step(STEP)
    out["reimport"] = topology(shape)
    print("topology", out["reimport"]); sys.stdout.flush()
    info, pf = classify(shape)
    parents = [x for x in info if x.get("parent")]
    caps = [x for x in info if "surf" in x and not x.get("parent")]
    out["parent_preserved"] = {
        "parent_bspline_faces_in_parent_step": len(pf),
        "matched_identically_in_wing_step": len(parents),
        "max_pole_difference_mm": max(x["parent"][1] for x in parents) if parents else None,
        "planar_end_faces_remaining": int(sum(1 for x in info if x.get("kind") == "non-bspline")),
        "note": "each parent B-spline face's poles compared with parent.step; the two planar attachment faces are replaced by the caps"}
    # identify faces
    def meanY(x): return float(x["surf"].P[..., 1].mean())
    out["parent_preserved"]["cap_faces"] = int(len(caps))
    def meanZ(x): return float(x["surf"].P[..., 2].mean())
    tips = {}
    for side, sg in (("right", 1), ("left", -1)):
        cf = [x for x in caps if np.sign(meanY(x)) == sg]
        lip = [x for x in cf if x["surf"].P.shape[1] == 61][0]
        fl = sorted([x for x in cf if x["surf"].P.shape[0] == 511], key=meanZ)
        for x in fl:
            x["trim"] = band_polyline(x["face"])
        pp = [x for x in parents if np.sign(meanY(x)) == sg]
        te = [x for x in pp if min(x["surf"].P.shape[:2]) <= 3]
        sk = sorted([x for x in pp if min(x["surf"].P.shape[:2]) > 3], key=meanZ)
        tips[side] = dict(lip=lip, up=fl[1], lo=fl[0], te=te[0] if te else None, pup=sk[1], plo=sk[0])
    print("faces identified", time.time() - t0); sys.stdout.flush()

    # ---------------- locked dimensions and the outline
    bb = out["reimport"]["bounding_box_mm"]
    rimg = np.linspace(0, 1, 8001)
    rim = cm.pa.tip_outline(rimg)
    T = tips["right"]
    Lip = T["lip"]["surf"]
    u, v, dist = project(Lip, rim, guess_grid=(800, 120))
    ii = int(np.argmax(dist))
    # plan-view silhouette of the cap against the locked outline (right tip)
    ug = param_grid(Lip.tu, Lip.du, 3); vg = param_grid(Lip.tv, Lip.dv, 3)
    Pl = Lip.grid(ug, vg)["S"].reshape(-1, 3)
    from scipy.spatial import cKDTree
    rimd = cm.pa.tip_outline(np.linspace(0, 1, 200001))
    tree = cKDTree(rimd[:, :2])
    dd, idx = tree.query(Pl[:, :2])
    tng = np.gradient(rimd[:, :2], axis=0); tng /= np.linalg.norm(tng, axis=-1, keepdims=True)
    nout = np.stack([-tng[:, 1], tng[:, 0]], -1)
    signed = np.sum((Pl[:, :2] - rimd[idx, :2]) * nout[idx], -1)
    io = int(np.argmax(signed))
    pio = float(np.linspace(0, 1, 200001)[idx[io]])
    rp = np.linspace(0, 1, 200001)[idx]
    le_half = rp < 0.5
    # reference: how far the parent's own section at the attachment sticks out ahead of the rim's start in plan
    sg_ = np.linspace(0, 1, 20001)
    secU = cm.pa.skin(sg_, cm.CUT, +1); secL = cm.pa.skin(sg_, cm.CUT, -1)
    parent_le_ahead = float(rim[0, 0] - min(secU[:, 0].min(), secL[:, 0].min()))
    out["locked_dimensions"] = {
        "plan_silhouette_outside_outline_leading_edge_half_max_mm": float(signed[le_half].max()),
        "plan_silhouette_outside_outline_trailing_edge_half_max_mm": float(signed[~le_half].max()),
        "parent_attachment_section_ahead_of_rim_start_in_plan_mm": parent_le_ahead,
        "span_extent_Y_mm": [bb[2], bb[3]], "finished_projected_span_mm": bb[3] - bb[2],
        "spec_projected_span_mm": 1000.0,
        "rim_to_cap_surface_max_mm": float(dist.max()), "rim_worst_p": float(rimg[ii]),
        "plan_silhouette_outside_locked_outline_max_mm": float(signed.max()), "silhouette_worst_at": loc(Pl[io]),
        "silhouette_worst_rim_p": pio,
        "silhouette_note": "the lip is symmetric about the parent's mid surface; near the trailing edge the locked rim lies up to 0.091 mm below that mid surface, so the lip apex (plan silhouette) sits slightly outboard of the rim"}
    print("outline", out["locked_dimensions"]); sys.stdout.flush()

    # ---------------- joins
    joins = {}
    for side in ("right", "left"):
        T = tips[side]
        res = {}
        for skin, capf in (("upper", "up"), ("lower", "lo")):
            Ps = T["p" + skin[:2]]["surf"]
            dirn, val = span_iso(Ps, cm.CUT)
            n = 6001
            uu = Ps.u0 + (Ps.u1 - Ps.u0) * (0.5 - 0.5 * np.cos(np.pi * np.linspace(0, 1, n))) if dirn == "v" else None
            if dirn == "v":
                ua, va = uu, np.full(n, val)
            else:
                va = Ps.v0 + (Ps.v1 - Ps.v0) * (0.5 - 0.5 * np.cos(np.pi * np.linspace(0, 1, n))); ua = np.full(n, val)
            Pa = Ps.pts(ua, va)["S"]
            Cf = T[capf]["surf"]
            cdir, cval = span_iso(Cf, cm.CUT)
            ub = match_on_iso(Cf, cval, Pa)
            j = join_check(Ps, ua, va, Cf, ub, np.full(n, cval), tangent_from="u" if dirn == "v" else "v")
            LEp = cm.pa.skin(0.0, cm.CUT, +1) * np.array([1, 1 if side == "right" else -1, 1])
            TEp = cm.pa.skin(1.0, cm.CUT, 1 if skin == "upper" else -1) * np.array([1, 1 if side == "right" else -1, 1])
            dend = np.minimum(np.linalg.norm(Pa - LEp, axis=-1), np.linalg.norm(Pa - TEp, axis=-1))
            excl = dend < 0.01
            res[skin] = summarize_join(j, skin, excl)
            res[skin]["all_samples_including_end_points"] = summarize_join(j, skin)
            res[skin]["exclusion"] = "samples within 0.01 mm of the leading-edge point or the trailing-edge corner (the corners of the cap face, where the band edge meets the attachment)"
            if side == "right":
                res[skin]["zebra"] = zebra_join(j, skin, excl=excl)
        # trailing-edge segment: parent TE face vs lip g = 1
        if T["te"] is not None:
            Pt = T["te"]["surf"]
            dirn, val = span_iso(Pt, cm.CUT)
            n = 801
            if dirn == "v":
                ua = np.linspace(Pt.u0, Pt.u1, n); va = np.full(n, val)
            else:
                va = np.linspace(Pt.v0, Pt.v1, n); ua = np.full(n, val)
            Pa = Pt.pts(ua, va)["S"]
            Lp = T["lip"]["surf"]
            # match along the lip's u = 1 boundary (v runs across the segment)
            vs_ = np.linspace(Lp.v0, Lp.v1, 20001)
            C = Lp.grid([Lp.u1], vs_)["S"][0]
            from scipy.spatial import cKDTree as KD
            _, idx = KD(C).query(Pa)
            j = join_check(Pt, ua, va, Lp, np.full(n, Lp.u1), vs_[idx], tangent_from="u" if dirn == "v" else "v")
            corners = np.minimum(np.linalg.norm(Pa - Pa[0], axis=-1), np.linalg.norm(Pa - Pa[-1], axis=-1)) < 0.01
            res["trailing_edge_segment"] = summarize_join(j, "te", corners)
            res["trailing_edge_segment"]["exclusion"] = "within 0.01 mm of the two trailing-edge corners (the parent's own creases end there)"
            if side == "right":
                res["trailing_edge_segment"]["zebra"] = zebra_join(j, "te", excl=corners)
        # internal lip / flat boundaries (band edges)
        Lp = T["lip"]["surf"]
        n = 4001
        ug = np.linspace(Lp.u0, Lp.u1, n)[1:-1]
        for nm, lv, capf in (("band_edge_upper", Lp.v1, "up"), ("band_edge_lower", Lp.v0, "lo")):
            Cf = T[capf]["surf"]
            Pb = Lp.pts(ug, np.full(len(ug), lv))["S"]
            uf, vf, dproj = project(Cf, Pb, guess_grid=(1200, 80))
            j = join_check(Lp, ug, np.full(len(ug), lv), Cf, uf, vf, tangent_from="u")
            pr = cm.p_of_g(ug)
            etaj = np.abs(j["P"][:, 1]) - cm.CUT
            corner_j = (etaj < cm.ETA_F) & (pr > 0.5)
            res[nm] = summarize_join(j, nm, corner_j)
            res[nm]["exclusion"] = "the trailing-edge corner zone (span < %.2f mm), where the band closes onto the blunt TE's sharp corner" % cm.ETA_F
            res[nm]["inside_te_corner_zone"] = summarize_join(j, nm, ~corner_j)
            # curvature jump split by design zone
            d22 = np.abs(j["ka"][2] - j["kb"][2])
            def bin_(m):
                m = m & ~corner_j
                if not np.any(m):
                    return None
                i = int(np.argmax(np.where(m, d22, -1)))
                return {"max_per_mm": float(d22[i]), "at_p": float(pr[i]), "at": loc(j["P"][i]),
                        "skin_curvature_there_per_mm": float(abs(j["kb"][2][i])), "lip_curvature_there_per_mm": float(abs(j["ka"][2][i]))}
            res[nm]["cross_join_curvature_jump_by_rim_parameter"] = {
                "p<0.01 (leading-edge corner, band < 0.2 mm)": bin_(pr < 0.01),
                "0.01<=p<0.78 (shifted nose / superelliptic m=3: designed curvature-continuous)": bin_((pr >= 0.01) & (pr < 0.78)),
                "0.78<=p<0.925 (blend to the inscribed circle)": bin_((pr >= 0.78) & (pr < 0.925)),
                "p>=0.925 outside the corner zone (inscribed circle = rule-mandated half-round: curvature jump expected)": bin_(pr >= 0.925)}
            if side == "right":
                res[nm]["zebra"] = zebra_join(j, nm, excl=corner_j)
        joins[side] = res
        print("joins", side, time.time() - t0); sys.stdout.flush()
    out["joins"] = joins

    # ---------------- curvature maps (larger principal curvature) on every cap face
    cmaps = {}
    keep = {}
    for side in ("right", "left"):
        T = tips[side]
        res = {}
        for nm in ("lip", "up", "lo"):
            S = T[nm]["surf"]
            lo_ = 1e-12
            hi_ = S.u1
            trim = T[nm].get("trim")
            r, k = curvature_map(S, nm, lo_, hi_, nper=4 if side == "right" else 2, refine=12 if side == "right" else 4,
                                 trim=trim)
            # how close the grid reaches the collapsed points
            first = S.grid([k["u"][0]], [0.5 * (S.v0 + S.v1)])["S"][0, 0]
            LEp = cm.pa.skin(0.0, cm.CUT, +1) * np.array([1, 1 if side == "right" else -1, 1])
            if nm == "lip":
                r["grid_closest_to_leading_edge_pole_mm"] = float(np.linalg.norm(k["P"][0] - LEp, axis=-1).max())
                r["note"] = "the lip's g = 0 row collapses to the leading-edge point (derivatives undefined): excluded"
            else:
                r["note"] = "trimmed continued-skin face: only grid points inside the band-edge trim are evaluated"
            res[nm] = r
            if side == "right":
                keep[nm] = k
            print("curvature", side, nm, json.dumps(r["outside_corner_zone"])[:400], time.time() - t0); sys.stdout.flush()
        cmaps[side] = res
    out["curvature_map"] = cmaps
    # overall summary
    allmain = [cmaps[s][f]["outside_corner_zone"] for s in cmaps for f in cmaps[s] if cmaps[s][f]["outside_corner_zone"]]
    allcorner = [cmaps[s][f]["te_corner_zone"] for s in cmaps for f in cmaps[s] if cmaps[s][f]["te_corner_zone"]]
    w = min(allmain, key=lambda x: x["ratio_to_rule"])
    wr = min(allmain, key=lambda x: x["smallest_radius_in_zone_mm"])
    out["printable_rule_summary"] = {
        "rule": "radius >= 0.4 mm where the attachment section is >= 0.8 mm thick at that chord fraction; else a half-round of radius t/2",
        "outside_trailing_edge_corner_zone": {"min_ratio_radius_to_rule": w["ratio_to_rule"], "at": w["at"], "radius_mm": w["radius_mm"],
                                              "chord_fraction": w["chord_fraction"], "attachment_thickness_mm": w["attachment_thickness_mm"],
                                              "rule_radius_mm": w["rule_radius_mm"],
                                              "smallest_radius_mm": wr["smallest_radius_in_zone_mm"], "smallest_radius_at": wr["smallest_radius_at"],
                                              "smallest_radius_chord_fraction": wr["smallest_radius_chord_fraction"],
                                              "smallest_radius_attachment_thickness_mm": wr["smallest_radius_attachment_thickness_mm"]},
        "te_corner_zone": {"definition": "span < %.2f mm beyond the attachment and chord fraction > 0.99: the blunt TE's sharp corners open into the half-round" % cm.ETA_F,
                           "min_radius_mm": min(x["smallest_radius_in_zone_mm"] for x in allcorner) if allcorner else None,
                           "min_ratio": min(x["ratio_to_rule"] for x in allcorner) if allcorner else None}}
    print("rule summary", out["printable_rule_summary"]); sys.stdout.flush()

    # ---------------- sections through the cap (perpendicular to the rim) and spanwise chord-fraction cuts
    T = tips["right"]
    cands = [(T[k]["surf"], outward_orient(T[k]["face"], T[k]["surf"]), T[k].get("trim")) for k in ("lip", "up", "lo", "pup", "plo")]
    if T["te"] is not None:
        cands.append((T["te"]["surf"], outward_orient(T["te"]["face"], T["te"]["surf"]), None))
    dzs = cm.design(np.array(SECTION_P))
    secs = []
    for i, p in enumerate(SECTION_P):
        Rp, dRp = cm.rim(np.array([p]))
        Rp = Rp[0]; Tn = dRp[0] / np.linalg.norm(dRp[0])
        nu_ = dzs["nu"][i]; nu_ = nu_ - (nu_ @ Tn) * Tn; nu_ /= np.linalg.norm(nu_)
        W = max(3.0 * dzs["bend"][i], 6.0)
        center = Rp + min(W / 2, 3.0) * nu_
        keep_fn = lambda P_, Rp=Rp, nu_=nu_, W=W: (P_[:, 1] > 300) & (((P_ - Rp) @ nu_) < W + 15.0) & (((P_ - Rp) @ nu_) > -2.0)
        s = run_section(shape, cands, Rp, Tn, keep_fn, center, nu_)
        if s is None:
            continue
        capm = s["P"][:, 1] >= cm.CUT - 1e-6
        bend = s["bend"]
        rad = 1 / np.abs(bend)
        chi, t, req, eta = chord_info(s["P"])
        ok = capm & np.isfinite(rad)
        imin = int(np.argmin(np.where(ok, rad, np.inf)))
        irim = int(np.argmin(np.linalg.norm(s["P"] - Rp, axis=-1)))
        iap = int(np.argmax(np.where(capm, -(s["P"] - Rp) @ nu_, -np.inf)))
        sgn = np.sign(np.where(np.abs(bend) > 1e-3, bend, 0))
        nz = np.where(sgn != 0)[0]
        rev = [loc(s["P"][nz[k]]) for k in range(1, len(nz)) if sgn[nz[k]] != sgn[nz[k - 1]]]
        conc = np.where(np.isfinite(bend), bend, 0.0)
        iconc = int(np.argmin(conc))
        chir, tr, rr, _ = chord_info(Rp[None, :])
        secs.append({"kind": "perpendicular to the rim", "rim_p": p, "plane_point": loc(Rp), "plane_normal": loc(Tn),
                     "points": int(len(rad)), "cap_points": int(capm.sum()),
                     "includes_parent_skin": bool(np.any(~capm)),
                     "min_radius_on_cap_mm": float(rad[imin]), "min_radius_at": loc(s["P"][imin]),
                     "min_radius_chord_fraction": float(chi[imin]), "rule_at_min_point_mm": float(req[imin]),
                     "radius_at_rim_point_mm": float(rad[irim]), "radius_at_lip_apex_mm": float(rad[iap]),
                     "attachment_thickness_at_rim_chord_fraction_mm": float(tr[0]), "rule_radius_mm": float(rr[0]),
                     "passes_rule_section_thickness": bool(rad[imin] >= rr[0] * 0.999),
                     "passes_rule_pointwise": bool(np.all((rad[ok] / req[ok]) >= 0.999)),
                     "min_ratio_pointwise": float(np.min(rad[ok] / req[ok])),
                     "bending_reversals": len(rev), "reversal_locations": rev[:10],
                     "strongest_concave_curvature_per_mm": float(min(conc[iconc], 0.0)),
                     "strongest_concave_at": loc(s["P"][iconc]),
                     "strongest_concave_on_cap": bool(s["P"][iconc, 1] >= cm.CUT)})
        print("section p", p, secs[-1]["min_radius_on_cap_mm"], secs[-1]["radius_at_lip_apex_mm"], len(rev)); sys.stdout.flush()
    # chord-fraction cuts following the sweep
    for chi0 in SECTION_CHI:
        o, cd, ud, c = cm.section_frame_v(np.array([cm.CUT, cm.CUT + 1.0]))
        A = o[0] + chi0 * c[0] * cd[0]; Bq = o[1] + chi0 * c[1] * cd[1]
        dsp = Bq - A; dsp /= np.linalg.norm(dsp)
        nrm = np.cross(dsp, ud[0]); nrm /= np.linalg.norm(nrm)
        nu_ = -dsp
        center = A + 30 * dsp
        keep_fn = lambda P_: (P_[:, 1] > 420)
        s = run_section(shape, cands, A, nrm, keep_fn, A + 70 * dsp, nu_)
        if s is None:
            continue
        capm = s["P"][:, 1] >= cm.CUT - 1e-6
        rad = 1 / np.abs(s["bend"])
        chi, t, req, eta = chord_info(s["P"])
        ok = capm & np.isfinite(rad)
        if not np.any(ok):
            continue
        imin = int(np.argmin(np.where(ok, rad, np.inf)))
        bend = s["bend"]
        sgn = np.sign(np.where(np.abs(bend) > 1e-3, bend, 0)); nz = np.where(sgn != 0)[0]
        rev = [loc(s["P"][nz[k]]) for k in range(1, len(nz)) if sgn[nz[k]] != sgn[nz[k - 1]]]
        secs.append({"kind": "constant chord fraction, following the sweep", "chord_fraction": chi0,
                     "points": int(len(rad)), "cap_points": int(capm.sum()),
                     "min_radius_on_cap_mm": float(rad[imin]), "min_radius_at": loc(s["P"][imin]),
                     "attachment_thickness_mm": float(cm.t_att(chi0)), "rule_radius_mm": float(cm.rule_radius(chi0)),
                     "min_ratio_pointwise": float(np.min(rad[ok] / req[ok])),
                     "bending_reversals": len(rev), "reversal_locations": rev[:10],
                     "note": "oblique to the rim where the rim is not perpendicular to the cut: by Meusnier such cuts show a tighter curve than the surface's true principal radius"})
        print("chord cut", chi0, secs[-1]["min_radius_on_cap_mm"], len(rev)); sys.stdout.flush()
    out["sections"] = secs
    out["seconds"] = time.time() - t0
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print("wrote", OUT, "%.1f s" % (time.time() - t0))


SECTION_P = [0.001, 0.005, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.92, 0.94, 0.96, 0.98,
             0.99, 0.995, 0.998, 0.9993, 0.9997]
SECTION_CHI = [0.01, 0.1, 0.3, 0.5, 0.7, 0.85, 0.95, 0.98]

if __name__ == "__main__":
    main()
