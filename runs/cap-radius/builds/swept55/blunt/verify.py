"""Verify the reimported wing.step against the skill 'close-rounded-wingtip' and write readback.json.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/swept55/blunt/verify.py
Everything below is measured on the STEP read back from disk (OCCT), not on the build objects.
"""
import hashlib
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cap_model as cm  # noqa: E402
import radius_map as rm  # noqa: E402

from OCP.BRep import BRep_Tool  # noqa: E402
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Curve2d, BRepAdaptor_Surface  # noqa: E402
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section  # noqa: E402
from OCP.BRepBndLib import BRepBndLib  # noqa: E402
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace  # noqa: E402
from OCP.BRepCheck import BRepCheck_Analyzer  # noqa: E402
from OCP.BRepGProp import BRepGProp  # noqa: E402
from OCP.Bnd import Bnd_Box  # noqa: E402
from OCP.GProp import GProp_GProps  # noqa: E402
from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf  # noqa: E402
from OCP.GeomAbs import GeomAbs_Plane, GeomAbs_BSplineSurface  # noqa: E402
from OCP.GeomLProp import GeomLProp_SLProps  # noqa: E402
from OCP.IFSelect import IFSelect_RetDone  # noqa: E402
from OCP.STEPControl import STEPControl_Reader  # noqa: E402
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID  # noqa: E402
from OCP.TopExp import TopExp, TopExp_Explorer  # noqa: E402
from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape  # noqa: E402
from OCP.TopoDS import TopoDS  # noqa: E402
from OCP.gp import gp_Dir, gp_Pln, gp_Pnt  # noqa: E402

STEP = os.path.join(HERE, "wing.step")
OUT = os.path.join(HERE, "readback.json")


def read(path):
    r = STEPControl_Reader()
    assert r.ReadFile(path) == IFSelect_RetDone
    r.TransferRoots()
    return r.OneShape()


def explore(shape, kind):
    out = []
    ex = TopExp_Explorer(shape, kind)
    while ex.More():
        out.append(ex.Current())
        ex.Next()
    return out


def bbox(shape):
    b = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape, b, False, False)
    return list(b.Get())


def srf_of(face):
    return BRep_Tool.Surface_s(face)


# ----------------------------------------------------------------------------- topology
def topology(shape):
    solids = explore(shape, TopAbs_SOLID)
    shells = explore(shape, TopAbs_SHELL)
    faces = [TopoDS.Face_s(f) for f in explore(shape, TopAbs_FACE)]
    emap = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(shape, TopAbs_EDGE, TopAbs_FACE, emap)
    counts = {}
    for i in range(1, emap.Extent() + 1):
        e = TopoDS.Edge_s(emap.FindKey(i))
        if BRep_Tool.Degenerated_s(e):
            counts["degenerated"] = counts.get("degenerated", 0) + 1
            continue
        n = emap.FindFromIndex(i).Size()
        counts[n] = counts.get(n, 0) + 1
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props, 1e-9)
    return dict(solids=len(solids), shells=len(shells), faces=len(faces),
                brepcheck_valid=bool(BRepCheck_Analyzer(shape).IsValid()),
                edges_by_face_count={str(k): v for k, v in counts.items()},
                closed=all(k in (2, "degenerated") for k in counts),
                volume_mm3=props.Mass(), bounding_box_mm=[round(x, 7) for x in bbox(shape)]), faces


def classify(faces):
    out = {"parent": [], "cap": {}}
    for f in faces:
        b = bbox(f)
        a = BRepAdaptor_Surface(f)
        if b[4] > cm.ATT + 1e-3 and b[1] > cm.ATT - 1e-3:
            side = "R"
        elif b[1] < -cm.ATT - 1e-3 and b[4] < -cm.ATT + 1e-3:
            side = "L"
        else:
            out["parent"].append(f)
            continue
        if a.GetType() == GeomAbs_Plane:
            key = side + "_disc"
        else:
            key = side + ("_upper" if b[5] > 1e-3 else "_lower")
        out["cap"][key] = f
    return out


# ----------------------------------------------------------------------------- parent preserved
def parent_preserved(parent_faces):
    """Retained parent faces vs parent.step faces (poles) and vs parent.py (sampled)."""
    ref = [TopoDS.Face_s(f) for f in explore(read(os.path.join(cm.PARENT_DIR, "parent.step")), TopAbs_FACE)]
    ref = [f for f in ref if BRepAdaptor_Surface(f).GetType() != GeomAbs_Plane]

    def poles(f):
        s = BRepAdaptor_Surface(f).BSpline()
        return np.array([[(s.Pole(i, j).X(), s.Pole(i, j).Y(), s.Pole(i, j).Z())
                          for j in range(1, s.NbVPoles() + 1)] for i in range(1, s.NbUPoles() + 1)])

    worst = 0.0
    matched = 0
    for f in parent_faces:
        P = poles(f)
        best = min((np.abs(P - poles(g)).max() for g in ref if poles(g).shape == P.shape), default=np.inf)
        worst = max(worst, best)
        matched += best < 1e-9
    # skins vs parent.py
    dev = 0.0
    for f in parent_faces:
        a = BRepAdaptor_Surface(f)
        s = srf_of(f)
        b = bbox(f)
        if b[5] - b[2] < 0.5:      # TE face
            continue
        side = 1 if b[5] > 1e-3 else -1
        sgn = 1 if b[4] > 1e-3 else -1
        U = np.linspace(a.FirstUParameter(), a.LastUParameter(), 41)
        V = np.linspace(a.FirstVParameter(), a.LastVParameter(), 21)
        for v in V:
            ref_pts = cm.parent.skin(U, sgn * v, side)
            for u, q in zip(U, ref_pts):
                p = s.Value(float(u), float(v))
                dev = max(dev, float(np.linalg.norm(np.array([p.X(), p.Y(), p.Z()]) - q)))
    return dict(faces_retained=len(parent_faces), faces_identical_to_parent_step=int(matched),
                max_pole_difference_mm=float(worst), max_distance_to_parent_py_mm=dev)


# ----------------------------------------------------------------------------- lip / outline
def lip_points(face, n=4001):
    """Points of the face's z = 0 boundary (the lip), by edge sampling."""
    pts = []
    for e in explore(face, TopAbs_EDGE):
        e = TopoDS.Edge_s(e)
        if BRep_Tool.Degenerated_s(e):
            continue
        c = BRepAdaptor_Curve(e)
        t = np.linspace(c.FirstParameter(), c.LastParameter(), n)
        P = np.array([[c.Value(float(x)).X(), c.Value(float(x)).Y(), c.Value(float(x)).Z()] for x in t])
        if np.abs(P[:, 2]).max() < 1e-6 and np.abs(P[:, 1]).max() > cm.ATT + 1:
            pts.append(P)
    return np.vstack(pts)


def outline_deviation(P):
    """Plan distance of lip points to the locked outline (both tips via |Y|)."""
    pp = np.linspace(0, 1, 200001)
    O = cm.parent.tip_outline(pp)
    from scipy.spatial import cKDTree
    tree = cKDTree(O[:, :2])
    Q = np.c_[P[:, 0], np.abs(P[:, 1])]
    d, i = tree.query(Q)
    # refine with the segment to the neighbour
    out = []
    for q, k in zip(Q, i):
        best = np.inf
        for k2 in (k - 1, k + 1):
            if 0 <= k2 < len(O):
                a, b = O[k, :2], O[k2, :2]
                t = np.clip(np.dot(q - a, b - a) / np.dot(b - a, b - a), 0, 1)
                best = min(best, np.linalg.norm(q - (a + t * (b - a))))
        out.append(best)
    out = np.array(out)
    j = int(np.argmax(out))
    return dict(max_plan_distance_mm=float(out.max()), at=P[j].round(5).tolist(),
                p_at=float(pp[i[j]]), max_abs_z_mm=float(np.abs(P[:, 2]).max()),
                n_points=int(len(P)))


# ----------------------------------------------------------------------------- joins
def props_at(srf, u, v):
    pr = GeomLProp_SLProps(srf, float(u), float(v), 2, 1e-12)
    p = pr.Value()
    out = dict(p=np.array([p.X(), p.Y(), p.Z()]))
    if pr.IsNormalDefined():
        n = pr.Normal(); out["n"] = np.array([n.X(), n.Y(), n.Z()])
    if pr.IsCurvatureDefined():
        d1, d2 = gp_Dir(), gp_Dir()
        pr.CurvatureDirections(d1, d2)
        out["k"] = (pr.MaxCurvature(), pr.MinCurvature())
        out["d"] = (np.array([d1.X(), d1.Y(), d1.Z()]), np.array([d2.X(), d2.Y(), d2.Z()]))
    return out


def project(srf, p):
    pj = GeomAPI_ProjectPointOnSurf(gp_Pnt(*map(float, p)), srf)
    if pj.NbPoints() == 0:
        return None
    u, v = pj.LowerDistanceParameters()
    return u, v, pj.LowerDistance()


def second_form(pr, n_ref, t):
    """Normal curvature along t, across (n x t) and the twist, w.r.t. the reference normal n_ref."""
    if "k" not in pr:
        return None
    k1, k2 = pr["k"]
    d1, d2 = pr["d"]
    sgn = 1.0 if np.dot(pr["n"], n_ref) > 0 else -1.0
    b = np.cross(n_ref, t)
    def kn(x):
        return sgn * (k1 * np.dot(x, d1) ** 2 + k2 * np.dot(x, d2) ** 2)
    tw = sgn * ((k1 - k2) * np.dot(t, d1) * np.dot(b, d1))
    return kn(t), kn(b), tw


def shared_edges(fa, fb):
    ea = {TopoDS.Edge_s(e).__hash__(): TopoDS.Edge_s(e) for e in explore(fa, TopAbs_EDGE)}
    out = []
    for e in explore(fb, TopAbs_EDGE):
        e = TopoDS.Edge_s(e)
        for e2 in ea.values():
            if e.IsSame(e2) and not BRep_Tool.Degenerated_s(e):
                out.append(e)
    return out


def join_check(fa, fb, n=600, end_skip=0.0, label=""):
    """Separation, tangent-plane angle and curvature in a common frame along the shared edge(s)."""
    sa, sb = srf_of(fa), srf_of(fb)
    rows = []
    for e in shared_edges(fa, fb):
        c = BRepAdaptor_Curve(e)
        t0, t1 = c.FirstParameter(), c.LastParameter()
        ts = np.linspace(t0, t1, n)
        for t in ts:
            p = c.Value(float(t))
            p = np.array([p.X(), p.Y(), p.Z()])
            v = c.DN(float(t), 1)
            tan = np.array([v.X(), v.Y(), v.Z()])
            if np.linalg.norm(tan) == 0:
                continue
            tan /= np.linalg.norm(tan)
            ua = project(sa, p); ub = project(sb, p)
            if ua is None or ub is None:
                continue
            A = props_at(sa, ua[0], ua[1]); Bp = props_at(sb, ub[0], ub[1])
            if "n" not in A or "n" not in Bp:
                continue
            nA = A["n"]; nB = Bp["n"] * np.sign(np.dot(Bp["n"], nA))
            ang = np.degrees(np.arccos(np.clip(np.dot(nA, nB), -1, 1)))
            fa_ = second_form(A, nA, tan); fb_ = second_form(Bp, nA, tan)
            rows.append(dict(p=p, sep=float(np.linalg.norm(A["p"] - Bp["p"])), ang=float(ang),
                             ka=fa_, kb=fb_))
    return rows


def summarize_join(rows, exclude_near=(), excl_r=0.0):
    keep = [r for r in rows if all(np.linalg.norm(r["p"][:2] - np.array(q)) > excl_r for q in exclude_near)]
    sep = np.array([r["sep"] for r in keep]); ang = np.array([r["ang"] for r in keep])
    out = dict(n=len(keep), max_separation_mm=float(sep.max()), max_angle_deg=float(ang.max()),
               at_max_angle=keep[int(np.argmax(ang))]["p"].round(4).tolist())
    dk = []
    for r in keep:
        if r["ka"] is None or r["kb"] is None:
            continue
        a, b = np.array(r["ka"]), np.array(r["kb"])
        dk.append(np.abs(a - b))
    if dk:
        dk = np.array(dk)
        out["max_curvature_difference_1_per_mm"] = dict(along_edge=float(dk[:, 0].max()),
                                                       across_edge=float(dk[:, 1].max()),
                                                       twist=float(dk[:, 2].max()))
        j = int(np.argmax(dk[:, 1]))
        out["at_max_across_difference"] = keep[j]["p"].round(4).tolist()
    return out


# ----------------------------------------------------------------------------- radius map
def radius_map_face(face, nu=900, nv=260):
    srf = srf_of(face)
    m = rm.map_face(srf, nu=nu, nv=nv)
    return srf, m


def refine_min(srf, m, mask, n_seeds=6, zone=None):
    """Refine local minima of the principal radius inside mask (grid seeds -> Nelder-Mead)."""
    R = np.where(mask, m["R"], np.nan)
    flat = np.argsort(np.nan_to_num(R, nan=np.inf).ravel())[: n_seeds * 20]
    seeds = []
    for f in flat:
        i = np.unravel_index(f, R.shape)
        if all(abs(i[0] - s[0]) + abs(i[1] - s[1]) > 6 for s in seeds):
            seeds.append(i)
        if len(seeds) >= n_seeds:
            break
    u0, u1, v0, v1 = srf.Bounds()
    best = None
    for i in seeds:
        x0 = np.array([m["U"][i], m["V"][i]])

        def f(x):
            u = np.clip(x[0], u0, u1); v = np.clip(x[1], v0, v1)
            pr = props_at(srf, u, v)
            if "k" not in pr:
                return 1e9
            if zone is not None and not zone(pr["p"]):
                return 1e9
            return 1 / max(abs(pr["k"][0]), abs(pr["k"][1]))
        r = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-7, maxiter=400))
        u = np.clip(r.x[0], u0, u1); v = np.clip(r.x[1], v0, v1)
        p = props_at(srf, u, v)["p"]
        if best is None or r.fun < best[0]:
            best = (float(r.fun), u, v, p)
    return best


def judge_face(face, label):
    srf, m = radius_map_face(face)
    R, xi, t, thick, P = m["R"], m["xi"], m["t"], m["thick"], m["P"]
    LE = np.array([cm.SWEEP * cm.ATT, cm.ATT])
    dLE = np.hypot(P[..., 0] - LE[0], np.abs(P[..., 1]) - LE[1])
    dTE = np.hypot(P[..., 0] - cm.P0_XY[0], np.abs(P[..., 1]) - cm.P0_XY[1])
    u0 = srf.Bounds()[0]
    collapsed = m["U"] <= u0 + 1e-12          # the collapsed LE edge itself
    res = {"grid": [int(R.shape[0]), int(R.shape[1])]}
    ok = ~collapsed & np.isfinite(R)
    res["closest_grid_point_to_collapsed_LE_corner_mm"] = float(dLE[ok].min())
    # thick zone
    mask = ok & thick
    in_thick = lambda p: rm.t_att(np.clip(cm.chord_fraction_of(p[0], p[1]), 0, 1)) >= 0.8
    b = refine_min(srf, m, mask, zone=in_thick)
    xb = cm.chord_fraction_of(b[3][0], b[3][1])
    res["thick"] = dict(rule="t_att >= 0.8 mm: every radius >= 0.4 mm", R_min=b[0], at=b[3].round(4).tolist(),
                        chord_fraction=float(xb), t_att=float(rm.t_att(np.clip(xb, 0, 1))),
                        grid_R_min=float(np.nanmin(np.where(mask, R, np.nan))))
    # thin zones
    for name, sel in (("thin_LE", ok & ~thick & (xi < 0.5)), ("thin_TE", ok & ~thick & (xi >= 0.5))):
        if not sel.any():
            continue
        lip = sel & (np.abs(P[..., 2]) < 1e-9)
        entry = {}
        if lip.any():
            ratio = R / (t / 2)
            j = np.unravel_index(np.nanargmin(np.where(lip, ratio, np.nan)), R.shape)
            entry["lip"] = dict(rule="t_att < 0.8 mm: lip half-round of radius t_att/2",
                                min_ratio_R_over_half_t=float(ratio[j]), R=float(R[j]), t_att=float(t[j]),
                                at=P[j].round(4).tolist())
            k = np.unravel_index(np.nanargmax(np.where(lip, ratio, np.nan)), R.shape)
            entry["lip"]["max_ratio_R_over_half_t"] = float(ratio[k])
        bands = {}
        for r_ex in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0):
            dd = dLE if name == "thin_LE" else dTE
            s2 = sel & (dd > r_ex)
            if s2.any():
                j = np.unravel_index(np.nanargmin(np.where(s2, R, np.nan)), R.shape)
                bands[f"beyond_{r_ex}mm_of_corner"] = dict(R_min=float(R[j]), t_att=float(t[j]),
                                                           ratio_R_over_half_t=float(R[j] / (t[j] / 2)),
                                                           at=P[j].round(4).tolist(),
                                                           distance_to_corner_mm=float(dd[j]))
        entry["interior_min_by_distance_from_corner"] = bands
        if name == "thin_TE" and lip.any():
            apex = np.array([cm.P0_XY[0] + cm.TE_HALF * cm.E_A[0], cm.P0_XY[1] + cm.TE_HALF * cm.E_A[1]])
            da = np.hypot(P[..., 0] - apex[0], np.abs(P[..., 1]) - apex[1])
            ratio = R / (t / 2)
            prof = {}
            for lo_, hi_ in ((0, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.5), (0.5, 1), (1, 3), (3, 30)):
                s3 = lip & (da >= lo_) & (da < hi_)
                if s3.any():
                    prof[f"{lo_}-{hi_}mm"] = [float(np.nanmin(ratio[s3])), float(np.nanmax(ratio[s3]))]
            entry["lip_ratio_R_over_half_t_by_distance_from_apex"] = prof
        res[name] = entry
    return res, m


# ----------------------------------------------------------------------------- sections
def section_points(faces, plane_pt, plane_n, step=0.02):
    """(face, u, v, xyz) samples of the plane section, (u, v) from the section's pcurves."""
    pln = gp_Pln(gp_Pnt(*map(float, plane_pt)), gp_Dir(*map(float, plane_n)))
    pts = []
    for f in faces:
        pf = BRepBuilderAPI_MakeFace(pln).Face()
        sec = BRepAlgoAPI_Section(f, pf, False)
        sec.ComputePCurveOn1(True)
        sec.Approximation(True)
        sec.Build()
        srf = srf_of(f)
        for e in explore(sec.Shape(), TopAbs_EDGE):
            e = TopoDS.Edge_s(e)
            try:
                c2 = BRepAdaptor_Curve2d(e, f)
            except Exception:
                continue
            a, b = c2.FirstParameter(), c2.LastParameter()
            p0 = srf.Value(*c2.Value(a).Coord()); p1 = srf.Value(*c2.Value(b).Coord())
            L = p0.Distance(p1)
            n = int(np.clip(300 + L / step, 300, 4000))
            for t in np.linspace(a, b, n):
                u, v = c2.Value(float(t)).Coord()
                pts.append((f, u, v))
    return pts

def section_curvature(pts, plane_n):
    """Signed curvature of the planar section, from the surface's second form (Meusnier)."""
    m = np.asarray(plane_n, float); m /= np.linalg.norm(m)
    rows = []
    for f, u, v in pts:
        s = srf_of(f)
        pr = props_at(s, u, v)
        p = pr["p"]
        if "k" not in pr or "n" not in pr:
            continue
        n = pr["n"] * (1 if f.Orientation().name == "TopAbs_FORWARD" else -1)   # outward
        t = np.cross(n, m); t /= np.linalg.norm(t)
        N = np.cross(m, t)
        k1, k2 = pr["k"]; d1, d2 = pr["d"]
        sg = 1 if f.Orientation().name == "TopAbs_FORWARD" else -1
        kn = sg * (k1 * np.dot(t, d1) ** 2 + k2 * np.dot(t, d2) ** 2)    # w.r.t. outward normal
        denom = np.dot(n, N)
        if abs(denom) < 1e-6:
            continue
        rows.append((p, -kn / denom))     # convex (bending away from outward normal) > 0
    return rows


def analyse_section(rows, label, xi_ref=None):
    if not rows:
        return dict(label=label, n=0)
    P = np.array([r[0] for r in rows]); k = np.array([r[1] for r in rows])
    R = 1 / np.maximum(np.abs(k), 1e-12)
    i = int(np.argmin(R))
    lip = np.abs(P[:, 2]) < 2e-3
    out = dict(label=label, n=int(len(P)), R_min=float(R[i]), at=P[i].round(4).tolist(),
               concave_points=int((k < -1e-3).sum()), min_curvature=float(k.min()))
    if (k < -1e-3).any():
        j = int(np.argmin(k))
        out["most_concave"] = dict(R=float(1 / abs(k[j])), at=P[j].round(4).tolist())
    if lip.any():
        out["lip_R"] = [float(x) for x in np.round(R[lip][np.argsort(P[lip, 0])][[0, -1]], 4)]
    xi_min = cm.chord_fraction_of(P[i, 0], P[i, 1])
    out["chord_fraction_at_R_min"] = float(xi_min)
    out["t_att_at_R_min"] = float(rm.t_att(np.clip(xi_min, 0, 1)))
    if xi_ref is not None:
        out["t_att_section"] = float(rm.t_att(np.clip(xi_ref, 0, 1)))
    return out


def sections(capR):
    faces = [capR["R_upper"], capR["R_lower"]]
    out = {"station_cuts": [], "sweep_cuts": [], "lip_normal_cuts": []}
    for Y in (390.5, 392, 395, 400, 410, 420, 430, 440, 445, 448, 449.5):
        rows = section_curvature(section_points(faces, (0, Y, 0), (0, 1, 0)), (0, 1, 0))
        out["station_cuts"].append(analyse_section(rows, f"Y = {Y}"))
    for xi in (0.0005, 0.002, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.97, 0.98, 0.99):
        d = np.array([cm.SWEEP - cm.TAPER * xi, 1.0, 0.0])
        nrm = np.array([1.0, -(cm.SWEEP - cm.TAPER * xi), 0.0])
        p0 = np.array([cm.SWEEP * cm.ATT + cm.C_ATT * xi, cm.ATT, 0.0])
        rows = section_curvature(section_points(faces, p0, nrm), nrm)
        out["sweep_cuts"].append(analyse_section(rows, f"constant chord fraction {xi}", xi))
    for p in (0.01, 0.03, 0.06, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97, 0.98, 0.99, 0.995, 0.998):
        xl, rl = cm.lip(p)
        X, Y = cm.plan_xy(xl, rl)
        nin = cm.lip_inward_normal(p)
        tan = np.array([-nin[1], nin[0], 0.0])        # plane normal = lip tangent
        rows = section_curvature(section_points(faces, (X, Y, 0), tan), tan)
        # keep the branch near the lip point (the plane may cut the cap elsewhere too)
        rows = [r for r in rows if np.hypot(r[0][0] - X, r[0][1] - Y) < 25]
        out["lip_normal_cuts"].append(analyse_section(rows, f"normal to the lip at p = {p}", float(xl)))
    return out


# ----------------------------------------------------------------------------- zebra
def zebra_join(fa, fb, views, stripe_axes, n=300, ds=0.01):
    """Numerical zebra across a join: stripe phase theta = asin(r . e) of the reflected view ray
    r = d - 2 (d.n) n (CAD normals of the reimported faces).  At each join point theta and its
    gradient (along the join, and across it into each face) are compared on both faces; at the
    actual stripe boundaries (theta = k * w) on the join, the stripe direction change is measured."""
    sa, sb = srf_of(fa), srf_of(fb)
    width = np.radians(5.0)
    res = []
    for e in shared_edges(fa, fb):
        c = BRepAdaptor_Curve(e)
        ts = np.linspace(c.FirstParameter(), c.LastParameter(), n)[1:-1]
        for t in ts:
            p = np.array(c.Value(float(t)).Coord())
            v = c.DN(float(t), 1); tan = np.array([v.X(), v.Y(), v.Z()]); tan /= np.linalg.norm(tan)
            ua, ub = project(sa, p), project(sb, p)
            if ua is None or ub is None:
                continue
            nA = props_at(sa, ua[0], ua[1]).get("n")
            if nA is None:
                continue
            stencils = []
            for srf, uv in ((sa, ua), (sb, ub)):
                pr = props_at(srf, uv[0], uv[1])
                n0 = pr["n"] * np.sign(np.dot(pr["n"], nA))
                b = np.cross(n0, tan); b /= np.linalg.norm(b)
                # which way is "into this face"?  (the side whose offset point lies closer to it)
                pp_, pm_ = project(srf, p + ds * b), project(srf, p - ds * b)
                sign = 1.0 if (pp_ is not None and (pm_ is None or pp_[2] < pm_[2])) else -1.0
                ns = {}
                for key, off in (("0", 0 * tan), ("t", ds * tan), ("b", ds * b * sign), ("bb", 2 * ds * b * sign)):
                    q = project(srf, p + off)
                    nn = props_at(srf, q[0], q[1])["n"]
                    ns[key] = nn * np.sign(np.dot(nn, nA))
                stencils.append((ns, sign))
            for d in views:
                for e_ax in stripe_axes:
                    th = []
                    for ns, sign in stencils:
                        val = {}
                        for key, nn in ns.items():
                            r = d - 2 * np.dot(d, nn) * nn
                            val[key] = np.arcsin(np.clip(np.dot(r, e_ax), -1, 1))
                        g_t = (val["t"] - val["0"]) / ds
                        g_b = (-3 * val["0"] + 4 * val["b"] - val["bb"]) / (2 * ds) * sign
                        th.append((val["0"], g_t, g_b))
                    res.append(dict(p=p, d=d, e=e_ax, thA=th[0][0], thB=th[1][0], gA=th[0][1:], gB=th[1][1:]))
    if not res:
        return {}
    zebra_join.last = (res, width)
    return summarize_zebra(res, width)


def summarize_zebra(res, width, exclude=(), excl_r=0.0):
    res = [r for r in res if all(np.hypot(r["p"][0] - q[0], abs(r["p"][1]) - q[1]) > excl_r for q in exclude)]
    if not res:
        return {}
    dth = np.array([abs(r["thA"] - r["thB"]) for r in res])
    dgb = np.array([abs(r["gA"][1] - r["gB"][1]) for r in res])
    gmag = np.array([np.hypot(*r["gA"]) for r in res])
    dir_change = []
    for r in res:
        # stripe boundary crossing near this join point?
        k = np.round(r["thA"] / width)
        if abs(r["thA"] - k * width) < 0.15 * width:
            a = np.degrees(np.arctan2(r["gA"][1], r["gA"][0])); b = np.degrees(np.arctan2(r["gB"][1], r["gB"][0]))
            dd = (a - b + 180) % 360 - 180
            dir_change.append((abs(dd), r))
    out = dict(samples=len(res), max_phase_jump_deg=float(np.degrees(dth.max())),
               at_max_phase_jump=res[int(np.argmax(dth))]["p"].round(4).tolist(),
               max_cross_gradient_jump_deg_per_mm=float(np.degrees(dgb.max())),
               median_gradient_deg_per_mm=float(np.degrees(np.median(gmag))),
               at_max_gradient_jump=res[int(np.argmax(dgb))]["p"].round(4).tolist())
    if dir_change:
        dmax = max(dir_change, key=lambda x: x[0])
        out["stripe_boundary_crossings"] = len(dir_change)
        out["max_stripe_direction_change_deg"] = float(dmax[0])
        out["at_max_direction_change"] = dmax[1]["p"].round(4).tolist()
    out["p95_cross_gradient_jump_deg_per_mm"] = float(np.degrees(np.percentile(dgb, 95)))
    return out


def zebra_join_all(fa, fb, views, axes, n):
    """Zebra summary over the whole join and away from the LE point / TE corner."""
    global _ZEB_RES
    full = zebra_join(fa, fb, views, axes, n)
    return full


# ----------------------------------------------------------------------------- main
def log(*a):
    print(*a, flush=True)


def main():
    t0 = time.time()
    data = open(STEP, "rb").read()
    rb = {"file": "wing.step", "step_sha256": hashlib.sha256(data).hexdigest(), "step_bytes": len(data)}
    shape = read(STEP)
    topo, faces = topology(shape)
    rb["reimport"] = topo
    cls = classify(faces)
    capR, capL = ({k: v for k, v in cls["cap"].items() if k.startswith(s)} for s in "RL")
    rb["faces"] = dict(parent=len(cls["parent"]), cap=sorted(cls["cap"]))
    log("topology", round(time.time() - t0, 1))
    rb["parent_preserved"] = parent_preserved(cls["parent"])
    log("parent", round(time.time() - t0, 1))
    # locked dimensions
    bb = topo["bounding_box_mm"]
    lipR = lip_points(capR["R_upper"]); lipL = lip_points(capL["L_upper"])
    rb["locked"] = dict(projected_span_mm=bb[4] - bb[1], y_extent_mm=[bb[1], bb[4]],
                        attachment_station_mm=[float(min(bbox(f)[1] for f in capR.values())),
                                               float(max(bbox(f)[4] for f in capL.values()))],
                        outline_right=outline_deviation(lipR), outline_left=outline_deviation(lipL),
                        tip_point=lipR[np.argmax(lipR[:, 1])].round(5).tolist())
    log("locked", round(time.time() - t0, 1))
    # joins
    par = cls["parent"]

    def find_parent_for(f):
        return [g for g in par if shared_edges(f, g)]
    joins = {}
    LE_pts = [(cm.SWEEP * cm.ATT, cm.ATT), (cm.SWEEP * cm.ATT, -cm.ATT)]
    TE_pts = [tuple(cm.P0_XY), (cm.P0_XY[0], -cm.P0_XY[1])]
    for key in ("R_upper", "R_lower", "L_upper", "L_lower"):
        f = cls["cap"][key]
        for g in find_parent_for(f):
            rows = join_check(f, g, n=800)
            joins[f"{key} / parent skin"] = dict(
                all_points=summarize_join(rows),
                excluding_0p05mm_of_LE_point_and_TE_corner=summarize_join(rows, LE_pts + TE_pts, 0.05))
    for s in "RL":
        rows = join_check(cls["cap"][s + "_upper"], cls["cap"][s + "_lower"], n=1500)
        joins[f"{s}_upper / {s}_lower (lip)"] = summarize_join(rows)
        for g in find_parent_for(cls["cap"][s + "_disc"]):
            rows = join_check(cls["cap"][s + "_disc"], g, n=100)
            joins[f"{s}_disc / parent TE face"] = summarize_join(rows)
        rows = join_check(cls["cap"][s + "_upper"], cls["cap"][s + "_disc"], n=200)
        angs = [(round(float(np.degrees(np.arctan2(r["p"][2], np.dot(r["p"][:2] - np.array(TE_pts[0 if s == 'R' else 1]),
                                                                       np.array([cm.E_A[0], cm.E_A[1] * (1 if s == 'R' else -1)]))))), 1),
                 round(r["ang"], 2)) for r in rows[::20]]
        joins[f"{s}_upper / {s}_disc (designed crease along the half-round)"] = dict(
            summarize_join(rows), crease_deg_vs_angle_around_half_round=angs)
    rb["joins"] = joins
    log("joins", round(time.time() - t0, 1))
    # radius maps
    radii = {}
    mR = None
    for key in ("R_upper", "R_lower", "L_upper", "L_lower"):
        nu, nv = (900, 260) if key == "R_upper" else (300, 120)
        srf_ = srf_of(cls["cap"][key])
        res, m = judge_face(cls["cap"][key], key) if key == "R_upper" else (None, None)
        if key != "R_upper":
            # the other three patches are mirror images: compare their poles with the mirrored R_upper
            ref = BRepAdaptor_Surface(cls["cap"]["R_upper"]).BSpline()
            s2 = BRepAdaptor_Surface(cls["cap"][key]).BSpline()
            sy = -1 if key.startswith("L") else 1
            sz = -1 if key.endswith("lower") else 1
            dmax = 0.0
            for i in range(1, ref.NbUPoles() + 1):
                for j in range(1, ref.NbVPoles() + 1):
                    a = ref.Pole(i, j); b = s2.Pole(i, j)
                    dmax = max(dmax, abs(a.X() - b.X()), abs(sy * a.Y() - b.Y()), abs(sz * a.Z() - b.Z()))
            res = dict(mirror_of="R_upper", max_pole_difference_mm=dmax,
                       note="radii identical to R_upper by symmetry")
        else:
            mR = m
        radii[key] = res
    rb["radius_map"] = radii
    log("radii", round(time.time() - t0, 1))
    rb["sections"] = sections(capR)
    log("sections", round(time.time() - t0, 1))
    # zebra
    views = [np.array(v, float) / np.linalg.norm(v) for v in
             ((0, 0, -1), (0.5, 0.3, -0.8), (-0.6, 0.2, -0.75), (0.2, -0.6, -0.75), (0.7, 0.0, -0.3))]
    axes = [np.array(a, float) / np.linalg.norm(a) for a in ((1, 0, 0), (0, 1, 0), (1, 1, 0), (1, -1, 0.3))]
    zeb = {}
    corners = [(cm.SWEEP * cm.ATT, cm.ATT), tuple(cm.P0_XY)]
    for key in ("R_upper", "R_lower"):
        f = cls["cap"][key]
        for g in find_parent_for(f):
            zeb[f"{key} / parent skin"] = dict(all=zebra_join(f, g, views, axes, n=160))
            res_, w_ = zebra_join.last
            for r_ex in (0.5, 1.0, 3.0):
                zeb[f"{key} / parent skin"][f"beyond_{r_ex}mm_of_LE_point_and_TE_corner"] = summarize_zebra(res_, w_, corners, r_ex)
    zeb["R_upper / R_lower (lip)"] = zebra_join(cls["cap"]["R_upper"], cls["cap"]["R_lower"],
                                                [-v for v in views[:3]] + views[:2], axes, n=160)
    rb["zebra"] = dict(stripe_width_deg=5.0, views=[v.round(3).tolist() for v in views],
                       stripe_axes=[a.round(3).tolist() for a in axes], joins=zeb)
    log("zebra", round(time.time() - t0, 1))
    np.savez_compressed(os.path.join(HERE, "radius_map_R_upper.npz"), U=mR["U"], V=mR["V"], P=mR["P"].astype(np.float32),
                        R=mR["R"].astype(np.float32), xi=mR["xi"].astype(np.float32))
    rb["runtime_s"] = round(time.time() - t0, 1)
    def clean(o):
        if isinstance(o, dict):
            return {k: clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if hasattr(o, "tolist"):
            return clean(o.tolist())
        if isinstance(o, float) and not np.isfinite(o):
            return None if np.isnan(o) else ("inf" if o > 0 else "-inf")
        return o
    json.dump(clean(rb), open(OUT, "w"), indent=1, allow_nan=False)
    print(json.dumps({k: rb[k] for k in ("step_sha256", "reimport", "parent_preserved", "locked")}, indent=1))


if __name__ == "__main__":
    main()
