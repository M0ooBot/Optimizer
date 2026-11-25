"""Verify wing.step for fin40/blunt on the reimported STEP and write readback.json.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/fin40/blunt/verify_cap.py

Follows the verification part of .agents/skills/close-rounded-wingtip/SKILL.md:
1. topology, validity, preserved parent skin, locked dimensions, joins
   (separation, tangent-plane angle, curvature in a common frame) and the
   shared perimeter between the two cap patches;
2. plane sections through the cap (constant chord position and cuts normal to
   the outline) with signed curvature along each trace;
3. a dense principal-curvature map of every cap face with local refinement;
4. a numerical zebra check across both wing-to-cap joins (and the lip).
All measurements use CAD surfaces and derivatives of the reimported file.
"""
import json
import pathlib
import sys
import time

import numpy as np
from scipy.optimize import minimize

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cap_def as cd  # noqa: E402
import occ_util as ou  # noqa: E402

from OCP.BRepCheck import BRepCheck_Analyzer  # noqa: E402
from OCP.TopExp import TopExp_Explorer  # noqa: E402
from OCP.TopAbs import TopAbs_SOLID, TopAbs_SHELL, TopAbs_FACE, TopAbs_EDGE  # noqa: E402
from OCP.TopoDS import TopoDS  # noqa: E402
from OCP.BRep import BRep_Tool  # noqa: E402
from OCP.ShapeAnalysis import ShapeAnalysis_FreeBounds  # noqa: E402
from OCP.GProp import GProp_GProps  # noqa: E402
from OCP.BRepGProp import BRepGProp  # noqa: E402
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section  # noqa: E402
from OCP.BRepAdaptor import BRepAdaptor_Curve  # noqa: E402
from OCP.gp import gp_Pln, gp_Pnt, gp_Dir  # noqa: E402
from OCP.GeomAbs import GeomAbs_SurfaceType  # noqa: E402

P = cd.parent
STEP = HERE / "wing.step"
OUT = HERE / "readback.json"
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:6.1f}s]", *a, flush=True)


def count(shape, kind):
    n, e = 0, TopExp_Explorer(shape, kind)
    while e.More():
        n += 1
        e.Next()
    return n


def thickness_at(x):
    """Parent attachment thickness at chord position x (attachment frame: LE 65, chord 85)."""
    xi = np.clip((np.asarray(x, float) - 65.0) / 85.0, 0.0, 1.0)
    t = P.section_thickness(xi, cd.Z0)
    return np.nan_to_num(t, nan=P.TE_THICKNESS)


def required_radius(x):
    t = thickness_at(x)
    return np.where(t >= 0.8, 0.4, t / 2.0), t


# ---------------------------------------------------------------- 1. topology
def topology(shape, F):
    ana = BRepCheck_Analyzer(shape, True)
    ana_exact = BRepCheck_Analyzer(shape, True, False, True)
    fb = ShapeAnalysis_FreeBounds(shape)
    closed = count(fb.GetClosedWires(), TopAbs_EDGE) == 0 and count(fb.GetOpenWires(), TopAbs_EDGE) == 0
    shells_closed = []
    e = TopExp_Explorer(shape, TopAbs_SHELL)
    while e.More():
        shells_closed.append(bool(BRep_Tool.IsClosed_s(e.Current())))
        e.Next()
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props, 1e-9)
    tol_max = 0.0
    e = TopExp_Explorer(shape, TopAbs_EDGE)
    while e.More():
        tol_max = max(tol_max, BRep_Tool.Tolerance_s(TopoDS.Edge_s(e.Current())))
        e.Next()
    return {
        "solids": count(shape, TopAbs_SOLID), "shells": count(shape, TopAbs_SHELL),
        "faces": count(shape, TopAbs_FACE),
        "face_names_and_types": {k: str(v.ad.GetType()).split(".")[-1] for k, v in F.items()},
        "brepcheck_valid": bool(ana.IsValid()), "brepcheck_exact_valid": bool(ana_exact.IsValid()),
        "free_edges_none": bool(closed), "shells_closed": shells_closed,
        "volume_mm3": props.Mass(), "max_edge_tolerance_mm": tol_max,
        "bounding_box_mm": ou.bbox(shape),
    }


# ------------------------------------------------------- parent preservation
def parent_preserved(F):
    out = {}
    s = np.linspace(0, 1, 81)
    z = np.linspace(0, 240, 49)
    S, Z = np.meshgrid(s, z)
    for side, name in ((+1, "parent_up"), (-1, "parent_lo")):
        pts = P.skin(S.ravel(), Z.ravel(), side)
        gaps = [F[name].uv(p)[2] for p in pts]
        out[name + "_max_distance_to_parent_py_mm"] = float(np.max(gaps))
    tt = np.linspace(0, 1, 9)
    T, Z2 = np.meshgrid(tt, z)
    pts = P.te_face(T.ravel(), Z2.ravel())
    out["te_face_max_distance_to_parent_py_mm"] = float(max(F["te"].uv(p)[2] for p in pts))
    return out


def parent_step_identity(F):
    """Distance from parent.step's skin and trailing-edge surfaces to the wing's corresponding faces."""
    shp = ou.read_step(cd.PARENT_DIR / "parent.step")
    res = {}
    for f in ou.faces(shp):
        G = ou.Face(f)
        x0, y0, z0, x1, y1, z1 = G.bbox
        if z1 < 1e-3 or z0 > 240 - 1e-3:
            continue                          # root face and parent.step's flat closing face
        if y1 - y0 < 1.0:
            name = "te"                       # the 0.4 mm trailing-edge face (a B-spline in parent.step)
        else:
            name = "parent_up" if (y0 + y1) > 0 else "parent_lo"
        d = 0.0
        for u in np.linspace(G.u0, G.u1, 61):
            for v in np.linspace(G.v0, G.v1, 61):
                p = G.d2(u, v)[0]
                if name == "te" and not G.inside(u, v):
                    continue
                d = max(d, F[name].uv(p)[2])
        res[name] = d
    return res


# ---------------------------------------------------------- locked dimensions
def section_points(shape, pln, n_per_edge=400):
    sec = BRepAlgoAPI_Section(shape, pln, False)
    sec.Approximation(True)
    sec.Build()
    pts = []
    e = TopExp_Explorer(sec.Shape(), TopAbs_EDGE)
    while e.More():
        c = BRepAdaptor_Curve(TopoDS.Edge_s(e.Current()))
        for t in np.linspace(c.FirstParameter(), c.LastParameter(), n_per_edge):
            q = c.Value(t)
            pts.append((q.X(), q.Y(), q.Z()))
        e.Next()
    return np.array(pts)


def polyline_distance(q, poly):
    """Distance from points q (n,2) to a dense polyline poly (m,2)."""
    from scipy.spatial import cKDTree
    tree = cKDTree(poly)
    _, i = tree.query(q)
    best = np.full(len(q), np.inf)
    for off in (-1, 0):
        a = np.clip(i + off, 0, len(poly) - 2)
        p0, p1 = poly[a], poly[a + 1]
        d = p1 - p0
        t = np.clip(np.einsum("ij,ij->i", q - p0, d) / np.einsum("ij,ij->i", d, d), 0, 1)
        best = np.minimum(best, np.linalg.norm(q - (p0 + t[:, None] * d), axis=1))
    return best


def locked_dimensions(shape, F):
    out = {}
    bb = ou.bbox(shape)
    crown = F["cap_up"].d2(1.0, 1.0)[0]
    out["crown_point_mm"] = crown.tolist()
    out["z_min_mm"], out["z_max_mm"] = bb[2], bb[5]
    out["finished_height_from_crown_point_mm"] = float(crown[2])
    # outline: section of the solid with Y = 0, compared with the locked outline (both directions)
    sec = section_points(shape, gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0)), 3000)
    cap_sec = sec[sec[:, 2] > cd.Z0 + 1e-9][:, [0, 2]]
    h = np.unique(np.concatenate([np.linspace(0, 1, 40001), 1 - np.geomspace(1e-12, 0.05, 20000),
                                  np.geomspace(1e-9, 0.05, 20000)]))
    le = P.tip_outline(h, -1)[:, [0, 2]]
    te = P.tip_outline(h, +1)[:, [0, 2]]
    poly = np.concatenate([le, te[::-1]])
    keep = np.concatenate([[True], np.linalg.norm(np.diff(poly, axis=0), axis=1) > 1e-12])
    poly = poly[keep]
    d_sec = polyline_distance(cap_sec, poly)
    i = int(np.argmax(d_sec))
    out["outline_max_deviation_section_to_outline_mm"] = float(d_sec[i])
    out["outline_max_deviation_at_xz"] = cap_sec[i].tolist()
    wc = cd.te_closure_height()
    lipreg = cap_sec[:, 1] > cd.Z0 + wc + 1e-3
    j = int(np.argmax(np.where(lipreg, d_sec, -1)))
    out["outline_max_deviation_where_lip_forms_silhouette_mm (z > closure)"] = float(d_sec[j])
    out["outline_max_deviation_where_lip_forms_silhouette_at_xz"] = cap_sec[j].tolist()
    out["outline_deviation_note"] = ("below the trailing-edge face closure the silhouette is the flat trailing-edge "
                                     "face, which lies on the parent's trailing-edge line; the locked outline departs "
                                     "from that line by a h^3/2, at most 7.1e-4 mm at the closure height")
    # outline -> model: exact distance from outline points to the Y = 0 section edges
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    sec_shape = BRepAlgoAPI_Section(shape, gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0)), False)
    sec_shape.Approximation(True)
    sec_shape.Build()
    hq = np.unique(np.concatenate([np.linspace(0, 1, 801), 1 - np.geomspace(1e-8, 0.05, 100),
                                   np.geomspace(1e-6, 0.06, 150)]))
    worst = (0.0, None)
    for branch in (-1, 1):
        for p in P.tip_outline(hq, branch):
            ds = BRepExtrema_DistShapeShape(BRepBuilderAPI_MakeVertex(gp_Pnt(*map(float, p))).Vertex(), sec_shape.Shape())
            if ds.Value() > worst[0]:
                worst = (ds.Value(), [float(p[0]), float(p[2])])
    out["outline_max_deviation_outline_to_section_mm"] = worst[0]
    out["outline_max_deviation_outline_to_section_at_xz"] = worst[1]
    # where the trailing-edge face is part of the silhouette the deviation is the outline's own
    # departure from the TE line, sigma_o(w) = TE(z) - x_T(z)
    wc = cd.te_closure_height()
    xL, XT, TE = cd.ends_values(cd.Z0 + wc)
    out["te_face_closure_height_above_attachment_mm"] = wc
    out["outline_departure_from_te_line_at_closure_mm"] = float(TE - P.tip_outline(np.array([wc / 40.0]), +1)[0, 0])
    w = np.linspace(wc, cd.WE, 400)
    xL, XT, TE = cd.ends_values(cd.Z0 + w)
    xt = P.tip_outline(w / 40.0, +1)[:, 0]
    out["design_lip_offset_from_outline_above_closure_max_mm"] = float(np.max(XT - xt))
    # chords and trailing-edge thickness
    out["attachment_section"] = {"LE_x": float(P.skin(0.0, 240.0, 1)[0]), "TE_x": float(P.skin(1.0, 240.0, 1)[0]),
                                 "chord_mm": float(P.skin(1.0, 240.0, 1)[0] - P.skin(0.0, 240.0, 1)[0])}
    te = F["te"]
    widths = {}
    for zz in (0.0, 60.0, 120.0, 180.0, 239.0, 240.0):
        x = 160.0 - 10.0 * zz / 240.0
        u, v, g = te.uv((x, 0.0, zz))
        # width of the TE face at this height: march along the face's V (=Y) direction
        ys = section_points(te.face, gp_Pln(gp_Pnt(0, 0, zz), gp_Dir(0, 0, 1)), 50)
        widths[str(zz)] = float(ys[:, 1].max() - ys[:, 1].min()) if len(ys) else None
    out["te_face_width_mm_at_z"] = widths
    rs = section_points(shape, gp_Pln(gp_Pnt(0, 0, 0.0), gp_Dir(0, 0, 1)), 400)
    out["root_chord_mm"] = float(rs[:, 0].max() - rs[:, 0].min()) if len(rs) else None
    return out


# -------------------------------------------------------------------- joins
def frame_compare(A, la, Bf, lb, e1):
    """Separation, normal angle and II components in a common frame (e1 along the join)."""
    n = la["n"]
    e1 = e1 - (e1 @ n) * n
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    comp = lambda W: np.array([e1 @ W @ e1, e1 @ W @ e2, e2 @ W @ e2])
    ca, cb = comp(la["W"]), comp(lb["W"])
    ang = np.degrees(np.arccos(np.clip(la["n"] @ lb["n"], -1, 1)))
    return float(np.linalg.norm(la["p"] - lb["p"])), float(ang), ca, cb


def join_report(F, face_a, face_b, pts, tangents, labels):
    A, Bf = F[face_a], F[face_b]
    rows = []
    for p, t, lab in zip(pts, tangents, labels):
        ua, va, ga = A.uv(p)
        ub, vb, gb = Bf.uv(p)
        la, lb = A.local(ua, va), Bf.local(ub, vb)
        if la is None or lb is None:
            continue
        sep, ang, ca, cb = frame_compare(A, la, Bf, lb, t)
        rows.append((lab, sep, ang, *ca, *cb, ga, gb))
    R = np.array(rows, float)
    dk = np.abs(R[:, 3:6] - R[:, 6:9])
    imax = lambda col: int(np.argmax(col))
    rep = {
        "samples": len(R),
        "max_separation_mm": float(R[:, 1].max()), "at": float(R[imax(R[:, 1]), 0]),
        "max_normal_angle_deg": float(R[:, 2].max()), "angle_at": float(R[imax(R[:, 2]), 0]),
        "max_abs_curvature_difference_per_mm": {
            "along_join_II(e1,e1)": float(dk[:, 0].max()),
            "twist_II(e1,e2)": float(dk[:, 1].max()),
            "across_join_II(e2,e2)": float(dk[:, 2].max())},
        "across_join_curvature_at_worst": {"label": float(R[imax(dk[:, 2]), 0]),
                                           face_a: float(R[imax(dk[:, 2]), 5]),
                                           face_b: float(R[imax(dk[:, 2]), 8])},
        "max_projection_gap_mm": float(max(R[:, 9].max(), R[:, 10].max())),
    }
    return rep, R


def joins(F):
    out = {}
    s = np.unique(np.concatenate([np.linspace(0, 1, 2001), np.geomspace(1e-6, 0.01, 200),
                                  1 - np.geomspace(1e-6, 0.01, 200)]))
    s = s[(s > 0) & (s < 1)]
    for side, pa, ca in ((+1, "parent_up", "cap_up"), (-1, "parent_lo", "cap_lo")):
        jet = P.skin_jet(s, np.full_like(s, 240.0), side)
        rep, R = join_report(F, pa, ca, jet["point"], jet["d_chord"], s)
        rep["label"] = "parent chord parameter s (0 = LE, 1 = TE corner)"
        out[f"{pa}__{ca}"] = rep
    # shared perimeter: the lip, cap_up vs cap_lo
    A = F["cap_up"]
    lip_pts, lip_t, lab = [], [], []
    for u in (A.u0, A.u1):
        for v in np.unique(np.concatenate([np.linspace(0, 1, 1001)[1:-1], 1 - np.geomspace(1e-5, 0.01, 50)])):
            if not A.inside(u, v):
                continue
            p, su, sv, *_ = A.d2(u, v)
            lip_pts.append(p)
            lip_t.append(sv)
            lab.append(v + (0 if u == A.u0 else 10))
    rep, R = join_report(F, "cap_up", "cap_lo", lip_pts, lip_t, lab)
    rep["label"] = "v on the cap patch (+10 for the trailing-edge branch); v = 1 is the crown"
    out["cap_up__cap_lo (lip)"] = rep
    # crease between cap skin and the flat trailing-edge face (continues the parent's corner)
    te = F["te"]
    nte = te.local(*te.uv((150.0, 0.0, 240.0))[:2])["n"]
    wc = cd.te_closure_height()
    crease = []
    for w in np.concatenate([np.linspace(0, wc, 60)[:-1], wc - np.geomspace(1e-6, 0.05, 20)]):
        z = cd.Z0 + w
        hw = float(cd.te_arch_halfwidth(z))
        x = 160.0 - 10.0 * z / 240.0
        u, v, g = F["cap_up"].uv((x, hw, z))
        l = F["cap_up"].local(u, v)
        crease.append((w, hw, float(np.degrees(np.arccos(np.clip(l["n"] @ nte, -1, 1))))))
    crease = np.array(crease)
    par = []
    for zz in (0.0, 120.0, 239.999):
        u, v, g = F["parent_up"].uv(P.skin(1.0, zz, 1))
        par.append(float(np.degrees(np.arccos(np.clip(F["parent_up"].local(u, v)["n"] @ nte, -1, 1)))))
    out["te_face_crease"] = {
        "note": "angle between the cap skin and the flat trailing-edge face along their shared edge; "
                "the parent has the same corner along its whole span",
        "parent_corner_angle_deg_at_z0_120_240": par,
        "cap_corner_angle_deg_at_attachment": float(crease[0, 2]),
        "cap_corner_angle_deg_near_closure (w, halfwidth, deg)": crease[np.argsort(crease[:, 0])][-5:].tolist(),
        "closure_height_mm": wc,
    }
    return out


# ------------------------------------------------------------ 2. wrap traces
def section_trace(shape, F, origin, normal, inward, q_keep, n_per_edge=1500):
    """Trace of the plane (origin, normal) through the cap. Returns ordered samples.

    inward: unit vector in the plane (and in XZ) pointing from the lip into the fin; q = (p - origin) . inward.
    Only points with q <= q_keep are kept.
    """
    sec = BRepAlgoAPI_Section(shape, gp_Pln(gp_Pnt(*origin), gp_Dir(*normal)), False)
    sec.Approximation(True)
    sec.Build()
    cands = [F[k] for k in ("cap_up", "cap_lo", "parent_up", "parent_lo", "te")]
    rows = []
    e = TopExp_Explorer(sec.Shape(), TopAbs_EDGE)
    while e.More():
        c = BRepAdaptor_Curve(TopoDS.Edge_s(e.Current()))
        for t in np.linspace(c.FirstParameter(), c.LastParameter(), n_per_edge):
            q3 = c.Value(t)
            p = np.array([q3.X(), q3.Y(), q3.Z()])
            q = (p - origin) @ inward
            if q > q_keep:
                continue
            best = None
            for G in cands:
                bb = G.bbox
                if not (bb[0] - 1e-3 <= p[0] <= bb[3] + 1e-3 and bb[1] - 1e-3 <= p[1] <= bb[4] + 1e-3
                        and bb[2] - 1e-3 <= p[2] <= bb[5] + 1e-3):
                    continue
                u, v, g = G.uv(p)
                if best is None or g < best[3]:
                    best = (G, u, v, g)
            G, u, v, g = best
            if G.name.startswith("cap") and v > 1 - 1e-5:
                continue    # collapsed crown point of the cap patches: derivatives undefined there
            if G.is_plane:
                rows.append((q, p[1], 0.0, *p, G.name, g, 0.0, 1.0))
                continue
            l = G.local(u, v)
            if l is None:
                continue
            T = np.cross(l["n"], normal)
            T /= np.linalg.norm(T)
            N0 = np.cross(normal, T)
            kn = T @ l["W"] @ T
            ks = -kn / max(abs(N0 @ l["n"]), 1e-12)
            rows.append((q, p[1], ks, *l["p"], G.name, g, -kn, abs(N0 @ l["n"])))
        e.Next()
    # order: +Y side from far to lip, then flat TE face (Y descending), then -Y side back out
    def key(r):
        q, y, ks, x, yy, z, name, g, kn, inc = r
        if name == "te":
            return (1, -y)
        return (0, -q) if y > 0 else (2, q)
    rows = sorted(rows, key=key)
    # remove duplicates at edge ends
    out, last = [], None
    for r in rows:
        pt = np.array(r[3:6])
        if last is not None and np.linalg.norm(pt - last) < 1e-9:
            continue
        out.append(r)
        last = pt
    return out


def analyse_trace(rows, label, chord_x_for_rule=None):
    q = np.array([r[0] for r in rows])
    ks = np.array([r[2] for r in rows])
    xyz = np.array([r[3:6] for r in rows])
    names = [r[6] for r in rows]
    kn = np.array([r[8] for r in rows])
    inc = np.array([r[9] for r in rows])
    on_cap = np.array([n.startswith("cap") for n in names])
    curved = np.array([n != "te" for n in names])
    # arc length
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(xyz, axis=0), axis=1))])
    R = np.where(np.abs(ks) > 0, 1 / np.maximum(np.abs(ks), 1e-300), np.inf)
    sel = curved & on_cap
    i = int(np.argmin(np.where(sel, R, np.inf)))
    lipi = int(np.argmin(np.where(on_cap, np.abs(xyz[:, 1]), np.inf)))
    # reversals: sign changes of the signed curvature on curved faces, ignoring |k| < 1e-4 /mm
    sig = np.where(curved & (np.abs(ks) > 1e-4), np.sign(ks), 0)
    nz = np.nonzero(sig)[0]
    rev = [float(s[nz[k + 1]]) for k in range(len(nz) - 1) if sig[nz[k]] != sig[nz[k + 1]]]
    concave = ks[curved] < -1e-4
    xr = xyz[i, 0] if chord_x_for_rule is None else chord_x_for_rule
    req, t = required_radius(xr)
    Rn = 1 / np.maximum(np.abs(kn), 1e-300)
    jn = int(np.argmin(np.where(sel, Rn, np.inf)))
    reqn, tn = required_radius(xyz[jn, 0])
    # progressive turn: on each side, radius should fall steadily toward the lip on the cap
    prog = {}
    for side, mask in (("+Y", xyz[:, 1] > 0), ("-Y", xyz[:, 1] < 0)):
        m = mask & on_cap & curved
        idx = np.nonzero(m)[0]
        if len(idx) < 5:
            continue
        d_lip = np.abs(s[idx] - s[lipi])
        order = np.argsort(-d_lip)
        r = R[idx][order]
        # interior local minima of radius (tight spot before the lip), ignoring the last 2 % near the lip
        mins = [k for k in range(1, len(r) - 1) if r[k] < r[k - 1] and r[k] < r[k + 1] and d_lip[order][k] > 0.02 * d_lip.max()]
        prog[side] = {"interior_radius_minima": [{"radius_mm": float(r[k]), "distance_to_lip_mm": float(d_lip[order][k])}
                                                 for k in mins][:5]}
    # monotonic tightening toward the lip on each side (cap part): count local minima of R before the lip
    return {
        "label": label, "samples": len(rows),
        "min_normal_section_radius_mm": float(Rn[jn]), "min_normal_section_at_xyz": xyz[jn].tolist(),
        "normal_section_required_mm": float(reqn), "normal_section_meets_rule": bool(Rn[jn] >= reqn * (1 - 1e-6)),
        "cut_incidence_at_min": float(inc[i]),
        "progression": prog,
        "min_radius_mm": float(R[i]), "min_radius_at_xyz": xyz[i].tolist(),
        "min_radius_face": names[i],
        "lip_radius_mm": float(R[lipi]), "lip_at_xyz": xyz[lipi].tolist(),
        "attachment_thickness_used_mm": float(t), "required_radius_mm": float(req),
        "meets_rule": bool(R[i] >= req * (1 - 1e-6)),
        "curvature_reversals_at_arclength_mm": rev[:20], "n_reversals": len(rev),
        "most_concave_curvature_per_mm": float(ks[curved].min()) if concave.any() else 0.0,
        "_profile": (s, ks, xyz, names),
    }


def wraps(shape, F):
    res = {"chord_cuts": [], "normal_cuts": []}
    # constant chord position cuts (vertical planes X = const), traced from Z = 230 around the lip
    xis = [0.0005, 0.002, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.551, 0.6, 0.7, 0.8, 0.9,
           0.95, 0.98, 0.99, 0.995, 0.998, 0.999, 0.9995, 0.9999]
    for xi in xis:
        x = 65.0 + 85.0 * xi
        rows = section_trace(shape, F, np.array([x, 0.0, 240.0]), np.array([1.0, 0, 0]),
                             np.array([0, 0, -1.0]), 10.0)
        a = analyse_trace(rows, f"X = {x:.4f} (chord fraction {xi})", chord_x_for_rule=x)
        prof = a.pop("_profile")
        a["chord_fraction"] = xi
        res["chord_cuts"].append(a)
    # cuts normal to the outline (follow the sweep): plane through the outline point, normal = outline tangent
    for branch in (-1, +1):
        for hh in [0.002, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]:
            jet = P.tip_outline_jet(np.array([hh]), branch)
            o = jet["point"][0]
            tan = jet["d_p"][0] / np.linalg.norm(jet["d_p"][0])
            inward = np.cross(np.array([0, 1.0, 0]), tan)
            # pick the in-plane direction pointing into the fin (toward the chord midline)
            if inward @ (np.array([112.0, 0, 250.0]) - o) < 0:
                inward = -inward
            rows = section_trace(shape, F, o, tan, inward, 12.0)
            a = analyse_trace(rows, f"normal to outline, {'LE' if branch < 0 else 'TE'} branch h = {hh}")
            a.pop("_profile")
            res["normal_cuts"].append(a)
    # crown: vertical plane through the crown (normal to the outline there is X)
    rows = section_trace(shape, F, np.array([112.08333333, 0, 280.0]), np.array([1.0, 0, 0]),
                         np.array([0, 0, -1.0]), 45.0)
    a = analyse_trace(rows, "crown, X = 112.0833")
    a.pop("_profile")
    res["normal_cuts"].append(a)
    return res


# ------------------------------------------------------ 3. curvature map
def curvature_map(F, nu=401, nv=801):
    out = {}
    vgrid = np.unique(np.concatenate([np.linspace(0, 1, nv), np.linspace(0, 0.03, 301),
                                      1 - np.geomspace(1e-5, 0.05, 120)]))
    for name in ("cap_up", "cap_lo"):
        G = F[name]
        ua, ub = G.u0, G.u1
        ugrid = np.unique(np.concatenate([np.linspace(ua, ub, nu), ua + np.geomspace(1e-6, 0.05, 40),
                                          ub - np.geomspace(1e-6, 0.05, 40), ub - np.linspace(0, 0.4, 120)]))
        best = []
        K = np.full((len(ugrid), len(vgrid)), np.nan)
        X = np.full_like(K, np.nan)
        for i, u in enumerate(ugrid):
            for j, v in enumerate(vgrid):
                if v > 1 - 1e-5 or not G.inside(u, v):
                    continue
                l = G.local(u, v)
                if l is None:
                    continue
                K[i, j] = l["kmax_abs"]
                X[i, j] = l["p"][0]
        req, _ = required_radius(np.nan_to_num(X, nan=100.0))
        R = 1 / K
        ratio = np.where(req > 0, R / np.where(req > 0, req, 1.0), np.inf)
        ratio = np.where(np.isfinite(K), ratio, np.nan)
        # refine the 12 smallest local minima of R/req and of R
        flat = np.argsort(np.nan_to_num(ratio, nan=np.inf).ravel())
        seeds, seen = [], []
        for f in flat:
            i, j = np.unravel_index(f, ratio.shape)
            if any(abs(i - a) < 6 and abs(j - b) < 6 for a, b in seen):
                continue
            seen.append((i, j))
            seeds.append((ugrid[i], vgrid[j]))
            if len(seeds) >= 12:
                break

        def obj(x):
            u, v = x
            if not (ua <= u <= ub and 0 <= v <= 1 - 1e-5) or not G.inside(u, v):
                return 1e3
            l = G.local(u, v)
            r = 1 / l["kmax_abs"]
            rq = float(required_radius(l["p"][0])[0])
            return r / rq if rq > 0 else 1e3

        refined = []
        for u0, v0 in seeds:
            m = minimize(obj, [u0, v0], method="Nelder-Mead",
                         options={"xatol": 1e-9, "fatol": 1e-12, "maxiter": 400})
            u, v = m.x
            l = G.local(u, v)
            r = 1 / l["kmax_abs"]
            rq, t = required_radius(l["p"][0])
            refined.append({"u": float(u), "v": float(v), "xyz": l["p"].tolist(), "radius_mm": float(r),
                            "required_mm": float(rq), "attachment_thickness_mm": float(t),
                            "ratio": float(r / rq) if rq > 0 else None})
        refined.sort(key=lambda d: d["radius_mm"])
        i, j = np.unravel_index(np.nanargmin(R), R.shape)
        crown = G.d2(1.0, 1.0)[0]
        reach = max(np.linalg.norm(G.d2(u, 1 - 1e-5)[0] - crown) for u in np.linspace(ua, ub, 41))
        out[name] = {
            "grid": [len(ugrid), len(vgrid)], "evaluated_points": int(np.isfinite(K).sum()),
            "grid_min_radius_mm": float(R[i, j]), "grid_min_at_uv": [float(ugrid[i]), float(vgrid[j])],
            "grid_min_ratio_to_required": float(np.nanmin(ratio)),
            "refined_minima": refined[:6],
            "smallest_radius_mm": refined[0]["radius_mm"] if refined[0]["radius_mm"] < R[i, j] else float(R[i, j]),
            "crown_excluded_within_mm": float(reach),
            "note": "v > 1 - 1e-5 excluded: the patch collapses to the crown point there (derivatives undefined)",
        }
        out[name + "_ratio_quantiles"] = np.nanquantile(ratio, [0, 0.001, 0.01, 0.5]).tolist()
    return out


# ---------------------------------------------------------- 4. zebra
def zebra_frame(axis):
    a = axis / np.linalg.norm(axis)
    b = np.cross(a, [0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.cross(a, [1.0, 0, 0])
    b /= np.linalg.norm(b)
    c = np.cross(a, b)
    return a, b, c


def zebra_phase(n, d, b, c, K):
    r = d - 2 * (d @ n) * n
    return K * np.arctan2(r @ c, r @ b) / (2 * np.pi), r


def zebra_grad(G, l, d, b, c, K):
    """Surface gradient (3D, in the tangent plane) of the stripe phase from CAD derivatives."""
    n = l["n"]
    r = d - 2 * (d @ n) * n
    rb, rc = r @ b, r @ c
    g = np.zeros(3)
    basis = []
    t1 = l["su"] / np.linalg.norm(l["su"])
    t2 = np.cross(n, t1)
    for t in (t1, t2):
        dn = G.dn(l, t)
        dr = -2 * ((d @ dn) * n + (d @ n) * dn)
        dphi = (rb * (dr @ c) - rc * (dr @ b)) / (rb * rb + rc * rc)
        basis.append(K * dphi / (2 * np.pi))
    return basis[0] * t1 + basis[1] * t2


def zebra(F, K=36):
    out = {}
    s = np.unique(np.concatenate([np.linspace(0, 1, 1601), np.geomspace(1e-5, 0.01, 60),
                                  1 - np.geomspace(1e-5, 0.01, 60)]))
    s = s[(s > 0) & (s < 1)]
    views, axes = [], []
    for side in (+1, -1):
        for az in (-60, -30, 0, 30, 60):
            for el in (-45, -20, 0, 20, 45):
                a, e = np.radians(az), np.radians(el)
                eye = np.array([np.sin(a) * np.cos(e), side * np.cos(a) * np.cos(e), np.sin(e)])
                views.append((side, az, el, -eye))
    for ang in (0, 30, 60, 90, 120, 150):
        axes.append(ang)
    for side, pa, ca in ((+1, "parent_up", "cap_up"), (-1, "parent_lo", "cap_lo")):
        A, C = F[pa], F[ca]
        jet = P.skin_jet(s, np.full_like(s, 240.0), side)
        locA, locC = [], []
        for p in jet["point"]:
            ua, va, _ = A.uv(p)
            uc, vc, _ = C.uv(p)
            locA.append(A.local(ua, va))
            locC.append(C.local(uc, vc))
        worst = {"phase_mismatch_stripes": (0, None), "gradient_abs_mismatch_stripes_per_mm": (0, None),
                 "gradient_rel_mismatch": (0, None),
                 "gradient_angle_deg": (0, None), "crossing_direction_change_deg": (0, None)}
        n_cross = 0
        n_eval = 0
        for side_v, az, el, d in views:
            if side_v != side:
                continue
            # stripe axes perpendicular to the view direction
            u1 = np.cross(d, [0, 0, 1.0])
            u1 /= np.linalg.norm(u1)
            u2 = np.cross(d, u1)
            for ang in axes:
                th = np.radians(ang)
                a, b, c = zebra_frame(np.cos(th) * u1 + np.sin(th) * u2)
                psiA, psiC, gA, gC, vis = [], [], [], [], []
                for la, lc in zip(locA, locC):
                    pA, _ = zebra_phase(la["n"], d, b, c, K)
                    pC, _ = zebra_phase(lc["n"], d, b, c, K)
                    psiA.append(pA)
                    psiC.append(pC)
                    gA.append(zebra_grad(A, la, d, b, c, K))
                    gC.append(zebra_grad(C, lc, d, b, c, K))
                    vis.append((d @ la["n"]) < -0.05)
                psiA, psiC = np.array(psiA), np.array(psiC)
                gA, gC = np.array(gA), np.array(gC)
                vis = np.array(vis)
                if not vis.any():
                    continue
                dpsi = np.abs((psiA - psiC + K / 2) % K - K / 2)
                gn = np.maximum(np.linalg.norm(gA, axis=1), 1e-12)
                absd = np.linalg.norm(gA - gC, axis=1)
                gmed = float(np.median(gn[vis]))
                rel = np.where(gn > 0.1 * gmed, absd / gn, 0.0)
                cosang = np.einsum("ij,ij->i", gA, gC) / (gn * np.maximum(np.linalg.norm(gC, axis=1), 1e-12))
                gang = np.degrees(np.arccos(np.clip(cosang, -1, 1)))
                n_eval += int(vis.sum())
                for key, arr in (("phase_mismatch_stripes", dpsi), ("gradient_abs_mismatch_stripes_per_mm", absd),
                                 ("gradient_rel_mismatch", rel),
                                 ("gradient_angle_deg", gang)):
                    m = np.where(vis, arr, -1)
                    k = int(np.argmax(m))
                    if m[k] > worst[key][0]:
                        worst[key] = (float(m[k]), {"s": float(s[k]), "view_az_el": [az, el], "stripe_axis_deg": ang,
                                                    "xyz": jet["point"][k].tolist(), "gradient_magnitude": float(gn[k]),
                                                    "median_gradient_this_view": gmed})
                # light/dark crossings on the join: psi crosses a multiple of 0.5
                ph = np.unwrap(psiA * 2 * np.pi / K) * K / (2 * np.pi)
                cells = np.floor(2 * ph)
                for k in np.nonzero(np.diff(cells) != 0)[0]:
                    if not (vis[k] and vis[k + 1]):
                        continue
                    n_cross += 1
                    # stripe direction on each face at the crossing: perpendicular to the phase gradient
                    for_a = np.cross(locA[k]["n"], gA[k])
                    for_c = np.cross(locC[k]["n"], gC[k])
                    cang = np.degrees(np.arccos(np.clip(abs(for_a @ for_c) /
                                                        (np.linalg.norm(for_a) * np.linalg.norm(for_c) + 1e-300), 0, 1)))
                    if cang > worst["crossing_direction_change_deg"][0]:
                        worst["crossing_direction_change_deg"] = (float(cang), {
                            "s": float(s[k]), "view_az_el": [az, el], "stripe_axis_deg": ang,
                            "xyz": jet["point"][k].tolist()})
        out[f"{pa}__{ca}"] = {"views": len([v for v in views if v[0] == side]), "stripe_axes": len(axes),
                              "stripes_per_turn": K, "visible_samples": n_eval, "crossings_checked": n_cross,
                              "worst": {k: {"value": v[0], "where": v[1]} for k, v in worst.items()},
                              "sampled_region": "attachment join, parent chord parameter s in [1e-5, 1 - 1e-5]",
                              "undefined_endpoints": "none on this join; s = 0 and 1 (LE point, TE corner) not sampled"}
    return out


def summarise(rb):
    """Put a short requirement-by-requirement summary first."""
    m = rb["curvature_map"]
    jn = rb["joins"]
    ws = rb["wrap_sections"]
    ld = rb["locked_dimensions"]
    worst_map = min((m[k]["refined_minima"][0] for k in ("cap_up", "cap_lo")), key=lambda d: d["radius_mm"])
    normal_min = min(ws["normal_cuts"], key=lambda a: a["min_radius_mm"])
    chord_nmin = min(ws["chord_cuts"], key=lambda a: a["min_normal_section_radius_mm"])
    chord_smin = min(ws["chord_cuts"], key=lambda a: a["min_radius_mm"])
    summary = {
        "closed_valid_single_solid": rb["reimport"]["solids"] == 1 and rb["reimport"]["brepcheck_valid"]
        and rb["reimport"]["brepcheck_exact_valid"] and rb["reimport"]["free_edges_none"],
        "parent_skin_preserved_max_mm": max(rb["parent_preserved"]["parent_up_max_distance_to_parent_py_mm"],
                                            rb["parent_preserved"]["parent_lo_max_distance_to_parent_py_mm"]),
        "finished_height_mm": ld["finished_height_from_crown_point_mm"],
        "outline_max_deviation_mm": max(ld["outline_max_deviation_section_to_outline_mm"],
                                        ld["outline_max_deviation_outline_to_section_mm"]),
        "outline_max_deviation_where_lip_is_silhouette_mm": ld["outline_max_deviation_where_lip_forms_silhouette_mm (z > closure)"],
        "attachment_join_max_separation_mm": max(jn[k]["max_separation_mm"] for k in ("parent_up__cap_up", "parent_lo__cap_lo")),
        "attachment_join_max_normal_angle_deg": max(jn[k]["max_normal_angle_deg"] for k in ("parent_up__cap_up", "parent_lo__cap_lo")),
        "attachment_join_max_curvature_difference_per_mm": max(max(jn[k]["max_abs_curvature_difference_per_mm"].values())
                                                               for k in ("parent_up__cap_up", "parent_lo__cap_lo")),
        "printable_rule": {
            "smallest_principal_radius_mm": worst_map["radius_mm"], "at_xyz": worst_map["xyz"],
            "attachment_thickness_there_mm": worst_map["attachment_thickness_mm"],
            "required_there_mm": worst_map["required_mm"],
            "smallest_ratio_radius_over_required": min(m[k]["grid_min_ratio_to_required"] for k in ("cap_up", "cap_lo")),
            "met": worst_map["radius_mm"] >= worst_map["required_mm"],
            "smallest_normal_cut_section_radius_mm": normal_min["min_radius_mm"], "normal_cut": normal_min["label"],
            "smallest_normal_section_radius_on_chord_cuts_mm": chord_nmin["min_normal_section_radius_mm"],
            "smallest_plane_section_radius_on_chord_cuts_mm": chord_smin["min_radius_mm"],
            "chord_cut_note": "constant-X cuts near the LE and TE cross the lip at a grazing angle (incidence "
                              f"{chord_smin['cut_incidence_at_min']:.3f} at {chord_smin['label']}); their plane-section radius "
                              "is the oblique slice of a ~0.6 mm round, not a surface radius",
        },
        "zebra_attachment_worst": {k: v["value"] for k, v in rb["zebra"]["parent_up__cap_up"]["worst"].items()},
    }
    return {"summary": summary, **rb}


def clean(o):
    """JSON-safe copy: NaN and infinities become null, numpy scalars become floats."""
    if isinstance(o, dict):
        return {str(k): clean(x) for k, x in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(x) for x in o]
    if isinstance(o, (np.floating, float)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def main():
    shape = ou.read_step(STEP)
    F = ou.classify(shape)
    rb = {"file": "runs/cap-radius/builds/fin40/blunt/wing.step", "step_sha256": ou.sha256(STEP),
          "step_bytes": STEP.stat().st_size, "units": "mm"}
    log("topology")
    rb["reimport"] = topology(shape, F)
    log("parent")
    rb["parent_preserved"] = parent_preserved(F)
    rb["parent_preserved"]["wing_vs_parent_step_skin_max_mm"] = parent_step_identity(F)
    log("locked")
    rb["locked_dimensions"] = locked_dimensions(shape, F)
    log("joins")
    rb["joins"] = joins(F)
    log("wraps")
    rb["wrap_sections"] = wraps(shape, F)
    log("map")
    rb["curvature_map"] = curvature_map(F)
    log("zebra")
    rb["zebra"] = zebra(F)
    rb = summarise(rb)
    OUT.write_text(json.dumps(clean(rb), indent=1))
    log("written", OUT)


if __name__ == "__main__":
    main()
