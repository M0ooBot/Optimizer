"""Verify the reimported wing.step (Wing B sharp, rounded tip cap) and write readback.json.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/wingb/sharp/verify_cap.py [--quick]

Checks (on the reimported STEP only, CAD normals and derivatives):
  1. topology / validity / volume / SHA-256; parent skin and root face preserved;
     locked span, attachment station and outline;
  2. wing-to-cap join: separation, tangent-plane angle, curvature in a common frame;
  3. sections (guide planes, between guides, spanwise, sweep-following, chordwise):
     signed curvature, smallest radius, lip radius, bending reversals, printable rule;
  4. map of the larger principal curvature over the whole cap with refinement;
  5. numerical zebra across the join (phase, gradient, stripe direction at crossings).
"""
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize

sys.path.insert(0, str(Path(__file__).resolve().parent))
import capgeom as g  # noqa: E402
import surfeval as se  # noqa: E402
from build_cap import p_of_t  # noqa: E402

QUICK = "--quick" in sys.argv
STEP = g.BUILD_DIR / "wing.step"
RULE_MIN = 0.4
RULE_T = 0.8


def required_radius(xi):
    return np.minimum(RULE_MIN, g.thickness_att(xi) / 2.0)


# ----------------------------------------------------------------------------
def read_step(path):
    from OCP.BRep import BRep_Tool
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.BRepGProp import BRepGProp
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    from OCP.GProp import GProp_GProps
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPControl import STEPControl_Reader
    from OCP.ShapeAnalysis import ShapeAnalysis_FreeBounds
    from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID, TopAbs_EDGE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    from OCP.TopTools import TopTools_HSequenceOfShape

    r = STEPControl_Reader()
    assert r.ReadFile(str(path)) == IFSelect_RetDone
    r.TransferRoots()
    shape = r.OneShape()

    def count(kind):
        e = TopExp_Explorer(shape, kind)
        n = 0
        while e.More():
            n += 1
            e.Next()
        return n
    faces = []
    e = TopExp_Explorer(shape, TopAbs_FACE)
    while e.More():
        f = TopoDS.Face_s(e.Current())
        a = BRepAdaptor_Surface(f)
        bb = Bnd_Box()
        BRepBndLib.Add_s(f, bb)
        faces.append(dict(face=f, type=str(a.GetType()).split(".")[-1], surf=BRep_Tool.Surface_s(f),
                          box=bb.Get(), adaptor=a))
        e.Next()
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props, 1e-9)
    bb = Bnd_Box()
    BRepBndLib.Add_s(shape, bb, False)
    fb = ShapeAnalysis_FreeBounds(shape)
    free = fb.GetClosedWires()
    free_open = fb.GetOpenWires()
    nfree = 0
    for w in (free, free_open):
        ee = TopExp_Explorer(w, TopAbs_EDGE)
        while ee.More():
            nfree += 1
            ee.Next()
    info = dict(solids=count(TopAbs_SOLID), shells=count(TopAbs_SHELL), faces=len(faces),
                face_types=[f["type"] for f in faces], brepcheck_valid=bool(BRepCheck_Analyzer(shape).IsValid()),
                free_boundary_edges=nfree, volume_mm3=props.Mass(), bounding_box_mm=list(bb.Get()))
    return shape, faces, info


def classify(faces):
    cap = skin = root = None
    for f in faces:
        x0, y0, z0, x1, y1, z1 = f["box"]
        if f["type"] == "GeomAbs_BSplineSurface" and x0 < -590:
            cap = f
        elif f["type"] == "GeomAbs_BSplineSurface":
            skin = f
        elif f["type"] == "GeomAbs_Plane":
            root = f
    return cap, skin, root


# ----------------------------------------------------------------------------
def check_parent_preserved(skin, root):
    par = g.parent()
    ps = par.surfaces()[0]
    s = skin["surf"]
    same_shape = (s.NbUPoles() == ps.NbUPoles() and s.NbVPoles() == ps.NbVPoles())
    dmax = 0.0
    if same_shape:
        for i in range(1, s.NbUPoles() + 1):
            for j in range(1, s.NbVPoles() + 1):
                dmax = max(dmax, s.Pole(i, j).Distance(ps.Pole(i, j)))
    uu = np.linspace(0, 1, 97)
    vv = np.linspace(0, par.SURFACE_DATA["v_attachment"], 41)
    dev = max(s.Value(float(u), float(v)).Distance(ps.Value(float(u), float(v))) for u in uu for v in vv)
    b = skin["adaptor"]
    pl = root["adaptor"].Plane()
    return dict(skin_poles_identical_shape=bool(same_shape), skin_max_pole_difference_mm=dmax,
                skin_max_sample_deviation_mm=dev,
                skin_parameter_range=[b.FirstUParameter(), b.LastUParameter(), b.FirstVParameter(), b.LastVParameter()],
                skin_v_attachment=par.SURFACE_DATA["v_attachment"],
                root_plane_origin=list(pl.Location().Coord()), root_plane_normal=list(pl.Axis().Direction().Coord()))


# ----------------------------------------------------------------------------
def tabulate(cap, nu, nv):
    t = np.linspace(0, 1, nu)
    t[0], t[-1] = 1e-4, 1 - 1e-4
    v = np.linspace(0, 2, nv)
    U, V = np.meshgrid(t, v, indexing="ij")
    D = se.d2(cap, U, V)
    C = se.curvature(D)
    # outward orientation: at the lip of the middle guide the normal points to -X (outboard)
    k = nu // 2
    j = nv // 2
    sgn = -np.sign(C["n"][k, j, 0])
    C["n"] = C["n"] * sgn
    C["k1"], C["k2"] = C["k1"] * sgn, C["k2"] * sgn
    C["LMN"] = tuple(x * sgn for x in C["LMN"])
    C["sign"] = sgn
    return U, V, D, C


def check_locked(cap, U, V, P):
    par = g.parent()
    X = P[..., 0]
    # max projected span of the cap (refine on the middle rows)
    i, j = np.unravel_index(np.argmin(X), X.shape)
    res = minimize(lambda x: cap.Value(float(np.clip(x[0], 0, 1)), float(np.clip(x[1], 0, 2))).X(),
                   [U[i, j], V[i, j]], method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-12))
    span = -res.fun
    # plan silhouette vs locked outline: signed distance of every sample from the outline
    pp = np.linspace(0, 1, 20001)
    O = par.tip_outline(pp)
    from scipy.spatial import cKDTree
    tree = cKDTree(O)
    pts = P[..., [0, 2]].reshape(-1, 2)
    dist, idx = tree.query(pts)
    T = np.gradient(O, axis=0)
    T /= np.linalg.norm(T, axis=1, keepdims=True)
    Nout = np.stack([T[:, 1], -T[:, 0]], -1)          # rotate tangent; fix orientation below
    mid = len(pp) // 2
    if Nout[mid, 0] > 0:                                 # outward at the tip is -X
        Nout = -Nout
    signed = np.sum((pts - O[idx]) * Nout[idx], 1)
    # lip trace: per guide row, the plan-outermost point along the outline normal
    lip_dev = []
    for r in range(1, U.shape[0] - 1, 4):
        pr = P[r][:, [0, 2]]
        d_, id_ = tree.query(pr)
        s_ = np.sum((pr - O[id_]) * Nout[id_], 1)
        lip_dev.append(np.max(s_))
    return dict(max_projected_span_mm=span, span_location=list(cap.Value(*map(float, res.x)).Coord()),
                max_outside_outline_mm=float(signed.max()),
                lip_to_outline_max_abs_mm=float(np.max(np.abs(lip_dev))),
                attachment_plane_x_spread_mm=float(np.ptp(np.r_[P[:, 0, 0], P[:, -1, 0]]) if False else
                                                   max(abs(P[:, 0, 0] + 540).max(), abs(P[:, -1, 0] + 540).max())))


# ----------------------------------------------------------------------------
def project(surface, pts, guess=None):
    from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
    from OCP.gp import gp_Pnt
    out = []
    for i, p in enumerate(pts):
        pr = GeomAPI_ProjectPointOnSurf(gp_Pnt(*map(float, p)), surface)
        if pr.NbPoints() == 0:
            out.append((np.nan, np.nan, np.nan))
            continue
        u, v = pr.LowerDistanceParameters()
        out.append((u, v, pr.LowerDistance()))
    return np.array(out)


def frame_curvatures(C, e1, n_ref):
    """Normal curvatures (II11, II22, II12) in the frame e1 (along the join), e2 = n x e1, n oriented like n_ref."""
    n = C["n"] * np.sign(np.sum(C["n"] * n_ref, -1))[..., None]
    sgn = np.sign(np.sum(C["n"] * n_ref, -1))
    e1 = e1 - np.sum(e1 * n, -1)[..., None] * n
    e1 /= np.linalg.norm(e1, axis=-1, keepdims=True)
    e2 = np.cross(n, e1)
    C2 = dict(C)
    C2["LMN"] = tuple(x * sgn for x in C["LMN"])
    k11 = se.normal_curvature(C2, e1)
    k22 = se.normal_curvature(C2, e2)
    kp = se.normal_curvature(C2, (e1 + e2) / np.sqrt(2))
    k12 = kp - 0.5 * (k11 + k22)
    return n, e1, e2, np.stack([k11, k22, k12], -1)


def check_join(cap, skin, nsamp):
    tab = g.attachment_tables()
    tt = np.linspace(0, 1, nsamp)
    xi = tab["xi_nose"] + (1 - tab["xi_nose"]) * (1 - np.cos(np.pi * tt)) / 2
    xi = xi[(xi > tab["xi_nose"] + 2e-4) & (xi < 1 - 2e-4)]
    jet = g.attachment_jet(xi)
    X, Z = g.plan_xz(np.full_like(xi, 540.0), xi)
    res = {}
    for side, name in ((1, "upper"), (-1, "lower")):
        Pp = np.stack([X, jet[side][0], Z], -1)
        e1 = np.gradient(Pp, axis=0)
        sk = project(skin["surf"], Pp)
        cp = project(cap, Pp)
        Dp = se.d2(skin["surf"], sk[:, 0], sk[:, 1])
        Dc = se.d2(cap, cp[:, 0], cp[:, 1])
        Cp, Cc = se.curvature(Dp), se.curvature(Dc)
        n_ref = np.zeros_like(Pp)
        n_ref[:, 1] = side
        npar, e1p, e2p, Kp = frame_curvatures(Cp, e1, n_ref)
        ncap, e1c, e2c, Kc = frame_curvatures(Cc, e1, npar)
        ang = np.degrees(np.arccos(np.clip(np.sum(npar * ncap, -1), -1, 1)))
        sep = np.linalg.norm(Dp[:, 0] - Dc[:, 0], axis=1)
        dK = np.abs(Kp - Kc)
        scale = np.maximum(np.abs(Kp).max(-1), 1e-3)
        i_a, i_k = int(np.argmax(ang)), int(np.argmax(dK[:, 1]))
        res[name] = dict(samples=int(len(xi)), xi_range=[float(xi[0]), float(xi[-1])],
                         max_separation_mm=float(sep.max()), max_skin_to_edge_mm=float(sk[:, 2].max()),
                         max_cap_to_edge_mm=float(cp[:, 2].max()),
                         max_tangent_angle_deg=float(ang.max()), at_xi_tangent=float(xi[i_a]),
                         max_dk_cross_per_mm=float(dK[:, 1].max()), at_xi_dk_cross=float(xi[i_k]),
                         parent_k_cross_there=float(Kp[i_k, 1]), cap_k_cross_there=float(Kc[i_k, 1]),
                         max_dk_along_per_mm=float(dK[:, 0].max()), max_dk_twist_per_mm=float(dK[:, 2].max()),
                         max_rel_dk=float(np.max(dK.max(-1) / scale)),
                         rms_dk_cross_per_mm=float(np.sqrt(np.mean(dK[:, 1] ** 2))))
        res[name]["_arrays"] = dict(xi=xi, P=Dp[:, 0], npar=npar, ncap=ncap, e1=e1p, e2=e2p, Kp=Kp, Kc=Kc)
    return res


# ----------------------------------------------------------------------------
def zebra_join(join, nstripes=24, exclude_mm=0.0, LE=None, TE=None):
    """Stripe phase phi = nstripes * acos(r.a)/pi of the reflected view ray, on both sides of the join."""
    views = []
    for el in (-60, -30, 0, 30, 60):
        for az in range(0, 360, 30):
            a, e = np.radians(az), np.radians(el)
            views.append(np.array([np.cos(e) * np.cos(a), np.sin(e), np.cos(e) * np.sin(a)]))
    axes = [np.array(x, float) / np.linalg.norm(x) for x in
            ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [0, 1, 1], [1, 0, 1], [1, -1, 0], [0, 1, -1])]
    worst = dict(phase=(0, None), grad_rel=(0, None), dir_deg=(0, None))
    ncross = 0
    for side, J in join.items():
        A = J["_arrays"]
        for d in views:
            for a in axes:
                ph, gr = [], []
                for n, K in ((A["npar"], A["Kp"]), (A["ncap"], A["Kc"])):
                    dn_ = np.sum(d * n, -1)
                    r = d - 2 * dn_[:, None] * n
                    ra = np.clip(np.sum(r * a, -1), -1, 1)
                    phi = nstripes * np.arccos(ra) / np.pi
                    # gradient: dn(t) = -(II(t,e1) e1 + II(t,e2) e2)
                    grads = []
                    for (c1, c2) in ((1.0, 0.0), (0.0, 1.0)):
                        s1 = K[:, 0] * c1 + K[:, 2] * c2
                        s2 = K[:, 2] * c1 + K[:, 1] * c2
                        dn = -(s1[:, None] * A["e1"] + s2[:, None] * A["e2"])
                        dr = -2 * (np.sum(d * dn, -1)[:, None] * n + dn_[:, None] * dn)
                        grads.append(-nstripes / np.pi * np.sum(a * dr, -1) / np.sqrt(np.maximum(1 - ra**2, 1e-12)))
                    ph.append(phi)
                    gr.append(np.stack(grads, -1))
                vis = (np.sum(d * A["npar"], -1) < -0.1) & (np.abs(np.sum(np.stack([d] * len(A["npar"])) *
                                                                   (d[None] - 0), -1)) > 0)
                ok = vis & (np.abs(np.sum((d - 2 * np.sum(d * A["npar"], -1)[:, None] * A["npar"]) * a, -1)) < 0.98)
                if not ok.any():
                    continue
                if exclude_mm > 0:
                    far = (np.linalg.norm(A["P"] - LE, axis=-1) > exclude_mm) & (np.linalg.norm(A["P"] - TE, axis=-1) > exclude_mm)
                    ok = ok & far
                    if not ok.any():
                        continue
                dph = np.abs(ph[0] - ph[1])
                i = np.argmax(np.where(ok, dph, -1))
                if dph[i] > worst["phase"][0]:
                    worst["phase"] = (float(dph[i]), dict(side=side, xi=float(A["xi"][i]), view=d.tolist(), axis=a.tolist()))
                # stripe edges crossing the join: phi passes a half-integer between neighbours
                f = ph[0] - 0.5
                cross = np.where(ok[:-1] & ok[1:] & (np.floor(f[:-1]) != np.floor(f[1:])))[0]
                for c in cross:
                    ncross += 1
                    ga, gb = gr[0][c], gr[1][c]
                    absd = float(np.linalg.norm(ga - gb))
                    if absd > worst.setdefault("abs", (0, None))[0]:
                        worst["abs"] = (absd, dict(side=side, xi=float(A["xi"][c]), view=[round(x, 4) for x in d],
                                                   axis=[round(x, 4) for x in a], grad_parent_per_mm=float(np.linalg.norm(ga)),
                                                   grad_cap_per_mm=float(np.linalg.norm(gb))))
                    if min(np.linalg.norm(ga), np.linalg.norm(gb)) < 0.2:
                        continue  # bands wider than 5 mm: stripe direction is not meaningful
                    ang = np.degrees(np.arccos(np.clip(np.dot(ga, gb) / (np.linalg.norm(ga) * np.linalg.norm(gb) + 1e-30), -1, 1)))
                    # relative to the larger gradient; gradients below 0.05 stripe/mm (bands wider than 20 mm) are
                    # floored so that near-flat reflections do not dominate
                    rel = np.linalg.norm(ga - gb) / max(np.linalg.norm(ga), np.linalg.norm(gb), 0.05)
                    loc = dict(side=side, xi=float(A["xi"][c]), view=[round(x, 4) for x in d], axis=[round(x, 4) for x in a],
                               grad_parent_per_mm=float(np.linalg.norm(ga)), grad_cap_per_mm=float(np.linalg.norm(gb)))
                    if ang > worst["dir_deg"][0]:
                        worst["dir_deg"] = (float(ang), loc)
                    if rel > worst["grad_rel"][0]:
                        worst["grad_rel"] = (float(rel), loc)
    return dict(stripes_per_180deg=nstripes, views=len(views), stripe_axes=len(axes), crossings_checked=ncross,
                excluded_within_mm_of_collapsed_points=exclude_mm,
                max_phase_mismatch_stripe_widths=worst["phase"][0], at_phase=worst["phase"][1],
                max_direction_change_at_crossing_deg=worst["dir_deg"][0], at_direction=worst["dir_deg"][1],
                max_relative_gradient_change_at_crossing=worst["grad_rel"][0], at_gradient=worst["grad_rel"][1],
                max_abs_gradient_change_at_crossing_stripes_per_mm=worst.get("abs", (0, None))[0],
                at_abs_gradient=worst.get("abs", (0, None))[1],
                note="direction and relative-gradient changes are evaluated at crossings where both stripe gradients "
                     "are at least 0.2 stripe/mm (bands narrower than 5 mm)")


# ----------------------------------------------------------------------------
def principal_map(cap, U, V, C):
    P = C["P"]
    xi = g.xi_of(-P[..., 0], P[..., 2])
    req = required_radius(xi)
    R = 1.0 / np.maximum(C["kmax"], 1e-300)
    margin = R / req
    LE = np.array(cap.Value(0.0, 1.0).Coord())
    TE = np.array(cap.Value(1.0, 1.0).Coord())
    dLE = np.linalg.norm(P - LE, axis=-1)
    dTE = np.linalg.norm(P - TE, axis=-1)

    def ratio_at(x):
        u, v = float(np.clip(x[0], 1e-6, 1 - 1e-6)), float(np.clip(x[1], 0, 2))
        D = se.d2(cap, np.array([u]), np.array([v]))
        Cx = se.curvature(D)
        Px = Cx["P"][0]
        r = 1 / Cx["kmax"][0]
        q = required_radius(g.xi_of(-Px[0], Px[2]))
        return r / q, r, q, Px
    # refine the worst local minima of R/req
    from scipy.ndimage import minimum_filter
    loc = (margin == minimum_filter(margin, size=7, mode="nearest"))
    cand = np.argwhere(loc)
    order = np.argsort(margin[loc])[:40]
    refined = []
    for idx in cand[order]:
        i, j = idx
        x0 = [U[i, j], V[i, j]]
        r = minimize(lambda x: ratio_at(x)[0], x0, method="Nelder-Mead",
                     options=dict(xatol=1e-9, fatol=1e-10, maxiter=400,
                                  initial_simplex=[x0, [x0[0] + 2e-4, x0[1]], [x0[0], x0[1] + 2e-3]]))
        rr, Rr, q, Px = ratio_at(r.x)
        refined.append(dict(u=float(r.x[0]), v=float(r.x[1]), p=float(p_of_t(r.x[0])), ratio=float(rr),
                            radius_mm=float(Rr), required_mm=float(q),
                            xi=float(g.xi_of(-Px[0], Px[2])), t_att_mm=float(g.thickness_att(g.xi_of(-Px[0], Px[2]))),
                            point=[float(c) for c in Px], d_LE_mm=float(np.linalg.norm(Px - LE)),
                            d_TE_mm=float(np.linalg.norm(Px - TE))))
    refined.sort(key=lambda d: d["ratio"])
    # smallest radius among points where the strict rule applies (t_att >= 0.8)
    strict = g.thickness_att(xi) >= RULE_T
    i = np.unravel_index(np.argmin(np.where(strict, R, np.inf)), R.shape)
    strict_min = dict(radius_mm=float(R[i]), p=float(p_of_t(U[i])), v=float(V[i]), xi=float(xi[i]),
                      t_att_mm=float(g.thickness_att(xi[i])), point=P[i].tolist())
    rs = minimize(lambda x: ratio_at(x)[1] if g.thickness_att(g.xi_of(-ratio_at(x)[3][0], ratio_at(x)[3][2])) >= RULE_T else 9,
                  [U[i], V[i]], method="Nelder-Mead", options=dict(xatol=1e-9, fatol=1e-10))
    _, Rr, q, Px = ratio_at(rs.x)
    strict_min_refined = dict(radius_mm=float(Rr), p=float(p_of_t(rs.x[0])), v=float(rs.x[1]),
                              xi=float(g.xi_of(-Px[0], Px[2])), t_att_mm=float(g.thickness_att(g.xi_of(-Px[0], Px[2]))),
                              point=[float(c) for c in Px])
    bands = []
    for lo, hi in ((0.8, 99), (0.4, 0.8), (0.1, 0.4), (0.02, 0.1), (0.005, 0.02), (0.0, 0.005)):
        t = g.thickness_att(xi)
        m = (t >= lo) & (t < hi)
        if m.any():
            k = np.unravel_index(np.argmin(np.where(m, margin, np.inf)), margin.shape)
            bands.append(dict(t_att_band_mm=[lo, hi], samples=int(m.sum()), min_ratio_R_over_required=float(margin[k]),
                              radius_mm=float(R[k]), required_mm=float(req[k]), p=float(p_of_t(U[k])), v=float(V[k]),
                              d_LE_mm=float(dLE[k]), d_TE_mm=float(dTE[k])))
    excl = []
    for lim in (0.0, 0.01, 0.02, 0.05, 0.1, 0.5):
        m = (dLE > lim) & (dTE > lim)
        k = np.unravel_index(np.argmin(np.where(m, margin, np.inf)), margin.shape)
        excl.append(dict(exclude_within_mm_of_collapsed_points=lim, min_ratio=float(margin[k]),
                         radius_mm=float(R[k]), required_mm=float(req[k]), xi=float(xi[k]), d_LE_mm=float(dLE[k]),
                         d_TE_mm=float(dTE[k])))
    fail = margin < 1.0
    tt = g.thickness_att(xi)
    nearLE = fail & (dLE < 1.0)
    rest = fail & ~nearLE
    fail_zone = dict(failing_samples=int(fail.sum()), of=int(fail.size),
                     near_LE_point_failing=int(nearLE.sum()),
                     near_LE_max_distance_mm=float(dLE[nearLE].max()) if nearLE.any() else None,
                     near_LE_min_ratio=float(margin[nearLE].min()) if nearLE.any() else None,
                     elsewhere_failing=int(rest.sum()),
                     elsewhere_max_t_att_mm=float(tt[rest].max()) if rest.any() else None,
                     elsewhere_max_required_mm=float(req[rest].max()) if rest.any() else None,
                     elsewhere_min_ratio=float(margin[rest].min()) if rest.any() else None,
                     elsewhere_max_d_TE_mm=float(dTE[rest].max()) if rest.any() else None)
    smallest = np.unravel_index(np.argmin(R), R.shape)
    return dict(grid=list(U.shape), u_range=[float(U.min()), float(U.max())],
                grid_closest_to_LE_point_mm=float(dLE.min()), grid_closest_to_TE_point_mm=float(dTE.min()),
                smallest_radius_anywhere_mm=float(R[smallest]),
                smallest_radius_anywhere_at=dict(p=float(p_of_t(U[smallest])), v=float(V[smallest]),
                                                 d_TE_mm=float(dTE[smallest]), d_LE_mm=float(dLE[smallest])),
                strict_zone_smallest_radius=strict_min, strict_zone_smallest_radius_refined=strict_min_refined,
                by_attachment_thickness=bands, exclusion_study=excl, failing_zone=fail_zone, worst_refined_minima=refined[:12])


# ----------------------------------------------------------------------------
def _root_refine(cap, n_p, p0, u, v, along):
    """Refine a plane crossing along u (along='u', v fixed) or along v (u fixed)."""
    def f(x):
        if along == "u":
            P = cap.Value(float(x), float(v))
        else:
            P = cap.Value(float(u), float(x))
        return n_p[0] * (P.X() - p0[0]) + n_p[1] * (P.Y() - p0[1]) + n_p[2] * (P.Z() - p0[2])
    return f


def trace_plane(cap, U, V, Pgrid, n_p, p0, both=False):
    """Points of the cap on the plane n_p.(x - p0) = 0, as (u, v) pairs (sign changes refined by brentq)."""
    F = np.tensordot(Pgrid - p0, n_p, axes=([-1], [0]))
    uv = []
    # crossings along u for every v column
    for j in range(F.shape[1]):
        col = F[:, j]
        idx = np.where(np.sign(col[:-1]) * np.sign(col[1:]) < 0)[0]
        for i in idx:
            f = _root_refine(cap, n_p, p0, None, V[0, j], "u")
            uv.append((brentq(f, U[i, j], U[i + 1, j], xtol=1e-14), V[0, j]))
    if both:
        for i in range(F.shape[0]):
            row = F[i, :]
            idx = np.where(np.sign(row[:-1]) * np.sign(row[1:]) < 0)[0]
            for j in idx:
                f = _root_refine(cap, n_p, p0, U[i, 0], None, "v")
                uv.append((U[i, 0], brentq(f, V[i, j], V[i, j + 1], xtol=1e-14)))
    return np.array(uv)


def section_curvature(cap, uv, n_p, sgn):
    D = se.d2(cap, uv[:, 0], uv[:, 1])
    C = se.curvature(D)
    n = C["n"] * sgn
    LMN = tuple(x * sgn for x in C["LMN"])
    C2 = dict(C, n=n, LMN=LMN)
    t = np.cross(n, n_p)
    t /= np.linalg.norm(t, axis=-1, keepdims=True)
    b = np.cross(n_p, t)
    kn = se.normal_curvature(C2, t)
    bn = np.sum(b * n, -1)
    c = -kn / np.maximum(np.abs(bn), 1e-12)          # signed plane-curve curvature, + convex
    return C["P"], n, c, -kn                         # -kn: surface normal curvature along the trace, + convex


def analyse_trace(name, P, n, c, cn, lip_dir=None, lip_pts=None):
    """c: plane-curve signed curvature of the trace; cn: surface normal curvature along the trace.

    Radii and the printable rule use the surface normal curvature along the trace (the radius the
    surface turns with in that direction). The plane-curve radius of an oblique cut through a lip is
    smaller than the surface radius by the cut's obliquity, so it is recorded separately.
    """
    xi = g.xi_of(-P[:, 0], P[:, 2])
    t_att = g.thickness_att(xi)
    req = required_radius(xi)
    R = 1 / np.maximum(np.abs(cn), 1e-300)
    Rplane = 1 / np.maximum(np.abs(c), 1e-300)
    ratio = R / req
    i = int(np.argmin(R))
    k = int(np.argmin(ratio))
    out = dict(name=name, points=int(len(P)), smallest_radius_mm=float(R[i]),
               smallest_at=dict(point=[round(float(x), 4) for x in P[i]], xi=float(xi[i]), t_att_mm=float(t_att[i]),
                                required_mm=float(req[i]), side="upper" if n[i, 1] >= 0 else "lower",
                                span_mm=float(-P[i, 0])),
               smallest_plane_curve_radius_mm=float(Rplane.min()),
               worst_ratio=float(ratio[k]),
               worst_at=dict(point=[round(float(x), 4) for x in P[k]], radius_mm=float(R[k]), xi=float(xi[k]),
                             t_att_mm=float(t_att[k]), required_mm=float(req[k])))
    if lip_dir is not None:
        j = int(np.argmax(n @ lip_dir))
        out["lip_radius_mm"] = float(R[j])
        out["lip_point"] = [round(float(x), 4) for x in P[j]]
        out["lip_xi"] = float(xi[j])
        out["lip_t_att_mm"] = float(t_att[j])
    if lip_pts is not None:
        out["lips"] = []
        for j in lip_pts:
            out["lips"].append(dict(point=[round(float(x), 4) for x in P[j]], radius_mm=float(R[j]), xi=float(xi[j]),
                                    t_att_mm=float(t_att[j])))
    # bending reversals: sign changes of the curvature, ignoring |c| < 2e-3 /mm (radius > 500 mm)
    # and requiring each signed run to last at least 3 samples
    s = np.where(cn > 2e-3, 1, np.where(cn < -2e-3, -1, 0))
    runs = []
    cur, start = 0, 0
    for i_, si in enumerate(s):
        if si != 0 and si != cur:
            if cur != 0:
                runs.append((cur, start, i_ - 1))
            cur, start = si, i_
    if cur != 0:
        runs.append((cur, start, len(s) - 1))
    runs = [r for r in runs if np.sum(s[r[1]:r[2] + 1] == r[0]) >= 3]
    rev = []
    for r0, r1 in zip(runs[:-1], runs[1:]):
        if r0[0] != r1[0]:
            m = (r0[2] + r1[1]) // 2
            rev.append(dict(span_mm=round(float(-P[m, 0]), 3), xi=round(float(xi[m]), 5),
                            side="upper" if n[m, 1] >= 0 else "lower", y=round(float(P[m, 1]), 4),
                            to="concave" if r1[0] < 0 else "convex"))
    out["bending_reversals"] = len(rev)
    out["reversal_locations"] = rev[:8]
    out["concave_min_radius_mm"] = float(1 / max(np.max(-cn), 1e-300)) if np.any(cn < 0) else None
    out["pass"] = bool(ratio[k] >= 1.0 - 1e-9)
    return out


def sections(cap, U, V, Pgrid, sgn):
    par = g.parent()
    res = []
    Y = np.array([0.0, 1.0, 0.0])

    def run(name, n_p, p0, lip_dir=None, chordwise=False):
        n_p = n_p / np.linalg.norm(n_p)
        uv = trace_plane(cap, U, V, Pgrid, n_p, p0, both=chordwise)
        if len(uv) < 5:
            return
        P, n, c, cn = section_curvature(cap, uv, n_p, sgn)
        if chordwise:
            ctr = P.mean(0)
            ang = np.arctan2(P[:, 1] - ctr[1], P[:, 2] - ctr[2])
            o = np.argsort(ang)
            P, n, c, cn = P[o], n[o], c[o], cn[o]
            lips = [int(np.argmax(P[:, 2])), int(np.argmin(P[:, 2]))]
            res.append(analyse_trace(name, P, n, c, cn, lip_pts=lips))
        else:
            o = np.lexsort((uv[:, 0], uv[:, 1]))
            res.append(analyse_trace(name, P[o], n[o], c[o], cn[o], lip_dir=lip_dir))
    tnodes = np.linspace(0, 1, 241)
    step = 8 if not QUICK else 24
    for k in range(step // 2, 240, step):
        for lab, tk in (("guide", tnodes[k]), ("between", 0.5 * (tnodes[k] + tnodes[k + 1]))):
            p = float(p_of_t(tk))
            xa = float(g.xi_attach(p))
            Xa, Za = g.plan_xz(540.0, xa)
            O = par.tip_outline(np.array([p]))[0]
            dvec = np.array([O[0] - Xa, 0.0, O[1] - Za])
            n_p = np.cross(dvec, Y)
            run(f"{lab}-plane p={p:.4f}", n_p, np.array([Xa, 0.0, Za]), lip_dir=dvec / np.linalg.norm(dvec))
    for xi in (0.0005, 0.002, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.93, 0.95, 0.97, 0.99, 0.997):
        Xa, Za = g.plan_xz(540.0, xi)
        run(f"spanwise Z=const xi_att={xi}", np.array([0.0, 0.0, 1.0]), np.array([Xa, 0.0, Za]),
            lip_dir=np.array([-1.0, 0.0, 0.0]))
        _, zle1, _, _, c1 = g.frame(540.0)
        dvec = np.array([-1.0, 0.0, -(zle1 - xi * c1)])
        run(f"sweep-following xi={xi}", np.cross(dvec, Y), np.array([Xa, 0.0, Za]),
            lip_dir=dvec / np.linalg.norm(dvec))
    for S in (540.5, 542, 545, 550, 560, 570, 575, 580, 585, 590, 595, 598, 599.5, 599.95):
        run(f"chordwise X=-{S}", np.array([1.0, 0.0, 0.0]), np.array([-S, 0.0, 0.0]), chordwise=True)
    return res


# ----------------------------------------------------------------------------
def main():
    t0 = time.time()
    data = STEP.read_bytes()
    shape, faces, info = read_step(STEP)
    cap_f, skin_f, root_f = classify(faces)
    cap = cap_f["surf"]
    out = dict(step_file=str(STEP.relative_to(g.ROOT)), step_sha256=hashlib.sha256(data).hexdigest(),
               step_bytes=len(data), reimport=info)
    out["parent_preserved"] = check_parent_preserved(skin_f, root_f)
    nu, nv = (241, 401) if QUICK else (721, 1201)
    U, V, D, C = tabulate(cap, nu, nv)
    out["locked"] = check_locked(cap, U, V, C["P"])
    join = check_join(cap, skin_f, 241 if QUICK else 801)
    LE = np.array(cap.Value(0.0, 1.0).Coord())
    TE = np.array(cap.Value(1.0, 1.0).Coord())
    out["zebra_join"] = [zebra_join(join, exclude_mm=e, LE=LE, TE=TE) for e in (0.0, 0.1, 0.5)]
    out["join"] = {k: {kk: vv for kk, vv in J.items() if not kk.startswith("_")} for k, J in join.items()}
    out["principal_curvature_map"] = principal_map(cap, U, V, C)
    sec = sections(cap, U, V, C["P"], C["sign"])
    out["sections"] = dict(count=len(sec), all_pass=all(s["pass"] for s in sec),
                           worst=sorted(sec, key=lambda s: s["worst_ratio"])[:5], traces=sec)
    out["seconds"] = time.time() - t0
    (g.BUILD_DIR / "readback.json").write_text(json.dumps(out, indent=1))
    summary = dict(sha=out["step_sha256"], reimport=info, parent=out["parent_preserved"], locked=out["locked"],
                   join={k: {kk: v for kk, v in J.items() if kk.startswith("max")} for k, J in out["join"].items()},
                   zebra=out["zebra_join"],
                   pmap={k: out["principal_curvature_map"][k] for k in
                         ("smallest_radius_anywhere_mm", "strict_zone_smallest_radius_refined", "by_attachment_thickness", "failing_zone",
                          "exclusion_study", "grid_closest_to_LE_point_mm", "grid_closest_to_TE_point_mm")},
                   sections_all_pass=out["sections"]["all_pass"],
                   sections_worst=[(s["name"], s["worst_ratio"], s["smallest_radius_mm"]) for s in out["sections"]["worst"]],
                   seconds=out["seconds"])
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
