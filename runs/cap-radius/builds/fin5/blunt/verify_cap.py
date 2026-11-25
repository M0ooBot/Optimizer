"""Verify the reimported wing.step (fin5/blunt rounded tip cap).

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/fin5/blunt/verify_cap.py

Writes readback.json and figures (PNG) next to this file. All geometry is read back
from wing.step; surface derivatives come from the STEP B-spline definitions
(NumPy basis derivatives, cross-checked against Open CASCADE's D2 evaluator).
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cap_analytic as ca  # noqa: E402
import surfeval as se  # noqa: E402
from build_cap import read_step, faces_of  # noqa: E402

from OCP.BRep import BRep_Tool  # noqa: E402
from OCP.GeomAdaptor import GeomAdaptor_Surface  # noqa: E402
from OCP.GeomAbs import GeomAbs_BSplineSurface, GeomAbs_Plane  # noqa: E402
from OCP.TopExp import TopExp_Explorer  # noqa: E402
from OCP.TopAbs import TopAbs_SOLID, TopAbs_SHELL, TopAbs_FACE, TopAbs_EDGE  # noqa: E402
from OCP.BRepCheck import BRepCheck_Analyzer  # noqa: E402
from OCP.GProp import GProp_GProps  # noqa: E402
from OCP.BRepGProp import BRepGProp  # noqa: E402
from OCP.Bnd import Bnd_Box  # noqa: E402
from OCP.BRepBndLib import BRepBndLib  # noqa: E402
from OCP.ShapeAnalysis import ShapeAnalysis_FreeBounds  # noqa: E402
from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape  # noqa: E402
from OCP.TopExp import TopExp  # noqa: E402

STEP = HERE / "wing.step"
OUT = HERE / "readback.json"
LE_POINT = np.array([65.0, 0.0, 240.0])
TE_CORNERS = [np.array([150.0, -0.2, 240.0]), np.array([150.0, 0.2, 240.0])]


def count(shape, kind):
    n, e = 0, TopExp_Explorer(shape, kind)
    while e.More():
        n += 1
        e.Next()
    return n


def identify(shape):
    out = {}
    for f in faces_of(shape):
        a = GeomAdaptor_Surface(BRep_Tool.Surface_s(f))
        if a.GetType() == GeomAbs_Plane:
            out["root"] = f
            continue
        b = a.BSpline()
        if b.NbVPoles() == 2:
            out["te"] = f
        elif b.Value(0.5, b.VKnot(1)).Z() >= ca.Z0 - 1e-9 and b.NbUPoles() > 40:
            out["cap"] = f
        else:
            out["skin+" if b.Value(0.5, 120.0).Y() > 0 else "skin-"] = f
    assert set(out) == {"root", "te", "cap", "skin+", "skin-"}, out.keys()
    return out


def required_radius(X):
    return ca.required_radius(X)


def edge_frame(d, sign, edge_dir_key):
    """Normal and second fundamental form in a frame whose e1 follows the given param direction."""
    if edge_dir_key == "Sv":
        d = {"S": d["S"], "Su": d["Sv"], "Sv": d["Su"], "Suu": d["Svv"], "Suv": d["Suv"], "Svv": d["Suu"]}
        sign = -sign  # swapping u and v flips Su x Sv
    return se.frame(d, sign)


def join_metrics(cap, par, cap_uv, par_uv, cap_dir, par_dir, label):
    dc = cap.points(*cap_uv); dp = par.points(*par_uv)
    nc, e1c, e2c, Bc = edge_frame(dc, cap.sign, cap_dir)
    npar, e1p, e2p, Bp = edge_frame(dp, par.sign, par_dir)
    # align e1 directions
    flip = np.sign(np.einsum("ij,ij->i", e1c, e1p))
    sep = np.linalg.norm(dc["S"] - dp["S"], axis=-1)
    ang = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", nc, npar), -1, 1)))
    # curvature in the parent's frame: k_n along e1 (edge) and e2 (across), twist
    def comps(B, e1, e2, n, E1, E2):
        # express B (in frame e1,e2) in the common frame (E1, E2) projected
        R = np.stack([np.stack([np.einsum("ij,ij->i", E1, e1), np.einsum("ij,ij->i", E1, e2)], -1),
                      np.stack([np.einsum("ij,ij->i", E2, e1), np.einsum("ij,ij->i", E2, e2)], -1)], -2)
        Bc_ = R @ B @ np.swapaxes(R, -1, -2)
        return Bc_[:, 0, 0], Bc_[:, 1, 1], Bc_[:, 0, 1]
    E1, E2 = e1p, e2p
    kc = comps(Bc, e1c, e2c, nc, E1, E2)
    kp = comps(Bp, e1p, e2p, npar, E1, E2)
    across_diff = np.abs(kc[1] - kp[1]); along_diff = np.abs(kc[0] - kp[0]); twist_diff = np.abs(kc[2] - kp[2])
    i = int(np.nanargmax(ang)); j = int(np.nanargmax(across_diff))
    return {
        "label": label,
        "samples": int(sep.size),
        "max_separation_mm": float(sep.max()),
        "max_normal_angle_deg": float(np.nanmax(ang)),
        "nan_samples": int(np.isnan(ang).sum()),
        "max_normal_angle_at_XYZ": dp["S"][i].round(5).tolist(),
        "max_across_normal_curvature_diff_per_mm": float(np.nanmax(across_diff)),
        "max_across_diff_at_XYZ": dp["S"][j].round(5).tolist(),
        "across_k_parent_there": float(kp[1][j]), "across_k_cap_there": float(kc[1][j]),
        "max_along_normal_curvature_diff_per_mm": float(np.nanmax(along_diff)),
        "max_twist_diff_per_mm": float(np.nanmax(twist_diff)),
        "_arrays": {"S": dp["S"], "ang": ang, "across": across_diff, "kc": kc[1], "kp": kp[1]},
    }


def main():
    res = {}
    data = STEP.read_bytes()
    res["step_file"] = STEP.name
    res["step_sha256"] = hashlib.sha256(data).hexdigest()
    res["step_bytes"] = len(data)
    shape = read_step(STEP)
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props, 1e-9)
    box = Bnd_Box(); box.SetGap(0.0); BRepBndLib.AddOptimal_s(shape, box, False, False)
    fb = ShapeAnalysis_FreeBounds(shape)
    free_edges = count(fb.GetOpenWires(), TopAbs_EDGE)
    res["reimport"] = {
        "solids": count(shape, TopAbs_SOLID), "shells": count(shape, TopAbs_SHELL), "faces": count(shape, TopAbs_FACE),
        "brepcheck_valid": bool(BRepCheck_Analyzer(shape).IsValid()),
        "free_edges": int(free_edges),
        "volume_mm3": props.Mass(),
        "bounding_box_mm": [round(x, 7) for x in box.Get()],
    }
    F = identify(shape)
    faces = {k: se.Face(F[k]) for k in ("cap", "skin+", "skin-", "te")}
    # orient normals outward from a known direction at the middle of each face
    outward = {"cap": [0, 0, 1], "skin+": [0, 1, 0], "skin-": [0, -1, 0], "te": [1, 0, 0]}
    for k, f in faces.items():
        um, vm = np.mean(f.urange), np.mean(f.vrange)
        d = f.points(np.array([um]), np.array([vm]))
        n = np.cross(d["Su"][0], d["Sv"][0])
        f.sign = 1.0 if np.dot(n, outward[k]) > 0 else -1.0
    res["reimport"]["normal_orientation_note"] = "normals oriented outward by checking each face at its parameter midpoint"
    res["reimport"]["cap_face"] = {"degree": [faces["cap"].pu, faces["cap"].pv], "poles": list(faces["cap"].P.shape[:2])}
    res["derivative_crosscheck_max_rel_err_vs_occ"] = {k: f.occ_check() for k, f in faces.items()}

    # ---------------- parent preserved
    pshape = read_step(ca.PARENT_DIR / "parent.step")
    PF = identify_parent(pshape)
    pres = {}
    for k in ("skin+", "skin-", "te"):
        a = se.Face(PF[k]); b = faces[k]
        pres[k] = {"same_pole_grid": a.P.shape == b.P.shape,
                   "max_pole_diff_mm": float(np.abs(a.P - b.P).max()) if a.P.shape == b.P.shape else None,
                   "same_knots": bool(np.allclose(a.ku, b.ku) and np.allclose(a.kv, b.kv))}
    # analytic parent vs reimported skins
    s = np.linspace(0, 1, 301); z = np.linspace(0, 240, 241)
    dmax = 0.0
    for side, key in ((+1, "skin+"), (-1, "skin-")):
        g = faces[key].grid(s, z)["S"]
        ref = ca.parent.skin(s[:, None], z[None, :], side)
        dmax = max(dmax, float(np.linalg.norm(g - ref, axis=-1).max()))
    pres["parent_py_to_reimported_skins_max_mm"] = dmax
    pv = GProp_GProps(); BRepGProp.VolumeProperties_s(pshape, pv, 1e-9)
    pres["parent_volume_mm3"] = pv.Mass()
    res["parent_preserved"] = pres

    cap = faces["cap"]
    # ---------------- locked outline / dimensions
    u = np.unique(np.concatenate([np.linspace(0, 1, 4001), np.geomspace(1e-7, 1e-2, 400), 1 - np.geomspace(1e-7, 1e-2, 400)]))
    lip = cap.grid(u, np.array([0.5]))["S"][:, 0]
    xi, h = ca.xi_h(lip[:, 0], lip[:, 2])
    resid = xi**2 + np.clip(h, 0, None) ** 3 - 1
    # distance to outline: Newton-like via gradient of F(X,Z)=xi^2+h^3-1
    a_ = 80 - 37.5 * lip[:, 2] / 240
    dFdX = 2 * xi / a_
    dFdZ = 3 * h**2 / 5 + 2 * xi * (-(27.5 / 240) * a_ - (lip[:, 0] - (80 + 27.5 * lip[:, 2] / 240)) * (-37.5 / 240)) / a_**2
    dist = np.abs(resid) / np.hypot(dFdX, dFdZ)
    # every cap point inside the outline (side view)
    tt = np.linspace(0, 1, 801)
    G = cap.grid(u[::4], tt)["S"]
    xi2, h2 = ca.xi_h(G[..., 0], G[..., 2])
    outside = xi2**2 + np.clip(h2, 0, None) ** 3 - 1
    dist_out = outside / np.hypot(2 * xi2 / (80 - 37.5 * G[..., 2] / 240), 3 * np.clip(h2, 0, None) ** 2 / 5 + 1e-12)
    crown_i = int(np.argmax(lip[:, 2]))
    res["locked"] = {
        "max_Z_mm": float(G[..., 2].max()),
        "crown_XYZ": lip[crown_i].tolist(),
        "crown_reference_XYZ": [108.07291666666667, 0.0, 245.0],
        "lip_max_distance_to_outline_mm": float(dist.max()),
        "lip_max_abs_Y_mm": float(np.abs(lip[:, 1]).max()),
        "cap_max_excursion_outside_outline_mm": float(max(0.0, dist_out.max())),
        "attachment_Z_mm": ca.Z0,
        "cap_min_Z_mm": float(G[..., 2].min()),
        "cap_min_Z_at_XYZ": G.reshape(-1, 3)[np.argmin(G[..., 2])].round(6).tolist(),
        "cap_max_outside_at_XYZ": G.reshape(-1, 3)[np.argmax(dist_out)].round(6).tolist(),
    }

    # ---------------- joins
    uu = np.unique(np.concatenate([np.linspace(0, 1, 2001), np.geomspace(1e-6, 1e-2, 300), 1 - np.geomspace(1e-6, 1e-2, 300)]))
    J = {}
    J["+Y_skin"] = join_metrics(cap, faces["skin+"], (uu, np.ones_like(uu)), (uu, np.full_like(uu, 240.0)), "Su", "Su", "+Y skin / cap at Z=240")
    J["-Y_skin"] = join_metrics(cap, faces["skin-"], (uu, np.zeros_like(uu)), (uu, np.full_like(uu, 240.0)), "Su", "Su", "-Y skin / cap at Z=240")
    tv = np.unique(np.concatenate([np.linspace(0, 1, 801), np.geomspace(1e-6, 1e-2, 100), 1 - np.geomspace(1e-6, 1e-2, 100)]))
    J["TE_face"] = join_metrics(cap, faces["te"], (np.ones_like(tv), tv), (np.full_like(tv, 240.0), tv), "Sv", "Sv", "TE face / cap along the TE segment")
    for k, j in J.items():
        arr = j.pop("_arrays")
        S = arr["S"]
        dc = np.min([np.linalg.norm(S - c, axis=-1) for c in TE_CORNERS], axis=0)
        dle = np.linalg.norm(S - LE_POINT, axis=-1)
        band = {}
        for lim in (0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0):
            m = (dle <= lim)
            if m.any() and k != "TE_face":
                band[f"within_{lim}mm_of_LE_point"] = {
                    "max_normal_angle_deg": float(np.nanmax(arr["ang"][m])),
                    "max_across_curvature_diff_per_mm": float(np.nanmax(arr["across"][m])),
                    "nan_samples(collapsed LE point)": int(np.isnan(arr["ang"][m]).sum())}
            m = (dc <= lim)
            if m.any():
                band[f"within_{lim}mm_of_TE_corner"] = {
                    "max_normal_angle_deg": float(np.nanmax(arr["ang"][m])),
                    "max_across_curvature_diff_per_mm": float(np.nanmax(arr["across"][m]))}
        if k == "TE_face":
            Yv = S[:, 1]
            for lim in (0.19, 0.18, 0.15, 0.1, 0.05):
                m = np.abs(Yv) <= lim
                band[f"|Y|<={lim}mm (corner at 0.2)"] = {"max_normal_angle_deg": float(arr["ang"][m].max()),
                                                    "max_across_curvature_diff_per_mm": float(arr["across"][m].max())}
        for lim in (0.01, 0.05, 0.1, 0.2, 0.5, 1.0):
            m = (dc > lim) & (dle > lim)
            band[f"beyond_{lim}mm_of_TE_corners_and_LE"] = {
                "max_normal_angle_deg": float(np.nanmax(arr["ang"][m])) if m.any() else None,
                "max_across_curvature_diff_per_mm": float(np.nanmax(arr["across"][m])) if m.any() else None}
        j["by_distance"] = band
    res["joins"] = J

    # ---------------- curvature map over the whole cap
    res["curvature_map"] = curvature_map(cap)
    globals()["_FACES"] = faces
    globals()["_RES"] = res
    return res, faces


def identify_parent(shape):
    out = {}
    for f in faces_of(shape):
        a = GeomAdaptor_Surface(BRep_Tool.Surface_s(f))
        if a.GetType() == GeomAbs_Plane:
            continue
        b = a.BSpline()
        if b.NbVPoles() == 2:
            out["te"] = f
        else:
            out["skin+" if b.Value(0.5, 120.0).Y() > 0 else "skin-"] = f
    return out


def rmap(cap, u, t):
    d = cap.grid(u, t)
    n, e1, e2, B = se.frame(d, cap.sign)
    k = se.principal(B)
    kmax = np.abs(k).max(-1)
    return d["S"], 1.0 / kmax, k


def curvature_map(cap):
    from scipy.optimize import minimize
    u = np.unique(np.concatenate([np.linspace(0, 1, 1201), np.geomspace(1e-6, 2e-2, 200), 1 - np.geomspace(1e-7, 2e-2, 300)]))
    t = np.unique(np.concatenate([np.linspace(0, 1, 1001), np.geomspace(1e-7, 2e-2, 150), 1 - np.geomspace(1e-7, 2e-2, 150)]))
    u = u[u >= 1e-6]
    S, r, k = rmap(cap, u, t)
    req = required_radius(S[..., 0])
    ratio = r / req
    dc = np.min([np.linalg.norm(S - c, axis=-1) for c in TE_CORNERS], axis=0)
    dle = np.linalg.norm(S - LE_POINT, axis=-1)
    out = {"grid": [int(u.size), int(t.size)], "u_min": float(u.min()),
           "grid_closest_to_LE_point_mm": float(dle.min()),
           "grid_closest_to_TE_corner_mm": float(dc.min())}

    def point_info(i):
        X = S[i][0]
        return {"r_mm": float(r[i]), "required_mm": float(req[i]), "ratio": float(ratio[i]),
                "XYZ": S[i].round(5).tolist(), "attachment_thickness_mm": float(ca.attachment_thickness(X)),
                "chord_fraction": float((X - 65) / 85), "dist_TE_corner_mm": float(dc[i]), "dist_LE_mm": float(dle[i])}
    i = np.unravel_index(np.argmin(r), r.shape)
    out["min_radius_anywhere_on_grid"] = point_info(i)
    for lim in (0.05, 0.1, 0.2, 0.25, 0.3, 0.5):
        m = dc > lim
        j = np.unravel_index(np.argmin(np.where(m, ratio, np.inf)), r.shape)
        jr = np.unravel_index(np.argmin(np.where(m, r, np.inf)), r.shape)
        out[f"beyond_{lim}mm_of_TE_corners"] = {"min_ratio": point_info(j), "min_radius": point_info(jr)}
    bad = ratio < 1
    out["points_below_rule"] = int(bad.sum())
    out["points_below_rule_max_dist_to_TE_corner_mm"] = float(dc[bad].max()) if bad.any() else 0.0
    out["points_below_rule_max_dist_to_LE_mm_if_not_near_TE"] = float(dle[bad & (dc > 1)].max()) if (bad & (dc > 1)).any() else None
    # refine local minima of ratio away from the corners (d > 0.25)
    from scipy.ndimage import minimum_filter
    loc = (ratio == minimum_filter(ratio, size=7)) & (dc > 0.25) & (dle > 1e-3)
    idx = np.argwhere(loc)
    order = np.argsort(ratio[loc])[:12]
    refined = []
    for (a, b) in idx[order]:
        def f(x):
            d = cap.points(np.array([x[0]]), np.array([x[1]]))
            n, e1, e2, B = se.frame(d, cap.sign)
            kk = np.abs(se.principal(B)).max()
            return float((1 / kk) / required_radius(d["S"][0, 0]))
        x0 = np.array([u[a], t[b]])
        sol = minimize(f, x0, method="L-BFGS-B", bounds=[(max(1e-6, u[max(a - 3, 0)]), u[min(a + 3, u.size - 1)]), (t[max(b - 3, 0)], t[min(b + 3, t.size - 1)])])
        d = cap.points(np.array([sol.x[0]]), np.array([sol.x[1]]))
        n, e1, e2, B = se.frame(d, cap.sign)
        kk = np.abs(se.principal(B)).max()
        X = d["S"][0]
        refined.append({"u": float(sol.x[0]), "t": float(sol.x[1]), "XYZ": X.round(5).tolist(), "r_mm": float(1 / kk),
                        "required_mm": float(required_radius(X[0])), "ratio": float(sol.fun),
                        "attachment_thickness_mm": float(ca.attachment_thickness(X[0])),
                        "dist_TE_corner_mm": float(min(np.linalg.norm(X - c) for c in TE_CORNERS))})
    refined.sort(key=lambda d: d["ratio"])
    out["refined_local_minima_of_ratio_beyond_0p25mm_of_TE_corners"] = refined[:8]
    # radius vs distance from the TE corners
    prof = []
    for lo, hi in ((0, 0.02), (0.02, 0.05), (0.05, 0.1), (0.1, 0.15), (0.15, 0.2), (0.2, 0.25), (0.25, 0.3), (0.3, 0.5), (0.5, 1.0)):
        m = (dc >= lo) & (dc < hi)
        if m.any():
            prof.append({"dist_mm": [lo, hi], "min_r_mm": float(r[m].min()), "min_ratio": float(ratio[m].min())})
    out["near_TE_corners_radius_by_distance"] = prof
    globals()["_MAP"] = (u, t, S, r, req, dc, dle)
    return out


# =========================================================================== sections
def _refine(face, uv, n_pl, p0, iters=6):
    u, v = uv[:, 0].copy(), uv[:, 1].copy()
    for _ in range(iters):
        d = face.points(u, v)
        f = (d["S"] - p0) @ n_pl
        fu = d["Su"] @ n_pl; fv = d["Sv"] @ n_pl
        g2 = fu**2 + fv**2 + 1e-300
        u = np.clip(u - f * fu / g2, *face.urange); v = np.clip(v - f * fv / g2, *face.vrange)
    return u, v


def plane_trace(face, n_pl, p0, ugrid, vgrid):
    """Ordered points of face ∩ plane: list of (u, v) arrays (one per connected piece)."""
    import contourpy
    S = face.grid(ugrid, vgrid)["S"]
    f = (S - p0) @ n_pl
    gen = contourpy.contour_generator(x=vgrid, y=ugrid, z=f)
    out = []
    for line in gen.lines(0.0):
        uv = np.stack([line[:, 1], line[:, 0]], -1)
        if len(uv) < 5:
            continue
        u, v = _refine(face, uv, n_pl, p0)
        out.append(np.stack([u, v], -1))
    return out


def trace_curvature(face, uv, n_pl):
    d = face.points(uv[:, 0], uv[:, 1])
    n, e1, e2, B = se.frame(d, face.sign)
    T = np.cross(n_pl, n); T /= np.linalg.norm(T, axis=-1, keepdims=True)
    Nc = np.cross(n_pl, T)
    Tf = np.stack([np.einsum("ij,ij->i", T, e1), np.einsum("ij,ij->i", T, e2)], -1)
    kn = np.einsum("ni,nij,nj->n", Tf, B, Tf)
    cosang = np.abs(np.einsum("ij,ij->i", n, Nc))
    ks = -kn / np.maximum(cosang, 1e-12)            # convex positive
    kp = se.principal(B)
    rp = 1.0 / np.abs(kp).max(-1)
    return d["S"], ks, rp


def _arclen(P):
    return np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=-1))])


def analyse_trace(name, P, ks, rp, kind, plane_desc, parent_join=None):
    s = _arclen(P)
    r = 1.0 / np.maximum(np.abs(ks), 1e-15)
    thr = 1e-3
    sig = np.where(np.abs(ks) > thr, np.sign(ks), 0)
    nz = np.nonzero(sig)[0]
    rev = [int(nz[k + 1]) for k in range(len(nz) - 1) if sig[nz[k]] != sig[nz[k + 1]]]
    i = int(np.argmin(r))
    X = P[i, 0]
    req = float(required_radius(X))
    out = {"name": name, "kind": kind, "plane": plane_desc, "points": int(len(P)), "length_mm": float(s[-1]),
           "min_section_radius_mm": float(r[i]), "min_section_radius_at_XYZ": P[i].round(4).tolist(),
           "surface_min_principal_radius_there_mm": float(rp[i]),
           "attachment_thickness_there_mm": float(ca.attachment_thickness(X)), "required_radius_mm": req,
           "surface_radius_meets_rule": bool(rp[i] >= req * (1 - 1e-9)),
           "min_surface_principal_radius_on_trace_mm": float(rp.min()),
           "min_surface_principal_radius_at_XYZ": P[int(np.argmin(rp))].round(4).tolist(),
           "bending_reversals(|k|>1e-3/mm)": [{"at_XYZ": P[j].round(4).tolist(), "arc_mm": float(s[j])} for j in rev],
           "most_concave_curvature_per_mm": float(ks.min()), "most_concave_at_XYZ": P[int(np.argmin(ks))].round(5).tolist()}
    lips = np.nonzero(np.diff(np.sign(P[:, 1])) != 0)[0]
    if len(lips):
        L = []
        for j in lips:
            L.append({"at_XYZ": P[j].round(4).tolist(), "lip_section_radius_mm": float(r[j]),
                      "lip_surface_min_principal_radius_mm": float(rp[j]),
                      "arc_from_min_radius_to_lip_mm": float(abs(s[j] - s[i]))})
        out["lip"] = L
        jl = lips[0]
        out["min_over_lip_section_radius_ratio"] = float(r[i] / r[jl])
    return out


def sections(faces, res):
    cap = faces["cap"]
    ug = np.unique(np.concatenate([np.linspace(0, 1, 1601), np.geomspace(1e-6, 0.02, 200), 1 - np.geomspace(1e-7, 0.02, 300)]))
    vg = np.unique(np.concatenate([np.linspace(0, 1, 1201), np.geomspace(1e-6, 0.02, 100), 1 - np.geomspace(1e-6, 0.02, 100)]))
    planes = []
    for s in (0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995):
        A = ca.parent.skin(np.array([s]), ca.Z0, +1)[0]; Pl = ca.lip_point(np.array([s]))[0]
        dirv = Pl - A; dirv[1] = 0; dirv /= np.linalg.norm(dirv)
        n = np.cross(dirv, [0, 1, 0]); n /= np.linalg.norm(n)
        planes.append((f"guide s={s}", "guide plane (contains Y and the side-view line A(s)-P(s))", n, np.array([A[0], 0, A[2]])))
    for X0 in (65.5, 66.0, 68.0, 70.0, 80.0, 90.0, 100.0, 108.07, 120.0, 130.0, 140.0, 145.0, 148.0, 149.0, 149.5, 149.8, 149.9, 149.95):
        planes.append((f"X={X0}", "constant-X plane", np.array([1.0, 0, 0]), np.array([X0, 0, 0])))
    for k in (0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98):
        dirv = np.array([65 - 75 * k, 0, 240.0]); dirv /= np.linalg.norm(dirv)
        n = np.cross(dirv, [0, 1, 0]); n /= np.linalg.norm(n)
        planes.append((f"sweep k={k}", "sweep-following plane (contains Y and the parent's iso-chord line at chord fraction k)", n, np.array([65 + 85 * k, 0, 240.0])))
    traces = []
    for name, desc, n, p0 in planes:
        pcs = plane_trace(cap, n, p0, ug, vg)
        if not pcs:
            continue
        uv = max(pcs, key=len)
        P, ks, rp = trace_curvature(cap, uv, n)
        if P[0, 1] > P[-1, 1]:
            P, ks, rp = P[::-1], ks[::-1], rp[::-1]
        a = analyse_trace(name, P, ks, rp, "spanwise wrap", desc)
        # parent side of the join: curvature of the same plane section on each skin at the attachment
        jinfo = {}
        for key in ("skin-", "skin+"):
            sk = faces[key]
            pp = plane_trace(sk, n, p0, np.linspace(0, 1, 2001), np.linspace(236.0, 240.0, 81))
            if pp:
                uvp = max(pp, key=len)
                Pp, ksp, rpp = trace_curvature(sk, uvp, n)
                j = int(np.argmax(Pp[:, 2]))
                end = P[0] if key == "skin-" else P[-1]
                kc = ks[0] if key == "skin-" else ks[-1]
                jinfo[key] = {"parent_section_curvature_at_join_per_mm": float(ksp[j]), "cap_section_curvature_at_join_per_mm": float(kc),
                              "join_point_gap_mm": float(np.linalg.norm(Pp[j] - end))}
        a["join_curvature"] = jinfo
        traces.append(a)
    # horizontal sections (thin ends of each section, near the LE and TE)
    hz = []
    for Z0 in (240.02, 240.05, 240.1, 240.2, 240.3, 240.5, 241.0, 242.0, 243.0, 244.0, 244.8):
        n = np.array([0, 0, 1.0]); p0 = np.array([0, 0, Z0])
        pcs = plane_trace(cap, n, p0, ug, vg)
        if not pcs:
            continue
        uv = max(pcs, key=len)
        P, ks, rp = trace_curvature(cap, uv, n)
        info = {"Z": Z0, "points": int(len(P))}
        for end, sel in (("LE end", P[:, 0] < P[:, 0].min() + 3.0), ("TE end", P[:, 0] > P[:, 0].max() - 3.0)):
            if not sel.any():
                continue
            r = 1 / np.maximum(np.abs(ks[sel]), 1e-15)
            Ps = P[sel]; i = int(np.argmin(r)); X = Ps[i, 0]
            tip = int(np.argmin(Ps[:, 0])) if end == "LE end" else int(np.argmax(Ps[:, 0]))
            req = float(required_radius(X))
            info[end] = {"min_section_radius_mm": float(r[i]), "at_XYZ": Ps[i].round(5).tolist(),
                         "surface_min_principal_radius_there_mm": float(rp[sel][i]),
                         "radius_at_section_tip_mm": float(r[tip]), "tip_XYZ": Ps[tip].round(5).tolist(),
                         "attachment_thickness_there_mm": float(ca.attachment_thickness(X)), "required_mm": req,
                         "dist_to_TE_corner_mm": float(min(np.linalg.norm(Ps[i] - c) for c in TE_CORNERS)),
                         "surface_meets_rule": bool(rp[sel][i] >= req * (1 - 1e-9))}
        hz.append(info)
    return {"wrap_traces": traces, "horizontal_sections": hz}


# =========================================================================== zebra
def _zebra_phase_grad(d, sign, eye_dir, axis, period):
    n, e1, e2, B = se.frame(d, sign)
    dv = -np.asarray(eye_dir, float)
    dn = n @ dv
    r = dv[None, :] - 2 * dn[:, None] * n
    ra = r @ axis
    ok = np.abs(ra) < 0.999
    phase = np.arcsin(np.clip(ra, -1, 1)) / period
    an = n @ axis
    d_t = dv[None, :] - dn[:, None] * n
    a_t = axis[None, :] - an[:, None] * n
    w = an[:, None] * d_t + dn[:, None] * a_t
    wf = np.stack([np.einsum("ij,ij->i", w, e1), np.einsum("ij,ij->i", w, e2)], -1)
    gf = 2 * np.einsum("nij,nj->ni", B, wf) / (period * np.sqrt(np.maximum(1 - ra**2, 1e-12)))[:, None]
    g = gf[:, :1] * e1 + gf[:, 1:] * e2
    return phase, g, ok, n


def zebra(faces, res):
    cap = faces["cap"]
    period = np.radians(10.0)
    u = np.unique(np.concatenate([np.linspace(0, 1, 4001), np.geomspace(1e-5, 0.01, 200), 1 - np.geomspace(1e-6, 0.01, 300)]))
    u = u[(u > 1e-5)]
    out = {"stripe_period_deg": 10.0, "phase_definition": "phase = asin(r.a)/period, r = reflected view ray, a = stripe axis (in the image plane); light/dark changes at half-integer... integer multiples of 0.5",
           "views_az_el_deg": [], "stripe_axes_deg": [0, 30, 60, 90, 120, 150]}
    summary = {}
    for side, key, tv in ((+1, "skin+", 1.0), (-1, "skin-", 0.0)):
        sk = faces[key]
        dc_ = cap.points(u, np.full_like(u, tv)); dp_ = sk.points(u, np.full_like(u, 240.0))
        S = dp_["S"]
        dco = np.min([np.linalg.norm(S - c, axis=-1) for c in TE_CORNERS], axis=0)
        dle = np.linalg.norm(S - LE_POINT, axis=-1)
        rows = []
        for az in (-60, -30, 0, 30, 60):
            for el in (-30, 0, 30, 60):
                e = np.array([np.sin(np.radians(az)) * np.cos(np.radians(el)), side * np.cos(np.radians(az)) * np.cos(np.radians(el)), np.sin(np.radians(el))])
                up = np.array([0, 0, 1.0]) if abs(el) < 80 else np.array([1.0, 0, 0])
                u1 = np.cross(up, e); u1 /= np.linalg.norm(u1); u2 = np.cross(e, u1)
                for al in out["stripe_axes_deg"]:
                    a = np.cos(np.radians(al)) * u1 + np.sin(np.radians(al)) * u2
                    pc, gc, okc, nc = _zebra_phase_grad(dc_, cap.sign, e, a, period)
                    pp, gp, okp, npar = _zebra_phase_grad(dp_, sk.sign, e, a, period)
                    vis = (nc @ e > 0.05) & okc & okp
                    dphi = np.abs(pc - pp)
                    gmag_c = np.linalg.norm(gc, axis=-1); gmag_p = np.linalg.norm(gp, axis=-1)
                    cosg = np.einsum("ij,ij->i", gc, gp) / np.maximum(gmag_c * gmag_p, 1e-300)
                    dang = np.degrees(np.arccos(np.clip(cosg, -1, 1)))
                    spacing = np.abs(gmag_c - gmag_p) / np.maximum(gmag_p, 1e-12)
                    # light/dark crossings along the join (parent side)
                    k2 = np.floor(2 * pp)
                    cr = np.nonzero((np.diff(k2) != 0) & vis[:-1] & vis[1:])[0]
                    rows.append((az, el, al, vis, dphi, dang, spacing, cr, gmag_p))
        def band(lim):
            m0 = (dco > lim) & (dle > lim)
            best = {"max_phase_diff_periods": 0.0, "max_stripe_dir_change_deg_at_crossings": 0.0,
                    "max_band_spacing_change_rel_at_crossings": 0.0, "crossings": 0}
            loc = {}
            for az, el, al, vis, dphi, dang, spacing, cr, gm in rows:
                m = vis & m0
                if m.any():
                    i = int(np.argmax(np.where(m, dphi, -1)))
                    if dphi[i] > best["max_phase_diff_periods"]:
                        best["max_phase_diff_periods"] = float(dphi[i]); loc["phase"] = {"XYZ": S[i].round(4).tolist(), "view_az_el": [az, el], "axis_deg": al}
                crm = [c for c in cr if m0[c] and m0[c + 1]]
                best["crossings"] += len(crm)
                for c in crm:
                    if dang[c] > best["max_stripe_dir_change_deg_at_crossings"]:
                        best["max_stripe_dir_change_deg_at_crossings"] = float(dang[c]); loc["direction"] = {"XYZ": S[c].round(4).tolist(), "view_az_el": [az, el], "axis_deg": al}
                    if spacing[c] > best["max_band_spacing_change_rel_at_crossings"]:
                        best["max_band_spacing_change_rel_at_crossings"] = float(spacing[c]); loc["spacing"] = {"XYZ": S[c].round(4).tolist(), "view_az_el": [az, el], "axis_deg": al}
            best["locations"] = loc
            return best
        summary[f"{key}/cap join"] = {f"beyond_{lim}mm_of_TE_corners_and_LE": band(lim) for lim in (0.01, 0.05, 0.1, 0.2, 0.5, 1.0)}
    out["views_az_el_deg"] = [[az, el] for az in (-60, -30, 0, 30, 60) for el in (-30, 0, 30, 60)]
    out["joins"] = summary
    out["sampled_region"] = "both cut edges from the LE point to the TE corners (u in [1e-5, 1], dense near both ends); visible points only (n.e > 0.05) and |r.a| < 0.999"
    return out


# =========================================================================== figures
def figures(faces):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    def zebra_colors(face, U, V, eye, axis, period):
        d = face.grid(U, V)
        shp = d["S"].shape[:2]
        flat = {k: v.reshape(-1, 3) for k, v in d.items()}
        ph, g, ok, n = _zebra_phase_grad(flat, face.sign, eye, axis, period)
        stripe = (np.floor(2 * ph) % 2).reshape(shp)
        shade = 0.55 + 0.45 * np.clip((n @ eye).reshape(shp), 0, 1)
        col = np.where(stripe[..., None] > 0, 0.92, 0.12) * shade[..., None]
        return d["S"], np.concatenate([np.repeat(col, 3, -1), np.ones(shp + (1,))], -1)

    def render(fname, parts, eye, axis_deg, xlim, ylim, zlim, title, period_deg=6.0):
        e = eye / np.linalg.norm(eye)
        up = np.array([0, 0, 1.0])
        u1 = np.cross(up, e); u1 /= np.linalg.norm(u1); u2 = np.cross(e, u1)
        a = np.cos(np.radians(axis_deg)) * u1 + np.sin(np.radians(axis_deg)) * u2
        fig = plt.figure(figsize=(9, 7), dpi=130)
        ax = fig.add_subplot(projection="3d", proj_type="ortho")
        for face, U, V in parts:
            S, C = zebra_colors(face, U, V, e, a, np.radians(period_deg))
            ax.plot_surface(S[..., 0], S[..., 1], S[..., 2], facecolors=C, linewidth=0, antialiased=False, shade=False, rstride=1, cstride=1)
        ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_zlim(*zlim)
        ax.set_box_aspect((xlim[1] - xlim[0], ylim[1] - ylim[0], zlim[1] - zlim[0]))
        el = np.degrees(np.arcsin(e[2])); az = np.degrees(np.arctan2(e[1], e[0]))
        ax.view_init(elev=el, azim=az)
        ax.set_axis_off()
        ax.set_title(title, fontsize=9)
        fig.savefig(HERE / fname, bbox_inches="tight")
        plt.close(fig)

    cap, sp, sm = faces["cap"], faces["skin+"], faces["skin-"]
    uc = np.unique(np.concatenate([np.linspace(0, 1, 700), np.geomspace(1e-5, 0.02, 60), 1 - np.geomspace(1e-6, 0.02, 80)]))
    tc = np.linspace(0, 1, 500)
    names = []
    # radius / requirement map over the cap (side view of the +Y half), from the curvature map grid
    u, t, S, r, req, dc, dle = globals()["_MAP"]
    half = t >= 0.5
    ratio = (r / req)[:, half]; P = S[:, half]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), dpi=130, gridspec_kw={"width_ratios": [3, 1.2]})
    for ax, (xl, zl) in zip(axs, (((63, 152), (239.8, 245.3)), ((149.4, 150.05), (239.98, 240.6)))):
        m = (P[..., 0] > xl[0]) & (P[..., 0] < xl[1]) & (P[..., 2] < zl[1])
        sc = ax.scatter(P[..., 0][m], P[..., 2][m], c=np.log10(np.clip(ratio[m], 1e-3, 1e2)), s=1.2, cmap="viridis", vmin=-1, vmax=1.5, rasterized=True)
        bad = m & (ratio < 1)
        ax.scatter(P[..., 0][bad], P[..., 2][bad], c="red", s=1.5, rasterized=True, label="r < required")
        ax.set_xlim(*xl); ax.set_ylim(*zl); ax.set_xlabel("X (mm)"); ax.set_ylabel("Z (mm)")
    axs[0].set_title("log10(min principal radius / required radius), +Y half, side view", fontsize=9)
    axs[1].set_title("TE corner close-up", fontsize=9)
    axs[1].legend(fontsize=7, loc="upper left")
    fig.colorbar(sc, ax=axs, shrink=0.8)
    fig.savefig(HERE / "fig_radius_ratio_map.png", bbox_inches="tight"); plt.close(fig)
    names.append("fig_radius_ratio_map.png")
    # whole fin, +Y side
    zs = np.linspace(0, 240, 200); us = np.linspace(0, 1, 300)
    render("fig_zebra_whole_fin.png", [(sp, us, zs), (cap, uc, tc)], np.array([0.25, 1.0, 0.25]), 30,
           (0, 160), (-20, 20), (0, 250), "Whole fin with cap, +Y side, zebra (fixed camera)")
    names.append("fig_zebra_whole_fin.png")
    # tip close-up, +Y side, two stripe angles
    zs2 = np.linspace(228, 240, 100)
    for ang in (30, 100):
        f = f"fig_zebra_tip_axis{ang}.png"
        render(f, [(sp, np.linspace(0, 1, 700), zs2), (cap, uc, tc)], np.array([0.2, 1.0, 0.8]), ang,
               (60, 155), (-12, 12), (228, 247), f"Tip, +Y side from above, stripe axis {ang} deg (attachment at Z = 240)", period_deg=8.0)
        names.append(f)
    # LE and TE close-ups
    uL = np.linspace(0, 0.08, 300)
    render("fig_zebra_LE.png", [(sp, uL, np.linspace(237, 240, 80)), (sm, uL, np.linspace(237, 240, 80)), (cap, np.unique(np.concatenate([np.geomspace(1e-5, 0.08, 400)])), np.linspace(0, 1, 700))],
           np.array([-1.0, 0.5, 0.35]), 45, (63, 70), (-3.5, 3.5), (237, 243), "LE end of the cap, zebra", period_deg=12.0)
    names.append("fig_zebra_LE.png")
    uT = np.unique(np.concatenate([1 - np.geomspace(1e-7, 0.06, 500), np.linspace(0.94, 1, 200)]))
    render("fig_zebra_TE.png", [(sp, uT, np.linspace(238.5, 240, 80)), (sm, uT, np.linspace(238.5, 240, 80)), (cap, uT, np.linspace(0, 1, 800)), (faces["te"], np.linspace(238.5, 240, 40), np.linspace(0, 1, 40))],
           np.array([1.0, 0.6, 0.35]), 45, (147.8, 150.3), (-1.25, 1.25), (238.5, 241.3), "TE end of the cap with the blunt TE face, zebra", period_deg=15.0)
    names.append("fig_zebra_TE.png")
    return names


def cap_volume_check(res):
    """Volume added by the cap: STEP (whole minus parent) against the analytic design integral."""
    def area(Z):
        h = (Z - 240) / 5; a = 80 - 37.5 * Z / 240; m = 80 + 27.5 * Z / 240
        xo = np.sqrt(max(1 - h**3, 0)); X = m + a * xo * (-np.cos(np.linspace(0, np.pi, 20001)))
        Y = np.sqrt(np.clip(ca.G_value(X, np.full_like(X, Z)), 0, None))
        return np.trapezoid(2 * Y, X)
    zs = 240 + 5 * (1 - np.cos(np.linspace(0, np.pi / 2, 801)))
    analytic = float(np.trapezoid([area(z) for z in zs], zs))
    step = res["reimport"]["volume_mm3"] - res["parent_preserved"]["parent_volume_mm3"]
    return {"cap_volume_step_mm3": step, "cap_volume_analytic_design_mm3": analytic, "difference_mm3": step - analytic}


def run_all():
    res, faces = main()
    res["cap_volume"] = cap_volume_check(res)
    res["sections"] = sections(faces, res)
    res["zebra"] = zebra(faces, res)
    res["figures"] = figures(faces)
    OUT.write_text(json.dumps(res, indent=1, default=float))
    return res


if __name__ == "__main__":
    if "--figures-only" in sys.argv:
        _res, _faces = main()
        print(figures(_faces))
        sys.exit(0)
    r = run_all()
    print(json.dumps({k: r[k] for k in ("step_sha256", "reimport", "locked")}, indent=1))
