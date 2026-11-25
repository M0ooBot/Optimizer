"""Principal-radius map of a cap face (OCCT surface derivatives) judged against the printable rule.

Used by verify.py on the reimported STEP; can also be run on the fitted patch (cap_patch.npz)
during design:   .venv/bin/python runs/cap-radius/builds/swept55/blunt/radius_map.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cap_model as cm  # noqa: E402

from OCP.GeomLProp import GeomLProp_SLProps  # noqa: E402
from OCP.gp import gp_Dir  # noqa: E402


def t_att(xi):
    """Parent attachment thickness (mm) at chord fraction xi (NaN outside [0, 1])."""
    xi = np.asarray(xi, float)
    out = np.full(xi.shape, np.nan)
    ok = (xi >= 0) & (xi <= 1)
    out[ok] = 2 * cm.Zp(xi[ok], cm.ATT)
    return out


def surface_props(srf, U, V):
    """Points, unit normals and principal curvatures (k1 >= k2, signed w.r.t. the normal),
    principal directions, on parameter arrays U, V (same shape)."""
    U = np.asarray(U, float); V = np.asarray(V, float)
    shp = U.shape
    P = np.zeros(shp + (3,)); N = np.full(shp + (3,), np.nan)
    K1 = np.full(shp, np.nan); K2 = np.full(shp, np.nan)
    D1 = np.full(shp + (3,), np.nan)
    props = GeomLProp_SLProps(srf, 2, 1e-12)
    d = gp_Dir()
    for idx in np.ndindex(shp):
        props.SetParameters(float(U[idx]), float(V[idx]))
        p = props.Value()
        P[idx] = (p.X(), p.Y(), p.Z())
        if props.IsNormalDefined():
            n = props.Normal(); N[idx] = (n.X(), n.Y(), n.Z())
        if props.IsCurvatureDefined():
            K1[idx] = props.MaxCurvature(); K2[idx] = props.MinCurvature()
            props.CurvatureDirections(d, gp_Dir())
            D1[idx] = (d.X(), d.Y(), d.Z())
    return P, N, K1, K2, D1


def judge(P, K1, K2):
    """Radius, chord fraction, attachment thickness and the rule margin at each point.

    thick (t >= 0.8): need R >= 0.4  -> margin = R / 0.4
    thin  (t <  0.8): lip half-round R = t / 2 -> ratio = R / (t / 2) (reported, see README)"""
    R = 1 / np.maximum(np.abs(K1), np.abs(K2))
    xi = cm.chord_fraction_of(P[..., 0], P[..., 1])
    t = t_att(np.clip(xi, 0, 1))
    thick = t >= 0.8
    return R, xi, t, thick


def map_face(srf, nu=800, nv=240, u_skip=0.0):
    u0, u1, v0, v1 = srf.Bounds()
    s = np.linspace(0, 1, nu)
    # cluster at both ends of U (LE corner, TE half-round) and at the lip
    Us = u0 + (u1 - u0) * (0.5 - 0.5 * np.cos(np.pi * s))
    Us = Us[(Us - u0) >= u_skip]
    t = np.linspace(0, 1, nv)
    Vs = v0 + (v1 - v0) * (1 - (1 - t) ** 1.5)
    UU, VV = np.meshgrid(Us, Vs, indexing="ij")
    P, N, K1, K2, D1 = surface_props(srf, UU, VV)
    R, xi, tt, thick = judge(P, K1, K2)
    return dict(U=UU, V=VV, P=P, N=N, K1=K1, K2=K2, D1=D1, R=R, xi=xi, t=tt, thick=thick)


def summarize(m, label=""):
    R, thick, t = m["R"], m["thick"], m["t"]
    out = {}
    Rt = np.where(thick, R, np.nan)
    if np.isfinite(Rt).any():
        i = np.unravel_index(np.nanargmin(Rt), R.shape)
        out["thick"] = dict(R_min=float(R[i]), at=m["P"][i].round(4).tolist(), xi=float(m["xi"][i]),
                            t_att=float(t[i]), U=float(m["U"][i]), V=float(m["V"][i]))
    ratio = np.where(~thick & np.isfinite(t), R / (t / 2), np.nan)
    if np.isfinite(ratio).any():
        j = np.unravel_index(np.nanargmin(ratio), R.shape)
        out["thin"] = dict(R=float(R[j]), ratio_to_half_t=float(ratio[j]), at=m["P"][j].round(4).tolist(),
                           xi=float(m["xi"][j]), t_att=float(t[j]), U=float(m["U"][j]), V=float(m["V"][j]))
        Rthin = np.where(~thick & np.isfinite(t), R, np.nan)
        k = np.unravel_index(np.nanargmin(Rthin), R.shape)
        out["thin_R_min"] = dict(R=float(R[k]), at=m["P"][k].round(4).tolist(), xi=float(m["xi"][k]),
                                 t_att=float(t[k]))
    return out


if __name__ == "__main__":
    import json
    sys.path.insert(0, HERE)
    from build_cap import to_occ_surface
    z = np.load(os.path.join(HERE, "cap_patch.npz"))
    srf = to_occ_surface(z["poles"], z["ku"], z["kv"])
    m = map_face(srf, nu=int(sys.argv[1]) if len(sys.argv) > 1 else 600, nv=200)
    print(json.dumps(summarize(m), indent=1))
