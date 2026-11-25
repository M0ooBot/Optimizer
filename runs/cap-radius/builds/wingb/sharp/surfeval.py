"""Vectorised evaluation helpers for OCP surfaces (points, derivatives, principal curvatures)."""
import numpy as np


def d2(surface, u, v):
    """Point and derivatives (P, Pu, Pv, Puu, Puv, Pvv) of an OCP Geom surface on arrays u, v."""
    from OCP.gp import gp_Pnt, gp_Vec
    u, v = np.broadcast_arrays(np.asarray(u, float), np.asarray(v, float))
    out = np.zeros(u.shape + (6, 3))
    p, du, dv, duu, dvv, duv = gp_Pnt(), gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec()
    for idx in np.ndindex(u.shape):
        surface.D2(float(u[idx]), float(v[idx]), p, du, dv, duu, dvv, duv)
        out[idx] = [p.Coord(), du.Coord(), dv.Coord(), duu.Coord(), duv.Coord(), dvv.Coord()]
    return out


def curvature(D):
    """Unit normal, principal curvatures (k1 >= k2, signed w.r.t. the normal Pu x Pv) and directions."""
    P, Pu, Pv, Puu, Puv, Pvv = [D[..., i, :] for i in range(6)]
    n = np.cross(Pu, Pv)
    nn = np.linalg.norm(n, axis=-1, keepdims=True)
    n = n / np.where(nn > 0, nn, 1)
    E, F, G = (Pu * Pu).sum(-1), (Pu * Pv).sum(-1), (Pv * Pv).sum(-1)
    L, M, N = (Puu * n).sum(-1), (Puv * n).sum(-1), (Pvv * n).sum(-1)
    den = E * G - F * F
    H = (E * N - 2 * F * M + G * L) / (2 * den)
    K = (L * N - M * M) / den
    disc = np.sqrt(np.maximum(H * H - K, 0))
    k1, k2 = H + disc, H - disc
    return dict(P=P, n=n, k1=k1, k2=k2, kmax=np.maximum(abs(k1), abs(k2)), H=H, K=K,
                EFG=(E, F, G), LMN=(L, M, N), Pu=Pu, Pv=Pv)


def normal_curvature(C, t):
    """Normal curvature of the surface in tangent direction t (3-vector array, any length)."""
    Pu, Pv = C["Pu"], C["Pv"]
    E, F, G = C["EFG"]
    L, M, N = C["LMN"]
    # solve t = a Pu + b Pv in the least-squares sense
    tu, tv = (t * Pu).sum(-1), (t * Pv).sum(-1)
    den = E * G - F * F
    a = (G * tu - F * tv) / den
    b = (E * tv - F * tu) / den
    I = E * a * a + 2 * F * a * b + G * b * b
    II = L * a * a + 2 * M * a * b + N * b * b
    return II / I
