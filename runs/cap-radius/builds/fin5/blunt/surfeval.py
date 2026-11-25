"""NumPy evaluation of B-spline faces read back from STEP (exact basis derivatives).

Used by verify_cap.py. Values are cross-checked against Open CASCADE's own D2
evaluator at random parameters (reported in readback.json).
"""
import numpy as np
from scipy.interpolate import BSpline
from OCP.BRep import BRep_Tool
from OCP.GeomAdaptor import GeomAdaptor_Surface
from OCP.TopAbs import TopAbs_REVERSED
from OCP.gp import gp_Pnt, gp_Vec


class Face:
    def __init__(self, face):
        self.face = face
        self.geom = BRep_Tool.Surface_s(face)
        b = GeomAdaptor_Surface(self.geom).BSpline()
        self.b = b
        self.pu, self.pv = b.UDegree(), b.VDegree()
        nu, nv = b.NbUPoles(), b.NbVPoles()
        assert not b.IsURational() and not b.IsVRational()
        self.P = np.array([[[b.Pole(i, j).X(), b.Pole(i, j).Y(), b.Pole(i, j).Z()] for j in range(1, nv + 1)]
                           for i in range(1, nu + 1)])
        ku, kv = [], []
        for i in range(1, b.NbUKnots() + 1):
            ku += [b.UKnot(i)] * b.UMultiplicity(i)
        for i in range(1, b.NbVKnots() + 1):
            kv += [b.VKnot(i)] * b.VMultiplicity(i)
        self.ku, self.kv = np.array(ku), np.array(kv)
        self.urange = (self.ku[0], self.ku[-1])
        self.vrange = (self.kv[0], self.kv[-1])
        self.sign = -1.0 if face.Orientation() == TopAbs_REVERSED else 1.0
        self._bu = [BSpline(self.ku, np.eye(nu), self.pu)]
        self._bv = [BSpline(self.kv, np.eye(nv), self.pv)]
        for k in (1, 2):
            self._bu.append(self._deriv(self._bu[0], k, nu))
            self._bv.append(self._deriv(self._bv[0], k, nv))

    @staticmethod
    def _deriv(b, k, n):
        if k <= b.k:
            return b.derivative(k)
        return lambda x: np.zeros(np.shape(x) + (n,))

    def basis(self, u, v):
        u = np.clip(np.asarray(u, float), *self.urange)
        v = np.clip(np.asarray(v, float), *self.vrange)
        return [b(u) for b in self._bu], [b(v) for b in self._bv]

    def grid(self, u, v):
        """Derivatives on the tensor grid u x v: dict of arrays (len(u), len(v), 3)."""
        Bu, Bv = self.basis(u, v)
        nu, nv, _ = self.P.shape
        T = [(Bu[a] @ self.P.reshape(nu, nv * 3)).reshape(-1, nv, 3) for a in range(3)]
        f = lambda a, b: np.einsum("abk,cb->ack", T[a], Bv[b], optimize=True)
        return {"S": f(0, 0), "Su": f(1, 0), "Sv": f(0, 1), "Suu": f(2, 0), "Suv": f(1, 1), "Svv": f(0, 2)}

    def points(self, u, v):
        """Derivatives at scattered (u, v): dict of arrays (n, 3)."""
        u = np.atleast_1d(np.asarray(u, float)); v = np.atleast_1d(np.asarray(v, float))
        Bu, Bv = self.basis(u, v)
        f = lambda a, b: np.einsum("ni,ijk,nj->nk", Bu[a], self.P, Bv[b], optimize=True)
        return {"S": f(0, 0), "Su": f(1, 0), "Sv": f(0, 1), "Suu": f(2, 0), "Suv": f(1, 1), "Svv": f(0, 2)}

    def occ_check(self, n=200, seed=1):
        rng = np.random.default_rng(seed)
        u = rng.uniform(*self.urange, n); v = rng.uniform(*self.vrange, n)
        d = self.points(u, v)
        err = 0.0
        for i in range(n):
            p = gp_Pnt(); a, b, c, e, f, g = gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec()
            self.geom.D2(float(u[i]), float(v[i]), p, a, b, c, e, f)
            occ = [p.Coord(), a.Coord(), b.Coord(), c.Coord(), e.Coord(), f.Coord()]
            mine = [d["S"][i], d["Su"][i], d["Sv"][i], d["Suu"][i], d["Svv"][i], d["Suv"][i]]
            for o, m in zip(occ, mine):
                err = max(err, np.abs(np.array(o) - m).max() / max(1.0, np.abs(m).max()))
        return err


def frame(d, sign=1.0):
    """Outward unit normal, orthonormal tangent frame (e1 along Su), and the second
    fundamental form B in that frame (normal curvature k_n(T) = T^T B T, convex < 0)."""
    Su, Sv = d["Su"], d["Sv"]
    n = np.cross(Su, Sv)
    n = sign * n / np.linalg.norm(n, axis=-1, keepdims=True)
    e1 = Su / np.linalg.norm(Su, axis=-1, keepdims=True)
    e2 = np.cross(n, e1)
    J = np.stack([np.stack([np.einsum("...k,...k", e1, Su), np.einsum("...k,...k", e1, Sv)], -1),
                  np.stack([np.einsum("...k,...k", e2, Su), np.einsum("...k,...k", e2, Sv)], -1)], -2)
    L = np.einsum("...k,...k", d["Suu"], n); M = np.einsum("...k,...k", d["Suv"], n); N = np.einsum("...k,...k", d["Svv"], n)
    II = np.stack([np.stack([L, M], -1), np.stack([M, N], -1)], -2)
    Ji = np.linalg.inv(J)
    B = np.swapaxes(Ji, -1, -2) @ II @ Ji
    return n, e1, e2, B


def principal(B):
    ev = np.linalg.eigvalsh(0.5 * (B + np.swapaxes(B, -1, -2)))
    return ev  # ascending; convex surfaces have negative values
