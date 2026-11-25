"""Shared helpers: read the STEP, evaluate its B-spline faces with NumPy, curvature, zebra phase.

Every quantity is computed from the B-spline data stored in the reimported STEP
(poles and knots read through OCP); verify_cap.py cross-checks this evaluator
against Open CASCADE's own GeomLProp evaluation.
"""
import hashlib
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline

from OCP.STEPControl import STEPControl_Reader
from OCP.IFSelect import IFSelect_RetDone
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_REVERSED
from OCP.TopoDS import TopoDS
from OCP.BRep import BRep_Tool
from OCP.GeomAdaptor import GeomAdaptor_Surface
from OCP.GeomAbs import GeomAbs_BSplineSurface, GeomAbs_Plane


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_step(path):
    r = STEPControl_Reader()
    if r.ReadFile(str(path)) != IFSelect_RetDone:
        raise RuntimeError(f"cannot read {path}")
    r.TransferRoots()
    return r.OneShape()


def faces_of(shape):
    out = []
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        out.append(TopoDS.Face_s(ex.Current()))
        ex.Next()
    return out


class BSurf:
    """Non-rational B-spline surface evaluated with NumPy; normals are outward (face orientation applied)."""

    def __init__(self, face):
        self.face = face
        self.geom = BRep_Tool.Surface_s(face)
        bs = GeomAdaptor_Surface(self.geom).BSpline()
        assert not bs.IsURational() and not bs.IsVRational()
        self.pu, self.pv = bs.UDegree(), bs.VDegree()
        self.nu, self.nv = bs.NbUPoles(), bs.NbVPoles()
        self.C = np.array([[[bs.Pole(i, j).X(), bs.Pole(i, j).Y(), bs.Pole(i, j).Z()]
                            for j in range(1, self.nv + 1)] for i in range(1, self.nu + 1)])
        self.tu = np.concatenate([[bs.UKnot(i)] * bs.UMultiplicity(i) for i in range(1, bs.NbUKnots() + 1)])
        self.tv = np.concatenate([[bs.VKnot(i)] * bs.VMultiplicity(i) for i in range(1, bs.NbVKnots() + 1)])
        self.urange = (self.tu[0], self.tu[-1])
        self.vrange = (self.tv[0], self.tv[-1])
        self._bu = BSpline(self.tu, np.eye(self.nu), self.pu, extrapolate=True)
        self._bv = BSpline(self.tv, np.eye(self.nv), self.pv, extrapolate=True)
        self.sign = -1.0 if face.Orientation() == TopAbs_REVERSED else 1.0

    def basis(self, u, v, du, dv):
        u = np.clip(np.asarray(u, float), *self.urange)
        v = np.clip(np.asarray(v, float), *self.vrange)
        return self._bu(u, nu=du), self._bv(v, nu=dv)

    def grid(self, us, vs, du=0, dv=0):
        Bu, Bv = self.basis(us, vs, du, dv)
        return np.einsum("ia,abc,jb->ijc", Bu, self.C, Bv, optimize=True)

    def at(self, u, v, du=0, dv=0):
        u, v = np.broadcast_arrays(np.asarray(u, float), np.asarray(v, float))
        sh = u.shape
        Bu, Bv = self.basis(u.ravel(), v.ravel(), du, dv)
        return np.einsum("na,abc,nb->nc", Bu, self.C, Bv, optimize=True).reshape(sh + (3,))

    def jet(self, u, v, grid=False):
        f = self.grid if grid else self.at
        d = {k: f(u, v, *k) for k in [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]}
        return d


def unit(a):
    return a / np.linalg.norm(a, axis=-1, keepdims=True)


def differential(jet, sign):
    """Outward unit normal and fundamental forms (II with convex = positive)."""
    Su, Sv = jet[(1, 0)], jet[(0, 1)]
    cr = np.cross(Su, Sv)
    area = np.linalg.norm(cr, axis=-1)
    n = sign * cr / area[..., None]
    E = np.sum(Su * Su, -1)
    F = np.sum(Su * Sv, -1)
    G = np.sum(Sv * Sv, -1)
    L = -np.sum(jet[(2, 0)] * n, -1)
    M = -np.sum(jet[(1, 1)] * n, -1)
    N = -np.sum(jet[(0, 2)] * n, -1)
    return dict(P=jet[(0, 0)], Su=Su, Sv=Sv, n=n, area=area, E=E, F=F, G=G, L=L, M=M, N=N)


def principal(d):
    det = d["E"] * d["G"] - d["F"] ** 2
    K = (d["L"] * d["N"] - d["M"] ** 2) / det
    H = (d["E"] * d["N"] - 2 * d["F"] * d["M"] + d["G"] * d["L"]) / (2 * det)
    disc = np.sqrt(np.clip(H * H - K, 0, None))
    return H + disc, H - disc, K, H     # k1 >= k2


def tangent_coords(d, t):
    """Coefficients (a, b) with t = a Su + b Sv for tangent vectors t."""
    E, F, G = d["E"], d["F"], d["G"]
    tu = np.sum(t * d["Su"], -1)
    tv = np.sum(t * d["Sv"], -1)
    det = E * G - F * F
    return (G * tu - F * tv) / det, (E * tv - F * tu) / det


def second_form(d, t1, t2):
    a1, b1 = tangent_coords(d, t1)
    a2, b2 = tangent_coords(d, t2)
    return d["L"] * a1 * a2 + d["M"] * (a1 * b2 + b1 * a2) + d["N"] * b1 * b2


def shape_operator(d, e1, e2):
    """Shape operator (convex positive) in an orthonormal tangent frame e1, e2."""
    w11 = second_form(d, e1, e1)
    w12 = second_form(d, e1, e2)
    w22 = second_form(d, e2, e2)
    return np.stack([np.stack([w11, w12], -1), np.stack([w12, w22], -1)], -2)


def zebra_phase(n, V, dvec, width):
    """Stripe phase of the reflection of view direction V in normal n; stripes of
    'width' radians of reflected angle across direction dvec. Light when floor(phase) is even."""
    r = V - 2 * np.sum(V * n, -1, keepdims=True) * n
    s = np.clip(np.sum(r * dvec, -1), -1, 1)
    return np.arcsin(s) / width, r


def zebra_gradient(n, W, e1, e2, V, dvec, width):
    """Surface gradient (3D tangent vector, per mm) of the zebra phase, from the shape operator W in (e1, e2)."""
    r = V - 2 * np.sum(V * n, -1, keepdims=True) * n
    s = np.clip(np.sum(r * dvec, -1), -1 + 1e-12, 1 - 1e-12)
    dp_dr = dvec / (width * np.sqrt(1 - s * s))[..., None]
    Vn = np.sum(V * n, -1)
    grads = []
    for i in range(2):
        dn = W[..., i, 0, None] * e1 + W[..., i, 1, None] * e2     # derivative of n along e_i
        dr = -2 * (np.sum(V * dn, -1)[..., None] * n + Vn[..., None] * dn)
        grads.append(np.sum(dp_dr * dr, -1))
    return grads[0][..., None] * e1 + grads[1][..., None] * e2
