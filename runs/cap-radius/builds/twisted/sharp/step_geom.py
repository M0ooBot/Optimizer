"""Helpers to read wing.step / parent.step and evaluate their faces with CAD derivatives (OCP)."""
import hashlib
import pathlib

import numpy as np
from OCP.STEPControl import STEPControl_Reader
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID, TopAbs_SHELL, TopAbs_REVERSED
from OCP.TopoDS import TopoDS
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.GeomAbs import GeomAbs_BSplineSurface
from OCP.gp import gp_Pnt, gp_Vec


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_step(path):
    r = STEPControl_Reader()
    if r.ReadFile(str(path)) != 1:
        raise RuntimeError(f"cannot read {path}")
    r.TransferRoots()
    return r.OneShape()


def explore(shape, kind):
    out = []
    ex = TopExp_Explorer(shape, kind)
    while ex.More():
        out.append(ex.Current())
        ex.Next()
    return out


class Face:
    """A B-spline face with outward-oriented CAD evaluation."""

    def __init__(self, face):
        self.face = TopoDS.Face_s(face)
        self.adaptor = BRepAdaptor_Surface(self.face)
        self.surface = self.adaptor.BSpline()
        self.sign = -1.0 if self.face.Orientation() == TopAbs_REVERSED else 1.0
        self.nu, self.nv = self.surface.NbUPoles(), self.surface.NbVPoles()
        self.v0 = self.surface.VKnot(1)
        self.v1 = self.surface.VKnot(self.surface.NbVKnots())

    def poles(self):
        s = self.surface
        P = np.empty((self.nu, self.nv, 3))
        for i in range(self.nu):
            for j in range(self.nv):
                p = s.Pole(i + 1, j + 1)
                P[i, j] = (p.X(), p.Y(), p.Z())
        return P

    def value(self, u, v):
        p = self.surface.Value(float(u), float(v))
        return np.array([p.X(), p.Y(), p.Z()])

    def d2(self, u, v):
        P = gp_Pnt()
        a, b, c, d, e = gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec()
        self.surface.D2(float(u), float(v), P, a, b, c, d, e)
        f = lambda q: np.array([q.X(), q.Y(), q.Z()])
        return f(P), f(a), f(b), f(c), f(d), f(e)

    def frame(self, u, v):
        """Point, derivatives, outward unit normal and second fundamental form coefficients."""
        X, Su, Sv, Suu, Svv, Suv = self.d2(u, v)
        N = np.cross(Su, Sv) * self.sign
        N /= np.linalg.norm(N)
        return dict(X=X, Su=Su, Sv=Sv, N=N, E=Su @ Su, F=Su @ Sv, G=Sv @ Sv,
                    L=Suu @ N, M=Suv @ N, Nn=Svv @ N)


def principal_from_frame(fr):
    E, F, G, L, M, Nn = fr["E"], fr["F"], fr["G"], fr["L"], fr["M"], fr["Nn"]
    det = E * G - F * F
    K = (L * Nn - M * M) / det
    H = (E * Nn - 2 * F * M + G * L) / (2 * det)
    disc = np.sqrt(max(H * H - K, 0.0))
    k1, k2 = H + disc, H - disc
    return max(abs(k1), abs(k2)), k1, k2


def second_form_in(fr, a, b):
    """II(a, b) for 3D tangent vectors a, b (outward normal; convex => negative)."""
    Su, Sv = fr["Su"], fr["Sv"]
    J = np.array([[fr["E"], fr["F"]], [fr["F"], fr["G"]]])
    ca = np.linalg.solve(J, np.array([a @ Su, a @ Sv]))
    cb = np.linalg.solve(J, np.array([b @ Su, b @ Sv]))
    II = np.array([[fr["L"], fr["M"]], [fr["M"], fr["Nn"]]])
    return ca @ II @ cb


def wing_faces(shape):
    """Classify the six B-spline faces of wing.step."""
    faces = {}
    for f in explore(shape, TopAbs_FACE):
        F_ = Face(f)
        if F_.adaptor.GetType() != GeomAbs_BSplineSurface:
            raise RuntimeError("unexpected non-B-spline face")
        mid = F_.value(0.5, 0.5 * (F_.v0 + F_.v1))
        side = "right" if mid[1] > 0 else "left"
        if F_.v1 > 100:      # parent skin: v = |Y| in [0, 440]
            up = F_.value(0.5, 0.0)[2] > 3.0
            faces[(side, "upper" if up else "lower")] = F_
        else:                 # cap: v = w in [0, 2]
            faces[(side, "cap")] = F_
    assert len(faces) == 6, sorted(faces)
    return faces
