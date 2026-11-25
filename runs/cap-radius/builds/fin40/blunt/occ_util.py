"""Small OCP helpers shared by the verification and rendering scripts."""
import hashlib
import pathlib

import numpy as np
from OCP.STEPControl import STEPControl_Reader
from OCP.IFSelect import IFSelect_RetDone
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_REVERSED, TopAbs_IN, TopAbs_ON
from OCP.TopoDS import TopoDS
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.GeomAbs import GeomAbs_Plane
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.BRep import BRep_Tool
from OCP.ShapeAnalysis import ShapeAnalysis_Surface
from OCP.BRepTopAdaptor import BRepTopAdaptor_FClass2d
from OCP.gp import gp_Pnt, gp_Vec, gp_Pnt2d


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def read_step(path):
    r = STEPControl_Reader()
    if r.ReadFile(str(path)) != IFSelect_RetDone:
        raise RuntimeError(f"cannot read {path}")
    r.TransferRoots()
    return r.OneShape()


def bbox(shape):
    b = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape, b, False, False)
    return list(b.Get())


class Face:
    """A face with fast access to its surface, derivatives and outward normal."""

    def __init__(self, face, name=""):
        self.face = face
        self.name = name
        self.ad = BRepAdaptor_Surface(face)
        self.is_plane = self.ad.GetType() == GeomAbs_Plane
        self.surf = BRep_Tool.Surface_s(face)
        self.sas = ShapeAnalysis_Surface(self.surf)
        self.sign = -1.0 if face.Orientation() == TopAbs_REVERSED else 1.0
        self.u0, self.u1 = self.ad.FirstUParameter(), self.ad.LastUParameter()
        self.v0, self.v1 = self.ad.FirstVParameter(), self.ad.LastVParameter()
        self.bbox = bbox(face)
        self._cls = None

    def inside(self, u, v, tol=1e-9):
        if self._cls is None:
            self._cls = BRepTopAdaptor_FClass2d(self.face, tol)
        st = self._cls.Perform(gp_Pnt2d(u, v))
        return st in (TopAbs_IN, TopAbs_ON)

    def uv(self, p, tol=1e-9):
        q = self.sas.ValueOfUV(gp_Pnt(*map(float, p)), tol)
        return q.X(), q.Y(), self.sas.Gap()

    def d2(self, u, v):
        P, Su, Sv, Suu, Svv, Suv = gp_Pnt(), gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec()
        self.surf.D2(u, v, P, Su, Sv, Suu, Svv, Suv)
        a = lambda w: np.array([w.X(), w.Y(), w.Z()])
        return a(P), a(Su), a(Sv), a(Suu), a(Suv), a(Svv)

    def local(self, u, v):
        """Point, outward unit normal, tangent basis and the second fundamental form as a 3x3 tensor.

        W(t1, t2) = t1 . W . t2 is II for tangent vectors (outward normal; convex => negative).
        Also returns the Weingarten derivatives n_u, n_v.
        """
        p, su, sv, suu, suv, svv = self.d2(u, v)
        n = np.cross(su, sv)
        nn = np.linalg.norm(n)
        if nn < 1e-14:
            return None
        n = self.sign * n / nn
        I = np.array([[su @ su, su @ sv], [su @ sv, sv @ sv]])
        II = np.array([[suu @ n, suv @ n], [suv @ n, svv @ n]])
        J = np.stack([su, sv], 1)              # 3x2
        Iinv = np.linalg.inv(I)
        Jp = J @ Iinv                          # maps uv-covectors
        W = Jp @ II @ Jp.T                     # 3x3 bilinear form on tangent vectors
        A = -Iinv @ II                         # Weingarten: [n_u n_v] = J A
        nu, nv = J @ A[:, 0], J @ A[:, 1]
        S = Iinv @ II
        k = np.linalg.eigvals(S).real
        return {"p": p, "n": n, "su": su, "sv": sv, "W": W, "nu": nu, "nv": nv,
                "k": k, "kmax_abs": float(np.max(np.abs(k)))}

    def dn(self, loc, t):
        """Derivative of the unit normal along a 3D tangent vector t."""
        J = np.stack([loc["su"], loc["sv"]], 1)
        ab, *_ = np.linalg.lstsq(J, t, rcond=None)
        return ab[0] * loc["nu"] + ab[1] * loc["nv"]


def faces(shape):
    out = []
    exp = TopExp_Explorer(shape, TopAbs_FACE)
    while exp.More():
        out.append(TopoDS.Face_s(exp.Current()))
        exp.Next()
    return out


def classify(shape):
    """Name the faces of wing.step: root, te, parent_up/lo, cap_up/lo."""
    named = {}
    for f in faces(shape):
        F = Face(f)
        x0, y0, z0, x1, y1, z1 = F.bbox
        if F.is_plane:
            name = "root" if z1 < 1e-3 else "te"
        elif z1 < 240.0 + 1e-3:
            name = "parent_up" if (y0 + y1) > 0 else "parent_lo"
        else:
            name = "cap_up" if (y0 + y1) > 0 else "cap_lo"
        if name in named:
            raise RuntimeError(f"two faces classified as {name}")
        F.name = name
        named[name] = F
    return named
