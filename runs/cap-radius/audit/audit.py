"""Audit a whole rounded-tip model STEP against its locked parent and the 0.4 mm printable-radius rule.

    python runs/cap-radius/audit/audit.py SHAPE VARIANT MODEL.step [--out FILE.json]

SHAPE is one of flat, twisted, swept55, fin40, fin5, wingb and VARIANT is sharp or blunt; the parent
is read from runs/cap-radius/parents/SHAPE/VARIANT/ (parent.py, parent.step, spec.json). Every number
comes from the reimported STEP through CAD surface derivatives (BRepAdaptor_Surface.D2); no mesh is
used. The JSON holds every measurement; a short text summary is printed and written next to it.
See README.md in this folder for the definitions, the cut families and the known limits.
"""

import argparse
import hashlib
import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq, minimize, minimize_scalar

from OCP.Bnd import Bnd_Box
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Curve2d, BRepAdaptor_Surface
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepCheck import BRepCheck_Analyzer, BRepCheck_NoError, BRepCheck_Shell
from OCP.BRepGProp import BRepGProp
from OCP.BRepTools import BRepTools
from OCP.BRepTopAdaptor import BRepTopAdaptor_FClass2d
from OCP.GeomAbs import GeomAbs_Plane
from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.ShapeAnalysis import ShapeAnalysis_Surface
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_IN, TopAbs_ON, TopAbs_REVERSED, TopAbs_SHELL, TopAbs_SOLID, TopAbs_VERTEX
from OCP.TopExp import TopExp, TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape, TopTools_IndexedMapOfShape
from OCP.gp import gp_Dir, gp_Pln, gp_Pnt, gp_Pnt2d, gp_Vec

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PARENTS = REPO / "runs" / "cap-radius" / "parents"
SHAPES = ("flat", "twisted", "swept55", "fin40", "fin5", "wingb")
VARIANTS = ("sharp", "blunt")

# The rule under test (skill step 5) and the audit's own acceptance tolerances.
RULE = {"min_radius_mm": 0.4, "thick_gate_mm": 0.8, "thin_ratio_pass": 0.9}
VERDICT_RTOL = 1e-6            # radius and ratio verdicts allow this relative round-off (a half-round at t = 0.8 is 0.4)
TOL = {
    "extent_mm": 1e-3,          # projected span / height against spec.json
    "retained_skin_mm": 1e-4,   # retained skin against the parent, both directions
    "outline_mm": 1e-3,         # outermost lip against the locked outline
    "join_gap_mm": 1e-4,        # G0 at the parent-to-cap joins
    "join_angle_deg": 0.01,     # G1 at the parent-to-cap joins
    "te_point_mm": 1e-3,        # sharp TE: the cap reaches the parent's TE point
}
SECTION_FRACTIONS = (0.02, 0.05, 0.1, 0.25, 0.4, 0.5, 0.65, 0.8, 0.9, 0.95, 0.97, 0.99)
THIN_TARGETS_TE = {"sharp": (0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05), "blunt": (0.7, 0.6, 0.5, 0.45)}
THIN_TARGETS_LE = (0.6, 0.4)
TE_APPROACH_MM = (1.0, 0.5, 0.2, 0.1, 0.05, 0.02)
CURVATURE_DEADBAND = 1e-5      # 1/mm; |kappa| below this counts as straight when looking for reversals
TURN_STEP_DEG = 1.5            # adaptive trace sampling: largest tangent turn between neighbouring samples
DEGENERATE_SIN = 1e-6          # |Su x Sv| / (|Su| |Sv|) below this: derivatives treated as undefined
DEGENERATE_AREA = 1e-9         # ... or |Su x Sv| below this fraction of the face's typical |Su| |Sv| (collapsed sides)
COLLAPSE_EXCLUSION_MM = 0.05   # verdicts on the map and joins skip points this close to a collapsed face corner (reported separately)
SWEEP55 = math.tan(math.radians(55.0))


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def clean(x):
    """Make a structure JSON-safe (numpy scalars and arrays, NaN/inf as None)."""
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        x = float(x)
        return x if math.isfinite(x) else None
    return x


# ---------------------------------------------------------------------------
# Parent and cut family
# ---------------------------------------------------------------------------

def load_parent(shape, variant):
    path = PARENTS / shape / variant / "parent.py"
    name = f"audit_parent_{shape}_{variant}"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class Tip:
    """One cap: outboard direction, attachment station and the cut family through it."""

    def __init__(self, case, name, span, direction, outboard):
        self.case, self.name, self.span, self.direction = case, name, span, direction
        self.outboard = np.asarray(outboard, float)   # unit vector; outboard coordinate = P . outboard
        self.station = case.station
        m = case.parent
        h = 0.5
        frames = [m.section_frame(span + direction * d) for d in (-h, 0.0, h)]
        A = [np.asarray(f["origin"], float) for f in frames]
        B = [float(f["chord"]) * np.asarray(f["chord_direction"], float) for f in frames]
        self.origin0, self.chordvec0 = A[1], B[1]
        self.chord = float(frames[1]["chord"])
        self.chord_dir = unit(B[1])
        self.up = unit(np.asarray(frames[1]["up_direction"], float))
        self.Sa = (A[2] - A[0]) / (2 * h)
        self.Sb = (B[2] - B[0]) / (2 * h)
        self.family = case.family

    def planes(self, f):
        """Vectorised cut planes at chord fractions f: point O, unit normal n, thickness dir T, outboard w."""
        f = np.atleast_1d(np.asarray(f, float))[:, None]
        T = np.broadcast_to(self.up, (f.shape[0], 3))
        if self.family == "legacy":
            return self._legacy(f[:, 0])
        O = self.origin0 + f * self.chordvec0
        S = self.Sa + f * self.Sb    # tangent of the parent's constant-chord-fraction line: follows its sweep, taper, dihedral
        n = unit(np.cross(S, T))
        n = n * np.sign(n @ self.chord_dir)[:, None]
        w = unit(S - np.sum(S * T, axis=1)[:, None] * T)
        w = w * np.sign(w @ self.outboard)[:, None]
        return O, n, T, w

    def _legacy(self, f):
        """The planes of runs/cap-trials/retest/audit.py (calibration only)."""
        shape, s = self.case.shape, self.direction
        k = f.shape[0]
        if shape.startswith("fin"):
            slope = (65 + 85 * f - 160 * f) / 240
            O = np.stack([160 * f, 0 * f, 0 * f], 1)
            n = unit(np.stack([np.ones(k), 0 * f, -slope], 1))
            T = np.broadcast_to([0.0, 1.0, 0.0], (k, 3))
            S = np.stack([slope, 0 * f, np.ones(k)], 1)
        else:
            if shape == "flat":
                slope, root = 0 * f, 120 * f
            elif shape == "swept55":
                slope, root = SWEEP55 + 0 * f, 95 * f
            elif shape == "twisted":
                chord = 100 * math.cos(math.radians(3))
                slope, root = (90 + chord * f - 180 * f) / 440, 180 * f
            else:
                raise ValueError("no legacy planes for " + shape)
            O = np.stack([root, 0 * f, 0 * f], 1)
            n = unit(np.stack([np.ones(k), -s * slope, 0 * f], 1))
            T = np.broadcast_to([0.0, 0.0, 1.0], (k, 3))
            S = np.stack([slope, s * np.ones(k), 0 * f], 1)
        w = unit(S - np.sum(S * T, axis=1)[:, None] * T)
        w = w * np.sign(w @ self.outboard)[:, None]
        return O, n, T, w

    def plane(self, f):
        O, n, T, w = self.planes([f])
        return O[0], n[0], T[0], w[0]

    def span_direction(self, f):
        """Unit outboard direction of the cut family inside the plane, keeping its rise (dihedral) for 3D rims."""
        O, n, T, w = self.plane(f)
        if self.family == "legacy" or self.case.outline_kind != "rim":
            return w
        S = self.Sa + f * self.Sb
        S = S - (S @ n) * n
        return unit(S) * np.sign(unit(S) @ self.outboard)

    def fraction_of(self, P):
        """Chord fraction of the family plane through each point (NaN outside -0.2..1.2)."""
        P = np.atleast_2d(P)
        grid = np.linspace(-0.2, 1.2, 141)
        O, n, _, _ = self.planes(grid)
        out = np.full(len(P), np.nan)
        for start in range(0, len(P), 20000):
            Q = P[start:start + 20000]
            D = np.einsum("mkj,kj->mk", Q[:, None, :] - O[None], n)
            sign = D[:, :-1] * D[:, 1:] <= 0
            has = sign.any(1)
            k = np.argmax(sign, axis=1)
            d0 = D[np.arange(len(Q)), k]
            d1 = D[np.arange(len(Q)), k + 1]
            lo, hi = grid[k], grid[k + 1]
            f = np.where(np.abs(d1 - d0) > 0, lo + (hi - lo) * d0 / (d0 - d1), lo)
            for _ in range(4):   # secant polish on the exact family
                Of, nf, _, _ = self.planes(f)
                g = np.sum((Q - Of) * nf, 1)
                Of2, nf2, _, _ = self.planes(f + 1e-7)
                g2 = np.sum((Q - Of2) * nf2, 1)
                slope = (g2 - g) / 1e-7
                f = f - np.where(np.abs(slope) > 0, g / slope, 0)
            out[start:start + len(Q)] = np.where(has, f, np.nan)
        return out


class Case:
    def __init__(self, shape, variant, family="native"):
        if family == "legacy" and shape == "wingb":
            family = "native"
        self.shape, self.variant, self.family = shape, variant, family
        self.parent = m = load_parent(shape, variant)
        self.spec = json.loads((PARENTS / shape / variant / "spec.json").read_text())
        self.blunt = m.TE_THICKNESS > 0
        self.station = float(m.ATTACHMENT)
        if shape in ("flat", "twisted", "swept55"):
            self.tips = [Tip(self, "right (+Y)", +self.station, +1, [0, 1, 0]),
                         Tip(self, "left (-Y)", -self.station, -1, [0, -1, 0])]
        elif shape.startswith("fin"):
            self.tips = [Tip(self, "top (+Z)", self.station, +1, [0, 0, 1])]
        else:
            self.tips = [Tip(self, "tip (-X)", self.station, +1, [-1, 0, 0])]
        sec = self.spec["attachment_section"]
        self.f_te = float(sec["trailing_edge_chord_fraction"])
        self.f_le = max(0.0, float(sec["leading_point_chord_fraction"]))
        self.crossings = [float(x) for x in self.spec["attachment_thickness"]["chord_fractions_where_thickness_is_0p8_mm"]]
        self._thickness_table()

    # -- attachment thickness -------------------------------------------------
    def thickness_exact(self, f):
        f = np.atleast_1d(np.asarray(f, float))
        return np.asarray(self.parent.section_thickness(f, self.tips[0].span), float).reshape(f.shape)

    def _thickness_table(self):
        k = np.linspace(0, np.pi, 1601)
        f = self.f_le + (self.f_te - self.f_le) * 0.5 * (1 - np.cos(k))
        f = np.unique(np.concatenate([f, self.crossings, [self.f_le, self.f_te]]))
        t = self.thickness_exact(f)
        good = np.isfinite(t)
        if self.variant == "sharp":
            t = np.where(good, t, 0.0)
            good = np.ones_like(good)
        self._tf = PchipInterpolator(f[good], t[good], extrapolate=False)
        self._tf_range = (f[good][0], f[good][-1])

    def thickness(self, f):
        f = np.asarray(f, float)
        lo, hi = self._tf_range
        return self._tf(np.clip(f, lo, hi))

    def fraction_at_thickness(self, target, side):
        lo, hi = self._tf_range
        if side == "te":
            a, b = self.crossings[-1], hi
        else:
            a, b = lo, self.crossings[0]
        g = lambda x: float(self.thickness_exact(x)[0]) - target
        try:
            ga, gb = g(a), g(b)
            if not (np.isfinite(ga) and np.isfinite(gb)) or ga * gb > 0:
                return None
            return brentq(g, a, b, xtol=1e-15, rtol=1e-15)
        except ValueError:
            return None

    def required_radius(self, t):
        t = np.asarray(t, float)
        return np.where(t >= RULE["thick_gate_mm"], RULE["min_radius_mm"], RULE["thin_ratio_pass"] * 0.5 * t)

    # -- locked outline -----------------------------------------------------------
    def outline(self, p, tip):
        """Locked outline points (N, 3); for Wing B a plan curve, so Y is NaN."""
        m, p = self.parent, np.atleast_1d(np.asarray(p, float))
        if self.shape in ("flat", "twisted", "swept55"):
            return np.asarray(m.tip_outline(p, tip.direction), float).reshape(-1, 3)
        if self.shape.startswith("fin"):
            h = np.where(p <= 0.5, 2 * p, 2 * (1 - p))
            a = np.asarray(m.tip_outline(h, -1), float)
            b = np.asarray(m.tip_outline(h, +1), float)
            return np.where((p <= 0.5)[:, None], a, b)
        xz = np.asarray(m.tip_outline(p), float).reshape(-1, 2)
        return np.stack([xz[:, 0], np.full(len(p), np.nan), xz[:, 1]], 1)

    @property
    def outline_is_plan(self):
        return self.shape == "wingb"

    @property
    def outline_kind(self):
        """plan: plan curve (Wing B); rim: 3D rim on the cambered mid-surface (twisted);
        midplane: curve in the symmetry plane of a symmetric section (flat, swept55, fins)."""
        return {"wingb": "plan", "twisted": "rim"}.get(self.shape, "midplane")

    def outline_in_plane(self, tip, O, n, near=None):
        p = np.linspace(0, 1, 4001)
        Q = self.outline(p, tip)
        Qz = np.nan_to_num(Q)
        d = (Qz - O) @ n
        idx = np.flatnonzero(d[:-1] * d[1:] <= 0)
        roots = []
        for i in idx:
            g = lambda x: float((np.nan_to_num(self.outline([x], tip)[0]) - O) @ n)
            try:
                roots.append(brentq(g, p[i], p[i + 1], xtol=1e-15))
            except ValueError:
                pass
        if not roots:
            return None, None
        pts = self.outline(roots, tip)
        if near is not None and len(roots) > 1:
            dist = np.linalg.norm(np.nan_to_num(pts - near) * (0 if self.outline_is_plan else 1)
                                  + np.nan_to_num(pts - near) * (np.array([1, 0, 1]) if self.outline_is_plan else 0), axis=1)
            k = int(np.argmin(dist))
        else:
            k = 0
        return roots[k], pts[k]

    # -- parent samples -----------------------------------------------------------
    def parent_skin_points(self):
        m = self.parent
        lo, hi = m.SPAN_LIMITS
        s = 0.5 * (1 - np.cos(np.linspace(0, np.pi, 61)))
        spans = np.unique(np.concatenate([np.linspace(lo, hi, 41), [lo + 1e-3, hi - 1e-3]]))
        pts = []
        if self.shape == "wingb":
            spans = spans[(spans > 0.2) & (spans <= 540)]
            for y in spans:
                for side in (+1, -1):
                    a, b = m.side_range(y, side)
                    u = float(a) + (float(b) - float(a)) * s
                    pts.append(np.asarray(m.skin(u, np.full_like(u, y)), float))
        else:
            S, Y = np.meshgrid(s, spans, indexing="ij")
            for side in (+1, -1):
                pts.append(np.asarray(m.skin(S, Y, side), float).reshape(-1, 3))
        if self.blunt:
            tt = np.linspace(0, 1, 9)
            T_, Y_ = np.meshgrid(tt, spans, indexing="ij")
            pts.append(np.asarray(m.te_face(T_, Y_), float).reshape(-1, 3))
        return np.concatenate(pts)


# ---------------------------------------------------------------------------
# STEP model
# ---------------------------------------------------------------------------

def read_step(path):
    reader = STEPControl_Reader()
    if reader.ReadFile(str(path)) != IFSelect_RetDone:
        raise IOError(f"cannot read {path}")
    reader.TransferRoots()
    return reader.OneShape()


def shapes_of(shape, kind):
    m = TopTools_IndexedMapOfShape()
    TopExp.MapShapes_s(shape, kind, m)
    return [m.FindKey(i) for i in range(1, m.Extent() + 1)], m


class Face:
    def __init__(self, index, face):
        self.index, self.face = index, face
        self.ad = BRepAdaptor_Surface(face)
        self.sign = -1.0 if face.Orientation() == TopAbs_REVERSED else 1.0
        self.type = str(self.ad.GetType()).split(".")[-1].replace("GeomAbs_", "")
        self.uv = BRepTools.UVBounds_s(face)
        self.cls = BRepTopAdaptor_FClass2d(face, 1e-9)
        self._p, self._v = gp_Pnt(), [gp_Vec() for _ in range(5)]
        self.role, self.tip = None, None
        self._sas = None
        u0, u1, v0, v1 = self.uv
        prods = []
        for u in np.linspace(u0, u1, 7)[1:-1]:
            for v in np.linspace(v0, v1, 7)[1:-1]:
                _, Su, Sv = self.d2(u, v)[:3]
                prods.append(np.linalg.norm(Su) * np.linalg.norm(Sv))
        self.scale = float(np.median(prods)) if prods else 1.0
        self.planar = self.ad.GetType() == GeomAbs_Plane or self._flat()

    def _flat(self):
        """A B-spline (or other) face that is geometrically a plane: one normal and coplanar points."""
        u0, u1, v0, v1 = self.uv
        pts, nrm = [], []
        for u in np.linspace(u0, u1, 7):
            for v in np.linspace(v0, v1, 7):
                P, Su, Sv = self.d2(u, v)[:3]
                c = np.cross(Su, Sv)
                pts.append(P)
                if np.linalg.norm(c) > 0:
                    nrm.append(c / np.linalg.norm(c))
        if len(nrm) < 3:
            return False
        nrm, pts = np.array(nrm), np.array(pts)
        n0 = nrm[0]
        if np.max(np.linalg.norm(np.cross(nrm, n0), axis=1)) > 1e-9:
            return False
        return bool(np.max(np.abs((pts - pts[0]) @ n0)) < 1e-9)

    def inside(self, u, v):
        return self.cls.Perform(gp_Pnt2d(u, v)) in (TopAbs_IN, TopAbs_ON)

    def d2(self, u, v):
        a, b, c, d, e = self._v
        self.ad.D2(float(u), float(v), self._p, a, b, c, d, e)
        return (np.array(self._p.Coord()), np.array(a.Coord()), np.array(b.Coord()),
                np.array(c.Coord()), np.array(e.Coord()), np.array(d.Coord()))   # P, Su, Sv, Suu, Suv, Svv

    def point(self, u, v):
        return np.array(self.ad.Value(float(u), float(v)).Coord())

    def jet(self, u, v):
        """Point, outward unit normal, ambient second-fundamental tensor H (3x3, n . d2r convention),
        principal curvatures (convex positive) and the derivative condition |Su x Sv| / (|Su||Sv|)."""
        P, Su, Sv, Suu, Suv, Svv = self.d2(u, v)
        cr = np.cross(Su, Sv)
        area = np.linalg.norm(cr)
        cond = area / max(np.linalg.norm(Su) * np.linalg.norm(Sv), 1e-300)
        if area == 0 or cond < DEGENERATE_SIN or area < DEGENERATE_AREA * self.scale:
            return dict(P=P, N=None, H=None, k=None, cond=cond)
        N = self.sign * cr / area
        I = np.array([[Su @ Su, Su @ Sv], [Su @ Sv, Sv @ Sv]])
        II = np.array([[N @ Suu, N @ Suv], [N @ Suv, N @ Svv]])
        J = np.column_stack((Su, Sv))
        Jp = np.linalg.solve(I, J.T)          # 2x3 pseudo-inverse
        H = Jp.T @ II @ Jp
        W = np.linalg.solve(I, II)
        k = -np.real(np.linalg.eigvals(W))    # convex (bending away from the outward normal) positive
        return dict(P=P, N=N, H=H, k=np.sort(k), cond=cond, Su=Su, Sv=Sv)

    def surface(self):
        return BRep_Tool.Surface_s(self.face)

    def nearest_distance(self, p):
        """Distance from p to the face's surface: nearest points of a cached (u, v) grid, then a bounded
        foot-point minimisation from each, plus OCCT's global projector as one more candidate."""
        from scipy.spatial import cKDTree
        u0, u1, v0, v1 = self.uv
        if not hasattr(self, "_tree"):
            us, vs = self.grid(121, 61)
            uv = np.array([(u, v) for u in us for v in vs])
            self._tree_uv = uv
            self._tree = cKDTree(np.array([self.point(u, v) for u, v in uv]))
        p = np.asarray(p, float)
        best = np.inf
        _, idx = self._tree.query(p, k=4)
        su, sv = (u1 - u0) or 1.0, (v1 - v0) or 1.0
        for k in np.atleast_1d(idx):
            x0 = [(self._tree_uv[k][0] - u0) / su, (self._tree_uv[k][1] - v0) / sv]

            def obj(x):
                P, Su, Sv = self.d2(u0 + x[0] * su, v0 + x[1] * sv)[:3]
                r = P - p
                return r @ r, np.array([2 * r @ Su * su, 2 * r @ Sv * sv])
            res = minimize(obj, x0, jac=True, bounds=[(0, 1), (0, 1)], method="L-BFGS-B",
                           options={"ftol": 1e-30, "gtol": 1e-20, "maxiter": 500})
            best = min(best, math.sqrt(max(res.fun, 0.0)))
        pr = GeomAPI_ProjectPointOnSurf(gp_Pnt(*map(float, p)), self.surface())
        if pr.NbPoints():
            best = min(best, pr.LowerDistance())
        return best

    def analysis(self):
        if self._sas is None:
            self._sas = ShapeAnalysis_Surface(self.surface())
        return self._sas

    def collapsed(self):
        if not hasattr(self, "_collapsed"):
            self._collapsed = collapsed_points(self)
        return self._collapsed

    def distance_to_collapse(self, P):
        pts = [c for _, c in self.collapsed()]
        if not pts:
            return np.full(len(np.atleast_2d(P)), np.inf)
        P = np.atleast_2d(P)
        return np.min(np.linalg.norm(P[:, None, :] - np.array(pts)[None], axis=2), axis=1)

    def grid(self, n_uniform, n_cos):
        u0, u1, v0, v1 = self.uv
        k = np.linspace(0, np.pi, n_cos)
        us = np.unique(np.concatenate([np.linspace(u0, u1, n_uniform), u0 + (u1 - u0) * 0.5 * (1 - np.cos(k))]))
        vs = np.unique(np.concatenate([np.linspace(v0, v1, n_uniform), v0 + (v1 - v0) * 0.5 * (1 - np.cos(k))]))
        return us, vs


def project_distance(points, faces, margin=5.0):
    """Smallest distance from each point to the (untrimmed) surfaces of the given faces."""
    best = np.full(len(points), np.inf)
    for fr in faces:
        box = Bnd_Box()
        BRepBndLib.Add_s(fr.face, box)
        lo = np.array(box.CornerMin().Coord()) - margin
        hi = np.array(box.CornerMax().Coord()) + margin
        inside = np.flatnonzero(np.all((points >= lo) & (points <= hi), axis=1))
        sas = fr.analysis()
        for i in inside:
            p = gp_Pnt(*map(float, points[i]))
            sas.ValueOfUV(p, 1e-10)
            d = sas.Gap()
            if not np.isfinite(d) or d > 1e-6:   # the projector can stall where a skin closes to a thin TE
                d = min(d, fr.nearest_distance(points[i]))
            best[i] = min(best[i], d)
    for i in np.flatnonzero(~np.isfinite(best)):   # far from every face's box
        for fr in faces:
            best[i] = min(best[i], fr.nearest_distance(points[i]))
    return best


class Model:
    def __init__(self, path, case):
        self.path, self.case = Path(path), case
        self.shape = read_step(path)
        self.solids, _ = shapes_of(self.shape, TopAbs_SOLID)
        self.shells, _ = shapes_of(self.shape, TopAbs_SHELL)
        faces, self.face_map = shapes_of(self.shape, TopAbs_FACE)
        self.faces = [Face(i, TopoDS.Face_s(f)) for i, f in enumerate(faces)]
        self.classify()

    def classify(self):
        case = self.case
        for fr in self.faces:
            u0, u1, v0, v1 = fr.uv
            pts = [fr.point(u, v) for u in np.linspace(u0, u1, 9) for v in np.linspace(v0, v1, 9) if fr.inside(u, v)]
            if not pts:
                pts = [fr.point(0.5 * (u0 + u1), 0.5 * (v0 + v1))]
            pts = np.array(pts)
            fr.role, fr.tip = "retained", None
            for tip in case.tips:
                o = pts @ tip.outboard - tip.station
                if o.max() > 1e-3:
                    fr.tip = tip
                    fr.role = "cap" if o.min() > -1e-4 else "straddles attachment"
            fr.centre = pts.mean(0)

    def cap_faces(self, tip=None):
        return [f for f in self.faces if f.role != "retained" and (tip is None or f.tip is tip)]

    def retained_faces(self):
        return [f for f in self.faces if f.role == "retained"]


# ---------------------------------------------------------------------------
# 1. Validity, extents, volume
# ---------------------------------------------------------------------------

def volume(shape, tolerance=None):
    props = GProp_GProps()
    if tolerance is None:
        BRepGProp.VolumeProperties_s(shape, props)
    else:
        BRepGProp.VolumeProperties_s(shape, props, tolerance, True, False)
    return props.Mass()


def surface_extreme(model, direction, faces=None):
    """Largest P . direction over the faces, from a seeded grid refined by bounded optimisation in (u, v)."""
    d = np.asarray(direction, float)
    best, where = -np.inf, None
    for fr in (faces or model.faces):
        u0, u1, v0, v1 = fr.uv
        seeds = []
        for u in np.linspace(u0, u1, 25):
            for v in np.linspace(v0, v1, 25):
                if fr.inside(u, v):
                    seeds.append((fr.point(u, v) @ d, u, v))
        seeds.sort(reverse=True)
        for val, u, v in seeds[:3]:
            obj = lambda x: -(fr.point(u0 + x[0] * (u1 - u0), v0 + x[1] * (v1 - v0)) @ d)
            x0 = [(u - u0) / (u1 - u0), (v - v0) / (v1 - v0)]
            r = minimize(obj, x0, bounds=[(0, 1), (0, 1)], method="L-BFGS-B", options={"ftol": 1e-16, "gtol": 1e-12})
            uu, vv = u0 + r.x[0] * (u1 - u0), v0 + r.x[1] * (v1 - v0)
            cand = -r.fun if fr.inside(uu, vv) else val
            if cand > best:
                best, where = cand, fr.point(uu, vv) if fr.inside(uu, vv) else fr.point(u, v)
    return best, where


def validity(model, case):
    shape = model.shape
    solid = model.solids[0] if model.solids else shape
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape, box, False, False)
    xmin, ymin, zmin, xmax, ymax, zmax = box.Get()
    spec_extent = case.spec["extent"]
    ends = {}
    for tip in case.tips:
        ends[tip.name] = surface_extreme(model, tip.outboard, model.cap_faces(tip) or None)
    if case.shape in ("flat", "twisted", "swept55"):
        (a, pa), (b, pb) = ends.values()
        measured, key, label = a + b, "finished_projected_span_mm", "projected span (Y)"
    elif case.shape.startswith("fin"):
        measured, key, label = list(ends.values())[0][0] - zmin, "finished_height_mm", "height (Z) above the root"
    else:
        measured, key, label = list(ends.values())[0][0] - (-xmax), "finished_projected_half_span_mm", "projected half-span from the root plane (-X)"
    target = float(spec_extent[key])
    out = {
        "solids": len(model.solids),
        "shells": len(model.shells),
        "faces": len(model.faces),
        "face_types": [f.type for f in model.faces],
        "face_roles": [f"{f.index}:{f.role}" + (f"@{f.tip.name}" if f.tip else "") for f in model.faces],
        "shell_closed": bool(model.shells) and all(BRepCheck_Shell(TopoDS.Shell_s(s)).Closed() == BRepCheck_NoError for s in model.shells),
        "brepcheck_valid": bool(BRepCheck_Analyzer(solid).IsValid()),
        "brepcheck_strict_valid": bool(BRepCheck_Analyzer(solid, True, False, True).IsValid()),
        "bounding_box_mm": [xmin, xmax, ymin, ymax, zmin, zmax],
        "extent": {"what": label, "measured_mm": measured, "spec_mm": target, "difference_mm": measured - target,
                   "outermost_points": {k: v[1] for k, v in ends.items()}, "bounding_box_note": "AddOptimal box shown for reference; extents come from surface optimisation"},
        "volume_mm3": volume(solid, 1e-9),
        "volume_occt_default_mm3": volume(solid),
        "straddling_faces": [f.index for f in model.faces if f.role == "straddles attachment"],
    }
    check = json.loads((PARENTS / case.shape / case.variant / "check.json").read_text())
    out["parent_volume_mm3"] = check["reimport"]["volume_mm3"]
    out["cap_volume_mm3"] = out["volume_mm3"] - out["parent_volume_mm3"]
    out["pass"] = (out["solids"] == 1 and out["shells"] == 1 and out["shell_closed"] and out["brepcheck_valid"]
                   and out["brepcheck_strict_valid"] and not out["straddling_faces"])
    out["extent"]["pass"] = abs(measured - target) <= TOL["extent_mm"]
    return out


# ---------------------------------------------------------------------------
# 2. Parent preservation and locked outline
# ---------------------------------------------------------------------------

def parent_faces(case):
    shape = read_step(PARENTS / case.shape / case.variant / "parent.step")
    faces, _ = shapes_of(shape, TopAbs_FACE)
    out = []
    for i, f in enumerate(faces):
        fr = Face(i, TopoDS.Face_s(f))
        u0, u1, v0, v1 = fr.uv
        c = fr.point(0.5 * (u0 + u1), 0.5 * (v0 + v1))
        if fr.planar and any(abs(c @ t.outboard - t.station) < 1e-6 for t in case.tips):
            continue    # the parent's flat attachment faces are replaced by the cap
        out.append(fr)
    return out


def preservation(model, case):
    retained = model.retained_faces()
    pts = case.parent_skin_points()
    d = project_distance(pts, retained)
    k = int(np.argmax(d))
    forward = {"samples": len(pts), "max_mm": float(d.max()), "at_xyz": pts[k],
               "what": "parent.py skin (and blunt TE face) points to the model's retained faces"}
    # Reverse: points on the model's retained faces to the parent's own surfaces.
    rev_pts = []
    for fr in retained:
        us, vs = fr.grid(31, 0)
        for u in us:
            for v in vs:
                if fr.inside(u, v):
                    rev_pts.append(fr.point(u, v))
    rev_pts = np.array(rev_pts)
    pf = parent_faces(case)
    d2 = project_distance(rev_pts, pf)
    k2 = int(np.argmax(d2))
    reverse = {"samples": len(rev_pts), "max_mm": float(d2.max()), "at_xyz": rev_pts[k2],
               "what": "points on the model's retained faces to the parent.step surfaces (attachment end faces excluded)"}
    worst = max(forward["max_mm"], reverse["max_mm"])
    return {"parent_to_model": forward, "model_to_parent": reverse, "max_mm": worst,
            "pass": worst <= TOL["retained_skin_mm"]}


def outline_surface_distance(model, case):
    """3D locked outlines only: distance from dense outline samples to the cap faces."""
    if case.outline_is_plan:
        return None
    out = []
    for tip in case.tips:
        p = np.linspace(0, 1, 401)
        Q = case.outline(p, tip)
        d = project_distance(Q, model.cap_faces(tip))
        k = int(np.argmax(d))
        out.append({"tip": tip.name, "samples": len(p), "max_mm": float(d.max()), "at_p": float(p[k]), "at_xyz": Q[k]})
    return out


# ---------------------------------------------------------------------------
# 3. Joins
# ---------------------------------------------------------------------------

def edge_faces(model):
    m = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(model.shape, TopAbs_EDGE, TopAbs_FACE, m)
    pairs = []
    for i in range(1, m.Extent() + 1):
        edge = TopoDS.Edge_s(m.FindKey(i))
        idx = []
        for face in m.FindFromIndex(i):
            j = model.face_map.FindIndex(face) - 1
            if j not in idx:
                idx.append(j)
        pairs.append((edge, idx))
    return pairs


def join_kind(a, b):
    capa, capb = a.role != "retained", b.role != "retained"
    if capa and capb:
        return "cap seam with a planar cap face" if (a.planar or b.planar) else "cap seam"
    if not (capa or capb):
        return None
    parent = a if not capa else b
    return "parent planar face to cap" if parent.planar else "parent skin to cap"


JOIN_FIELDS = ("gap_mm", "angle_deg", "curvature_mismatch_per_mm", "cross_join_normal_curvature_mismatch_per_mm")


def sample_join(model, case, edge, fa, fb, n=401):
    """Dense samples along one shared edge: gap, tangent-plane angle and the second fundamental forms of both
    faces in one tangent frame (e1 along the edge, e2 = N x e1). Samples within COLLAPSE_EXCLUSION_MM of a
    collapsed corner are kept in the 'all' maxima but left out of the judged maxima."""
    if BRep_Tool.Degenerated_s(edge):
        return None
    curve = BRepAdaptor_Curve(edge)
    t0, t1 = curve.FirstParameter(), curve.LastParameter()
    pa, pb = BRepAdaptor_Curve2d(edge, fa.face), BRepAdaptor_Curve2d(edge, fb.face)
    k = np.arange(1, n - 1)
    s = 0.5 * (1 - np.cos(np.pi * k / (n - 1)))
    rows, skipped = [], 0
    end_a, end_b = np.array(curve.Value(t0).Coord()), np.array(curve.Value(t1).Coord())
    tip = fa.tip or fb.tip
    for x in s:
        t = t0 + (t1 - t0) * x
        p3, v3 = gp_Pnt(), gp_Vec()
        curve.D1(t, p3, v3)
        ja = fa.jet(*pa.Value(t).Coord())
        jb = fb.jet(*pb.Value(t).Coord())
        if ja["N"] is None or jb["N"] is None or v3.Magnitude() == 0:
            skipped += 1
            continue
        e1 = unit(np.array(v3.Coord()))
        e2 = unit(np.cross(ja["N"], e1))
        IIa = np.array([[e1 @ ja["H"] @ e1, e1 @ ja["H"] @ e2], [e2 @ ja["H"] @ e1, e2 @ ja["H"] @ e2]])
        IIb = np.array([[e1 @ jb["H"] @ e1, e1 @ jb["H"] @ e2], [e2 @ jb["H"] @ e1, e2 @ jb["H"] @ e2]])
        ang = math.degrees(math.atan2(np.linalg.norm(np.cross(ja["N"], jb["N"])), np.clip(ja["N"] @ jb["N"], -1, 1)))
        rows.append((x, ja["P"], np.linalg.norm(ja["P"] - jb["P"]), ang, np.linalg.norm(IIa - IIb),
                     abs(IIa[1, 1] - IIb[1, 1]), IIa, IIb))
    if not rows:
        return {"samples": 0, "skipped_degenerate": skipped}
    P = np.array([r[1] for r in rows])
    near = np.minimum(fa.distance_to_collapse(P), fb.distance_to_collapse(P)) < COLLAPSE_EXCLUSION_MM
    cols = {name: np.array([r[i] for r in rows]) for name, i in zip(JOIN_FIELDS, (2, 3, 4, 5))}
    frac = tip.fraction_of(P) if tip is not None else np.full(len(P), np.nan)
    res = {"samples": len(rows), "skipped_degenerate": skipped, "samples_near_collapse": int(near.sum()),
           "closest_sample_to_ends_mm": [float(np.linalg.norm(P[0] - end_a)), float(np.linalg.norm(P[-1] - end_b))],
           "chord_fraction_range": [float(np.nanmin(frac)), float(np.nanmax(frac))] if np.isfinite(frac).any() else None,
           "length_mm": float(sum(np.linalg.norm(np.diff(P, axis=0), axis=1)))}
    judged = ~near if (~near).any() else np.ones(len(P), bool)
    for name, v in cols.items():
        j = int(np.argmax(np.where(judged, v, -np.inf)))
        ja_ = int(np.argmax(v))
        res[name] = {"max": float(v[j]), "at_xyz": P[j], "at_chord_fraction": frac[j], "median": float(np.median(v[judged])),
                     "max_including_near_collapse": float(v[ja_]), "at_xyz_including_near_collapse": P[ja_]}
    j = int(np.argmax(np.where(judged, cols["curvature_mismatch_per_mm"], -np.inf)))
    res["curvature_mismatch_per_mm"]["II_parent_side_at_max"] = rows[j][6]
    res["curvature_mismatch_per_mm"]["II_cap_side_at_max"] = rows[j][7]
    return res


def joins(model, case):
    out = []
    for edge, idx in edge_faces(model):
        if len(idx) != 2:
            if len(idx) == 1 and model.faces[idx[0]].role != "retained" and not BRep_Tool.Degenerated_s(edge):
                out.append({"kind": "cap face closing seam (not sampled)", "faces": idx})
            continue
        fa, fb = model.faces[idx[0]], model.faces[idx[1]]
        kind = join_kind(fa, fb)
        if kind is None:
            continue
        if fb.role == "retained":
            fa, fb = fb, fa
        r = sample_join(model, case, edge, fa, fb)
        if r is None:
            continue
        r.update(kind=kind, faces=[fa.index, fb.index], tip=(fa.tip or fb.tip).name if (fa.tip or fb.tip) else None)
        out.append(r)
    summary = {}
    for kind in ("parent skin to cap", "parent planar face to cap", "cap seam", "cap seam with a planar cap face"):
        rows = [r for r in out if r.get("kind") == kind and r.get("samples")]
        if not rows:
            continue
        summary[kind] = {"edges": len(rows), "samples": sum(r["samples"] for r in rows),
                         "skipped_degenerate": sum(r["skipped_degenerate"] for r in rows),
                         "samples_near_collapse": sum(r["samples_near_collapse"] for r in rows),
                         "closest_sample_to_edge_ends_mm": min(min(r["closest_sample_to_ends_mm"]) for r in rows)}
        for name in JOIN_FIELDS:
            r = max(rows, key=lambda r: r[name]["max"])
            summary[kind][name] = {"max": r[name]["max"], "at_xyz": r[name]["at_xyz"], "at_chord_fraction": r[name]["at_chord_fraction"],
                                   "max_including_near_collapse": max(q[name]["max_including_near_collapse"] for q in rows)}
    skin = summary.get("parent skin to cap")
    ok = bool(skin) and skin["gap_mm"]["max"] <= TOL["join_gap_mm"] and skin["angle_deg"]["max"] <= TOL["join_angle_deg"]
    seam = summary.get("cap seam")
    seam_ok = (not seam) or (seam["gap_mm"]["max"] <= TOL["join_gap_mm"] and seam["angle_deg"]["max"] <= TOL["join_angle_deg"])
    return {"edges": out, "summary": summary, "pass_G1": ok, "cap_seams_pass_G1": seam_ok}


# ---------------------------------------------------------------------------
# 4. Section traces
# ---------------------------------------------------------------------------

class Segment:
    """One section edge on one cap face, sampled adaptively on the exact plane."""

    def __init__(self, fr, edge, O, n):
        self.fr, self.edge, self.O, self.n = fr, edge, O, n
        try:
            self.pc = BRepAdaptor_Curve2d(edge, fr.face)
            self.t0, self.t1 = self.pc.FirstParameter(), self.pc.LastParameter()
            self.c3 = None
        except Exception:
            self.pc = None
            self.c3 = BRepAdaptor_Curve(edge)
            self.t0, self.t1 = self.c3.FirstParameter(), self.c3.LastParameter()
        self.cache = {}

    def uv(self, t):
        if self.pc is not None:
            u, v = self.pc.Value(t).Coord()
        else:
            sas = self.fr.analysis()
            u, v = sas.ValueOfUV(self.c3.Value(t), 1e-10).Coord()
        u0, u1, v0, v1 = self.fr.uv
        for _ in range(6):   # Newton onto the exact cut plane (the approximated pcurve is only close)
            P, Su, Sv = self.fr.d2(u, v)[:3]
            d = (P - self.O) @ self.n
            g = np.array([self.n @ Su, self.n @ Sv])
            gg = g @ g
            if abs(d) < 1e-13 or gg == 0:
                break
            u, v = np.clip(u - d * g[0] / gg, u0, u1), np.clip(v - d * g[1] / gg, v0, v1)
        return u, v

    def at(self, t):
        if t in self.cache:
            return self.cache[t]
        j = self.fr.jet(*self.uv(t))
        if j["N"] is None:
            r = dict(t=t, P=j["P"], kappa=np.nan, tau=None, ok=False)
        else:
            cr = np.cross(j["N"], self.n)
            s = np.linalg.norm(cr)
            tau = cr / s if s > 0 else None
            kn = -(tau @ j["H"] @ tau) if tau is not None else np.nan
            kappa = kn / s if tau is not None else np.nan
            kp = j["k"][np.argmax(np.abs(j["k"]))]
            r = dict(t=t, P=j["P"], kappa=kappa, kn=kn, kp=kp, tau=tau, ok=tau is not None, N=j["N"], H=j["H"])
        self.cache[t] = r
        return r

    def sample(self, max_samples=8000):
        ts = list(np.linspace(self.t0, self.t1, 65))
        span = self.t1 - self.t0
        for _ in range(60):
            rows = [self.at(t) for t in ts]
            new = []
            for a, b in zip(rows, rows[1:]):
                if b["t"] - a["t"] <= 1e-13 * max(span, 1e-300):
                    continue
                turn = 180.0 if (a["tau"] is None or b["tau"] is None) else math.degrees(
                    math.acos(np.clip(abs(a["tau"] @ b["tau"]), -1, 1)))
                if turn > TURN_STEP_DEG or np.linalg.norm(b["P"] - a["P"]) > 0.25:
                    new.append(0.5 * (a["t"] + b["t"]))
            if not new or len(ts) + len(new) > max_samples:
                break
            ts = sorted(set(ts) | set(new))
        self.rows = [self.at(t) for t in ts]
        return self.rows

    def refine_max(self, i):
        lo = self.rows[max(0, i - 1)]["t"]
        hi = self.rows[min(len(self.rows) - 1, i + 1)]["t"]
        if hi <= lo:
            return self.rows[i]
        f = lambda t: -abs(self.at(float(t))["kappa"]) if self.at(float(t))["ok"] else 0.0
        r = minimize_scalar(f, bounds=(lo, hi), method="bounded", options={"xatol": 1e-14 * max(1.0, abs(hi))})
        best = self.at(float(r.x))
        return best if best["ok"] and abs(best["kappa"]) >= abs(self.rows[i]["kappa"]) else self.rows[i]

    def refine_extreme(self, i, w):
        lo = self.rows[max(0, i - 1)]["t"]
        hi = self.rows[min(len(self.rows) - 1, i + 1)]["t"]
        if hi <= lo:
            return self.rows[i]
        r = minimize_scalar(lambda t: -(self.at(float(t))["P"] @ w), bounds=(lo, hi), method="bounded",
                            options={"xatol": 1e-14 * max(1.0, abs(hi))})
        best = self.at(float(r.x))
        return best if best["P"] @ w >= self.rows[i]["P"] @ w else self.rows[i]

    def root(self, ta, tb):
        g = lambda t: self.at(float(t))["kappa"]
        try:
            return self.at(float(brentq(g, ta, tb, xtol=1e-15 * max(1.0, abs(tb)))))
        except ValueError:
            return self.at(0.5 * (ta + tb))


def section_segments(model, tip, O, n):
    pln = gp_Pln(gp_Pnt(*map(float, O)), gp_Dir(*map(float, n)))
    segs = []
    for fr in model.cap_faces(tip):
        sec = BRepAlgoAPI_Section(fr.face, pln, False)
        sec.ComputePCurveOn1(True)
        sec.Approximation(True)
        sec.Build()
        if not sec.IsDone():
            raise RuntimeError(f"section failed on face {fr.index}")
        edges, _ = shapes_of(sec.Shape(), TopAbs_EDGE)
        for e in edges:
            segs.append(Segment(fr, TopoDS.Edge_s(e), O, n))
    return segs


def chain(segs, tip, T, tol=1e-3):
    """Order the sampled segments from the upper attachment, round the lip, to the lower attachment."""
    ends = []
    for k, s in enumerate(segs):
        ends.append((k, 0, s.rows[0]["P"]))
        ends.append((k, 1, s.rows[-1]["P"]))
    att = [e for e in ends if abs(e[2] @ tip.outboard - tip.station) < 1e-3]
    if not att:
        return None, "no trace end lies on the attachment"
    start = max(att, key=lambda e: e[2] @ T)
    order, used, cur, gaps = [], set(), start, []
    while True:
        k, side = cur[0], cur[1]
        used.add(k)
        order.append((k, side == 1))   # reversed if we enter at its end
        exit_p = segs[k].rows[0 if side == 1 else -1]["P"]
        cand = [e for e in ends if e[0] not in used]
        if not cand:
            break
        nxt = min(cand, key=lambda e: np.linalg.norm(e[2] - exit_p))
        d = np.linalg.norm(nxt[2] - exit_p)
        if d > tol:
            break
        gaps.append(d)
        cur = nxt
    return order, {"pieces_unused": len(segs) - len(used), "max_junction_gap_mm": max(gaps) if gaps else 0.0}


def trace(model, case, tip, f, label, light=False):
    O, n, T, w = tip.plane(f)
    t_att = float(case.thickness_exact(f)[0]) if np.isfinite(case.thickness_exact(f)[0]) else float(case.thickness(f))
    res = {"fraction": f, "label": label, "attachment_thickness_mm": t_att,
           "zone": "thick" if t_att >= RULE["thick_gate_mm"] else "thin", "plane_point": O, "plane_normal": n}
    try:
        segs = section_segments(model, tip, O, n)
        if not segs:
            res["error"] = "plane misses the cap"
            return res
        for s in segs:
            s.sample(max_samples=600 if light else 8000)
        order, info = chain(segs, tip, T)
        if order is None:
            res["error"] = info
            return res
        res.update(info)
    except Exception as exc:   # noqa: BLE001 - record and move on; the summary marks the cut as failed
        res["error"] = f"{type(exc).__name__}: {exc}"
        return res
    rows, owner = [], []
    for k, rev in order:
        r = segs[k].rows[::-1] if rev else segs[k].rows
        for i, row in enumerate(r):
            rows.append(row)
            owner.append((k, (len(r) - 1 - i) if rev else i))
    P = np.array([r["P"] for r in rows])
    kap = np.array([r["kappa"] for r in rows])
    arc = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
    length = arc[-1]
    # Outermost point along the in-plane outboard direction (for a 3D rim, the direction keeps the rim's rise).
    wl = tip.span_direction(f)
    i = int(np.argmax(P @ wl))
    k, j = owner[i]
    outer = segs[k].refine_extreme(j, wl)
    outer_arc = arc[i] / length if length else np.nan
    q_p, Q = case.outline_in_plane(tip, O, n, near=outer["P"])
    kind = case.outline_kind
    # The lip: where the locked perimeter crosses the cut (rim and mid-plane outlines); for a plan outline the
    # perimeter fixes only the plan position, so the lip is the outermost point in plan.
    lip, lip_arc, trace_distance = outer, outer_arc, None
    if Q is not None and kind != "plan":
        ii = int(np.argmin(np.linalg.norm(P - Q, axis=1)))
        kk, jj = owner[ii]
        seg = segs[kk]
        lo_, hi_ = seg.rows[max(0, jj - 1)]["t"], seg.rows[min(len(seg.rows) - 1, jj + 1)]["t"]
        lip = rows[ii]
        if hi_ > lo_:
            rr = minimize_scalar(lambda t: np.linalg.norm(seg.at(float(t))["P"] - Q), bounds=(lo_, hi_), method="bounded",
                                 options={"xatol": 1e-14 * max(1.0, abs(hi_))})
            cand = seg.at(float(rr.x))
            if np.linalg.norm(cand["P"] - Q) < np.linalg.norm(lip["P"] - Q):
                lip = cand
        trace_distance = float(np.linalg.norm(lip["P"] - Q))
        lip_arc = arc[ii] / length if length else np.nan

    def radius(row, key="kappa"):
        return 1 / abs(row[key]) if row.get("ok") and row[key] else None

    res.update(length_mm=length, samples=len(rows), plane_error_mm=float(np.max(np.abs((P - O) @ n))),
               trace_ends_thickness_mm=float(abs((P[0] - P[-1]) @ T)),
               lip={"defined_as": "outermost point in plan" if kind == "plan" or Q is None else "trace point at the locked outline",
                    "xyz": lip["P"], "radius_mm": radius(lip), "principal_radius_mm": radius(lip, "kp"),
                    "signed_curvature_per_mm": lip["kappa"], "arclength_fraction": lip_arc},
               outermost={"xyz": outer["P"], "radius_mm": radius(outer), "arclength_fraction": outer_arc})
    if Q is not None:
        if kind == "plan":
            dev = float((np.nan_to_num(outer["P"] - Q) * np.array([1, 0, 1])) @ wl)
        else:
            dev = float((outer["P"] - Q) @ wl)
        res["outline"] = {"kind": kind, "p": q_p, "xyz": Q, "outermost_offset_mm": dev}
        if trace_distance is not None:
            res["outline"]["trace_distance_mm"] = trace_distance
        # Judged deviation: plan outline = outermost point in plan; mid-plane outline = on the cap and outermost;
        # cambered 3D rim = on the cap (its outermost point sits off the rim by construction).
        res["outline"]["deviation_mm"] = {"plan": abs(dev), "midplane": max(abs(dev), trace_distance or 0.0),
                                          "rim": trace_distance}[kind]
        tan = np.nan_to_num(case.outline([min(q_p + 1e-6, 1)], tip)[0] - case.outline([max(q_p - 1e-6, 0)], tip)[0])
        if np.linalg.norm(tan) > 0:
            ell = unit(tan)
            res["lip"]["cut_obliquity_deg"] = math.degrees(math.acos(min(1.0, abs(ell @ n))))
            if lip.get("ok"):
                Nl = lip["N"]
                lt = ell - (ell @ Nl) * Nl
                if np.linalg.norm(lt) > 0:
                    b = unit(np.cross(Nl, unit(lt)))
                    kperp = -(b @ lip["H"] @ b)
                    res["lip"]["lip_normal_radius_mm"] = 1 / abs(kperp) if kperp else np.inf
    if light:
        return res
    # Smallest radius with local refinement of the sampled peak.
    finite = np.isfinite(kap)
    i = int(np.nanargmax(np.where(finite, np.abs(kap), -1)))
    k, j = owner[i]
    best = segs[k].refine_max(j)
    for seg in segs:   # a refined peak on another branch may beat the sampled one
        pass
    r_min = 1 / abs(best["kappa"]) if best["kappa"] else np.inf
    res["min_radius"] = {"radius_mm": r_min, "xyz": best["P"], "signed_curvature_per_mm": best["kappa"],
                         "normal_section_radius_mm": 1 / abs(best["kn"]) if best["kn"] else np.inf,
                         "principal_radius_mm": 1 / abs(best["kp"]) if best["kp"] else np.inf,
                         "arclength_fraction": arc[i] / length if length else np.nan,
                         "distance_from_lip_mm": float(np.linalg.norm(best["P"] - lip["P"])),
                         "convex": bool(best["kappa"] > 0)}
    # Bending reversals between curvature samples outside the deadband, refined by root finding.
    nz = np.flatnonzero(finite & (np.abs(kap) > CURVATURE_DEADBAND))
    revs = []
    for a, b in zip(nz, nz[1:]):
        if kap[a] * kap[b] < 0:
            (ka, ja), (kb, jb) = owner[a], owner[b]
            if ka == kb:
                ta, tb = sorted((segs[ka].rows[ja]["t"], segs[kb].rows[jb]["t"]))
                r = segs[ka].root(ta, tb)
                revs.append({"xyz": r["P"], "arclength_fraction": arc[a] / length, "on": "face interior"})
            else:
                revs.append({"xyz": P[b], "arclength_fraction": arc[b] / length, "on": "face boundary"})
    res["reversals"] = revs
    # Tangent breaks where the trace crosses from one cap face to the next (a crease is a zero radius).
    turns = []
    for a in range(len(rows) - 1):
        if owner[a][0] != owner[a + 1][0] and rows[a]["tau"] is not None and rows[a + 1]["tau"] is not None:
            turns.append(math.degrees(math.acos(min(1.0, abs(rows[a]["tau"] @ rows[a + 1]["tau"])))))
    res["max_face_junction_turn_deg"] = max(turns) if turns else 0.0
    if res["max_face_junction_turn_deg"] > TOL["join_angle_deg"]:
        r_min = 0.0
        res["min_radius"]["radius_mm"] = 0.0
        res["min_radius"]["crease"] = True
    # Peak magnitudes on either side of the lip (a tight bend before the lip shows up here).
    side_peaks = {}
    for name, m in (("upper", np.arange(len(P)) < np.argmax(P @ wl)), ("lower", np.arange(len(P)) > np.argmax(P @ wl))):
        if m.any() and finite[m].any():
            q = np.flatnonzero(m & finite)
            ii = q[np.argmax(np.abs(kap[q]))]
            side_peaks[name] = {"radius_mm": 1 / abs(kap[ii]), "arclength_fraction": arc[ii] / length}
    res["side_peaks"] = side_peaks
    req = float(case.required_radius(t_att))
    res["required_radius_mm"] = req
    if res["zone"] == "thick":
        res["margin_mm"] = r_min - RULE["min_radius_mm"]
        res["pass"] = bool(r_min >= RULE["min_radius_mm"] * (1 - VERDICT_RTOL))
    else:
        res["ratio_to_half_thickness"] = r_min / (0.5 * t_att) if t_att > 0 else None
        res["lip_ratio_to_half_thickness"] = (res["lip"]["radius_mm"] / (0.5 * t_att)) if t_att > 0 and res["lip"]["radius_mm"] else None
        res["pass"] = bool(t_att > 0 and r_min / (0.5 * t_att) >= RULE["thin_ratio_pass"] * (1 - VERDICT_RTOL))
    return res


def section_fractions(case):
    rows = [(f, "listed") for f in SECTION_FRACTIONS]
    rows += [(f, "0.8 mm thickness crossing") for f in case.crossings]
    for t in THIN_TARGETS_TE[case.variant]:
        f = case.fraction_at_thickness(t, "te")
        if f is not None:
            rows.append((f, f"thin zone, t = {t} mm (TE side)"))
    for t in THIN_TARGETS_LE:
        f = case.fraction_at_thickness(t, "le")
        if f is not None:
            rows.append((f, f"thin zone, t = {t} mm (LE side)"))
    return sorted(rows)


def sections(model, case):
    out = {}
    for tip in case.tips:
        out[tip.name] = [trace(model, case, tip, f, label) for f, label in section_fractions(case)]
    return out


def outline_lip_scan(model, case, n=97):
    """Outermost lip (or, for a 3D rim, the cap surface) against the locked outline on a dense set of cut planes."""
    fr = np.concatenate([[0.005, 0.01], np.linspace(0.02, 0.98, n), [0.99, 0.995]])
    out = {}
    for tip in case.tips:
        rows, errors = [], 0
        for f in fr:
            r = trace(model, case, tip, float(f), "outline scan", light=True)
            if "outline" in r:
                o = r["outline"]
                rows.append((f, o["deviation_mm"], o["outermost_offset_mm"], o.get("trace_distance_mm", np.nan)))
            else:
                errors += 1
        if not rows:
            out[tip.name] = {"planes": 0, "errors": errors}
            continue
        a = np.array(rows, float)
        k = int(np.nanargmax(a[:, 1]))
        out[tip.name] = {"kind": case.outline_kind, "planes": len(rows), "errors": errors,
                         "max_deviation_mm": float(a[k, 1]), "at_fraction": float(a[k, 0]),
                         "max_abs_lip_outboard_offset_mm": float(np.nanmax(np.abs(a[:, 2]))),
                         "max_outline_to_trace_mm": float(np.nanmax(a[:, 3])) if np.isfinite(a[:, 3]).any() else None}
    return out


# ---------------------------------------------------------------------------
# 5. Whole-cap principal-curvature map
# ---------------------------------------------------------------------------

def collapsed_points(fr):
    u0, u1, v0, v1 = fr.uv
    out = []
    for name, pts in (("u=min", [fr.point(u0, v) for v in np.linspace(v0, v1, 21)]),
                      ("u=max", [fr.point(u1, v) for v in np.linspace(v0, v1, 21)]),
                      ("v=min", [fr.point(u, v0) for u in np.linspace(u0, u1, 21)]),
                      ("v=max", [fr.point(u, v1) for u in np.linspace(u0, u1, 21)])):
        pts = np.array(pts)
        if np.max(np.linalg.norm(pts - pts.mean(0), axis=1)) < 1e-6:
            out.append((f"collapsed side {name}", pts.mean(0)))
    for e, _ in [(e, None) for e in shapes_of(fr.face, TopAbs_EDGE)[0]]:
        e = TopoDS.Edge_s(e)
        if BRep_Tool.Degenerated_s(e):
            v = shapes_of(e, TopAbs_VERTEX)[0][0]
            out.append(("degenerated edge", np.array(BRep_Tool.Pnt_s(TopoDS.Vertex_s(v)).Coord())))
    return out


def cap_map(model, case, n_uniform=241, n_cos=81, n_refine=12):
    """Larger principal curvature over every cap face on a dense (u, v) grid, judged point by point against the
    rule at the chord fraction of the family plane through the point; local extrema refined in (u, v)."""
    faces_out, raw = [], []
    excl = COLLAPSE_EXCLUSION_MM
    for fr in model.cap_faces():
        tip = fr.tip
        us, vs = fr.grid(n_uniform, n_cos)
        P = np.full((len(us), len(vs), 3), np.nan)
        K = np.full((len(us), len(vs)), np.nan)      # larger-magnitude principal curvature, signed (convex +)
        excluded, outside = 0, 0
        for a, u in enumerate(us):
            for b, v in enumerate(vs):
                if not fr.inside(u, v):
                    outside += 1
                    continue
                j = fr.jet(u, v)
                P[a, b] = j["P"]
                if j["N"] is None:
                    excluded += 1
                    continue
                K[a, b] = j["k"][np.argmax(np.abs(j["k"]))]
        valid = np.isfinite(K)
        flatP = P.reshape(-1, 3)
        have = np.isfinite(flatP).all(1)
        frac = np.full(len(flatP), np.nan)
        frac[have] = tip.fraction_of(flatP[have])
        frac = frac.reshape(K.shape)
        dcol = np.full(len(flatP), np.inf)
        dcol[have] = fr.distance_to_collapse(flatP[have])
        dcol = dcol.reshape(K.shape)
        t = case.thickness(frac)
        req = case.required_radius(np.nan_to_num(t))
        with np.errstate(divide="ignore", invalid="ignore"):
            R = 1 / np.abs(K)
            score = np.where(valid & np.isfinite(t) & (req > 0), R / req, np.nan)
        judged = valid & (dcol >= excl)
        u0, u1, v0, v1 = fr.uv
        su, sv = (u1 - u0), (v1 - v0)

        def evaluate(u, v, honour_exclusion=True):
            if not (u0 <= u <= u1 and v0 <= v <= v1) or not fr.inside(u, v):
                return None
            j = fr.jet(u, v)
            if j["N"] is None:
                return None
            dc = float(fr.distance_to_collapse(j["P"])[0])
            if honour_exclusion and dc < excl:
                return None
            k = j["k"][np.argmax(np.abs(j["k"]))]
            f = float(tip.fraction_of(j["P"][None])[0])
            tt = float(case.thickness(f))
            rq = float(case.required_radius(tt))
            return dict(u=u, v=v, P=j["P"], k=k, R=1 / abs(k) if k else np.inf, f=f, t=tt, req=rq, dc=dc,
                        score=(1 / abs(k)) / rq if rq > 0 and k else np.inf)

        def at_grid(a, b):
            return dict(u=us[a], v=vs[b], P=P[a, b], k=K[a, b], R=R[a, b], f=frac[a, b], t=t[a, b], req=req[a, b],
                        score=score[a, b], dc=dcol[a, b])

        # Seeds: lowest local minima of the score and highest local maxima of |K| among judged grid points.
        seeds = set()
        for field, sense in ((score, 1), (np.abs(K), -1)):
            Fm = np.where(judged & np.isfinite(field), sense * field, np.inf)
            pad = np.pad(Fm, 1, constant_values=np.inf)
            local = np.ones_like(Fm, bool)
            for da in (-1, 0, 1):
                for db in (-1, 0, 1):
                    if da or db:
                        local &= Fm <= pad[1 + da:1 + da + Fm.shape[0], 1 + db:1 + db + Fm.shape[1]]
            cand = np.argwhere(local & np.isfinite(Fm))
            seeds.update(map(tuple, sorted(cand.tolist(), key=lambda ab: Fm[ab[0], ab[1]])[:n_refine]))
        refined = []
        for a, b in seeds:
            x0 = np.array([(us[a] - u0) / su, (vs[b] - v0) / sv])
            simplex = [x0, x0 + [1e-3, 0], x0 + [0, 1e-3]]
            for key, sign in (("score", 1), ("k", -1)):
                def obj(x):
                    r = evaluate(u0 + x[0] * su, v0 + x[1] * sv)
                    if r is None:
                        return 1e30
                    return r["score"] if key == "score" else -abs(r["k"])
                try:
                    r = minimize(obj, x0, method="Nelder-Mead",
                                 options={"xatol": 1e-11, "fatol": 1e-14, "maxiter": 400, "initial_simplex": simplex})
                    best = evaluate(u0 + r.x[0] * su, v0 + r.x[1] * sv)
                    if best is not None:
                        refined.append(best)
                except Exception:   # noqa: BLE001
                    pass
        cands = refined + [at_grid(a, b) for a, b in seeds]
        if not cands and judged.any():
            gi = np.nanargmin(np.where(judged, score, np.inf))
            cands.append(at_grid(*np.unravel_index(gi, score.shape)))
        worst = min(cands, key=lambda r: r["score"]) if cands else None
        tightest = max(cands, key=lambda r: abs(r["k"])) if cands else None
        zone_worst = {}
        for name in ("thick", "thin"):
            zc = [c for c in cands if np.isfinite(c["t"]) and (c["t"] >= RULE["thick_gate_mm"]) == (name == "thick")]
            mask = judged & np.isfinite(t) & ((np.nan_to_num(t) >= RULE["thick_gate_mm"]) == (name == "thick"))
            if mask.any():
                zc.append(at_grid(*np.unravel_index(np.nanargmin(np.where(mask, score, np.inf)), score.shape)))
            if zc:
                zone_worst[name] = describe_point(min(zc, key=lambda r: r["score"]))
        # Near the collapsed corners: what the grid shows inside the exclusion radius.
        bands = []
        for lo, hi in ((0, 1e-3), (1e-3, 1e-2), (1e-2, excl)):
            m = valid & (dcol >= lo) & (dcol < hi)
            if m.any():
                gi = np.nanargmin(np.where(m, score, np.inf))
                bands.append({"distance_mm": [lo, hi], "points": int(m.sum()), "min_radius_mm": float(np.nanmin(np.where(m, R, np.inf))),
                              "worst": describe_point(at_grid(*np.unravel_index(gi, score.shape)))})
        cps = fr.collapsed()
        vp = P[valid]
        rec = {"face": fr.index, "tip": tip.name, "surface": fr.type, "grid": [len(us), len(vs)],
               "points_evaluated": int(np.isfinite(P).all(-1).sum()), "points_outside_trim": outside,
               "points_excluded_undefined_derivatives": excluded,
               "points_within_collapse_exclusion": int((valid & (dcol < excl)).sum()),
               "points_outside_chord_range": int(np.sum(valid & ~((frac >= case._tf_range[0]) & (frac <= case._tf_range[1])))),
               "collapsed_points": [{"what": w_, "xyz": c,
                                     "nearest_grid_point_with_derivatives_mm": float(np.min(np.linalg.norm(vp - c, axis=1))) if len(vp) else None}
                                    for w_, c in cps],
               "near_collapse": bands,
               "worst": describe_point(worst) if worst else None, "tightest": describe_point(tightest) if tightest else None,
               "zone_worst": zone_worst}
        faces_out.append(rec)
        raw.append(dict(face=fr, P=P[valid], R=R[valid], f=frac[valid], t=t[valid], score=score[valid], dc=dcol[valid]))
    model.map_raw = raw
    judged_faces = [r for r in faces_out if r["worst"]]
    if not judged_faces:
        return {"faces": faces_out, "pass": False}
    overall = min(judged_faces, key=lambda r: r["worst"]["score_radius_over_required"])
    tight = max(judged_faces, key=lambda r: abs(r["tightest"]["signed_curvature_per_mm"]))
    near = [b for r in faces_out for b in r["near_collapse"]]
    return {"collapse_exclusion_mm": excl, "faces": faces_out, "worst": overall["worst"], "worst_face": overall["face"],
            "tightest": tight["tightest"], "tightest_face": tight["face"],
            "worst_within_collapse_exclusion": min((b["worst"] for b in near), key=lambda w: w["score_radius_over_required"], default=None),
            "closest_grid_approach_to_collapse_mm": min((c["nearest_grid_point_with_derivatives_mm"] for r in faces_out
                                                        for c in r["collapsed_points"] if c["nearest_grid_point_with_derivatives_mm"] is not None), default=None),
            "pass": bool(overall["worst"]["score_radius_over_required"] >= 1.0 - VERDICT_RTOL)}


def describe_point(r):
    t = float(r["t"]) if r["t"] is not None and np.isfinite(r["t"]) else None
    d = {"xyz": r["P"], "uv": [r["u"], r["v"]], "radius_mm": float(r["R"]), "signed_curvature_per_mm": float(r["k"]),
         "chord_fraction": float(r["f"]) if np.isfinite(r["f"]) else None, "attachment_thickness_mm": t,
         "required_radius_mm": float(r["req"]), "score_radius_over_required": float(r["score"]),
         "distance_to_collapsed_corner_mm": float(r.get("dc", np.inf))}
    if t is not None and t >= RULE["thick_gate_mm"]:
        d["zone"], d["margin_mm"] = "thick", float(r["R"]) - RULE["min_radius_mm"]
    elif t is not None:
        d["zone"], d["ratio_to_half_thickness"] = "thin", float(r["R"]) / (0.5 * t) if t > 0 else None
    return d


# ---------------------------------------------------------------------------
# 6. Trailing-edge closure
# ---------------------------------------------------------------------------

def te_closure(model, case, joins_out):
    """How the cap meets the trailing edge, and the radius there."""
    m = case.parent
    out = {"kind": "blunt face end" if case.blunt else "sharp point", "tips": []}
    raw = getattr(model, "map_raw", [])
    for tip in case.tips:
        rec = {"tip": tip.name}
        approach = []
        for d in TE_APPROACH_MM:
            f = case.f_te - d / tip.chord
            r = trace(model, case, tip, f, f"{d} mm before the TE")
            approach.append({"distance_to_te_mm": d, "fraction": f, "attachment_thickness_mm": r["attachment_thickness_mm"],
                             "min_radius_mm": r.get("min_radius", {}).get("radius_mm"),
                             "lip_radius_mm": (r.get("lip") or {}).get("radius_mm"),
                             "lip_normal_radius_mm": (r.get("lip") or {}).get("lip_normal_radius_mm"),
                             "cut_obliquity_deg": (r.get("lip") or {}).get("cut_obliquity_deg"),
                             "ratio_to_half_thickness": r.get("ratio_to_half_thickness"),
                             "pass": r.get("pass"), "error": r.get("error")})
        rec["approach_sections"] = approach
        if case.blunt:
            seg = np.asarray(m.te_face(np.array([0.0, 1.0]), np.array([tip.span, tip.span])), float).reshape(2, 3)
            rec["parent_te_segment"] = seg
            te_joins = [e for e in joins_out["edges"] if e.get("kind") == "parent planar face to cap" and e.get("tip") == tip.name]
            rec["te_face_joins"] = [{k: e[k] for k in ("faces", "length_mm", "gap_mm", "angle_deg", "samples", "skipped_degenerate")}
                                    for e in te_joins if e.get("samples")]
            planar = [e for e in joins_out["edges"] if e.get("kind") == "cap seam with a planar cap face" and e.get("tip") == tip.name and e.get("samples")]
            rec["seams_with_planar_cap_faces"] = [{k: e[k] for k in ("faces", "length_mm", "angle_deg")} for e in planar]
            a, b = seg
            def seg_dist(P):
                ab = b - a
                s_ = np.clip(((P - a) @ ab) / (ab @ ab), 0, 1)
                return np.linalg.norm(P - (a + s_[:, None] * ab), axis=1)
            near = []
            for r_ in raw:
                if r_["face"].tip is tip and not r_["face"].planar and len(r_["P"]):
                    msk = seg_dist(r_["P"]) <= 0.1
                    near += list(r_["R"][msk])
            if near:
                rec["curved_cap_principal_radius_within_0p1mm_of_te_segment"] = {
                    "points": len(near), "min_mm": float(np.min(near)), "median_mm": float(np.median(np.minimum(near, 1e6))),
                    "half_te_thickness_mm": 0.5 * m.TE_THICKNESS}
            rec["closes"] = bool(te_joins) and all(e["gap_mm"]["max"] <= TOL["join_gap_mm"] for e in te_joins if e.get("samples"))
        else:
            if case.shape == "wingb":
                te = np.asarray(m.skin(np.array([1.0]), np.array([tip.span])), float).reshape(3)
            else:
                te = np.asarray(m.skin(1.0, tip.span, +1), float).reshape(3)
            verts = []
            for fr in model.cap_faces(tip):
                for v in shapes_of(fr.face, TopAbs_VERTEX)[0]:
                    verts.append(np.array(BRep_Tool.Pnt_s(TopoDS.Vertex_s(v)).Coord()))
                verts += [c for _, c in fr.collapsed()]
            dist = min(np.linalg.norm(v - te) for v in verts) if verts else np.inf
            rec["parent_te_point"] = te
            rec["nearest_cap_vertex_or_collapse_mm"] = float(dist)
            near = []
            for r_ in raw:
                if r_["face"].tip is tip and len(r_["P"]):
                    msk = np.linalg.norm(r_["P"] - te, axis=1) <= 0.1
                    near += list(zip(r_["R"][msk], np.linalg.norm(r_["P"][msk] - te, axis=1)))
            if near:
                rr = np.array(near)
                rec["cap_principal_radius_within_0p1mm_of_te_point"] = {"points": len(rr), "min_mm": float(rr[:, 0].min()),
                                                                       "closest_point_mm": float(rr[:, 1].min())}
            rec["closes"] = bool(dist <= TOL["te_point_mm"])
        ok_sections = [a_ for a_ in approach if a_["error"] is None]
        rec["pass"] = bool(rec["closes"] and ok_sections and all(a_["pass"] for a_ in ok_sections))
        out["tips"].append(rec)
    out["pass"] = all(t["pass"] for t in out["tips"])
    return out


# ---------------------------------------------------------------------------
# 7. Summary
# ---------------------------------------------------------------------------

def summarize(report):
    req = []
    v = report["validity"]
    req.append(("1 validity", v["pass"],
                f"solids {v['solids']}, shells {v['shells']}, closed {v['shell_closed']}, valid {v['brepcheck_valid']}, "
                f"strict {v['brepcheck_strict_valid']}, faces {v['faces']}" + (f", faces straddling the attachment {v['straddling_faces']}" if v["straddling_faces"] else "")))
    e = v["extent"]
    req.append(("1 extent", e["pass"], f"{e['what']} {e['measured_mm']:.6f} mm vs spec {e['spec_mm']:.3f} (diff {e['difference_mm']:.2e}); "
                f"volume {v['volume_mm3']:.3f} mm3 (cap adds {v['cap_volume_mm3']:.3f}; OCCT default integration gives {v['volume_occt_default_mm3']:.3f})"))
    p = report["preservation"]
    req.append(("2 retained skin", p["pass"], f"max deviation {p['max_mm']:.2e} mm (parent->model {p['parent_to_model']['max_mm']:.2e}, "
                f"model->parent {p['model_to_parent']['max_mm']:.2e})"))
    o = report["outline"]
    scan = o["lip_scan"]
    worst = max((s.get("max_deviation_mm", np.inf) for s in scan.values()), default=np.inf)
    offs = max((s.get("max_abs_lip_outboard_offset_mm", np.nan) for s in scan.values()), default=np.nan)
    surf = o["outline_to_surface"]
    sd = max((s["max_mm"] for s in surf), default=0.0) if surf else None
    ok = worst <= TOL["outline_mm"] and (sd is None or sd <= TOL["outline_mm"])
    kind = next(iter(scan.values())).get("kind")
    what = {"plan": "outermost lip in plan", "midplane": "outermost lip", "rim": "outline-to-cap distance in the cut (3D rim)"}.get(kind, "")
    req.append(("2 locked outline", ok, f"max {what} deviation {worst:.2e} mm over {sum(s.get('planes', 0) for s in scan.values())} cut planes"
                + (f"; outline-to-surface max {sd:.2e} mm" if sd is not None else "")
                + (f"; outermost point of each cut sits up to {offs:.3g} mm from the rim along the rising span direction (informational)" if kind == "rim" else "")))
    js = report["joins"]["summary"]
    j = js.get("parent skin to cap")
    if j:
        req.append(("3 joins G1", report["joins"]["pass_G1"], f"{j['edges']} edges, {j['samples']} samples: gap {j['gap_mm']['max']:.2e} mm, "
                    f"angle {j['angle_deg']['max']:.2e} deg; curvature mismatch {j['curvature_mismatch_per_mm']['max']:.3e} /mm "
                    f"(cross-join normal curvature {j['cross_join_normal_curvature_mismatch_per_mm']['max']:.3e} /mm); "
                    f"closest sample to an edge end {j['closest_sample_to_edge_ends_mm']:.2g} mm, {j['samples_near_collapse']} samples near collapsed corners not judged"))
    else:
        req.append(("3 joins G1", False, "no parent-skin-to-cap edge found"))
    sm = js.get("cap seam")
    if sm:
        req.append(("3 cap seams G1", report["joins"]["cap_seams_pass_G1"], f"{sm['edges']} edges between curved cap faces: gap {sm['gap_mm']['max']:.2e} mm, "
                    f"angle {sm['angle_deg']['max']:.2e} deg, curvature mismatch {sm['curvature_mismatch_per_mm']['max']:.3e} /mm"))
    cuts = [c for rows in report["sections"].values() for c in rows]
    errs = [c for c in cuts if "error" in c]
    thick = [c for c in cuts if "error" not in c and c["zone"] == "thick"]
    thin = [c for c in cuts if "error" not in c and c["zone"] == "thin"]
    s_ok = all(c["pass"] for c in thick + thin) and not errs
    wt = min(thick, key=lambda c: c["min_radius"]["radius_mm"]) if thick else None
    wn = min(thin, key=lambda c: c["ratio_to_half_thickness"] or 0) if thin else None
    txt = f"{len(cuts)} cuts ({len(errs)} failed to trace)"
    if wt:
        txt += f"; thick zone min radius {wt['min_radius']['radius_mm']:.6g} mm at x/c {wt['fraction']:.6g} ({sum(not c['pass'] for c in thick)} of {len(thick)} below 0.4)"
    if wn:
        txt += f"; thin zone min ratio {wn['ratio_to_half_thickness']:.4g} at x/c {wn['fraction']:.8g}, t {wn['attachment_thickness_mm']:.4g} mm ({sum(not c['pass'] for c in thin)} of {len(thin)} below 0.9)"
    oblique = [c for c in thick + thin if not c["pass"] and (c["lip"].get("cut_obliquity_deg") or 0) > 60]
    if oblique:
        txt += f"; {len(oblique)} of the failing cuts cross the lip more than 60 deg from square"
    nrev = sum(len(c.get("reversals", [])) for c in cuts)
    txt += f"; bending reversals {nrev} on {sum(bool(c.get('reversals')) for c in cuts)} cuts"
    req.append(("4 section rule", s_ok, txt))
    ln_thick = [(c["lip"].get("lip_normal_radius_mm"), c["fraction"]) for c in thick if c["lip"].get("lip_normal_radius_mm")]
    ln_thin = [(c["lip"]["lip_normal_radius_mm"] / (0.5 * c["attachment_thickness_mm"]), c["fraction"]) for c in thin
               if c["lip"].get("lip_normal_radius_mm") and c["attachment_thickness_mm"] > 0]
    info = "radius across the lip, perpendicular to the locked outline, at each cut's lip point"
    if ln_thick:
        r_, f_ = min(ln_thick)
        info += f"; thick zone min {r_:.4g} mm at x/c {f_:.6g} ({sum(r < RULE['min_radius_mm'] for r, _ in ln_thick)} below 0.4)"
    if ln_thin:
        r_, f_ = min(ln_thin)
        info += f"; thin zone min ratio {r_:.4g} at x/c {f_:.8g} ({sum(r < RULE['thin_ratio_pass'] for r, _ in ln_thin)} below 0.9)"
    req.append(("4b lip-normal radius (informational)", None, info))
    cm = report["cap_map"]
    if cm.get("worst"):
        w = cm["worst"]
        zone = w.get("zone")
        wc = cm.get("worst_within_collapse_exclusion")
        req.append(("5 whole-cap map", cm["pass"], f"worst point radius {w['radius_mm']:.6g} mm at x/c {w['chord_fraction']:.6g}, t {w['attachment_thickness_mm']:.4g} mm, "
                    + (f"margin {w['margin_mm']:.4g} mm" if zone == "thick" else f"ratio {w.get('ratio_to_half_thickness'):.4g}")
                    + f"; tightest {cm['tightest']['radius_mm']:.6g} mm at x/c {cm['tightest']['chord_fraction']:.6g}"
                    + f"; points within {cm['collapse_exclusion_mm']} mm of collapsed corners not judged (grid reaches {cm['closest_grid_approach_to_collapse_mm']:.2g} mm"
                    + (f"; worst there {wc['radius_mm']:.3g} mm, score {wc['score_radius_over_required']:.3g})" if wc else ")")))
    else:
        req.append(("5 whole-cap map", False, "no cap point evaluated"))
    te = report["te_closure"]
    tt = []
    for t in te["tips"]:
        a = [x for x in t["approach_sections"] if x["error"] is None]
        last = a[-1] if a else None
        s = f"{t['tip']}: closes {t['closes']}"
        if last:
            s += (f", {last['distance_to_te_mm']} mm before TE: t {last['attachment_thickness_mm']:.4g}, section min R {last['min_radius_mm']:.4g}, "
                  f"ratio {last['ratio_to_half_thickness'] if last['ratio_to_half_thickness'] is None else round(last['ratio_to_half_thickness'], 4)}")
        if "curved_cap_principal_radius_within_0p1mm_of_te_segment" in t:
            r = t["curved_cap_principal_radius_within_0p1mm_of_te_segment"]
            s += f"; principal R within 0.1 mm of the TE face end min {r['min_mm']:.4g} median {r['median_mm']:.4g}"
        if "cap_principal_radius_within_0p1mm_of_te_point" in t:
            r = t["cap_principal_radius_within_0p1mm_of_te_point"]
            s += f"; principal R within 0.1 mm of the TE point min {r['min_mm']:.4g}"
        tt.append(s)
    req.append(("6 TE closure", te["pass"], te["kind"] + "; " + "; ".join(tt)))
    overall = all(r[1] for r in req if r[1] is not None)
    lines = [f"{report['shape']} {report['variant']}  {report['model']}  sha256 {report['step_sha256'][:12]}  family {report['family']}",
             f"OVERALL {'PASS' if overall else 'FAIL'}"]
    lines += [f"  [{'INFO' if ok is None else 'PASS' if ok else 'FAIL'}] {name}: {text}" for name, ok, text in req]
    return {"requirements": [{"name": n, "pass": ok, "detail": t} for n, ok, t in req], "pass": overall}, "\n".join(lines)


def audit(shape, variant, step, family="native", grid=241, quick=False):
    t0 = time.time()
    case = Case(shape, variant, family)
    model = Model(step, case)
    report = {"tool": "runs/cap-radius/audit/audit.py", "shape": shape, "variant": variant, "model": str(step),
              "step_sha256": hashlib.sha256(Path(step).read_bytes()).hexdigest(), "family": case.family,
              "family_definition": family_text(case), "rule": RULE, "tolerances": TOL}
    timing = {}

    def timed(name, fn):
        s = time.time()
        r = fn()
        timing[name] = round(time.time() - s, 1)
        return r

    report["validity"] = timed("validity", lambda: validity(model, case))
    report["preservation"] = timed("preservation", lambda: preservation(model, case))
    report["outline"] = {"lip_scan": timed("outline_scan", lambda: outline_lip_scan(model, case, n=25 if quick else 97)),
                         "outline_to_surface": timed("outline_surface", lambda: outline_surface_distance(model, case))}
    report["joins"] = timed("joins", lambda: joins(model, case))
    report["sections"] = timed("sections", lambda: sections(model, case))
    report["cap_map"] = timed("cap_map", lambda: cap_map(model, case, n_uniform=81 if quick else grid, n_cos=31 if quick else 81))
    report["te_closure"] = timed("te_closure", lambda: te_closure(model, case, report["joins"]))
    report["summary"], text = summarize(report)
    timing["total"] = round(time.time() - t0, 1)
    report["timing_s"] = timing
    return clean(report), text


def family_text(case):
    if case.family == "legacy":
        return "planes of runs/cap-trials/retest/audit.py (calibration)"
    return ("plane through the attachment chord-line point at chord fraction f, containing the attachment section's "
            "thickness direction and the tangent of the parent's constant-chord-fraction line d(chord point at f)/d(span)")


def section_table(report):
    lines = []
    for tip, rows in report["sections"].items():
        lines.append(f"  sections, {tip}:")
        lines.append("      x/c          t(mm)   zone  minR(mm)  at-arc  lipR(mm)   ratio  rev outline-dev(mm) obliq(deg) lipNormalR(mm) verdict  note")
        for c in rows:
            if "error" in c:
                lines.append(f"    {c['fraction']:<12.8g} {c['attachment_thickness_mm']:8.4g}  {c['zone']:5}  ERROR {c['error']}")
                continue
            ratio = c.get("ratio_to_half_thickness")
            od = (c.get("outline") or {}).get("deviation_mm")
            obl = c["lip"].get("cut_obliquity_deg")
            lines.append(f"    {c['fraction']:<12.8g} {c['attachment_thickness_mm']:8.4g}  {c['zone']:5} {c['min_radius']['radius_mm']:9.5g}"
                         f" {c['min_radius']['arclength_fraction']:6.3f} {c['lip']['radius_mm'] or float('nan'):9.5g}"
                         f" {ratio if ratio is not None else float('nan'):7.4g} {len(c['reversals']):3d} {od if od is not None else float('nan'):14.3e}"
                         f" {obl if obl is not None else float('nan'):9.2f} {c['lip'].get('lip_normal_radius_mm') or float('nan'):13.5g}"
                         f"   {'pass' if c['pass'] else 'FAIL'}   {c['label']}")
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("shape", choices=SHAPES)
    ap.add_argument("variant", choices=VARIANTS)
    ap.add_argument("step", type=Path)
    ap.add_argument("--out", type=Path, help="JSON output (default: results/SHAPE-VARIANT.json in this folder)")
    ap.add_argument("--family", choices=("native", "legacy"), default="native",
                    help="cut family; 'legacy' reproduces runs/cap-trials/retest/audit.py planes for calibration")
    ap.add_argument("--grid", type=int, default=241, help="uniform grid lines per face direction in the cap map")
    ap.add_argument("--quick", action="store_true", help="coarser outline scan and cap map (smoke tests only)")
    a = ap.parse_args(argv)
    out = a.out or HERE / "results" / f"{a.shape}-{a.variant}.json"
    report, text = audit(a.shape, a.variant, a.step, a.family, a.grid, a.quick)
    text = text + "\n" + section_table(report)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    out.with_suffix(".txt").write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0 if report["summary"]["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
