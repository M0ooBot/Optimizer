"""Export each analytic parent as one closed STEP solid.

Every skin is interpolated by a quintic tensor-product B-spline through the
exact parent.py surface. Node counts double until the surface matches the
exact definition to FIT_TOLERANCE at the midpoints between nodes. A blunt
trailing edge is the ruled face between the two skins' trailing-edge curves.
The attachment (and, for fins, the root) is closed with a planar face bounded
by the skins' own isocurves.
"""

import json
import sys

import cadquery as cq
import numpy as np
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepLib import BRepLib
from OCP.Geom import Geom_BSplineSurface
from OCP.GeomFill import GeomFill
from OCP.TColgp import TColgp_Array2OfPnt
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCP.gp import gp_Pnt
from scipy.interpolate import BSpline, make_interp_spline

from common import PARENTS, load_parent

FIT_TOLERANCE = 2e-7  # mm
SEW_TOLERANCE = 1e-6
KNOT_TOLERANCE = 1e-7
DEGREE = 5


def chord_nodes(n):
    # Cluster toward both ends of the chord parameter.
    return (1 - np.cos(np.linspace(0, np.pi, n))) / 2


def occ_knots(spline):
    values, counts = np.unique(spline.t, return_counts=True)
    knots = TColStd_Array1OfReal(1, len(values))
    mults = TColStd_Array1OfInteger(1, len(values))
    for i, (value, count) in enumerate(zip(values, counts), 1):
        knots.SetValue(i, float(value))
        mults.SetValue(i, int(count))
    return knots, mults


class Fit:
    """A tensor-product interpolant of one skin over chord parameter s and span station.

    Where the shape lists CHORD_BREAKS (lines across which the skin is only
    C1), each chord segment is interpolated separately and the segments are
    joined with full-multiplicity knots, so no quintic piece spans a break.
    """

    def __init__(self, module, side, span_lo, span_hi, n_s, n_y):
        breaks = [0.0, *getattr(module, "CHORD_BREAKS", ()), 1.0]
        segments = []
        for a, b in zip(breaks, breaks[1:]):
            n = max(9, int(round(n_s * (b - a))) | 1)
            segments.append(a + (b - a) * chord_nodes(n))
        self.s_nodes = np.unique(np.concatenate(segments))
        self.y_nodes = np.linspace(span_lo, span_hi, n_y)
        k_y = min(DEGREE, n_y - 1)
        knots, poles = [], []
        for index, nodes in enumerate(segments):
            S, Y = np.meshgrid(nodes, self.y_nodes, indexing="ij")
            along_y = make_interp_spline(self.y_nodes, module.skin(S, Y, side), k=k_y, axis=1)
            along_s = make_interp_spline(nodes, along_y.c.transpose(1, 0, 2), k=DEGREE, axis=0)
            t = along_s.t
            first = index == 0
            last = index == len(segments) - 1
            if first:
                knots.append(t[: len(t) - (DEGREE + 1)])
            else:
                knots.append(np.full(DEGREE, nodes[0]))
                knots.append(t[DEGREE + 1: len(t) - (DEGREE + 1)])
            if last:
                knots.append(t[len(t) - (DEGREE + 1):])
            poles.append(along_s.c if first else along_s.c[1:])
        self.y_spline = along_y
        self.s_spline = BSpline(np.concatenate(knots), np.concatenate(poles), DEGREE)
        self.poles = self.s_spline.c  # (poles in s, poles in y, 3)

    def __call__(self, s, y):
        bs = BSpline.design_matrix(np.clip(s, 0, 1), self.s_spline.t, self.s_spline.k).toarray()
        by = BSpline.design_matrix(np.clip(y, self.y_nodes[0], self.y_nodes[-1]),
                                   self.y_spline.t, self.y_spline.k).toarray()
        return np.einsum("pi,ijc,pj->pc", bs, self.poles, by)

    def surface(self):
        ns, ny, _ = self.poles.shape
        poles = TColgp_Array2OfPnt(1, ns, 1, ny)
        for i in range(ns):
            for j in range(ny):
                poles.SetValue(i + 1, j + 1, gp_Pnt(*map(float, self.poles[i, j])))
        uk, um = occ_knots(self.s_spline)
        vk, vm = occ_knots(self.y_spline)
        surface = Geom_BSplineSurface(poles, uk, vk, um, vm, self.s_spline.k, self.y_spline.k)
        # Segments meet with matching position; lower the joint multiplicity to
        # C1 where the tangent planes agree to within KNOT_TOLERANCE.
        self.c1_joints = []
        for index in range(1, surface.NbUKnots() + 1):
            if surface.UMultiplicity(index) == DEGREE:
                self.c1_joints.append(bool(surface.RemoveUKnot(index, DEGREE - 1, KNOT_TOLERANCE)))
        return surface


def fit_error(fit, module, side):
    """Largest error between chord nodes (at span nodes) and between span nodes (at chord nodes)."""
    s_mid = (fit.s_nodes[1:] + fit.s_nodes[:-1]) / 2
    y_mid = (fit.y_nodes[1:] + fit.y_nodes[:-1]) / 2
    ys = np.unique(np.r_[fit.y_nodes[::max(1, len(fit.y_nodes) // 4)], fit.y_nodes[-1]])
    S, Y = np.meshgrid(s_mid, ys, indexing="ij")
    error_s = np.linalg.norm(fit(S.ravel(), Y.ravel()) - module.skin(S, Y, side).reshape(-1, 3), axis=1).max()
    ss = np.unique(np.r_[fit.s_nodes[::2], 1.0])
    S, Y = np.meshgrid(ss, y_mid, indexing="ij")
    error_y = np.linalg.norm(fit(S.ravel(), Y.ravel()) - module.skin(S, Y, side).reshape(-1, 3), axis=1).max()
    return error_s, error_y


def fit_skin(module, side, span_lo, span_hi):
    n_s, n_y = 33, 5
    while True:
        fit = Fit(module, side, span_lo, span_hi, n_s, n_y)
        error_s, error_y = fit_error(fit, module, side)
        done_s, done_y = error_s <= FIT_TOLERANCE or n_s >= 1025, error_y <= FIT_TOLERANCE or n_y >= 65
        if done_s and done_y:
            return fit, max(error_s, error_y)
        if not done_s:
            n_s = 2 * n_s - 1
        if not done_y:
            n_y = 2 * n_y - 1


def edge(curve):
    return cq.Edge(BRepBuilderAPI_MakeEdge(curve).Edge())


def planar_face(edges):
    wire = cq.Wire.assembleEdges(edges)
    return cq.Face(BRepBuilderAPI_MakeFace(wire.wrapped, True).Face())


def half_faces(module, span_lo, span_hi, blunt, close):
    """Skins (and TE face) between two stations, plus planar faces at the stations in close."""
    fits, faces, errors = {}, [], {}
    for side in (+1, -1):
        fits[side], errors[side] = fit_skin(module, side, span_lo, span_hi)
        faces.append(cq.Face(BRepBuilderAPI_MakeFace(fits[side].surface(), SEW_TOLERANCE).Face()))
    surfaces = {side: fits[side].surface() for side in fits}
    if blunt:
        ruled = GeomFill.Surface_s(surfaces[-1].UIso(1.0), surfaces[+1].UIso(1.0))
        faces.append(cq.Face(BRepBuilderAPI_MakeFace(ruled, SEW_TOLERANCE).Face()))
    for station in close:
        edges = [edge(surfaces[side].VIso(station)) for side in (+1, -1)]
        if blunt:
            upper = surfaces[+1].Value(1.0, station)
            lower = surfaces[-1].Value(1.0, station)
            edges.append(cq.Edge.makeLine(cq.Vector(lower.X(), lower.Y(), lower.Z()),
                                          cq.Vector(upper.X(), upper.Y(), upper.Z())))
        faces.append(planar_face(edges))
    return faces, fits, errors


def solid_from(faces):
    sewing = BRepBuilderAPI_Sewing(SEW_TOLERANCE)
    for face in faces:
        sewing.Add(face.wrapped)
    sewing.Perform()
    shell = cq.Shape.cast(sewing.SewedShape())
    if not isinstance(shell, cq.Shell) or sewing.NbFreeEdges() or sewing.NbMultipleEdges():
        raise RuntimeError(f"Sewing gave {type(shell).__name__}, free edges {sewing.NbFreeEdges()}, "
                           f"multiple edges {sewing.NbMultipleEdges()}")
    shell.wrapped.Closed(True)
    solid = cq.Solid.makeSolid(shell)
    if not BRepLib.OrientClosedSolid_s(solid.wrapped):
        raise RuntimeError("Could not orient the closed solid")
    if not BRepCheck_Analyzer(solid.wrapped, True, False, True).IsValid():
        bad = [(face.geomType(), [str(s) for s in BRepCheck_Analyzer(face.wrapped, True, False, True)
                                  .Result(face.wrapped).Status()])
               for face in solid.Faces() if not BRepCheck_Analyzer(face.wrapped, True, False, True).IsValid()]
        raise RuntimeError(f"Solid failed strict BRepCheck: {bad}")
    return solid


def precise_volume(solid):
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(solid.wrapped, props, 1e-9, True, False)
    return props.Mass()


def build(shape, variant):
    module = load_parent(shape, variant)
    blunt = module.TE_THICKNESS > 0
    lo, hi = module.SPAN_LIMITS
    if module.FULL_SPAN and shape == "flat":
        faces, fits, errors = half_faces(module, lo, hi, blunt, close=[lo, hi])
        pieces = [(lo, hi, errors)]
    elif module.FULL_SPAN:
        right, fits, errors = half_faces(module, 0.0, hi, blunt, close=[hi])
        faces = right + [face.mirror("XZ") for face in right]
        pieces = [(0.0, hi, errors)]
    else:
        faces, fits, errors = half_faces(module, lo, hi, blunt, close=[lo, hi])
        pieces = [(lo, hi, errors)]
    solid = solid_from(faces)
    out = PARENTS / shape / variant / "parent.step"
    cq.exporters.export(solid, str(out))
    report = {
        "faces": len(solid.Faces()),
        "volume_mm3": precise_volume(solid),
        "fit_tolerance_mm": FIT_TOLERANCE,
        "skin_fit_max_midpoint_error_mm": max(max(e.values()) for _, _, e in pieces),
        "chord_break_joints_reduced_to_C1": {("upper" if side > 0 else "lower"): getattr(fits[side], "c1_joints", [])
                                              for side in fits},
        "skin_poles": {("upper" if side > 0 else "lower"): list(fits[side].poles.shape[:2]) for side in fits},
        "step_bytes": out.stat().st_size,
    }
    (PARENTS / shape / variant / "build.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    targets = sys.argv[1:] or ["flat", "twisted", "swept55", "fin40", "fin5"]
    for shape in targets:
        for variant in ("sharp", "blunt"):
            print(shape, variant, build(shape, variant), flush=True)
