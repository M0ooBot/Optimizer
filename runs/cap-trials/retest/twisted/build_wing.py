
"""Build the wing as six sewn tensor-product B-spline faces."""
from pathlib import Path

import cadquery as cq
import numpy as np
from scipy.interpolate import make_interp_spline
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing
from OCP.Geom import Geom_BSplineSurface
from OCP.TColgp import TColgp_Array2OfPnt
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCP.gp import gp_Pnt

from parent_definition import CUT, reach, section, span_jet

HERE = Path(__file__).resolve().parent
CHORD_NODES = np.sort(np.r_[np.linspace(0, 1, 301), 2/np.pi*np.arcsin(np.sqrt([.36, .44]))])
PARENT_NODES = np.linspace(0, CUT, 13)
WRAP_NODES = np.linspace(-1, 1, 65)


def jet_boundary(first, second):
    return [(1, first), (2, second)]


def fit_surface(nodes, points, start, end):
    span_fit = make_interp_spline(nodes, points, k=5, bc_type=(start, end))
    chord_fit = make_interp_spline(CHORD_NODES, span_fit.c.transpose(1, 0, 2), k=5)
    poles = TColgp_Array2OfPnt(1, chord_fit.c.shape[0], 1, chord_fit.c.shape[1])
    for i, row in enumerate(chord_fit.c, 1):
        for j, xyz in enumerate(row, 1):
            poles.SetValue(i, j, gp_Pnt(*xyz))

    def occt_knots(spline):
        values, counts = np.unique(spline.t, return_counts=True)
        knots = TColStd_Array1OfReal(1, len(values))
        multiplicities = TColStd_Array1OfInteger(1, len(values))
        for index, (value, count) in enumerate(zip(values, counts), 1):
            knots.SetValue(index, value)
            multiplicities.SetValue(index, int(count))
        return knots, multiplicities

    chord_knots, chord_mult = occt_knots(chord_fit)
    span_knots, span_mult = occt_knots(span_fit)
    surface = Geom_BSplineSurface(poles, chord_knots, span_knots, chord_mult, span_mult, 5, 5)
    return cq.Face(BRepBuilderAPI_MakeFace(surface, 1e-8).Face())


def parent_face(sign):
    u = CHORD_NODES
    points = np.array([section(y, u, sign) for y in PARENT_NODES])
    start = jet_boundary(*span_jet(0, u, sign))
    end = jet_boundary(*span_jet(CUT, u, sign))
    return fit_surface(PARENT_NODES, points, start, end)


def thickness_factor(q):
    return q * (15 + q*q*(-10 + 3*q*q))/8


def cap_face():
    u = CHORD_NODES
    distance = reach(u)
    points = np.array([
        section(CUT + distance*(1-q*q), u, thickness_factor(q))
        for q in WRAP_NODES
    ])

    def boundary(q):
        tangent, curvature = span_jet(CUT, u, q)
        first = (-2*q*distance)[:, None]
        second = (-2*distance)[:, None]
        return jet_boundary(tangent*first, curvature*first**2 + tangent*second)

    return fit_surface(WRAP_NODES, points, boundary(-1), boundary(1))


def build():
    half_faces = [parent_face(1), parent_face(-1), cap_face()]
    sewing = BRepBuilderAPI_Sewing(1e-6)
    for face in half_faces:
        sewing.Add(face.wrapped)
        sewing.Add(face.mirror("XZ").wrapped)
    sewing.Perform()
    stitched = cq.Shape.cast(sewing.SewedShape())
    if stitched.ShapeType() != "Shell" or sewing.NbFreeEdges():
        raise RuntimeError(f"Expected closed shell: {stitched.ShapeType()}, {sewing.NbFreeEdges()} free edges")
    solid = cq.Solid.makeSolid(stitched)
    if solid.Volume() < 0:
        solid = cq.Solid(solid.wrapped.Reversed())
    if not solid.isValid() or solid.Volume() <= 0:
        raise RuntimeError("Invalid wing solid")
    cq.exporters.export(solid, str(HERE / "wing.step"))
    print(f"Valid solid with {len(solid.Faces())} faces, volume {solid.Volume():.6f} mm³")


if __name__ == "__main__":
    build()
