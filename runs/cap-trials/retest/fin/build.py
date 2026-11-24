from pathlib import Path

import cadquery as cq
import numpy as np
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing
from OCP.Geom import Geom_BSplineSurface
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal, TColStd_Array2OfReal
from OCP.TColgp import TColgp_Array2OfPnt
from OCP.gp import gp_Pnt
from numpy.polynomial import Polynomial
from scipy.interpolate import make_interp_spline

HERE = Path(__file__).resolve().parent
ATTACHMENT = 240.0
M_SLOPE = 27.5 / ATTACHMENT
A_SLOPE = -37.5 / ATTACHMENT
SECTION_RATIO = 0.12


def cap_span(reach):
    h = Polynomial([0, 2, -1])
    z = ATTACHMENT + reach * h
    m = 80 + M_SLOPE * z
    a = 80 + A_SLOPE * z
    correction = {40: 0.0, 5: 0.50}[reach]
    q = 1 + h + h**2
    x = a * Polynomial([1, -1])
    y = SECTION_RATIO * x * (1 - correction * h**3)

    def eval_span(t, order=0):
        root = np.sqrt(q(t))
        if order == 0:
            return np.array([m(t), x(t) * root, y(t) * root, z(t)]).T
        root_d = q.deriv()(t) / (2 * root)
        if order == 1:
            return np.array([m.deriv()(t), x.deriv()(t) * root + x(t) * root_d,
                             y.deriv()(t) * root + y(t) * root_d, z.deriv()(t)]).T
        root_dd = q.deriv(2)(t) / (2 * root) - root_d**2 / root
        return np.array([m.deriv(2)(t),
                         x.deriv(2)(t) * root + 2 * x.deriv()(t) * root_d + x(t) * root_dd,
                         y.deriv(2)(t) * root + 2 * y.deriv()(t) * root_d + y(t) * root_dd,
                         z.deriv(2)(t)]).T

    nodes = np.linspace(0, 1, 33)
    spline = make_interp_spline(nodes, eval_span(nodes), k=5, bc_type=(
        [(1, eval_span(0, 1)), (2, eval_span(0, 2))],
        [(1, eval_span(1, 1)), (2, eval_span(1, 2))]
    ))
    knots, multiplicities = np.unique(spline.t, return_counts=True)
    return spline.c, knots, multiplicities, spline.k


def occ_array(values, cls):
    result = cls(1, len(values))
    for i, value in enumerate(values, 1):
        result.SetValue(i, value.item() if isinstance(value, np.generic) else value)
    return result


def side(span, sign):
    controls, knots, multiplicities, degree = span
    angular_controls = [(1, 0), (1, sign), (0, sign), (-1, sign), (-1, 0)]
    angular_weights = [1, 2**-0.5, 1, 2**-0.5, 1]
    poles = TColgp_Array2OfPnt(1, 5, 1, len(controls))
    weights = TColStd_Array2OfReal(1, 5, 1, len(controls))
    for i, ((cosine, sine), weight) in enumerate(zip(angular_controls, angular_weights), 1):
        for j, (m, a, b, z) in enumerate(controls, 1):
            poles.SetValue(i, j, gp_Pnt(m + a * cosine, b * sine, z))
            weights.SetValue(i, j, weight)
    surface = Geom_BSplineSurface(
        poles, weights, occ_array([0.0, 0.5, 1.0], TColStd_Array1OfReal),
        occ_array(knots, TColStd_Array1OfReal),
        occ_array([3, 2, 3], TColStd_Array1OfInteger),
        occ_array(multiplicities, TColStd_Array1OfInteger), 2, degree)
    return BRepBuilderAPI_MakeFace(surface, 1e-7).Face()


def build(reach):
    parent = (np.array([[80, 80, 9.6, 0.0], [107.5, 42.5, 5.1, ATTACHMENT]]),
              [0.0, 1.0], [2, 2], 1)
    cap = cap_span(reach)
    ellipse = cq.Workplane("XY", origin=(80, 0, 0)).ellipse(80, 9.6).val()
    root = cq.Face.makeFromWires(ellipse)
    sew = BRepBuilderAPI_Sewing(1e-6)
    sew.Add(root.wrapped)
    for sign in [-1, 1]:
        sew.Add(side(parent, sign))
        sew.Add(side(cap, sign))
    sew.Perform()
    shell = cq.Shape.cast(sew.SewedShape())
    assert isinstance(shell, cq.Shell), shell.ShapeType()
    assert shell.isValid() and sew.NbFreeEdges() == sew.NbMultipleEdges() == 0
    solid = cq.Solid.makeSolid(shell)
    assert solid.isValid() and solid.Volume() > 0
    return solid


if __name__ == "__main__":
    for reach in [40, 5]:
        folder = HERE / str(reach)
        folder.mkdir(exist_ok=True)
        solid = build(reach)
        cq.exporters.export(solid, str(folder / "wing.step"))
        print(reach, len(solid.Faces()), solid.Volume())
