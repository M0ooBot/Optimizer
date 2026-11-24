"""Build exact polynomial B-spline surfaces and export a sewn whole-wing solid."""
from math import comb, radians, tan
from pathlib import Path

import cadquery as cq
import numpy as np
from numpy.polynomial import Polynomial
from OCP.BRep import BRep_Tool
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.Geom import Geom_BezierSurface
from OCP.TColgp import TColgp_Array2OfPnt
from OCP.gp import gp_Pnt

HERE = Path(__file__).resolve().parent
STATION = 390.0
SWEEP = tan(radians(55))
TAPER = -115.0 / STATION
IDENTITY = Polynomial([0, 1])
ONE = Polynomial([1])
Q = IDENTITY**2
THICKNESS = 0.45 * (0.2969 * IDENTITY - 0.1260 * Q - 0.3516 * Q**2
                    + 0.2843 * Q**3 - 0.1036 * Q**4)


def product(u, t):
    """Tensor power coefficients for the product of two univariate polynomials."""
    return np.outer(u.coef, t.coef)


def combine(*terms):
    out = np.zeros(tuple(max(term.shape[i] for term in terms) for i in range(2)))
    for term in terms:
        out[:term.shape[0], :term.shape[1]] += term
    return out


def multiply(left, right):
    out = np.zeros(tuple(a + b - 1 for a, b in zip(left.shape, right.shape)))
    for (i, j), value in np.ndenumerate(left):
        out[i:i+right.shape[0], j:j+right.shape[1]] += value * right
    return out


def controls(power, degree):
    """Express a power-basis polynomial in the degree-elevated Bernstein basis."""
    n, m = degree
    power_to_bernstein = [
        np.array([[comb(i, j) / comb(k, j) if j <= i else 0
                   for j in range(size)] for i in range(k + 1)])
        for k, size in zip(degree, power.shape)
    ]
    return power_to_bernstein[0] @ power @ power_to_bernstein[1].T


def surface(side, skin, cap):
    if cap:
        rise = 240.0 * Q * (ONE - Q)
        y = combine(product(ONE, STATION * ONE), product(rise, 2*IDENTITY - IDENTITY**2))
        chord = combine(product(ONE, 210.0 * ONE), TAPER * y)
        fairing = 3*Q**2 - 2*Q**3 - Q
        ease = 10*IDENTITY**3 - 15*IDENTITY**4 + 6*IDENTITY**5
        chord_position = combine(product(Q, ONE), product(fairing, ease))
        # Ease the forward thickness shoulder without changing the lip tangent or radius.
        roll = ONE - IDENTITY**3 - 4.4*IDENTITY**3*(ONE - IDENTITY)**2
        z = multiply(chord, product(THICKNESS, roll))
    else:
        y = product(ONE, STATION * IDENTITY)
        chord = combine(product(ONE, 210.0 * ONE), TAPER * y)
        chord_position = product(Q, ONE)
        z = multiply(chord, product(THICKNESS, ONE))
    x = combine(SWEEP * y, multiply(chord, chord_position))
    coordinates = [x, side * y, skin * z]
    degree = tuple(max(poly.shape[d] for poly in coordinates) - 1 for d in range(2))
    xx, yy, zz = [controls(poly, degree) for poly in coordinates]
    if cap:
        # Enforce exact collapsed terminal edges and the common planar lip.
        for row in (0, -1):
            xx[row, :] = xx[row, 0]
            yy[row, :] = yy[row, 0]
            zz[row, :] = 0
        zz[:, -1] = 0
    else:
        zz[[0, -1], :] = 0
    poles = TColgp_Array2OfPnt(1, degree[0] + 1, 1, degree[1] + 1)
    for i, j in np.ndindex(xx.shape):
        poles.SetValue(i + 1, j + 1, gp_Pnt(float(xx[i, j]), float(yy[i, j]), float(zz[i, j])))
    return Geom_BezierSurface(poles)


def build():
    sewing = BRepBuilderAPI_Sewing(1e-6)
    for side in (-1, 1):
        for skin in (-1, 1):
            for cap in (False, True):
                patch = surface(side, skin, cap)
                face = BRepBuilderAPI_MakeFace(patch, 1e-8).Face()
                sewing.Add(face)
    sewing.Perform()
    assert sewing.NbFreeEdges() == 0, sewing.NbFreeEdges()
    shells = cq.Shape.cast(sewing.SewedShape()).Shells()
    assert len(shells) == 1 and BRep_Tool.IsClosed_s(shells[0].wrapped)
    solid = cq.Solid.makeSolid(shells[0])
    assert BRepCheck_Analyzer(solid.wrapped, True, False, True).IsValid()
    assert solid.Volume() > 0
    cq.exporters.export(solid, str(HERE / "swept55-wing.step"))


if __name__ == "__main__":
    build()
