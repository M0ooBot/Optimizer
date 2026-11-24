"""Build a 120 mm chord, 800 mm span wing with exact polynomial tip surfaces."""
from math import comb
from pathlib import Path

import cadquery as cq
import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.Geom import Geom_BezierSurface
from OCP.TColgp import TColgp_Array2OfPnt
from OCP.gp import gp_Pnt

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "flat-wing.step"
CHORD = 120.0
CUT = 350.0
CAP = 50.0

# Increasing u goes aft. Both airfoil edges have a common vertical tangent.
X = np.array([0, 0, 0, 10, -15, 6]) * CHORD
HEIGHT = np.array([0, 34, -54, 20, 0, 0])
RUNOUT = np.array([0, 4, -4, 0, 0, 0]) * CAP


def bernstein(power):
    degree = len(power) - 1
    return np.array([
        sum(power[j] * comb(i, j) / comb(degree, j) for j in range(i + 1))
        for i in range(degree + 1)
    ])


def make_face(poles):
    nu, nv, _ = poles.shape
    grid = TColgp_Array2OfPnt(1, nu, 1, nv)
    for i, row in enumerate(poles, 1):
        for j, point in enumerate(row, 1):
            grid.SetValue(i, j, gp_Pnt(*point))
    surface = Geom_BezierSurface(grid)
    return BRepBuilderAPI_MakeFace(surface, 1e-7).Face()


def construct():
    faces = []
    x = bernstein(X)
    runout = bernstein(RUNOUT)
    for skin_sign in (1, -1):
        h = skin_sign * bernstein(HEIGHT)
        parent = np.empty((len(x), 2, 3))
        parent[:, 0, :] = np.column_stack((x, np.full(len(x), -CUT), h))
        parent[:, 1, :] = np.column_stack((x, np.full(len(x), CUT), h))
        faces.append(make_face(parent))

        # At v=0 the cap control row is the actual parent end-pole row.
        # At v=1 upper and lower use the same perimeter and form a vertical tangent plane.
        for column, side in ((0, -1), (1, 1)):
            cut_poles = parent[:, column, :]
            roll = np.empty((len(x), 6, 3))
            roll[:, :, 0] = cut_poles[:, 0, None]
            roll[:, :, 1] = cut_poles[:, 1, None] + side * np.outer(
                runout, bernstein([0, 2, -1, 0, 0, 0])
            )
            roll[:, :, 2] = np.outer(cut_poles[:, 2], bernstein([1, 0, 0, -2.5, 1.875, -0.375]))
            faces.append(make_face(roll))

    sew = BRepBuilderAPI_Sewing(1e-7)
    for face in faces:
        sew.Add(face)
    sew.Perform()
    shells = cq.Shape.cast(sew.SewedShape()).Shells()
    if len(shells) != 1 or not BRep_Tool.IsClosed_s(shells[0].wrapped):
        raise RuntimeError("Expected a single sewn, closed shell")
    solid = cq.Solid.makeSolid(shells[0])
    if not BRepCheck_Analyzer(solid.wrapped, True, False, True).IsValid():
        raise RuntimeError("Strict OCCT validity check failed")
    return solid


if __name__ == "__main__":
    cq.exporters.export(construct(), str(OUTPUT))
