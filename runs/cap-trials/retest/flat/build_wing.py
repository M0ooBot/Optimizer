"""Construct exact polynomial skins and roll surfaces for the 800 mm flat wing."""
from pathlib import Path

import cadquery as cq
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing
from OCP.BRepLib import BRepLib
from OCP.Geom import Geom_BezierSurface
from OCP.TColgp import TColgp_Array2OfPnt
from OCP.gp import gp_Pnt

OUT = Path(__file__).resolve().parent
# Degree-five Bernstein ordinates of the locked chordwise polynomials.
X = [0, 0, 0, 120, 120, 120]
D = [0, 40, 60, 60, 40, 0]
H = [0, 6.8, 8.2, 6.2, 2.8, 0]


def x_at(u):
    return 120*u**3*(10-15*u+6*u*u)


def h_at(u):
    return u*(1-u)*(34-20*u)


def d_at(u):
    return 200*u*(1-u)


def face(side=None, skin=1):
    """Return a central skin or one half of the cap at the indicated end."""
    if side is None:
        span_poles = [-350, 350]
        n = 2
    else:
        # f(t)=2t-t² and g(t)=1-2t³+t⁴ in a shared quartic basis.
        span_poles = [0, 1/2, 5/6, 1, 1]
        height_poles = [1, 1, 1, 1/2, 0]
        n = 5
    poles = TColgp_Array2OfPnt(1, 6, 1, n)
    for i, (x, d, h) in enumerate(zip(X, D, H), 1):
        for j, span in enumerate(span_poles, 1):
            y = span if side is None else side * (350 + d * span)
            z = skin * h * (1 if side is None else height_poles[j-1])
            poles.SetValue(i, j, gp_Pnt(x, y, z))
    return BRepBuilderAPI_MakeFace(Geom_BezierSurface(poles), 1e-9).Face()


def main():
    sewing = BRepBuilderAPI_Sewing(1e-7)
    for side in (None, -1, 1):
        for skin in (-1, 1):
            sewing.Add(face(side, skin))
    sewing.Perform()
    shell = cq.Shape.cast(sewing.SewedShape())
    assert isinstance(shell, cq.Shell) and sewing.NbFreeEdges() == 0, 'Sewing failed to close the shell'
    shell.wrapped.Closed(True)
    solid = cq.Solid.makeSolid(shell)
    assert BRepLib.OrientClosedSolid_s(solid.wrapped) and solid.isValid()
    cq.exporters.export(solid, str(OUT / 'wing.step'))


if __name__ == '__main__':
    main()
