"""Build the swept fin and its rolled cap as four rational B-spline skin faces."""

from math import sqrt
from pathlib import Path

import cadquery as cq
import numpy as np
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing
from OCP.BRepCheck import BRepCheck_Analyzer, BRepCheck_Shell, BRepCheck_NoError
from OCP.Geom import Geom_BSplineSurface
from OCP.TColgp import TColgp_Array2OfPnt
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal, TColStd_Array2OfReal
from OCP.gp import gp_Pnt
from scipy.interpolate import make_interp_spline

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "vertical-fin.step"
CUT = 240.0
HEIGHT = 280.0
ROOT_CHORD = 160.0
CUT_CHORD = 85.0
LE_SWEEP = 65.0
THICKNESS_RATIO = 0.12
CENTER_SLOPE = (LE_SWEEP + CUT_CHORD / 2 - ROOT_CHORD / 2) / CUT
RADIUS_SLOPE = (CUT_CHORD - ROOT_CHORD) / (2 * CUT)
SEW_TOL = 1e-6


def center(z):
    return ROOT_CHORD / 2 + CENTER_SLOPE * z


def radius(z):
    return ROOT_CHORD / 2 + RADIUS_SLOPE * z


def arrays(values, kind):
    result = kind(1, len(values))
    for index, value in enumerate(values, 1):
        result.SetValue(index, value)
    return result


def skin(stations, knots, mults, degree, side):
    # Exact rational half-ellipse in the chordwise direction. The two patches
    # use the same sweep, taper and section scale in the height direction.
    ellipse_poles = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0)]
    ellipse_weights = [1, sqrt(0.5), 1, sqrt(0.5), 1]
    poles = TColgp_Array2OfPnt(1, 5, 1, len(stations))
    weights = TColStd_Array2OfReal(1, 5, 1, len(stations))
    for u, ((x, y), weight) in enumerate(zip(ellipse_poles, ellipse_weights), 1):
        for v, (mid, z, scale) in enumerate(stations, 1):
            poles.SetValue(u, v, gp_Pnt(mid + scale*x, side*THICKNESS_RATIO*scale*y, z))
            weights.SetValue(u, v, weight)
    return Geom_BSplineSurface(
        poles, weights,
        arrays([0.0, 0.5, 1.0], TColStd_Array1OfReal),
        arrays(list(map(float, knots)), TColStd_Array1OfReal),
        arrays([3, 2, 3], TColStd_Array1OfInteger),
        arrays(list(map(int, mults)), TColStd_Array1OfInteger),
        2, degree, False, False,
    )


def cap_stations():
    # R^2 = 1 - h^3 in normalized physical height h. This gives the same
    # position, slope and curvature as the parent at h=0. Near the apex,
    # w = sqrt(1-h) regularizes the radius for spline construction.
    parameter = np.linspace(0, 1, 33)
    w = 1 - parameter
    z = CUT + (HEIGHT - CUT)*(1 - w*w)
    scale = radius(z)*w*np.sqrt(3 - 3*w*w + w**4)
    along = 2*(HEIGHT-CUT)*np.array([CENTER_SLOPE, 1, RADIUS_SLOPE])
    at_join = [(1, along), (2, -along)]
    at_apex = [
        (1, [0, 0, -radius(HEIGHT)*sqrt(3)]),
        (2, -2*(HEIGHT-CUT)*np.array([CENTER_SLOPE, 1, 0])),
    ]
    spline = make_interp_spline(
        parameter, np.column_stack((center(z), z, scale)), k=5,
        bc_type=(at_join, at_apex),
    )
    knots, mults = np.unique(spline.t, return_counts=True)
    return spline.c, knots, mults, spline.k


def make_fin():
    root_edges = []
    sew = BRepBuilderAPI_Sewing(SEW_TOL)
    root_and_cut = np.array([[center(z), z, radius(z)] for z in [0, CUT]])
    cap = cap_stations()
    for side in (-1, 1):
        parent = skin(root_and_cut, [0, 1], [2, 2], 1, side)
        cap_surface = skin(*cap, side)
        root_edges.append(cq.Edge(BRepBuilderAPI_MakeEdge(parent.VIso(0)).Edge()))
        for surface in (parent, cap_surface):
            sew.Add(BRepBuilderAPI_MakeFace(surface, SEW_TOL).Face())
    root = cq.Face.makeFromWires(cq.Wire.assembleEdges(root_edges))
    sew.Add(root.wrapped)
    sew.Perform()
    shell = cq.Shape.cast(sew.SewedShape())
    if not isinstance(shell, cq.Shell) or BRepCheck_Shell(shell.wrapped).Closed() != BRepCheck_NoError:
        raise ValueError(f"Sewing gave {type(shell)} with {sew.NbFreeEdges()} free edges and {sew.NbDegeneratedShapes()} degeneracies, closed {shell.Closed()}")
    shell.wrapped.Closed(True)
    solid = cq.Solid.makeSolid(shell)
    if solid.Volume() < 0:
        solid = cq.Solid(solid.wrapped.Reversed())
    if not BRepCheck_Analyzer(solid.wrapped, True, False, True).IsValid():
        raise ValueError("The fin failed strict Open CASCADE validation")
    return solid


if __name__ == "__main__":
    solid = make_fin()
    solid.exportStep(str(OUTPUT), write_pcurves=True)
    print(f"Exported {OUTPUT.name}")
