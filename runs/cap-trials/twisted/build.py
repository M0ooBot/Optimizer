"""Build a swept and washed-out wing with curvature-matched rolled tips."""

from pathlib import Path

import cadquery as cq
import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepLib import BRepLib
from OCP.Geom import Geom_BSplineSurface
from OCP.gp import gp_Pnt
from OCP.TColgp import TColgp_Array2OfPnt
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from scipy.interpolate import make_interp_spline


HERE = Path(__file__).resolve().parent
CUT, TIP = 440.0, 500.0
DIHEDRAL = np.tan(np.deg2rad(5))
WASHOUT = np.deg2rad(3)
ANGLE = np.pi / 2
U = np.linspace(0., 1., 121)


def airfoil(u):
    """Closed-TE NACA 2412. Smooth the two camber parabolas over 0.36–0.44c."""
    x = np.sin(np.pi * u / 2)**2
    forward = 0.02 / 0.4**2 * (0.8 * x - x*x)
    aft = 0.02 / 0.6**2 * (0.2 + 0.8*x - x*x)
    fraction = np.clip((x - 0.36) / 0.08, 0, 1)
    blend = fraction**3 * (10 - 15*fraction + 6*fraction**2)
    blend_derivative = 30*fraction**2 * (1-fraction)**2 / 0.08
    camber = forward + (aft-forward)*blend
    slope = (0.8 - 2*x) * (0.02/0.4**2 * (1-blend) + 0.02/0.6**2 * blend)
    slope += (aft-forward)*blend_derivative
    thickness = 0.6 * (0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x*x + 0.2843*x**3 - 0.1036*x**4)
    thickness = np.where((u == 0) | (u == 1), 0, thickness)
    camber = np.where((u == 0) | (u == 1), 0, camber)
    normal = np.stack((-slope, np.ones_like(slope)), axis=-1) / np.sqrt(1 + slope*slope)[..., None]
    return np.stack((x, camber), axis=-1), thickness[..., None]*normal


def section(y, u, thickness_factor):
    """Apply taper and true twist about the rising leading-edge reference line."""
    xy, normal_offset = airfoil(u)
    x, z = (xy + np.asarray(thickness_factor)[..., None]*normal_offset).T
    theta = WASHOUT * y/CUT
    chord = 180 - 80*y/CUT
    return np.stack((
        90*y/CUT + chord*(x*np.cos(theta) - z*np.sin(theta)),
        np.broadcast_to(y, np.shape(x)),
        y*DIHEDRAL + chord*(x*np.sin(theta) + z*np.cos(theta)),
    ), axis=-1)


def span_jet(y, u, thickness_factor):
    """The first and second span derivatives of the unchanged parent skin."""
    position = section(y, u, thickness_factor)
    chord = 180 - 80*y/CUT
    theta_rate = WASHOUT/CUT
    chord_rate = -80/CUT
    baseline = np.array([90/CUT, 1., DIHEDRAL])
    offset = position - np.asarray(y)[..., None]*baseline
    rotated = np.stack((-offset[:, 2], np.zeros(len(u)), offset[:, 0]), axis=-1)
    return (
        baseline + np.asarray(chord_rate/chord)[..., None]*offset + theta_rate*rotated,
        2*theta_rate*np.asarray(chord_rate/chord)[..., None]*rotated - theta_rate**2*offset,
    )


def reach(u):
    """Outboard distance to the common 3D perimeter for each chord guide."""
    return np.where((u == 0) | (u == 1), 0, (TIP-CUT)*np.sin(np.pi*u))


def cap_profiles(v):
    radius_progress = np.cos(ANGLE*v)
    y = CUT + reach(U)*radius_progress
    # One guide crosses the perimeter and connects both skin attachments.
    roll = np.sin(ANGLE*v)*np.sqrt(1+radius_progress**2)*(1+0.25*radius_progress**4)
    return section(y, U, roll)


def axis_arrays(spline):
    knots, mults = np.unique(spline.t, return_counts=True)
    values = TColStd_Array1OfReal(1, len(knots))
    multiplicities = TColStd_Array1OfInteger(1, len(knots))
    for index, (knot, mult) in enumerate(zip(knots, mults), start=1):
        values.SetValue(index, float(knot))
        multiplicities.SetValue(index, int(mult))
    return values, multiplicities


def face_from_guides(spline, mirror, is_cap):
    along = make_interp_spline(U, spline.c.transpose(1, 0, 2), k=5)
    net = along.c * [1, mirror, 1]
    # Make the constant-span boundary coordinates exact, and the cap's two
    # collapsed guide endpoints exactly coincident with their cut vertices.
    net[:, 0, 1] = mirror*CUT if is_cap else 0.
    net[:, -1, 1] = mirror*CUT
    if is_cap:
        net[0, :, :] = net[0, 0, :]
        net[-1, :, :] = net[-1, 0, :]
    poles = TColgp_Array2OfPnt(1, net.shape[0], 1, net.shape[1])
    for i, row in enumerate(net, start=1):
        for j, xyz in enumerate(row, start=1):
            poles.SetValue(i, j, gp_Pnt(*map(float, xyz)))
    u_knots, u_mults = axis_arrays(along)
    v_knots, v_mults = axis_arrays(spline)
    surface = Geom_BSplineSurface(poles, u_knots, v_knots, u_mults, v_mults, along.k, spline.k)
    # Keep both cap faces on the same underlying support at the perimeter.
    ranges = [(-1, 0), (0, 1)] if is_cap else [(0, CUT)]
    return [BRepBuilderAPI_MakeFace(surface, 0, 1, low, high, 1e-8).Face() for low, high in ranges]


def parent(side):
    v = np.linspace(0., CUT, 17)
    return make_interp_spline(
        v, np.array([section(y, U, side) for y in v]), k=5,
        bc_type=tuple(list(enumerate(span_jet(y, U, side), start=1)) for y in [0., CUT]),
    )


def cap():
    v = np.linspace(-1., 1., 161)
    d = reach(U)
    lower_first, lower_second = span_jet(CUT, U, -1)
    upper_first, upper_second = span_jet(CUT, U, 1)
    return make_interp_spline(
        v, np.array([cap_profiles(t) for t in v]), k=5,
        bc_type=(
            [(1, ANGLE*d[:, None]*lower_first), (2, (ANGLE*d[:, None])**2*lower_second)],
            [(1, -ANGLE*d[:, None]*upper_first), (2, (ANGLE*d[:, None])**2*upper_second)],
        ),
    )


def read_solid(path):
    solids = cq.importers.importStep(str(path)).val().Solids()
    if len(solids) != 1:
        raise RuntimeError(f"Expected one solid, found {len(solids)}")
    solid = solids[0]
    if len(solid.Shells()) != 1 or not BRep_Tool.IsClosed_s(solid.Shells()[0].wrapped):
        raise RuntimeError("Expected one closed shell")
    if not BRepCheck_Analyzer(solid.wrapped, True, False, True).IsValid():
        raise RuntimeError("Strict OCCT validity check failed")
    return solid


def main():
    sewing = BRepBuilderAPI_Sewing(2e-6)
    for spline, is_cap in [(parent(1), False), (parent(-1), False), (cap(), True)]:
        for mirror in (1, -1):
            for face in face_from_guides(spline, mirror, is_cap):
                sewing.Add(face)
    sewing.Perform()
    shell = cq.Shape.cast(sewing.SewedShape())
    if not isinstance(shell, cq.Shell) or sewing.NbFreeEdges() or sewing.NbMultipleEdges():
        raise RuntimeError(f"Sewing failed with {sewing.NbFreeEdges()} free and {sewing.NbMultipleEdges()} multiple edges")
    solid = cq.Solid.makeSolid(shell)
    if not BRepLib.OrientClosedSolid_s(solid.wrapped):
        raise RuntimeError("Could not orient the closed wing")
    path = HERE / "twisted-wing.step"
    cq.exporters.export(solid, str(path))
    solid = read_solid(path)
    print(f"Reimported {len(solid.Faces())} faces and {solid.Volume():,.3f} mm³")


if __name__ == "__main__":
    main()
