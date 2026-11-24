"""Build a boundary-surface tip on the saved, untwisted Wing B STEP."""

from dataclasses import dataclass
from pathlib import Path
from math import comb

import cadquery as cq
import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.BRepLib import BRepLib
from OCP.Geom import Geom_BSplineSurface
from OCP.GProp import GProp_GProps
from OCP.gp import gp_Dir, gp_Pln, gp_Pnt, gp_Vec
from OCP.TColgp import TColgp_Array2OfPnt
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from scipy.interpolate import make_interp_spline
from scipy.optimize import brentq


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "paired-wing-trial/a/measured/01-wing-occt.step"
OUTPUT = HERE / "wing-b-rounded.step"
CUT = 540.0
SPAN = 600.0
BOW_COEFFICIENT = 60.0 / SPAN**2
PLAN_CONTROL_ARM = (SPAN - CUT) * 4 / 3
SEW_TOLERANCE = 0.0002
AXES = np.array([-1.0, 1.0, -1.0])  # STEP XYZ to span, height, aft


def read_solid(path):
    solids = cq.importers.importStep(str(path)).val().Solids()
    if len(solids) != 1 or not BRepCheck_Analyzer(solids[0].wrapped).IsValid():
        raise ValueError(f"Expected one valid solid in {path}")
    return solids[0]


@dataclass
class Cut:
    skin: object

    def __post_init__(self):
        self.v = brentq(lambda v: self.point(0, v)[0] - CUT, 0, 1, xtol=1e-14)
        self.nose_u = brentq(lambda u: self.differential(u)[0][2], 0.47, 0.53, xtol=1e-14)
        self.nose = self.point(self.nose_u, self.v)
        self.tail = self.point(0, self.v)
        self.span_slopes = [self.span_slope(u) for u in (self.nose_u, 0.0)]
        self.nose[0] = self.tail[0] = CUT
        _, nose_span = self.differential(self.nose_u)
        tail_lower, tail_span = self.differential(0.0)
        tail_upper, _ = self.differential(1.0)
        bow_slope = 2 * BOW_COEFFICIENT * CUT
        tail_height_per_chord = (tail_lower[1] / tail_lower[2] + tail_upper[1] / tail_upper[2]) / 2
        tail_chord_per_u = (self.tail[2] - self.nose[2]) * np.pi / 2
        self.rim_slope_correction = (
            3 * PLAN_CONTROL_ARM * (nose_span[1] / nose_span[0] - bow_slope),
            -3 * PLAN_CONTROL_ARM * (tail_span[1] / tail_span[0] - bow_slope)
            - tail_height_per_chord * tail_chord_per_u,
        )

    def point(self, u, v):
        return np.array(self.skin.Value(u, v).Coord()) * AXES

    def differential(self, u):
        point, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
        self.skin.D1(u, self.v, point, du, dv)
        return np.array(du.Coord()) * AXES, np.array(dv.Coord()) * AXES

    def span_slope(self, u):
        _, dv = self.differential(u)
        return dv[2] / dv[0]

    def profiles(self, chord):
        if chord <= self.nose[2]:
            return self.nose[None, :].repeat(2, axis=0), [self.nose_u, self.nose_u]
        if chord >= self.tail[2]:
            return self.tail[None, :].repeat(2, axis=0), [1.0, 0.0]

        def chord_error(u):
            return self.point(u, self.v)[2] - chord

        upper = brentq(chord_error, self.nose_u, 1, xtol=1e-14)
        lower = brentq(chord_error, 0, self.nose_u, xtol=1e-14)
        positions = np.array([self.point(u, self.v) for u in (upper, lower)])
        positions[:, 0] = CUT
        return positions, [upper, lower]

    def height_change(self, u, span, chord):
        du, dv = self.differential(u)
        along_span = span / dv[0]
        along_profile = (chord - along_span * dv[2]) / du[2]
        return along_span * dv[1] + along_profile * du[1]


def plan_perimeter(cut, u):
    """Control arms lie along the original plan tangents at both ends."""
    arm = PLAN_CONTROL_ARM
    start = cut.nose[[0, 2]]
    stop = cut.tail[[0, 2]]
    nose_slope, tail_slope = cut.span_slopes
    controls = np.array(
        [start, start + [arm, arm * nose_slope], stop + [arm, arm * tail_slope], stop]
    )
    return sum(
        comb(3, i) * u**i * (1 - u) ** (3 - i) * point
        for i, point in enumerate(controls)
    )


def rim_height_adjustment(cut, u):
    """Continue the existing 3D edge tangents through the rim endpoints."""
    def taper(distance, extent):
        x = min(distance / extent, 1.0)
        return 1 - x**3 * (10 - 15 * x + 6 * x**2)

    nose, tail = cut.rim_slope_correction
    return nose * u * taper(u, 0.25) + tail * (u - 1) * taper(1 - u, 0.35)


def wrap_guide(cut, u, t):
    """Two branches of one wrap from upper cut through perimeter to lower cut."""
    chord = cut.nose[2] + (cut.tail[2] - cut.nose[2]) * (1 - np.cos(np.pi * u / 2))
    cut_points, parameters = cut.profiles(chord)
    if u in (0.0, 1.0):
        points = cut_points[:, None, :].repeat(len(t), axis=1)
        return points, np.zeros((2, 3)), np.zeros((2, 3))

    origin = cut_points[0, [0, 2]]
    direction = plan_perimeter(cut, u) - origin
    span, chord_shift = direction
    cut_slopes = np.array([cut.height_change(p, span, chord_shift) for p in parameters])
    half_thickness = (cut_points[0, 1] - cut_points[1, 1]) / 2
    height = cut_points[:, 1].mean()
    if half_thickness <= 0:
        raise ValueError(f"Invalid cut thickness at perimeter parameter {u}")

    # The middle of the rim follows the saved quadratic front-view bow.
    # Corrections vanish toward the rim and impose the existing skin tangents at the cut.
    remaining = 1 - t
    extension = 1 - remaining**2
    plan = origin + extension[:, None] * direction
    bow = BOW_COEFFICIENT * (plan[:, 0] ** 2 - CUT**2)
    bow_slope = 2 * BOW_COEFFICIENT * CUT * span
    fairing = extension * (1 - extension) ** 2
    center = (
        height + bow + (cut_slopes.mean() - bow_slope) * fairing
        + rim_height_adjustment(cut, u) * extension**2
    )

    # The terminal radius follows the local lip height and eases down near the TE.
    # This form gives finite curvature at the outer perimeter and zero height there.
    tip_thinning = 0.35
    radius = 0.7 * half_thickness * (1 - 0.3 * u**4)
    lip_height = half_thickness * np.exp(-tip_thinning)
    rounding = lip_height**2 / (2 * np.linalg.norm(direction) * radius)
    thickness_slope = (cut_slopes[0] - cut_slopes[1]) / (2 * half_thickness) + rounding / 2
    thickness = (
        half_thickness
        * np.exp(thickness_slope * fairing - tip_thinning * extension**2)
        * remaining
        / np.sqrt(remaining**2 + rounding * (1 - remaining**2))
    )

    points = np.array(
        [np.column_stack((plan[:, 0], center + side * thickness, plan[:, 1])) for side in (1, -1)]
    )
    cut_tangent = 2 * np.array([[span, slope, chord_shift] for slope in cut_slopes])
    rim_speed = lip_height / np.sqrt(rounding)
    rim_tangent = np.array([[0, -rim_speed, 0], [0, rim_speed, 0]])
    return points, cut_tangent, rim_tangent


def knot_arrays(spline):
    knots, mults = np.unique(spline.t, return_counts=True)
    occt_knots = TColStd_Array1OfReal(1, len(knots))
    occt_mults = TColStd_Array1OfInteger(1, len(knots))
    for index, (knot, mult) in enumerate(zip(knots, mults), start=1):
        occt_knots.SetValue(index, float(knot))
        occt_mults.SetValue(index, int(mult))
    return occt_knots, occt_mults


def boundary_faces(cut):
    # Sample the continuous wrap guides in their natural boundary parameters.
    # Denser sampling near the endpoints keeps the cut and perimeter fits within sewing tolerance.
    u = (1 - np.cos(np.linspace(0, np.pi, 201))) / 2
    u = np.sort(np.concatenate((u, np.linspace(0, 0.018, 41)[1:])))
    t = (1 - np.cos(np.linspace(0, np.pi, 49))) / 2
    guides, cut_tangents, rim_tangents = [np.array(a) for a in zip(*(wrap_guide(cut, p, t) for p in u))]

    faces = []
    for branch in range(2):
        across = make_interp_spline(
            t,
            guides[:, branch].transpose(1, 0, 2),
            axis=0,
            bc_type=([(1, cut_tangents[:, branch])], [(1, rim_tangents[:, branch])]),
        )
        along = make_interp_spline(u, across.c.transpose(1, 0, 2), axis=0)
        u_knots, u_mults = knot_arrays(along)
        v_knots, v_mults = knot_arrays(across)
        net = along.c * AXES
        poles = TColgp_Array2OfPnt(1, net.shape[0], 1, net.shape[1])
        for i, row in enumerate(net, start=1):
            for j, xyz in enumerate(row, start=1):
                poles.SetValue(i, j, gp_Pnt(*map(float, xyz)))
        surface = Geom_BSplineSurface(poles, u_knots, v_knots, u_mults, v_mults, along.k, across.k)
        faces.append(cq.Face(BRepBuilderAPI_MakeFace(surface, 1e-8).Face()))
    return faces


def section(solid, span):
    cutter = BRepAlgoAPI_Section(solid.wrapped, gp_Pln(gp_Pnt(-span, 0, 0), gp_Dir(1, 0, 0)))
    cutter.Approximation(True)
    cutter.Build()
    return cq.Shape.cast(cutter.Shape())


def assemble(inboard, tip):
    sewing = BRepBuilderAPI_Sewing(SEW_TOLERANCE)
    for face in [*inboard, *tip]:
        sewing.Add(face.wrapped)
    sewing.Perform()
    shell = cq.Shape.cast(sewing.SewedShape())
    if not isinstance(shell, cq.Shell) or sewing.NbFreeEdges() != 0 or sewing.NbMultipleEdges() != 0:
        raise RuntimeError(
            f"Tip did not sew to a single closed shell. "
            f"Free edges {sewing.NbFreeEdges()}, multiple edges {sewing.NbMultipleEdges()}, "
            f"shells {len(shell.Shells())}."
        )
    solid = cq.Solid.makeSolid(shell)
    if not BRepLib.OrientClosedSolid_s(solid.wrapped) or not BRepCheck_Analyzer(solid.wrapped).IsValid():
        raise RuntimeError("The sewn wing is not a valid oriented solid")
    return solid


def volume(solid):
    properties = GProp_GProps()
    BRepGProp.VolumeProperties_s(solid.wrapped, properties, 1e-7, True, False)
    return properties.Mass()


def main():
    original = read_solid(SOURCE)
    skin_face = next(face for face in original.Faces() if face.geomType() == "BSPLINE")
    root = next(face for face in original.Faces() if abs(face.Center().x) < 1e-8)
    cut = Cut(BRep_Tool.Surface_s(skin_face.wrapped))
    print(f"Cut chord from {cut.nose[2]:.6f} to {cut.tail[2]:.6f} mm aft")
    trimmed = cq.Face(BRepBuilderAPI_MakeFace(cut.skin, 0, 1, 0, cut.v, 1e-8).Face())
    tip = boundary_faces(cut)
    solid = assemble([root, trimmed], tip)
    cq.exporters.export(solid, str(OUTPUT))
    rebuilt = read_solid(OUTPUT)
    bounds = rebuilt.BoundingBox()
    if abs(bounds.xlen - SPAN) > 0.001:
        raise RuntimeError(f"Projected half-span {bounds.xlen:.6f} mm")
    print(
        f"Read back one solid with {len(rebuilt.Faces())} faces, "
        f"{bounds.xlen:.6f} mm span and {volume(rebuilt):.2f} mm³ volume "
        f"versus {volume(original):.2f} mm³ before rounding"
    )
    for station in (100.0, 300.0, 500.0, CUT - 0.1):
        prior = section(original, station)
        after = section(rebuilt, station)
        if len(prior.Edges()) == 0 or len(after.Edges()) == 0:
            raise RuntimeError(f"Missing readback seam check at {station:.1f} mm")
        displacement = max(
            cq.Vertex.makeVertex(*edge.positionAt(fraction).toTuple()).distance(prior)
            for edge in after.Edges()
            for fraction in np.linspace(0, 1, 31)
        )
        if displacement > SEW_TOLERANCE:
            raise RuntimeError(f"Inboard profile changed at {station:.1f} mm by {displacement:.6f} mm")
        print(f"Inboard cut at {station:.1f} mm moved at most {displacement:.6f} mm")


if __name__ == "__main__":
    main()
