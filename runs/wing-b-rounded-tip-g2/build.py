"""Build the curvature-matched boundary-surface tip on the saved Wing B STEP."""

from dataclasses import dataclass
from math import comb
from pathlib import Path
from tempfile import TemporaryDirectory

import cadquery as cq
import numpy as np
from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.BRepLib import BRepLib
from OCP.Geom import Geom_BSplineSurface
from OCP.Geom2d import Geom2d_BSplineCurve, Geom2d_Line
from OCP.GProp import GProp_GProps
from OCP.gp import gp_Dir, gp_Dir2d, gp_Pln, gp_Pnt, gp_Pnt2d, gp_Vec
from OCP.TColgp import TColgp_Array1OfPnt2d, TColgp_Array2OfPnt
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
EDGE_FAIRING_FRACTION = 0.18
SKIN_BLEND_FRACTION = 0.62
AFT_REMAP_START = 0.78
AFT_REMAP_STRENGTH = 0.045
SEW_TOLERANCE = 0.0002
AXES = np.array([-1.0, 1.0, -1.0])  # STEP XYZ to span, height, aft


def read_solid(path):
    solids = cq.importers.importStep(str(path)).val().Solids()
    if len(solids) != 1 or not BRepCheck_Analyzer(solids[0].wrapped, True, False, True).IsValid():
        raise ValueError(f"Expected one valid solid in {path}")
    return solids[0]


@dataclass
class Cut:
    skin: object

    def __post_init__(self):
        self.v = self.span_parameter(CUT)
        self.nose_u = self.nose_parameter(self.v)
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
        # Taylor data of the true source outlines in span, including motion of
        # the leading-edge chord extremum in the native surface parameter.
        offsets = np.arange(-2.0, 3.0)
        edge_samples = []
        for offset in offsets:
            v = self.span_parameter(CUT + offset)
            edge_samples.append([self.point(self.nose_parameter(v), v)[1:], self.point(0, v)[1:]])
        self.edge_derivatives = np.array([
            np.polynomial.polynomial.polyfit(offsets, np.array(edge_samples)[:, side], 4)[1:3]
            for side in range(2)
        ])

    def span_parameter(self, span):
        return brentq(lambda v: self.point(0, v)[0] - span, 0, 1, xtol=1e-14)

    def nose_parameter(self, v):
        return brentq(lambda u: self.differential(u, v)[0][2], 0.47, 0.53, xtol=1e-14)

    def point(self, u, v):
        return np.array(self.skin.Value(u, v).Coord()) * AXES

    def differential(self, u, v=None):
        point, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
        self.skin.D1(u, self.v if v is None else v, point, du, dv)
        return np.array(du.Coord()) * AXES, np.array(dv.Coord()) * AXES

    def span_slope(self, u):
        _, dv = self.differential(u)
        return dv[2] / dv[0]

    def profiles(self, chord, span=CUT):
        if span == CUT:
            v, nose_u, nose, tail = self.v, self.nose_u, self.nose, self.tail
        else:
            v = self.span_parameter(span)
            nose_u = self.nose_parameter(v)
            nose = self.point(nose_u, v)
            tail = self.point(0, v)
        if chord <= nose[2]:
            if chord < nose[2] - 1e-6:
                raise ValueError(f"Guide is forward of parent wing at {span:.6f}, {chord:.6f} mm")
            return nose[None, :].repeat(2, axis=0), [nose_u, nose_u]
        if chord >= tail[2]:
            if chord > tail[2] + 1e-6:
                raise ValueError(f"Guide is aft of parent wing at {span:.6f}, {chord:.6f} mm")
            return tail[None, :].repeat(2, axis=0), [1.0, 0.0]

        def chord_error(u):
            return self.point(u, v)[2] - chord

        upper = brentq(chord_error, nose_u, 1, xtol=1e-14)
        lower = brentq(chord_error, 0, nose_u, xtol=1e-14)
        positions = np.array([self.point(u, v) for u in (upper, lower)])
        positions[:, 0] = span
        return positions, [upper, lower]

    def height_jet(self, u, direction):
        """Return height derivatives over the straight projected wrap guide."""
        point = gp_Pnt()
        derivatives = [gp_Vec() for _ in range(5)]
        self.skin.D2(u, self.v, point, *derivatives)
        du, dv, duu, dvv, duv = [np.array(vector.Coord()) * AXES for vector in derivatives]
        jacobian = np.column_stack((du[[0, 2]], dv[[0, 2]]))
        rate_u, rate_v = np.linalg.solve(jacobian, direction)
        acceleration = rate_u**2 * duu + 2 * rate_u * rate_v * duv + rate_v**2 * dvv
        # Cancel projected acceleration because the chosen guide is straight
        # in plan. Its remaining vertical acceleration matches parent curvature.
        accel_u, accel_v = np.linalg.solve(jacobian, -acceleration[[0, 2]])
        h1 = du[1] * rate_u + dv[1] * rate_v
        h2 = acceleration[1] + du[1] * accel_u + dv[1] * accel_v
        return h1, h2


def transition_weight(distance, extent):
    x = min(distance / extent, 1.0)
    return 1 - x**3 * (10 - 15 * x + 6 * x**2)


def continue_outline(cut, u, span, value, axis):
    for side, distance in enumerate((u, 1 - u)):
        if distance < EDGE_FAIRING_FRACTION:
            root = (cut.nose, cut.tail)[side][axis + 1]
            slope, second_order = cut.edge_derivatives[side, :, axis]
            delta = span - CUT
            continued = root + delta * (slope + delta * second_order)
            value += (continued - value) * transition_weight(distance, EDGE_FAIRING_FRACTION)
    return value


def plan_perimeter(cut, u):
    """Control arms lie along the original plan tangents at both ends."""
    arm = PLAN_CONTROL_ARM
    start = cut.nose[[0, 2]]
    stop = cut.tail[[0, 2]]
    nose_slope, tail_slope = cut.span_slopes
    controls = np.array(
        [start, start + [arm, arm * nose_slope], stop + [arm, arm * tail_slope], stop]
    )
    point = sum(
        comb(3, i) * u**i * (1 - u) ** (3 - i) * point
        for i, point in enumerate(controls)
    )
    point[1] = continue_outline(cut, u, point[0], point[1], 1)
    return point


def rim_height_adjustment(cut, u):
    """Retain the first rounded tip's rim elevation as a baseline."""
    nose, tail = cut.rim_slope_correction
    return nose * u * transition_weight(u, 0.25) + tail * (u - 1) * transition_weight(1 - u, 0.35)


def reference_cut_chord(cut, u):
    return cut.nose[2] + (cut.tail[2] - cut.nose[2]) * (1 - np.cos(np.pi * u / 2))


def cut_chord(cut, u):
    # Let the aft guides approach zero cut thickness quadratically, spreading
    # the turnover aft of the joint while retaining the same perimeter.
    distance = 1 - u
    forward_fraction = np.cos(np.pi * u / 2)
    fade = transition_weight(distance, 1 - AFT_REMAP_START)
    scaled = forward_fraction * distance / (distance + AFT_REMAP_STRENGTH * fade)
    return cut.tail[2] - (cut.tail[2] - cut.nose[2]) * scaled


def wrap_guide(cut, u, t):
    """Two branches of one wrap from upper cut through perimeter to lower cut."""
    chord = cut_chord(cut, u)
    cut_points, parameters = cut.profiles(chord)
    if u in (0.0, 1.0):
        points = cut_points[:, None, :].repeat(len(t), axis=1)
        zero = np.zeros((2, 3))
        return points, zero, zero, zero, zero

    origin = cut_points[0, [0, 2]]
    direction = plan_perimeter(cut, u) - origin
    span, chord_shift = direction
    cut_slopes, cut_curvatures = np.array([cut.height_jet(p, direction) for p in parameters]).T
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
    rim_fairing = rim_height_adjustment(cut, u)
    rim_bow = BOW_COEFFICIENT * ((CUT + span)**2 - CUT**2)
    reference_points, _ = cut.profiles(reference_cut_chord(cut, u))
    reference_height = reference_points[:, 1].mean()
    baseline_rim_height = reference_height + rim_bow + rim_fairing
    rim_fairing += reference_height - height
    rim_fairing += continue_outline(cut, u, CUT + span, baseline_rim_height, 0) - baseline_rim_height
    center_slope_correction = cut_slopes.mean() - bow_slope
    center = height + bow + center_slope_correction * fairing + rim_fairing * extension**2

    # The terminal radius follows the local lip height and eases down near the TE.
    # This form gives finite curvature at the outer perimeter and zero height there.
    tip_thinning = 0.35
    radius = 0.7 * half_thickness * (1 - 0.3 * u**4)
    lip_height = half_thickness * np.exp(-tip_thinning)
    rounding = lip_height**2 / (2 * np.linalg.norm(direction) * radius)
    relative_thickness_slope = (cut_slopes[0] - cut_slopes[1]) / (2 * half_thickness)
    thickness_slope = relative_thickness_slope + rounding / 2
    thickness = (
        half_thickness
        * np.exp(thickness_slope * fairing - tip_thinning * extension**2)
        * remaining
        / np.sqrt(remaining**2 + rounding * (1 - remaining**2))
    )

    side = np.array([1, -1])
    heights = center + side[:, None] * thickness
    # Continue the actual parent skin in the same plan coordinates and fade its
    # height into the rolled guide. Blending position instead of extrapolating
    # curvature avoids amplifying small variations in the parent section fit.
    for index in np.flatnonzero(extension < SKIN_BLEND_FRACTION):
        s, c = plan[index]
        parent, _ = cut.profiles(c, s)
        weight = transition_weight(extension[index], SKIN_BLEND_FRACTION)
        heights[:, index] += (parent[:, 1] - heights[:, index]) * weight
    points = np.array([np.column_stack((plan[:, 0], y, plan[:, 1])) for y in heights])
    cut_tangent = 2 * np.array([[span, slope, chord_shift] for slope in cut_slopes])
    cut_acceleration = -cut_tangent
    cut_acceleration[:, 1] += 4 * cut_curvatures
    rim_speed = lip_height / np.sqrt(rounding)
    rim_tangent = np.array([[0, -rim_speed, 0], [0, rim_speed, 0]])
    rim_slope = 2 * (BOW_COEFFICIENT * (CUT + span) * span + rim_fairing)
    rim_acceleration = -2 * np.array([[span, rim_slope, chord_shift]]).repeat(2, axis=0)
    return points, cut_tangent, cut_acceleration, rim_tangent, rim_acceleration


def knot_arrays(spline):
    knots, mults = np.unique(spline.t, return_counts=True)
    occt_knots = TColStd_Array1OfReal(1, len(knots))
    occt_mults = TColStd_Array1OfInteger(1, len(knots))
    for index, (knot, mult) in enumerate(zip(knots, mults), start=1):
        occt_knots.SetValue(index, float(knot))
        occt_mults.SetValue(index, int(mult))
    return occt_knots, occt_mults


def cut_samples(cut):
    # Fit the edge more tightly than OCCT's reconstructed STEP edge tolerance.
    # Refine where the native profile changes fastest, especially around the nose.
    u = (1 - np.cos(np.linspace(0, np.pi, 601))) / 2
    u = np.sort(np.concatenate((u[u >= 0.008], np.linspace(0, 0.008, 401), np.linspace(0.008, 0.025, 81)[1:])))
    values = {float(p): cut.profiles(cut_chord(cut, p))[0] for p in u}
    for _ in range(8):
        spline = make_interp_spline(u, np.array([values[p] for p in u]))
        tests = (u[:-1, None] + np.diff(u)[:, None] * np.array([0.21, 0.5, 0.79])).ravel()
        exact = np.array([cut.profiles(cut_chord(cut, p))[0] for p in tests])
        errors = np.linalg.norm(spline(tests) - exact, axis=2).max(axis=1).reshape(-1, 3).max(axis=1)
        refine = (errors > 1.8e-8) & (np.diff(u) > 1e-5)
        if not refine.any():
            break
        midpoint = tests[1::3][refine]
        values.update(zip(midpoint, exact[1::3][refine]))
        u = np.sort(np.concatenate((u, midpoint)))
    return u


def boundary_faces(cut):
    # Sample the continuous wrap guides in their natural boundary parameters.
    # Denser sampling near the endpoints keeps the cut and perimeter fits within sewing tolerance.
    u = cut_samples(cut)
    t = 1 - (1 - np.linspace(0, 1, 69))**2
    guides, cut_tangents, cut_accels, rim_tangents, rim_accels = [
        np.array(a) for a in zip(*(wrap_guide(cut, p, t) for p in u))
    ]

    faces = []
    for branch in range(2):
        across = make_interp_spline(
            t,
            guides[:, branch].transpose(1, 0, 2),
            axis=0,
            k=5,
            bc_type=(
                [(1, cut_tangents[:, branch]), (2, cut_accels[:, branch])],
                [(1, rim_tangents[:, branch]), (2, rim_accels[:, branch])],
            ),
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


def set_cut_parameter_curves(shell, cut):
    """Give each cap seam the common 3D edge's native parameterization."""
    builder = BRep_Builder()
    parent = next(face for face in shell.Faces() if -CUT < face.Center().x < -1)
    native_cut = cut.skin.VIso(cut.v)
    native_pcurve = Geom2d_Line(gp_Pnt2d(0, cut.v), gp_Dir2d(1, 0))
    remap_threshold = cut_chord(cut, AFT_REMAP_START)
    for face in shell.Faces():
        if face.Center().x >= -CUT:
            continue
        for edge in face.Edges():
            if edge.BoundingBox().xlen > SEW_TOLERANCE:
                continue
            curve = BRepAdaptor_Curve(edge.wrapped)
            first, last = curve.FirstParameter(), curve.LastParameter()
            if abs(first) < 1e-7 and abs(last - 1) < 1e-7:
                # Sewing can choose the cap spline as the shared 3D edge.
                # Use the exact native isocurve and its parent-side pcurve.
                if curve.Value(first).Z() < curve.Value(last).Z():
                    first, last = 0.0, cut.nose_u
                else:
                    first, last = cut.nose_u, 1.0
                builder.UpdateEdge(edge.wrapped, native_cut, 5e-6)
                builder.Range(edge.wrapped, first, last)
                builder.UpdateEdge(edge.wrapped, native_pcurve, parent.wrapped, 5e-6)
                builder.Range(edge.wrapped, parent.wrapped, first, last)
                for vertex in edge.Vertices():
                    point = BRep_Tool.Pnt_s(vertex.wrapped)
                    bound = min((first, last), key=lambda p: point.Distance(native_cut.Value(p)))
                    builder.UpdateVertex(vertex.wrapped, bound, edge.wrapped, 5e-6)
                curve = BRepAdaptor_Curve(edge.wrapped)
            phase = np.linspace(0, np.pi, 2401)
            parameters = first + (last - first) * (1 - np.cos(phase)) / 2
            chord = np.array([-curve.Value(parameter).Z() for parameter in parameters])
            fraction = np.clip((chord - cut.nose[2]) / (cut.tail[2] - cut.nose[2]), 0, 1)
            u = (4 / np.pi) * np.arcsin(np.sqrt(fraction / 2))
            for index in np.flatnonzero((chord > remap_threshold) & (chord < cut.tail[2])):
                u[index] = brentq(lambda p: cut_chord(cut, p) - chord[index], AFT_REMAP_START, 1, xtol=1e-14)
            u[0], u[-1] = (0.0, 1.0) if u[0] < u[-1] else (1.0, 0.0)
            spline = make_interp_spline(parameters, u)
            knots, mults = knot_arrays(spline)
            poles = TColgp_Array1OfPnt2d(1, len(spline.c))
            for index, value in enumerate(spline.c, start=1):
                poles.SetValue(index, gp_Pnt2d(float(value), 0))
            pcurve = Geom2d_BSplineCurve(poles, knots, mults, spline.k)
            builder.UpdateEdge(edge.wrapped, pcurve, face.wrapped, 5e-6)
            builder.Range(edge.wrapped, face.wrapped, first, last)


def assemble(inboard, tip, cut):
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
    set_cut_parameter_curves(shell, cut)
    solid = cq.Solid.makeSolid(shell)
    if not BRepLib.OrientClosedSolid_s(solid.wrapped) or not BRepCheck_Analyzer(
        solid.wrapped, True, False, True
    ).IsValid():
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
    solid = assemble([root, trimmed], tip, cut)
    with TemporaryDirectory(dir=HERE, prefix="step-") as folder:
        pending = Path(folder) / OUTPUT.name
        cq.exporters.export(solid, str(pending))
        rebuilt = read_solid(pending)
        verify_readback(original, rebuilt)
        pending.replace(OUTPUT)


def verify_readback(original, rebuilt):
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
