"""Independently measure the exported Wing B STEP surfaces in millimetres."""

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import cadquery as cq
import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepCheck import BRepCheck_Analyzer, BRepCheck_NoError
from OCP.BRepClass import BRepClass_FaceClassifier
from OCP.BRepTools import BRepTools
from OCP.Bnd import Bnd_Box
from OCP.TopAbs import TopAbs_IN, TopAbs_ON, TopAbs_REVERSED
from OCP.gp import gp_Dir, gp_Pln, gp_Pnt, gp_Pnt2d, gp_Vec
from scipy.optimize import brentq, minimize_scalar, root


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "paired-wing-trial/a/measured/01-wing-occt.step"
CUT = 540.0
AXES = np.array([-1.0, 1.0, -1.0])  # World XYZ to span, height, aft
SPAN = np.array([1.0, 0.0, 0.0])
AFT = np.array([0.0, 0.0, 1.0])
CHORD_FRACTIONS = np.unique(np.r_[np.geomspace(1e-8, 1e-4, 17), .0005, .001, .005, .01, .02, .03,
                                  np.linspace(.04, .96, 47), .97, .98, .99, .995, .999, .9995,
                                  1 - np.geomspace(1e-6, 1e-4, 9)])
TRACE_FRACTIONS = (.001, .005, .01, .02, .05, .15, .25, .35, .5, .65, .75,
                   .85, .95, .98, .99, .995, .999)
TRACE_OFFSETS = (.01, .02, .05, .1, .2, .35, .5, .75, 1, 1.5, 2, 2.5, 3, 4,
                 5, 6, 7, 8, 9, 10, 12, 15, 18, 20, 25, 30, 40, 50, 55)


def norm(x):
    return float(np.linalg.norm(x))


def unit(x):
    length = norm(x)
    if length < 1e-13:
        raise ValueError("Undefined tangent direction")
    return x / length


def to_design(vector):
    return np.array(vector.Coord()) * AXES


def angle_deg(a, b):
    return float(np.degrees(np.arctan2(norm(np.cross(a, b)), a @ b)))


def box_of(shape):
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape.wrapped, box, False, False)
    return box.Get()


def digest_of(path):
    with path.open("rb") as file:
        return hashlib.file_digest(file, "sha256").hexdigest()


@dataclass
class Jet:
    position: np.ndarray
    first: np.ndarray
    second: np.ndarray
    normal: np.ndarray

    def frame(self, direction=SPAN):
        tangent = unit(direction - self.normal * (self.normal @ direction))
        return np.column_stack((tangent, np.cross(self.normal, tangent)))

    def shape_operator(self, frame):
        # The pseudoinverse converts unit physical tangent vectors to UV speeds.
        # The second fundamental form is invariant to UV scale and skew.
        speeds = np.linalg.lstsq(self.first, frame, rcond=None)[0]
        ii = self.second @ self.normal
        return speeds.T @ ii @ speeds

    def fixed_aft_normal_curvature(self):
        tangent = unit(np.cross(AFT, self.normal))
        return float(self.shape_operator(tangent[:, None])[0, 0])


class Surface:
    def __init__(self, face):
        self.face = face.wrapped
        self.adaptor = BRepAdaptor_Surface(self.face)
        self.sign = -1 if self.face.Orientation() == TopAbs_REVERSED else 1
        u0, u1, v0, v1 = BRepTools.UVBounds_s(self.face)
        self.lower, self.upper = np.array([u0, v0]), np.array([u1, v1])

    def point(self, uv):
        return to_design(self.adaptor.Value(*uv))

    def d1(self, uv):
        values = gp_Pnt(), gp_Vec(), gp_Vec()
        self.adaptor.D1(*uv, *values)
        point, du, dv = map(to_design, values)
        return point, np.column_stack((du, dv))

    def jet(self, uv):
        values = (gp_Pnt(), *(gp_Vec() for _ in range(5)))
        self.adaptor.D2(*uv, *values)
        point, du, dv, duu, dvv, duv = map(to_design, values)
        first = np.column_stack((du, dv))
        second = np.array([[duu, duv], [duv, dvv]])
        return Jet(point, first, second, self.sign * unit(np.cross(du, dv)))

    def condition_ratio(self, uv):
        scales = np.linalg.svd(self.d1(uv)[1], compute_uv=False)
        return float(scales[-1] / scales[0]) if scales[0] else 0.0

    def invert_plan(self, target, guess):
        target = np.asarray(target)

        def residual(uv):
            point, tangent = self.d1(np.clip(uv, self.lower, self.upper))
            uv_delta = uv - np.clip(uv, self.lower, self.upper)
            return (point + tangent @ uv_delta)[[0, 2]] - target

        def jacobian(uv):
            return self.d1(np.clip(uv, self.lower, self.upper))[1][[0, 2]]

        fit = root(residual, guess, jac=jacobian, tol=1e-10)
        uv = fit.x
        if norm(residual(uv)) > 1e-7 or np.any(uv < self.lower - 1e-9) or np.any(uv > self.upper + 1e-9):
            return None
        uv = np.clip(uv, self.lower, self.upper)
        where = BRepClass_FaceClassifier(self.face, gp_Pnt2d(*uv), 1e-9).State()
        return uv if where in (TopAbs_IN, TopAbs_ON) else None


class Boundary:
    def __init__(self, edge, surface):
        self.surface = surface
        self.start, self.end = BRep_Tool.Range_s(edge.wrapped, surface.face)
        self.pcurve = BRep_Tool.CurveOnSurface_s(edge.wrapped, surface.face, self.start, self.end)

    def uv(self, parameter):
        return np.array(self.pcurve.Value(parameter).Coord())

    def at_aft(self, aft):
        def error(parameter):
            return self.surface.point(self.uv(parameter))[2] - aft

        # Some lip pcurves turn infinitesimally near their ends. The join curves
        # are monotone in aft away from the leading nose and trailing vertex.
        parameter = brentq(error, self.start, self.end, xtol=1e-14)
        return self.uv(parameter)


def compare(first, second, direction=SPAN):
    frame = first.frame(direction)
    ii_first = first.shape_operator(frame)
    ii_second = second.shape_operator(frame)
    error = ii_second - ii_first
    return {
        "gap_mm": norm(first.position - second.position),
        "tangent_angle_deg": angle_deg(first.normal, second.normal),
        "inboard_ii_per_mm": ii_first.tolist(),
        "outboard_ii_per_mm": ii_second.tolist(),
        "ii_difference_per_mm": error.tolist(),
        "ii_difference_norm_per_mm": float(np.linalg.norm(error, ord=2)),
    }


def edge_curve(edge):
    curve = BRepAdaptor_Curve(edge.wrapped)
    return curve, curve.FirstParameter(), curve.LastParameter()


def sample_join(edge, inner_face, outer_face):
    curve, lo, hi = edge_curve(edge)
    nose, tail = sorted(to_design(curve.Value(v))[2] for v in (lo, hi))
    inside, outside = Surface(inner_face), Surface(outer_face)
    inner_edge, outer_edge = Boundary(edge, inside), Boundary(edge, outside)
    points, traces, corner_queries = [], [], []
    for fraction in CHORD_FRACTIONS:
        aft = nose + fraction * (tail - nose)
        uv_inner, uv_outer = inner_edge.at_aft(aft), outer_edge.at_aft(aft)
        inner_jet, outer_jet = inside.jet(uv_inner), outside.jet(uv_outer)
        points.append(dict(chord_fraction=fraction, aft_mm=aft, **compare(inner_jet, outer_jet)))
    aft_inputs = [nose + fraction * (tail - nose) for fraction in TRACE_FRACTIONS]
    if nose < 150.0 < tail:
        aft_inputs.append(150.0)
    for aft in aft_inputs:
        uv = outer_edge.at_aft(aft)
        inner_jet = inside.jet(inner_edge.at_aft(aft))
        outer_jet = outside.jet(uv)
        trace = []
        offsets = [0, *TRACE_OFFSETS]
        for offset in offsets:
            if offset:
                solved = outside.invert_plan((CUT + offset, aft), uv)
                if solved is None:
                    break
                uv = solved
            jet = outside.jet(uv)
            shape = jet.shape_operator(jet.frame())
            trace.append(dict(span_mm=jet.position[0], height_mm=jet.position[1],
                              normal=jet.normal.tolist(), ii_per_mm=shape.tolist(),
                              fixed_aft_normal_curvature_per_mm=jet.fixed_aft_normal_curvature()))
        traces.append(dict(chord_fraction=(aft - nose) / (tail - nose), aft_mm=aft,
                           join=compare(inner_jet, outer_jet), values=trace,
                           inboard_fixed_aft_normal_curvature_per_mm=inner_jet.fixed_aft_normal_curvature()))
    for label, inner_p, outer_p in zip(("first", "last"), (inner_edge.start, inner_edge.end),
                                     (outer_edge.start, outer_edge.end)):
        corner = dict(edge_end=label, position_mm=inside.point(inner_edge.uv(inner_p)).tolist())
        condition = outside.condition_ratio(outer_edge.uv(outer_p))
        corner["outboard_uv_condition_ratio"] = condition
        if condition < 1e-7:
            corner["derivatives"] = "Degenerate outboard UV at boundary vertex; see interior limits"
            corner_queries.append(corner)
            continue
        try:
            corner.update(compare(inside.jet(inner_edge.uv(inner_p)), outside.jet(outer_edge.uv(outer_p))))
        except (RuntimeError, ValueError) as error:
            corner["derivatives"] = str(error)
        corner_queries.append(corner)
    return dict(side="upper" if inside.jet(inner_edge.at_aft((nose + tail) / 2)).normal[1] > 0 else "lower",
                edge_tolerance_mm=BRep_Tool.Tolerance_s(edge.wrapped), points=points,
                fixed_aft_traces=traces, endpoints=corner_queries)


def vector_curvature(first, second):
    return (second - unit(first) * (unit(first) @ second)) / (first @ first)


def edge_jet(curve, parameter):
    values = gp_Pnt(), gp_Vec(), gp_Vec()
    curve.D2(parameter, *values)
    return tuple(map(to_design, values))


def sample_lip(edge, faces):
    curve, lo, hi = edge_curve(edge)
    boundaries = [Boundary(edge, Surface(face)) for face in faces]
    rows = []
    fractions = np.unique(np.r_[.00001, .0001, .0005, .001, .002, .005, .01, .02,
                                np.linspace(.025, .975, 191), .98, .99, .995, .998,
                                .999, .9995, .9999, .99999])
    for fraction in fractions:
        parameter = lo + fraction * (hi - lo)
        position, first, second = edge_jet(curve, parameter)
        plan_curvature = vector_curvature(first[[0, 2]], second[[0, 2]])
        row = dict(edge_fraction=float(fraction), position_mm=position.tolist(),
                   span_tangent=unit(first).tolist(), plan_curvature_per_mm=norm(plan_curvature),
                   space_curvature_per_mm=norm(vector_curvature(first, second)))
        if len(boundaries) == 2 and BRep_Tool.SameParameter_s(edge.wrapped):
            jets = [boundary.surface.jet(boundary.uv(parameter)) for boundary in boundaries]
            row.update(compare(*jets, direction=first))
        rows.append(row)
    result = minimize_scalar(lambda parameter: curve.Value(parameter).X(), bounds=(lo, hi),
                             method="bounded", options={"xatol": 1e-15})
    wrap_radii = []
    if len(boundaries) == 2 and BRep_Tool.SameParameter_s(edge.wrapped):
        front, back = sorted(to_design(curve.Value(parameter))[2] for parameter in (lo, hi))
        for aft in (90, 110, 135, 150):
            if not front < aft < back:
                continue
            parameter = brentq(lambda p: -curve.Value(p).Z() - aft, lo, hi, xtol=1e-14)
            jets = [boundary.surface.jet(boundary.uv(parameter)) for boundary in boundaries]
            radii = [norm(jet.normal[:2]) / abs(jet.fixed_aft_normal_curvature()) for jet in jets]
            wrap_radii.append(dict(aft_mm=aft, position_mm=to_design(curve.Value(parameter)).tolist(),
                                   fixed_aft_radii_mm=radii))
    silhouette = []
    for span in (540.01, 540.1, 540.25, 540.5, 541, 542, 545, 550, 555, 560, 565,
                 570, 575, 580, 585, 590, 595, 598, 599, 599.5, 599.9, 599.99):
        def error(parameter):
            return -curve.Value(parameter).X() - span

        sections = [edge_jet(curve, brentq(error, start, stop, xtol=1e-14))
                    for start, stop in ((lo, result.x), (result.x, hi))]
        row = {"span_mm": span}
        for label, (position, first, second) in zip(("leading", "trailing"),
                                                    sorted(sections, key=lambda item: item[0][2])):
            row[f"{label}_height_mm"], row[f"{label}_aft_mm"] = position[1:]
            row[f"{label}_height_slope"], row[f"{label}_plan_slope"] = first[1:] / first[0]
            row[f"{label}_height_second_per_mm"], row[f"{label}_plan_second_per_mm"] = (
                (second[1:] * first[0] - first[1:] * second[0]) / first[0]**3
            )
        silhouette.append(row)
    return dict(edge_tolerance_mm=BRep_Tool.Tolerance_s(edge.wrapped),
                max_span_position_mm=to_design(curve.Value(result.x)).tolist(), points=rows,
                front_and_top_outline=silhouette, wrap_radii=wrap_radii)


def section(shape, span):
    plane = gp_Pln(gp_Pnt(-span, 0, 0), gp_Dir(1, 0, 0))
    operation = BRepAlgoAPI_Section(shape.wrapped, plane, False)
    operation.Approximation(True)
    operation.Build()
    result = cq.Shape.cast(operation.Shape())
    if not result.Edges():
        raise ValueError(f"No section at {span:g} mm")
    return result


def extreme_aft(shape, sense):
    candidates = []
    for edge in shape.Edges():
        curve, lo, hi = edge_curve(edge)
        parameters = np.linspace(lo, hi, 61)
        values = np.array([-sense * curve.Value(parameter).Z() for parameter in parameters])
        index = np.argmin(values)
        lower = parameters[max(0, index - 1)]
        upper = parameters[min(len(parameters) - 1, index + 1)]
        fit = minimize_scalar(lambda p: -sense * curve.Value(p).Z(), bounds=(lower, upper),
                              method="bounded", options={"xatol": 1e-14})
        candidates.extend(to_design(curve.Value(p)) for p in (lower, fit.x, upper))
    return min(candidates, key=lambda point: sense * point[2])


def outline(shape):
    stations = (50, 300, 535, 539.9, 540.1, 545, 565, 580, 590, 599, 599.9)
    rows = []
    for span in stations:
        profile = section(shape, span)
        front, back = extreme_aft(profile, 1), extreme_aft(profile, -1)
        _, y_min, _, _, y_max, _ = box_of(profile)
        rows.append(dict(span_mm=span, leading_height_mm=front[1], leading_aft_mm=front[2],
                         trailing_height_mm=back[1], trailing_aft_mm=back[2],
                         section_height_mm=y_max - y_min, chord_mm=back[2] - front[2]))
    return rows


def inboard_distances(shape, original):
    rows = []
    for span in (50, 300, 510, 539, 539.9):
        before, after = section(original, span), section(shape, span)
        distances = [target.distance(cq.Vertex.makeVertex(*edge.positionAt(fraction).toTuple()))
                     for source, target in ((before, after), (after, before))
                     for edge in source.Edges() for fraction in np.linspace(0, 1, 31)]
        rows.append(dict(span_mm=span, sampled_difference_mm=max(distances)))
    return rows


def inboard_support(inner_face, original):
    source_face = next(face for face in original.Faces() if face.geomType() == "BSPLINE")
    source, trimmed = (BRep_Tool.Surface_s(face.wrapped) for face in (source_face, inner_face))

    def knot_vectors(surface, axis):
        if axis == "u":
            return tuple((surface.UKnot(i), surface.UMultiplicity(i)) for i in range(1, surface.NbUKnots() + 1))
        return tuple((surface.VKnot(i), surface.VMultiplicity(i)) for i in range(1, surface.NbVKnots() + 1))

    def counts(surface):
        return surface.NbUPoles(), surface.NbVPoles(), surface.UDegree(), surface.VDegree()

    report = dict(source_uv_bounds=BRepTools.UVBounds_s(source_face.wrapped),
                  trimmed_uv_bounds=BRepTools.UVBounds_s(inner_face.wrapped),
                  same_degree_and_pole_counts=counts(source) == counts(trimmed),
                  same_knots_and_multiplicities=all(knot_vectors(source, axis) == knot_vectors(trimmed, axis)
                                                   for axis in ("u", "v")))
    if report["same_degree_and_pole_counts"]:
        poles = [(i, j) for i in range(1, source.NbUPoles() + 1) for j in range(1, source.NbVPoles() + 1)]
        report["max_pole_distance_mm"] = max(source.Pole(i, j).Distance(trimmed.Pole(i, j)) for i, j in poles)
        report["max_weight_difference"] = max(abs(source.Weight(i, j) - trimmed.Weight(i, j)) for i, j in poles)
    return report


def attached_faces(edge, faces):
    return [face for face in faces if any(edge.wrapped.IsSame(other.wrapped) for other in face.Edges())]


def exact_edge_check(edge, face):
    curve, lo, hi = edge_curve(edge)
    if not BRep_Tool.SameParameter_s(edge.wrapped):
        return {"same_parameter": False}
    boundary = Boundary(edge, Surface(face))

    def deviation(parameter):
        return norm(boundary.surface.point(boundary.uv(parameter)) - to_design(curve.Value(parameter)))

    end_cluster = np.geomspace(1e-9, .009, 100)
    fractions = np.unique(np.r_[end_cluster, np.linspace(0, 1, 3001), 1 - end_cluster])
    largest, at_fraction = max((deviation(lo + fraction * (hi - lo)), float(fraction)) for fraction in fractions)
    return dict(same_parameter=True, max_sampled_deviation_mm=largest,
                edge_parameter_fraction=at_fraction, edge_tolerance_mm=BRep_Tool.Tolerance_s(edge.wrapped))


def review(path, *, with_outline=True):
    digest = digest_of(path)
    shape = cq.importers.importStep(str(path)).val()
    if digest_of(path) != digest:
        raise RuntimeError(f"STEP changed during import: {path}")
    original = cq.importers.importStep(str(SOURCE)).val()
    faces = shape.Faces()
    inner_faces = [face for face in faces if box_of(face)[0] > -CUT - .001 and face.geomType() != "PLANE"]
    outer_faces = [face for face in faces if box_of(face)[0] < -CUT - .001]
    joins, lips = [], []
    for edge in shape.Edges():
        interior, exterior = attached_faces(edge, inner_faces), attached_faces(edge, outer_faces)
        edge_box = box_of(edge)
        if len(interior) == len(exterior) == 1 and abs(edge_box[0] + CUT) < .001 and abs(edge_box[3] + CUT) < .001:
            row = sample_join(edge, *interior, *exterior)
            row["outboard_pcurve_check"] = exact_edge_check(edge, exterior[0])
            joins.append(row)
        if len(exterior) == 2 and edge_box[0] < -CUT - .001:
            lips.append(sample_lip(edge, exterior))

    analyzer = BRepCheck_Analyzer(shape.wrapped, True, False, True)
    exact_face_status = [
        [str(status).split(".")[-1] for status in analyzer.Result(face.wrapped).Status() if status != BRepCheck_NoError]
        for face in faces
    ]
    report = dict(model=path.name, sha256=digest, solid_count=len(shape.Solids()), face_count=len(faces),
                  valid=shape.isValid(), exact_valid=analyzer.IsValid(), exact_face_status=exact_face_status,
                  closed_shells=[shell.Closed() for shell in shape.Shells()], bounds_world_xyz_mm=box_of(shape),
                  joins=joins, lips=lips)
    if len(inner_faces) == 1 and inner_faces[0].geomType() == "BSPLINE":
        report["inboard_bspline_support"] = inboard_support(inner_faces[0], original)
    preserved = inner_faces[0] if len(inner_faces) == 1 else shape
    report["inboard_section_distances"] = inboard_distances(preserved, original)
    if with_outline:
        report["outline"] = outline(shape)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", nargs="?", type=Path, default=HERE / "wing-b-rounded.step")
    parser.add_argument("--output", type=Path, default=HERE / "g2-readback.json")
    parser.add_argument("--skip-outline", action="store_true", help="Skip the full-solid section sweep")
    args = parser.parse_args()
    report = review(args.step, with_outline=not args.skip_outline)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"{args.step.name}: {report['solid_count']} solid, {report['face_count']} faces, "
          f"valid={report['valid']}, exact_valid={report['exact_valid']}")
    for branch in report["joins"]:
        interior = [row for row in branch["points"] if .001 <= row["chord_fraction"] <= .999]
        print(branch["side"], "0.1–99.9% chord", {
            key: max(row[key] for row in interior)
            for key in ("gap_mm", "tangent_angle_deg", "ii_difference_norm_per_mm")
        })
    print(f"{len(report['joins'])} join edges and {len(report['lips'])} lip edges")


if __name__ == "__main__":
    main()
