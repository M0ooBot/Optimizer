"""Build the Wing B parents (sharp and blunt), their parent.py modules and build reports.

sharp: the skin surface of the accepted parent STEP
       runs/paired-wing-trial/a/measured/01-wing-occt.step, cut exactly at the
       540 mm attachment by B-spline segmentation (no refit).
blunt: the same Phase 0 sections, regenerated through the Airfoil Converter's
       blunt_trailing_edge(points, 0.4, keep_chord=True) at every loft station,
       closed with a straight trailing-edge line and lofted with the same
       cq.Solid.makeLoft call as the original (planeopt.phase0_evidence).
"""

import json
import pprint
import sys

import cadquery as cq
import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace
from OCP.Geom import Geom_BSplineCurve, Geom_BSplineSurface
from OCP.GeomAbs import GeomAbs_BSplineSurface
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.TColgp import TColgp_Array1OfPnt
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCP.gp import gp_Pnt, gp_Vec
from scipy.interpolate import make_lsq_spline

from build_analytic import precise_volume, solid_from, SEW_TOLERANCE
from common import PARENTS, REPO, converter

sys.path.insert(0, str(REPO / "src"))
from planeopt.phase0_evidence import _profile_wire, _stations  # noqa: E402
from planeopt.phase0_geometry import _profile, load_definition, section_frame  # noqa: E402

SOURCE = REPO / "runs/paired-wing-trial/a/measured/01-wing-occt.step"
DEFINITION = REPO / "runs/paired-wing-trial/a/bowed/definition.json"
FIXTURE = REPO / "aircraft/phase0_sd7037"
CUT = 540.0
TE = 0.4
AXES = np.array([-1.0, 1.0, -1.0])  # STEP XYZ to (span, height, aft), as in the accepted cap build
PLAN_CONTROL_ARM = (600.0 - CUT) * 4 / 3
EDGE_FAIRING_FRACTION = 0.18
# Common-knot representation of each Converter section spline. The raw splines have
# different knot vectors (the cut drops a different number of points per station),
# which ThruSections would merge into a ~76 MB surface; a least-squares refit on one
# shared knot vector stays within this tolerance of each raw spline.
REFIT_TOLERANCE = 2e-5  # mm


def read_solid(path):
    solids = cq.importers.importStep(str(path)).val().Solids()
    assert len(solids) == 1
    return solids[0]


def bspline_faces(solid):
    faces = [f for f in solid.Faces() if BRepAdaptor_Surface(f.wrapped).GetType() == GeomAbs_BSplineSurface]
    faces.sort(key=lambda f: -f.Area())
    return faces


def as_bspline(face):
    return BRepAdaptor_Surface(face.wrapped).BSpline()


def bisect(function, lo, hi, iterations=100):
    f_lo = function(lo)
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        f_mid = function(mid)
        if np.sign(f_mid) == np.sign(f_lo):
            lo, f_lo = mid, f_mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def v_at_span(surface, span, u=0.0):
    return bisect(lambda v: -surface.Value(u, v).X() - span, 0.0, 1.0)


def segment(surface, v_end):
    copy = surface.Copy()  # OCP returns the derived Geom_BSplineSurface
    u0, u1, _, _ = copy.Bounds()
    copy.Segment(u0, u1, 0.0, v_end)
    return copy


def edge(curve):
    return cq.Edge(BRepBuilderAPI_MakeEdge(curve).Edge())


def planar(edges):
    wire = cq.Wire.assembleEdges(edges)
    return cq.Face(BRepBuilderAPI_MakeFace(wire.wrapped, True).Face())


def x_spread(surface, v):
    us = np.linspace(*surface.Bounds()[:2], 41)
    xs = [surface.Value(float(u), v).X() for u in us]
    return float(np.ptp(xs))


def outline_data(original_skin):
    """Freeze the plan perimeter inputs exactly as runs/wing-b-rounded-tip-g2/build.py derives them."""
    def point(u, v):
        return np.array(original_skin.Value(u, v).Coord()) * AXES

    def d1(u, v):
        p, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
        original_skin.D1(u, v, p, du, dv)
        return np.array(du.Coord()) * AXES, np.array(dv.Coord()) * AXES

    span_parameter = lambda span: bisect(lambda v: point(0, v)[0] - span, 0.0, 1.0)
    nose_parameter = lambda v: bisect(lambda u: d1(u, v)[0][2], 0.47, 0.53)
    v = span_parameter(CUT)
    nose_u = nose_parameter(v)
    nose, tail = point(nose_u, v), point(0.0, v)
    nose[0] = tail[0] = CUT
    slopes = [d1(u, v)[1][2] / d1(u, v)[1][0] for u in (nose_u, 0.0)]
    offsets = np.arange(-2.0, 3.0)
    samples = []
    for offset in offsets:
        w = span_parameter(CUT + offset)
        samples.append([point(nose_parameter(w), w)[1:], point(0, w)[1:]])
    samples = np.array(samples)
    taylor = [np.polynomial.polynomial.polyfit(offsets, samples[:, side], 4)[1:3, 1] for side in range(2)]
    return {
        "nose_span_aft": [float(nose[0]), float(nose[2])],
        "tail_span_aft": [float(tail[0]), float(tail[2])],
        "nose_plan_slope": float(slopes[0]),
        "tail_plan_slope": float(slopes[1]),
        "control_arm_mm": PLAN_CONTROL_ARM,
        "edge_fairing_fraction": EDGE_FAIRING_FRACTION,
        "edge_aft_taylor": [[float(a), float(b)] for a, b in taylor],
        "source": "runs/paired-wing-trial/a/measured/01-wing-occt.step at 540 mm (sharp parent)",
    }


def converter_sections():
    """Phase 0 loft stations blunted through the Converter: [(fraction, 3D points, 2D points)]."""
    geometry = converter()
    wing = load_definition(DEFINITION).wing
    profile = _profile(FIXTURE)
    sections = []
    for fraction in _stations(wing, 4):
        frame = section_frame(wing, float(fraction))
        chord = frame.chord_mm
        flat = geometry.blunt_trailing_edge([(x * chord, y * chord) for x, y in profile], TE, keep_chord=True)
        le, ca, ua = (np.array(a) for a in (frame.le, frame.chord_axis, frame.up_axis))
        points = [tuple(le + ca * x + ua * y) for x, y in flat]
        sections.append((float(fraction), points, flat))
    return sections


def common_knots(curves):
    """Normalised knots dense enough for every section: the first section's knots, halved."""
    first = curves[0]
    t0, t1 = first.FirstParameter(), first.LastParameter()
    knots = np.array([(first.Knot(i) - t0) / (t1 - t0) for i in range(1, first.NbKnots() + 1)])
    knots = np.unique(np.r_[knots, (knots[1:] + knots[:-1]) / 2])
    return np.r_[[0.0] * 3, knots, [1.0] * 3]


def refit(curve, knots):
    """Least-squares cubic on the common knots with exact end points; returns (Geom curve, max error)."""
    t0, t1 = curve.FirstParameter(), curve.LastParameter()
    own = np.array([(curve.Knot(i) - t0) / (t1 - t0) for i in range(1, curve.NbKnots() + 1)])
    dense = np.unique(np.concatenate([np.linspace(a, b, 24) for a, b in zip(own, own[1:])]))
    dense = np.unique(np.r_[dense, knots[3:-3]])
    value = lambda ts: np.array([curve.Value(t0 + t * (t1 - t0)).Coord() for t in ts])
    points = value(dense)
    weights = np.ones(len(dense))
    weights[[0, -1]] = 1e6
    fit = make_lsq_spline(dense, points, knots, k=3, w=weights)
    coefficients = fit.c.copy()
    coefficients[0], coefficients[-1] = points[0], points[-1]
    fit.c = coefficients
    check = np.unique(np.concatenate([np.linspace(a, b, 37) for a, b in zip(own, own[1:])]))
    error = float(np.linalg.norm(fit(check) - value(check), axis=1).max())
    poles = TColgp_Array1OfPnt(1, len(coefficients))
    for i, c in enumerate(coefficients, 1):
        poles.SetValue(i, gp_Pnt(*map(float, c)))
    values, counts = np.unique(knots, return_counts=True)
    k = TColStd_Array1OfReal(1, len(values))
    m = TColStd_Array1OfInteger(1, len(values))
    for i, (a, b) in enumerate(zip(values, counts), 1):
        k.SetValue(i, float(a))
        m.SetValue(i, int(b))
    return Geom_BSplineCurve(poles, k, m, 3), error


def build_sharp():
    original = read_solid(SOURCE)
    skin = as_bspline(bspline_faces(original)[0])
    v = v_at_span(skin, CUT)
    cut = segment(skin, v)
    faces = [cq.Face(BRepBuilderAPI_MakeFace(cut, SEW_TOLERANCE).Face()),
             planar([edge(cut.VIso(0.0))]), planar([edge(cut.VIso(v))])]
    solid = solid_from(faces)
    return solid, skin, {"v_attachment_in_source": v, "attachment_plane_x_spread_mm": x_spread(skin, v),
                         "source_skin_poles": [skin.NbUPoles(), skin.NbVPoles()],
                         "source_skin_degrees": [skin.UDegree(), skin.VDegree()]}


def build_blunt():
    sections = converter_sections()
    exact = [BRep_Tool.Curve_s(cq.Edge.makeSpline([cq.Vector(*p) for p in pts]).wrapped, 0, 1)
             for _, pts, _ in sections]
    knots = common_knots(exact)
    wires, errors = [], []
    for (fraction, pts, flat), curve in zip(sections, exact):
        fitted, error = refit(curve, knots)
        errors.append(error)
        start, end = cq.Vector(*pts[0]), cq.Vector(*pts[-1])
        wires.append(cq.Wire.assembleEdges([edge(fitted), cq.Edge.makeLine(end, start)]))
    loft = cq.Solid.makeLoft(wires)
    skin_face, te_face = bspline_faces(loft)[:2]
    skin, te = as_bspline(skin_face), as_bspline(te_face)
    v = v_at_span(skin, CUT)
    cut_skin, cut_te = segment(skin, v), segment(te, v)
    faces = [cq.Face(BRepBuilderAPI_MakeFace(cut_skin, SEW_TOLERANCE).Face()),
             cq.Face(BRepBuilderAPI_MakeFace(cut_te, SEW_TOLERANCE).Face())]
    for station in (0.0, v):
        faces.append(planar([edge(cut_skin.VIso(station)), edge(cut_te.VIso(station))]))
    solid = solid_from(faces)
    stations_json = [{"fraction": f, "converter_points": len(flat),
                      "te_gap_mm": float(abs(flat[0][1] - flat[-1][1])), "te_x_mm": float(flat[0][0])}
                     for f, _, flat in sections]
    if max(errors) > REFIT_TOLERANCE:
        raise RuntimeError(f"Common-knot refit error {max(errors):.3g} mm exceeds {REFIT_TOLERANCE} mm")
    return solid, skin, {"v_attachment_in_loft": v, "attachment_plane_x_spread_mm": x_spread(skin, v),
                         "loft_skin_poles": [skin.NbUPoles(), skin.NbVPoles()],
                         "loft_skin_degrees": [skin.UDegree(), skin.VDegree()],
                         "loft_stations": len(sections),
                         "common_knot_refit_max_error_mm": max(errors),
                         "refit_tolerance_mm": REFIT_TOLERANCE,
                         "converter_stations": stations_json}


def surface_data(skin_surface, v_attachment):
    """Orientation facts parent.py needs: which u end is upper, nose bracket, TE face direction."""
    u0, u1 = skin_surface.Bounds()[:2]
    mid_v = v_attachment / 2
    zs = np.array([skin_surface.Value(float(u), mid_v).Z() for u in np.linspace(u0, u1, 2001)])
    nose = np.linspace(u0, u1, 2001)[zs.argmax()]
    quarter = lambda end: skin_surface.Value(float(nose + 0.5 * (end - nose)), mid_v).Y()
    upper_end = u1 if quarter(u1) > quarter(u0) else u0
    lower_end = u0 if upper_end == u1 else u1
    return {"u_reference": float(0.5 * (u0 + u1)), "v_attachment": float(v_attachment),
            "nose_bracket": [float(nose - 0.08), float(nose + 0.08)],
            "nose_sign": 1.0 if upper_end > lower_end else -1.0,
            "u_end_upper": float(upper_end), "u_end_lower": float(lower_end)}


def compose(variant, data, outline):
    source = (PARENTS / "tools/shapes/wingb.py").read_text(encoding="utf-8")
    doc_end = source.index('"""', 3) + 3
    doc, body = source[:doc_end], source[doc_end:]
    body = body.replace('VARIANT = "sharp"', f'VARIANT = "{variant}"', 1)
    body = body.replace("TE_THICKNESS = 0.0", f"TE_THICKNESS = {0.4 if variant == 'blunt' else 0.0}", 1)
    body = body.replace("SURFACE_DATA = {}", f"SURFACE_DATA = {pprint.pformat(data, sort_dicts=False)}", 1)
    body = body.replace("OUTLINE_DATA = {}", f"OUTLINE_DATA = {pprint.pformat(outline, sort_dicts=False)}", 1)
    note = ("\nVariant: sharp. The accepted parent surface, segmented at the attachment.\n"
            if variant == "sharp" else
            "\nVariant: blunt. Sections regenerated through the Airfoil Converter with TE thickness\n"
            "0.4 mm and 'Keep chord after the cut', lofted as the original; the trailing edge is a\n"
            "flat (ruled) face 0.4 mm thick.\n")
    text = doc[:-3] + note + "\nGenerated by runs/cap-radius/parents/tools/build_wingb.py; do not edit by hand.\n" \
        + '"""\n\nimport numpy as np\n' + body
    (PARENTS / "wingb" / variant / "parent.py").write_text(text, encoding="utf-8")


def main():
    sharp, source_skin, sharp_report = build_sharp()
    outline = outline_data(source_skin)
    blunt, blunt_skin, blunt_report = build_blunt()
    for variant, solid, report, skin in (("sharp", sharp, sharp_report, source_skin),
                                         ("blunt", blunt, blunt_report, blunt_skin)):
        folder = PARENTS / "wingb" / variant
        folder.mkdir(parents=True, exist_ok=True)
        step = folder / "parent.step"
        cq.exporters.export(solid, str(step))
        v = report.get("v_attachment_in_source", report.get("v_attachment_in_loft"))
        # parent.step holds the segmented surface, whose v already ends at the attachment.
        faces = bspline_faces(read_solid(step))
        reread = faces[0]
        data = surface_data(as_bspline(reread), as_bspline(reread).Bounds()[3])
        if variant == "blunt":
            te = as_bspline(faces[1])
            u0, u1 = te.Bounds()[:2]
            mid = data["v_attachment"] / 2
            data["te_u_range"] = [float(u0), float(u1)]
            data["te_u_upward"] = bool(te.Value(u1, mid).Y() > te.Value(u0, mid).Y())
        compose(variant, data, outline)
        report.update({"faces": len(solid.Faces()), "volume_mm3": precise_volume(solid), "step_bytes": step.stat().st_size,
                       "surface_data": data})
        (folder / "build.json").write_text(json.dumps(report, indent=2) + "\n")
        print(variant, {k: v for k, v in report.items() if k != "converter_stations"})


if __name__ == "__main__":
    main()
