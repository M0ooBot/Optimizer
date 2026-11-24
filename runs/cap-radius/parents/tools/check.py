"""Write check.json for every parent: reimport validity, fidelity, trailing edge and Converter agreement."""

import hashlib
import json
import sys

import cadquery as cq
import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.BRepCheck import BRepCheck_Analyzer, BRepCheck_NoError, BRepCheck_Shell
from OCP.GeomAbs import GeomAbs_BSplineSurface, GeomAbs_Plane
from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCP.gp import gp_Dir, gp_Pln, gp_Pnt
from scipy.spatial import cKDTree

from build_analytic import precise_volume
from common import PARENTS, REPO, SHAPES, VARIANTS, converter, load_parent

ACCEPTED = {
    "flat": REPO / "runs/cap-trials/retest/flat/wing.step",
    "twisted": REPO / "runs/cap-trials/retest/twisted/wing.step",
    "swept55": REPO / "runs/cap-trials/retest/swept55/wing.step",
    "fin40": REPO / "runs/cap-trials/retest/fin/40/wing.step",
    "fin5": REPO / "runs/cap-trials/retest/fin/5/wing.step",
    "wingb": REPO / "runs/wing-b-rounded-tip-g2/wing-b-rounded.step",
}
CONVERTER_REVISION = "d77bbdf45a8c0aa92efbc2a3d3007cfe4f904105"


def read(path):
    shape = cq.importers.importStep(str(path)).val()
    return shape


def surface_type(face):
    return BRepAdaptor_Surface(face.wrapped).GetType()


def project(points, surfaces, margin=1.0):
    """Smallest distance from each point to any of the given Geom surfaces.

    Each surface is only tried for points inside its bounding box grown by
    margin; ShapeAnalysis_Surface finds the nearest parameters (grid start
    plus Newton), with a global GeomAPI projection as the fallback.
    """
    from OCP.Bnd import Bnd_Box
    from OCP.BndLib import BndLib_AddSurface
    from OCP.GeomAdaptor import GeomAdaptor_Surface
    from OCP.ShapeAnalysis import ShapeAnalysis_Surface
    best = np.full(len(points), np.inf)
    for surface in surfaces:
        box = Bnd_Box()
        BndLib_AddSurface.Add_s(GeomAdaptor_Surface(surface), 0.0, box)
        lo = np.array(box.CornerMin().Coord()) - margin
        hi = np.array(box.CornerMax().Coord()) + margin
        inside = np.flatnonzero(np.all((points >= lo) & (points <= hi), axis=1))
        analysis = ShapeAnalysis_Surface(surface)
        for i in inside:
            p = gp_Pnt(*map(float, points[i]))
            uv = analysis.ValueOfUV(p, 1e-10)
            d = analysis.Gap()
            if not np.isfinite(d) or d > 1e-3:
                projector = GeomAPI_ProjectPointOnSurf(p, surface)
                if projector.NbPoints():
                    d = min(d, projector.LowerDistance())
            best[i] = min(best[i], d)
    return best


def retained_faces(shape, solid):
    """Accepted-model faces that belong to the parent (inside the attachment)."""
    keep = []
    for face in solid.Faces():
        if surface_type(face) != GeomAbs_BSplineSurface:
            continue
        box = face.BoundingBox()
        if shape in ("flat", "twisted", "swept55"):
            limit = {"flat": 350.0, "twisted": 440.0, "swept55": 390.0}[shape]
            ok = max(abs(box.ymin), abs(box.ymax)) <= limit + 1e-3
        elif shape.startswith("fin"):
            ok = box.zmax <= 240.0 + 1e-3
        else:
            ok = box.xmin >= -540.0 - 1e-3
        if ok:
            keep.append(BRep_Tool.Surface_s(face.wrapped))
    return keep


def sample_grid(module, n_s=41, n_y=21):
    lo, hi = module.SPAN_LIMITS
    if module.SHAPE == "wingb":
        lo, hi = 0.5, 539.5
    s = (1 - np.cos(np.linspace(0, np.pi, n_s))) / 2
    y = np.linspace(lo, hi, n_y)
    return np.meshgrid(s, y, indexing="ij")


def exact_skin_points(module):
    S, Y = sample_grid(module)
    if module.SHAPE == "wingb":
        pts = []
        for side in (+1, -1):
            for y in np.unique(Y):
                a, b = module.side_range(y, side)
                u = float(a) + (float(b) - float(a)) * np.unique(S)
                pts.append(module.skin(u, np.full_like(u, y)))
        return np.concatenate(pts)
    return np.concatenate([module.skin(S, Y, side).reshape(-1, 3) for side in (+1, -1)])


def step_skin_points(solid, n=41):
    """Points on every curved face of the STEP over its own parameter grid."""
    pts = []
    for face in solid.Faces():
        if surface_type(face) != GeomAbs_BSplineSurface:
            continue
        adaptor = BRepAdaptor_Surface(face.wrapped)
        u = np.linspace(adaptor.FirstUParameter(), adaptor.LastUParameter(), n)
        v = np.linspace(adaptor.FirstVParameter(), adaptor.LastVParameter(), 21)
        for a in u:
            for b in v:
                pts.append(adaptor.Value(a, b).Coord())
    return np.array(pts)


def validity(path):
    shape = read(path)
    solids = shape.Solids()
    solid = solids[0]
    shells = solid.Shells()
    return solid, {
        "solids": len(solids),
        "shells": len(shells),
        "faces": len(solid.Faces()),
        "face_types": sorted(str(surface_type(f)).split(".")[-1] for f in solid.Faces()),
        "brepcheck_valid": bool(BRepCheck_Analyzer(solid.wrapped).IsValid()),
        "brepcheck_strict_valid": bool(BRepCheck_Analyzer(solid.wrapped, True, False, True).IsValid()),
        "shell_closed": all(BRepCheck_Shell(s.wrapped).Closed() == BRepCheck_NoError for s in shells),
        "volume_mm3": precise_volume(solid),
        "bounding_box_mm": [round(v, 9) for v in (lambda b: (b.xmin, b.xmax, b.ymin, b.ymax, b.zmin, b.zmax))(solid.BoundingBox())],
    }


def section_edges(solid, origin, normal):
    plane = gp_Pln(gp_Pnt(*map(float, origin)), gp_Dir(*map(float, normal)))
    cutter = BRepAlgoAPI_Section(solid.wrapped, plane)
    cutter.Build()
    return cq.Shape.cast(cutter.Shape()).Edges()


def te_thickness_on_step(module, solid, station):
    """Length of the straight trailing-edge segment of the STEP section at a station."""
    frame = module.section_frame(station)
    normal = np.cross(frame["chord_direction"], frame["up_direction"])
    origin = frame["origin"]
    best = None
    for e in section_edges(solid, origin, normal):
        a, b = np.array(e.startPoint().toTuple()), np.array(e.endPoint().toTuple())
        xa = (a - origin) @ frame["chord_direction"] / frame["chord"]
        xb = (b - origin) @ frame["chord_direction"] / frame["chord"]
        if min(xa, xb) > 0.95 and e.geomType() == "LINE" or (min(xa, xb) > 0.95 and abs(xa - xb) < 1e-6):
            length = float(np.linalg.norm(a - b))
            if best is None or abs(length - 0.4) < abs(best - 0.4):
                best = length
    return best


def te_on_boundary_face(module, solid, station):
    """Trailing-edge edge length on a planar root/attachment face."""
    frame = module.section_frame(station)
    normal = np.cross(frame["chord_direction"], frame["up_direction"])
    for face in solid.Faces():
        if surface_type(face) != GeomAbs_Plane:
            continue
        c = np.array(face.Center().toTuple())
        if abs((c - frame["origin"]) @ normal) > 1e-6:
            continue
        lengths = []
        for e in face.Edges():
            a, b = np.array(e.startPoint().toTuple()), np.array(e.endPoint().toTuple())
            xa = (a - frame["origin"]) @ frame["chord_direction"] / frame["chord"]
            xb = (b - frame["origin"]) @ frame["chord_direction"] / frame["chord"]
            chord = float(np.linalg.norm(a - b))
            if min(xa, xb) > 0.95 and chord > 0 and abs(e.Length() - chord) < 1e-7:
                lengths.append(chord)
        return min(lengths) if lengths else None
    return None


def to_2d(points, frame):
    d = points - frame["origin"]
    return np.stack([d @ frame["chord_direction"], d @ frame["up_direction"]], -1)


def nearest_on_section(module, span, query):
    """Distance from 2D points to the finished section of an analytic parent (Newton refined)."""
    frame = module.section_frame(span)
    s = np.unique(np.r_[(1 - np.cos(np.linspace(0, np.pi, 40001))) / 2])
    best = np.full(len(query), np.inf)
    for side in (+1, -1):
        dense = to_2d(module.skin(s, np.full_like(s, span), side), frame)
        _, index = cKDTree(dense).query(query)
        t = s[index]
        for _ in range(8):
            j = module.skin_jet(t, np.full_like(t, span), side)
            p = to_2d(j["point"], frame)
            d1 = np.stack([j["d_chord"] @ frame["chord_direction"], j["d_chord"] @ frame["up_direction"]], -1)
            d2 = np.stack([j["d_chord_chord"] @ frame["chord_direction"], j["d_chord_chord"] @ frame["up_direction"]], -1)
            r = p - query
            g = (r * d1).sum(-1)
            h = (d1 * d1).sum(-1) + (r * d2).sum(-1)
            t = np.clip(t - g / h, 0, 1)
        p = to_2d(module.skin(t, np.full_like(t, span), side), frame)
        best = np.minimum(best, np.linalg.norm(p - query, axis=1))
    return best


def converter_agreement(shape, span):
    """Run the Converter on the sharp section and compare with the analytic blunt section."""
    geometry = converter()
    sharp, blunt = load_parent(shape, "sharp"), load_parent(shape, "blunt")
    frame = sharp.section_frame(span)
    # The Converter's thickest-station search is quadratic in the point count.
    s = (1 - np.cos(np.linspace(0, np.pi, 2501))) / 2
    upper = to_2d(sharp.skin(s, np.full_like(s, span), +1), frame)
    lower = to_2d(sharp.skin(s, np.full_like(s, span), -1), frame)
    loop = np.concatenate([upper[::-1], lower[1:]])
    out = np.array(geometry.blunt_trailing_edge([tuple(p) for p in loop], 0.4, keep_chord=True))
    scale_converter = out[:, 1].max() / loop[:, 1].max()
    cut = {k: float(np.ravel(v)[0]) for k, v in blunt.cut(np.array([span])).items()}
    corners = to_2d(np.array([blunt.skin(1.0, span, side) for side in (+1, -1)]), frame)
    distance = nearest_on_section(blunt, span, out)
    return {
        "station_mm": span,
        "sampled_points": len(loop),
        "converter_points": len(out),
        "converter_te_thickness_mm": float(abs(out[0, 1] - out[-1, 1])),
        "converter_te_chord_position_mm": float(out[0, 0]),
        "analytic_te_chord_position_mm": float(corners[:, 0].mean()),
        "converter_scale": float(scale_converter),
        "analytic_scale": cut["scale"],
        "scale_difference": float(abs(scale_converter - cut["scale"])),
        "converter_cut_chord_fraction": float(out[0, 0] / scale_converter / frame["chord"]),
        "analytic_cut_chord_fraction": cut["cut_chord_fraction"],
        "max_distance_converter_points_to_analytic_section_mm": float(distance.max()),
    }


def wingb_converter_agreement(solid):
    """Converter output at loft stations and at the attachment against the blunt Wing B STEP skin."""
    sys.path.insert(0, str(REPO / "src"))
    from planeopt.phase0_geometry import _profile, load_definition, section_frame
    geometry = converter()
    wing = load_definition(REPO / "runs/paired-wing-trial/a/bowed/definition.json").wing
    profile = _profile(REPO / "aircraft/phase0_sd7037")
    surfaces = [BRep_Tool.Surface_s(f.wrapped) for f in solid.Faces() if surface_type(f) == GeomAbs_BSplineSurface]
    rows = []
    for fraction in (0.0, 0.25, 0.5, 0.75, 0.875, 0.9):
        frame = section_frame(wing, fraction)
        flat = geometry.blunt_trailing_edge([(x * frame.chord_mm, y * frame.chord_mm) for x, y in profile], 0.4, keep_chord=True)
        le, ca, ua = (np.array(a) for a in (frame.le, frame.chord_axis, frame.up_axis))
        points = np.array([le + ca * x + ua * y for x, y in flat])
        d = project(points, surfaces)
        rows.append({"span_fraction": fraction, "projected_span_mm": 600 * fraction,
                     "is_loft_station": bool(abs(fraction * 48 - round(fraction * 48)) < 1e-9),
                     "converter_points": len(points),
                     "converter_te_thickness_mm": float(abs(flat[0][1] - flat[-1][1])),
                     "max_distance_to_step_mm": float(d.max())})
    return rows


def segment_distance(point, a, b):
    t = np.clip((point - a) @ (b - a) / ((b - a) @ (b - a)), 0, 1) if np.any(a != b) else 0.0
    return float(np.linalg.norm(point - (a + t * (b - a))))


def outline_endpoints(module):
    """Gap between the locked outline's ends and the attachment's leading and trailing edges."""
    att = module.ATTACHMENT
    if module.SHAPE == "wingb":
        ends = module.tip_outline(np.array([0.0, 1.0]))
        le = module.skin(module.nose_u(att), att)[[0, 2]]
        if module.TE_THICKNESS > 0:
            corners = module.te_face(np.array([0.0, 1.0]), att)[:, [0, 2]]
        else:
            corners = np.repeat(module.skin(0.0, att)[None, [0, 2]], 2, axis=0)
        return {"plane": "plan (STEP X, Z)", "leading_edge_gap_mm": float(np.linalg.norm(ends[0] - le)),
                "trailing_edge_gap_mm": segment_distance(ends[1], *corners)}
    if module.SHAPE.startswith("fin"):
        start, end = module.tip_outline(0.0, -1), module.tip_outline(0.0, +1)
    else:
        start, end = module.tip_outline(0.0), module.tip_outline(1.0)
    le = module.skin(0.0, att, +1)
    corners = [module.skin(1.0, att, side) for side in (+1, -1)]
    return {"plane": "3D", "leading_edge_gap_mm": float(np.linalg.norm(start - le)),
            "trailing_edge_gap_mm": segment_distance(end, *corners),
            "note": "the trailing-edge gap is the distance to the attachment's trailing edge (point or blunt segment)"}


def gprop_volume(solid, tolerance=None):
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    props = GProp_GProps()
    if tolerance is None:
        BRepGProp.VolumeProperties_s(solid.wrapped, props)
    else:
        BRepGProp.VolumeProperties_s(solid.wrapped, props, tolerance, True, False)
    return props.Mass()


def face_area(face):
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    props = GProp_GProps()
    BRepGProp.SurfaceProperties_s(face.wrapped, props, 1e-9, True)
    return props.Mass()


def step_section_area(solid, origin, normal):
    wire = cq.Wire.assembleEdges(section_edges(solid, origin, normal))
    return face_area(cq.Face.makeFromWires(wire))


def shoelace(points):
    x, y = np.asarray(points).T
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def polyline_thickness(points, x):
    ys = []
    for (x0, y0), (x1, y1) in zip(points, list(points[1:]) + [points[0]]):
        if (x0 <= x < x1) or (x1 <= x < x0):
            ys.append(y0 + (y1 - y0) * (x - x0) / (x1 - x0))
    return max(ys) - min(ys)


def section_area_check(shape, variant, module, solid):
    """Cross-section area of the reimported STEP against the section it was built from.

    Analytic parents: shoelace area of 40,001 parent.py points around the
    section. Wing B: area of the cq spline face through the section points
    (Converter output for blunt, the Converter-resampled SD7037 for sharp) at
    the definition frame, i.e. exactly what each loft station is built from.
    Wing B also compares the thickness distribution with the Converter's points.
    """
    rows = []
    if shape == "wingb":
        sys.path.insert(0, str(REPO / "src"))
        from planeopt.phase0_evidence import _profile_wire
        from planeopt.phase0_geometry import _profile, load_definition, section_frame
        geometry = converter()
        wing = load_definition(REPO / "runs/paired-wing-trial/a/bowed/definition.json").wing
        profile = _profile(REPO / "aircraft/phase0_sd7037")
        fractions = np.array([0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95])
        for station in (1.0, 100.0, 300.0, 500.0, 539.0):
            frame = section_frame(wing, station / 600)
            c = frame.chord_mm
            flat = [(x * c, y * c) for x, y in profile]
            if variant == "blunt":
                flat = geometry.blunt_trailing_edge(flat, 0.4, keep_chord=True)
            le, ca, ua = (np.array(a) for a in (frame.le, frame.chord_axis, frame.up_axis))
            points = [tuple(le + ca * x + ua * y) for x, y in flat]
            reference = face_area(cq.Face.makeFromWires(_profile_wire(cq, points)))
            area = step_section_area(solid, (-station, 0, 0), (1, 0, 0))
            converter_thickness = np.array([polyline_thickness(flat, f * c) for f in fractions])
            parent_thickness = module.section_thickness(fractions, station)
            rows.append({"station_mm": station, "step_area_mm2": area, "reference_area_mm2": reference,
                         "relative_difference": float(area / reference - 1),
                         "chord_fractions": fractions.tolist(),
                         "step_thickness_mm": parent_thickness.tolist(),
                         "converter_points_thickness_mm": converter_thickness.tolist(),
                         "max_thickness_difference_mm": float(np.abs(parent_thickness - converter_thickness).max())})
        note = ("reference = spline face through the loft-station points; converter thickness is read off the "
                "Converter's point polyline, so it differs from the smooth section by up to a few micrometres near the nose")
    else:
        lo, hi = module.SPAN_LIMITS
        span = hi - lo
        stations = [lo + 1.0, lo + 0.25 * span, lo + 0.5 * span + (1.0 if module.FULL_SPAN else 0.0),
                    lo + 0.75 * span, hi - 1.0]
        s = (1 - np.cos(np.linspace(0, np.pi, 20001))) / 2
        for station in stations:
            frame = module.section_frame(station)
            upper = to_2d(module.skin(s, np.full_like(s, station), +1), frame)
            lower = to_2d(module.skin(s, np.full_like(s, station), -1), frame)
            reference = shoelace(np.concatenate([upper[::-1], lower[1:]]))
            normal = np.cross(frame["chord_direction"], frame["up_direction"])
            area = step_section_area(solid, frame["origin"], normal)
            rows.append({"station_mm": station, "step_area_mm2": area, "reference_area_mm2": float(reference),
                         "relative_difference": float(area / reference - 1)})
        note = "reference = shoelace area of 40,001 parent.py points around the section"
    return {"note": note, "stations": rows,
            "max_abs_relative_difference": max(abs(r["relative_difference"]) for r in rows)}


def derivative_check(module):
    """Largest relative mismatch between the analytic jets and central differences."""
    rng = np.random.default_rng(1)
    lo, hi = module.SPAN_LIMITS
    s = rng.uniform(0.02, 0.98, 40)
    y = rng.uniform(lo + 2, hi - 2, 40)
    worst = {}
    sides = (+1, -1)
    if module.SHAPE == "wingb":
        s, y, sides = s[:12], y[:12], (None,)
    for side in sides:
        jet = (lambda a, b: module.skin_jet(a, b)) if side is None else (lambda a, b: module.skin_jet(a, b, side))
        j = jet(s, y)
        h, k = 1e-6, 1e-3
        p = lambda a, b: jet(a, b)["point"]
        q = lambda a, b, key: jet(a, b)[key]
        fd = {
            "d_chord": (p(s + h, y) - p(s - h, y)) / (2 * h),
            "d_span": (p(s, y + k) - p(s, y - k)) / (2 * k),
            "d_chord_chord": (q(s + h, y, "d_chord") - q(s - h, y, "d_chord")) / (2 * h),
            "d_chord_span": (q(s, y + k, "d_chord") - q(s, y - k, "d_chord")) / (2 * k),
            "d_span_span": (q(s, y + k, "d_span") - q(s, y - k, "d_span")) / (2 * k),
        }
        for key, value in fd.items():
            scale = np.abs(j[key]).max() + 1e-12
            worst[key] = max(worst.get(key, 0.0), float(np.abs(j[key] - value).max() / scale))
    return worst


def check(shape, variant):
    folder = PARENTS / shape / variant
    module = load_parent(shape, variant)
    step = folder / "parent.step"
    solid, report = validity(step)
    report = {"step_sha256": hashlib.sha256(step.read_bytes()).hexdigest(), "step_bytes": step.stat().st_size,
              "reimport": report}
    step_surfaces = [BRep_Tool.Surface_s(f.wrapped) for f in solid.Faces() if surface_type(f) == GeomAbs_BSplineSurface]
    exact = exact_skin_points(module)
    if shape != "wingb":
        report["parent_py_to_parent_step_max_mm"] = float(project(exact, step_surfaces).max())
    report["derivative_check_max_relative_error"] = derivative_check(module)
    report["tip_outline_endpoints"] = outline_endpoints(module)
    report["section_areas"] = section_area_check(shape, variant, module, solid)
    report["volume_note"] = {
        "precise_mm3": report["reimport"]["volume_mm3"],
        "occt_default_mm3": gprop_volume(solid),
        "explanation": "volume_mm3 uses BRepGProp with a 1e-9 tolerance; OCCT's default integration (cadquery "
                       "Solid.Volume()) can be badly wrong on these B-spline skins (Wing B sharp: 17% low)",
    }
    att = module.ATTACHMENT
    if variant == "sharp":
        accepted = read(ACCEPTED[shape]).Solids()[0]
        kept = retained_faces(shape, accepted)
        # Identify the accepted whole model by hash only; README.md maps it to its file.
        report["accepted_model_sha256"] = hashlib.sha256(ACCEPTED[shape].read_bytes()).hexdigest()
        report["accepted_retained_faces"] = len(kept)
        report["max_deviation_parent_py_to_accepted_retained_skin_mm"] = float(project(exact, kept).max())
        on_step = step_skin_points(solid)
        report["max_deviation_parent_step_to_accepted_retained_skin_mm"] = float(project(on_step, kept).max())
        report["sample_points"] = {"parent_py": len(exact), "parent_step": len(on_step)}
    else:
        if shape == "wingb":
            stations = [0.0, 135.0, 270.0, 405.0, 540.0]
        elif module.FULL_SPAN:
            stations = [-att, -att / 2, 0.0, att / 2, att]
        else:
            stations = [0.0, att / 4, att / 2, 3 * att / 4, att]
        rows = []
        for station in stations:
            boundary = station in (module.SPAN_LIMITS[0], module.SPAN_LIMITS[1]) and not (
                module.FULL_SPAN and station == 0.0)
            measured = (te_on_boundary_face(module, solid, station) if boundary
                        else te_thickness_on_step(module, solid, station))
            if shape == "wingb":
                ends = module.te_face(np.array([0.0, 1.0]), station)
            else:
                ends = np.array([module.skin(1.0, station, side) for side in (-1, +1)])
            rows.append({"station_mm": station, "measured_on": "planar end face edge" if boundary else "plane section of the STEP",
                         "step_te_thickness_mm": measured,
                         "parent_py_te_thickness_mm": float(np.linalg.norm(ends[1] - ends[0]))})
        report["te_thickness"] = rows
        if shape == "wingb":
            report["converter_agreement"] = {"revision": CONVERTER_REVISION, "sections": wingb_converter_agreement(solid)}
        else:
            report["converter_agreement"] = {"revision": CONVERTER_REVISION,
                                             "root": converter_agreement(shape, 0.0),
                                             "attachment": converter_agreement(shape, att)}
    if shape == "wingb":
        report["note"] = ("parent.py evaluates the B-spline faces of parent.step, so parent.py and parent.step "
                          "are the same surface by construction")
    return report


def main():
    targets = sys.argv[1:] or SHAPES
    for shape in targets:
        for variant in VARIANTS:
            data = check(shape, variant)
            (PARENTS / shape / variant / "check.json").write_text(json.dumps(data, indent=1) + "\n")
            print(shape, variant, json.dumps({k: v for k, v in data.items() if k not in ("step_sha256",)})[:1500], flush=True)


if __name__ == "__main__":
    main()
