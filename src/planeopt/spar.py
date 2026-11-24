"""Straight round spar search and B-rep check on a half-wing STEP, in millimetres."""
import argparse
import hashlib
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory

import cadquery as cq
import numpy as np
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from scipy.optimize import LinearConstraint, brentq, linprog, minimize
from scipy.spatial import ConvexHull

RADIUS = 3.2
WALL = 1.0
TOL = 1e-7

def sections(solid, stations, count=481):
    faces = [face for face in solid.Faces() if face.geomType() == "BSPLINE"]
    if len(faces) != 1:
        raise RuntimeError("Expected the single closed Wing B loft skin face")
    surface = faces[0]._geomAdaptor()
    a, b, c, d = surface.Bounds()
    u_parameters = np.linspace(a, b, count)
    rings = []
    for span in stations:
        target = min(max(span, -surface.Value(a,c).X()), -surface.Value(a,d).X())
        v = brentq(lambda t: -surface.Value(a, t).X()-target, c, d)
        ring = surface.VIso(v)
        points = np.asarray([ring.Value(u).Coord() for u in u_parameters])
        if np.max(np.abs(points[:, 0]+span)) > 1e-5:
            raise RuntimeError("The skin iso-curve does not have a constant physical span")
        rings.append(points)
    return rings

def support_lp(stations, profiles, bounds, scale):
    rows, rhs = [], []
    for span, profile in zip(stations, profiles):
        for a, b, c in ConvexHull(profile[:, 1:]).equations:
            rows.append((a, b, a*span/scale, b*span/scale, 1))
            rhs.append(-c)
    result = linprog([0,0,0,0,-1], A_ub=np.asarray(rows), b_ub=np.asarray(rhs),
                     bounds=[*bounds, (None, None)], method="highs")
    if not result.success:
        raise RuntimeError(result.message)
    return result.x

def axis(params, scale, last_joint, engagement):
    y, z, dy, dz = params[:4]
    root = np.array([0., y, z])
    direction = np.array([-1., dy/scale, dz/scale])
    direction /= np.linalg.norm(direction)
    radial_span = RADIUS * np.linalg.norm(direction[1:])
    full_start = (last_joint + radial_span) / -direction[0]
    end_length = full_start + engagement
    return root, direction, end_length, radial_span


def cylinder_sdf(points, params, scale, last_joint, engagement):
    root, direction, end_length, _ = axis(params, scale, last_joint, engagement)
    relatives = points - root
    axial = np.sum(relatives * direction, axis=1)
    radial = np.linalg.norm(relatives - axial[:, None] * direction, axis=1) - RADIUS
    upper = axial - end_length
    # The root trim is treated separately. The infinite inboard extension
    # provides conservative distances near the entry for interior points.
    return np.hypot(np.maximum(radial, 0), np.maximum(upper, 0)) + np.minimum(np.maximum(radial, upper), 0)


def search(seed, stations, profiles, bounds, scale, last_joint, engagement):
    points = np.concatenate(profiles)
    guess = seed.copy()
    guess[4] -= RADIUS
    best = guess.copy()
    best[4] = cylinder_sdf(points, guess, scale, last_joint, engagement).min()
    ring_count = len(profiles[0])
    nearest_by_ring = np.argmin(cylinder_sdf(points, guess, scale, last_joint, engagement).reshape(-1, ring_count), axis=1)
    active = np.arange(len(profiles))*ring_count + nearest_by_ring
    angles = np.linspace(0, 2*np.pi, 25)[:-1]
    normals = np.column_stack((np.cos(angles), np.sin(angles)))
    rows, upper = [], []
    center_steps = np.unique(np.r_[np.arange(0,len(stations),10), np.searchsorted(stations, last_joint)])
    for i in center_steps:
        if stations[i] <= scale:
            rows.extend(np.column_stack((normals, normals*stations[i]/scale, np.zeros(len(normals)))))
            upper.extend((profiles[i][:,1:] @ normals.T).max(axis=0))
    interior_corridor = LinearConstraint(np.asarray(rows), -np.inf, upper)
    for _ in range(12):
        values = cylinder_sdf(points, guess, scale, last_joint, engagement)
        active = np.unique(np.r_[active, np.argpartition(values, 31)[:32]])
        candidate = minimize(lambda p: -p[4], guess, method="SLSQP",
                             bounds=[*bounds, (-10, 25)],
                             constraints=[interior_corridor, {"type":"ineq",
                                           "fun":lambda p: cylinder_sdf(points[active], p, scale, last_joint, engagement)-p[4]}],
                             options={"ftol":1e-9, "maxiter":100})
        guess = candidate.x
        actual = cylinder_sdf(points, guess, scale, last_joint, engagement).min()
        if actual > best[4]:
            best = np.r_[guess[:4], actual]
        if guess[4]-actual <= 1e-6:
            break
    # CAD containment separately verifies the axis selected by the sampled corridor.
    return best


def slab(shape, inboard, outboard):
    box = shape.BoundingBox()
    return cq.Solid.makeBox(outboard - inboard, box.ylen + 2, box.zlen + 2,
                            cq.Vector(-outboard, box.ymin - 1, box.zmin - 1))


def extrema(a, b):
    result = BRepExtrema_DistShapeShape(a.wrapped, b.wrapped, 1e-8)
    if not result.IsDone() or result.NbSolution() < 1:
        raise RuntimeError("OCCT could not evaluate the B-rep distance")
    return result


def on_channel(point, params, scale, last_joint, engagement):
    root, direction, end_length, _ = axis(params, scale, last_joint, engagement)
    axial = np.dot(point-root, direction)
    radial = point-root-axial*direction
    radial_length = np.linalg.norm(radial)
    if axial < end_length and axial-end_length < radial_length-RADIUS:
        return root + axial*direction + RADIUS*radial/radial_length
    cap = root + end_length*direction
    if radial_length <= RADIUS:
        return cap + radial
    return cap + RADIUS*radial/radial_length


def refine_clearance(skin_face, params, scale, last_joint, engagement):
    surface = skin_face._geomAdaptor()
    u0, u1, v0, v1 = surface.Bounds()
    uv_bounds = ((u0,u1),(v0,v1))

    def clearance(uv):
        point = np.asarray(surface.Value(*uv).Coord())
        return cylinder_sdf(point[None,:], params, scale, last_joint, engagement)[0]

    u_grid = np.linspace(u0, u1, 181)
    v_grid = np.linspace(v0, v1, 151)
    grid = np.asarray([(u,v) for v in v_grid for u in u_grid])
    points = np.asarray([surface.Value(*uv).Coord() for uv in grid])
    values = cylinder_sdf(points, params, scale, last_joint, engagement)
    minima = []
    for i in range(0, len(grid), 5*len(u_grid)):
        stop = min(i+5*len(u_grid), len(grid))
        start = grid[i+np.argmin(values[i:stop])]
        candidate = minimize(clearance, start, bounds=uv_bounds, method="Nelder-Mead",
                             options={"xatol":1e-11,"fatol":1e-11})
        point = np.asarray(surface.Value(*candidate.x).Coord())
        witness = on_channel(point, params, scale, last_joint, engagement)
        minima.append((float(candidate.fun), point, witness))
    return minima


def verify_and_export(solid, seed, stations, profiles, bounds, source, out, span,
                      joints, engagement, height, requested_tip, scale):
    last_joint = joints[-1]
    params = search(seed, stations, profiles, bounds, scale, last_joint, engagement)
    root, direction, end_length, radial_span = axis(params, scale, last_joint, engagement)
    end = root + end_length*direction
    entry = root + last_joint/-direction[0]*direction
    root_setback = (RADIUS + WALL + 1.) / -direction[0]
    tool = cq.Solid.makeCylinder(RADIUS, end_length+root_setback,
                                 cq.Vector(*(root-root_setback*direction)), cq.Vector(*direction))
    channel = tool.intersect(slab(tool, 0, max(span, -end[0]) + RADIUS + WALL + 1)).Solids()[0]
    faces = [face for face in solid.Faces()
             if not (face.BoundingBox().xlen < 1e-5 and abs(face.Center().x) < 1e-5)]
    skin = cq.Compound.makeCompound(faces)
    dist = extrema(channel, skin)
    cad_wall = dist.Value()
    channel_witness = dist.PointOnShape1(1).Coord()
    skin_witness = dist.PointOnShape2(1).Coord()
    points = np.concatenate(profiles)
    values = cylinder_sdf(points, params, scale, last_joint, engagement)
    worst = int(np.argmin(values))
    continuous = refine_clearance(next(f for f in faces if f.geomType() == "BSPLINE"),
                                  params, scale, last_joint, engagement)
    continuous.sort(key=lambda item: item[0])
    if continuous[0][0] < 0:
        signed_wall, skin_point, channel_point = continuous[0]
        skin_witness, channel_witness = skin_point.tolist(), channel_point.tolist()
    else:
        signed_wall = cad_wall
        for wall, skin_point, channel_point in continuous:
            # The radial projection only describes the root-clipped cylinder if
            # the projected point lies on the outboard side of the root plane.
            if channel_point[0] <= 0 and wall < signed_wall:
                signed_wall = wall
                skin_witness, channel_witness = skin_point.tolist(), channel_point.tolist()
    # An interior witness and a disjoint connected channel establish containment
    # apart from the deliberate opening on the excluded root mounting face.
    center_station = root + direction*(end_length/2)
    center_inside = solid.isInside(cq.Vector(*center_station), TOL)
    # This in-wing expanded cylinder encloses the one-millimetre offset
    # of the clipped channel and checks the wall floor independently of extrema.
    clearance_envelope_contained = None
    channel_farthest_span = -end[0] + radial_span
    if center_inside and cad_wall >= WALL and channel_farthest_span < span:
        expanded = cq.Solid.makeCylinder(RADIUS + WALL, end_length + root_setback + WALL,
                                         cq.Vector(*(root-root_setback*direction)), cq.Vector(*direction))
        expanded = expanded.intersect(slab(expanded, 0, max(span, -end[0]) + RADIUS + WALL + 1))
        clearance_envelope_contained = not expanded.cut(solid).Solids()
    reported_wall = float(signed_wall) if center_inside else None
    passing = bool(center_inside and reported_wall >= WALL and clearance_envelope_contained)
    failure_reason = None
    if channel_farthest_span >= span:
        failure_reason = "Channel extends to or beyond the outer tip face"
    elif not center_inside:
        failure_reason = "Axis midpoint is outside the wing"
    elif reported_wall < WALL:
        failure_reason = "Minimum exterior wall is below the requirement"
    elif not clearance_envelope_contained:
        failure_reason = "Minimum wall could not be certified by expanded channel containment"
    result = {
        "input_step":source.as_posix(), "input_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
        "coordinates":"CAD xyz, mm. Span grows along -X, up +Y, aft -Z.",
        "cut_plane_x_mm":[-j for j in joints],
        "stock_diameter_mm":6., "channel_diameter_mm":2*RADIUS, "radial_installation_allowance_mm":.2,
        "required_wall_mm":WALL, "required_full_diameter_engagement_mm":engagement,
        "tip_piece_length_mm":requested_tip, "max_print_height_mm":height,
        "projected_span_mm":span, "split_span_mm":list(joints),
        "axis_root_mm":root.tolist(), "axis_direction":direction.tolist(),
        "axis_channel_end_mm":end.tolist(), "axis_length_mm":end_length,
        "tool_start_mm":(root-root_setback*direction).tolist(),
        "axis_last_joint_mm":entry.tolist(), "axis_last_joint_to_end_mm":float(np.linalg.norm(end-entry)),
        "joint_obliquity_deg":float(np.degrees(np.arccos(-direction[0]))),
        "radial_span_projection_mm":radial_span,
        "channel_farthest_span_mm":float(channel_farthest_span),
        "full_diameter_engagement_mm":end_length-(last_joint+radial_span)/-direction[0],
        "full_diameter_stock_engagement_mm":end_length-(last_joint+3.*(radial_span/RADIUS))/-direction[0],
        "min_wall_mm":reported_wall,
        "cad_brep_boundary_distance_mm":float(cad_wall),
        "signed_sampled_wall_mm":float(values[worst]),
        "continuous_skin_clearance_mm":float(continuous[0][0]),
        "worst_channel_point_mm":channel_witness, "worst_skin_point_mm":skin_witness,
        "min_wall_minus_requirement_mm":None if reported_wall is None else reported_wall-WALL,
        "center_inside_cad":center_inside, "geometry_pass":passing,
        "clearance_envelope_contained":clearance_envelope_contained,
        "failure_reason":failure_reason,
        "search": {
            "root_yz_bounds_mm":bounds[:2], "yz_delta_at_search_reference_span_mm_bounds_mm":bounds[2:],
            "span_station_range_mm":[float(stations[0]),float(stations[-1])],
            "span_station_count":len(stations), "points_per_section":len(profiles[0]),
            "search_reference_span_mm":scale, "seed_planar_axis_radius_mm":float(seed[4]), "params_yz_root_delta_at_search_reference_span_mm":params[:4].tolist(),
            "method":"Sampled section convex hull linear program, then sequential constrained maximization of signed 3D clearance to the finite cylinder. Continuous skin optimization, B-rep distance and clearance-envelope containment evaluated separately.",
            "signed_wall_when_intersecting":"Negative continuous clearance is a breach witness. CAD boundary distance is zero when the bore touches or intersects the skin."
        }
    }
    if passing:
        bored = solid.cut(tool)
        if len(bored.Solids()) != 1 or not bored.isValid():
            raise RuntimeError("Bored wing is not one valid solid")
        result["print_pieces"] = []
        sections = []
        for index, (start, stop) in enumerate(zip((0., *joints), (*joints, span+1e-5)), 1):
            blank = solid.intersect(slab(solid, start, stop))
            piece = blank.cut(tool)
            if len(piece.Solids()) != 1 or not piece.isValid():
                raise RuntimeError(f"Print section {index} is not one valid solid")
            measured_height = piece.BoundingBox().xlen
            if measured_height > height+1e-5:
                raise RuntimeError(f"Print section {index} exceeds {height} mm")
            filename = f"print-{index:02d}.step"
            cq.exporters.export(piece, str(out/filename))
            check = cq.importers.importStep(str(out/filename)).solids().vals()
            if len(check) != 1 or not check[0].isValid():
                raise RuntimeError(f"STEP export {index} did not reimport as one valid solid")
            cylinders = [face._geomAdaptor() for face in check[0].Faces() if face.geomType() == "CYLINDER"]
            if not cylinders or any(
                abs(cylinder.Radius()-RADIUS) > TOL
                or np.linalg.norm(np.cross(np.asarray(cylinder.Axis().Location().Coord())-root, direction)) > TOL
                or abs(abs(np.dot(cylinder.Axis().Direction().Coord(), direction))-1) > TOL
                for cylinder in cylinders
            ):
                raise RuntimeError(f"Print section {index} did not retain the 6.4 mm bore")
            sections.append(piece)
            result["print_pieces"].append({
                "file":filename, "span_mm":[start, min(stop,span)], "span_print_height_mm":measured_height,
                "reimported_cylindrical_bore_faces":len(cylinders)
            })
        assembly_path = out / "print-assembly.step"
        cq.exporters.export(cq.Compound.makeCompound(sections), str(assembly_path))
        assembly = cq.importers.importStep(str(assembly_path)).solids().vals()
        if len(assembly) != len(sections) or not all(piece.isValid() for piece in assembly):
            raise RuntimeError("Print assembly did not reimport with all valid pieces")
    return result

def solve_spar(step: Path, out: Path, *, tip_piece_length_mm: float = 250.,
               min_tip_engagement_mm: float = 20., max_print_height_mm: float = 250.,
               split_span_mm: tuple[float, ...] | None = None) -> dict:
    """Search a Wing B loft STEP and export bored print pieces for a passing fit."""
    step, out = Path(step), Path(out)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError(f"Choose a fresh or empty output directory: {out}")
    tip, engagement, height = (tip_piece_length_mm, min_tip_engagement_mm, max_print_height_mm)
    if not all(math.isfinite(x) and x > 0 for x in (tip, engagement, height)):
        raise ValueError("Tip piece length, engagement, and maximum print height must be positive")
    solids = cq.importers.importStep(str(step)).solids().vals()
    if len(solids) != 1 or not solids[0].isValid():
        raise ValueError("Expected one valid input solid")
    solid = solids[0]
    box = solid.BoundingBox()
    caps = sorted(face.Center().x for face in solid.Faces() if face.BoundingBox().xlen < 1e-5)
    if len(caps) != 2 or abs(caps[1]) > 1e-5:
        raise ValueError("Expected planar end faces with the root at X=0")
    span = -caps[0]
    if split_span_mm is None:
        if not 0 < tip < span:
            raise ValueError("Tip piece length must be shorter than the projected half span")
        remaining = span - tip
        inboard_piece_count = math.ceil(remaining / height - 1e-10)
        joints = tuple(remaining*i/inboard_piece_count for i in range(1, inboard_piece_count+1))
    else:
        joints = tuple(float(s) for s in split_span_mm)
        if not joints or not all(math.isfinite(j) for j in joints):
            raise ValueError("Supply one or more finite split spans")
        tip = span - joints[-1]
    lengths = np.diff([0., *joints, span])
    if np.any(lengths <= TOL) or np.any(lengths > height+1e-5):
        raise ValueError("Split spans must be increasing and every piece must fit the maximum print height")
    last_joint = joints[-1]
    scale = last_joint + engagement
    # Bounding-box differences cap the possible axis slope because each
    # candidate must be contained between the root and the final joint.
    bounds = [(box.ymin, box.ymax), (box.zmin, box.zmax),
              (-box.ylen*scale/last_joint, box.ylen*scale/last_joint),
              (-box.zlen*scale/last_joint, box.zlen*scale/last_joint)]
    # A skin point far outboard of the cap cannot control a practical bore.
    # Include the rest of the span if the required engagement approaches the tip.
    stop = min(span, scale+20)
    stations = np.unique(np.r_[np.arange(0., stop, 1.), stop, *joints])
    profiles = sections(solid, stations)
    seed = support_lp(stations[stations <= scale],
                      [profile for station, profile in zip(stations, profiles) if station <= scale],
                      bounds, scale)
    out.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".spar-", dir=out.parent) as temporary:
        stage = Path(temporary) / "result"
        stage.mkdir()
        np.savez(stage / "search-sections.npz", stations=stations,
                 points=np.concatenate(profiles), counts=[len(p) for p in profiles])
        result = verify_and_export(solid, seed, stations, profiles, bounds, step, stage, span,
                                   joints, engagement, height, tip, scale)
        (stage / "result.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
        if out.exists():
            out.rmdir()
        stage.rename(out)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--tip-length", type=float, default=250.)
    parser.add_argument("--engagement", type=float, default=20.)
    parser.add_argument("--max-height", type=float, default=250.)
    parser.add_argument("--split-span", nargs="+", type=float)
    args = parser.parse_args()
    solve_spar(args.step, args.out, tip_piece_length_mm=args.tip_length,
               min_tip_engagement_mm=args.engagement, max_print_height_mm=args.max_height,
               split_span_mm=args.split_span)
