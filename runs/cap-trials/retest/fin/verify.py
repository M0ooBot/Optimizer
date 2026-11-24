import hashlib
import json
from itertools import product
from pathlib import Path

import cadquery as cq
import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.gp import gp_Dir, gp_Pln, gp_Pnt, gp_Vec
from scipy.optimize import brentq
from scipy.signal import find_peaks

from build import ATTACHMENT, A_SLOPE, M_SLOPE

HERE = Path(__file__).resolve().parent
CAMERAS = [(0.15, 1, 0.25), (0.35, 1, 0.65), (-0.55, 1, -0.25), (0, 1, 0)]
ANGLES = np.deg2rad([0, 27, 53, 90, 123, 158])
FREQUENCY = 7 * np.pi
FRACTIONS = [-.99, -.97, -.9, -.7, -.35, 0, .35, .7, .9, .97, .99]


def arr(vec):
    return np.array([vec.X(), vec.Y(), vec.Z()])


def unit(vec):
    return vec / np.linalg.norm(vec)


def surface_geometry(surface, u, v, side):
    point = gp_Pnt()
    du, dv, duu, dvv, duv = [gp_Vec() for _ in range(5)]
    surface.D2(u, v, point, du, dv, duu, dvv, duv)
    du, dv, duu, dvv, duv = [arr(vector) for vector in [du, dv, duu, dvv, duv]]
    cross = np.cross(du, dv)
    normal = side * unit(cross)
    normal_derivatives = side * (np.eye(3) - np.outer(normal, normal)) @ np.column_stack([
        np.cross(duu, dv) + np.cross(du, duv),
        np.cross(duv, dv) + np.cross(du, dvv)
    ]) / np.linalg.norm(cross)
    tangents = np.column_stack([du, dv])
    inverse = np.linalg.pinv(tangents)
    return arr(point), normal, normal_derivatives @ inverse, inverse


def reflection(geometry, view, stripe):
    _, normal, normal_jacobian, _ = geometry
    projection = np.dot(view, normal)
    phase = FREQUENCY * (2 * projection * np.dot(normal, stripe) - np.dot(view, stripe))
    phase_gradient = 2 * FREQUENCY * (
        projection * stripe + np.dot(stripe, normal) * view
    ) @ normal_jacobian
    return phase, phase_gradient


def angle(a, b):
    return np.degrees(np.arctan2(np.linalg.norm(np.cross(unit(a), unit(b))), np.dot(unit(a), unit(b))))


def max_record(records, value, **location):
    if value > records.get("value", -1):
        records.clear()
        records.update(value=float(value), **location)


def check_joins(surfaces):
    max_values = {name: {} for name in (
        "gap_mm", "normal_angle_deg", "shape_operator_gap_per_mm",
        "zebra_phase_rad", "zebra_gradient_rad_per_mm", "zebra_crossing_direction_deg"
    )}
    grid = np.linspace(0, 1, 361)
    crossing_count = 0
    for side in (-1, 1):
        parent, cap = surfaces[side]
        geom_parent = [surface_geometry(parent, u, 1, side) for u in grid]
        geom_cap = [surface_geometry(cap, u, 0, side) for u in grid]
        for u, first, second in zip(grid, geom_parent, geom_cap):
            location = dict(side=side, u=float(u), point_mm=first[0].tolist())
            max_record(max_values["gap_mm"], np.linalg.norm(first[0]-second[0]), **location)
            max_record(max_values["normal_angle_deg"], angle(first[1], second[1]), **location)
            max_record(max_values["shape_operator_gap_per_mm"], np.linalg.norm(first[2]-second[2], ord=2),
                       **location)
        for vi, ai in product(range(len(CAMERAS)), range(len(ANGLES))):
            view = unit(np.array(CAMERAS[vi]) * [1, side, 1])
            right = unit(np.cross(view, [0, 0, 1]))
            up = unit(np.cross(right, view))
            stripe = np.cos(ANGLES[ai]) * right + np.sin(ANGLES[ai]) * up
            phase_first = []
            for u, first, second in zip(grid, geom_parent, geom_cap):
                p1, g1 = reflection(first, view, stripe)
                p2, g2 = reflection(second, view, stripe)
                phase_first.append(p1)
                location = dict(side=side, u=float(u), camera=vi, stripe_angle_deg=float(np.degrees(ANGLES[ai])),
                                point_mm=first[0].tolist())
                max_record(max_values["zebra_phase_rad"], abs(p1-p2), **location)
                max_record(max_values["zebra_gradient_rad_per_mm"], np.linalg.norm(g1-g2), **location)
            levels = np.arange(-7, 7) * np.pi + np.pi/2
            crossings = set()
            for low, high, p_low, p_high in zip(grid[:-1], grid[1:], phase_first[:-1], phase_first[1:]):
                for level in levels[(levels >= min(p_low, p_high)) & (levels <= max(p_low, p_high))]:
                    u = brentq(lambda u: reflection(surface_geometry(parent, u, 1, side), view, stripe)[0] - level,
                               low, high, xtol=5e-15)
                    crossings.add(round(u, 13))
            for u in sorted(crossings):
                first = surface_geometry(parent, u, 1, side)
                second = surface_geometry(cap, u, 0, side)
                _, g1 = reflection(first, view, stripe)
                _, g2 = reflection(second, view, stripe)
                if min(np.linalg.norm(g1), np.linalg.norm(g2)) < 1e-10:
                    continue
                crossing_count += 1
                max_record(max_values["zebra_crossing_direction_deg"],
                           angle(np.cross(first[1], g1), np.cross(second[1], g2)), side=side, u=u,
                           camera=vi, stripe_angle_deg=float(np.degrees(ANGLES[ai])), point_mm=first[0].tolist())
    return dict(maxima=max_values, crossing_count=crossing_count, camera_vectors=CAMERAS,
                stripe_angles_deg=np.degrees(ANGLES).tolist(), stripe_wave_number=float(FREQUENCY),
                samples_per_side=len(grid), note="Angles and derivatives compared in common physical frame.")


def perimeter(surfaces, reach):
    max_values = {name: {} for name in (
        "locked_outline_error_mm", "side_gap_mm", "normal_angle_deg", "shape_operator_gap_per_mm")}
    for u in (0.0, 1.0):
        chord_sign = 1 if u == 0 else -1
        for v in np.linspace(0, 0.9999, 501):
            first = surface_geometry(surfaces[1][1], u, v, 1)
            second = surface_geometry(surfaces[-1][1], u, v, -1)
            point = first[0]
            z = point[2]
            h = (z-ATTACHMENT)/reach
            expected_x = 80 + M_SLOPE*z + chord_sign * (80 + A_SLOPE*z) * np.sqrt(max(0,1-h**3))
            location = dict(u=u, v=float(v), point_mm=point.tolist())
            max_record(max_values["locked_outline_error_mm"], np.hypot(point[0]-expected_x, point[1]), **location)
            max_record(max_values["side_gap_mm"], np.linalg.norm(first[0]-second[0]), **location)
            max_record(max_values["normal_angle_deg"], angle(first[1], second[1]), **location)
            max_record(max_values["shape_operator_gap_per_mm"], np.linalg.norm(first[2]-second[2], ord=2),
                       **location)
    tip = arr(surfaces[1][1].Value(.5, 1))
    return dict(maxima=max_values, crown_mm=tip.tolist(), samples_per_branch=501,
                excluded="Crown is a collapsed UV boundary. Ordinary surface derivatives are undefined exactly there.")


def section(solid, surfaces, f):
    slope = M_SLOPE + f * A_SLOPE
    offset = 80 * (1 + f)
    normal = unit(np.array([1, 0, -slope]))
    cutting_plane = gp_Pln(gp_Pnt(offset, 0, 0), gp_Dir(*normal))
    operation = BRepAlgoAPI_Section(solid.wrapped, cutting_plane, False)
    operation.Approximation(True)
    operation.Build()
    assert operation.IsDone()
    edges = cq.Shape.cast(operation.Shape()).Edges()
    cap_edges = [edge for edge in edges if max(edge.startPoint().z, edge.endPoint().z) > ATTACHMENT + 1e-5]
    assert len(cap_edges) == 2, (f, [(ed.startPoint(), ed.endPoint()) for ed in edges])
    halves = []
    for side in (-1, 1):
        cap = surfaces[side][1]

        def residual(u, v):
            p = cap.Value(u, v)
            return p.X() - slope * p.Z() - offset

        if f == 0:
            lip_v = 1.0
        else:
            lip_u = float(f < 0)
            lip_v = brentq(lambda v: residual(lip_u, v), 0, 1, xtol=1e-14)
        signed = []
        xyz = []
        for i, s in enumerate(np.linspace(0, 1, 2001)):
            v = lip_v * (2*s - s*s)
            if f == 0:
                u = 0.5
            elif i == 2000:
                u = lip_u
            else:
                u = brentq(lambda u: residual(u, v), 0, 1, xtol=1e-14)
            p, du, dv, duu, dvv, duv = gp_Pnt(), *(gp_Vec() for _ in range(5))
            cap.D2(u, v, p, du, dv, duu, dvv, duv)
            du, dv, duu, dvv, duv = [arr(vector) for vector in [du, dv, duu, dvv, duv]]
            if f == 0:
                tangent, acceleration = dv, dvv
            else:
                frame = np.column_stack([du, dv])
                gradient = normal @ frame
                direction = unit(np.array([gradient[1], -gradient[0]]))
                tangent = frame @ direction
                raw_second = duu*direction[0]**2 + dvv*direction[1]**2 + 2*duv*direction[0]*direction[1]
                param_acceleration = -(normal @ raw_second) * gradient / np.dot(gradient, gradient)
                acceleration = raw_second + frame @ param_acceleration
            signed.append(float(side * np.dot(np.cross(tangent, acceleration), normal) / np.linalg.norm(tangent)**3))
            xyz.append(arr(p).tolist())
        signed = np.array(signed)
        signs = np.sign(signed[np.abs(signed) > 2e-5])
        i = int(np.argmax(np.abs(signed)))
        curvature_peaks = []
        for peak in find_peaks(signed, prominence=0.002)[0]:
            curvature_peaks.append(dict(radius_mm=float(1/abs(signed[peak])), point_mm=xyz[peak]))
        edge = next(edge for edge in cap_edges if side * max((edge.startPoint().y, edge.endPoint().y), key=abs) > 0)
        cut_lip = max(edge.startPoint(), edge.endPoint(), key=lambda p: p.z)
        actual_lip_gap = np.linalg.norm(np.array(xyz[-1])-np.array(cut_lip.toTuple()))
        halves.append(dict(side=side, join_curvature_per_mm=float(signed[0]),
                           lip_point_mm=xyz[-1], lip_radius_mm=float(1/abs(signed[-1])),
                           min_radius_mm=float(1/abs(signed[i])), min_radius_point_mm=xyz[i],
                           min_radius_sample=i, signed_curvature_reversals=int(np.count_nonzero(np.diff(signs))),
                           minimum_signed_curvature_per_mm=float(np.min(signed)),
                           interior_curvature_maxima=curvature_peaks, occ_section_lip_gap_mm=float(actual_lip_gap)))
    return dict(fraction=f, cutting_plane=f"x = {offset:.9f} + {slope:.12f} z", total_edges=len(edges),
                cap_samples_per_side=2001, reversal_threshold_per_mm=2e-5, halves=halves,
                method="OCC solid-plane section with curvature evaluated by implicit differentiation on the original STEP faces. The center section uses the stable span derivative at the collapsed UV crown.")


def verify(reach):
    folder = HERE / str(reach)
    path = folder / "wing.step"
    shape = cq.importers.importStep(str(path)).val()
    solids = shape.Solids()
    assert len(solids) == 1
    solid = solids[0]
    assert solid.isValid() and all(shell.Closed() for shell in solid.Shells())
    surfaces = {}
    for side in (-1, 1):
        faces = sorted((face for face in solid.Faces() if side*face.Center().y > 0.1), key=lambda face: face.Center().z)
        assert len(faces) == 2
        surfaces[side] = [BRep_Tool.Surface_s(face.wrapped) for face in faces]
    maximum_ellipse_error = 0.0
    for parent, _ in surfaces.values():
        for u, v in product(np.linspace(0, 1, 101), np.linspace(0, 1, 21)):
            x, y, z = arr(parent.Value(u, v))
            a = 80 + A_SLOPE*z
            maximum_ellipse_error = max(maximum_ellipse_error, abs(((x-80-M_SLOPE*z)/a)**2+(y/(.12*a))**2-1))
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(solid.wrapped, props, 1e-9, True)
    result = dict(step_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), reach_mm=reach,
                  solid_count=len(solids), valid=True, shells_closed=True, faces=len(solid.Faces()),
                  volume_mm3=props.Mass(), parent_ellipse_equation_max_residual=maximum_ellipse_error,
                  joins=check_joins(surfaces), perimeter=perimeter(surfaces, reach),
                  sections=[section(solid, surfaces, f) for f in FRACTIONS])
    (folder / "readback.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(reach, result["step_sha256"], result["joins"]["crossing_count"], flush=True)


if __name__ == "__main__":
    for case in (40, 5):
        verify(case)
