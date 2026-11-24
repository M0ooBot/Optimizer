
"""Quantify continuity and planar wrap cuts on the reimported STEP."""
import hashlib
import json
from pathlib import Path

import cadquery as cq
import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepCheck import BRepCheck_NoError, BRepCheck_Shell
from OCP.gp import gp_Pnt, gp_Vec
from OCP.TopAbs import TopAbs_REVERSED
from scipy.optimize import brentq, newton

from parent_definition import CUT, reach, section, span_jet

HERE = Path(__file__).resolve().parent
STEP = HERE / "wing.step"


def unit(v):
    return v / np.linalg.norm(v)


def angle(a, b):
    return float(np.degrees(np.arctan2(np.linalg.norm(np.cross(a, b)), np.dot(a,b))))


def xyz(p):
    return np.array(p.Coord())


def derivatives(surface, u, v):
    p = gp_Pnt()
    vectors = [gp_Vec() for _ in range(5)]
    surface.D2(float(u), float(v), p, *vectors)
    return xyz(p), *(xyz(vec) for vec in vectors)


def frame(surface, u, v, orientation):
    p, pu, pv, puu, pvv, puv = derivatives(surface, u, v)
    tangent = np.column_stack([pu, pv])
    inverse = np.linalg.pinv(tangent)
    normal = orientation * unit(np.cross(pu, pv))
    second = np.array([[normal @ puu, normal @ puv],
                       [normal @ puv, normal @ pvv]])
    return p, normal, -inverse.T @ second @ inverse


def stripe(frame_data, camera, direction):
    point, normal, normal_gradient = frame_data
    view = unit(camera - point)
    view_gradient = -(np.eye(3) - np.outer(view, view)) / np.linalg.norm(camera - point)
    reflection = 2 * (view @ normal) * normal - view
    derivative = 2 * np.outer(normal, view @ normal_gradient + normal @ view_gradient)
    derivative += 2 * (view @ normal) * normal_gradient - view_gradient
    direction = 4 * direction
    return direction @ reflection, (direction @ derivative) @ (np.eye(3) - np.outer(normal, normal))


def maximum(record, name, value, **location):
    if name not in record or value > record[name]["value"]:
        record[name] = dict(value=float(value), **location)


def classify(solid):
    records = {}
    for face in solid.Faces():
        surface = BRep_Tool.Surface_s(face.wrapped)
        _, _, lo, hi = surface.Bounds()
        middle = surface.Value(.5, (lo + hi)/2)
        side = 1 if middle.Y() > 0 else -1
        name = "cap" if lo < 0 else ("upper" if middle.Z() > 22 else "lower")
        orientation = -1 if face.wrapped.Orientation() == TopAbs_REVERSED else 1
        records[side, name] = surface, orientation
    return records


def verify_joins(faces):
    u_samples = np.linspace(.001, .999, 999)
    cameras = np.array([[125., 0., 1100.], [660., 0., 570.], [-550., 0., 390.],
                        [200., 0., -1000.]])
    directions = [np.array([np.cos(a), np.sin(a), 0.]) for a in np.deg2rad([0, 30, 60, 90, 120, 150])]
    report = {"sample_u_range": [.001, .999], "views_xyz_mm": cameras.tolist(),
              "stripe_angles_degrees": [0, 30, 60, 90, 120, 150], "crossings": 0}
    for side in [-1, 1]:
        cap, cap_sign = faces[side, "cap"]
        for name, q in [("upper", 1), ("lower", -1)]:
            skin, skin_sign = faces[side, name]

            def paired(u):
                return frame(skin, u, CUT, skin_sign), frame(cap, u, q, cap_sign)

            sampled = [paired(u) for u in u_samples]
            for u, (a, b) in zip(u_samples, sampled):
                location = dict(side=side, skin=name, u=u)
                maximum(report, "separation_mm", np.linalg.norm(a[0] - b[0]), **location)
                maximum(report, "tangent_angle_deg", angle(a[1], b[1]), **location)
                maximum(report, "shape_operator_difference_per_mm", np.linalg.norm(a[2] - b[2], ord=2), **location)
            for view_id, camera in enumerate(cameras):
                for stripe_id, direction in enumerate(directions):
                    direction = direction * [1, side, 1]
                    phases = []
                    for u, (a, b) in zip(u_samples, sampled):
                        pa, ga = stripe(a, camera, direction)
                        pb, gb = stripe(b, camera, direction)
                        phases.append(pa)
                        location = dict(side=side, skin=name, u=u, view=view_id, stripe=stripe_id)
                        maximum(report, "stripe_phase_cycles", abs(pa - pb), **location)
                        maximum(report, "stripe_gradient_cycles_per_mm", np.linalg.norm(ga - gb), **location)
                    for i, (p0, p1) in enumerate(zip(phases[:-1], phases[1:])):
                        for level in np.arange(np.ceil(2*min(p0,p1))/2, max(p0,p1), .5):
                            u = brentq(lambda t: stripe(paired(t)[0], camera, direction)[0]-level,
                                       u_samples[i], u_samples[i+1])
                            a, b = paired(u)
                            _, ga = stripe(a, camera, direction)
                            _, gb = stripe(b, camera, direction)
                            if min(np.linalg.norm(ga), np.linalg.norm(gb)) > 1e-9:
                                report["crossings"] += 1
                                maximum(report, "crossing_stripe_angle_deg", angle(ga, gb),
                                        side=side, skin=name, u=u, view=view_id, stripe=stripe_id)
    return report


def planar_sections(faces):
    cap, _ = faces[1, "cap"]
    records = []
    for station in [.02, .035, .055, .085, .12, .18, .24, .32, .41, .5,
                    .59, .68, .76, .82, .88, .915, .945, .965, .98]:
        span_direction, _ = span_jet(CUT, np.array([station]), 0)
        upper = section(CUT, station, 1)
        lower = section(CUT, station, -1)
        origin = section(CUT, station, 0)
        normal = unit(np.cross(span_direction[0], upper-lower))
        q_samples = np.unique(np.r_[np.linspace(-1, 1, 1201), np.linspace(-.04, .04, 401), np.linspace(-1, -.985, 121)])
        positions, curvature = [], []
        u = station
        for q in q_samples:
            u = newton(lambda t: normal @ (xyz(cap.Value(t, q))-origin), u, tol=1e-12)
            if not 0 < u < 1:
                raise RuntimeError(f"Plane section branch leaves cap at {station=}, {q=}, {u=}")
            p, pu, pv, puu, pvv, puv = derivatives(cap, u, q)
            du = -(normal @ pv) / (normal @ pu)
            acceleration = puu*du*du + 2*puv*du + pvv
            tangent = pv + pu*du
            acceleration -= pu * (normal @ acceleration) / (normal @ pu)
            positions.append(p)
            curvature.append((normal @ np.cross(tangent, acceleration)) / np.linalg.norm(tangent)**3)
        p = np.asarray(positions)
        k = np.asarray(curvature)
        i = np.argmax(abs(k))
        lip = np.argmin(abs(q_samples))
        mask = np.flatnonzero(abs(k) > 1e-7)
        crossings = []
        for a, b in zip(mask[:-1], mask[1:]):
            if k[a]*k[b] < 0:
                crossings.append([float(q_samples[a]), float(q_samples[b])])
        records.append({"station_u_at_upper_seam": station, "plane_normal": normal.tolist(),
                        "samples": len(k), "minimum_radius_mm": float(1/abs(k[i])),
                        "minimum_radius_xyz_mm": p[i].tolist(), "minimum_radius_q": float(q_samples[i]),
                        "lip_radius_mm": float(1/abs(k[lip])), "lip_xyz_mm": p[lip].tolist(),
                        "minimum_to_lip_radius_ratio": float(abs(k[lip]/k[i])),
                        "signed_curvature_at_ends_per_mm": [float(k[0]), float(k[-1])],
                        "sign_change_q_brackets": crossings})
    return records


def geometry(faces):
    result = {}
    for side in [-1, 1]:
        for name, sign in [("upper", 1), ("lower", -1)]:
            skin, _ = faces[side, name]
            for u in np.linspace(0, 1, 501):
                for y in np.linspace(0, CUT, 19):
                    expected = section(y, u, sign) * [1, side, 1]
                    maximum(result, "parent_sampling_error_mm", np.linalg.norm(xyz(skin.Value(u,y))-expected),
                            side=side, skin=name, u=u, y_mm=y)
        cap, _ = faces[side, "cap"]
        for u in np.linspace(0, 1, 1201):
            expected = section(CUT+reach(u), u, 0) * [1, side, 1]
            maximum(result, "rim_sampling_error_mm", np.linalg.norm(xyz(cap.Value(u,0))-expected), side=side, u=u)
        result[f"tip_y_side_{side}_mm"] = cap.Value(.5,0).Y()
    result["projected_span_mm"] = result["tip_y_side_1_mm"] - result["tip_y_side_-1_mm"]
    result["root_chord_mm"] = faces[1, "upper"][0].Value(0,0).Distance(faces[1, "upper"][0].Value(1,0))
    result["cut_chord_mm"] = faces[1, "upper"][0].Value(0,CUT).Distance(faces[1, "upper"][0].Value(1,CUT))
    cut_le = xyz(faces[1, "upper"][0].Value(0,CUT))
    cut_te = xyz(faces[1, "upper"][0].Value(1,CUT))
    result["cut_le_xyz_mm"] = cut_le.tolist()
    result["washout_degrees"] = float(np.degrees(np.arctan2(cut_te[2]-cut_le[2], cut_te[0]-cut_le[0])))
    result["dihedral_degrees"] = float(np.degrees(np.arctan2(cut_le[2], cut_le[1])))
    return result


def main():
    imported = cq.importers.importStep(str(STEP)).val()
    shells = imported.Shells()
    closed = len(shells) == 1 and BRepCheck_Shell(shells[0].wrapped).Closed() == BRepCheck_NoError
    assert imported.isValid() and closed and len(imported.Solids()) == 1
    faces = classify(imported)
    data = {
        "step_sha256": hashlib.sha256(STEP.read_bytes()).hexdigest(),
        "valid": imported.isValid(), "closed_shell": closed,
        "solid_count": len(imported.Solids()), "face_count": len(imported.Faces()),
        "volume_mm3": imported.Volume(),
        "geometry": geometry(faces),
        "join_and_zebra": verify_joins(faces),
        "right_cap_planar_sections": planar_sections(faces),
    }
    output = HERE / "readback.json"
    output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
