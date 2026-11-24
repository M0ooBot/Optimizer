"""Measure signed curvature along oblique sections of the exported 55° wing."""

import hashlib
import json

import numpy as np
from OCP.gp import gp_Pnt, gp_Vec
from scipy.optimize import brentq

from build import STATION, SWEEP
from check import FRACTIONS, HERE, STEP, load, unit, xyz


PLANE_NORMAL = np.array([1.0, -SWEEP, 0.0])
OUTBOARD = unit(np.array([SWEEP, 1.0, 0.0]))
PARENT_SAMPLES = np.linspace(0.94, 1.0, 81)
CAP_SAMPLES = np.linspace(0.0, 1.0, 1001)
SIGN_TOLERANCE = 1e-9


def section_point(patch, offset, v):
    """Recover a plane–surface intersection point and its first two derivatives."""
    u = brentq(lambda u: PLANE_NORMAL @ xyz(patch, u, v) - offset, 0, 1)
    point = gp_Pnt()
    derivatives = [gp_Vec() for _ in range(5)]
    patch.D2(u, float(v), point, *derivatives)
    du, dv, duu, dvv, duv = [np.array(vector.Coord()) for vector in derivatives]
    # Differentiate the constraint PLANE_NORMAL @ S(u(v), v) = offset twice.
    chord_slope = PLANE_NORMAL @ du
    if chord_slope <= 0:
        raise ValueError(f"Expected an increasing chord coordinate at u={u}, v={v}")
    u_first = -(PLANE_NORMAL @ dv) / chord_slope
    first = dv + du * u_first
    second = dvv + 2 * duv * u_first + duu * u_first**2
    second -= du * (PLANE_NORMAL @ second) / chord_slope
    curvature = ((OUTBOARD @ first) * second[2] - (OUTBOARD @ second) * first[2])
    curvature /= np.linalg.norm(first)**3
    return dict(u=float(u), v=float(v), xyz_mm=list(point.Coord()),
                signed_curvature_per_mm=float(curvature))


def extrema(samples):
    def curvature(sample):
        return sample["signed_curvature_per_mm"]

    return dict(minimum=min(samples, key=curvature), maximum=max(samples, key=curvature))


def audit_section(parent, cap, fraction, cut_chord):
    offset = cut_chord * fraction
    parent_trace = [section_point(parent, offset, v) for v in PARENT_SAMPLES]
    cap_trace = [section_point(cap, offset, v) for v in CAP_SAMPLES]
    signed_extrema = extrema(cap_trace)
    tightest = max(signed_extrema.values(), key=lambda sample: abs(sample["signed_curvature_per_mm"]))
    lip = cap_trace[-1]
    signs = [sample for sample in cap_trace if abs(sample["signed_curvature_per_mm"]) > SIGN_TOLERANCE]
    inflections = [
        section_point(cap, offset, brentq(
            lambda v: section_point(cap, offset, v)["signed_curvature_per_mm"],
            before["v"], after["v"]))
        for before, after in zip(signs, signs[1:])
        if before["signed_curvature_per_mm"] * after["signed_curvature_per_mm"] < 0
    ]
    return dict(
        cut_chord_fraction=fraction,
        plane_offset_mm=offset,
        path_start=parent_trace[0],
        attachment=cap_trace[0],
        parent_signed_curvature_extrema=extrema(parent_trace),
        cap_signed_curvature_extrema=signed_extrema,
        smallest_sampled_cap_radius=dict(tightest, radius_mm=1 / abs(tightest["signed_curvature_per_mm"])),
        lip=dict(lip, radius_mm=1 / abs(lip["signed_curvature_per_mm"])),
        resolved_cap_inflections=inflections,
    )


def main():
    _, patches = load()
    parent, cap = patches[1, 1, False], patches[1, 1, True]
    cut_chord = float(np.linalg.norm(xyz(parent, 1, 1) - xyz(parent, 0, 1)))
    report = dict(
        step_file=STEP.name,
        step_sha256=hashlib.sha256(STEP.read_bytes()).hexdigest(),
        geometry="Right upper parent and cap. Reflections have the same radius, with reversed "
                 "signed curvature on the lower face in the corresponding outboard/up frame.",
        section_paths=dict(
            equation="x - tan(55°)*y = cut_chord_mm * cut_chord_fraction",
            plane_normal=PLANE_NORMAL.tolist(),
            cut_chord_mm=cut_chord,
            signed_curvature_frame=dict(outboard=OUTBOARD.tolist(), up=[0, 0, 1]),
        ),
        coverage=dict(
            attachment_y_mm=STATION,
            parent_start_y_mm=STATION * float(PARENT_SAMPLES[0]),
            parent_v_interval=[float(PARENT_SAMPLES[0]), float(PARENT_SAMPLES[-1])],
            cap_v_interval=[float(CAP_SAMPLES[0]), float(CAP_SAMPLES[-1])],
            face_trim="Rectangular [0, 1] × [0, 1] on each STEP surface. Solve u in [0, 1] "
                      "at each v, from the retained parent through the cap to its lip at v=1.",
            parent_uniform_v_samples=len(PARENT_SAMPLES),
            cap_uniform_v_samples=len(CAP_SAMPLES),
            inflection_bracket_curvature_tolerance_per_mm=SIGN_TOLERANCE,
            scope="Extrema use the sampled section points. Inflections use bracketed sign "
                  "changes. The section planes omit the collapsed end vertices.",
        ),
        sections=[audit_section(parent, cap, f, cut_chord) for f in FRACTIONS],
    )
    (HERE / "section-audit.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
