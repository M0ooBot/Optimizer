"""Measure the round-trip STEP without importing its geometry generator."""
import json
from pathlib import Path

import cadquery as cq
import numpy as np
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.BRepBndLib import BRepBndLib
from OCP.Bnd import Bnd_Box
from OCP.IntCurvesFace import IntCurvesFace_ShapeIntersector
from OCP.gp import gp_Dir, gp_Lin, gp_Pln, gp_Pnt
from scipy.optimize import minimize_scalar


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "paired-wing-trial/a/measured/01-wing-occt.step"


def section(solid, point, normal):
    plane = gp_Pln(gp_Pnt(*point), gp_Dir(*normal))
    cut = BRepAlgoAPI_Section(solid.wrapped, plane, True)
    if not cut.IsDone():
        raise ValueError(f"Section failed at {point}")
    result = cq.Shape.cast(cut.Shape())
    if not result.Edges():
        raise ValueError(f"Empty section at {point}")
    return result


def bounds(shape):
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape.wrapped, box, False, False)
    return box.Get()


def sampled_delta(first, second):
    return max(second.distance(cq.Vertex.makeVertex(*p.toTuple()))
               for edge in first.Edges() for p in edge.sample(61)[0])


def lip_measurement(solid, probe, aft):
    # Chordwise position is held fixed, producing a true span-height trace.
    cut = section(solid, (0, 0, -aft), (0, 0, 1))
    candidates = []
    for edge in cut.Edges():
        curve = BRepAdaptor_Curve(edge.wrapped)
        grid = np.linspace(curve.FirstParameter(), curve.LastParameter(), 129)
        index = min(range(len(grid)), key=lambda i: curve.Value(grid[i]).X())
        lo, hi = grid[max(0, index - 1)], grid[min(len(grid) - 1, index + 1)]
        opt = minimize_scalar(lambda t: curve.Value(t).X(), bounds=(lo, hi), method="bounded")
        t = min([lo, hi, opt.x, grid[index]], key=lambda t: curve.Value(t).X())
        point = curve.Value(t)
        candidates.append((point.X(), point.Y()))
    x, y = min(candidates)
    # Section curves can be polylines. Intersect the actual STEP faces on both
    # sides of the rim to recover x(y) and estimate its local bending radius.
    heights = np.array([-.12, -.08, -.04, .04, .08, .12])
    xs = []
    for dy in heights:
        ray = gp_Lin(gp_Pnt(-605, y + dy, -aft), gp_Dir(1, 0, 0))
        probe.Perform(ray, 0, 75)
        if not probe.IsDone() or not probe.NbPnt():
            raise ValueError(f"No skin at {aft} mm aft, {y + dy} mm high")
        xs.append(min(probe.Pnt(i).X() for i in range(1, probe.NbPnt() + 1)) - x)
    a, slope, intercept = np.polyfit(heights, xs, 2)
    radius = (1 + slope**2)**1.5 / (2 * a) if a > 0 else None
    residual = np.max(np.abs(np.polyval([a, slope, intercept], heights) - xs))
    return dict(aft_mm=aft, outboard_mm=-x, height_mm=y, fitted_section_radius_mm=radius,
                fit_residual_mm=float(residual),
                tangent_departure_from_vertical_deg=float(np.degrees(np.arctan(abs(slope)))))


def main():
    original = cq.importers.importStep(str(PARENT)).val()
    model = cq.importers.importStep(str(HERE / "wing-b-rounded.step")).val()
    result = {"solid_count": len(model.Solids()), "valid": model.isValid(),
              "closed_shells": [shell.Closed() for shell in model.Shells()],
              "bounds_xyz_mm": bounds(model),
              "inboard_sections": [], "outer_sections": []}
    for span in (50, 300, 510, 539):
        plane = ((-span, 0, 0), (1, 0, 0))
        before, after = section(original, *plane), section(model, *plane)
        error = max(sampled_delta(before, after), sampled_delta(after, before))
        result["inboard_sections"].append(dict(span_mm=span, sampled_error_mm=error))
    for span in (540.1, 550, 565, 580, 590, 595, 599, 599.9):
        box = bounds(section(model, (-span, 0, 0), (1, 0, 0)))
        result["outer_sections"].append(dict(span_mm=span, chord_mm=box[5] - box[2],
                                             section_height_mm=box[4] - box[1]))
    probe = IntCurvesFace_ShapeIntersector()
    probe.Load(model.wrapped, 1e-8)
    result["rim_cuts"] = [lip_measurement(model, probe, aft) for aft in (90, 110, 135)]
    print(json.dumps(result, indent=2, allow_nan=False))
    assert result["valid"] and result["solid_count"] == 1 and all(result["closed_shells"])
    assert abs(result["bounds_xyz_mm"][0] + 600) < .25
    assert max(row["sampled_error_mm"] for row in result["inboard_sections"]) < .25
    chords = [row["chord_mm"] for row in result["outer_sections"]]
    assert all(a > b > 0 for a, b in zip(chords, chords[1:]))
    for row in result["rim_cuts"]:
        assert row["fitted_section_radius_mm"] is not None and row["fitted_section_radius_mm"] > 0
        assert row["fit_residual_mm"] < .005
        assert row["tangent_departure_from_vertical_deg"] < 3
    (HERE / "readback.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
