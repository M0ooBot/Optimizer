"""Calibrate audit.py on the accepted whole models and compare with the earlier paired section audit.

    python runs/cap-radius/audit/calibration/calibrate.py [--skip-audits]

1. Runs audit.py (native cut family) on the six accepted whole models, whose sharp parents are the sharp
   parents of this trial, writing calibration/accepted/SHAPE-sharp.{json,txt}.
2. Re-traces the exact planes and chord fractions of runs/cap-trials/retest/audit.py with audit.py's section
   code (--family legacy) and compares smallest radius, lip radius and reversals with
   runs/cap-trials/retest/comparison.json ('rebuilt' = the accepted models). For the twisted wing it also
   traces the native family at the same fractions to show what the plane definition changes.
3. Compares the Wing B joins with runs/wing-b-rounded-tip-g2/g2-readback.json.
Writes calibration/legacy_comparison.json and calibration/legacy_comparison.md.
"""

import argparse
import hashlib
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import audit as A  # noqa: E402

ACCEPTED = {
    "flat": "runs/cap-trials/retest/flat/wing.step",
    "twisted": "runs/cap-trials/retest/twisted/wing.step",
    "swept55": "runs/cap-trials/retest/swept55/wing.step",
    "fin40": "runs/cap-trials/retest/fin/40/wing.step",
    "fin5": "runs/cap-trials/retest/fin/5/wing.step",
    "wingb": "runs/wing-b-rounded-tip-g2/wing-b-rounded.step",
}
LEGACY_KEY = {"flat": "flat", "twisted": "twisted", "swept55": "swept55", "fin40": "fin/40", "fin5": "fin/5"}


def run_audits():
    out = HERE / "accepted"
    out.mkdir(exist_ok=True)

    def one(shape):
        cmd = [sys.executable, str(HERE.parent / "audit.py"), shape, "sharp", str(A.REPO / ACCEPTED[shape]),
               "--out", str(out / f"{shape}-sharp.json")]
        r = subprocess.run(cmd, capture_output=True, text=True)
        return shape, r.returncode, r.stdout[-400:] + r.stderr[-2000:]

    with ThreadPoolExecutor(6) as pool:
        for shape, code, tail in pool.map(one, ACCEPTED):
            print(f"audit {shape}: exit {code}")
            if code not in (0, 1):
                print(tail)


def legacy_comparison():
    reference = json.loads((A.REPO / "runs/cap-trials/retest/comparison.json").read_text())
    rows = {}
    for shape, key in LEGACY_KEY.items():
        ref = reference[key]
        path = A.REPO / ACCEPTED[shape]
        if hashlib.sha256(path.read_bytes()).hexdigest() != ref["step_sha256"][1]:
            raise ValueError(f"{path} is not the model comparison.json measured")
        legacy = A.Case(shape, "sharp", "legacy")
        native = A.Case(shape, "sharp", "native")
        model_l = A.Model(path, legacy)
        model_n = A.Model(path, native)
        tip_l, tip_n = legacy.tips[0], native.tips[0]
        cut_rows = []
        for c in ref["sections"]:
            f = c["chord_fraction"]
            r = A.trace(model_l, legacy, tip_l, f, "legacy plane")
            old = c["rebuilt"]
            row = {"fraction": f, "previous_min_mm": old["minimum_mm"], "previous_lip_mm": old["lip_mm"],
                   "previous_branches_with_reversal": old["branches_with_reversal"], "previous_min_xyz": old["min_xyz_mm"]}
            if "error" in r:
                row["error"] = r["error"]
            else:
                row.update(min_mm=r["min_radius"]["radius_mm"], lip_mm=r["lip"]["radius_mm"], reversals=len(r["reversals"]),
                           min_xyz=r["min_radius"]["xyz"],
                           min_relative_difference=r["min_radius"]["radius_mm"] / old["minimum_mm"] - 1,
                           min_location_difference_mm=float(np.linalg.norm(np.array(r["min_radius"]["xyz"]) - old["min_xyz_mm"])),
                           lip_relative_difference=(r["lip"]["radius_mm"] / old["lip_mm"] - 1) if r["lip"]["radius_mm"] else None)
            if shape == "twisted":
                rn = A.trace(model_n, native, tip_n, f, "native plane")
                if "error" not in rn:
                    row.update(native_min_mm=rn["min_radius"]["radius_mm"], native_lip_mm=rn["lip"]["radius_mm"],
                               native_reversals=len(rn["reversals"]))
            cut_rows.append(row)
        O, n, _, _ = tip_l.plane(0.5)
        On, nn, _, _ = tip_n.plane(0.5)
        rows[shape] = {"cuts": cut_rows, "plane_normal_at_0.5": {"legacy": n, "native": nn},
                       "native_equals_legacy_planes": bool(np.allclose(n, nn, atol=1e-12) and abs((On - O) @ n) < 1e-9)}
        print(shape, "done", flush=True)
    return rows


def wingb_joins():
    g2 = json.loads((A.REPO / "runs/wing-b-rounded-tip-g2/g2-readback.json").read_text())
    mine = json.loads((HERE / "accepted" / "wingb-sharp.json").read_text())
    prev = {"gap_mm": 0.0, "angle_deg": 0.0, "ii_difference_norm_per_mm": 0.0}
    for j in g2["joins"]:
        for p in j["points"]:
            prev["gap_mm"] = max(prev["gap_mm"], p["gap_mm"])
            prev["angle_deg"] = max(prev["angle_deg"], p["tangent_angle_deg"])
            prev["ii_difference_norm_per_mm"] = max(prev["ii_difference_norm_per_mm"], float(np.linalg.norm(p["ii_difference_per_mm"])))
    s = mine["joins"]["summary"]["parent skin to cap"]
    return {"previous (check_g2.py, 85 points per side, chord fractions 1e-8 to 1)": prev,
            "audit.py (399 points per edge)": {"gap_mm": s["gap_mm"]["max"], "angle_deg": s["angle_deg"]["max"],
                                              "curvature_mismatch_per_mm": s["curvature_mismatch_per_mm"]["max"],
                                              "curvature_mismatch_including_near_collapse_per_mm": s["curvature_mismatch_per_mm"]["max_including_near_collapse"]}}


def markdown(rows, wingb):
    lines = ["# Section radii on the earlier audit's planes", "",
             "Smallest radius and lip radius on the cuts of runs/cap-trials/retest/audit.py, retraced with audit.py "
             "(--family legacy) on the same accepted STEP files. 'Previous' is the 'rebuilt' column of "
             "runs/cap-trials/retest/comparison.json. Regenerate with calibration/calibrate.py.", ""]
    for shape, r in rows.items():
        head = "| x/c | previous min (mm) | audit.py min (mm) | rel. diff | previous lip (mm) | audit.py lip (mm) | reversals prev / now |"
        rule = "|---:|---:|---:|---:|---:|---:|---:|"
        if shape == "twisted":
            head += " native planes: min / lip (mm) |"
            rule += "---:|"
        lines += [f"## {shape}", "", f"Native planes identical to the earlier planes: {r['native_equals_legacy_planes']}.", "", head, rule]
        for c in r["cuts"]:
            if "error" in c:
                lines.append(f"| {c['fraction']} | {c['previous_min_mm']:.6g} | error: {c['error']} | | | | |")
                continue
            extra = ""
            if shape == "twisted":
                extra = f" {c.get('native_min_mm', float('nan')):.6g} / {c.get('native_lip_mm') or float('nan'):.6g} |"
            lines.append(f"| {c['fraction']} | {c['previous_min_mm']:.6g} | {c['min_mm']:.6g} | {c['min_relative_difference']:+.1e} | "
                         f"{c['previous_lip_mm']:.6g} | {c['lip_mm']:.6g} | {c['previous_branches_with_reversal']} / {c['reversals']} |" + extra)
        lines.append("")
    lines += ["## Wing B joins", "", "~~~", json.dumps(A.clean(wingb), indent=1), "~~~", ""]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--skip-audits", action="store_true")
    a = ap.parse_args()
    if not a.skip_audits:
        run_audits()
    rows = legacy_comparison()
    wb = wingb_joins()
    (HERE / "legacy_comparison.json").write_text(json.dumps(A.clean({"sections": rows, "wingb_joins": wb}), indent=1) + "\n")
    (HERE / "legacy_comparison.md").write_text(markdown(A.clean(rows), wb))
    print((HERE / "legacy_comparison.md").read_text())


if __name__ == "__main__":
    main()
