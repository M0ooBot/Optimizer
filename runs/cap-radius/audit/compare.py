"""Run audit.py over several whole models and write one summary table (Markdown and JSON).

    python runs/cap-radius/audit/compare.py --entry SHAPE VARIANT MODEL.step [--entry ...]
    python runs/cap-radius/audit/compare.py --list entries.json      # [{"shape", "variant", "step", "label"?}, ...]
    python runs/cap-radius/audit/compare.py --builds                 # every runs/cap-radius/builds/S/V/wing.step

Options: --out-dir DIR (default runs/cap-radius/audit/results), --jobs N (parallel audits, default 3),
--reuse (keep an existing audit JSON whose STEP hash matches), --name NAME (table file stem, default summary).
Each audit is written to DIR/LABEL.json and DIR/LABEL.txt; the table goes to DIR/NAME.md and DIR/NAME.json.
"""

import argparse
import hashlib
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
BUILDS = REPO / "runs" / "cap-radius" / "builds"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_one(entry, out_dir, reuse):
    label = entry.get("label") or f"{entry['shape']}-{entry['variant']}"
    out = out_dir / f"{label}.json"
    if reuse and out.exists():
        try:
            if json.loads(out.read_text())["step_sha256"] == sha(entry["step"]):
                return label, out, "reused"
        except (KeyError, ValueError):
            pass
    cmd = [sys.executable, str(HERE / "audit.py"), entry["shape"], entry["variant"], str(entry["step"]), "--out", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode not in (0, 1) or not out.exists():
        return label, None, (r.stdout[-1000:] + r.stderr[-3000:])
    return label, out, "ran"


def fmt(x, spec=".3g"):
    return "" if x is None else format(x, spec)


def row(label, entry, report):
    req = {r["name"]: r for r in report["summary"]["requirements"]}
    v, p, j = report["validity"], report["preservation"], report["joins"]["summary"].get("parent skin to cap") or {}
    cuts = [c for rows in report["sections"].values() for c in rows if "error" not in c]
    thick = [c for c in cuts if c["zone"] == "thick"]
    thin = [c for c in cuts if c["zone"] == "thin"]
    wt = min(thick, key=lambda c: c["min_radius"]["radius_mm"]) if thick else None
    wn = min(thin, key=lambda c: c["ratio_to_half_thickness"] or 0) if thin else None
    scan = report["outline"]["lip_scan"]
    cm = report["cap_map"]
    w = cm.get("worst") or {}
    te = report["te_closure"]
    return {
        "label": label, "shape": report["shape"], "variant": report["variant"], "step": str(entry["step"]),
        "step_sha256": report["step_sha256"], "overall_pass": report["summary"]["pass"],
        "requirements": {k: r["pass"] for k, r in req.items()},
        "valid": req["1 validity"]["pass"], "extent_difference_mm": v["extent"]["difference_mm"],
        "volume_mm3": v["volume_mm3"], "cap_volume_mm3": v["cap_volume_mm3"],
        "retained_skin_max_mm": p["max_mm"],
        "outline_max_deviation_mm": max((s.get("max_deviation_mm") or 0) for s in scan.values()),
        "join_gap_mm": (j.get("gap_mm") or {}).get("max"), "join_angle_deg": (j.get("angle_deg") or {}).get("max"),
        "join_curvature_mismatch_per_mm": (j.get("curvature_mismatch_per_mm") or {}).get("max"),
        "sections_failed_to_trace": sum("error" in c for rows in report["sections"].values() for c in rows),
        "section_thick_min_radius_mm": wt["min_radius"]["radius_mm"] if wt else None,
        "section_thick_min_at": wt["fraction"] if wt else None,
        "section_thick_cuts_below_0p4": sum(not c["pass"] for c in thick),
        "section_thin_min_ratio": wn["ratio_to_half_thickness"] if wn else None,
        "section_thin_min_at": wn["fraction"] if wn else None,
        "section_thin_cuts_below_0p9": sum(not c["pass"] for c in thin),
        "section_cuts": len(cuts),
        "reversals": sum(len(c["reversals"]) for c in cuts),
        "map_worst_radius_mm": w.get("radius_mm"), "map_worst_fraction": w.get("chord_fraction"), "map_worst_zone": w.get("zone"),
        "map_worst_margin_mm": w.get("margin_mm"), "map_worst_ratio": w.get("ratio_to_half_thickness"),
        "map_tightest_radius_mm": (cm.get("tightest") or {}).get("radius_mm"),
        "te_closure_pass": te["pass"], "te_closes": all(t["closes"] for t in te["tips"]),
    }


def markdown(rows):
    head = ["model", "overall", "valid", "extent Δ (mm)", "skin dev (mm)", "outline dev (mm)", "join gap (mm) / angle (deg)",
            "join curv. mismatch (1/mm)", "sections: thick min R (mm) @ x/c [cuts < 0.4]", "sections: thin min R/(t/2) @ x/c [cuts < 0.9]",
            "reversals", "map worst R (mm) @ x/c, zone: margin or ratio", "map tightest R (mm)", "TE closure"]
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    pf = lambda b: "PASS" if b else "FAIL"
    for r in rows:
        if "error" in r:
            out.append(f"| {r['label']} | ERROR | " + " | ".join([""] * (len(head) - 2)) + " |")
            continue
        mw = (f"{fmt(r['map_worst_radius_mm'], '.4g')} @ {fmt(r['map_worst_fraction'], '.5g')}, {r['map_worst_zone']}: "
              + (f"{fmt(r['map_worst_margin_mm'], '+.3g')} mm" if r["map_worst_zone"] == "thick" else f"{fmt(r['map_worst_ratio'], '.3g')}"))
        out.append("| " + " | ".join([
            f"{r['shape']} {r['variant']}" + (f" ({r['label']})" if r["label"] != f"{r['shape']}-{r['variant']}" else ""),
            pf(r["overall_pass"]), pf(r["valid"]), fmt(r["extent_difference_mm"], ".2e"), fmt(r["retained_skin_max_mm"], ".2e"),
            fmt(r["outline_max_deviation_mm"], ".2e"), f"{fmt(r['join_gap_mm'], '.1e')} / {fmt(r['join_angle_deg'], '.1e')}",
            fmt(r["join_curvature_mismatch_per_mm"], ".3g"),
            f"{fmt(r['section_thick_min_radius_mm'], '.4g')} @ {fmt(r['section_thick_min_at'], '.4g')} [{r['section_thick_cuts_below_0p4']}]",
            f"{fmt(r['section_thin_min_ratio'], '.3g')} @ {fmt(r['section_thin_min_at'], '.6g')} [{r['section_thin_cuts_below_0p9']}]",
            str(r["reversals"]), mw, fmt(r["map_tightest_radius_mm"], ".4g"), pf(r["te_closure_pass"])]) + " |")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--entry", nargs=3, action="append", metavar=("SHAPE", "VARIANT", "STEP"), default=[])
    ap.add_argument("--list", type=Path)
    ap.add_argument("--builds", action="store_true", help="audit every runs/cap-radius/builds/S/V/wing.step")
    ap.add_argument("--out-dir", type=Path, default=HERE / "results")
    ap.add_argument("--name", default="summary")
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--reuse", action="store_true")
    a = ap.parse_args()
    entries = [{"shape": s, "variant": v, "step": Path(p)} for s, v, p in a.entry]
    if a.list:
        entries += [dict(e, step=Path(e["step"])) for e in json.loads(a.list.read_text())]
    if a.builds:
        for step in sorted(BUILDS.glob("*/*/wing.step")):
            entries.append({"shape": step.parent.parent.name, "variant": step.parent.name, "step": step})
    if not entries:
        ap.error("no models given")
    a.out_dir.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(a.jobs) as pool:
        results = list(pool.map(lambda e: (e, *run_one(e, a.out_dir, a.reuse)), entries))
    rows = []
    for entry, label, out, status in results:
        if out is None:
            rows.append({"label": label, "error": status})
            print(f"{label}: audit failed\n{status}")
            continue
        report = json.loads(out.read_text())
        rows.append(row(label, entry, report))
        print(f"{label}: {status}, {'PASS' if report['summary']['pass'] else 'FAIL'}")
    table = markdown(rows)
    (a.out_dir / f"{a.name}.md").write_text(f"# Cap audit summary\n\nGenerated by runs/cap-radius/audit/compare.py. "
                                            f"Per-model detail: LABEL.txt and LABEL.json in this folder.\n\n{table}")
    (a.out_dir / f"{a.name}.json").write_text(json.dumps(rows, indent=1, default=str) + "\n")
    print(table)


if __name__ == "__main__":
    main()
