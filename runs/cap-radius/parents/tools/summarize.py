"""Print the README summary table from spec.json, build.json and check.json."""

import json

from common import PARENTS, SHAPES, VARIANTS


def fmt(x, digits=3):
    if x is None:
        return "n/a"
    if x == 0:
        return "0"
    if abs(x) < 1e-3 or abs(x) >= 1e5:
        return f"{x:.{digits - 1}e}"
    return f"{x:.{digits}g}" if abs(x) < 1 else f"{x:.{digits + 3}g}"


def rows():
    out = []
    for shape in SHAPES:
        for variant in VARIANTS:
            folder = PARENTS / shape / variant
            spec = json.loads((folder / "spec.json").read_text())
            check = json.loads((folder / "check.json").read_text())
            build = json.loads((folder / "build.json").read_text())
            a = spec["attachment_section"]
            r = check["reimport"]
            crossings = spec["attachment_thickness"]["chord_fractions_where_thickness_is_0p8_mm"]
            if variant == "sharp":
                fidelity = check["max_deviation_parent_py_to_accepted_retained_skin_mm"]
                fidelity_step = check["max_deviation_parent_step_to_accepted_retained_skin_mm"]
                te = "sharp"
                conv = ""
            else:
                fidelity = None
                fidelity_step = None
                te_values = [row["step_te_thickness_mm"] for row in check["te_thickness"]]
                te = f"{min(te_values):.6f}-{max(te_values):.6f}"
                agreement = check["converter_agreement"]
                if shape == "wingb":
                    conv = fmt(max(s["max_distance_to_step_mm"] for s in agreement["sections"]))
                else:
                    conv = fmt(max(agreement[k]["max_distance_converter_points_to_analytic_section_mm"]
                                   for k in ("root", "attachment")))
            out.append({
                "shape": shape, "variant": variant,
                "attachment": spec["attachment_stations"],
                "chord": a["chord_mm"], "te_att": a["trailing_edge_thickness_mm"],
                "cross": crossings, "reach": spec["cap_reach_mm"],
                "faces": r["faces"], "valid": r["brepcheck_strict_valid"] and r["shell_closed"] and r["solids"] == 1,
                "volume": r["volume_mm3"], "bytes": check["step_bytes"],
                "fit": check.get("parent_py_to_parent_step_max_mm"),
                "fidelity": fidelity, "fidelity_step": fidelity_step, "te": te, "conv": conv,
                "cut": spec.get("keep_chord_cut"),
            })
    return out


def table():
    lines = ["| Parent | Attachment | Chord at attachment (mm) | TE at attachment (mm) | 0.8 mm thickness at x/c | Cap reach (mm) | STEP faces, strict valid, closed | parent.py vs parent.step (mm) | Sharp: deviation from accepted skin, parent.py / parent.step (mm) | Blunt: TE on STEP (mm) | Blunt: Converter agreement (mm) |",
             "|---|---|---:|---:|---|---:|---|---:|---|---|---:|"]
    for r in rows():
        att = r["attachment"]
        where = ("Y = ±%g" % att["Y_mm"][1]) if "Y_mm" in att else ("Z = %g" % att["Z_mm"][0]) if "Z_mm" in att else "X = -540"
        fid = f"{fmt(r['fidelity'])} / {fmt(r['fidelity_step'])}" if r["variant"] == "sharp" else ""
        lines.append(
            f"| {r['shape']} {r['variant']} | {where} | {r['chord']:.6g} | {r['te_att']:.6f} | "
            f"{', '.join(f'{c:.6f}' for c in r['cross'])} | {r['reach']:g} | {r['faces']}, {'yes' if r['valid'] else 'NO'} | "
            f"{fmt(r['fit']) if r['fit'] is not None else 'same surface'} | {fid} | {r['te'] if r['variant'] == 'blunt' else ''} | {r['conv']} |")
    return "\n".join(lines)


if __name__ == "__main__":
    print(table())
    for r in rows():
        print(r["shape"], r["variant"], "volume", round(r["volume"], 3), "bytes", r["bytes"], "cut", r["cut"])
