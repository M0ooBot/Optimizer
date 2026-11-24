"""Build and verify a fin with a user-specified cap reach in an isolated folder."""

import argparse
import json
from math import isfinite

import build_fin


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cap-reach-mm", type=float, required=True)
    args = parser.parse_args()
    reach = args.cap_reach_mm
    if not isfinite(reach) or reach <= 0:
        parser.error("--cap-reach-mm must be a finite positive number")
    if build_fin.radius(build_fin.CUT + reach) <= 0:
        parser.error("--cap-reach-mm extends beyond the positive radius of the parent taper")
    destination = build_fin.HERE / f"cap-{reach}mm"
    destination.mkdir(exist_ok=True)
    build_fin.HEIGHT = build_fin.CUT + reach
    build_fin.HERE = destination
    build_fin.OUTPUT = destination / "vertical-fin.step"

    # Import the analyzers after setting this run's geometry and paths.
    import check_fin
    import render_fin

    solid = build_fin.make_fin()
    solid.exportStep(str(build_fin.OUTPUT), write_pcurves=True)
    evidence = check_fin.check(build_fin.OUTPUT)
    (destination / "readback.json").write_text(json.dumps(evidence, indent=2) + "\n")
    render_fin.render_overview(check_fin.patches(build_fin.OUTPUT)[0])
    render_fin.render_sections()
    print(f"Verified STEP and evidence in {destination}")


if __name__ == "__main__":
    main()
