"""Local stdio MCP tools for the SD7037 Phase 0 study."""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
import math
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory
from typing import Any

from mcp.server.fastmcp import FastMCP


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "aircraft" / "phase0_sd7037"
BASELINE = REPO / "aircraft" / "phase0_shape_study" / "baseline.json"
CONVERTER_REVISION = "d77bbdf45a8c0aa92efbc2a3d3007cfe4f904105"
LABEL = re.compile(r"[a-z][a-z0-9_-]{0,63}", re.ASCII)


def checked_label(value: str) -> str:
    if not LABEL.fullmatch(value) or Path(value).is_reserved() or value == "baseline":
        raise ValueError("Use a lowercase revision label starting with a letter, other than 'baseline'.")
    return value


def make_server(study_dir: Path, converter_source: Path | None = None) -> FastMCP:
    from planeopt.phase0_geometry import (
        apply_edit,
        build_curves,
        calibrate_projected_span,
        hinge_fit,
        load_definition,
        projected_skin_span,
        save_definition,
        section_frame,
        write_curves,
    )

    server = FastMCP("SD7037 Phase 0")

    def definition_path(label: str) -> Path:
        if label == "baseline":
            return BASELINE
        path = (study_dir / checked_label(label) / "definition.json").resolve()
        path.relative_to(study_dir)
        return path

    @server.tool()
    def phase0_reference() -> dict[str, Any]:
        """Read the straight study definition and reference dimensions in millimetres."""
        return {
            "definition": json.loads(BASELINE.read_text(encoding="utf-8")),
            "pinned_converter_revision": CONVERTER_REVISION,
            "axes": {"span": "-X", "up": "+Y", "nose": "+Z"},
            "section_twist": {
                "tip_twist_deg": "Tip incidence relative to the fixed zero-degree root. Positive lowers the trailing edge.",
                "midspan_twist_peak_deg": (
                    "Optional incidence at half-span. Two quintic ramps join there with zero slope and curvature. "
                    "Omit on a new definition or set null in an edit for linear root-to-tip twist."
                ),
                "pivot": "Local leading edge",
            },
            "wing_locks": {
                "root": "root section and placement",
                "projected_span": "maximum sampled outboard skin X distance from the root plane",
                "tip_height": "leading-edge Y coordinate when requested",
                "tip_aft": "leading-edge aft displacement (sweep_mm + rake_mm) when requested",
            },
            "fixture": str(FIXTURE),
            "study_dir": str(study_dir),
        }

    @server.tool()
    def phase0_verify_converter_fixture() -> dict[str, Any]:
        """Regenerate the fixture with Airfoil Converter and compare every generated curve."""
        try:
            from airfoil_converter import store, wing_build, writer
        except ModuleNotFoundError as exc:
            if exc.name != "airfoil_converter":
                raise
            raise ValueError("Start the server with --converter-source for fixture verification.") from exc

        sidecar = store.load(str(FIXTURE / "SD7037-Baseline.airfoils.json"))
        for rib in sidecar.exports:
            rib.source = str(FIXTURE / "sd7037-il.csv")
        record, = sidecar.wings
        result = wing_build.build_wing(
            record.spec, sidecar, record.stem, own_names=record.feature_names(), check=True
        )
        expected = {curve.feature for curve in record.curves if not curve.retired}
        regenerated = {curve.feature: writer.curve_bytes(curve.points) for curve in result.curves}
        mismatches = sorted(expected.symmetric_difference(regenerated))
        for feature in expected.intersection(regenerated):
            stored = FIXTURE / "curves" / f"{feature}.sldcrv"
            if not stored.is_file() or stored.read_bytes() != regenerated[feature]:
                mismatches.append(feature)
        checkout_revision = None
        if converter_source is not None:
            try:
                checkout = subprocess.run(
                    ["git", "-C", str(converter_source), "rev-parse", "HEAD"],
                    capture_output=True, text=True, check=False, timeout=5,
                )
                if checkout.returncode == 0:
                    checkout_revision = checkout.stdout.strip()
            except (OSError, subprocess.TimeoutExpired):
                pass
        return {
            "matches": not mismatches,
            "mismatches": sorted(mismatches),
            "curves": len(regenerated),
            "sections": result.station_count,
            "pinned_converter_revision": CONVERTER_REVISION,
            "checkout_revision": checkout_revision,
            "checkout_matches_pin": (
                checkout_revision == CONVERTER_REVISION if checkout_revision else None
            ),
        }

    @server.tool()
    def phase0_shape(
        revision: str,
        parent: str = "baseline",
        wing: dict[str, Any] | None = None,
        tail: dict[str, Any] | None = None,
        lock_tip_height: bool = False,
        lock_tail: bool = False,
        lock_tip_aft: bool = False,
        calibrate_wing_span: bool = True,
        brief: str | None = None,
    ) -> dict[str, Any]:
        """Change named wing or tail parameters and save curves and an editable JSON definition.

        Use phase0_reference for field names and units. For a local wing follow-up,
        pass the preceding revision as parent and set lock_tip_height and lock_tail
        when the brief fixes them. Set lock_tip_aft to hold the leading-edge tip
        aft coordinate. Set wing.span_curve to full_span to distribute sweep
        and rise from the root, with no separate onset parameters.
        Set midspan_twist_peak_deg with tip_twist_deg for smooth twist through
        a prescribed half-span incidence, or set it to null for linear twist.
        Span and root locks apply to every edit.
        Span calibration trims the leading-edge guide to meet the skin envelope.
        """
        revision = checked_label(revision)
        destination = study_dir / revision
        if destination.exists():
            raise ValueError(f"Revision directory {destination} already exists.")
        previous = load_definition(definition_path(parent))
        definition = apply_edit(
            previous, name=revision, wing=wing, tail=tail, lock_tip_height=lock_tip_height,
            lock_tail=lock_tail, lock_tip_aft=lock_tip_aft,
        )
        if calibrate_wing_span:
            definition = replace(definition, wing=calibrate_projected_span(definition.wing, FIXTURE))
        span = projected_skin_span(definition.wing, FIXTURE)
        if abs(span - definition.wing.projected_span_mm) > definition.tolerance_mm:
            raise ValueError("Sampled wing skin span conflicts with the locked span. Enable span calibration.")
        curves = build_curves(definition, FIXTURE)
        measurements = {}
        for name, surface in (("wing", definition.wing), ("tail", definition.tail)):
            margin, station = hinge_fit(surface, FIXTURE)
            frame = section_frame(surface, 1.0)
            measurements[name] = {
                "leading_edge_projected_span_mm": surface.root_le_mm[0] - frame.le[0],
                "target_projected_span_mm": surface.projected_span_mm,
                "sampled_projected_skin_span_mm": span if name == "wing" else projected_skin_span(surface, FIXTURE),
                "tip_le_xyz_mm": frame.le,
                "tip_quarter_chord_xyz_mm": frame.quarter_chord,
                "minimum_sampled_hinge_skin_margin_mm": margin,
                "minimum_hinge_margin_station_fraction": station,
            }
        with TemporaryDirectory(prefix=".shape-", dir=study_dir) as temp:
            staged = Path(temp) / revision
            staged.mkdir()
            save_definition(definition, staged / "definition.json")
            curve_dir = staged / "curves"
            curve_dir.mkdir()
            write_curves(curves, curve_dir)
            (staged / "request.json").write_text(
                json.dumps(
                    {
                        "revision": revision, "parent": parent, "wing": wing, "tail": tail,
                        "lock_tip_height": lock_tip_height, "lock_tail": lock_tail,
                        "lock_tip_aft": lock_tip_aft,
                        "calibrate_wing_span": calibrate_wing_span,
                        "brief": brief,
                    },
                    indent=2, allow_nan=False,
                ) + "\n",
                encoding="utf-8",
            )
            staged.rename(destination)
        return {
            "definition": str(destination / "definition.json"),
            "curves": str(destination / "curves"),
            "curve_count": len(curves),
            "request": str(destination / "request.json"),
            "generated_geometry_measurements": measurements,
        }

    @server.tool()
    def phase0_report(
        report_name: str,
        revisions: list[str],
        step_files: dict[str, dict[str, str]] | None = None,
        fallback: bool = False,
        rebuild_of: dict[str, str] | None = None,
    ) -> dict[str, Any]:
        """Write comparison plots and measurements for ordered revisions.

        Optional step_files maps each revision to wing or tail STEP paths relative
        to the study directory. Set fallback to generate Open CASCADE STEP lofts
        for surfaces without supplied STEP files. Install the cad extra for STEP.
        CAD comparisons use the saved reference STEP for baseline wing unless
        you explicitly supply a baseline wing STEP of your own.
        Put baseline first to measure against its locked dimensions. Pass
        rebuild_of as a mapping from rebuilt revision to original revision.
        """
        from planeopt.phase0_evidence import write_evidence

        if not revisions or revisions[0] != "baseline" or len(set(revisions)) != len(revisions):
            raise ValueError("Supply distinct revision labels with 'baseline' first.")
        if any(new == original or new not in revisions or original not in revisions
               for new, original in (rebuild_of or {}).items()):
            raise ValueError("Each rebuild comparison needs two distinct listed revisions.")
        destination = study_dir / checked_label(report_name)
        if destination.exists():
            raise ValueError(f"Report directory {destination} already exists.")
        definitions = {label: load_definition(definition_path(label)) for label in revisions}
        step_paths = {}
        step_sources = {}
        if fallback or step_files:
            step_paths["baseline", "wing"] = FIXTURE / "SD7037-Baseline.step"
            step_sources["baseline", "wing"] = "reference_step"
        for label, surfaces in (step_files or {}).items():
            if label not in definitions or surfaces.keys() - {"wing", "tail"}:
                raise ValueError(f"Unknown revision or surface in step_files for {label}.")
            for surface, filename in surfaces.items():
                relative_path = Path(filename)
                if relative_path.is_absolute():
                    raise ValueError(f"Give a relative STEP path for {label} {surface}.")
                path = (study_dir / relative_path).resolve()
                path.relative_to(study_dir)
                if path.suffix.lower() not in {".step", ".stp"} or not path.is_file():
                    raise ValueError(f"Expected an existing STEP file for {label} {surface}.")
                step_paths[label, surface] = path
                step_sources[label, surface] = "supplied_step"
        with TemporaryDirectory(prefix=".report-", dir=study_dir) as temp:
            staged = Path(temp) / report_name
            report = write_evidence(
                staged, definitions, FIXTURE, step_paths=step_paths, fallback=fallback,
                step_sources=step_sources, rebuild_of=rebuild_of,
            )
            overlay = destination / Path(report["overlay"]).relative_to(staged)
            report["overlay"] = str(overlay)
            for surfaces in report["revisions"].values():
                for evidence in surfaces.values():
                    readback = evidence.get("step_readback", {})
                    if readback.get("source") == "occt_fallback":
                        relative_path = Path(readback["path"]).relative_to(staged)
                        readback["path"] = str(destination / relative_path)
            (staged / "evidence.json").write_text(
                json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8",
            )
            staged.rename(destination)
        return {
            "evidence": (destination / "evidence.json").as_uri(),
            "overlay": overlay.as_uri(),
            "files": sorted(path.name for path in destination.iterdir()),
            "measurements": report,
        }

    @server.tool()
    def phase0_spar(
        report_name: str,
        wing_step: str,
        tip_piece_length_mm: float = 250.,
        min_tip_engagement_mm: float = 20.,
        max_print_height_mm: float = 250.,
    ) -> dict[str, Any]:
        """Solve and verify a straight spar channel in a Wing B STEP solid.

        Give a STEP path relative to the study directory. Set the span length
        of the tip piece, the minimum full-diameter spar depth inside that
        piece, and the maximum spanwise height allowed for any print piece.
        The solver determines the other cuts. Install the cad extra.
        """
        settings = {
            "tip_piece_length_mm": tip_piece_length_mm,
            "min_tip_engagement_mm": min_tip_engagement_mm,
            "max_print_height_mm": max_print_height_mm,
        }
        for name, value in settings.items():
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive.")
        if tip_piece_length_mm > max_print_height_mm:
            raise ValueError("tip_piece_length_mm exceeds max_print_height_mm.")
        destination = study_dir / checked_label(report_name)
        if destination.exists():
            raise ValueError(f"Report directory {destination} already exists.")
        relative_path = Path(wing_step)
        if relative_path.is_absolute():
            raise ValueError("Give a relative STEP path for wing_step.")
        source = (study_dir / relative_path).resolve()
        if not source.is_relative_to(study_dir):
            raise ValueError("wing_step must stay inside the study directory.")
        if source.suffix.lower() not in {".step", ".stp"} or not source.is_file():
            raise ValueError("Expected an existing STEP file for wing_step.")

        from planeopt.spar import solve_spar

        with TemporaryDirectory(prefix=".spar-", dir=study_dir) as temp:
            staged = Path(temp) / report_name
            report = solve_spar(source, staged, **settings)
            staged.rename(destination)
        return {
            "evidence": (destination / "result.json").as_uri(),
            "files": {path.name: path.as_uri() for path in sorted(destination.iterdir())
                      if path.is_file()},
            "measurements": report,
        }

    return server


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-dir", required=True, type=Path, help="Directory for study outputs")
    parser.add_argument(
        "--converter-source", type=Path,
        help="Root of a pinned Airfoil Converter checkout for fixture verification",
    )
    args = parser.parse_args()
    study_dir = args.study_dir.resolve()
    if study_dir.is_relative_to(FIXTURE.resolve()):
        parser.error("--study-dir must be outside the baseline fixture")
    if args.converter_source is not None:
        source = args.converter_source.resolve() / "src"
        if not (source / "airfoil_converter" / "wing_build.py").is_file():
            parser.error("--converter-source must contain src/airfoil_converter/wing_build.py")
        sys.path.insert(0, str(source))
    study_dir.mkdir(parents=True, exist_ok=True)
    make_server(study_dir, args.converter_source).run()


if __name__ == "__main__":
    main()
