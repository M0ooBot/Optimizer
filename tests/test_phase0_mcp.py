"""Exercise the Phase 0 study through its stdio MCP entry point."""

import asyncio
import json
from pathlib import Path
import sys

import pytest


mcp = pytest.importorskip("mcp")
from mcp.client.stdio import stdio_client


REPO = Path(__file__).resolve().parents[1]


def test_study_through_stdio(tmp_path):
    async def exercise():
        server = mcp.StdioServerParameters(
            command=sys.executable,
            args=[str(REPO / "tools" / "phase0_mcp.py"), "--study-dir", str(tmp_path)],
        )
        async with stdio_client(server) as streams, mcp.ClientSession(*streams) as session:
            await session.initialize()
            tools = await session.list_tools()
            assert {tool.name for tool in tools.tools} >= {
                "phase0_reference", "phase0_shape", "phase0_report", "phase0_spar",
            }

            async def call(name, **arguments):
                result = await session.call_tool(name, arguments)
                assert not result.isError, result.content
                return result.structuredContent

            reference = await call("phase0_reference")
            assert reference["definition"]["wing"]["projected_span_mm"] == 600
            assert "midspan_twist_peak_deg" in reference["section_twist"]
            twist = await call(
                "phase0_shape", revision="reversing",
                wing={"midspan_twist_peak_deg": 4, "tip_twist_deg": -4}, lock_tail=True,
            )
            twist_definition = json.loads(Path(twist["definition"]).read_text(encoding="utf-8"))
            assert twist_definition["wing"]["midspan_twist_peak_deg"] == 4
            assert twist_definition["wing"]["tip_twist_deg"] == -4
            initial = await call(
                "phase0_shape", revision="swept",
                wing={
                    "sweep_mm": 105, "tip_chord_mm": 90, "rise_mm": 45,
                    "tip_twist_deg": -3, "section_orientation": "front_normal",
                    "hinge": {"start": 0.2, "end": 0.49},
                },
                brief="Sweep and taper the wing and raise its outer span.",
            )
            follow_up = await call(
                "phase0_shape", revision="rise-later", parent="swept",
                wing={"rise_start": 0.7}, lock_tip_height=True, lock_tail=True,
                brief="Start the rise farther outboard and preserve the tip height and tail.",
            )
            before = json.loads(Path(initial["definition"]).read_text(encoding="utf-8"))
            after = json.loads(Path(follow_up["definition"]).read_text(encoding="utf-8"))
            assert after["tail"] == before["tail"]
            assert after["wing"].pop("tip_span_trim_mm") > 0
            assert before["wing"].pop("tip_span_trim_mm") > 0
            assert after["wing"] == {**before["wing"], "rise_start": 0.7}
            for result in [initial, follow_up]:
                generated = result["generated_geometry_measurements"]["wing"]
                assert generated["tip_le_xyz_mm"][1] == 45
                assert generated["sampled_projected_skin_span_mm"] == pytest.approx(600, abs=0.000001)
                assert generated["minimum_sampled_hinge_skin_margin_mm"] > 0
            request = json.loads(Path(follow_up["request"]).read_text(encoding="utf-8"))
            assert request["brief"].startswith("Start the rise farther outboard")

            full_span = await call(
                "phase0_shape", revision="distributed", parent="rise-later",
                wing={"span_curve": "full_span", "sweep_mm": 85, "rake_mm": 20},
                lock_tip_aft=True, lock_tip_height=True, lock_tail=True,
            )
            full_definition = json.loads(Path(full_span["definition"]).read_text(encoding="utf-8"))
            assert full_definition["wing"]["rise_start"] is None
            assert full_definition["wing"]["rake_start"] is None
            assert full_definition["tail"] == after["tail"]
            full_measurement = full_span["generated_geometry_measurements"]["wing"]
            assert full_measurement["tip_le_xyz_mm"][1:] == [45, -105]
            assert full_measurement["sampled_projected_skin_span_mm"] == pytest.approx(600, abs=0.000001)

            conflict = await session.call_tool("phase0_shape", {
                "revision": "invalid", "parent": "rise-later",
                "wing": {"projected_span_mm": 605},
            })
            assert conflict.isError
            assert "wing.projected_span_mm" in conflict.content[0].text
            assert not (tmp_path / "invalid").exists()
            conflict = await session.call_tool("phase0_shape", {
                "revision": "invalid-aft", "parent": "distributed",
                "wing": {"sweep_mm": 90}, "lock_tip_aft": True,
            })
            assert conflict.isError
            assert "wing.sweep_mm + wing.rake_mm" in conflict.content[0].text
            report = await call(
                "phase0_report", report_name="comparison",
                revisions=["baseline", "swept", "rise-later", "distributed"],
            )
            assert report["measurements"]["baseline_revision"] == "baseline"
            generated = report["measurements"]["revisions"]["distributed"]["wing"]["generated"]
            assert generated["span_curve"] == "full_span"
            assert generated["rise_onset_fraction"] is None
            assert report["overlay"].startswith("file:")
            assert (tmp_path / "comparison" / "evidence.json").is_file()

    asyncio.run(exercise())


def test_spar_rejects_invalid_settings_and_step_paths(tmp_path):
    async def exercise():
        server = mcp.StdioServerParameters(
            command=sys.executable,
            args=[str(REPO / "tools" / "phase0_mcp.py"), "--study-dir", str(tmp_path)],
        )
        async with stdio_client(server) as streams, mcp.ClientSession(*streams) as session:
            await session.initialize()
            tools = await session.list_tools()
            properties = next(tool for tool in tools.tools
                              if tool.name == "phase0_spar").inputSchema["properties"]
            assert properties["tip_piece_length_mm"]["default"] == 250.
            assert properties["min_tip_engagement_mm"]["default"] == 20.
            assert properties["max_print_height_mm"]["default"] == 250.
            for options, message in [
                ({"tip_piece_length_mm": 251.}, "exceeds max_print_height_mm"),
                ({"min_tip_engagement_mm": -1.}, "finite and positive"),
                ({"min_tip_engagement_mm": 250.}, "existing STEP file"),
                ({"max_print_height_mm": 0.}, "finite and positive"),
                ({"wing_step": str(tmp_path / "external.step")}, "relative STEP path"),
                ({"wing_step": "../external.step"}, "inside the study directory"),
                ({}, "existing STEP file"),
            ]:
                arguments = {"report_name": "spar-study", "wing_step": "wing.step", **options}
                result = await session.call_tool("phase0_spar", arguments)
                assert result.isError
                assert message in result.content[0].text
            assert list(tmp_path.iterdir()) == []

    asyncio.run(exercise())
