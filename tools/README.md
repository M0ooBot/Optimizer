Offline pipelines (EXECUTION_PLAN.md section 5) — arrive at M1:

- `ingest_props.py` — APC performance tables -> data/props/ CSVs + smooth CT/CP fits
- `fit_profile.py`  — slicer results for scaled sections -> ConstructionProfile constants

## SD7037 Phase 0 study

Launch the study server as a stdio subprocess from a local MCP client. Set the
client's working directory to the Optimizer checkout. The server stores each
revision under the study directory, including its definition, generated curve
files and replayable tool arguments.

```powershell
uv run --locked --extra phase0 --extra cad python tools/phase0_mcp.py --study-dir C:/Studies/SD7037
```

Add `--converter-source C:/path/to/Airfoil-Converter` for fixture verification.
Use [the pinned Converter revision](https://github.com/Trying-to-Fly/Airfoil-Converter/tree/d77bbdf45a8c0aa92efbc2a3d3007cfe4f904105).
The local adapter uses the [MCP stdio transport](https://modelcontextprotocol.io/specification/2025-06-18/basic/transports).

1. Call `phase0_reference` to read the dimensions, axes and field names. Call
   `phase0_verify_converter_fixture` to regenerate the original wing fixture
   with the pinned Converter's Wing builder and compare the curve files.
   The study applies its sweep, rake, rise and twist to the saved Converter
   profiles and skin anchors.
2. Call `phase0_shape` with a new `revision` label and partial `wing` and `tail`
   parameters. Values are in millimetres and degrees. Span and chord fractions
   are in the range 0–1. The server validates the span and attachment locks and
   writes a complete geometry definition and curves to the new revision folder.
   Pass the designer's original wording in `brief` to record it with the tool
   arguments. The response includes the smallest sampled hinge skin margin
   in each surface. A negative value locates an axis that leaves the skin.
3. For a follow-up moving the start of the rise outboard, pass the previous
   revision as `parent`, set only `wing.rise_start`, and set both
   `lock_tip_height` and `lock_tail` to true. A shape call with the saved
   revision as `parent` and no geometry changes regenerates its curves.
   The projected span lock uses the skin envelope. The server calibrates the
   leading-edge guide for this span when section planes tilt. The tip height
   lock uses the leading edge, and the tool also reports quarter-chord tip height.
   To distribute the bend and sweep from the root, set
   `wing.span_curve` to `"full_span"` and omit `rise_start` and
   `rake_start`. Set `lock_tip_aft` to preserve the leading-edge tip's
   aft coordinate while changing how sweep and rake are distributed.
   The optional `wing.full_span_linear_rise_fraction` blends the quadratic
   rise at 0 with linear rise at 1 while preserving the leading-edge tip's
   height and aft coordinate. The default is 0. A nonzero blend with nonzero
   rise requires parallel sections to preserve the fixed root mounting plane.
   To prescribe section twist, set `wing.tip_twist_deg` and optionally
   `wing.midspan_twist_peak_deg`. Each section rotates about its local leading
   edge. The root stays at zero degrees and positive incidence lowers the
   trailing edge. With a midpoint angle, two quintic ramps meet at half-span
   with zero slope and curvature. Setting the midpoint angle to null restores
   the linear root-to-tip law.
4. Call `phase0_report` with `report_name` and an ordered `revisions` list
   such as `["baseline", "swept", "rise-later"]`. Use `fallback: true`
   for offline Open CASCADE STEP lofts and sampled readback. CAD comparisons
   use the supplied original baseline wing STEP as their reference. To measure a
   supplied export, pass `step_files` as a mapping such as
   `{"swept": {"wing": "swept/wing.step", "tail": "swept/tail.step"}}`.
   The files must be inside the study directory. Pass `rebuild_of` as a mapping
   such as `{"swept-rebuild": "swept"}` to compare a regenerated revision
   with its original. The tool returns the measurements and the PNG overlay
   location.
5. Call `phase0_spar` with a new `report_name` and a `wing_step` path inside
   the study directory. Set `tip_piece_length_mm` to the outboard piece's
   spanwise length, `min_tip_engagement_mm` to the minimum full-diameter
   spar depth measured from that piece's joint, and `max_print_height_mm`
   to the spanwise height limit for every piece. The defaults are 250 mm,
   20 mm, and 250 mm. The tool solves for a straight 6.4 mm bore with at
   least 1 mm of exterior wall, measures the joint geometry, and returns
   the evidence. Passing results include the bored STEP pieces and assembly.

Example `phase0_shape` arguments for a local edit of an existing `swept`
revision:

```json
{
  "revision": "rise-later",
  "parent": "swept",
  "wing": {"rise_start": 0.7},
  "lock_tip_height": true,
  "lock_tail": true,
  "brief": "Start the rise farther outboard and preserve the tip height and tail."
}
```

### SolidWorks handoff

SW-Agent's advertised write scope includes approved STEP export of an eligible
saved independent part and selected edits to existing closed parts. Its
`newDocument` route is inspection-only. An agent needs a supported bridge
operation that creates a study part, imports the generated profile and guide
curves as editable features, and builds and rebuilds the loft before it can
complete a native CAD round trip. The XYZ curve files and offline STEP solids are
available for geometric comparison while that bridge operation is absent.

If a native study part is available, start with `sw_status`, inspect the
part with `sw_feature_tree`, and use `scripts_list` to check for an eligible
packaged export. The guarded `sw_run_snippet` STEP route requires one eligible
saved part in `targets`, effect `writes-new-files`, an existing `outputFolder`
inside the study directory and one STEP filename in `params.fileName`. Use
the bridge approvals and poll `sw_change` when the bridge returns a pending
ID. Give the measured export to `phase0_report` as a supplied STEP.
