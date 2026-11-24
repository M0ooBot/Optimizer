# Wing B handoff

Start a new study from the [revised Wing B definition](../../runs/wing-b-spar/shape/constrained/redistributed/definition.json) and the matching [unbored half-wing STEP](../../runs/wing-b-spar/shape/constrained/redistributed/wing.step). This is the 600 mm SD7037 half-wing with a 200 mm root chord, 75 mm tip chord, 60 mm tip rise and installed twist of +4° at midspan and −4° at the tip. Its CAD coordinates put the root leading edge at the origin, with outboard along −X, up along +Y and forward along +Z.

For a new MCP study, copy the definition to `<study-dir>/wing-b/definition.json` and the STEP to `<study-dir>/wing-b/wing.step`. Use `wing-b` as the parent of a new `phase0_shape` revision. The [MCP instructions](../../tools/README.md#sd7037-phase-0-study) cover setup, rebuilds, readback and adjustable spar settings. Codex and Claude share the [wing-shaping skill](../../.agents/skills/shape-wing-edges/SKILL.md).

## Reference geometry and measurements

The [full-wing preview](../../runs/wing-b-spar/longest-spar/boundary-verified/full-wing.png) mirrors the modeled half at the root. Download and open the [offline 3D viewer](../../runs/wing-b-spar/longest-spar/boundary-verified/full-wing.html) to inspect it.

On this geometry the longest verified straight spar measures 430.4405 mm from root center to cap center. The 6.4 mm finished channel leaves 1.00010 mm of exterior wall at its closest point. The [verified measurements and axis coordinates](../../runs/wing-b-spar/longest-spar/boundary-verified/result.json), [bored three-piece half-wing](../../runs/wing-b-spar/longest-spar/boundary-verified/print-assembly.step) and [upper-bound certificate](../../runs/wing-b-spar/longest-spar/upper-bound.json) document this CAD limit. The bound excludes lengths of 431 mm or more for this shape. The test used a 250 mm tip piece, a 250 mm print-height limit and a minimum engagement of 20 mm. New geometry needs its own fit check.

## Possible next trial

The [rounded wingtip case study](../../docs/studies/WING_B_ROUNDED_TIP.md) uses the earlier untwisted Wing B. It includes the rolled-cap construction, two STEP revisions, a full-wing viewer and a numerical zebra comparison.

Explore whether an agent can add a winglet to this half-wing using the offline geometry workflow. First establish the intended winglet shape and which dimensions to hold fixed, especially whether the 600 mm projected half-span includes the winglet. Compare the new tip against the saved STEP and inspect the leading edge, trailing edge and surface through the junction from both the front and the top. SolidWorks work is deferred.
