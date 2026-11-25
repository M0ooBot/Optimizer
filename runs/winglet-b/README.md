# Wing B blended winglet

The twisted Wing B ([definition](../wing-b-spar/shape/constrained/redistributed/definition.json)) continued into a blended winglet. The 600 mm projected half-span includes the winglet. The trailing edge is blunt at 0.4 mm throughout, so the 0.4 mm minimum radius applies everywhere, including the winglet's leading edge.

| Stage | Status |
| --- | --- |
| 1. Law and outline: edge curves, section planes, checks on the law | **In tuning.** Bend: centred, chosen 2026-11-24. Edges redrawn in the side view. Configuration study done ([stage1](stage1/)) |
| 2. Surface with the tip open; locked blunt parent; STEP readback | Next |
| 3. Close the winglet tip with the [rounded wingtip skill](../../.agents/skills/close-rounded-wingtip/SKILL.md) | |
| 4. Whole-surface curvature map and zebra review | |
| 5. A follow-up edit, e.g. move the bend's peak and keep the height | |

## Decisions

- **Shared code for geometry, a skill for judgment** (2026-11-24). The winglet's laws, join continuity, solvers and checks live in one tested module, [phase0_winglet](../../src/planeopt/phase0_winglet.py). The skill to come will teach an agent which settings to tune, how to read the checks and the viewer, and how to extend the laws. A new winglet is a new `WingletDefinition`, not new geometry code. Insights found in review (draw both edges; change sections with the turn) become code and checks once, so no later agent has to rediscover them. Open: make the laws pluggable, so an agent with a new idea (another curvature shape, a split tip) writes only that piece and inherits the join and the checks. A parent adapter is also needed, so the law can continue a wing imported from CAD and not only a Phase 0 definition.
- **Suggestion for the rounded cap** (not yet discussed with whoever owns the cap work). The cap trials went the other way: each agent wrote its own geometry-specific construction from the [rounded wingtip skill](../../.agents/skills/close-rounded-wingtip/SKILL.md) ([cap trials](../cap-trials/README.md)). The cap could try this format too: one shared cap generator (attachment curves from the parent, the 3D perimeter, wrap guides, the thickness-gated 0.4 mm rule), the [independent cap audit](https://github.com/M0ooBot/Optimizer/tree/cap-radius-trials/runs/cap-radius/audit) as its standard check, and the skill reduced to judgment. Repository issues are disabled, so this note stands in for one.

## How the winglet is defined

The winglet is not a second body blended onto the wing. Wing B is generated from rules: a leading-edge path, a chord, a twist about the leading edge and the SD7037 section in spanwise-parallel planes ([phase0_geometry](../../src/planeopt/phase0_geometry.py)). [phase0_winglet](../../src/planeopt/phase0_winglet.py) continues each rule past the attachment, starting from the parent's value and its first and second derivatives along the leading edge, so the skin is curvature-continuous with Wing B by construction. The checks below measure it.

The whole winglet is the transition. From the attachment to the tip:

- **Bend.** In the front view the leading edge turns from Wing B's local dihedral (7.3°) to the 75° cant at the tip. Its curvature follows a flat-topped bell. It ramps up from Wing B's curvature over the first `ramp_in` of the winglet's length, holds, and ramps down to zero over the last `ramp_out`. Two settings cover the shape: their sum sets how spread out the turn is (0.5 / 0.5 is a single peak; short ramps spread it evenly) and their ratio moves the peak toward the winglet root or the tip. The peak curvature follows from the turn, so a narrower bell means a tighter bend.
- **Section planes** start spanwise-parallel like Wing B's and rotate until they are normal to the leading edge in the front view at the tip.
- **Leading and trailing edges are drawn in the winglet's planform**: looking straight at its face, normal to the plane of its tip, along its span direction against aft. Each edge's aft position is a function of its own distance u along that span direction. It leaves the attachment in Wing B's direction and curvature (needed for a seamless join), then bends aft progressively, like a parabola: its aft slope only grows. The bend eases in over `edge_ease` of the edge's length and builds as `(u / length) ** edge_progression` (0 is a parabola; larger values start gentler and bend harder near the tip). The setting is the edge's overall sweep (`sweep_deg`, `te_sweep_deg`): the straight line from start to tip, aft positive. The tip direction is not forced. The chord is the distance between the edges. Alternatively `tip_chord_mm` sets the tip chord and the trailing edge's overall sweep follows. See [Planform](#planform).
- **Twist** blends from Wing B's washout into 2° toe-out.
- **Sweep, taper, twist, morph and plane rotation advance with the turn**, not with distance: each is blended by how much of the 7°-to-75° turn the leading edge has already made. Each section then changes while the skin is already rising. Advancing them with distance instead made the winglet sag below Wing B's line in the front view before the bend lifted it (see [Front view](#front-view)).
- **Section morph.** Each section is split into a camber line and a thickness distribution. The SD7037 camber fades to zero and its thickness blends into a symmetric 10% section (NACA 0010 thickness). The trailing-edge thickness is put back as a ramp fixed at 0.4 mm, so it never scales with the chord. The Airfoil Converter ([d77bbdf](https://github.com/Trying-to-Fly/Airfoil-Converter/tree/d77bbdf45a8c0aa92efbc2a3d3007cfe4f904105)) cuts both end sections' blunt trailing edges with keep-chord, as the cap-radius parents do.

Solving order: the winglet's length meets the height; the attachment station meets the 600 mm span, with the tip extended 12 mm along its tangent as an allowance for the cap. The attachment may not come within 10 mm of the aileron hinge's 510 mm end; if it would have to, the bend fraction is solved instead (see [Bend distribution](#bend-distribution)).

## Bend distribution

**Chosen: centred (35% / 35%), the default in `WingletDefinition`.** The ramps are continuous settings, not modes. The flat top's width is `1 − ramp_in − ramp_out`: zero is a pure bell, larger values approach an even turn. The peak sits at `ramp_in + width / 2` along the winglet. A third setting, `bend_fraction`, is the share of the winglet the bell occupies. Below 1, the turn finishes early and a straight leg runs to the tip (still tapering and sweeping, because both edges are drawn). That saves room sideways at the price of a tighter bend. With `bend_fraction=None` the law keeps the bend over the whole winglet when that fits the span. Otherwise it starts the winglet at the hinge clearance and picks the longest bend that fits.

## Tuning

Tuning now works on one configuration, the baseline, until it is right; then the study is rebuilt from it. `build.py` with no flags builds only the baseline into [stage1/tuning/](stage1/tuning/) and `viewer.py` writes [stage1/tuning.html](stage1/tuning.html) (about 30 seconds together). `--all` builds and views the whole study (about 20 minutes).

## Configuration study

One setting at a time around the baseline: **60 mm, 75° cant, overall sweep 35° LE / 10° TE aft, 2° toe-out, centred bend**. It was last built with the planform edge rule; it is not rebuilt during tuning. The baseline is 60 mm, not 65 mm, because at 65 mm the winglet already starts within 0.4 mm of the hinge clearance and almost any change pushed it past the span. [Open the study viewer](stage1/viewer.html) (download it and open it locally) and pick a configuration from the grouped list. *Winglet planform* is the view the edges are drawn in, and *Baseline ghost* overlays the baseline. Each group's front view, plan view and curvature are overlaid in [stage1/configs/](stage1/configs/) (`group-*.png`).

Not built: leading-edge sweep, overall 50° (the tip chord would fall below the airfoil family's 18 mm). Edge angles are in the winglet's planform, against its span direction, aft positive, at the start and at the tip:

| Configuration | Starts at (mm) | Bend | Tightest bend radius (mm) | Leading edge | Trailing edge | Tip chord (mm) | Smallest LE radius (mm) | Every guide bends aft progressively |
| --- | ---: | --- | ---: | --- | --- | ---: | ---: | --- |
| Bend distribution (ramps in / out): Even 20 / 20 | 526.2 | whole | 67.0 | 26° → 43° | -3° → 23° | 49.7 | 0.73 | yes |
| Bend distribution (ramps in / out): Centred 35 / 35 (baseline) | 526.2 | whole | 55.3 | 26° → 43° | -3° → 23° | 49.7 | 0.73 | yes |
| Bend distribution (ramps in / out): Bell 50 / 50 | 526.2 | whole | 43.1 | 26° → 43° | -3° → 23° | 49.7 | 0.72 | **no** (2.4°) |
| Bend distribution (ramps in / out): Early peak 30 / 60 | 536.7 | whole | 43.1 | 27° → 43° | -2° → 23° | 48.9 | 0.71 | **no** (2.1°) |
| Bend distribution (ramps in / out): Late peak 50 / 25 | 520.0 | 96% | 53.8 | 26° → 43° | -3° → 24° | 50.2 | 0.73 | yes |
| Trailing-edge sweep, overall: 10° forward | 527.2 | whole | 55.3 | 27° → 43° | -3° → -18° | 24.5 | 0.39 | **no** (15.7°) |
| Trailing-edge sweep, overall: 0° | 526.7 | whole | 55.3 | 27° → 43° | -3° → 3° | 37.1 | 0.58 | **no** (1.2°) |
| Trailing-edge sweep, overall: 15° aft | 525.9 | whole | 55.3 | 26° → 43° | -3° → 32° | 56.3 | 0.73 | yes |
| Trailing-edge sweep, overall: 20° aft | 525.7 | whole | 55.3 | 26° → 43° | -3° → 40° | 63.2 | 0.73 | yes |
| Trailing-edge sweep, overall: 25° aft | 525.4 | whole | 55.3 | 26° → 43° | -3° → 47° | 70.5 | 0.73 | yes, **chord grows near the tip** |
| Leading-edge sweep, overall: 20° | 525.2 | whole | 55.3 | 26° → 11° | -3° → 23° | 75.4 | 0.73 | **no** (15.3°), **chord grows near the tip** |
| Leading-edge sweep, overall: 25° | 525.5 | whole | 55.3 | 26° → 23° | -3° → 23° | 67.6 | 0.73 | **no** (4.7°), **chord grows near the tip** |
| Leading-edge sweep, overall: 35° | 526.2 | whole | 55.3 | 26° → 43° | -3° → 23° | 49.7 | 0.73 | yes |
| Leading-edge sweep, overall: 45° | 527.1 | whole | 55.3 | 27° → 58° | -3° → 23° | 26.8 | 0.58 | yes |
| Height: 45 mm | 543.2 | whole | 41.5 | 27° → 43° | -2° → 23° | 56.0 | 0.70 | **no** (1.9°) |
| Height: 50 mm | 537.5 | whole | 46.1 | 27° → 43° | -2° → 23° | 53.9 | 0.71 | **no** (1.0°) |
| Height: 55 mm | 531.9 | whole | 50.7 | 27° → 43° | -3° → 23° | 51.8 | 0.72 | **no** (0.5°) |
| Height: 65 mm | 520.5 | whole | 59.9 | 26° → 43° | -3° → 23° | 47.6 | 0.73 | yes |
| Height: 70 mm | 520.0 | 95% | 59.1 | 26° → 43° | -3° → 23° | 45.1 | 0.70 | yes |
| Cant at the tip: 60° | 520.0 | 78% | 55.0 | 17° → 49° | -2° → 23° | 44.7 | 0.74 | yes |
| Cant at the tip: 65° | 520.0 | 89% | 58.7 | 19° → 48° | -2° → 23° | 46.1 | 0.74 | yes |
| Cant at the tip: 70° | 520.0 | 99% | 60.7 | 22° → 46° | -3° → 23° | 48.0 | 0.74 | yes |
| Cant at the tip: 75° | 526.2 | whole | 55.3 | 26° → 43° | -3° → 23° | 49.7 | 0.73 | yes |
| Cant at the tip: 80° | 532.9 | whole | 49.8 | 33° → 37° | -3° → 24° | 51.3 | 0.72 | **no** (2.8°) |
| Cant at the tip: 85° | 539.2 | whole | 45.4 | 42° → 24° | -4° → 24° | 52.7 | 0.66 | **no** (18.1°) |
| Toe-out: 0° | 525.7 | whole | 55.3 | 26° → 43° | -3° → 23° | 49.8 | 0.73 | yes |
| Toe-out: 2° | 526.2 | whole | 55.3 | 26° → 43° | -3° → 23° | 49.7 | 0.73 | yes |
| Toe-out: 4° | 526.5 | whole | 55.3 | 27° → 43° | -3° → 23° | 49.7 | 0.73 | **no** (0.3°) |
| Toe-out: 6° | 526.8 | whole | 55.3 | 27° → 43° | -3° → 23° | 49.8 | 0.72 | **no** (0.5°) |
| Combinations: 70 mm and 70° cant | 520.0 | 86% | 56.7 | 22° → 46° | -3° → 23° | 43.0 | 0.71 | yes |

What the study shows:

- **An edge can only bend aft if its overall sweep is above its starting angle.** In this view the leading edge leaves Wing B at about 26° and the trailing edge at about 3° forward. Asking for less (leading edge 20° or 25°, trailing edge 10° forward), or a cant that raises the starting angle (85° starts the leading edge at 42°), makes that edge bend forward instead. Those configurations fail the check.
- **Tighter bends wobble the skin guides.** The bell and early-peak bends, short winglets (45–55 mm) and 80° cant keep the edges and the chord-line guides right. Skin guides near the surface counter-turn by 0.3–2.8°, from the section's thickness swinging as it rotates through a tight bend.
- **The overall sweeps set the tip chord.** A trailing edge from 0° to 25° aft gives 37–71 mm. A 45° leading edge gives 27 mm.

## Planform

Two earlier versions drew S-shaped edges. The first forced each edge to a set sweep at the tip. The second was drawn seen from outboard, where Wing B at the join is seen nearly end-on and its raked leading edge appears to leave at 56°. From there no edge can bend aft progressively, so it had to steepen first, which read as a bend at one spot. The rule from review: every guide bends one way, aft and progressively (aft slope and its rate both positive, like a parabola), with the tip direction free and similar overall sweeps.

Drawn in the winglet's own planform, the leading edge leaves at about 26° and the trailing edge at about 3° forward. Both bend aft progressively to their overall sweeps. On the baseline the leading edge goes 26° → 43° and the trailing edge −3° → 23°; the guides between them do the same, with no counter-turn on the chord line and at most 0.11° on the skin. Each configuration's report records every guide under measured.planform.

## Checks, every configuration

| Check | Result |
| --- | --- |
| Projected half-span with the 12 mm cap allowance | 600 mm |
| Height, tip cant, overall sweeps, toe-out | as requested, measured on the generated sections |
| Trailing edge | 0.4 mm at every station |
| Join to Wing B (LE and three skin points) | largest second-derivative jump 7.2e-8 per mm against values of about 3e-4 to 9e-4 |
| Every guide bends aft progressively in the winglet planform (edges, chord-line guides, skin guides within 0.25°) | baseline and most configurations; see the table for the rest |
| Bend curvature never negative; dihedral rises throughout; every section has positive thickness | pass |
| Chord decreases throughout | pass, except where flagged in the table |
| Leading-edge radius ≥ 0.4 mm | pass; smallest 0.46 mm across the edge |
| Droop below Wing B's front-view tangent | upper outline at most 0.012 mm, skin lines at most 0.062 mm (see [Front view](#front-view)) |

## Front view

A review of the bell variant saw the winglet first droop below Wing B and then curve up. It did. With the sections changing in step with distance along the winglet, the shrinking chord and the fading camber lowered the upper outline before the bend had lifted the skin. The upper outline fell up to 0.12 mm below Wing B's tangent line and skin lines near the trailing edge up to 0.22 mm, about 20 mm along.

Tying those changes to the turn instead fixes it. What remains is Wing B's own shape carried through the join for curvature continuity. Wing B's washout makes its trailing edge curve down in the front view (curvature −8.3e-4 per mm; it would fall 0.09 mm below its tangent over 15 mm), and its upper outline too (−2.7e-4 per mm). The winglet lifts both sooner than Wing B alone would. Each configuration's report records these measurements under measured.front_view_sag.

## Findings outside the winglet

- **SD7037's resampled nose has a micron-scale tight spot.** A spline through the Converter's 321-point SD7037 has a radius of about 0.00086 chord right at the foremost point (0.08 mm at 87.5 mm chord), against about 0.008 chord (0.7 mm) over the nose as a whole. The points there are unevenly spaced. The bump is a few microns high, far below what a printer resolves, but a curvature map of Wing B's own CAD would flag it against the 0.4 mm rule. This study resamples every section on one cosine grid, which smooths it out.

## Regenerate

From the repository root, with the pinned Converter extracted somewhere outside the checkout:

```sh
git -C <Airfoil-Converter> archive d77bbdf src | tar -x -C <dir>
uv run --locked --extra phase0 --extra cad python runs/winglet-b/stage1/build.py --converter-src <dir>/src
uv run --locked --extra phase0 --extra cad python runs/winglet-b/stage1/viewer.py --converter-src <dir>/src
uv run --locked --extra phase0 --extra cad python -m pytest tests/test_phase0_winglet.py
```

The tests use a NACA stand-in for the Converter so they run without it.

[Airfoil Converter tools this work needs](CONVERTER_TOOLS.md), now served by the [Converter MCP server](CONVERTER_TOOLS.md#converter-mcp-server)
