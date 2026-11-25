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
- **Leading and trailing edges are drawn in the side view** (seen from outboard: up against aft). Each edge's aft position is a function of its own height above its own start. It leaves the attachment in Wing B's direction and curvature (needed for a seamless join), eases into a one-way turn over `edge_ease` of its height, and keeps turning to the tip. The turn fades toward the tip as `(1 − y / height) ** edge_fade`. The setting is the edge's overall sweep (`sweep_deg`, `te_sweep_deg`): the straight line from start to tip, aft positive. The tip direction is not forced, so nothing curves back in at the tip. The chord is the distance between the edges. Alternatively `tip_chord_mm` sets the tip chord and the trailing edge's overall sweep follows. See [Side view](#side-view).
- **Twist** blends from Wing B's washout into 2° toe-out.
- **Sweep, taper, twist, morph and plane rotation advance with the turn**, not with distance: each is blended by how much of the 7°-to-75° turn the leading edge has already made. Each section then changes while the skin is already rising. Advancing them with distance instead made the winglet sag below Wing B's line in the front view before the bend lifted it (see [Front view](#front-view)).
- **Section morph.** Each section is split into a camber line and a thickness distribution. The SD7037 camber fades to zero and its thickness blends into a symmetric 10% section (NACA 0010 thickness). The trailing-edge thickness is put back as a ramp fixed at 0.4 mm, so it never scales with the chord. The Airfoil Converter ([d77bbdf](https://github.com/Trying-to-Fly/Airfoil-Converter/tree/d77bbdf45a8c0aa92efbc2a3d3007cfe4f904105)) cuts both end sections' blunt trailing edges with keep-chord, as the cap-radius parents do.

Solving order: the winglet's length meets the height; the attachment station meets the 600 mm span, with the tip extended 12 mm along its tangent as an allowance for the cap. The attachment may not come within 10 mm of the aileron hinge's 510 mm end; if it would have to, the bend fraction is solved instead (see [Bend distribution](#bend-distribution)).

## Bend distribution

**Chosen: centred (35% / 35%), the default in `WingletDefinition`.** The ramps are continuous settings, not modes. The flat top's width is `1 − ramp_in − ramp_out`: zero is a pure bell, larger values approach an even turn. The peak sits at `ramp_in + width / 2` along the winglet. A third setting, `bend_fraction`, is the share of the winglet the bell occupies. Below 1, the turn finishes early and a straight leg runs to the tip (still tapering and sweeping, because both edges are drawn). That saves room sideways at the price of a tighter bend. With `bend_fraction=None` the law keeps the bend over the whole winglet when that fits the span. Otherwise it starts the winglet at the hinge clearance and picks the longest bend that fits.

## Configuration study

One setting at a time around the baseline: **60 mm, 75° cant, overall sweep 35° LE / 10° TE aft, 2° toe-out, centred bend**. The baseline is 60 mm, not 65 mm, because at 65 mm the winglet already starts within 0.4 mm of the hinge clearance and almost any change pushed it past the span. [Open the viewer](stage1/viewer.html) (download it and open it locally) and pick a configuration from the grouped list. The camera stays put as you switch, *Baseline ghost* overlays the baseline, and *From outboard* is the side view the edge rules are drawn in. Each group's front view, plan view and curvature are overlaid in [stage1/configs/](stage1/configs/) (`group-*.png`), and every configuration has its own report and views there.

Every configuration fits. The tall and low-cant ones fit only because the bend finishes early. Edge angles are seen from outboard, from vertical, aft positive, at the start and at the tip:

| Configuration | Starts at (mm) | Bend | Tightest bend radius (mm) | Leading edge | Trailing edge | Tip chord (mm) | Smallest LE radius (mm) | |
| --- | ---: | --- | ---: | --- | --- | ---: | ---: | --- |
| Bend distribution (ramps in / out): Even 20 / 20 | 525.8 | whole | 67.0 | 56° → 19° | -7° → 17° | 58.1 | 0.69 |  |
| Bend distribution (ramps in / out): Centred 35 / 35 (baseline) | 525.9 | whole | 55.3 | 56° → 19° | -7° → 17° | 58.1 | 0.69 |  |
| Bend distribution (ramps in / out): Bell 50 / 50 | 525.9 | whole | 43.1 | 56° → 19° | -7° → 17° | 58.1 | 0.69 |  |
| Bend distribution (ramps in / out): Early peak 30 / 60 | 536.5 | whole | 43.1 | 56° → 19° | -7° → 16° | 55.8 | 0.67 |  |
| Bend distribution (ramps in / out): Late peak 50 / 25 | 520.0 | 95% | 53.5 | 56° → 19° | -8° → 17° | 59.4 | 0.71 |  |
| Trailing-edge sweep, overall: 10° forward | 526.6 | whole | 55.3 | 56° → 19° | -7° → -12° | 38.4 | 0.46 |  |
| Trailing-edge sweep, overall: 0° | 526.3 | whole | 55.3 | 56° → 19° | -7° → 3° | 48.2 | 0.58 |  |
| Trailing-edge sweep, overall: 15° aft | 525.7 | whole | 55.3 | 56° → 19° | -7° → 23° | 63.3 | 0.73 | **chord grows near the tip** |
| Trailing-edge sweep, overall: 20° aft | 525.4 | whole | 55.3 | 56° → 19° | -7° → 30° | 68.7 | 0.73 | **chord grows near the tip** |
| Trailing-edge sweep, overall: 25° aft | 525.2 | whole | 55.3 | 56° → 19° | -7° → 36° | 74.4 | 0.73 | **chord grows near the tip** |
| Leading-edge sweep, overall: 20° | 525.1 | whole | 55.3 | 56° → -8° | -7° → 17° | 78.5 | 0.73 | **chord grows near the tip** |
| Leading-edge sweep, overall: 25° | 525.3 | whole | 55.3 | 56° → 0° | -7° → 17° | 72.3 | 0.73 | **chord grows near the tip** |
| Leading-edge sweep, overall: 35° | 525.9 | whole | 55.3 | 56° → 19° | -7° → 17° | 58.1 | 0.69 |  |
| Leading-edge sweep, overall: 45° | 526.6 | whole | 55.3 | 56° → 38° | -7° → 17° | 39.9 | 0.57 |  |
| Leading-edge sweep, overall: 50° | 527.1 | whole | 55.3 | 56° → 46° | -7° → 17° | 28.3 | 0.47 |  |
| Height: 45 mm | 542.9 | whole | 41.5 | 57° → 19° | -6° → 16° | 62.3 | 0.71 |  |
| Height: 50 mm | 537.3 | whole | 46.1 | 56° → 19° | -7° → 16° | 60.9 | 0.72 |  |
| Height: 55 mm | 531.6 | whole | 50.7 | 56° → 19° | -7° → 17° | 59.5 | 0.71 |  |
| Height: 65 mm | 520.2 | whole | 59.9 | 56° → 19° | -8° → 17° | 56.7 | 0.68 |  |
| Height: 70 mm | 520.0 | 94% | 58.8 | 56° → 19° | -8° → 17° | 54.1 | 0.65 |  |
| Cant at the tip: 60° | 520.0 | 77% | 54.2 | 56° → 19° | -8° → 17° | 59.5 | 0.70 |  |
| Cant at the tip: 65° | 520.0 | 89% | 58.0 | 56° → 19° | -8° → 17° | 59.4 | 0.70 |  |
| Cant at the tip: 70° | 520.0 | 98% | 60.2 | 56° → 19° | -8° → 17° | 59.4 | 0.71 |  |
| Cant at the tip: 75° | 525.9 | whole | 55.3 | 56° → 19° | -7° → 17° | 58.1 | 0.69 |  |
| Cant at the tip: 80° | 532.7 | whole | 49.8 | 56° → 19° | -7° → 16° | 56.6 | 0.68 |  |
| Cant at the tip: 85° | 539.1 | whole | 45.4 | 56° → 19° | -7° → 16° | 55.2 | 0.66 |  |
| Toe-out: 0° | 525.3 | whole | 55.3 | 56° → 19° | -7° → 17° | 58.1 | 0.69 |  |
| Toe-out: 2° | 525.9 | whole | 55.3 | 56° → 19° | -7° → 17° | 58.1 | 0.69 |  |
| Toe-out: 4° | 526.3 | whole | 55.3 | 56° → 19° | -7° → 17° | 58.2 | 0.69 |  |
| Toe-out: 6° | 526.6 | whole | 55.3 | 56° → 19° | -7° → 17° | 58.4 | 0.70 |  |
| Combinations: 70 mm and 70° cant | 520.0 | 85% | 56.2 | 56° → 19° | -8° → 17° | 54.2 | 0.64 |  |

What the study shows:

- **The overall sweeps set the tip chord.** A trailing edge from 10° forward to 25° aft gives 38–74 mm. A leading edge from 20° to 50° gives 79–28 mm.
- **The chord grows again near the tip** when the trailing edge ends up swept further back than the leading edge. That happens at trailing-edge sweep, overall 15° aft, trailing-edge sweep, overall 20° aft, trailing-edge sweep, overall 25° aft, leading-edge sweep, overall 20°, leading-edge sweep, overall 25°: the edges' overall sweeps are too far apart in the wrong direction. These configurations are flagged, not hidden.
- **Height and cant set the bend.** A taller or less steep winglet needs more room sideways. It starts farther inboard until it reaches the hinge clearance, then the bend shortens. The tightest radius runs from 41.5 mm (45 mm tall) to 60 mm (70° cant).
- **Toe-out changes almost nothing in these views.**
- **The leading-edge rule binds at small tip chords.** At 50° leading-edge sweep the tip chord is 28 mm. The nose radius is 0.34 mm streamwise and passes only because sweep makes the cross-section across the edge rounder (0.48 mm).

## Side view

Seen from outboard, the first edge laws drew an S. The leading edge went 56° → 34° → 37° and the trailing edge swept back past its target before curving back in (to 11.4°, then 10.4° at the tip). Each edge was forced to arrive at an exact sweep at the tip, and its sweep was defined against distance in the front view, which the side view compresses as the winglet stands up. The rule from review:

1. Every guide curve bends one way only, from the attachment to the tip. Leading edges bow forward and trailing edges bow aft, so the side-view outline is convex.
2. Nothing forces the tip direction, so no guide curves back in.
3. The guides share a similar overall sweep.

At the attachment each edge must still leave in Wing B's direction. Seen from outboard that is 56° for the leading edge (Wing B's raked tip seen nearly end-on) and about 7° forward for the trailing edge. So the leading edge steepens and the trailing edge swings aft, each one way.

Now both edges turn one way, with no counter-turn in any configuration (measured where the view sees them, past the stretch that is end-on). The baseline's leading edge goes 56° → 19° and its trailing edge −7° → 17°. With `edge_fade` = 2 the turn happens mostly low on the winglet, and both edges end nearly parallel. With a constant turn rate (`edge_fade` = 0), averaging 35° from a 56° start ended the leading edge swept forward.

**Still open: guides between the edges.** Chord-line guides at 25%, 50% and 75% of the chord counter-turn by up to 2.0° on the baseline and 6.0° at 60° cant. Skin guides are similar. Their aft positions are blends of the two edges, which turn in opposite senses. Their heights also carry the section's twist, camber and thickness, which swing from vertical to sideways as the section plane rotates. Both bend the guide's side-view path. Each configuration's report records the counter-turn of every guide under measured.side_view.

## Checks, every configuration

| Check | Result |
| --- | --- |
| Projected half-span with the 12 mm cap allowance | 600 mm |
| Height, tip cant, overall sweeps, toe-out | as requested, measured on the generated sections |
| Trailing edge | 0.4 mm at every station |
| Join to Wing B (LE and three skin points) | largest second-derivative jump 7.2e-8 per mm against values of about 3e-4 to 9e-4 |
| Edges turn one way seen from outboard; bend curvature never negative; dihedral rises throughout; every section has positive thickness | pass |
| Chord decreases throughout | pass, except the 5 configurations flagged above |
| Leading-edge radius ≥ 0.4 mm | pass; smallest 0.46 mm across the edge |
| Droop below Wing B's front-view tangent | upper outline at most 0.012 mm, skin lines at most 0.062 mm (see [Front view](#front-view)) |
| Guides between the edges, counter-turn seen from outboard | up to 6.0° (open; see [Side view](#side-view)) |

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

[Airfoil Converter tools this work needs](CONVERTER_TOOLS.md)
