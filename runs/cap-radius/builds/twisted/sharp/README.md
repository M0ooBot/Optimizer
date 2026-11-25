# Twisted / sharp parent: rounded tip caps

`wing.step` is the full 1000 mm wing: the locked parent from
`runs/cap-radius/parents/twisted/sharp/parent.step` with its two flat attachment faces replaced by a
rounded cap at each tip. It is one closed, valid solid with 6 B-spline faces: the 4 parent skin faces,
unchanged, and one cap face per tip. Units are mm, axes as in `spec.json` (X aft, Y span, +Y right
tip, Z up).

- File: `wing.step` (AP214, 16.1 MB). SHA-256 `9db6d83f6f6e4e965525d4946c4651f35f161e180c068d3e03536f35ad94973f`
- Volume: 1 514 991.8 mm³ (BRepGProp, 1e-9 tolerance). The parent alone is 1 448 878.3 mm³.
- Projected span: 1000.0000 mm. Bounding box: X −0.014…189.881, Z −7.628…50.760.

## Construction

The method follows `.agents/skills/close-rounded-wingtip/SKILL.md`. `cap_model.py` defines it
analytically.

- **Profiles.** The cap starts on the upper attachment cut edge (w = 0), runs through the locked 3D rim
  (w = 1) and ends on the lower cut edge (w = 2). The cut edges are the parent skin at |Y| = 440. The
  rim is `parent.tip_outline`: it leaves the leading-edge point tangent to the leading edge and
  reaches the trailing-edge point tangent to the trailing edge.
- **Guides.** There is a wrap guide at every chord parameter p, which is the parent's own native chord
  parameter (u = s = p). Each guide runs from the upper skin point at p, through the rim point p, to
  the lower skin point at p. The cap point is

  `S(p, w) = place(Y, airfoil(p, sigma))`, with Y = 440 + 60 sin(πp)·(1 − v²), v = |1 − w|, and sigma = ±F.

  `place` and `airfoil` are the parent's own functions, with taper, twist, dihedral and sweep
  continued outboard. The thickness side (±1) is replaced by a signed thickness fraction sigma. At
  sigma = ±1 the formula gives the parent skin exactly; at sigma = 0 it gives the camber surface, which
  is where the locked rim lies.
- **Roll profile.** With x = v²/E(p), the thickness fraction is F = sqrt(G(x)) for x < 1 and 1 beyond,
  where G(x) = 1 − (1 − x)⁴(1 + 1.4631x + 1.3569x² + 1.9858x³).
  - Near the rim, F ∝ v, so the lip is round and the parametrisation stays smooth through it.
  - At x = 1, F meets the unreduced thickness with C3 continuity.
  - At the attachment, the cap continues the parent skin to second order (exact G2 in the analytic
    model).
  - G's coefficients were optimised so that, in the thin aft region, the lip is a near-circular
    half-round with no tighter shoulder.
- **Roll share E(p).** E(p) is the fraction of the reach that the roll uses. It is 1 forward of p = 0.32,
  so the roll is spread over the whole reach. Between p = 0.32 and 0.88 it falls on a quintic smoothstep
  to 0.0837. That value makes the aft lip radius about 1.03–1.09 × half the attachment thickness. The
  schedule was set from the analytic principal curvatures so that the minimum radius is ≥ 0.4 mm
  wherever the attachment section is ≥ 0.8 mm thick, with about 35 % margin over the mid chord.
- **Fit.** Each cap is one degree 5×5 B-spline face with 511 × 157 poles (`build_wing.py`):
  - In u it uses the parent skin's own basis (same knots and degree). Its first and last pole rows are
    copied from the parent's attachment rows, so the cap and the parent share the attachment edges
    exactly.
  - In w it is a constrained least-squares fit, with value, first and second derivative imposed at
    both attachments. The knots are dense on the rim half.
  - The u direction is interpolated at the Greville abscissae.
  - The LE and TE pole columns collapse to the parent's LE and TE points.
  - Fit error against the analytic cap is 5.9e-5 mm at most. The worst principal-curvature error is
    0.17 % of the guide's peak curvature, and the guide minimum radius is within 1.8e-5 (relative).
- **Left tip.** The poles are mirrored in XZ. The faces are sewn at 1e-6 mm (0 free edges) and made
  into a solid.

## Measured on the reimported STEP (`verify_wing.py` → `readback.json`)

| Check | Result |
|---|---|
| Solid | 1 solid, 1 shell, 6 faces, BRepCheck valid, shell closed, 0 free edges, outward normals verified |
| Parent preserved | parent face poles identical to parent.step (max diff 0.0); parent skin vs parent.py 8.9e-8 mm (18 492 samples) |
| Span / attachment / outline | span 1000.0000 mm; cap edges at \|Y\| = 440.000; rim vs locked outline 4.6e-10 mm; rim is the outermost point of every guide |
| Join separation (4 joins × 4001 points) | ≤ 4.3e-14 mm |
| Join tangent-plane angle | ≤ 1.5e-6° |
| Join curvature, common frame (II along, across, twist) | max abs difference 4.1e-6 /mm; relative to local parent curvature ≤ 3.7e-6 for chord fraction 0.001–0.999, 0.47 % at the last sample next to the TE point |
| Zebra across joins (9 views × 6 stripe angles, 6794 light/dark crossings) | phase mismatch ≤ 1.1e-7 period; stripe-direction change at crossings ≤ 0.00025°; stripe-gradient mismatch ≤ 1.8 % (worst at the last sample, chord fraction 0.999997, 0.02 mm from the TE point) |
| **Printable radius, map (713 × 1569 grid per cap + refinement)** | **thick zone (t ≥ 0.8 mm): smallest radius 0.411 mm at chord fraction 0.9720, where t = 0.800 mm (point 187.60, ±459.55, 45.72), on the shoulder just inboard of the lip; passes 0.4 mm** |
| Thin aft zone (chord fraction > 0.9720) | lip radius 1.033–1.089 × t/2; smallest radius anywhere 1.028 × t/2 |
| Thin LE zone (chord fraction < 0.00043) | smallest radius 1.50 mm (t/2 ≤ 0.40 mm) |
| Rim-normal sections (62) | thick zone smallest 0.416 mm (chord fraction 0.9718); thin zone ≥ 1.034 × t/2; no bending reversals on the cap |
| Guide-plane sections (36 + 26 between) | no bending reversals on the cap; in-plane radii fall below 0.4 mm near the LE and TE (see below) |

Mid-chord radii are about 0.55–0.75 mm, and 1.5–1.9 mm near the leading edge (from the map's
`by_chord_position` table).

## Decisions and caveats

- **Rule judged on principal curvature.** The printable rule is judged on the larger principal
  curvature of the cap, at the chord position of the guide: the rim point's chord fraction,
  sin²(πp/2), with t taken from `parent.section_thickness` at Y = 440.
- **Guide-plane sections understate the radius near the LE and TE.** Near the LE and TE the locked rim
  runs almost spanwise, so a plane through a guide cuts the lip obliquely. Its in-plane radius is the
  true lip radius scaled by about |rim tangent · plane normal| (Meusnier); `readback.json`
  `sections_obliquity` lists the two side by side; at all 18 guide stations their ratio matches
  |rim tangent · plane normal| to within 1 %. For example, the guide plane at p = 0.85 shows
  0.29 mm where the true lip radius is 0.76 mm. The true lip cross-section is the rim-normal cut, and
  those cuts agree with the curvature map.
- **Thin LE zone.** For chord fraction below 0.00043 the cap is much rounder than the thin rule asks
  (1.5 mm), because there the lip is the parent's own nose (nose radius about 1.6 mm) swept around
  the outline. It is not a t/2 half-round.
- **Thin aft zone.** The lip is a G2 near-half-round slightly larger than t/2. The parent's trailing-edge wedge
  (about 16.5° included angle) leaves room for it. A G2 roll cannot be an exact circle of radius t/2.
- **Inherited C1 lines.** The parent is only C1 along its camber-blend lines (chord fraction 0.36 and
  0.44): its chordwise curvature jumps 3–8 % there. The cap inherits these lines, and they show as
  small stripe kinks in the zebra renders. The parent is locked, so they were left as they are.
- **Visible roll onset.** Aft of 32 % chord, the roll starts on the curve w = 1 − sqrt(E(p)). Before
  that curve the skin is almost straight spanwise; after it, curvature begins smoothly (C3). On a
  spanwise-straight skin this turns the zebra stripes quickly but continuously, which looks like a
  bend in `zebra/*_tip_top.png`. It is present in the analytic design too, not a fit artefact.
- **Collapsed corners excluded.** The cap collapses to the LE and TE points. The curvature map stops
  at p = 1e-4 and 1 − 1e-4, about 0.019 mm from those points.

## Files and how to run (from the project root, with `.venv/bin/python`)

- `cap_model.py`: the analytic cap, the design parameters and the rule. It imports the parent's `parent.py`.
- `build_wing.py`: fits the caps, sews the solid and writes `wing.step` and `build.json`. Takes about 40 s.
- `step_geom.py`: STEP reading and CAD-derivative helpers.
- `verify_wing.py`: all checks on the reimported STEP; writes `readback.json`. Takes about 12 min.
- `render_zebra.py`: zebra images in `zebra/`, shaded from CAD normals at each pixel's (u, v), with
  no edge lines, centred and mirrored cameras, and mirrored stripes. Takes about 2 min.
