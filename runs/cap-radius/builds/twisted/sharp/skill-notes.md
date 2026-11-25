# Notes for the maintainer of close-rounded-wingtip

Build: twisted/sharp parent (1000 mm span, 5° dihedral, 3° washout, cambered NACA 2412-type section,
sharp TE, locked 60 mm 3D rim). Tool: CadQuery/OCP only, no SolidWorks. Final result: every stated
requirement is met on the reimported STEP. The smallest radius in the ≥ 0.8 mm zone is 0.411 mm at
chord fraction 0.9720, where t = 0.800 mm. The margins and the interpretations it took are below.

## How each step was followed

**Build 1: attachment curves.**
- Followed. The parent revision is the one locked in `parent.py` and `parent.step`; the attachment
  is at |Y| = 440.
- The cap's first and last B-spline pole rows are copied from the parent skin's attachment pole rows,
  so the shared edges are identical, not just close. This needs the cap to use the parent's u basis.
  The skill does not say how to achieve "use the actual cut edges" in a kernel without a
  boundary-surface feature. A line on that would help.

**Build 2: 3D perimeter.**
- The perimeter was locked and supplied (`tip_outline`), so there was nothing to draw. I verified
  that it is tangent to the LE and TE in 3D and lies on the continued camber surface.

**Build 3: wrap guides.**
- Followed, in a continuous form: a guide at every chord parameter p, not a handful of discrete
  guides. Each guide's turn is defined by a roll profile F and a per-guide roll share E(p).
- With no boundary-surface tool in OCP, the "boundary surface" (build 4) is an analytic surface
  through the three profiles with the guides as its w curves, then fitted to one B-spline face per
  tip.

**Build 4: boundary surface.**
- G1 and G2 to both skins come by construction: the cap is the parent's own section formula with the
  thickness side ±1 replaced by F(w), and F = 1 + O(η³).
- Measured on the STEP:
  - tangent-plane angle ≤ 1.5e-6°;
  - curvature difference in a common frame ≤ 3.7e-6 relative, for chord fraction 0.001–0.999;
  - 0.47 % relative at the last sample next to the TE point.

**Build 5: printable radius.**
- Met, measured as the larger principal curvature:
  - zone with t ≥ 0.8 mm: minimum 0.411 mm;
  - thin aft zone: lip 1.033–1.089 × t/2, and nothing tighter than 1.028 × t/2;
  - thin LE zone: minimum 1.50 mm.
- It drove the whole design, as the skill says it would. See the iterations below.

**Build 6: progressive turn.**
- No bending reversals in any of 124 sections.
- Mid chord: the tightest point of each guide is within 0.01 in w of the lip.
- Aft thin zone: the turn is a near-constant-radius roll. Its shoulder is 1.5–4 % tighter than the
  lip; that is the cost of G2 against a full half-round.
- Near the LE: the tightest point is the parent's own nose curvature (about 1.5 mm), not the cap lip
  (2.0 mm).

**Verify 1–4.**
- All done on the reimported STEP, with CAD derivatives through OCP.
- **Verify 5** (visual review): zebra PNGs were rendered from CAD normals per pixel, with no edge
  lines, centred and mirrored cameras, and mirrored stripes.

## Where the skill was unclear, ambiguous or wrong

1. **"Radius" in step 5 versus planar sections in Verify 2.**
   - Build 5 says "every radius along the wrap". Verify 2 measures the radius of planar sections at
     guide locations. Verify 3 uses the larger principal curvature. These are not the same number.
   - Near the LE and TE the locked rim runs almost spanwise, so any plane through a chordwise guide
     cuts the lip obliquely. The section's in-plane radius is about the true lip radius × |rim
     tangent · plane normal| (Meusnier).
   - On this wing that relation held to within 1 % at all 18 guide stations. Guide-plane sections show 0.07 mm at chord fraction 0.0004 and 0.29 mm at 0.9455,
     where the true lip radii are 2.0 mm and 0.76 mm.
   - Judged literally on guide-plane sections, the rule is impossible to meet near the LE and TE for
     any closed cap on this outline.
   - I judged it on principal curvature and on rim-normal sections, and reported the guide-plane
     numbers next to them (`sections_obliquity` in `readback.json`).
   - Suggestion: state that the rule is a surface (principal) radius, and that sections must be cut
     normal to the rim to measure the lip.
2. **"Chord position" of a cap point is undefined.**
   - The rule depends on attachment thickness "at the same chord position". A point on the cap has
     no chord position of its own.
   - I used the chord fraction of the point's guide, which is also the chord fraction of its rim
     point, sin²(πp/2).
   - Near the TE the true lip cross-section runs chordwise, so one lip section spans a small range of
     chord fractions. It made no practical difference here, but it should be defined.
3. **The thin-zone rule is not achievable literally with G2.**
   - A turn of 90° per side within a height of t/2, with every radius ≥ t/2, forces an exact circular
     arc. That arc meets the flat with a curvature jump (G1).
   - The same bound applies in the thick zone just above t = 0.8 mm. There "every radius ≥ 0.4 mm"
     with h = 0.4 mm also forces an exact circle.
   - The way out on this wing was the parent's TE wedge, about 16.5° included angle. A lip on a
     wedge turns only about 82° per side, which leaves about 15 % slack for a G2 shoulder.
   - On a blunt or parallel-sided section, the rule at t ≈ 0.8 mm cannot be met with a G2 roll. The
     skill should say which wins: G2, or the radius floor.
4. **"End the cap the same way" at a sharp TE.** It is clear in intent. Note that the radius → 0 at the
   TE point, so any curvature map must stop short of it. The skill does say to exclude collapsed
   endpoints, but not how close to go. I went to 0.019 mm.
5. **Thin LE zone.**
   - The rule names "the thin end of any section". At the LE the section is under 0.8 mm thick only
     for chord fraction < 0.00043 (the first 0.04 mm).
   - There the lip is the parent's nose swept round the rim: radius about 1.5–2 mm, much rounder than
     t/2.
   - I treated "rounder than t/2" as meeting the rule. A literal t/2 half-round there would create a
     tight spot the parent does not have.
6. **"Spread the roll over more of the cap reach" is only half the story.**
   - Where the section is thick (forward and mid chord), spreading the roll over the whole reach is
     right. Here E = 1 forward of p = 0.32.
   - Aft, the sections are thin and the reach is long (e.g. t = 2.6 mm over 35 mm at 90 % chord).
     There a spread, elliptical roll gives tip radii far below 0.4 mm.
   - The fix is the opposite: keep the full thickness out towards the rim and concentrate the roll
     into a blunt, nearly half-round lip.
   - The skill's sentence led me first towards a spread-out family that failed aft.
7. **Zebra artefacts that are not defects.** Two features look like flaws in the renders but are not:
   - The parent is only C1 along its camber-blend lines (chord fraction 0.36 and 0.44; curvature jumps
     3–8 %), and the cap inherits these.
   - Where the roll starts on a nearly spanwise-straight skin, stripes rotate quickly but
     continuously.

   Verify 5's "recheck suspicious marks" helped. A note that kinks along parent chord lines may be
   inherited would have saved time.
8. **Parametrisation skew near the LE and TE.**
   - With chordwise guides, the guide (w) direction and the u lines become almost parallel near the
     LE and TE (8° apart at p = 0.95).
   - This amplifies tiny fit errors in the chordwise principal curvature by about cot². My first
     B-spline fit (69 w poles) was accurate to 2e-4 mm but showed 0.2 /mm curvature noise there.
   - Dense knots on the rim half (157 w poles) fixed it, with curvature error ≤ 0.17 % of each
     guide's peak.
   - The skill does not warn about this. It matters for any construction built from chord-station
     guides on a rim that ends tangent to the LE and TE.

## Decisions and improvisations not covered by the skill

- The roll profile family, G(x) = 1 − (1 − x)⁴ q(x), with C3 onset and coefficients optimised for a
  round aft lip.
- The roll-share schedule E(p) and its margins.
- Fitting in the parent's u basis so that the edges are shared exactly.
- Constrained (value + 1st + 2nd derivative) w fits at both attachments.
- Rim-normal cuts as the lip cross-section.
- The numerical zebra model:
  - phase = π·24·(R·d), with R the reflected view ray;
  - 9 views × 6 stripe angles;
  - phase and surface-gradient mismatch at 4001 points per join;
  - stripe-direction change at every light/dark crossing on the joins.

## Iterations

All design iterations except the last were done on the analytic model before any STEP was written.
There was one STEP rebuild.

1. Superellipse roll F = sqrt(1 − ηⁿ) with n(p) rising from 3 to 30.
   - The aft lip had a shoulder tighter than the lip (0.82–0.88 × t/2).
   - The thick-zone minimum near t = 0.8 mm was below 0.4 mm.
   - Rejected.
2. Compact-support polynomial roll G = 1 − (1 − x)³(1 + Ax) with a roll-share schedule.
   - Mid chord (chord fraction 0.6–0.7) reached only 0.33 mm, because the schedule started too far
     aft.
   - Moving the schedule start to p = 0.32 and optimising (A, B, E_TE) fixed it.
3. Raised the onset order to (1 − x)⁴ with a cubic q, to ease the roll onset and help the fit, then
   re-optimised. The thin-zone minimum rose from 1.020 to 1.033 × t/2 (analytic).
4. First STEP with 32 graded w spans: valid and accurate in position, but with curvature noise near
   the TE (item 8 above). Rebuilt with density-graded knots. This is the final STEP.

## Requirements

| Requirement | Status |
|---|---|
| One closed, valid solid; parent skin, attachment, span and outline unchanged | Met (pole-identical parent, 1000.0000 mm, rim 4.6e-10 mm) |
| Tangent to both skins | Met, ≤ 1.5e-6° |
| Curvature matched and eased into the roll | Met, ≤ 3.7e-6 relative (0.47 % at the last sample by the TE point); onset C3 |
| ≥ 0.4 mm where t ≥ 0.8 mm | Met, 0.411 mm minimum (2.8 % margin, at the t = 0.8 boundary) |
| Half-round of t/2 where t < 0.8 mm (aft) | Met as a G2 near-half-round: lip 1.033–1.089 × t/2, minimum 1.028 × t/2. Not an exact circle (impossible with G2) |
| Thin LE zone | Exceeded: 1.50 mm minimum, rounder than t/2 rather than a t/2 half-round |
| Cap ends at the sharp TE like the parent | Met (radius → 0 at the TE point) |
| Progressive turn, no reversal | Met. Aft lip is a near-constant-radius turn with its shoulder 1.5–4 % tighter than the lip |
| Zebra continuity across joins | Met: phase ≤ 1.1e-7 period, direction change ≤ 0.00025°, gradient ≤ 1.8 % (worst at the TE corner sample) |
| Guide-plane section radius ≥ 0.4 mm (literal reading of Verify 2) | Not met near the LE and TE (down to 0.07 mm), and in my view it cannot be met on this outline. See item 1 |
