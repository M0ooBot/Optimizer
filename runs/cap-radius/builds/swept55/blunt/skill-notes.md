# Skill notes: close-rounded-wingtip on swept55 / blunt

These notes are for the skill's maintainer. They cover what I did for each step, where the skill was unclear or wrong, what I had to invent, and what I met.

## Steps followed

### Build

1. **Attachment curves.**
   - Used the parent's cut edges at Y = ±390, unchanged. The poles of the retained parent faces are identical to `parent.step`.
   - On a blunt TE, the two cut edges end at the TE corners (Z = ±0.2). Between them is the flat-face segment, which the three-profile recipe does not cover. See ambiguity 1.
2. **Perimeter.**
   - Used the locked outline at Z = 0. It is tangent to the LE and TE lines at its ends, by construction.
   - The only deviation is in the last 0.5 mm before the TE end, where the lip is nudged onto the TE plane: ≤ 0.16 µm measured, ≤ 1.2 µm by design. This lets the half-round apex sit exactly in that plane.
3. **Spanwise guides.**
   - The guides follow the sweep: they run along lines of constant chord fraction of the continued planform. This is what the outline's own definition implies, and near the LE it makes the cap continue the parent's nose.
   - Near the TE the guides had to fan forward (see "Improvised").
   - The guides are defined analytically (`cap_model.py`). The lip radius is set per guide from a target.
4. **Boundary surface.**
   - There is no SolidWorks; OCCT has no equivalent of a curvature-matched boundary surface. So I defined the surface analytically (profiles plus guides plus a height law) and fitted one quintic B-spline per patch.
   - The fit is constrained to:
     - match the parent's first and second cross-derivatives at the cut, which gives G1 and G2,
     - put every lip pole at Z = 0, with the next pole row directly above, so the mirrored lower patch shares the tangent plane,
     - put every pole of the TE-end guide in the TE plane.
   - The upper and lower patches are separate; they are G1 at the lip by construction and G2 by symmetry.
5. **Printable radius rule.** Designed in from the start: the lip-radius target is part of the guide definition. Measured on the STEP, see below.
6. **Progressive turn.**
   - Every guide profile is a flat run (exactly the parent), then a C2 ramp, then a stretched arc that ends vertical at the lip. So the curvature rises toward the lip.
   - The lip-normal sections show no tight bend before the lip and no reversal, apart from numerical-level ones at R ≥ 29 mm.

### Verify

All checks were run on the reimported STEP (`verify.py`):

- topology and validity, parent preservation, and the locked dimensions;
- separation, angle and curvature in a common frame along all joins (four cap/parent joins, two lips, two half-disc/TE-face joins, and the designed crease);
- 48 sections: lip-normal, station, and sweep-following;
- a principal-curvature map: 900 × 260 grid plus Nelder–Mead refinement, judged against the rule by the plan chord fraction of each point;
- numerical zebra across the joins;
- renders, in `figures/`.

## Requirements: met, partly met, not met

| Requirement | Status |
|---|---|
| One closed valid solid; parent, span and outline unchanged | Met. Outline deviation ≤ 0.16 µm, at the TE end. |
| Tangency to both parent skins | Met. ≤ 3.6e-5° (≤ 5e-4° at the collapsed LE point). |
| Curvature matched at the attachment (smooth reflection) | Met more than 3 mm from the TE corner: zebra cross-gradient jump ≤ 3.5e-4 °/mm against a 0.28 °/mm median. **Partly met near the TE corner**: up to 2.9 °/mm at 0.5–3 mm and 24.5 °/mm at 0.26 mm from it. There the cap is only 0.2–1 mm long and the fit's curvature error dominates. At the corner vertex itself, G2 is impossible: the half-round starts with curvature 5 mm⁻¹ where the parent edge is straight. |
| Lip: finite radius, upper and lower share a tangent plane | Met. |
| t ≥ 0.8: every radius ≥ 0.4 | **Met: min 0.412 mm**, at chord fraction 0.97996, t = 0.800. |
| t < 0.8 at the LE: lip half-round with R = t/2 | Met in spirit, not literally. The lip there is the parent's nose (R ≈ 1.9 mm, ≥ 4.8 × t/2). Making it "t/2" would mean a lip sharper than the parent's own LE. See ambiguity 2. |
| t < 0.8 at the TE: lip half-round with R = t/2; on a blunt TE, the half-round closes the TE face with radius t_TE/2 | Met for the lip: R/(t/2) = 1.07–1.17 more than 0.1 mm beyond the apex. The TE face ends in a semicircle of radius exactly 0.2. Within 0.05 mm of the apex the lip is flatter than a half-round (the flat face fades into it). |
| "Every turn printable" inside the thin TE zone | **Not met within ~1 mm of each TE corner.** Interior radius:<br>- 0.067 mm at the half-round edge, 0.04 mm from the corner (0.34 × t/2);<br>- 0.083 beyond 0.1 mm;<br>- 0.107 beyond 0.2 mm;<br>- 0.153 beyond 0.3 mm;<br>- 0.186 beyond 0.5 mm;<br>- ≥ t/2 beyond 1.0 mm.<br>See ambiguity 3. |
| Zebra continuity across both joins | Met away from the TE corner. The stripes kink within about 1.5 mm of it: direction change up to 138° at stripe boundaries there, against ≤ 0.02° elsewhere. The renders also show the stripes bending within about 3 mm outboard of the join (taper let-out; see item 8). |

## Where the skill was unclear, ambiguous or wrong

1. **The blunt TE does not fit the three-profile recipe.**
   - The upper and lower cut edges end 0.4 mm apart. The perimeter ends at the outline's end point, which is the middle of the TE segment. Something must close the segment, and the skill only says "the half-round closes the end of the flat trailing-edge face, so its radius there is half the trailing-edge thickness".
   - My reading: the flat face continues as a coplanar half-disc of radius 0.2. The cap patches meet its semicircular edge, and the parent's TE-corner crease runs around that edge, fading to 0° at the apex where the lip takes over.
   - Any G1 alternative has a singular point at each TE-corner vertex. The cap's normal there would have to equal both the skin's normal and the TE face's normal.
   - **Suggest:** state this construction explicitly, or say which alternative is intended.
2. **"Thickness at the attachment at the same chord position" is ambiguous on a swept outline that is tangent to the LE and TE.**
   - Near the corners, the roll section normal to the lip crosses many chord positions. A sweep-following section meets the lip at a grazing angle.
   - At the LE, lip points with chord fraction < 0.001 sit on the parent's nose, whose radius (0.85 streamwise, about 1.48 normal to the LE) is far above t/2.
   - I judged every surface point by its own plan chord fraction in the continued planform.
   - **Suggest:** define the chord position for the rule (the lip point's? the section's?), and say the LE corner inherits the nose.
3. **The rule cannot hold at the end of a blunt TE.**
   - The parent's TE corners are 84° creases that stop at the attachment. A rounded cap must dissolve them along the half-round's 0.31 mm edge: the surface normal turns about 81° along that arc. That forces a twist of roughly 5 rad/mm plus the half-round's own 5 mm⁻¹.
   - The principal radius next to the edge therefore cannot reach t/2 = 0.2. My estimate of the ideal is about 0.12–0.17. I reached 0.067 at the edge and ≥ 0.1 beyond 0.2 mm, after trying several fields.
   - The skill should say what is acceptable here: a crease continuing the parent's (as I did), a minimum radius, or a different TE treatment such as a longer flat-face taper.
4. **"Spreading the roll over more of the cap reach" helps the interior, not the lip.**
   - For an elliptic-type roll the lip radius is H²/W, so a wider roll makes the lip *sharper*.
   - What makes the 0.4 mm rule feasible at t ≈ 0.8:
     - concentrating the turn near the lip (lip R = H at best, with a flat run),
     - and using the section's forward thickening (the TE wedge) through an arc tangent to the parent's slope, which gives up to about 1.17 H.
   - With a G2 flat start, a quarter-circle of radius exactly H is impossible: the curvature must overshoot somewhere. So the threshold band t ∈ [0.8, 0.84] needs that wedge.
   - **Suggest:** state that the best achievable radius is about the local half-thickness (more with a thickening section), and that the margin at the threshold is thin.
5. **"Record the radius at the lip in that same section" can mislead.**
   - Sweep-following and station sections cross the lip obliquely near the LE, TE and tip apex. Their in-plane lip radius equals the surface radius × sin θ, which gave 0.02–0.17 mm where the surface radius is 0.2–1.9.
   - **Suggest:** require lip-normal sections for lip radii; use sweep cuts only for reversals.
6. **Collapsed endpoints.** The LE corner of each patch is collapsed (a degenerate edge). Excluded as the skill says; the grid reaches 1e-12 mm of it and radii are ≥ 1.34 mm beyond 0.05 mm. The half-round apex is not collapsed, but the lip radius there tends to infinity by design.
7. **Reversals.**
   - The deliberate taper let-out (thickness kept along the cap rather than continuing the chord taper) shows as a concave bend along the sweep with R ≈ 850–910 mm in the first 3 mm.
   - The station cut at Y = 390.5 shows a concave R of 2.5 mm in the TE-corner transition.
   - An intermediate design (scaling the roll by the height at the roll top to gain wedge thickness) met the 0.4 rule but produced a concave band with R ≈ 1.1 mm. It was caught by the sections and replaced by the tangent-arc form. The section check earned its place.

8. **The zebra numbers at the join do not capture "easing".**
   - The join is G2: the phase-gradient jump is ≤ 3.5e-4 °/mm more than 3 mm from the corners. Yet the renders show the stripes bending within about 3 mm outboard of the join.
   - The cause is my taper let-out, which is G2 at the cut but changes the spanwise slope over only 3 mm.
   - Step 4's "ease that curvature into the roll" has no measurable criterion.
   - **Suggest:** add one, such as the rate of change of the stripe gradient over the first few mm, or of curvature along the guides.

## Improvised (not covered by the skill)

- **Taper let-out.** The parent's half-thickness is continued outboard with the chord taper saturating over 3 mm (`r_eff = 390 + 3 tanh((r − 390)/3)`). This keeps G2 at the cut without thinning the section under the lip.
- **TE fan.** Within about 60 mm of the TE corner, the guides start up to 2 mm (plus 0.06 × the lip distance) forward of the constant-chord-fraction start. Pure sweep-following guides there meet the lip at 0–5°, and the roll collapses to micrometres: it effectively continues the TE crease.
- **TE-corner field.** A parabolic-coordinate roll field (`te_field`) that equals the half-round on the TE line inboard of the apex and the target bullnose outboard of it. Its coefficient fixes how the crease fades around the half-round (80.9° → 1.2°), and the fan slope is matched to it so the patch stays regular at the apex corner.
- **Gap form.** `z = Zp − H_lip·Δ(τ)`, where Δ is an arc tangent to the parent's own slope at the roll top. This keeps the roll convex and lets the TE wedge raise the lip radius to about 1.17 × t/2.
- **Fit constraints.**
  - G1 and G2 at the cut, lip poles at Z = 0 with the next row straight above them, the TE-end guide in the TE plane.
  - The patch parameter is graded toward the LE corner and the TE apex, but not too far: grading to nanometre scale made the fit's curvature error at the TE corner much worse.
  - An extra U-derivative constraint at the half-round guide made the lip ring and was dropped.

## Iterations

About 9 design iterations and 6 fitting iterations:

1. Sweep-following guides with a flat-start ellipse. Fine at mid-chord; micrometre radii near the TE corner, because the guides grazed the lip.
2. TE fan plus a curved bulge to even out the crease: about 0.05 at the half-round, plus a plan-mapping fold, a sign bug, and a stalled attachment parameter.
3. The parabolic TE-corner field. A τ normalisation factor-of-2 error was fixed. Result: about 0.07–0.09 next to the half-round, 0.2 on the lip.
4. Ramp and lip-target tuning. The thick-zone minimum was 0.364–0.386 at the threshold band: the flat-start overshoot.
5. Roll height taken from the roll top (wedge). Rule met (0.429), but with a concave band of R ≈ 1.1 mm.
6. Height shed earlier. Worse (0.13 and still concave).
7. Gap form with a tangent arc (final): 0.412, convex.
8. Fit iterations: knots, grading, cross-derivative constraints, and the end-derivative attempt (rejected).

About 60% of the effort went into the TE corner.
