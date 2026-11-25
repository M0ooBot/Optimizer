# Skill notes: close-rounded-wingtip on swept55 / sharp

This is an account for the skill's maintainer. Every number comes from the reimported `wing.step` (`readback.json`) unless it is marked "analytic".

## What I followed, step by step

### Build

1. **Attachment curves.** Followed exactly. The cap's cut edges are the parent STEP's own edge poles: the u degree is raised 5 → 7, which is exact. The parent skin is untouched (0.0 mm against `parent.step`). The span and the station hold.
2. **Perimeter curve.** It was locked, so I used it as given. Its tangency to the leading and trailing edges and its Z = 0 elevation (the wing is level and untwisted) already held.
3. **Wrap guides.** Followed, but as a continuous family rather than a few drawn curves. Every u = const line of the cap is a guide: a straight plan line at constant chord fraction on the parent's planform continued outboard, ending at outline point p = q. Each guide's profile is the parent's skin carried out to x0(q), then a superellipse roll.
4. **Boundary surface.**
   - CadQuery/OCP has no boundary surface with curvature-matched guides, so I wrote the B-spline directly: one patch through the lip, fixed boundary rows for C0, G1 and G2, and least squares for the interior.
   - Measured result: tangent angle ≤ 1.5e-6° and curvature difference ≤ 1.1e-5 /mm. The perimeter is interior to the face, so there is no perimeter seam.
5. **Printable radius.** Treated as the design driver from the start. The design-stage tool is `design.py`, which solves the roll start per chord position against the rule.
6. **Progressive turn.** Checked per guide on the analytic surface, then on the STEP. From u ≈ 0.09 aft, every wrap is tightest exactly at its lip. The leading-edge zone is the exception, described below.

### Verify

1. **Reimport.** Done: topology, volume, parent preserved, outline, span, station, joins in a common tangent frame.
2. **Sections.** Done.
   - 64 sweep-following cuts at the u knots, midway between them, and thin ones near both edges. Each is traced from the parent skin 20 mm inboard of the cut, around the lip, and back.
   - 34 cuts normal to the lip.
   - Curvature on all of them uses CAD derivatives and Meusnier's theorem, so it is exact, not finite-differenced.
3. **Whole-cap map.** Done: 618 × 501 grid per cap, with every local minimum of radius/rule refined by Nelder–Mead. The grid reaches 0.042 mm from the leading-edge collapsed point and 0.036 mm from the trailing-edge one.
4. **Numerical zebra.** Done.
   - Phase and surface gradient come from CAD normals with Weingarten derivatives.
   - 4 views × 4 stripe angles, stripe axis mirrored for Y < 0.
   - 1 937 light/dark crossings located on the joins.
5. **Visual review.** Done in `render.py`. It uses depth-buffered point splats of the exact CAD surface (poles read back from the STEP; normals agree with OCC's D1 to 1e-13), with no mesh and no edge lines. The views are a centred whole-wing view with mirrored stripes and tip close-ups.

## Where the skill was unclear, ambiguous or wrong

1. **"Chord position" for a cap point is undefined.** The rule depends on "the thickness … at the attachment, at the same chord position", but a point on the cap is not at any attachment chord position. I used the guide it lies on (u = the parent's chord parameter), which is natural here because the outline maps p ↔ chord fraction U(p). The skill should say how to assign it, especially for guide layouts that do not pair cut points and lip points one to one.

2. **On a swept tip, the advice to spread the roll is backwards for the aft half.**
   - From mid-chord aft, the locked outline runs back to the trailing edge, so the aft lip is in effect a trailing edge of the tip. Its radius is set by the chordwise thickness across neighbouring guides, not by the shape along any one guide.
   - An elliptic, fully spread roll (my first family) gave lip radii of 0.08–0.1 mm at u = 0.6–0.8, where the section is 2–6 mm thick.
   - Meeting 0.4 mm there needs the opposite of spreading: hold the parent skin out to 80–93 % of the reach and round only the end. That keeps ≥ 0.8 mm of thickness up to the aft lip.
   - The skill's "usually means spreading the roll over more of the cap reach" pointed the wrong way for this shape.

3. **Section radius depends on the cut direction.**
   - Verify step 2 asks for sweep-following cuts and says to judge their smallest radius against the rule. A planar cut that crosses the lip obliquely reads radius × cos(angle).
   - On this wing those cuts read 0.005–0.5 × the rule at the lip everywhere, even at u = 0.31 (8.4 mm thick: 0.21 mm read against 0.87 mm true).
   - Taken literally, no rounded tip on a 55° planform can pass. I judged printability on the surface's normal and principal radius, reported the raw section radii as well, and added lip-normal sections, which do see the true lip radius.
   - The skill should say which radius is meant.

4. **"A half-round whose radius is half the local thickness" is ambiguous and not exactly achievable.**
   - "Local thickness" could mean the attachment (as step 5 says) or the cap section itself. The parent's taper alone makes the cap section 6–8 % thinner than the attachment at the trailing-edge lips.
   - The aft lip sits on a wedge, not a slab. A tangent circle can be up to about 1.18 × t/2 there, but only with a curvature jump.
   - A G2-smooth roll whose curvature peaks at the lip reaches about 1.0–1.08 × t/2 with my families.
   - I targeted 1.03 × t/2 and got 1.018–1.034 × t/2.

5. **The switch at t = 0.8 mm has zero margin.** Just above 0.8 mm the rule is ≥ 0.4 mm, which equals t/2, so a design aimed at t/2 sits exactly on the limit. I needed the 3 % margin above, and the taper easing, to read ≥ 0.4 on the STEP: 0.4148 at t = 0.803.

6. **The thin zone at the leading edge does not behave like the trailing edge.** The rule calls the first 0.1 % of chord (t < 0.8 mm) a half-round zone. There, the cap continues the parent's leading-edge roll, whose true radius is 1.48 mm. So the radius is 1.33–2.2 mm, far above t/2. I counted that as meeting the rule ("as round as the section allows"), but the wording does not fit.

7. **The leading-edge zone breaks the "tightest at the lip" pattern.** For u < 0.09, the tightest turn is mid-wrap (1.33 mm), not at the lip (2.2 mm). It is inherited from the parent's leading-edge roll, and it is not the "tight bend just before the perimeter" the skill warns about. I left it and documented it.

8. **Reversal counting needs a tolerance.**
   - The parent is a ruled, tapered skin (saddle-shaped, K ≤ 0), so any cut slightly off its rulings shows reversals at the 1e-5 /mm level.
   - The taper easing I added is a deliberate, very slight concave bend: radius ≥ 7 m on sweep cuts, ≥ 1.35 m on lip-normal cuts.
   - Fit ripples at the roll-start knot on the thin trailing-edge sections are ≤ 1.1e-3 of the local lip curvature.
   - I counted reversals only above 1e-5 /mm and reported the largest concave curvature. The skill gives no zebra or curvature tolerances ("judge against the design's tolerances"), so I report numbers only.

9. **Collapsed endpoints are inevitable here, and the skill does not say so.** The perimeter starts at the cut's leading-edge point and ends at its trailing-edge point, so the boundary surface must collapse to a point at both ends. Keeping G1 there needs the cross-boundary derivative to vanish at the poles; I used D(u) = 240u(1−u) as a factor of the derivative row. This is worth a sentence in the skill.

10. **"Guide locations and between guides" assumes discrete guides.** With a continuous guide family, I used the B-spline u knots and the midpoints between them.

11. **I could not read the Wing B study the skill links to.** The ground rules forbid it, so whatever guidance it holds was unavailable.

## What I decided or improvised that the skill does not cover

- A one-knob profile family, x0(q), solved per chord position against the rule, with the shoulder-free limit computed for every q.
- The superellipse exponent drops to 2.5 toward the trailing edge, for a rounder lip at the same thickness.
- Taper easing aft of q = 0.6.
- A slope limit and smoothing on x0(q). Without them, the jump from a spread roll to a held-out roll at q ≈ 0.40 made a concave crease (analytic radius 12.5 mm).
- A wrap reparameterisation that puts the roll start on a knot line.
- A thickness-weighted fit.
- Lip-normal sections as the primary section check.

## Iterations

There were 6 iterations: 4 design families on the analytic surface, then 3 CAD builds, the first of which tested the fourth family.

1. **Analytic: z = √(1 − xᵃ).** Pushing a far enough for the aft lip made a shoulder 10–17 % tighter than the lip. The best shoulder-free result was about 0.92 × t/2. Dropped.
2. **Analytic: flat run plus cubic superellipse (x0).** Shoulder-free, but only 0.91–0.95 × t/2 near t = 0.8, so the 0.4 mm floor failed.
3. **Analytic: exponent 2.25 plus taper easing.** Met 1.03 × t/2.
4. **CAD build 1.**
   - The STEP map failed within about 1 mm of the trailing-edge point (ratio down to 0.02), because the fit lost the sub-micron shape there.
   - Concave fit ripples up to 4.6e-3 /mm appeared at the roll start.
5. **CAD build 2.**
   - Fixes: thickness-weighted fit, roll start aligned to a C2 knot, exponent eased to 2.5.
   - The map passed everywhere, but the analytic design still had a concave crease (12.5 mm) at q ≈ 0.40.
6. **CAD build 3 (final).** Slope-limited and smoothed x0(q).

## Requirements: met, partly met, not met

| Requirement | Status |
|---|---|
| One closed, valid solid; parent, station, span and outline unchanged | Met. 0 free edges; parent 0.0 mm; outline 8.8e-9 mm; span 900.0000002 mm |
| Tangent match to both skins at the attachment | Met. ≤ 1.5e-6° |
| Curvature match and eased into the roll | Met at the join (≤ 1.1e-5 /mm). Inside the cap, curvature rises from 0 at the roll start |
| Radius ≥ 0.4 mm where t ≥ 0.8 mm | Met on the whole-cap map and on lip-normal and surface-radius section checks. Smallest 0.4148 mm at ξ = 0.960 (t = 0.803 mm). Not met by raw sweep-cut section radius (see note 3) |
| Half-round of t/2 where t < 0.8 mm (trailing edge) | Met as lip radius 1.018–1.034 × t/2 with the lip tightest (by attachment thickness). Not an exact half-round; see note 4. Near the leading edge, met with radius ≫ t/2 (note 6) |
| Progressive turn, no shoulder before the lip | Met for u ≥ 0.09. Leading-edge zone: deliberate mid-wrap minimum (note 7) |
| No bending reversals | Partly met. Only very slight concave bends, all recorded with location and size (note 8) |
| Numerical zebra continuity across both joins | Met. Phase ≤ 2.9e-8 stripe; gradient ≤ 4.7e-4 relative; direction change ≤ 5.3e-5° |
| Visual review, centred camera, mirrored stripes, no edges | Done (`figures/`). No marks at the joins. A smooth stripe hook near the aft lip, which is the held-thick aft roll |
