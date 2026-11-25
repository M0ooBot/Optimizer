# Skill notes: close-rounded-wingtip on fin5 / sharp

These notes are for the skill's maintainer. I read only the skill, my parent folder (`parent.py`, `parent.step`, `spec.json`, `build.json`, `check.json`) and `src/`/`tools/` listings. I did not read the linked Wing B study, as the ground rules required, so anything the skill leaves to that study was unavailable to me.

## Steps followed, and how

**Build 1: attachment curves.** Followed. The cap's first pole row is the parent skin's top-edge pole row, copied from `parent.step`. The cap also reuses the skin's u knots and degree, so both joins are the parent's own edges. On a fin these are the +Y and -Y cut edges, as the skill says.

**Build 2: 3D perimeter.** The perimeter was locked, so I used it as given. It lies in Y = 0 and leaves the leading and trailing edges tangent to them in 3D. Its curvature is also zero there, because d²/dh² sqrt(1 - h³) = 0 at h = 0. There is no twist and the section is symmetric, so the elevation check in the skill was trivial.

**Build 3: spanwise wrap guides.** This is where I departed most. I did not draw individual guides. The locked outline is exactly the parent's elliptical chord shrunk by sqrt(1 - h³), so I built the cap as a family of elliptical horizontal sections. The chord follows the outline, and the half thickness follows one designed wrap profile Yp(h). That profile is the mid-chord guide. It was designed by its curvature distribution, not drawn. The guides at other chord positions follow from the section rule. Two consequences:
- My meridians (theta = const) run from the attachment to a collapsed crown pole. They do not "intersect all three profiles in order". The skill's guide, from skin to perimeter to skin at one chord location, corresponds to a cut through my surface. I measured those cuts (constant-X, sweep-following and perimeter-normal planes) but did not construct them.
- The skill's boundary surface "upper cut edge → perimeter → lower cut edge" is realised as two patches per side, mirrored. Each is bounded by its attachment edge, the two perimeter branches and the crown pole. The skill's fallback text ("If the CAD tool constructs upper and lower cap patches separately…") covers this.

**Build 4: boundary surface with tangency and curvature matching.** Met exactly by construction. G1 and G2 are imposed through the first three pole rows; the parent is ruled, so only S_vv ∥ S_Z is needed. The perimeter gets a common tangent plane and, by mirror symmetry, common curvature.

**Build 5: printable radius.** Met with margin everywhere. I treated it as a design requirement, as the skill asks. In practice my first design failed it and I only found out from the prototype's curvature map (iteration 1 below).

**Build 6: progressive turn.** Met only partly; see "Requirements" below. This geometry does not allow a strictly progressive turn with G2.

**Verify 1–5.** All five were done on the reimported STEP. Verify 5 used point-splatted zebra PNGs from exact CAD normals (`zebra_*.png`), because no CAD viewer was available.

## Where the skill was unclear, ambiguous or wrong

1. **"At the same chord position" is undefined for cap points.** A cap point does not sit at an attachment chord station, and the fin is swept and tapered. I used the sweep-following chord fraction (X − LE(Z))/c(Z), which equals the attachment chord fraction along the parent's rulings. The verdict here does not depend on that choice, because the smallest radius found (0.612 mm) is above 0.4 mm everywhere.
2. **The rule's rationale fails near the ends of a steep outline.** "A turn through the lip can be no wider than half the section's thickness" holds in a plane that contains the thickness direction and the surface normal. Near the leading and trailing edges, the outline climbs steeply (about 75° from horizontal at h = 0), so the turn across the lip lies mostly in the tangent plane. The principal radius there is set by the parent's nose (0.61–0.63 mm), not by t/2. The "< 0.8 mm: half-round of radius t/2" branch was therefore never binding, and the cap there is rounder than t/2. The skill could say that the nose radius, not t/2, governs where the outline runs along an edge.
3. **Verify 2 does not say which curvature of a section to judge.** Its "signed curvature along each trace" can mean the in-plane curvature of the cut curve or the surface's normal curvature along it. Near the leading and trailing edges, constant-X and sweep-following cuts graze the surface at the lip: the plane is 58–80° from the surface normal. Their in-plane radius is far below the surface's own:
   - 0.10–0.39 mm in-plane where t = 0.46–1.12 mm, against a surface radius of 0.62–0.82 mm.
   - Applied literally, the rule flags 19 of 80 traces (`readback.json` → `sections.traces_failing_rule_on_section_radius`).
   - The same cuts through the parent's own nose would read just as tight, and no smooth cap can avoid it.

   I judged printability on the principal curvature (Verify 3), which no trace fails. I also added cuts normal to the perimeter, which contain the surface normal at the lip and give the true lip radius. The skill should name those cuts, or say to judge sections by normal curvature.
4. **"On a sharp trailing edge the half-round shrinks to the parent's trailing-edge point."** On this parent, "sharp" means zero thickness at the TE, but the plan-view TE is still the ellipse's end, with a nose radius of 0.612 mm. The cap simply continues it.
5. **Tolerances are not given.** Verify 4 says to judge against "the design's tolerances", but the task gives none. I report raw numbers. The stripe period (0.1 in R·A) and the view set are my choices.
6. **Collapsed endpoints.** The skill anticipates them. My crown pole loses second derivatives to round-off within about 10 µm: the computed curvature rises from 0.1732 /mm to 1.8 /mm at 1e-6 in v. I stopped derivative-based sampling at v = 1 − 1e-4, which is 7 µm from the crown, and report the crown radii from the limit (5.76 mm and 522 mm).

## What I decided or improvised

- **Elliptic-section construction with the parent's own u parameter.** This gives exact G2 and the exact outline (within 5e-9 mm) without separate guide curves.
- **Wrap profile from a designed curvature distribution.** A smoothstep ramp over 12 % of the arc, then a smooth relaxation to the crown, solved so the curve lands horizontal at Z = 245. There is a double knot in v at the ramp end, where the profile is only C3.
- **The 12 % ramp, as a judgement call.**
  - Shorter ramps give a more even wrap but a more abrupt curvature onset just after the G2 join.
  - Scan at mid-chord (tightest / crown radius): 10 % → 4.39/5.63 mm, 12 % → 4.27/5.77 mm, 15 % → 4.10/6.02 mm, 20 % → 3.84/6.54 mm.
- **Chord-position definition** for the rule (item 1 above).
- **Zebra implementation.**
  - Phase: (reflected view ray · stripe axis)/0.1.
  - Gradient: from the exact shape operator.
  - Coverage: 40 views (both sides, azimuth −60..60°, elevation −30..60°) × 6 stripe axes, with visibility culling.
  - Crossings: located by sign change and linear interpolation.
- **Section families.** 27 constant-X, 27 sweep-following (the parent's constant-chord-fraction rulings) and 26 perimeter-normal planes.
- **Picture renderer.** Orthographic point-splatting with a z-buffer, used instead of a CAD viewer.

## Design iterations: 3, plus profile tuning

1. **Bezier wrap guides through a lip map h = sin(theta).** Each guide was a degree-6 Bezier starting along the skin ruling and turning in the plane of Y and the lip-to-root direction. Rejected at prototype stage: the smallest principal radius was 0.12–0.14 mm near the leading and trailing edges, and below 0.4 mm wherever t was up to about 2 mm. The lip map put each lip point well forward of its attachment point near the ends, by up to 0.4 mm at t = 2 mm. The guides therefore had to bend across the nose, in the direction of the surface normal, within about 1 mm.
2. **Elliptic sections with Yp = 5.1 sqrt(1 − h³).** The smallest radius was 0.612 mm, so the printable rule passed. But mid-chord had a shoulder: 3.37 mm at 40 % of the roll, opening to a flat 7.5 mm crown. That is exactly the pattern Build 6 warns about. Polynomial variants of the profile (1 − h³(A + Bh + Ch²)) gave wavy curvature and were dropped.
3. **Elliptic sections with the curvature-designed profile (final).** Mid-chord runs 4.27 → 5.77 mm. The smallest radius is unchanged at 0.612 mm, because it sits at the trailing-edge nose and is set by the parent.

## Requirements: met, partly met, not met

| requirement | status |
|---|---|
| one closed, valid solid; parent skin, attachment, span and outline unchanged | **met**: valid, 0 free edges; skin poles identical; outline within 4.8e-9 mm; height 245.000 mm |
| tangent and curvature match at both skin joins | **met**: ≤ 1.5e-6 deg; curvature difference 3.7e-8 relative |
| common tangent plane and smooth turn on the perimeter | **met**: ≤ 1.7e-6 deg; 4.7e-6 relative |
| printable radius: 0.4 mm where t ≥ 0.8 mm, t/2 below, judged on principal curvature | **met**: smallest 0.612 mm at (149.997, 0, 240.080), where t = 0.02 mm; smallest margin 1.68 (0.670 mm where t = 0.811 mm); 0 failing points |
| printable radius judged on the in-plane curvature of grazing constant-X or ruling cuts | **not met literally**: 0.10–0.39 mm on 19 traces near the edges; a measurement artefact (item 3) |
| end the cap the parent's way at the TE | **met**: the perimeter starts at the parent's TE point and continues its 0.61 mm nose |
| progressive turn with no tight bend before the lip | **partly met**: no reversals anywhere, but the tightest point sits below the lip at every station. Lip/tightest ratio is 1.35 at mid-chord (4.27 → 5.77 mm) and about 1.15–1.25 at 10–20 % chord. The locked 5 mm reach is below the 5.1 mm half thickness, so a G2 wrap cannot tighten all the way to the lip |
| zebra continuity across joins | **met**: phase jump ≤ 1.4e-9 periods and stripe direction ≤ 1.4e-5 deg on the skin joins; on the perimeter, 9.2e-8 periods and 5.6e-5 deg |
| review in context, fixed camera, no edge lines | **met, with my own renderer**: `zebra_*.png` |
