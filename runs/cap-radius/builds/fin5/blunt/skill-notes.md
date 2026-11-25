# Skill notes: close-rounded-wingtip on fin5 / blunt

For the skill maintainer: how each step went on a swept, tapered fin with 12% elliptical sections and a blunt 0.4 mm TE, plus a 5 mm cap on a locked outline. There is no SolidWorks here; everything was built with CadQuery/OCP and NumPy.

## Steps followed and how

**Build 1: attachment curves.** Used as stated. The cap's two long boundary rows are the parent skins' top pole rows. I inserted knots into them but did not move them, so the edges are shared exactly (gap 2.6e-11 mm). The skill does not say what to do with the third attachment curve that a blunt parent has: the TE segment. I made it the cap's u = 1 edge.

**Build 2: 3D perimeter.** The perimeter is fixed by the locked outline, which is already tangent to the LE and TE lines. The skill's plan view and elevation view collapse into the Y = 0 plane for a symmetric fin.

**Build 3–4: guides and boundary surface.** OCP has no boundary-surface feature with profile and guide tangency settings, so I could not follow the recipe literally. Instead:
- I wrote the cap as an analytic implicit surface Y² = G_par·ψ. It is G2 to the skins by construction, and its lip is exactly the outline.
- The guides are planar cuts of that surface. Each cut contains Y and the side-view line from the attachment point A(s) to its paired perimeter point P(s). The pairing uses equal arc-length fractions.
- I fitted one B-spline face to exact points, plus the analytic first and second cross-boundary derivatives.

"Match tangency/curvature to both skins" became explicit least-squares equations rather than CAD options. That worked well: the joins are tangent to 6.6e-5° and the curvature difference is 1.7e-5 /mm.

The skill offers two patches as an option. I used one face that runs across the lip, so the lip has no seam. The skill's "check the shared perimeter separately" therefore did not apply.

**Build 5: printable radius.** I used the spec's attachment-thickness table and mapped each cap point to chord fraction (X − 65)/85. The required radius was taken as 0.4 where t ≥ 0.8, and t/2 where t < 0.8.

**Build 6: progressive turn.** Followed: a design study over the profile parameter ε, described under the iterations below.

**Verify 1–5.** All five checks were done on the reimported STEP:
- solid, parent and locked-geometry checks;
- gap, tangent angle and curvature in a common frame at matching parameters along both skin joins and the TE segment;
- 41 plane traces plus 11 horizontal sections, with signed curvature, reversals, the smallest radius and the lip radius;
- a dense principal-curvature map with refined local minima;
- a numerical zebra check: 20 views × 6 stripe axes, phase, gradient, and light/dark crossings with stripe direction and spacing;
- matplotlib zebra renders with a fixed camera.

## Where the skill was unclear, ambiguous or wrong

1. **The blunt-TE ending is ambiguous.** "The half-round closes the end of the flat trailing-edge face, so its radius there is half the trailing-edge thickness" has two readings:
   - (a) In the TE-face plane, the face ends in a semicircle of radius 0.2. The skins then meet that face extension at the parent's 18° crease. That puts a crease, radius 0, inside the cap.
   - (b) In the lip's cross-section near the TE end, the lip is a 0.2 mm half-round. That cannot join the flat 0.4 mm face at Z = 240 without a singular transition, because the outline's TE branch follows the TE line to within a·h³/2.

   I used neither literally. My cap continues the flat face tangentially (normal angle 2.3e-4° in the middle of the segment) and rolls over it with radii of 0.2 mm or more, except at the corners. The lip near the TE end is therefore blunter than a 0.2 mm half-round, not tighter. The skill should say which reading it means.
2. **The TE corner vertices are unavoidable singular points, and the skill does not mention them.** On a blunt TE, the skin and the TE face meet at a crease (18° here). Any cap that is tangent to the skins and contains the TE segment needs two tangent planes at each corner, so its radius tends to 0 there. The radius rule therefore cannot hold right up to those points.

   The skill only allows excluding "collapsed endpoints where ordinary surface derivatives are undefined". These corners are not collapsed: on a B-spline patch they are regular corners with a horizontal tangent plane, 87.7° off the skin. The skill should name them, say how close to them the rule must hold, and ideally say how to shrink the zone.
3. **A progressive turn can conflict with curvature matching.** If curvature rises monotonically from 0 at a G2 attachment to its maximum at the lip, the guide always climbs more height than it runs inward. Here the crown climbs 5 mm while the half-thickness is 5.10 mm, so a monotonic profile is impossible: some opening-out toward the lip is forced. The cleanest result was a quick rise followed by an almost constant radius: the crown section is 4.44 mm at its tightest and 5.52 mm at the lip. The skill hints at this ("a flatter crown … can call for a different progression"). Stating the inequality would save a design loop.
4. **Section radius is not principal radius.** Oblique cuts, such as horizontal sections near the crown or constant-X planes near the LE and TE, report section radii far below the surface's own radii (Meusnier's theorem). For example, the Z = 244.8 section has a 0.15 mm section radius where the surface radius is 5.1 mm. Step 2 asks for the "smallest sampled radius" of sections and judges it against the rule. It should say to judge principal (or normal-section) radii, or to cut along the outline normal.
5. **No zebra formula and no tolerances.** I defined phase = asin(r·a)/10°. "Judge against the design's tolerances" had no tolerances to judge against, so I report the numbers without a pass/fail.
6. **Where the half-round rule applies.** It is unclear whether the "thinner than 0.8 mm" rule is a minimum radius (≥ t/2) or an exact shape (= t/2). I treated it as a minimum. Near the LE, t < 0.8 only within 0.13 mm of the leading edge, and the parent's own nose radius there is 0.61 mm.

## Decisions and improvisations not covered by the skill

- The analytic implicit cap Y² = G_par·ψ. The outline's own sqrt(1 − h³) form suggested ψ, and it gives G2 and the exact lip for free.
- Pairing of attachment points with perimeter points by equal arc-length fraction, and planar guides that contain Y.
- Local trailing-edge coordinates. Near the TE segment, the cap's closure is narrower than double precision in global X (1e-28 mm at s = 1 − 1e-10).
- Expressing the derivative equations in pole units. Without that scaling the joint least-squares solve was singular.
- Choosing ε = 0.25. This trades crown flatness against the LE lip radius: smaller ε would push the LE lip below 0.48 mm.

## Iterations

- **Design: 2.**
  - Pure cubic ψ: crown section 3.37 mm at its tightest, 7.45 mm at the lip. That is the "tight bend then flatten" pattern the skill warns about.
  - ε-blended ψ, studied at ε = ∞, 1, 0.5, 0.3, 0.2 and 0.1; chose 0.25.
- **Representation: 5.**
  1. A separable fit with samples stopped at s = 0.9999. Near the TE it produced radii of 5e-8 mm, the cap dipped 3 µm below Z = 240, and it was 1 µm outside the outline.
  2. Local TE coordinates with samples taken to s = 1 − 1e-10. Fit error 3.6e-4 mm, but curvature at the joins was poor: 80 /mm beyond 0.2 mm, caused by very small t-spans.
  3. Imposing the derivative rows one after another. The fit blew up near the TE because the unconstrained end poles went to the minimum-norm solution.
  4. A joint sparse least-squares solve with weighted derivative equations. It was singular until the equations were normalised.
  5. Relaxed t-knot grading at the edges to reduce ripple.

## Requirements: met, partly met, not met

| Requirement | Status |
|---|---|
| Closed valid solid; parent, attachment, span and outline unchanged | Met. Lip within 5.4e-6 mm of the outline; parent poles identical |
| Tangent to both skins | Met beyond 0.01 mm of the TE corners and the LE point (≤ 7.0e-4°; ≤ 6.6e-5° beyond 0.1 mm). At the corner points themselves: 87.7° |
| Curvature matched and eased into the roll | Met. ≤ 1.7e-5 /mm beyond 1 mm; ≤ 2.8e-4 /mm beyond 0.2 mm |
| Printable radius ≥ 0.4 mm (t ≥ 0.8) or ≥ t/2 | Met everywhere except within 0.20 mm of the two TE corner vertices. There: 0.19 mm against 0.23 mm at 0.15–0.20 mm, 0.03 mm at 0.10–0.15 mm, about 0 at the vertex. Elsewhere the minimum is 0.477 mm at (65.131, 0, 240.433), where t = 0.80 and 0.4 mm is required. At the LE: below the rule only within 0.26 µm of the collapsed point |
| End the cap like the parent's TE ("half-round, r = 0.2") | Partly met. The face continues tangentially and the lip is blunter than 0.2 mm, not a 0.2 mm half-round (ambiguity 1) |
| Progressive turn, no shoulder | Partly met. Off the crown the tightest point is the lip. At the crown the section opens from 4.44 mm to 5.52 mm (ratio 0.80); point 3 shows some opening-out is forced |
| No unintended bending reversals | Met beyond 0.4 mm of the TE corners. Closer in, concave ripples within 10 µm above the attachment: ≤ 0.008 /mm beyond 0.3 mm, up to 4.5 /mm inside 0.15 mm |
| Zebra continuity | Beyond 0.2 mm of the corners and LE: phase ≤ 1.3e-5 period, direction ≤ 0.033°, spacing ≤ 0.03%. Beyond 0.1 mm: direction ≤ 0.14°. Beyond 0.01 mm: direction ≤ 1.75° |
| Visual review with a CAD viewer | Partly met. Only matplotlib renders were possible; no CAD viewer |
