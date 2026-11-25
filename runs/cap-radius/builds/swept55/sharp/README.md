# swept55 / sharp: rounded tip caps

`wing.step` is the full symmetric wing: the parent from `parents/swept55/sharp/parent.step`, unchanged, with a rounded cap closing each tip. The parent is a level, untwisted 55° swept wing with a NACA 0009 section and a sharp trailing edge, 210 mm root chord, 95 mm chord at |Y| = 390 mm, and 900 mm projected span. It is one closed, valid solid with 6 faces:

- 4 parent skin faces, copied from `parent.step` with identical poles;
- 2 cap faces, the right one and its exact mirror in Y. Each is one B-spline face, degree 7 × 5 with 89 × 139 poles, running from the upper cut edge around the lip to the lower cut edge.

The parent's two flat end faces are dropped. STEP SHA-256: `5d511033d025d9084c9e6b5b5d0763e7ce853632b304e4ade9df16d91f6f2420` (2.36 MB, AP214, mm). A rebuild writes a new header timestamp, so the hash changes even when the geometry does not.

## How the cap is built

All of this is in `cap_model.py`.

**Guides.** Each wrap guide sits at constant chord fraction ξ = U(q) on the parent's own planform, continued outboard. In plan it is a straight line that follows the sweep. It starts on the cut at chord parameter q and ends on the locked outline at p = q. The outline formula puts lip point p at exactly that chord fraction.

**Wrap shape.** Along a guide, with x = sin θ the fraction of the reach (D = 240 q(1−q)):

`Z = t(q) · C(q,x) · ζ`

- ζ = 1 up to x0(q). This is the parent's own skin carried outboard.
- After x0 the roll is a superellipse, ζ = √(1 − ξᵇ) with ξ = (x − x0)/(1 − x0). It reaches the lip with a vertical tangent.
- At the cut the cap equals the extended parent plus an offset of order x³, so position, tangent plane and curvature match (G2).

**Knobs along the chord**, chosen by `design.py` and stored in `design_table.json`:

- **x0(q)**, where the roll starts:
  - 0 over the front 40 % of the chord, which gives a near-elliptic roll spread over the whole reach;
  - rising to about 0.93 aft, which holds the skin out and rounds only the last 7 % of the reach.
- **b(q)**, the superellipse exponent: 3, blended to 2.5 toward the trailing edge.
- **Taper easing k_u(q)**: from q = 0.6 aft, the parent's taper is eased out over the flat run, so the lip keeps the attachment thickness. Without it the cap thins by 6–8 % and the 0.4 mm floor fails at t ≈ 0.8 mm.

**Design targets.** x0 is solved per chord position so that the lip's principal radius is min(0.5 mm, 1.03 · t/2), where t is the attachment thickness at that chord position. It is capped where any larger value would put a tighter spot inboard of the lip, then slope-limited and smoothed along the chord.

**CAD construction.** All of this is in `build_cap.py`.

- **u direction.** u is the parent's own chord parameter on its own knots, with the degree raised from 5 to 7 (exact).
- **Fixed rows at each cut edge:**
  - row 0 is the parent's cut-edge poles, so position matches exactly;
  - row 1 is π·D(u)·∂P/∂v of the parent, so the tangent plane is exact;
  - row 2 sets the second derivative to π²D²·P_vv, so curvature matches.
- **Collapsed ends.** Columns 0 and 88 collapse to the leading-edge and trailing-edge points. D(u) vanishes there, so the cap closes as a two-pole patch.
- **Free poles** are a least-squares fit to the analytic cap. The fit is weighted by 1/(t + 0.01 mm) so the sub-micron trailing-edge rows keep their shape.
- **Wrap parameter.** It is shifted so that the start of the roll lies on a C2 knot line.

**Assembly.** The faces are sewn at 1e-6 mm. No free edges, BRepCheck valid.

## Measured on the reimported STEP (`readback.json`)

| Check | Result |
|---|---|
| Topology | 1 solid, 1 shell, 6 faces, 0 free edges, BRepCheck valid; volume 1 208 744.5 mm³ (parent 1 164 279.9; each cap +22 232 mm³) |
| Parent skin preserved | 0.0 mm against `parent.step`; 1.75e-7 mm against `parent.py` (the parent STEP's own fit error) |
| Locked dimensions | Bounding box Y = ±450.0000001; cut edges at |Y| = 390 exactly; lip to locked outline ≤ 8.8e-9 mm; lip normal horizontal (|n_z| ≤ 1e-7); no cap point outboard of the outline (≤ 6e-9 mm) |
| Wing-to-cap joins (4, u ∈ [1e-4, 1 − 1e-4]) | Gap ≤ 7.8e-12 mm; tangent-plane angle ≤ 1.5e-6°; curvature difference in a common frame ≤ 1.1e-5 /mm (≤ 6.4e-6 relative away from the ends) |
| Numerical zebra on the joins (4 views × 4 stripe angles) | Phase jump ≤ 2.9e-8 stripe; stripe-gradient jump ≤ 4.7e-4 relative; direction change at 1 937 light/dark crossings ≤ 5.3e-5° |
| Smallest principal radius where t ≥ 0.8 mm | 0.4148 mm on the lip at u = 0.880 (ξ = 0.960, t = 0.803 mm), (677.23, ±415.35, 0). Rule: ≥ 0.4 |
| Lip radius, thick part | 0.47–0.52 mm from mid-chord aft (target 0.5); 0.52–2.2 mm forward of it |
| Thin trailing-edge zone (t < 0.8 mm) | Lip radius 1.018–1.034 × t/2, and every wrap's tightest point is its lip. Rule: half-round of t/2 |
| Whole-cap map (618 × 501 grid per cap, local minima refined) | 0 points below the rule; smallest radius/rule ratio 1.018. The grid reaches 0.042 mm from the leading-edge collapsed point and 0.036 mm from the trailing-edge one |
| Lip-normal sections (34) | All meet the rule; worst 1.033 × the rule |
| Sweep-following sections (64, at and between u knots) | The surface's normal radius along every trace meets the rule (worst 1.032×). The raw planar-section radius does not (worst 0.005×), because these cuts cross the lip obliquely; see below |

## Deviations and deliberate features

- **Sweep-following sections read tighter than the surface.** A vertical cut along the sweep meets the lip at an angle, so the section curve's radius there is the surface radius × cos(angle). An example is 0.21 mm against a true 0.87 mm at u = 0.31. Printability is judged on the surface radius; both numbers are in `readback.json`.
- **Near the leading edge (u < 0.09) the tightest turn is mid-wrap, not at the lip.** It is 1.33–1.88 mm, against 1.9–2.2 mm at the lip. This is the parent's own leading-edge roll (true radius 1.478 mm at the cut) carried into the tip and opened up to the lip. It is 3.3× or more above the rule.
- **Bending reversals, all slight.**
  - Sweep cuts: concave curvature ≤ 1.4e-4 /mm (radius ≥ 7 m) in the flat run. It comes from the taper easing and the parent's ruled, slightly saddle-shaped skin. The saddle shape is inherited: 31 % of cap grid points have K < 0.
  - Lip-normal cuts at p = 0.44, 0.6 and 0.8: ≤ 7.4e-4 /mm (radius ≥ 1.35 m), where the aft cap is held thick.
  - Thin trailing-edge sections, at the start of the roll: ≤ 1.1e-3 × the local lip curvature. This is a fit ripple at the roll-start knot line. Its largest absolute value, 3.95 /mm, is on a section 0.6 µm thick.
- **Shape of the aft tip.** The aft half of the cap stays near full thickness until close to the lip. This is required to keep a ≥ 0.4 mm radius on the aft (trailing-edge-like) part of the outline. The zebra figures show it as a smooth hook in the stripes near the aft lip.

## Files and how to run (from the project root)

```
PY=.venv/bin/python; D=runs/cap-radius/builds/swept55/sharp
$PY $D/design.py        # x0(q) table -> design_table.json (about 45 s)
$PY $D/map_analytic.py  # optional: principal-radius check of the analytic cap
$PY $D/build_cap.py     # wing.step + build_report.json (about 20 s)
$PY $D/verify.py        # readback.json + map_right.npz (about 35 s)
$PY $D/render.py        # figures/*.png: zebra views from CAD normals, radius map
```

`cap_model.py` is the analytic cap, with exact derivatives through the parent's jet class. `skill-notes.md` is the account for the skill maintainer.
