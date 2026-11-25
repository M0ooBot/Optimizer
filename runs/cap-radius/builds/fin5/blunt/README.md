# fin5 / blunt: rounded tip cap

`wing.step` contains the whole fin (parent plus cap) as one closed solid: 1 solid, 1 shell, 5 faces. Units are mm, X runs aft, Y is thickness and Z is span. The parent is unchanged. The cap closes the fin between the attachment at Z = 240 and the crown at Z = 245, and its lip lies on the locked side-view outline.

STEP SHA-256: `cc8b5b6248eeae837dd693c4eb0dbedd32e916639017fdb223768da48b563318` (816,285 bytes). The STEP writer stamps the time into the file, so each rebuild gives a new hash; `readback.json` always records the current one.

## How it was built

1. **Analytic cap** (`cap_analytic.py`). The fin is symmetric, so each skin can be written as Y² = G(X, Z). For the parent, G_par = 0.0144·a²·(1+ξ)(ξ_f−ξ), with ξ = (X−m)/a the normalised chord position (m = mid-chord, a = half chord, ξ_f set by the keep-chord blunt cut). The cap uses
   G_cap = G_par · ψ, where ψ = 1 − (1+ε)·x³/(x+ε), x = h/(1−ξ²)^(1/3), h = (Z−240)/5 and ε = 0.25.
   - ψ = 1 + O(h³) at the attachment, so the cap is curvature-continuous (G2) with both skins.
   - ψ = 0 exactly on the locked outline ξ² + h³ = 1, so the lip *is* the outline.
   - Y² = G is smooth across Y = 0, so the two sides share a tangent plane and their curvature matches at the lip.
   - ε shapes the wrap. With ε = ∞ (a pure cubic) the crown section tightens to r = 3.4 mm and then opens out to 7.5 mm at the lip. With ε = 0.25 the crown section runs from 4.4 mm to 5.5 mm, and the LE lip radius stays at 0.48 mm or more.
2. **Perimeter and guides.** The perimeter is the locked outline, which leaves tangent to the LE and TE lines. Each wrap guide is the cut through the analytic cap by the plane that contains the Y axis and the side-view line from the attachment point A(s) to its perimeter point P(s). A and P are paired by equal arc-length fractions along the cut edge and along the outline.
3. **One B-spline face**, degree 5×5 with 89×93 poles (`build_cap.py`).
   - u is the parent's own chord parameter s: the parent's knots, refined at both ends.
   - t is the normalised arc length around each guide: −Y cut edge, then the lip, then the +Y cut edge.
   - The boundary poles are fixed: the parent's top pole rows (so the edges are shared exactly), the LE point (a collapsed u = 0 edge) and the straight TE segment (u = 1).
   - The interior poles come from one sparse least-squares solve. It fits 1.37 M exact analytic points and, with heavy weights, the analytic first and second cross-boundary derivatives at both cut edges. Near the TE the samples use local coordinates (distance forward of the TE line) because the thin closure there is below double precision in global X.
   - The fit error is 6.9e-5 mm or less for s ≤ 0.99 and 3.6e-4 mm at most, in the TE-corner zone.
4. The cap face replaces the parent's flat attachment face. The faces are sewn at 1e-6 mm, made into a solid and exported as AP214.

## Measured on the reimported STEP (`verify_cap.py` → `readback.json`)

| Check | Result |
|---|---|
| Solid | Valid (BRepCheck), 0 free edges; volume 352,563.395 mm³ |
| Cap volume | Whole minus parent is 2,402.964 mm³; the analytic design gives 2,402.962 mm³ |
| Parent | Retained skins and TE face have the same poles and knots as `parent.step` (difference 0); `parent.py` matches the skins to 8.6e-8 mm |
| Locked geometry | Lip within 5.4e-6 mm of the outline; crown Z = 245.000000 at X = 108.0728; cap never more than 3.8e-5 mm outside the outline; attachment at Z = 240; bounding box X 0–160, Y ±9.601, Z 0–245 |
| Skin joins (±Y) | Gap ≤ 2.6e-11 mm. Tangent angle ≤ 6.6e-5° beyond 0.1 mm of the TE corners and LE point, ≤ 7.0e-4° beyond 0.01 mm. Cross-edge curvature difference ≤ 1.7e-5 /mm beyond 1 mm and ≤ 2.8e-4 /mm beyond 0.2 mm |
| TE face join | Shared segment; tangent angle 2.3e-4° for |Y| ≤ 0.05 and 0.023° for |Y| ≤ 0.1, rising to 3.9° at 0.01 mm from a corner (see the limitation below) |
| Tightest turn (whole-cap map, refined) | Outside the TE-corner zones: r = 0.477 mm at (65.131, 0, 240.433), the lip near the LE. The attachment there is 0.80 mm thick, so 0.4 mm is required (ratio 1.19). Near the TE: 0.606 mm at (149.773, 0, 241.028) |
| Wrap sections (41 planes: guides, constant X, sweep-following) | Every section beyond 0.4 mm of the TE corners meets the rule, with no bending reversal. At the crown the section is tightest at 4.44 mm, about 6 mm before the lip, and is 5.52 mm at the lip (ratio 0.80). Off the crown the tightest point is the lip |
| Zebra, 20 views × 6 stripe axes, 10° period, both joins | Beyond 0.2 mm of the TE corners and LE: phase mismatch ≤ 1.3e-5 period, stripe-direction change at light/dark crossings ≤ 0.033°, band-spacing change ≤ 0.03% (2,681 crossings). Beyond 0.1 mm: ≤ 0.14° |

## Limitations

- **TE corner vertices.** The printable-radius rule fails within 0.20 mm of the two TE corner points (150, ±0.2, 240):
  - 0.15–0.20 mm from a corner: smallest radius 0.19 mm against a required 0.23 mm.
  - 0.10–0.15 mm: 0.030 mm.
  - At the corner the radius goes to about 0.

  This is where the parent's 18° skin/TE-face crease ends. A cap that is tangent to the skins and contains the TE segment cannot have a single tangent plane at those points, so its curvature has no upper bound there. The analytic design behaves the same way (radius below the rule within 0.22 mm). In the same zones the B-spline has small concave ripples within 10 µm above the attachment: at most 0.008 /mm beyond 0.3 mm of a corner, larger closer in.
- **LE point.** The u = 0 edge collapses to the LE point, where derivatives are undefined. The curvature grid reaches 1.4e-5 mm from it. Radii below the rule occur only within 0.26 µm of the point.
- **TE end of the lip.** Near the TE, the lip continues the flat TE face and rolls over it. It is blunter than a 0.2 mm half-round: the section tip radius is about 680 mm at Z = 240.05, 11 mm at 240.2 and 1.1 mm at 240.5. It is therefore not literally a "half-round of radius 0.2" (see `skill-notes.md`).

## Files and how to run (from the project root)

```
.venv/bin/python runs/cap-radius/builds/fin5/blunt/build_cap.py    # ~4 min: wing.step, build.json
.venv/bin/python runs/cap-radius/builds/fin5/blunt/verify_cap.py   # ~3 min: readback.json, fig_*.png
```

- `cap_analytic.py`: the cap definition, perimeter pairing and wrap sections.
- `build_cap.py`: the fit, sewing and export.
- `surfeval.py`: B-spline evaluation of the reimported faces, cross-checked against Open CASCADE D2 to a relative 3.7e-9.
- `verify_cap.py`: all of the checks above.
- Figures:
  - `fig_radius_ratio_map.png`: radius divided by required radius over the cap, with the TE-corner close-up.
  - `fig_zebra_*.png`: matplotlib zebra renders with a fixed camera and no edge lines.
