# Swept vertical fin with 5 mm rounded cap

`wing.step` contains the whole fin as one closed, valid solid. Axes are X aft, Y transverse, and Z up, in millimetres. The root chord is 160 mm. At Z = 240 mm the leading edge is at X = 65 mm and the chord is 85 mm. The parent skin uses exact 12% elliptical sections with linear sweep and taper.

The cap follows the locked XZ outline `m ± a sqrt(1 − h³)` at Y = 0, where `h = (Z − 240) / 5`, `m = 80 + (27.5/240) Z`, and `a = 80 − (37.5/240) Z`. Two mirrored tensor-product patches use exact rational half-ellipses across the chord and a quintic spline along the span. The cap span parameter `t` gives `h = 2t − t²`, which keeps the wrap regular along a physical center section through the tip. The cap matches value and the first two derivatives of the parent at the join.

The cap multiplies transverse thickness by `1 − 0.5 h³`. This preserves both attachment derivatives and the locked XZ silhouette. The uncorrected 12% cap had a center-section shoulder radius of approximately 3.4 mm and a larger 7.47 mm lip radius. With the correction, the center's interior curvature maximum has a radius of 3.007 mm at Z = 241.870 mm. The center lip turns more tightly at 1.867 mm. Smaller interior maxima remain around Z = 241.87–241.96 mm over the sampled central chord range. The curvature retains its sign through those local variations.

## STEP readback

- Crown at (108.072917, 0, 245) mm. Maximum sampled locked-outline error was 7.11e-09 mm.
- At both parent-to-cap joins, the maximum sampled gap was 0 mm, normal angle 3.06e-09°, and common-frame shape operator mismatch 5.95e-10 mm⁻¹. Along the two branches of the cap symmetry seam, the largest sampled normal angle was 2.53e-10°. The crown is a collapsed UV boundary and is excluded from the seam's ordinary surface derivative checks.
- Numerical reflection stripes sampled 361 points per join across four camera directions and six stripe angles. The largest phase difference was 1.23e-09 rad and the largest phase-gradient difference was 1.35e-08 rad/mm. The largest direction difference over 928 located light/dark crossings was 5.77e-06°.
- Eleven solid-plane sections follow the fin sweep and taper at constant uncontracted chord fraction. Curvatures use implicit derivatives of the reimported STEP faces, sampled 2,001 times per side per cut. Every cut had its smallest sampled radius at the lip and no signed-curvature reversals above 2 × 10⁻⁵ mm⁻¹. Lip radii span 0.5370–2.6839 mm. The center lip radius is 1.8675 mm.

`readback.json` records the STEP SHA-256, sample locations, maxima and section results. `zebra.png` shows a fixed lateral orthographic view and cap detail using STEP surface normals, with edge overlays hidden.

## Reproduce

From the project root, run `python fin/build.py`, then `python fin/verify.py`. These use Python, CadQuery 2.8.0 and its Open CASCADE bindings, NumPy and SciPy. Run `python fin/render_zebra.py` with Matplotlib to regenerate the zebra views. Re-export may change the STEP header timestamp and the matching hash.
