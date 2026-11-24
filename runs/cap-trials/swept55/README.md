# Rounded tips on a 55° swept wing

The [whole-wing STEP](swept55-wing.step) contains one closed eight-face solid. It measures 900 mm across the tips, with a 55° leading-edge sweep measured aft from the spanwise axis, cap cuts at ±390 mm, and 210 mm and 95 mm chords at the root and cuts. The chord lines are level and untwisted. View the [whole-wing plan and zebra image](whole-wing-review.png), [enlarged zebra views](zebra-detail.png), and [oblique sections](lip-sections.png).

## Geometry

The new parent wing uses a symmetric NACA 0009 thickness distribution and straight linear taper. The final thickness coefficient changes from −0.1015 to −0.1036 to close the trailing edge at zero thickness. At the root, this removes the nominal 0.397 mm full trailing-edge gap of the original formula. The two parent halves mirror at the root without a center fillet. The straight tapered parent and its sections were generated for this study. [NASA four-digit airfoil formulation](https://ntrs.nasa.gov/api/citations/19760003945/downloads/19760003945.pdf)

Each wingtip has a single planar lip at zero elevation. In plan, its spanwise extent outboard of the cut is 240q(1−q) mm, where q runs from zero at the leading edge to one at the trailing edge. Its chordwise coordinate within the extrapolated tapered band uses 3q²−2q³. This gives the lip the outgoing leading-edge tangent and incoming trailing-edge tangent, and places its greatest outboard extent exactly 60 mm beyond the cut.

The cap upper and lower faces share that lip. The fixed-q curves act as spanwise wrap guides from the exact parent cut through the lip to the opposite skin. At the cut, the chordwise shift and thickness roll begin at third order, allowing the cap to match the parent curvature. At the lip, the cap's first derivative is vertical, and the projected outboard distance turns quadratically. The mirrored lower and upper patches therefore meet with a common tangent plane and curvature. The extra easing term reduces the initial forward section shoulder while preserving the lip radius.

CadQuery's high-level loft does not expose the required second-derivative boundary conditions. The build converts the polynomial maps directly to OCCT Bézier faces, sews them as one shell, and checks the solid before exporting. The evidence uses an independent STEP read that checks the resulting B-spline faces. Zebra shading uses OCCT derivatives and a custom reflection model.

## STEP measurements

| Check | Measured result |
| --- | --- |
| Strict OCCT solid validity and closure | Pass |
| Preserved generated parent surface | Maximum sampled deviation 5.46 × 10⁻¹² mm |
| Left/right mirror symmetry | Maximum sampled deviation 6.06 × 10⁻¹¹ mm |
| Lip end tangents relative to parent leading and trailing edges | Maximum 1.55 × 10⁻¹¹° on the right-side lip |
| Skin-to-cap gap | Maximum 2.74 × 10⁻¹² mm |
| Skin-to-cap tangent-plane angle | Maximum 1.46 × 10⁻⁹° |
| Skin-to-cap curvature tensor difference | Maximum 1.52 × 10⁻¹⁰ mm⁻¹ |
| Sectional lip radius | 0.172–1.17 mm across nine fixed oblique planes from 2% to 97% chord at the cut |

The [measurements](readback.json) record the STEP hash and all six joins. Each join uses 701 locations from 0.1% to 99.9% in q and 144 combinations of view and stripe angle. The largest stripe phase mismatch was 1.94 × 10⁻¹⁰ cycles, and the largest gradient mismatch was 1.17 × 10⁻⁹ cycles/mm. Both occurred at 0.1% on the left upper cut. The largest direction change at an actual light/dark transition was 1.03 × 10⁻⁷° at 0.176% on that cut. Both upper-to-lower lip joins stayed within 3.15 × 10⁻¹⁶ cycles/mm of gradient continuity and 3.63 × 10⁻¹⁴° of stripe direction continuity. The reimport passes the numerical study limits listed in the JSON.

The cap patches collapse to vertices at the leading and trailing ends. Those vertices have no regular two-direction derivative and are outside the join sweep. The nine illustrated sections run from 23.4 mm inboard of the cuts to their respective lips. Their oblique planes follow the leading-edge sweep and the plots use an expanded vertical scale. Nearer the leading edge, beyond the first illustrated cut, changing local chord position can still produce a small thickness shoulder. The sampled radii describe the geometric roll and do not establish a manufacturing minimum. The root seam and sharp trailing edge are outside the cap continuity checks.

The [section curvature audit](section-audit.json) follows the same nine planes from the parent through the cap. It resolves curvature sign changes at 2%, 90%, and 97% chord. The smallest sampled radius is 0.0421 mm at 97% chord and (y = 398.439) mm. The radius at the shared lip in that section is 0.1719 mm.

## Reproduce

Run these commands from the repository root using the installed Python environment. The study used CadQuery 2.8.0, OCP 7.9.3.1, NumPy 2.5.4, SciPy 1.18.1, and Matplotlib 3.11.2.

```powershell
python -B runs/cap-trials/swept55/build.py
python -B runs/cap-trials/swept55/check.py
python -B runs/cap-trials/swept55/section_audit.py
python -B runs/cap-trials/swept55/zebra.py
```

The build writes the STEP. The checker reads the STEP and writes its fingerprint, measurements, and section plot. The zebra renderer reads the same STEP and writes the whole-wing and close-up PNGs. Running the build again can change the STEP header timestamp, so rerun the checker afterward.
