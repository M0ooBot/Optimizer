# Rounded vertical fin

[Whole-fin STEP](vertical-fin.step) · [Outline and zebra views](fin-overview.png) · [Wrap sections](fin-sections.png) · [STEP readback](readback.json)

| Dimension | Reimported STEP |
| --- | ---: |
| Overall height | 280 mm |
| Root chord | 160 mm |
| Cap attachment height | 240 mm |
| Chord at attachment | 85 mm |
| Leading-edge sweep aft at attachment | 65 mm |

This is a new geometry study. The parent has linear sweep and taper, symmetric elliptical sections with 12% thickness-to-chord, and a flat root. The center plane is at Y = 0, positive X is aft, and positive Z is up. The whole fin is a solid envelope with no hinge, mounting geometry, or aerodynamic section prescription.

## Cap construction

The cap uses the left and right parent edges at Z = 240 mm. Its perimeter follows the leading edge up to the crown and back down to the trailing edge. The STEP represents that outline as two branches that meet at the crown. Both cap patches share the same perimeter. Sections inclined with the local sweep trace the wrap from one cut edge, through the perimeter, and back to the other cut edge.

The CAD kernel constructs those boundaries and wraps through a common section law in two tensor-product surfaces. In the cap, the squared section scale is (1-h^3), where (h=(Z-240)/40). This preserves position and the first two height derivatives at the attachment. It yields a rounded crown at Z = 280 mm. Rational B-splines represent the elliptical sections exactly and a constrained spline represents the cap contraction in height. Both patches and the measured sections follow the resulting surface directly. The build uses CadQuery and Open CASCADE through Python. The editable construction is the source script, and the exchange artifact is STEP.

## STEP verification

The reimport contains one closed, strictly valid Open CASCADE solid with five faces. The two parent skins agree with the designed swept ellipse to a maximum dimensionless profile residual of 2.2 × 10⁻¹³. Each perimeter branch follows the associated parent edge tangent at the attachment within 8 × 10⁻¹⁰°. Both skins register within numerical precision over the sampled attachment.

| Maximum sampled mismatch | Parent-to-cap joins | Cap perimeter branches |
| --- | ---: | ---: |
| Position gap | 0 mm | 0 mm |
| Tangent-plane angle | 7.7 × 10⁻¹⁰° | 3.1 × 10⁻¹⁰° |
| Normal-derivative operator | 2.7 × 10⁻¹¹ mm⁻¹ | 8.7 × 10⁻¹¹ mm⁻¹ |
| Zebra phase | 8.0 × 10⁻¹¹ cycles | 3.2 × 10⁻¹¹ cycles |
| Zebra phase gradient | 1.6 × 10⁻¹⁰ cycles/mm | 5.2 × 10⁻¹⁰ cycles/mm |
| Zebra stripe direction at crossings | 1.2 × 10⁻⁷° | 3.3 × 10⁻¹⁰° |

The zebra sweep uses 501 samples on each attachment and perimeter branch and 168 camera and stripe-angle combinations. It evaluates the phase with the same shader formula as the rendered views and solves light/dark crossings on the CAD geometry. The largest attachment stripe turn occurred near X = 139.02 mm, Y = −3.42 mm, Z = 240 mm. The largest perimeter gradient difference was on the trailing branch near X = 112.15 mm, Y = 0, Z = 279.99996 mm. The exact crown is a collapsed surface parameter. The perimeter sweep reaches Z = 279.99996 mm, and separate near-crown samples converge to principal radii of about 0.710 mm and 49.28 mm. The shape has a small transverse radius and the stripe spacing tightens through that turn.

Ten oblique sections run from Z = 220 mm to the perimeter at chord positions from 5% to 95%. Their lip radii range from 0.262 mm at 5% chord to 0.702 mm at 55% chord. The traces show curvature building outboard of the attachment before tightening at the lip. The images show the complete fin in two fixed zebra poses with CAD edge lines hidden.

## Reproduce

Run these commands from the trial directory with Python and the project CAD dependencies. This run used CadQuery 2.8, Open CASCADE 7.9.3, and Matplotlib 3.11.

```powershell
python build_fin.py
python check_fin.py
python render_fin.py
```

The check records the SHA-256 of the STEP it measures. Rebuilding changes the STEP header timestamp.
