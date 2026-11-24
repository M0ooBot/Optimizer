# Swept and washed-out rounded-tip wing

[Whole-wing STEP](twisted-wing.step) · [Whole-wing zebra](zebra.png) · [Fixed-aft sections](sections.png) · [Guide curvature](guide-curvature.png) · [STEP measurements](checks.json)

## Geometry

The mirrored wing uses a cambered NACA 2412 section. The thickness is applied normal to the camber line, with the closed-trailing-edge coefficient −0.1036. A quintic blend between the two camber parabolas over 36–44% chord fairs the camber curvature transition. From the root to the 440 mm cap station, taper, leading-edge sweep and washout vary linearly. The dihedral reference is the leading-edge chord datum, and the section rotates about that datum. The wing has the centerline crease implied by the mirrored panels and a sharp airfoil trailing edge.

| Quantity | Reimported STEP |
| --- | ---: |
| Projected full span | 1,000.000 mm |
| Root chord | 180.000 mm |
| Cap station, each side | 440.000 mm |
| True chord at cap station | 100.000 mm |
| Projected chord at cap station | 99.863 mm |
| Leading-edge datum sweep at cap station | 90.000 mm |
| Foremost physical skin point sweep at cap station | 89.996 mm |
| Leading-edge datum dihedral | 5.000° |
| Washout at cap station | 3.000° |

The upper and lower cut profiles use the parent skin section and its first and second spanwise derivatives at 440 mm. The perimeter continues the leading and trailing directions in 3D and reaches 500 mm projected semispan. The wrap guides reduce the skin thickness smoothly to zero at the perimeter while carrying the taper, dihedral, twist and camber through the cap. The upper and lower faces share one spline surface through the perimeter. The retained parent spline agrees with the analytical parent to 0.00049 mm in the sampled region.

## Readback

The exported STEP reimported as one closed eight-face solid and passed the strict Open CASCADE validity check. Across the four skin attachments, the sampled gap was at most 5.6 × 10⁻¹² mm, the tangent-plane mismatch was 1.5 × 10⁻⁷°, and the curvature-tensor mismatch was 2.8 × 10⁻⁶ mm⁻¹. The paired cap faces returned identical positions, tangent planes and curvature tensors at the sampled perimeter coordinates.

The numerical zebra sweep used 1,001 coordinates per join, seven cameras and twelve stripe angles per camera. The maximum phase mismatch was 3.6 × 10⁻⁸ cycles and the maximum gradient mismatch was 4.1 × 10⁻⁵ cycles/mm. Across 7,456 solved light/dark crossings, the greatest stripe-direction mismatch was 0.00188° at the left lower, near-nose attachment at approximately (90.010, −440.000, 38.436) mm. The sweep excludes the cap's collapsed end vertices where ordinary two-direction surface derivatives are undefined. The zebra images show centered views with mirrored stripe lighting and hidden CAD edges.

The fixed-aft sections trace both skins into the lip, including the thin aft section at 187 mm. Signed guide curvature shows a gradual turn and a lower-side inflection just outboard of the cut. The lip radius is about 1.76 mm at 20% chord, 1.15 mm at midchord, 0.141 mm at 90% and 0.0052 mm at 99%. The aft radius decreases toward the closed trailing edge. This is a geometry study with no minimum manufacturing radius specified.

## Reproduce

Use Python 3.12 with CadQuery 2.8.0, cadquery-ocp 7.9.3.1.1, NumPy, SciPy and Matplotlib. From this directory run:

```text
python build.py
python check.py
```

The check script reads the STEP from disk, records its SHA-256, and regenerates the measurements and figures. The SolidWorks bridge reported that SolidWorks was closed. Construction, STEP validation and reflection analysis ran in CadQuery and Open CASCADE through WSL.
