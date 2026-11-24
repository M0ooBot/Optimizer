
# Swept and washed-out wing

[wing.step](wing.step) contains one valid closed solid with six B-spline faces. The projected span is 1000 mm. The caps begin at ±440 mm and the root and cut chords measure 180 and 100 mm. The leading edge moves 90 mm aft at the cut. The retained wing has 5° dihedral and 3° washout.

## Build and check

Run the following commands from this directory using Python with CadQuery 2.8.0, NumPy, and SciPy. The optional image command also uses Matplotlib.

```powershell
python build_wing.py
python verify_wing.py
python render_zebra.py
```

The builder reads the supplied `parent_definition.py` without modifying it. The verifier records the hash of the actual export in [readback.json](readback.json). The two centered views in [zebra.png](zebra.png) use mirrored lighting and CAD normals from the STEP.

## Construction

Each cap is one surface that runs from the lower cut edge through the 3D rim and back to the upper cut edge. The wrap parameter `q` increases from −1 on the lower skin to +1 on the upper skin. Outboard distance grows as `reach(u) * (1 − q²)`. The thickness factor is `(15q − 10q³ + 3q⁵)/8`. Its first two derivatives vanish at each cut edge, so the cap can inherit the parent skin's position and first two span derivatives through the change of parameters. The quadratic span term and nonzero thickness derivative at the rim produce a finite turn radius at interior chord stations. The nose and trailing edge collapse to points at the ends of the cap patch.

The same chord knots interpolate every face. The parent uses 13 span stations and the wrap uses 65, both with imposed endpoint derivatives. The left wing is an XZ reflection of the right wing. The caps and four retained skin faces are sewn into a shell and converted to a solid.

## STEP readback

| Check | Maximum or result |
|---|---:|
| Sampled error against locked parent skin | 0.0000563 mm |
| Sampled error against locked 3D rim | 0.000000332 mm |
| Wing-to-cap separation | 5.64 × 10⁻¹¹ mm |
| Wing-to-cap tangent-plane angle | 1.60 × 10⁻⁷° |
| Difference in the shape operator | 1.94 × 10⁻⁶ mm⁻¹ |
| Reflection stripe phase difference | 1.92 × 10⁻⁸ cycles |
| Reflection stripe gradient difference | 1.45 × 10⁻⁵ cycles/mm |
| Stripe direction difference at 988 actual crossings | 0.000189° |

Join and stripe results cover both wing halves, both skin joins and `0.001 ≤ u ≤ 0.999`. The cap has singular surface coordinates at the collapsed endpoints. Zebra testing uses four perspective views and six stripe angles. The rim lies on the interior of the cap surface, so there is no separate face seam there.

Nineteen physical cut planes run from the lower attachment over the rim to the upper attachment. Each plane follows the local mean span direction at its cut station and contains the upper and lower cut points. The minimum radius is at least 99.0% of the radius at the rim in every measured plane. The highest-radius section at `u=0.32` has a 0.978 mm rim radius. The aftmost sampled cut at `u=0.98` has a 0.0000943 mm rim radius as the section closes toward the trailing edge. At every sampled plane the signed curvature reverses close to the lower attachment, where the inherited lower-skin curvature transitions into the wrap. The JSON gives each radius, location, and sign-change bracket. The model is checked geometrically and has no specified tooling-radius limit.
