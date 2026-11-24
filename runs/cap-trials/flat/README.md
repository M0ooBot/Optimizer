# Straight, level wing with rounded tips

[Whole-wing STEP](flat-wing.step) · [Centered zebra](zebra.png) · [Lip sections](sections.png) · [Numerical STEP readback](readback.json)

The wing is a single mirrored solid spanning 800 mm. The central 700 mm has a 120 mm constant chord, with zero sweep, dihedral and twist. The outer 50 mm on each side is reserved for the cap. Coordinates are millimetres, with X aft, Y spanwise and Z up.

## Parent and cap

The study uses an original symmetric polynomial section with a rounded leading and trailing edge. Its thickness is 12.483 mm, or 10.4% of chord. The parameter `u` runs from zero at the leading edge to one at the trailing edge. The parent section has `X = 120u³(10 − 15u + 6u²)` and `Z = ±u(1 − u)(34 − 20u)`. The parent upper and lower surfaces extend these curves directly from Y = −350 to +350.

The build copies the actual parent control rows at the cut into the inner boundaries of the cap patches. The shared midplane perimeter on each side runs from leading edge to trailing edge. Its local outboard offset is `200u(1 − u)` mm, which peaks at 50 mm. Both perimeter endpoints align with the parent edge tangents. The plan curvature also tends to zero at those endpoints.

The tensor-product surface encodes a continuous family of fixed-aft wrap guides. For roll parameter `v` from zero at the cut to one at the rim, the guide scales the outboard run by `2v − v²` and the signed parent height by `1 − 2.5v³ + 1.875v⁴ − 0.375v⁵`. Height has zero first and second derivatives at the parent. The upper and lower patches have the same tangent plane at the midplane perimeter and turn smoothly around it.

## STEP verification

| Check | Result |
| --- | --- |
| Topology and size | One closed six-face solid passed the strict OCCT check. Chord and span were 120 and 800 mm within the 0.000001 mm study limit. |
| Parent-to-cap joins | Zero sampled separation, 9.54 × 10⁻¹⁵° maximum tangent-plane gap and 1.47 × 10⁻¹⁶/mm maximum curvature-tensor gap. |
| Shared perimeter | Zero sampled separation, tangent-plane gap and curvature-tensor gap. |
| Zebra | The worst phase gap was 2.22 × 10⁻¹⁵ cycles and the worst gradient gap was 1.99 × 10⁻¹⁵ cycles/mm, near the aft edge. The stripe-direction gap was zero at 3,924 solved light/dark crossings. |
| Lip sections | Curvature grew into the turn with its sampled peak at the perimeter. Rim radius was 1.266 mm at midchord and 0.089 mm at 119.88 mm aft. |

The numerical sweep covers 0.12 to 119.88 mm aft at 1,001 stations per join. It uses 84 camera and stripe combinations for each of the four skin-to-cap joins. The checker evaluates the actual endpoint direction of the outline. The exact leading- and trailing-edge vertices collapse the cap's surface parameters, so the two-direction derivative comparison covers the interior interval. The local fixed-aft lip radius tends to zero as the section ends.

The high-level loft interface provides less direct control over the boundary derivatives, so the build creates polynomial Bézier patches through OCCT and uses CadQuery for solid creation and STEP export. The static zebra views evaluate CAD normals directly in centered top and bottom projections. The numerical sweep supplies additional angles. The study evaluates surface geometry and imposes no manufacturing radius target.

## Reproduce

From the repository root, run with Python, CadQuery/OCP, NumPy, SciPy and Matplotlib installed.

```powershell
python -B runs/cap-trials/flat/build.py
python -B runs/cap-trials/flat/check.py
```

Run only `check.py` to verify the delivered STEP without regenerating it. The report records the hash of the STEP it inspected. A rebuild can change the hash because the exporter writes a timestamp.
