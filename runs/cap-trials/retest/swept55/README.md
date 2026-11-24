# 55° swept wing with rounded tips

The STEP contains the complete symmetric wing as one closed solid with eight faces. The projected span is 900 mm and its reimported volume is 1,206,170.742 mm³. The level, untwisted parent uses a closed-trailing-edge NACA 0009 section, with the chord linearly tapered from 210 mm at the root to 95 mm at each attachment station, 390 mm outboard. Coordinates are X aft, Y span, and Z up.

## Reproduce

Run with Python, CadQuery 2.8, NumPy, and SciPy. The optional preview also uses Matplotlib.

```powershell
python build.py
python verify.py
python plot_zebra.py
```

The scripts write the files beside themselves. The STEP exporter records its export date. Each verification run records the hash of the STEP it tested.

## Construction

The parent patches use the same chordwise thickness polynomial at every spanwise station. The polynomial approximates the analytic NACA profile in the regular parameter `u = 3q² − 2q³`. The leading-edge tangent and closed trailing-edge endpoint are imposed during fitting.

For each fixed `q`, the upper cap starts on the parent at `|Y| = 390 mm`. Its spanwise progression is `r = 390 + 240q(1−q) f(v)`, where `f(v) = 1.6v − 0.8v³ + 0.2v⁴` for `0 ≤ v ≤ 1`. The wrap multiplies the local airfoil height by `g(v) = 1 − 2v³ + v⁴`. The derivatives of `g` through second order match a constant at the attachment, and `f′(1) = 0` gives the lip a vertical tangent through the wrap. The lower cap mirrors this surface in the chord plane. The construction preserves the prescribed lip at `v = 1`, including its endpoint tangents in plan. The two exact endpoint guides collapse where the cap meets the parent leading and trailing edges.

## Reimported STEP checks

| Check | Result |
| --- | ---: |
| Closed, valid solid | 1 solid, 1 closed shell, 8 faces |
| Projected span | 900 mm |
| Largest sampled parent NACA height and position error | 3.64 × 10⁻⁹ mm |
| Largest sampled locked lip position error | 6.48 × 10⁻¹² mm |
| Worst parent-to-cap gap | 2.74 × 10⁻¹² mm |
| Worst parent-to-cap tangent-plane angle | 2.98 × 10⁻¹¹° |
| Worst parent-to-cap normal-derivative mismatch | 3.11 × 10⁻¹⁰ mm⁻¹ |
| Worst shared-lip gap and tangent-plane angle | 2.16 × 10⁻¹⁶ mm, 0° |
| Worst shared-lip normal-derivative mismatch | 1.41 × 10⁻¹⁶ mm⁻¹ |

The numerical zebra check uses an orthographic reflection phase evaluated from CAD face normals and its spatial gradient from the surface derivatives. It samples 501 stations per wing-to-cap join from `q = 0.002` to `0.998`, with three view directions and five stripe angles mirrored for each half and skin. The largest phase mismatch was 3.03 × 10⁻¹² cycles and the largest spatial phase-gradient mismatch was 1.58 × 10⁻⁹ cycles/mm. Across 534 light/dark crossings, the largest stripe tangent mismatch was 2.10 × 10⁻¹⁰°. The reported error magnitudes are numerical residuals. The exact `q = 0` and `q = 1` cap endpoints have a collapsed wrap direction, so ordinary surface derivatives there are outside this comparison. The edge-free [zebra preview](zebra.png) uses lighting symmetric about the centerline.

The sections follow sweep in planes of constant normalized chord, `X − (tan 55° − 115u/390)|Y| = 210u`. Each trace is evaluated through both CAD faces, from the upper attachment around the common lip to the lower attachment. The largest plane residual was 5.24 × 10⁻¹¹ mm. At each of 27 section locations, including thin regions near the leading and trailing edges, the signed curvature kept the same bending direction and the smallest sampled radius was at the lip. The section at `q = 0.5` has a lip radius of 0.16059 mm. Sampled lip radii range from 5.49 × 10⁻⁷ mm at `q = 0.995` to 0.23923 mm at `q = 0.3`. The very small radii near the trailing-edge endpoint reflect the rapidly vanishing thickness of the specified parent profile. The radius results describe the CAD surface geometry.

[readback.json](readback.json) contains the sampled errors, worst-case locations, section measurements, and STEP SHA-256.
