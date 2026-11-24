# Independent cap audit for the 0.4 mm radius trial

`audit.py` measures one whole-model STEP against its locked parent in `runs/cap-radius/parents/S/V/` and against the printable-radius rule of the close-rounded-wingtip skill (build step 5). Every number comes from the reimported STEP through CAD surface derivatives (`BRepAdaptor_Surface.D2`), never from a mesh. `compare.py` runs it over several models and writes one table.

~~~sh
.venv/bin/python runs/cap-radius/audit/audit.py SHAPE VARIANT MODEL.step [--out FILE.json] [--family native|legacy] [--grid N] [--quick]
.venv/bin/python runs/cap-radius/audit/compare.py --builds                       # every runs/cap-radius/builds/S/V/wing.step
.venv/bin/python runs/cap-radius/audit/compare.py --entry flat sharp a.step --entry wingb blunt b.step [--list entries.json] [--out-dir DIR] [--jobs 3] [--reuse]
~~~

SHAPE is flat, twisted, swept55, fin40, fin5 or wingb; VARIANT is sharp or blunt. `audit.py` writes FILE.json (default `results/SHAPE-VARIANT.json`) and FILE.txt (the summary below plus a per-cut table) and exits 0 only if every judged requirement passes. A full audit takes about 1.5 to 5 minutes on its own (82 s for the swept55 build); six run at once took 5 to 15 minutes each. `compare.py` writes `DIR/LABEL.{json,txt}` per model and `DIR/summary.{md,json}`; `--reuse` keeps an audit whose recorded STEP SHA-256 still matches.

## Conventions

- Millimetres. Normals point out of the solid (face orientation applied). Signed curvature is positive where the surface bends away from its outward normal (convex), negative where it is concave. Radius = 1/|curvature|.
- **Gate.** t is the parent's attachment-section thickness at the same chord fraction (`parent.py section_thickness`, perpendicular to the chord as the Converter measures it). Where t >= 0.8 mm the radius must be >= 0.4 mm (thick zone). Where t < 0.8 mm the score is radius / (t/2) and >= 0.9 passes (thin zone). Verdicts allow 1e-6 relative round-off, so an exact half-round at t = 0.8 mm (radius 0.4) passes.
- **Chord fraction of a cap point** is the fraction f of the cut-family plane that passes through it (solved on the exact family, not interpolated). The same assignment is used for the section cuts and the cap map.
- **Collapsed corners.** Caps built as boundary surfaces usually collapse one side of the (u, v) patch to the attachment's LE or TE point. Points where |Su x Sv| / (|Su||Sv|) < 1e-6, or |Su x Sv| < 1e-9 of the face's typical |Su||Sv|, have undefined derivatives and are excluded and counted. Points within 0.05 mm (`COLLAPSE_EXCLUSION_MM`) of a collapsed corner are measured and reported, in distance bands, but are left out of the map and join verdicts; the report gives how close the grid gets.

## Cut families

The native family is one plane per chord fraction f at each tip, built from `parent.py section_frame`: the plane passes through the attachment chord-line point at f, contains the attachment section's thickness direction, and contains the tangent of the parent's constant-chord-fraction line, d(chord point at f)/d(span), evaluated at the attachment. So each plane follows the local sweep, taper and dihedral of its own chord line, and the cap's continuation of that line.

| Shape | What the native plane at f is |
|---|---|
| flat | X = 120 f (identical to the earlier audit's planes) |
| twisted | through the twisted attachment chord point; contains the section's tilted up direction (-sin 3°, 0, cos 3°) and the chord line's sweep, taper and 5° dihedral rise |
| swept55 | vertical, through the chord point; swept at dX/d|Y| = tan 55° - (115/390) f, i.e. 55° at the LE, 48.6° at the TE. It contains the locked outline's point of chord fraction U(p) = f exactly, because the outline is the planform continued along these lines |
| fin40, fin5 | through (160 f, 0, 0) and (65 + 85 f, 0, 240): the fin's straight constant-fraction lines (identical to the earlier audit's planes) |
| wingb | vertical (contains Y, the section's up direction), through the definition-frame chord point at the 540 mm attachment, along the frame's constant-fraction line (sweep from the LE reference, taper, quadratic dihedral) |

`--family legacy` reproduces the planes of `runs/cap-trials/retest/audit.py` (calibration only): same as native for flat and the fins; parallel to the 55° LE for swept55 at every fraction; vertical planes through the root and attachment chord points in plan for twisted.

The task asked for planes "following the sweep" on swept55. The LE-parallel planes of the earlier audit follow the LE only; aft of mid-chord they run aft of the cap's own chord lines (the TE is swept 6.4° less), so near the TE a cut entering the cap at f = 0.96 leaves through the thinner trailing region and is judged against the wrong thickness. The native planes follow each chord line's own sweep. Both are available; the calibration below shows the difference.

Cut fractions per tip: 0.02, 0.05, 0.1, 0.25, 0.4, 0.5, 0.65, 0.8, 0.9, 0.95, 0.97, 0.99; both 0.8 mm crossings from `spec.json`; inside the thin zone the fractions where t = 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05 mm (sharp) or 0.7, 0.6, 0.5, 0.45 mm (blunt) on the TE side and t = 0.6, 0.4 mm on the LE side. Full-span wings are cut at both tips.

## What each measurement means

1. **Validity.** Solid and shell counts, shell closure (`BRepCheck_Shell`), `BRepCheck_Analyzer` standard and strict (exact). Extents come from bounded optimisation of the outboard coordinate over the cap faces (the `AddOptimal` box under-reads a single wrap face by 0.05 mm on the twisted model): projected span for the full wings, height above the root for the fins, half-span from the root plane for Wing B, against `spec.json` (tolerance 1e-3 mm). Volume uses `BRepGProp` with a 1e-9 tolerance; the OCCT default is reported beside it (on Wing B it reads 684,720 instead of 732,404 mm³).
2. **Parent preservation.** Retained skin: about 5,200 `parent.py` skin points (61 chord stations, cosine-clustered, by 43 span stations including both attachments; plus the blunt TE face) projected onto the model's non-cap faces, and about 3,800 points of those faces projected back onto `parent.step`; the larger maximum is judged (1e-4 mm). The projection falls back to a seeded foot-point minimisation, because OCCT's projectors stall where a skin closes to a thin TE (they returned 6e-3 mm on swept55 where the true distance is 3.6e-9 mm). Locked outline, on 101 extra cuts per tip (f = 0.005 to 0.995): for the mid-plane outlines (flat, swept55, fins) the outline must lie on the cut trace and be its outermost point along the plane's outboard direction; for the cambered twisted rim it must lie on the cut trace (its outermost point sits up to 0.13 mm off the rim by construction, reported); for Wing B's plan outline the outermost point in plan is compared. 3D outlines are also projected onto the cap faces (401 samples per tip). Tolerance 1e-3 mm.
3. **Joins.** Every edge shared by a parent face and a cap face, 399 cosine-clustered samples. At each, both faces are evaluated through their own pcurves at the same edge parameter: gap, angle between outward normals, and both second fundamental forms in one frame (e1 along the edge, e2 = N x e1); the curvature mismatch is the Frobenius norm of the 2x2 difference, and the cross-join normal-curvature difference (e2 e2 term) is listed separately. G1 is judged (gap <= 1e-4 mm, angle <= 0.01°); curvature is reported. Seams between curved cap faces get the same treatment and G1 judgement; seams with flat cap faces and parent-TE-face joins are reported under TE closure.
4. **Section traces.** Each cap face is cut by the plane (`BRepAlgoAPI_Section` with pcurves); every sample is pulled onto the exact plane by Newton steps in (u, v), so the trace is exact to 1e-13 mm. Edges are chained from the upper attachment, round the lip, to the lower attachment and sampled adaptively until the tangent turns less than 1.5° between samples (and samples are under 0.25 mm apart). Section curvature is the surface's normal curvature along the trace divided by |N x n_plane| (Meusnier), signed as above. Reported per cut: smallest radius (sampled peak refined by bounded 1D minimisation) and its position, arc-length fraction and principal radius; bending reversals (sign changes outside a 1e-5 /mm dead band, refined by root finding); tangent breaks at face junctions (a crease counts as radius 0); the lip radius; the outermost point; t; the verdict. **The lip** is the trace point at the locked outline (rim and mid-plane outlines) or the outermost point in plan (Wing B). Two context numbers are added for every cut: the **cut obliquity**, the angle between the plane normal and the outline tangent at the lip (0° = the cut crosses the lip square), and the **lip-normal radius**, the radius across the lip perpendicular to the outline at that point (summary line 4b, informational).
5. **Whole-cap map.** Each cap face on a (u, v) grid of 241 uniform plus 81 cosine-clustered lines per direction (about 100,000 points per face), trimmed by the face's own classifier. Larger-magnitude principal curvature from the shape operator I⁻¹II; chord fraction, t, required radius and score at every point. The twelve lowest local minima of score and twelve highest local maxima of curvature per face are refined by Nelder-Mead in (u, v). Reported: worst point (score), tightest point, worst per zone, all with location, chord fraction, t, margin or ratio, and distance to the nearest collapsed corner; then the near-corner bands.
6. **TE closure.** Sharp: the distance from the parent's TE point to the nearest cap vertex or collapsed corner (<= 1e-3 mm), and the principal radius within 0.1 mm of it. Blunt: the parent TE-face to cap joins (length, gap, angle), seams with flat closing faces, and the principal radius of curved cap faces within 0.1 mm of the 0.4 mm TE segment. Both: section cuts 1, 0.5, 0.2, 0.1, 0.05 and 0.02 mm ahead of the TE with their rule verdicts. Pass = closes and every approach cut passes.
7. **Summary.** One line per requirement with the governing numbers; OVERALL is the AND of the judged lines (4b is informational).

## Calibration

All files are in `calibration/`. Regenerate with `calibration/synthetic.py` and `calibration/calibrate.py`.

### Known-answer geometry (`synthetic.py`, `synthetic/synthetic_check.json`)

- **Sphere and torus** (NURBS-converted `BRepPrimAPI` solids). Principal curvatures against 1/R and against (1/r, cos φ/(R + r cos φ)) at about 500 points each: largest error 4e-15 /mm, with the concave inner torus correctly negative. Plane sections of the torus at z = 0, 1, 1.9 (outer circle convex +1/(R+w), inner concave -1/(R-w)) and a 40° oblique plane through the sphere: largest in-plane curvature error 8e-16 /mm.
- **Revolved caps on the flat parents.** Each attachment section's upper curve revolved 180° about the attachment chord line, sewn to the parent and NURBS-converted; on the blunt parent the 0.4 mm TE face is closed by a flat half disc of radius 0.2 mm. Every X = const cut is then a half circle of radius h = t/2, so the right answer is ratio 1 everywhere. Closed-form references come from `parent.py` jets.

| Check (sharp / blunt) | audit.py | Closed form |
|---|---|---|
| Section radius, 48 / 40 cuts | relative error <= 1.9e-9 / 2.8e-10; all ratios 1.0000; no reversals | h(X) |
| t at the cut vs 2h | <= 6.8e-10 / 3.6e-12 relative | |
| Principal curvatures, 3,042 points | relative error <= 4.2e-10 / 4.1e-10 | 1/(h√(1+h_x²)), -h_xx/(1+h_x²)^1.5 |
| Map worst point (refined) | 0.99132502132 / 0.99135282268 mm | minimum over the profile 0.99132502131 / 0.99135282264 mm (9e-12 / 3e-11 relative) |
| Join curvature mismatch, largest of 399 samples | 0.726078 / 0.726148 /mm | round-the-circle curvature at that sample, to 1e-10; its maximum between samples is 1.2e-4 higher |
| Lip vs locked outline | 44.03485 / 44.03464 mm | 44.03485 / 44.03464 mm (the revolved lip is not the locked outline) |
| Cap volume | relative error -5.8e-12 / -4.8e-12 | π∫h² dX |
| Span | 712.4826786386 / 712.4830287105 mm | 700 + 2 max h, to 7e-11 mm |
| TE | closes at the point; 0.02 mm ahead ratio 1.0000 / TE face joined to the half disc with gap 0 and angle 0, 0.02 mm ahead ratio 1.0000 | |

These models fail only the requirements they should (span and outline, which they do not try to meet).

### Accepted whole models against the earlier section audit (`legacy_comparison.md`)

On the earlier audit's own planes (`--family legacy`) and fractions, the smallest radius agrees on every cut of all five models to within 7e-11 relative (flat 0.387942 mm at 0.97, swept55 0.00075794 mm at 0.97, twisted 0.00288628 mm at 0.99, fin40 0.261767 mm at 0.05, fin5 1.88896 mm at 0.45), and the lip radius agrees to 9e-6 relative (twisted; 1e-13 or better elsewhere) once the lip is defined the same way, as the trace point at the locked outline. Differences:

- **Twisted reversals at 0.02 and 0.05** (earlier 1, now 0). The earlier sign change sits 0.0012 mm from the attachment, where |curvature| is 1.2e-6 /mm (radius 800 m). It is inside this audit's 1e-5 /mm dead band (the earlier one was 1e-8 /mm), so it is treated as straight. The other eight twisted reversals agree.
- **Twisted native planes** contain the section's 3°-tilted up direction and the chord line's dihedral rise instead of the vertical; the smallest radius changes by -0.4% to +1.0% (0.5257 vs 0.5210 mm at 0.02, 0.00291557 vs 0.00288628 mm at 0.99).
- **Swept55 native planes** follow each chord line's sweep instead of the LE. The smallest radius at 0.95 is 0.00817 mm (LE-parallel cut: 0.00206 mm); the map's worst point is 0.0385 mm at the 0.8 mm crossing (f = 0.960).
- **Wing B joins** against `runs/wing-b-rounded-tip-g2/check_g2.py` (85 points per side): gap 1.8e-8 vs 1.5e-8 mm, angle 2.2e-4° vs 3.5e-5°, curvature mismatch 1.2e-3 vs 1.7e-2 /mm. The earlier sampler went to chord fraction 1e-8 of the corners; this one stops 0.05 mm short of collapsed corners for the verdict (1.6e-3 /mm with those samples included), and compares the forms in the edge frame.

Accepted-model results with the native family (`accepted/summary.md`, per-model `accepted/*.txt`): all six are valid, at the spec extents, preserve the parent to 2.7e-7 mm and the outline to 7.7e-7 mm, and are G1 at the joins. Against the 0.4 mm rule (they predate it), the whole-cap map passes flat (worst 0.845 mm), fin40 (0.577 mm) and fin5 (0.668 mm) and fails twisted (0.048 mm at the TE-side crossing), swept55 (0.039 mm, same place) and Wing B (0.199 mm at x/c 0.0021 near the nose).

## Known limits

- **The section verdict depends on the cut, not only on the surface.** Near both ends of the chord the locked outline runs almost spanwise (it leaves the LE along the LE and returns along the TE), so constant-chord cuts cross the lip at 60° to 90° from square and show it as the tight end of an elongated section. The section radius then falls far below the surface's own radius across the lip: on the accepted swept55 cap at f = 0.02 the cut crosses the lip 86° from square and shows 0.130 mm where the lip-normal radius is 2.09 mm. Every accepted model, including the flat one whose map passes everywhere, fails item 4 in the thin zones for this reason; the synthetic revolved cap passes only because its cuts are exactly its circles. Read item 4 with the obliquity column, line 4b and the map. Item 4 is judged as specified.
- **The flat parent's TE and LE are vertical** (x ∝ u³ while z ∝ u), so its thin zones are within 0.026 mm (TE) and 0.002 mm (LE) of the ends, and cuts there are nearly tangent to the surface.
- **The map's gate is the attachment thickness at the point's family fraction**, not the cap's own local thickness. Near the TE of swept or tapered tips the cap is thinner than t(f) of its plane wherever the family plane leaves the planform, and the map then judges it against a larger allowance.
- **Points within 0.05 mm of a collapsed corner are not judged** in the map or the joins. They are measured and reported. Several accepted caps have a cone-like point there (principal radius down to 1e-7 mm within 1e-4 mm of the corner).
- **Reversals** use a 1e-5 /mm dead band; a sign change through a flatter region is not counted. Reversals exactly at a face junction are located at the junction.
- **Faces that straddle the attachment** are flagged and fail validity; their inboard part would otherwise be judged as cap.
- **Seams of a periodic cap face** (a face meeting itself) are listed, not sampled.
- **Retained-skin projection** is onto the untrimmed surfaces of the retained faces; closure and validity are what show that those faces cover the parent.
- **Wing B's outline is a plan curve**, so only the plan position of the lip is checked; its height is free.
- The parent `spec.json` extents and outlines are trusted as given.

## Files

| Path | Contents |
|---|---|
| `audit.py` | the audit |
| `compare.py` | batch runner and summary table |
| `calibration/synthetic.py`, `calibration/synthetic/` | known-answer models (two STEP files), their audits and `synthetic_check.json` |
| `calibration/calibrate.py`, `calibration/accepted/` | audits of the six accepted models and `summary.md` |
| `calibration/legacy_comparison.{md,json}` | earlier-plane comparison and the Wing B join comparison |
| `results/` | audits of build STEP files, keyed by SHA-256 (`builds-snapshot.md` lists what existed when it ran) |
