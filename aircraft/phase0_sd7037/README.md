# Phase 0 baseline wing

Start a candidate in a separate folder and keep this baseline for comparison. This is a single starboard half-wing with local origin at the root leading edge, span along −X, nose toward +Z and up along +Y. Root and tip are SD7037, 200 mm chord, 600 mm projected half-span, zero sweep, dihedral and twist. The planform area of this half is 0.120 m². It is a geometric test fixture with a continuous external skin; the spar channel and movable controls are later study features.

## Files

- [STEP solid](SD7037-Baseline.step) opens as a geometric reference. It has no native SolidWorks loft history.
- [Converter definition](SD7037-Baseline.airfoils.json) stores the two source-rib recipes, their placements, the wing recipe and the names of its generated sections and guides. The source is the Converter's bundled [SD7037 CSV](sd7037-il.csv), rescaled from 175 to 200 mm.
- [Section and guide curves](curves/) are the Converter's own millimetre XYZ exports. [Preview](preview.png) shows the half-wing and the common root/tip profile.

## Reuse in the next session

Use the Wing geometry routines from [Airfoil Converter](https://github.com/Trying-to-Fly/Airfoil-Converter/tree/d77bbdf45a8c0aa92efbc2a3d3007cfe4f904105). Load the sidecar with `store.load`, then call `wing_build.build_wing(record.spec, sidecar, record.stem, own_names=record.feature_names(), check=True)` for the wing record. For a working copy, point the rib records' `source` at the copied CSV and every record's `output_folder` at its own curves folder before exporting. Straight edges are expressed as blank edge-source fields. To explore curved edges, the Converter accepts guide files. Check its current geometric capabilities against the requested bend and section orientations.

The baseline has nine generated section curves at 75 mm intervals and 34 edge/surface guides. All 43 regenerated curves matched their saved files byte for byte. Each section is the same translated outline. For this constant-section case the STEP was constructed with Open CASCADE by interpolating the Converter's root curve and extruding it 600 mm. Reading that STEP back yielded one valid solid (volume 1,229,124.5 mm³). Sampled Converter section points on all nine stations were within 0.000001 mm of its boundary. A native Converter-linked SolidWorks loft and CAD readback are still pending.
