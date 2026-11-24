# Locked parents for the 0.4 mm cap-radius trial

Each folder `S/V/` holds one locked parent: the part of an accepted model that is not cap, with no cap geometry. Shapes S are `flat`, `twisted`, `swept55`, `fin40`, `fin5` and `wingb`. Variants V are `sharp`, the parent as accepted, and `blunt`, the same parent with its trailing edge cut to exactly 0.4 mm.

| File | Contents |
| --- | --- |
| `parent.py` | Exact definition. It gives skin points and exact first and second derivatives in (chord parameter, span station), the blunt TE face, the cut, section thickness and the locked tip outline. The analytic parents need NumPy only. `wingb` reads its B-spline faces from `parent.step` through OCP. |
| `parent.step` | The parent alone as one closed solid with flat faces at the attachment(s). Wings with full-span accepted models are mirrored between the attachments. Wing B runs from the root to 540 mm. The fins run from the root to Z = 240 mm. |
| `spec.json` | Axes, extent, attachment station, cap reach and the locked outline with its parameters and samples. It also has the attachment chord and TE thickness, and the attachment thickness table (177 rows, clustered at both ends) with the chord fractions where the thickness crosses 0.8 mm. |
| `check.json` | Reimport validity and `parent.py`-to-`parent.step` fidelity. Sharp: deviation from the accepted model's retained skin. Blunt: TE thickness on the STEP and agreement with the Airfoil Converter. Also a jet-versus-finite-difference derivative check and the outline end gaps. |
| `build.json` | Build facts: fit sizes, pole counts, Wing B loft data. |

## Summary

| Parent | Attachment | Chord at attachment (mm) | TE at attachment (mm) | 0.8 mm thickness at x/c | Cap reach (mm) | STEP faces, strict valid, closed | parent.py vs parent.step (mm) | Sharp: deviation from accepted skin, parent.py / parent.step (mm) | Blunt: TE on STEP (mm) | Blunt: Converter agreement (mm) |
|---|---|---:|---:|---|---:|---|---:|---|---|---:|
| flat sharp | Y = ±350 | 120 | 0.000000 | 0.000017, 0.999784 | 50 | 4, yes | 5.10e-12 | 1.27e-13 / 4.64e-12 |  |  |
| flat blunt | Y = ±350 | 120 | 0.400000 | 0.000017, 0.999812 | 50 | 5, yes | 5.57e-12 |  | 0.400000-0.400000 | 1.19e-07 |
| twisted sharp | Y = ±440 | 100 | 0.000000 | 0.000435, 0.972039 | 60 | 6, yes | 8.79e-08 | 1.22e-06 / 1.63e-07 |  |  |
| twisted blunt | Y = ±440 | 100 | 0.400000 | 0.000427, 0.985897 | 60 | 8, yes | 5.15e-09 |  | 0.400000-0.400000 | 3.98e-08 |
| swept55 sharp | Y = ±390 | 95 | 0.000000 | 0.001021, 0.960383 | 60 | 6, yes | 1.49e-07 | 3.58e-09 / 1.89e-07 |  |  |
| swept55 blunt | Y = ±390 | 95 | 0.400000 | 0.001001, 0.979955 | 60 | 8, yes | 1.15e-07 |  | 0.400000-0.400000 | 2.77e-08 |
| fin40 sharp | Z = 240 | 85 | 0.000000 | 0.001540, 0.998460 | 40 | 4, yes | 7.62e-09 | 1.62e-12 / 8.59e-09 |  |  |
| fin40 blunt | Z = 240 | 85 | 0.400000 | 0.001540, 0.998845 | 40 | 5, yes | 8.40e-09 |  | 0.400000-0.400000 | 7.63e-07 |
| fin5 sharp | Z = 240 | 85 | 0.000000 | 0.001540, 0.998460 | 5 | 4, yes | 7.62e-09 | 1.62e-12 / 8.59e-09 |  |  |
| fin5 blunt | Z = 240 | 85 | 0.400000 | 0.001540, 0.998845 | 5 | 5, yes | 8.40e-09 |  | 0.400000-0.400000 | 7.63e-07 |
| wingb sharp | X = -540 | 87.5 | 0.000000 | 0.001315, 0.933747 | 60 | 3, yes | same surface | 1.64e-11 / 1.28e-11 |  |  |
| wingb blunt | X = -540 | 87.5 | 0.400000 | 0.001266, 0.972328 | 60 | 4, yes | same surface |  | 0.400000-0.400000 | 9.09e-06 |

Chord fractions are measured from each section's 2D origin, the leading-edge point, along the chord. Thickness is taken perpendicular to the chord line, as the Converter measures it. Sharp fidelity is the largest distance from about 1,700 sampled `parent.py` points (and 1,700 to 3,400 points sampled on `parent.step`'s curved faces) to the accepted model's retained skin faces. Blunt TE thickness was measured on the reimported STEP at five stations, including both ends: on the planar end-face edge at the attachment and root, and by plane sections in between. Converter agreement is the largest distance from the Converter's own output points to the analytic blunt section, at the root and the attachment. For `wingb` it is the largest distance from the Converter's output points to the blunt STEP skin, at five loft stations and at the 540 mm attachment (the attachment gives the 9.1e-6 mm). Every `parent.step` reimports as one solid with one closed shell. It passes standard and strict BRepCheck.

### Volumes and section areas

All volumes here use `BRepGProp` with a 1e-9 tolerance (`check.json` `reimport.volume_mm3`). Don't use OCCT's default integration, which is what CadQuery's `Solid.Volume()` and slab-cut volumes use. On the Wing B sharp parent it reports 589,733 mm³ instead of 711,013 mm³, 17% low, and that once made the blunt Wing B look 27% larger than the sharp one. `check.json` records both values under `volume_note`.

`check.json` also compares each reimported STEP section with the section it was built from (`section_areas`), at five stations:

- **Analytic parents:** the reference is the shoelace area of 40,001 `parent.py` points around the section.
- **Wing B:** the reference is the spline face through that station's loft points: the Converter output for blunt, the resampled SD7037 for sharp. The check also compares the thickness distribution with the Converter's points.

| Parent | Volume sharp → blunt (mm³) | Increase | Blunt/sharp section area, root → attachment | Largest STEP area error vs reference, sharp / blunt |
| --- | --- | ---: | --- | --- |
| flat | 864,000.0 → 864,047.8 | +0.006% | 1.0001 | 1.5e-6 / 1.5e-6 |
| twisted | 1,448,878 → 1,476,610 | +1.91% | 1.0153 → 1.0275 | 1.6e-6 / 1.6e-6 |
| swept55 | 1,164,280 → 1,191,108 | +2.30% | 1.0175 → 1.0387 | 1.6e-6 / 1.5e-6 |
| fin40, fin5 | 350,036.3 → 350,160.4 | +0.035% | 1.0002 → 1.0008 | 1.6e-6 / 1.6e-6 |
| wingb | 711,013 → 747,481 | +5.13% | 1.0413 → 1.0787 | 3.0e-4 / 1.7e-5 |

Each area ratio is close to the square of the keep-chord scale below, less the small cut-off tip. At the attachment that gives 1.0395 against 1.0387 measured for swept55, 1.0279 against 1.0275 for twisted, and 1.0798 against 1.0787 for Wing B. Wing B grows most because SD7037 has a thin, cusped trailing edge. The 0.4 mm cut therefore falls 2.0% of chord from the TE at the root and 3.8% at the attachment, which makes the scale 1.021 → 1.039. The 3.0e-4 on Wing B sharp is the accepted loft's own approximation of its sections (degree 7, 8 rows); the surface was left as accepted. On Wing B blunt the STEP thickness agrees with the Converter's points to 0.0045 mm or better at x/c = 0.01–0.95. Most of that is the Converter polyline's chord error near the nose.

Keep-chord cut of the blunt variants (cut chord fraction of the original section → uniform scale):

| Parent | Root | Attachment |
| --- | --- | --- |
| flat | 0.9999720 → 1.000028 | same (constant chord) |
| twisted | 0.9923934 → 1.007664 | 0.9863378 → 1.013850 |
| swept55 | 0.9912907 → 1.008786 | 0.9808224 → 1.019553 |
| fin40, fin5 | 0.9998915 → 1.000109 | 0.9996157 → 1.000384 |
| wingb | 0.9797566 → 1.020664 | 0.9623400 → 1.039138 |

The Converter matched the analytic rule closely. The scale differed by at most 7.1e-9, the cut chord fraction by at most 7.0e-9, and the TE thickness by at most 5e-13 mm. The remaining difference comes from the Converter's polyline interpolation of the 5,001 sampled section points.

## How the parents were recovered

- **flat.** The accepted builder's degree-5 chordwise polynomial is x = 120u³(10 − 15u + 6u²) with half thickness u(1 − u)(34 − 20u). It is extruded over Y = ±350. The outline is the accepted rim, `|Y| = 350 + 200u(1 − u)` at Z = 0.
- **twisted.** The supplied `parent_definition.py` is transcribed term by term into `parent.py`, including its 3D rim, `reach(u)`.
- **swept55.** The section is the analytic NACA 0009 closed-TE thickness in the accepted parameter q (U = 3q² − 2q³) with linear taper and 55° sweep. The lip outline is the parent planform continued to r = 390 + 240q(1 − q) at chord fraction U(q).
- **fin40 and fin5.** These use the exact 12% elliptical sections with linear sweep and taper and a flat root. The outline is `m ± a sqrt(1 − h³)` with h = (Z − 240)/40 or /5.
- **wingb.** The sharp parent is the skin surface of `01-wing-occt.step`, cut at X = −540 by exact B-spline segmentation, with no refit. The outline's anchors are taken from that surface by the accepted build's own procedure and frozen in `parent.py`. The blunt parent uses the same 49 Phase 0 loft stations (`_stations(wing, 4)`) and the same definition frames. Each section goes through the Converter's `blunt_trailing_edge(points, 0.4, keep_chord=True)` (revision d77bbdf). It is then closed by a straight TE line, as `_profile_wire` does, and lofted with the same `cq.Solid.makeLoft` call. Both faces are cut at X = −540 by segmentation.

Accepted models used for the sharp comparison (`check.json` names them by hash):

| Parent | Accepted model | SHA-256 |
| --- | --- | --- |
| flat | runs/cap-trials/retest/flat/wing.step | bc9ee86b…01caae5 |
| twisted | runs/cap-trials/retest/twisted/wing.step | 5b6bb3ae…743160f9 |
| swept55 | runs/cap-trials/retest/swept55/wing.step | 57132969…2c23e9 |
| fin40 | runs/cap-trials/retest/fin/40/wing.step | 0a30d2d6…28aaf1 |
| fin5 | runs/cap-trials/retest/fin/5/wing.step | ddf994a7…eaf57e |
| wingb | runs/wing-b-rounded-tip-g2/wing-b-rounded.step | 9a948464…2b3001 |

## Judgment calls

1. **The swept55 formula.** The accepted patch is a degree-14 polynomial fit to NACA 0009, accurate to 1e-9 of chord. `parent.py` uses the NACA formula itself. The two agree to 3.6e-9 mm on the accepted skin.
2. **The blunt chord parameter.** On every blunt parent, s = 0 at the leading edge and s = 1 at the cut. On the twisted wing, s maps to the native u as u = s + (u_end − 1)·max(0, 2s − 1)⁵ rather than linearly. This keeps the C1 camber-blend lines (x/c = 0.36 and 0.44) at fixed s, so the STEP skin can carry a knot line there. On the twisted parents those lines are C1 joints in the B-spline (multiplicity 4), placed within 1e-7 mm by knot removal.
3. **Scaling origin.** As in the Converter, keep-chord scaling is about the section's 2D origin, the leading-edge point of the definition, not the most forward skin point. Chord is measured as the Converter measures it, from the most forward point of the loop. On the cambered twisted section that point is 0.014 mm ahead of the LE point at 180 mm chord. The blunt TE therefore sits about 0.0001 mm forward of the sharp TE. On Wing B the resampled SD7037 nose is 0.0095 mm aft of the frame origin at the attachment. The blunt nose and TE therefore move aft by about 0.0004 mm, and the blunt TE sits at x/c = 1.0000043.
4. **Outline end corrections.** These apply to blunt twisted and wingb only. Because of item 3, the locked outline's ends missed the blunt attachment's LE point and TE face by 0.00011 mm (twisted) and 0.00037 mm (wingb). That is more than a 1e-6 mm sewing tolerance. Each blunt `parent.py` adds the exact end offsets, faded out over 18% of the outline parameter with the quintic weight already used by the Wing B fairing. The outline moves by at most 0.00011 mm (twisted) and 0.00037 mm (wingb), and its ends now lie on the parent. The sharp outlines are unchanged.
5. **Thickness at the twisted leading edge.** At x/c = 0 on the twisted section, the thickness is that of the outermost upper and lower crossings (0.25 mm at 100 mm chord), because the upper skin reaches slightly ahead of the origin. The table starts at x/c = 0. It ends at the TE position of the finished section: 1 for the symmetric parents, 0.9999989 for blunt twisted, and the measured TE for Wing B (1.0000043 blunt).
6. **Wing B blunt loft.** Run directly, the original pipeline loft (one spline through each station's Converter points) gives sections with different knot vectors, because the cut drops a different number of points at each station. ThruSections merged them into a 14,016 × 51-pole surface and a 76 MB STEP. Instead, each station's spline is refit once by least squares onto one shared cubic knot vector, with its end points exact. The refit stays within 9.6e-6 mm of the raw spline. The loft call is unchanged. OCCT then built a degree-3, 51-row loft in v; the accepted sharp loft is degree 7 with 8 rows. With sections that are not exact scalings of one profile, the same call takes its interpolating path. The TE gap of a lofted pair of 0.4 mm lines stays exactly 0.4 mm along the span.
7. **Wing B parent.py.** `wingb/*/parent.py` evaluates the B-spline faces in `parent.step` with OCP rather than a closed-form formula. Its chord parameter is the native u, running once around the section (`nose_u`, `side_range`), and `skin_jet` converts the v derivatives to projected-span derivatives. The definition frame (LE reference, chord, up) is exact and is given in closed form.
8. **One-sided root derivatives.** On the full-span wings, |Y| enters the planform, so span derivatives at Y = 0 are those of the right half. The twisted and swept parents have a C0 crease at the root, as in the accepted models. Their STEP files keep separate left and right faces.
9. **STEP fits.** Analytic skins are quintic tensor-product interpolants refined until midpoint errors are below 2e-7 mm, measured against `parent.py`. Over the sampled grids, `parent.py` and `parent.step` agree to at most 1.5e-7 mm.

## Targets not met

- The Wing B blunt common-knot refit (judgment call 6) is accurate to 9.6e-6 mm, and the Converter-to-STEP agreement is 9.1e-6 mm. All analytic parents agree with the Converter to better than 8e-7 mm. The goal of 1e-7 mm was not reached here. Holding it would have needed about four times the knots or the 76 MB merged surface.
- The Wing B blunt loft does not have the original's degree-7, 8-row spanwise structure (judgment call 6).
- Everything else met its target: validity, closure, TE thickness of 0.4 mm within 2.1e-12 mm everywhere, sharp fidelity of 1.2e-6 mm or better, and outline end gaps of 0.

## Regenerate

Run the following commands from `runs/cap-radius/parents/tools` with the worktree's `.venv` Python (CadQuery 2.8.0, OCP 7.9.3, NumPy, SciPy). The Airfoil Converter must be checked out at revision d77bbdf at `/home/m0obot/repos/Airfoil-Converter`; see `common.py`.

```sh
python compose.py         # analytic parent.py files from jet.py + shapes/*.py
python build_analytic.py  # analytic parent.step files (+ build.json)
python build_wingb.py     # Wing B parent.step and parent.py (both variants)
python make_spec.py       # spec.json
python check.py           # check.json (about 7 minutes)
python summarize.py       # the summary table above
```

STEP headers carry a timestamp, so rerun `check.py` after rebuilding. Each `check.json` records the SHA-256 of the STEP it measured.
