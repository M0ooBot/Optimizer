# Wing B (sharp TE) rounded tip cap

Rounded closure of the bowed SD7037 half-wing (root 200 mm, 600 mm projected half-span) from the
attachment at 540 mm projected span out to the locked tip outline. Built with OCP / Open CASCADE 7.9.3
(no SolidWorks). Parent skin, attachment station, span and outline are unchanged.

## Files

| File | What it is |
|---|---|
| `wing.step` | Whole model: parent skin + parent root face + cap, one closed valid solid (4 950 559 bytes, SHA-256 `18881aabebb7e63d17ba7329739efaf24dfa03552ed98ec935e3ab84ff836a20`). |
| `build_cap.py` | Builds the cap and writes `wing.step` and `build.json`. |
| `capgeom.py` | Cap geometry: parent access, extended planform, base surface, guide design laws. |
| `surfeval.py` | Vectorised surface derivatives / principal curvatures. |
| `verify_cap.py` | All checks on the reimported STEP; writes `readback.json`. |
| `render_zebra.py` | Zebra views (`zebra_*.png`) and the radius map (`radius_map.png`) from CAD normals of the reimported STEP. |
| `readback.json` | Every measurement below, with locations. |

Run from the project root, in order (about 10 s, 6-11 min, 1 min):

    .venv/bin/python runs/cap-radius/builds/wingb/sharp/build_cap.py
    .venv/bin/python runs/cap-radius/builds/wingb/sharp/verify_cap.py      # --quick for a coarse run
    .venv/bin/python runs/cap-radius/builds/wingb/sharp/render_zebra.py

## Construction

Coordinates are the STEP's (X = -span, Y up, Z forward). The cap is laid out in the parent's
section frame continued past 540 mm: S = projected span, xi = chord fraction of the frame
(LE reference Z = -(35f + 45f^2), chord 200 - 125f, f = S/600). The locked outline's faired ends are
exactly this continuation, so its LE end is the extended LE and its TE end the extended TE.

1. **Profiles.** Upper and lower cut edges are the parent skin's isocurve at the attachment
   (v = 0.8986, X = -540). The perimeter is the locked outline in plan, at the height of the
   parent's camber surface continued to second order (rise 60 f^2 and section camber), so it
   follows the parent's rise.
2. **Base surface.** Camber Cm and squared half-thickness Q = h^2 are taken from the parent at 540 mm
   with their first and second span derivatives at fixed xi (exact chain rule through the B-spline).
   Q (not h) is continued, so the round nose and the sharp TE stay smooth. Cm is continued exactly;
   the thickness taper is levelled off with a tanh (O(d^3) from the parent, net -1.4 % half-thickness),
   so the cap is curvature-continuous with the parent at the cut.
3. **Guides.** One wrap guide per outline parameter p (239 interior guides, cosine spaced, plus the
   collapsed LE and TE points). Guide p starts at the attachment station xi_a(p) (about p; quadratic
   near the nose so the guides next to the LE leave the cut spanwise), runs straight in (S, xi) to the
   outline point O(p), turns through the lip and returns to the lower cut at the same station. Along
   it Q(v) = Qb(v) - Q_O (1 - v/a)^m for v < a (v = plan distance from the lip), so the closure
   starts curvature-continuously and the lip has a finite radius.
4. **Lip radius laws (printable rule built in).** For each guide the closure strength is solved so
   that the smallest radius of the lip section normal to the outline, including the shear from the
   camber slope (which tightens one side of an inclined lip by about (1 + m^2)^1.5), equals a target:
   0.6 mm where the section is thick enough, otherwise the roundest circle-like G2 lip the local
   thickness allows (exponent m >= 2.1, second-order coefficient <= circle). Where the parent's nose
   already gives a rounder lip the target blends into it. All laws are smooth in p.
5. **Lip bead (thin TE band only).** A curvature-continuous lip cannot reach a radius equal to the
   half-thickness (only a true circle can, and that is only tangent-continuous). To still give
   R >= 0.4 mm where t_att >= 0.8 mm and R >= t_att/2 where thinner, the guides whose lip lies where
   t_att < 1.3 mm swell within 3 closure lengths of the lip (half-thickness +12 %, ramped in between
   1.3 and 0.8 mm; +8 % more below 0.5 mm). This is a deliberate interior bend (a shallow concave
   band just inboard of the lip, see Sections).
6. **Surface.** The guide net (241 x 201 points; wrap parameter asinh-spaced about the lip) is
   interpolated by one B-spline face: cubic across guides, quintic along the wrap with the exact
   first and second wrap derivatives of the base surface imposed at both cut edges. The face is
   degenerate at the LE and TE points. It is sewn (tolerance 1e-4 mm) to the parent's skin and root
   faces in place of the parent's flat attachment face.

## Measured on the reimported STEP (readback.json)

All numbers from `verify_cap.py` (full run: principal-curvature grid 721 x 1201, 785 join stations per
side, 110 section traces).

**Solid.** 1 solid, 1 shell, 3 faces (parent skin, parent root plane, cap), BRepCheck valid, no free
boundary edges. Volume 731 173.83 mm^3 (parent alone 711 012.72 mm^3, cap 20 161 mm^3; 1e-9 integration
tolerance). STEP 4.95 MB.

**Parent preserved.** Skin face poles identical to parent.step (max pole difference 0, sample deviation 0);
root plane X = 0 unchanged.

**Locked dimensions.** Max projected span 600.0000000001 mm, at (-600.000, 62.225, -117.054) (the
outline's tip point). Cap plan silhouette outside the locked outline by at most 5.6e-6 mm; lip trace
on the outline within 4.9e-6 mm. Attachment edges at X = -540 within 7e-11 mm.

**Wing-to-cap join** (same physical points, frame: along the cut / across the cut / twist):

| | upper | lower |
|---|---|---|
| max separation | 1.8e-5 mm | 1.9e-5 mm |
| max tangent-plane angle | 0.032 deg (xi 0.0027) | 0.025 deg (xi 0.0009) |
| same, 0.01 < xi < 0.99 | 0.0045 deg | 0.013 deg |
| max cross-cut curvature difference | 0.0011 /mm | 0.0025 /mm (xi 0.0004) |
| same, 0.01 < xi < 0.99 | 0.00026 /mm | 0.00022 /mm |
| rms cross-cut curvature difference | 1.0e-4 /mm | 1.7e-4 /mm |
| max twist-curvature difference | 0.0056 /mm | 0.0066 /mm |

The largest differences are within 0.25 mm of the round nose, where the parent's curvature is up to
13 /mm.

**Zebra across the join** (24 stripes per 180 deg, 60 views x 8 stripe axes, CAD normals and shape
operators on both faces, 5029 light/dark crossings on the join): max stripe-phase mismatch 0.0084 stripe
widths (upper, xi 0.0027); max change of stripe direction at a crossing 4.9 deg and max relative change
of stripe gradient 0.15, both at the lower side xi 0.0093 (0.8 mm behind the nose, gradients 0.72 vs
0.63 stripe/mm); max absolute gradient change 0.39 stripe/mm where the gradient is 7.5 stripe/mm (nose).
Direction and relative gradient are evaluated only where both gradients are >= 0.2 stripe/mm. The same
numbers hold with 0.1 mm and 0.5 mm excluded around the collapsed points. Views: `zebra_top.png`,
`zebra_bottom.png`, `zebra_oblique.png`, `zebra_tip-on.png` (fixed orthographic cameras, no edges).

**Printable radius (larger principal curvature, whole cap, refined minima).**

| attachment thickness at the point's chord station | samples | min R / required | radius there |
|---|---|---|---|
| >= 0.8 mm (R >= 0.4 required) | 570 305 | 1.066 | 0.4265 mm |
| 0.4-0.8 mm (half-round) | 47 198 | 1.050 | 0.2437 mm |
| 0.1-0.4 mm | 61 955 | 1.056 | 0.2084 mm |
| 0.02-0.1 mm | 46 665 | 1.037 | 0.0111 mm |
| 0.005-0.02 mm | 32 925 | 1.014 | 0.0027 mm |
| < 0.005 mm | 106 873 | 0 (collapsed points) | - |

- Smallest radius where the 0.4 mm rule applies: **0.4242 mm** (refined) at (-573.079, 55.224, -149.750),
  outline parameter p = 0.835, chord station xi = 0.9337 where the attachment thickness is 0.800 mm,
  on the lower side 0.3 mm inboard of the lip. On the LE side the smallest such radius is 0.505 mm
  (xi 0.0014, S = 543.1).
- Failing samples: 1113 of 865 921. 963 lie within 0.001 mm of the collapsed LE point; the other 150
  lie within 2.2 mm of the TE point where t_att < 1.5e-6 mm (required radius < 0.8 um). The grid reaches
  5.6e-6 mm from the LE point and 2.1e-6 mm from the TE point.
- Lip radius along the outline: about 0.6-0.69 mm over the thick part (designed 0.6 mm minimum
  including shear), 0.8-2.4 mm on the LE side (0.05 < p < 0.25) where the parent's nose makes the lip
  rounder, falling to the half-round of the local thickness (0.37 mm at p = 0.854, 0.037 mm at
  p = 0.919) and to zero at the TE point.

**Sections (110 traces: 50 vertical planes along guides and halfway between guides, 18 spanwise
Z = const and 18 constant-xi (sweep-following) planes from xi 0.0005 to 0.997, 14 chordwise
X = const planes from 540.5 to 599.95 mm).** All pass the rule on the surface radius along the trace.
Examples (smallest radius / lip radius, mm): guide p = 0.47: 0.636 / 0.677; p = 0.73: 0.514 / 0.599;
p = 0.85: 0.298 / 0.370 (t_att 0.56 mm at the minimum); spanwise xi 0.3: 0.543 / 0.625; chordwise
X = -599.95: 0.648 (LE) / 0.669 (TE). The smallest radius on a guide is 5-20 % below its lip value and
sits next to the lip on one side; that is the camber/rise slope shearing the lip (accounted for in the
design), not a shoulder. Bending reversals: none on guides with p between 0.03 and 0.80 (their only
concavity is the parent's bow plus the levelling-off of the thickness taper, radius > 440 mm, below the
counting threshold of 500 mm); 4-8 on the thin-TE guides
from the lip bead (concave radius 14 mm at p = 0.85, 1.2 mm at p = 0.92); 2-7 on traces through the
LE micro-nose (xi < 0.0002) and the parent's own concave aft lower surface in chordwise cuts. See
`readback.json` `sections.traces[].reversal_locations`.

## Decisions and limits

- **Chord position** of any cap point is its xi in the continued planform frame; the printable rule
  compares the larger principal curvature there with min(0.4, t_att(xi)/2), t_att from the parent's
  attachment section (spec.json definition).
- **Section radii** are reported as the surface's normal-curvature radius along the cut. The plane-curve
  radius of a cut that crosses the lip obliquely is smaller by the obliquity and says nothing about the
  surface; it is recorded as `smallest_plane_curve_radius_mm`.
- **Collapsed points.** The face collapses at the LE and TE points. Below 1 um from the LE point the
  fitted surface's curvature is not meaningful; near the sharp TE the half-round radius shrinks with
  the section and is not resolved below about 1 um (t_att < 1.5 um, within 2.2 mm of the TE point).
- **Parent nose.** The parent's own nose has a 0.076 mm radius at the nose point (radius 0.7-0.9 mm
  0.001 mm away). The cap continues it along the LE end of the outline; by the chord-position rule it
  is within the thin (t < 0.8 mm) zone.
- **Bead** makes the cap's half-thickness up to 20 % larger than the parent's at the same chord
  station, within three closure lengths of the lip (about 0.1-2 mm), only where t_att < 1.3 mm.
