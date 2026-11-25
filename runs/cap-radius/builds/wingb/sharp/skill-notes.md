# Skill notes: close-rounded-wingtip on Wing B (sharp TE)

Honest account for the skill maintainer. Build: CadQuery/OCP only, no SolidWorks.

## What I followed, and how

**Build 1 (attachment curves).** Used the parent skin's isocurve at v = 0.8986 (X = -540) as the upper and
lower cut edges; the cap boundary interpolates exact points on it (join separation <= 1.9e-5 mm).

**Build 2 (3D perimeter).** The plan shape is the locked outline; its elevation is the parent's camber
surface continued to second order (rise 60 f^2 plus section camber). Tangency to the LE/TE holds by
construction (the outline's ends are the parent's LE/TE continuation).

**Build 3 (wrap guides).** 239 guides, one per outline parameter, each from an upper cut station through
the perimeter back to the lower cut station at the same chord station. Not drawn by hand: each is a
closed-form profile Q = Qb - Q_O (1 - v/a)^m on the parent's continued surface.

**Build 4 (boundary surface).** OCC has no SolidWorks-style boundary surface with guides, so I
interpolated the guide net myself (one face, cubic across guides, quintic along the wrap with the exact
first and second wrap derivatives clamped at the cut). Tangency and curvature at the attachment come
from the base surface being the parent's second-order continuation.

**Build 5 (printable rule).** Built into the guide laws from the start: each guide's closure strength
is solved so that the smallest lip-section radius reaches 0.6 mm where possible, otherwise the roundest
G2 lip the thickness allows.

**Build 6 (progressive turn).** Closure exponent m >= 2.1 and a circle-like bound on the second-order
term, so the curvature rises into the lip instead of forming a shoulder.

**Verify 1-5.** All done on the reimported STEP: join separation, angle and curvature in a common frame;
110 plane sections (guide, between-guide, spanwise, constant-chord "sweep-following", chordwise) with
signed curvature and reversals; a 721 x 1201 principal-curvature map with refined minima; a numerical
zebra across the join (60 views x 8 stripe axes, phase, gradient, direction at 5029 crossings); and
rendered zebra views with fixed cameras.

## Where the skill was unclear, ambiguous or wrong

1. **The thin-section rule is not achievable as written for a curvature-continuous cap.** "A full
   half-round whose radius is half the local thickness" can only be an exact circle, and a circle
   meets the flat skin with only tangent continuity. Any G2 roll that drops by h has some radius
   below h. The same applies to sections just over 0.8 mm, where R >= 0.4 mm also means a circle.
   With my closure family the best G2 lip is about 0.96 h. On top of that, an inclined section
   (camber slope plus the wing's 0.18 rise) cannot hold a circle of half the vertical thickness at
   all; its perpendicular thickness is smaller. I met the rule by locally thickening the thin-TE
   guides next to the lip (a "lip bead", up to +20 % half-thickness within 0.1-2 mm of the lip). The
   skill should say which to give up (G2 inside the cap, the exact half-round, or the section
   thickness) or allow a tolerance such as 0.95 h.
2. **"Same chord position" is undefined for points inside the cap.** Guides cross chord stations, and
   the outline's LE/TE ends run along the continued LE/TE. I used the chord fraction of the parent's
   section frame continued past the attachment. Judging pointwise this way in a wedge-shaped TE
   region raises the requirement just inboard of the lip (the section gets thicker), so a lip exactly
   at t/2 fails a few percent there. That cost another 8 % of bead below 0.5 mm.
3. **The section-trace radius depends on the cut angle.** A plane that crosses a lip obliquely gives a
   plane-curve radius smaller than the surface radius by the obliquity. Taken literally, verify step 2
   fails any lip cut obliquely (my constant-chord cuts in the TE band read 0.04-0.15 of the rule). I
   judged sections by the surface's normal curvature along the trace and recorded the plane-curve
   radius separately. The skill should say which one it means.
4. **"Spread the roll over more of the cap reach" points the wrong way.** For a given half-thickness h
   and reach L, a roll spread evenly over the reach (an elliptical closure) has a lip radius of h^2/L:
   4^2/60 = 0.27 mm at mid-chord here, below 0.4. A larger lip radius needs a blunter closure that
   keeps the thickness longer and turns more compactly (R = m h^2 / 2L for Q = h^2 (1 - s^m)). The
   skill's point is presumably "don't concentrate the turn at a point", but the wording invites the
   elliptical closure, which fails.
5. **Camber/rise shear.** The vertical-tangent lip radius is not the smallest radius once the section is
   inclined: one side tightens by about (1 + m^2)^1.5 (m = camber slope across the lip). The design has
   to aim above the target. The skill does not mention this.
6. **"Sweep-following cuts"**: ambiguous between the LE sweep and the local constant-chord lines (which
   fan from 0.19 at the LE to -0.015 at the TE). I used constant-chord lines.
7. **Locked outline coincides with the LE/TE continuation** for the first and last 18 % of the outline.
   The cap's LE end is therefore the parent's nose continued, and on a sharp TE the half-round shrinks
   to zero along about 20 mm of outline, not only at the TE point. The skill's "shrinks to the parent's
   trailing-edge point" does not cover how it behaves along the TE end of the outline.
8. **Zebra metrics are underspecified.** "Stripe phase and its surface gradient" needs a definition. I
   used phase = N acos(r . a)/pi of the reflected view ray r. Stripe direction and relative gradient
   only mean something where the gradient is not tiny, so I floored them at 0.2 stripe/mm.
9. **Parent artefact.** The parent's nose has a 0.076 mm radius spot within about 0.001 mm of the nose
   point (0.7-0.9 mm just outside it). A G2 cap inherits it along its LE end. The rule passes there only
   because the chord station is in the thin zone. Worth a note in the skill: check the parent's own
   radii before blaming the cap.

## What I had to decide or improvise (not covered)

- The whole surface construction (no boundary-surface tool), its parameterisation (cosine-spaced guides,
  asinh wrap spacing about the lip) and the interpolation end conditions.
- The base surface past the cut (second-order continuation of camber and of squared half-thickness, with
  the planform thickness taper levelled off), and Q = h^2 as the variable so noses and lips stay smooth.
- The lip-radius target (0.6 mm, 1.5 x the minimum), the smooth law blending, and the lip bead.
- Guide stations near the nose: guides that start on the cut edge and run around the nose make an
  ill-conditioned corner (0.02 mm radius artefact). They must leave the cut spanwise.
- Exclusions at the collapsed LE/TE points (failures only within 0.001 mm of the LE point, and within
  2.2 mm of the TE point where the required radius is below 1 um).

## Iterations (about 10, all on the same construction idea)

1. First net with an arc-length/turning-angle wrap parameter and bicubic interpolation: lip radius rang
   0.40-0.73 mm along guides, and TE rows folded (radius about 0). The wrap parameter jumped in speed.
2. Cosine guide spacing and a single wrap parameter through the lip: ringing reduced.
3. asinh wrap spacing about the lip: folds gone; thick-region radii 0.49-0.64 mm.
4. Shear-aware targets and a thickness bulge for the thin band; first bulge also hit the nose
   (radius 0.2 mm on the LE strip); restricted to the TE.
5. Join tangent error 0.16 deg from interpolation end conditions: clamped quintic wrap direction,
   now <= 0.032 deg.
6. LE corner artefact (0.02 mm radius on the first guide): guides made to leave the cut spanwise.
7. Taper return with an overshoot added a concave band: replaced by a monotone tanh level-off.
8. Zebra review showed seams where the closure law switched regime (exponent jumping 45 -> 2.3 over a
   few guides) and a kink in the lip-radius target: smooth laws. No radius check had flagged these;
   only the rendered zebra did.
9. The bulge band (16 mm wide) produced zebra fingers and islands: replaced by the local lip bead.
10. Final build and full verification.

## Requirements: met, partly met, not met

| Requirement | Status | Numbers |
|---|---|---|
| One closed valid solid; parent skin, station, span, outline unchanged | met | valid, 0 free edges; skin poles identical; span 600.0000000001 mm; silhouette within 5.6e-6 mm of the outline |
| Tangent match at both joins | met | separation <= 1.9e-5 mm; angle <= 0.032 deg (<= 0.013 deg for 0.01 < xi < 0.99) |
| Curvature match / eased into the roll | met | cross-cut curvature difference <= 0.0025 /mm (rms 1-2e-4 /mm); closure starts C2 |
| Printable rule, t_att >= 0.8 mm: R >= 0.4 mm | met | min 0.4242 mm at xi 0.9337 (t_att 0.800 mm), p 0.835, (-573.08, 55.22, -149.75) |
| Thin sections: half-round of half the thickness | met with a lip bead, down to t_att about 1.5 um | min R / (t_att/2) >= 1.014 for t_att >= 5 um; unresolved within 2.2 mm of the TE point (t_att < 1.5 um) and 0.001 mm of the LE point |
| Ends like the parent (sharp TE: half-round shrinks to the TE point) | met | lip radius goes to 0 at the TE point |
| Progressive turn, no shoulder | met | minimum on each guide is 5-20 % below the lip radius, next to the lip (shear), no reversal on thick guides |
| Zebra continuity across the joins | met (no tolerance stated) | phase mismatch <= 0.0084 stripe width; direction change at crossings <= 4.9 deg (near the nose) |
| Exact half-round in the strict sense (circle) | not met by design | the lip is G2, so about 0.96 of a circle's radius before the bead |
