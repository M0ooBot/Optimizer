---
name: shape-wing-edges
description: Shape or review a CAD wing's leading and trailing edges when the designer asks about sweep, taper, upsweep, a raked tip, or the smoothness of a wing silhouette.
---

# Shape wing edges

**Choose the silhouette before choosing the curve.** Trace both edges from root to tip in front and top views. For a brief that calls for a flowing wing, describe where each edge begins to turn and how its direction evolves. A trailing edge can remain nearly straight while the leading edge sweeps progressively aft. Preserve the designer's locked dimensions and interfaces, including any declared hinge axis. A deliberate panel break, straight inner panel, or localized tip bend can be the right response to a different brief.

**Choose geometry that expresses that intent.** A spline can be smooth at every join and still concentrate nearly all the bend near the tip. If the brief calls for a full-span curve, choose a construction that can turn across the full span. In the Optimizer Phase 0 tool, `span_curve="full_span"` is one such option. Extra sections can help a loft follow its guides between ribs. If the guide itself turns too late, revise its path.

**Use measurements when they clarify a decision.** If a bend looks abrupt, compare the direction or displacement of each affected edge at the root, one or more intermediate stations, and the tip. In plan view, departure from the initial tangent can help distinguish curvature from overall sweep. Choose stations that reveal the feature in question. Treat these measurements as diagnostic unless the designer specifies a numerical requirement.

**Inspect the rebuilt wing.** Compare both edges on the exported CAD geometry in true front and top views at a common scale. Inspect between defining sections and recheck the locks. Revise and recheck an unintended kink, ripple, or abrupt transition. If a lock prevents the requested shape, identify the conflicting feature. If only the curves or intended geometry are available, say which CAD checks remain open.

**Worked example.** Two SD7037 wings had the same 600 mm projected half-span and about 65 mm of tip rise. Wing 1 rose only about 3 mm by midspan, which made the outer bend look abrupt to the designer. For Wing 2 the agent kept the 15 section stations and tip position, chose a full-span guide, and checked the rebuilt leading edge at midspan. It rose about 16 mm there, giving the gradual silhouette the designer wanted. The midpoint was informative because the endpoints matched. Its value was specific to that design.
