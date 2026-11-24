---
name: shape-wingtip-caps
description: Build or revise a rounded CAD wingtip cap, control its planform and rolled lip, and check wing-to-cap continuity using reimported surfaces and zebra reflections.
---

# Shape wingtip caps

**Establish the attachment and finish.** Identify the exact parent wing revision and where the new surface may begin. Trace the intended outline from leading to trailing edge in plan, then follow the parent wing's rise and twist in front view. Treat planform rounding and the radius seen in a section across the lip as independent decisions. Preserve the specified span, interfaces and control clearance. A photo can specify a silhouette without specifying an airfoil or a lip radius.

**Build from the skin.** For a rolled closure, use the actual upper and lower cut edges, an open three-dimensional perimeter, and guides that run from the upper skin around the lip to the underside. Control how the roll tightens toward the trailing edge. Where a fair reflection is important, sample the parent skin and carry its local slope and curvature into the transition. Ease those constraints out across the cap. Inspect the aft guides as the parent section approaches zero thickness. If the perimeter needs local fairing, measure and disclose its displacement against the intended outline.

**Verify the exported shape.** Reimport the STEP. Check the solid, the retained skin and locked dimensions. At corresponding physical positions on both sides of the join, measure separation, tangent-plane angle and curvature in a common tangent frame. Inspect several sections beyond the join for a ridge displaced outboard. Check the upper-to-lower rim separately. At a collapsed endpoint, report one-sided behavior and the closest regular sample. If a sewn boundary fails strict CAD validation, inspect agreement between the shared 3D edge and each face's parameter-space edge curve.

**Use zebra for comparison.** Show a symmetric design at full span with a centered camera and mirrored stripe lighting. Preserve camera, zoom and light direction while switching revisions. Hide displayed CAD edge lines and use surface normals from the CAD faces. Recheck suspicious triangular islands with finer sampling. For quantitative review, compare the continuous reflection field and its surface gradient across the join. To report a kink in an actual stripe, locate the same light/dark transition on each face. Record the STEP revision, sampled region, views and exceptions. Use the design's tolerances when making an acceptance decision.

See the [Wing B cap study](../../../docs/studies/WING_B_ROUNDED_TIP.md) for the build, a failed tangent-only blend, and the measured effect of its revision. Its dimensions, guide mapping and sampling densities belong to that fixture.
