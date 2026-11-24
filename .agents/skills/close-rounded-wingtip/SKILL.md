---
name: close-rounded-wingtip
description: Build a rounded closure for a wing, horizontal stabilizer, or vertical fin using the existing skin edges, a 3D perimeter curve, and spanwise wrap guides. Verify the blend into both skins, the radius through the outer lip, and the continuity of zebra reflections on the exported geometry.
---

# Close a rounded wingtip or tail tip

**Continue the wing’s leading and trailing edges with one tangent 3D perimeter curve, then use spanwise wrap guides to shape the roll through that perimeter.**

## Build the cap

1. **Take the attachment curves from the parent skin.** Identify the parent revision and the station where the cap begins. Use the actual upper and lower cut edges as the first and last profiles of the boundary surface. On a vertical fin, use the left and right cut edges. Preserve the specified span, interfaces, and control-surface clearance.

2. **Draw one 3D perimeter curve for the outer lip.** Start at the leading-edge end of the cut and head outboard, tangent to the existing leading edge. Sweep around the tip and arrive at the trailing-edge end of the cut, tangent to the trailing edge in its inboard direction. Match both tangents in 3D. Use a plan view to shape the outline and a second view to check its elevation against the parent surface’s rise and twist. This curve is the middle profile of the boundary surface.

3. **Draw the spanwise wrap guides.** Pick a location along the chord. Start the guide on the upper cut edge, follow the upper skin’s slope outboard, turn smoothly through the perimeter, and return inboard to the lower cut edge following the lower skin’s slope. Give the turn through the perimeter a **finite lip radius** that meets the printable minimum in step 5. Each guide intersects all three profiles in order. Place additional guides at other chord locations where the roll needs more control, particularly as the section thins toward the trailing edge.

4. **Build the boundary surface.** Select **upper cut edge → 3D perimeter → lower cut edge** as the three profiles in the first direction. Use the spanwise wrap guides in the second direction. Match tangency to both parent skins at the attachment. Where the design calls for a smooth reflection, also match each skin’s local curvature and ease that curvature into the roll. If the CAD tool constructs upper and lower cap patches separately, give them a common tangent plane along the perimeter and a smooth, rounded turn through it.

5. **Keep every turn printable.** The minimum radius on the cap is **0.4 mm**, the width of the printer nozzle. A turn through the lip can be no wider than half the section's thickness, so the rule depends on the thickness of the parent section at the attachment, at the same chord position:

   - **At least 0.8 mm thick:** every radius along the wrap must be 0.4 mm or larger, at the lip and anywhere inside the cap.
   - **Thinner than 0.8 mm** (near a trailing edge, or the thin end of any section): make the lip a full half-round whose radius is half the local thickness, so it is as round as the section allows.

   Where the parent's trailing edge ends, end the cap the same way. On a sharp trailing edge the half-round shrinks to the parent's trailing-edge point. On a blunt trailing edge, the half-round closes the end of the flat trailing-edge face, so its radius there is half the trailing-edge thickness. Treat the 0.4 mm minimum as a design requirement from the start rather than a check after the build: it usually means spreading the roll over more of the cap reach.

6. **Shape the turn through the whole wrap.** For a simple rounded tip, shape each spanwise guide so the skin turns progressively from the attachment into the lip. Watch for a guide that makes a tight bend just before the perimeter and then flattens as it reaches the lip. If that happens, redistribute the bend along the guide while preserving the parent skin, the attachment, and the specified perimeter. A flatter crown or a cambered section can call for a different progression. Use the intended section shape to judge it.

## Verify the finished cap

1. **Check the reimported STEP.** Verify that the intended solid is closed, the retained parent skin is preserved, and the locked dimensions hold. At corresponding physical locations along each wing-to-cap join, measure separation and tangent-plane angle. Where curvature continuity is required, compare curvature in a common tangent frame. Check the shared perimeter separately if the cap uses two patches.

2. **Measure where the cap turns most tightly.** On the reimported STEP, take sections through the cap at guide locations and between guides. Include thin sections near the leading and trailing edges. On a swept wing, include cuts that follow the sweep so they trace the wrap from the parent skin out to the lip. Continue each trace around the lip and back to the opposite skin.

   Measure signed curvature along each trace. Record where the direction of bending reverses, the smallest sampled radius and its location, and the radius at the lip in that same section. Record the parent's attachment thickness at that section and judge the smallest radius against the printable rule in step 5 of the build. Sample more densely around any apparent reversal or tight spot to confirm it. If an intended simple roll tightens sharply before the lip and then opens out, revise the guides and repeat the section and zebra checks. Where the shape deliberately calls for an interior bend, document its location and purpose.

3. **Map the tightest turn over the whole cap.** Sections can pass between tight spots. Evaluate the larger principal curvature of every cap face on a dense grid of the reimported STEP, refine any local maximum, and compare each against the printable rule for its chord position. Report the smallest radius found, its location, and the attachment thickness used to judge it. Exclude collapsed endpoints where ordinary surface derivatives are undefined, and say how close to them the grid reaches.

4. **Run the numerical zebra check on the exported geometry.** Use CAD face normals and surface derivatives from the reimported STEP to compute the same reflection pattern as the viewer. Across both wing-to-cap joins, compare the stripe phase and its surface gradient. These show whether the bands meet and whether their spacing or direction changes abruptly. Locate actual light/dark crossings on the join and measure the change in stripe direction at the same crossing on both faces. Sweep multiple views and stripe angles. Report the largest mismatch and its location, the sampled region, and any endpoints where ordinary surface derivatives are undefined. Judge the results against the design’s tolerances.

5. **Review flagged locations in context.** Show the cap on the complete lifting surface. For a mirrored wing or tailplane, use a centered camera and mirrored stripe lighting. Hide displayed CAD edge lines and hold the camera, zoom, and lighting fixed when comparing revisions. Recheck suspicious triangular marks with finer CAD-normal sampling before attributing them to the surface.

The [Wing B rounded-tip study](../../../docs/studies/WING_B_ROUNDED_TIP.md) records the original build, the shoulder found after the first tangent-matched blend, and the measured effect of the curvature-matched revision. Its dimensions, guide placement, and sampling densities are specific to that wing.
