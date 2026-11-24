# Airfoil Converter tools the winglet work needs

A running list for the planned Converter MCP server, kept as the winglet study finds them. Each entry says what the winglet needed and whether the pinned Converter ([d77bbdf](https://github.com/Trying-to-Fly/Airfoil-Converter/tree/d77bbdf45a8c0aa92efbc2a3d3007cfe4f904105)) already has it. Building the server is a separate piece of work.

| Tool | Winglet need | In the Converter today |
| --- | --- | --- |
| Blunt trailing edge, keep chord | Every section ends in a 0.4 mm trailing edge | **Yes**: `geometry.blunt_trailing_edge(points, thickness, keep_chord=True)`. It cuts on a polyline, so its output as a function of chord is only piecewise smooth. The winglet study samples it at a few chords and splines between them. A server tool could return the cut station and scale too, so a caller can build a smooth family. |
| Place a section in any plane | Section planes rotate through the bend | **Partly**: the GUI places sections on main planes, 3-point, 2-point and normal-to-line planes, with in-plane pitch. An agent needs one call that takes a leading-edge point, chord direction and up direction directly. |
| Split a section into camber and thickness | The morph works on camber and thickness, not points | **No.** Needs the true nose (a sampled outline rarely has a vertex there) and a common chordwise grid. |
| Morph between two sections | SD7037 into a symmetric section through the bend, with the trailing edge held in mm | **No.** Blend camber and thickness by a weight, then re-add the trailing-edge ramp at a fixed thickness. |
| Analytic and UIUC profiles | The target section (NACA 0010 thickness here) | **Partly**: the Converter reads airfoil-plotter CSVs and curve files. NACA 4-digit generation and a UIUC coordinate fetch would let an agent pick sections without a separate tool. |
| Leading-edge radius of a section | The 0.4 mm rule applies to leading edges when the trailing edge is blunt | **No.** Report the nose radius from a spline through the points, and the radius normal to a swept leading edge given the sweep. |
| Sections and curves along a law | Dozens of sections through the bend, with leading-edge, trailing-edge and skin guide curves | **Partly**: the Wing export writes sections, edge curves and skin guides for a wing whose sections are placed by the Converter. It does not take sections in rotating planes from a caller. |
| Loft in SolidWorks through those sections and guides | Stage 2 onward | **Partly**: `swloft` lofts profiles and guides. Untested with sections whose planes rotate through 75°. |
