# Handoff: the Optimizer in the design loop

Status: **proposal, 2026-11-18. Nothing here is a decision in this repo yet.** It comes from a
design discussion with the user that spanned several projects, and from an independent review of
that design the same day. A session here reviews it against the code and `docs/MODEL_DETAILS.md`,
`docs/OPTIMIZATION_CONCEPT.md` and `docs/EQUIPMENT_PLAN.md`, decides, and records the result in
those documents. Until then they win. The user plans substantial changes to this tool, so this
file describes what the loop asks of it, not how the tool should be built.

This repo is now `MoZ200701/Optimizer`, a fork of `Trying-to-Fly/Optimizer` (2026-11-18), which
stays as it is.

Background, in other repos:

- [Engineering-Harness](https://github.com/MoZ200701/Engineering-Harness): `CHARTER.md`. The
  principles shared by all the projects, the architecture (Part 1b), the whole loop from
  specification to manufacture, the terms (Part 2a), and the review's questions and decisions
  (Part 2b, OQ-A1 to OQ-D4). Batch A is about this tool.
- [Plane-Aero](https://github.com/MoZ200701/Plane-Aero): `DESIGN.md`. The aerodynamics workflow:
  this tool, Plane-Builder and CFD.
- [Plane-Builder](https://github.com/MoZ200701/Plane-Builder) branch `handoff-terminal-first`, and
  the SW-Agent, SW-Agent-Sim and Plane Slicer documents linked from the charter.

---

## 1. What the user decided (direction, not design)

1. **Rough design, then verification** (OQ-A1). This tool produces the rough design. The user
   designs the detailed plane in Plane-Builder, which owns the shape and takes priority. That
   design comes **back here to be verified**: this tool checks the user's design rather than
   replacing it (charter P7). The imported-fuselage route (MODEL_DETAILS §7.5) is the starting
   point.
2. **New capabilities** (OQ-A2):
   - an **SAE Aero score objective** (the competition's scoring, under its weight limit);
   - **takeoff** and a **climb case** (MODEL_DETAILS §2 says "Climb, launch, and reserve are
     outside the model");
   - the **structural envelope**: g limits, speeds and weights the structure is checked against.
   **Loads per part do not come from here**; CFD and FEA produce them, *except wing structural
   loads (item 9, OQ-E6)*. **Mass budgets per part**
   are set by the user and tracked by the Systems tool, not produced here.
3. **The parts list is a shared data file that no tool owns** (OQ-A3). This tool and the
   planned Systems tool both read it. It starts from this repo's equipment-manifest format
   (MODEL_DETAILS §10, EQUIPMENT_PLAN), extended with electrical ratings, thermal limits and the
   source of each number. **Positions:** this tool's optimized positions (§10.3) are proposals
   until the user places a part in SolidWorks; from then the CAD position is the truth.
4. **The propulsion model is extended, as the one copy** (OQ-A4): battery internal resistance,
   wire losses, climb, and an honest ESC efficiency estimate instead of a constant. The Systems
   tool reuses it for the electrical budget and the heat that drives cooling studies.
   **There is no thrust stand,** and the design does not assume one; §2.3 becomes "no bench
   measurement" (item 13). A thrust stand may be considered in the future.
5. **One route for measured masses** (OQ-A5). Every mass goes into the Systems ledger (slicer
   grams, weighed parts, datasheets, allowances), and this tool calibrates from the ledger only.
   §1.4's reality factor would read from it.
6. **Corrections from CFD** (Plane-Aero §4): measured values replace this tool's assumptions
   where CFD can see them: the fuselage form factor, the interference margin, the prop-wash and
   tail-efficiency effects.
7. **Later, an MCP server** that Claude uses from the terminal, beside the others (charter A2).
8. **A "verify as-is" mode** (OQ-E3, second review): evaluate the design the user gives it,
   **wings and tails included**, without re-optimizing it. Today only a fuselage can be imported
   (MODEL_DETAILS §7.5, which also re-optimizes it), and the wing is one parametric family (§9.1),
   which may not represent a wing reshaped in Plane-Builder.
9. **Wing structural loads come from here** (OQ-E6): the span loading at each
   point of the structural envelope, with that point's trimmed state (attitude, elevator). CFD
   supplies local loads and spot-checks; near the stall it is not the source.
   *Corrected after the third review:* the in-loop model is a lifting line (MODEL_DETAILS §3.4),
   with spanwise lift but no chordwise pressure, and trim is solved at cruise only today, so trim
   per envelope point is a new capability. How the near-stall case is checked is open (OQ-H1).
10. **What a re-run produces** (OQ-F1): the verify-as-is verdict on the user's design, plus
    optional suggestions ("wing area +4% would gain 1.2 min"), which the user carries into
    Plane-Builder if they agree. It never overwrites the design.
11. **The parts list, the ledger and positions** (OQ-F2), refining item 3: the parts list says
    what each component is; Claude drafts it and every actual edit needs the user's OK. Masses
    live in a separate **ledger** of measurements (source and date, never overwritten), from which
    this tool calibrates; it keeps the weighed-vs-slicer pairs §1.4's reality factor needs.
    Positions come in 3D from the user's SolidWorks assembly once a part is placed. The manifest
    today is Python code (`Item(...)` objects with computed masses, x-only positions), so
    extracting it is a redesign of that part of this tool.
12. **The frame** (OQ-F3): the SolidWorks frame is the shared frame everywhere (mm; origin at the
    nose; +Z forward, +Y up, −X starboard). While this tool keeps its own frame, it converts at
    its own edge.
13. **Flight data** (OQ-F4): §2.3's "no measurement path" becomes "no bench measurement"; flight
    data is allowed once a plane flies. A motor-off glide gives drag; a powered cruise at the same
    speed then gives the propulsion chain's efficiency.
14. **After the third review** (OQ-H1, H2, H8): the spar this tool sizes (§1.3) goes into
    stage-1 FEA as a beam, at its size and position, until the user's detailed part declares the
    real one. The governing wing case is checked by a static load test of the printed wing, not
    by software, so this tool's span loading is checked for its *shape* only (against CFD at one
    moderate angle). Positions are measured from a fixed datum per plane, not the nose.

Charter-wide rules that apply here too:

- **Claude only, for now** (A1b); **agent instructions in `AGENTS.md`**, with `CLAUDE.md` reduced
  to the single line `@AGENTS.md` (A1c).
- **The architecture's terms are canonical** (Part 2a): "station" means Plane-Builder's
  definition sections; fuselage results are grouped by x-bands; "flight case" (aero) and "load
  case" (structure) are distinct. This tool's internal names are its own until it exchanges data.
- **Results say where their inputs came from** (C6): fingerprints of inputs, and what was not
  checked (C1), as this tool's *not checked by this model* column already does.

## 2. What it touches here

To decide, after reviewing the code:

- **§7 (the parametric loft) vs Plane-Builder's shape.** How the rough fuselage relates to the
  detailed one; what the brief sent to Plane-Builder contains; how the imported design is
  verified: *settled by item 8, checked as-is.*
- **§5 (objectives).** How the SAE score is written as an objective, and what takeoff and climb
  need from the aero and propulsion models.
- **The structural envelope.** What it contains and how it is exported (today the only
  structural case is the spar's root bending at about 5 g, §4).
- **§10 (the equipment manifest)** becoming the shared parts list: what moves out of this repo
  into the shared file, and how this tool reads it back.
- **§2 (propulsion).** The extensions in §1.4 above, and what a Systems tool needs from it.
- **§1.4 (calibration).** Reading measured masses from the Systems ledger.
- **Corrections from CFD.** How a measured form factor or interference value enters the model,
  and how it stays tied to the shape it was measured on.

## 3. What not to take from this file

- Nothing here changes current work or the open branches.
- The user plans to change this tool substantially; anything above can be met by a rebuilt tool
  as well as by this one.
- The direction in §1 is the user's; everything in §2 is for this repo to work out.
