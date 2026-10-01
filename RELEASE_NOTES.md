# V4.4.0

Pathway-level pool multiplicity without changing accepted circle geometry.

- Added opt-in counts from signed stoichiometric fluxes, normalized by an
  explicit positive scale. Canonical and local/shared carbon views agree;
  cofactor views retain their own chemical display names.
- Candidate 01 now shows 3 F6P, 2 G3P, 2 Xu5P, 3 Ru5P, 3 H6P and
  3 Formaldehyde. Other metabolites also show explicit counts, including 1.
  Counts are pool-throughput totals, not concentrations or branch amounts.
- Retained the original 12 reactions, stoichiometry, fluxes, all node positions,
  all reaction paths, the full 360-degree circle and pool codes. Only the shared
  G3P label offset moved 15 SVG units to avoid its longer label touching its node.
- Added independent browser recomputation of all visible numeric labels and
  label-versus-node clearance. A deliberately wrong count remains detectable
  even when copied consistently to every view. Failed acceptance exits nonzero.
- Added 16 regressions; full suite: 58 tests. Includes fractional/zero/reverse
  fluxes, nonunit stoichiometry, invalid data, normalization, idempotence,
  original chemistry, actual browser labels, deliberate numeric corruption,
  node/label collisions, and existing drag/download/reload regressions.

# V4.3.0

Connected regeneration and explicit same-pool relationships.

- Kept the eleven exact-radius outer arcs and unchanged biochemical inputs.
- Added validated shared actual metabolite nodes. TKT1 and TKT2 now feed one
  Xu5P node; RPE curves back to the actual Ru5P ring node instead of floating
  between unconnected local copies in the center.
- Added shared editable local carbon views. TAL output and TKT1 input reuse
  one G3P view; dragging it updates both reaction metadata entries and survives
  a real browser JSON download/reload without changing the outer G3P node.
- Added distinct `[A]` / `[B]` pool codes, with a visible key inside the figure.
  Candidate 01 now has three visible local carbon views, down from seven.
- Added opt-in consumption arrows terminating at the input reaction junction,
  including reverse-flux routing and marker-reference verification.
- Extended actual-SVG semantic QA to reject disconnected shared-pool clones
  and incorrect visible pool codes, not merely missing biochemical incidences.
- Side reactions reuse authoritative endpoint positions and reject attempts to
  silently move a backbone node off the circle. Legacy layout grammar remains
  available; the new relationship choices are explicit layout configuration.

# V4.2.0

Exact-circle, incidence-based rendering and fail-closed browser acceptance.

- Added an explicit, topologically validated `circular_backbone` layout mode.
  Candidate 01 uses a 475-unit radius, a full 360-degree closed backbone,
  eleven circular reaction arcs, and a central RPE reaction.
- Replaced paired substrate/product paths in this mode with one main trunk
  and one branch per additional endpoint. No duplicate shared trunks are
  painted. Consumption branches merge without misleading product arrows.
- Added explicitly marked, editable `local_instances` for repeated carbon
  metabolites; these remain distinct from nucleotide/water/proton cofactors.
  The biochemical model, all twelve reactions, and all input fluxes remain
  unchanged.
- Added browser-level incidence accounting and actual-SVG geometry QA.
  Checks cover screen scaling, radial error, crossings, long overlaps,
  labels, canvas clipping, node obstacles, and arrow semantics.
- Corrected an acceptance blind spot: the old audit ignored entire path
  pairs when they shared an endpoint. A regression now checks that the
  archived V4.1 spiderweb is detected as containing actual crossings.
  Earlier release notes' zero-crossing claims describe that old limited
  audit and must not be interpreted as visual acceptance.
- Added fail-closed CLI acceptance, real Chromium counterexample tests,
  and actual JSON download/reload tests.
- PNG export now uses the declared SVG size instead of the editor viewport.
- Retained the legacy layout mode and its regression suite.

# V4.1.0

Hyperedge-aware multi-metabolite routing update.

- Replaced visible pairwise routing for multi-substrate / multi-product
  reactions with a shared reaction bus: reactant branches converge locally,
  traverse one short shared trunk, then diverge to products.
- Preserved the existing circular / elliptical structure-first layout as the
  upstream geometry invariant; the new bus changes routing semantics rather
  than flattening the cycle.
- Prevented edge-channel and waypoint optimizers from splitting an arbitrary
  D3Flux substrate/product drawing pair away from its shared reaction bus.
- Extended reaction-hub optimization to repair obstacle-only defects, not only
  explicit edge crossings.
- Candidate 01 now renders five grouped reaction buses and passes browser
  geometry QA with zero sampled crossings and zero sampled obstacle hits.

# V4.0.0

Structure-first circular publication layout update.

- Added `ellipse` and `ellipse_arc` scientific-layout primitives so pathway
  geometry is generated as a cycle/arc skeleton before edge routing.
- Added exact angular-order search for small opt-in structures. Candidate 01
  automatically changes the PPP regeneration order from
  `Xu5P,E4P,S7P,R5P` to `Xu5P,S7P,E4P,R5P`, enabling the downstream edge
  optimizer to reach zero crossings without waypoint detours.
- Removed Candidate 01's horizontal `carbon_backbone` lane and replaced it
  with an outer carbon-cycle arc plus an inner regeneration ring.
- Cubic routes now receive continuous endpoint direction vectors rather than
  being restricted to N/S/E/W cardinal ports, reducing visible kinks on
  circular layouts.
- Publication QA now reports structure checks, near-edge clearance, and route
  stretch relative to the required reaction-hub path.
- Candidate 01 browser geometry QA passes with 0 sampled path crossings and
  0 sampled obstacle hits; 20 of 21 rendered paths are cubic.

# V3.2.0

Crossing-aware smooth routing update.

- Added exact bundle-channel search for small reaction bundles.
- Added per-edge corridor optimization when reaction-level channels are
  insufficient.
- Added optional local reaction-hub nudging as a constrained fallback.
- Added residual-crossing waypoint detours without moving biochemical nodes.
- Replaced piecewise independent cubic handles with a shared-knot smooth
  multi-segment cubic spline, including waypoint support.
- Route acceptance now minimizes a unified defect cost: crossings plus
  unrelated-node proximity, so an apparent crossing fix cannot simply create
  a new overlap elsewhere.
- Added opt-in flexible-node nudging as the final fallback for congested visual
  anchors.
- Added a zero-defect compaction pass: after crossings and node-proximity
  defects reach zero, the renderer contracts waypoints, corridor offsets, and
  flexible-node geometry only when all zero-defect invariants remain true.
- On Candidate 01 the compaction pass shortens the audited route network by
  about 6.8% (7947.97 to 7407.19 layout units) while keeping crossings and
  obstacle hits at zero.
- Candidate 01 route QA now reduces estimated inter-reaction crossings from
  2 to 0 and route obstacle hits to 0 while preserving stoichiometry and
  solved fluxes.
- Renderer manifests now distinguish cubic, quadratic, and linear paths.

# V3.1.0

Scientific-pathway publication layout update.

- Fixed layout authority: `layout.json` coordinates are no longer silently
  overwritten by the default saved Map JSON. Coordinate inheritance now
  requires explicit `--inherit-map-coordinates`.
- Added a thin `scientific_layout` compiler for semantic roles, ordered lanes,
  explicit anchors, local cofactors, reaction routing, and deterministic label
  placement without changing stoichiometry or fluxes.
- Added semantic edge routing (`straight`, `orthogonal`, `arc`) while retaining
  the native D3Flux quadratic renderer as a fallback.
- Added publication mode where flux magnitude is encoded by line width;
  low-flux reactions are no longer represented as dashed/uncertain links.
- Added publication typography, in-canvas legend placement, muted semantic
  colors, and compact reaction display labels.
- Added publication QA and regression tests for coordinate authority, semantic
  layout, cofactor localization, flux styling, and scene validity.
- Candidate 01 was rebuilt from the original `model.tsv` and `fit_result.json`
  as `candidate01_publication_v3_1_0.*`.
- Explicit UTF-8 JSON reading fixes Windows locale failures for labels such as
  `NADP⁺`, `H⁺`, and `H₂O`.

# V3.0.2

Final archived renderer for the DSM123 visualization test phase.

- Repacked from the latest renderer source and `test10.json` initial template

- Interactive flux, concentration, and network rendering
- Map JSON save, load, refresh recovery, and unsaved-change warning
- Pathway and metabolite palettes embedded in Map JSON
- Reaction flux labels and adjustable label sizes
- Delete, undo, restore, split, and merge tools
- Draggable circle guides with radius and opacity controls
- Metabolite and reaction snapping, movement along circles, and drag-to-release
- Offline HTML plus SVG and PNG export
- Optional saved Map JSON as the initial visual template
- Bundled `initial_template.json`, sourced from `test10.json`, is the default
