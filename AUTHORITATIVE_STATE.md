# Authoritative renderer state

Last consolidated: **2026-10-01**

## Status

- Base release: **v4.4.0**
- Current implementation: v4.4.0 plus validated maintenance fixes through
  **2026-09-30**.
- Verification on the authoritative Windows checkout:
  **128 tests passed** with the complete unittest suite.
- The source directory was historically a directory bundle rather than a Git
  repository. This consolidation makes the source tree itself the authority,
  so old chat conversations are no longer required to reconstruct the software.

## Canonical entry points

### Single pathway

~~~bash
python -X utf8 render_network.py FIT_RESULT.json MODEL.tsv OUTPUT.html --layout LAYOUT.json
~~~

Required inputs:

- model.tsv — reaction stoichiometry
- fit_result.json — solved fluxes and optional concentrations
- layout.json — visual and semantic layout contract

Generated beside OUTPUT.html:

- HTML
- SVG
- PNG
- Map JSON
- manifest with backend and browser acceptance information

### Module-conditioned batch

Use the stable wrapper:

~~~bash
python -X utf8 render_batch.py prepare  --source SOURCE.zip --output OUTPUT_DIR
python -X utf8 render_batch.py render   --output OUTPUT_DIR
python -X utf8 render_batch.py validate --output OUTPUT_DIR
python -X utf8 render_batch.py export   --output OUTPUT_DIR
~~~

The validated implementation remains in module47_batch_20260922.py for
research provenance.

Batch source contract:

- route_models.json
- ranking_all.tsv
- reaction_evidence.tsv
- metabolite_names.json
- manifest.json

The prepare action converts those frozen research results into standard
renderer inputs without re-solving chemistry or changing fluxes:

~~~text
source export
  -> model.tsv
  -> fit_result.json
  -> layout.json
  -> source_pathway.json
  -> layout_meta.json
  -> input_validation.json
  -> native renderer
  -> HTML / SVG / PNG / Map JSON / manifest
  -> strict browser visual QA
  -> vector PDF export
~~~

## Correctness invariants

1. Reaction stoichiometry and solved flux values are never silently changed by
   layout code.
2. Every visible biochemical incidence must correspond to a real model
   incidence.
3. Circular backbones must satisfy their declared topology and radius
   constraints.
4. Multi-substrate and multi-product reactions remain grouped reaction events,
   not arbitrary pairwise chords.
5. Repeated metabolite views preserve explicit pool, source, and split
   semantics.
6. A render is accepted only when both backend checks and actual-browser SVG
   QA pass.
7. Stale success manifests cannot mask a failed new render.

## Important code

- render_network.py — single-render CLI
- scripts_render_d3flux.py — renderer orchestration and Map JSON generation
- scripts_standalone_d3flux.py — offline HTML, SVG, and PNG export
- scientific_layout.py — publication/scientific layout layer
- circular_layout.py — circular backbone, multiplicity, and relationship logic
- open_side_layout.py — open side-chain placement
- branch_allocation.py — local-flow quantity allocation
- source_split.py — source/split relationship semantics
- label_clearance.py / label_clearance.js — label-only repair
- browser_circle_qa.js — browser-side strict SVG audit
- module47_batch_20260922.py — validated batch adapter
- render_batch.py — stable batch CLI alias

## Repository fixture policy

Historical render directories contain many generated HTML, PNG, SVG, PDF, and
preview files. They are not source code and are ignored by default. Only the
small subset required by the regression suite is committed as explicit
fixtures.

This distinction is intentional:

- source, contracts, tests, and minimal fixtures belong in Git;
- scratch renders, delivery PDFs, previews, and backups remain local generated
  artifacts.

## Historical delivery locations

These locations are delivery records, not repository dependencies:

- General paper figures: E:\论文路径
- ED 15-path delivery: E:\论文路径\ED15_20260930_latest

The repository must remain usable even when those folders do not exist.
