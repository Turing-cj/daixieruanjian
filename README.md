# Metabolic Pathway Renderer

Offline publication-oriented metabolic-network renderer built on the bundled
D3Flux / FreeFlux Studio stack.

**Current validated line:** v4.4.0 plus maintenance fixes through
**2026-09-30**. The authoritative checkout passed **128/128 regression tests**
on 2026-10-01.

For the exact consolidated project state, see AUTHORITATIVE_STATE.md.
Release history remains in RELEASE_NOTES.md.

## What this repository does

The renderer converts solved pathway chemistry and fluxes into an editable,
offline metabolic-network figure while preserving biochemical semantics.

Core input contract:

~~~text
model.tsv
fit_result.json
layout.json
~~~

Core output contract:

~~~text
HTML
SVG
PNG
Map JSON
manifest.json
~~~

The manifest records backend layout checks and the result of strict
actual-browser SVG QA.

For the module-conditioned research workflow, the upstream frozen export is:

~~~text
route_models.json
ranking_all.tsv
reaction_evidence.tsv
metabolite_names.json
manifest.json
~~~

The prepare stage converts that export into the standard renderer input
contract. It does not re-solve the pathway.

## Quick start

Python 3.11 is the currently verified runtime.

~~~bash
python -m pip install -r requirements.txt
python -m playwright install chromium
~~~

The repository vendors the modified d3flux Python package and browser-side
JavaScript needed for offline figures.

### Render one pathway

~~~bash
python -X utf8 render_network.py candidate01_demo/fit_result.json candidate01_demo/model.tsv candidate01_demo/output.html --layout candidate01_demo/layout.json
~~~

The command writes output.html, output.svg, output.png, output.json, and
output.manifest.json beside the requested HTML file.

### Run the module-conditioned batch

Use the stable CLI entry point:

~~~bash
python -X utf8 render_batch.py prepare  --source PATH_TO_EXPORT.zip --output PATH_TO_WORKDIR
python -X utf8 render_batch.py render   --output PATH_TO_WORKDIR
python -X utf8 render_batch.py validate --output PATH_TO_WORKDIR
python -X utf8 render_batch.py export   --output PATH_TO_WORKDIR
~~~

Focus on a subset with:

~~~bash
python -X utf8 render_batch.py render --output PATH_TO_WORKDIR --ids ID1,ID2
~~~

## Execution model

### Single-render path

~~~text
model.tsv + fit_result.json + layout.json
        |
        v
scripts_render_d3flux.py
        |
        +-- scientific_layout.py
        +-- circular_layout.py
        +-- branch_allocation.py
        +-- source_split.py
        +-- open_side_layout.py
        |
        v
standalone offline HTML
        |
        +-- actual SVG browser audit
        +-- label clearance when explicitly enabled
        |
        v
HTML + SVG + PNG + Map JSON + manifest
~~~

### Batch path

~~~text
frozen research export ZIP
        |
        v
render_batch.py prepare
        |
        v
standard renderer inputs + provenance JSON
        |
        v
render
        |
        v
strict visual QA
        |
        v
validate
        |
        v
vector PDF export + delivery validation
~~~

## Visual and biochemical invariants

The renderer is deliberately fail-closed around defects observed during
development:

- layout must not alter source stoichiometry or solved fluxes;
- every substrate/product incidence must be represented by the rendered scene;
- a declared circular backbone must remain a true equal-radius cycle;
- open output chains must not become artificial second rings;
- multi-endpoint reactions remain one grouped reaction event;
- shared pools must use real shared-node semantics, not disconnected look-alike
  copies;
- pathway quantity labels are independently recomputed from chemistry and flux;
- source/split labels distinguish real parent-child allocation from independent
  production;
- path continuity, crossings, overlaps, node hits, text collisions, clipping,
  arrow markers, and screen-axis distortion are checked on the actual SVG in
  Chromium;
- a stale success manifest cannot hide a new render failure.

## Important files

| File | Responsibility |
| --- | --- |
| render_network.py | Stable single-path render CLI |
| render_batch.py | Stable batch CLI |
| module47_batch_20260922.py | Validated module-conditioned batch implementation |
| scripts_render_d3flux.py | Main renderer orchestration |
| scripts_standalone_d3flux.py | Offline HTML and SVG/PNG export |
| scientific_layout.py | Publication layout and semantic roles |
| circular_layout.py | Circular backbone and pool/multiplicity logic |
| open_side_layout.py | Open side-chain geometry |
| branch_allocation.py | Local-flow amount allocation |
| source_split.py | Source/split relationship model |
| browser_circle_qa.js | Strict browser-side visual/semantic QA |
| label_clearance.py / label_clearance.js | Label-only clearance pass |
| BRANCH_QUANTITY_CONTRACT.md | Local-flow quantity contract |
| SOURCE_SPLIT_CONTRACT.md | Source/split semantic contract |

## Tests

Run the complete regression suite:

~~~bash
python -X utf8 -m unittest discover -v
~~~

The suite includes pure-Python invariants and real Chromium tests. Browser
tests render fresh scenes, inspect actual SVG geometry, deliberately inject
counterexamples, and verify Map JSON save/load behavior.

The repository commits only the small set of historical fixtures required by
those regressions. Bulk generated renders, PDFs, previews, scratch iterations,
and backups are intentionally excluded from Git.

## Version and provenance

VERSION remains 4.4.0 because the core release line is unchanged. The
post-release September 22-30 work strengthens batch layout, branch quantity
semantics, source/split semantics, connectivity, and visual QA without changing
the meaning of the v4.4 pathway-multiplicity release.

Dated filenames and fixture directories are retained where they encode
validated research provenance. Stable user-facing entry points are provided so
normal use no longer depends on remembering those historical dates.
