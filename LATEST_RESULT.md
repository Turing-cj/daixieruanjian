# Latest consolidated state

Updated: **2026-10-01**

The authoritative source directory is this repository.

## Verified source state

- VERSION: 4.4.0
- Includes maintenance changes through 2026-09-30 in:
  - module47_batch_20260922.py
  - source_split.py
  - browser_circle_qa.js
  - associated regression tests
- Full validation on the authoritative Windows checkout:
  - python -X utf8 -m unittest discover -v
  - **128 tests passed**

## Stable commands

Single pathway:

~~~bash
python -X utf8 render_network.py FIT_RESULT.json MODEL.tsv OUTPUT.html --layout LAYOUT.json
~~~

Module-conditioned batch:

~~~bash
python -X utf8 render_batch.py prepare  --source SOURCE.zip --output OUTPUT_DIR
python -X utf8 render_batch.py render   --output OUTPUT_DIR
python -X utf8 render_batch.py validate --output OUTPUT_DIR
python -X utf8 render_batch.py export   --output OUTPUT_DIR
~~~

## Historical delivery records

The prior nine-path source/split delivery remains locally under
nine_paths_source_split_20260922.

The latest ED 15-path paper delivery was recorded at:

~~~text
E:\论文路径\ED15_20260930_latest
~~~

Those generated delivery directories are not authoritative source and are not
required for a clean clone of the renderer.

See AUTHORITATIVE_STATE.md for the complete repository contract.
