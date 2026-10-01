"""Resolve enlarged-label collisions against actual SVG metrics, not case IDs."""
from __future__ import annotations
import copy
from pathlib import Path
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parent


def clear_publication_labels(layout, html_path):
    """Return a copied layout plus an audit of label-only position changes."""
    out = copy.deepcopy(layout)
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(viewport={'width': int(layout['canvas'][0]), 'height': int(layout['canvas'][1])})
        page.goto(Path(html_path).resolve().as_uri())
        page.wait_for_function('''() => {
          const ps=[...document.querySelectorAll('path.link')];
          if(!window.freefluxAuditSnapshot || !ps.length || ps.some(p=>!p.getAttribute('d'))) return false;
          const geometry=ps.map(p=>p.getAttribute('d')).join('');
          if(window.labelGeometry!==geometry) {window.labelGeometry=geometry;window.labelStable=performance.now();}
          return performance.now()-window.labelStable>300;
        }''', timeout=20000)
        page.evaluate('document.fonts.ready')
        result = page.evaluate((ROOT / 'label_clearance.js').read_text(encoding='utf-8'))
        browser.close()
    if not result['paths_unchanged']:
        raise RuntimeError('Label placement unexpectedly changed reaction geometry')
    spec = out['circular_backbone']
    overrides = spec.setdefault('view_overrides', {})
    for change in result['changes']:
        view = change['id']
        position = {key: change[key] for key in ('label_dx', 'label_dy')}
        if view in spec.get('shared_instances', {}):
            group = spec['shared_instances'][view]
            group['label_offset'] = [position['label_dx'], position['label_dy']]
            for rid in group['reactions']:
                overrides.setdefault(group['metabolite'] + '_' + rid, {}).update(position)
        else:
            overrides.setdefault(view, {}).update(position)
    return out, result
