"""Reusable publication sizing; keep all biochemical and node coordinates intact."""
from __future__ import annotations
import copy
import math


def enlarge_publication(layout):
    """Apply one idempotent readable profile to canonical and copied views."""
    out = copy.deepcopy(layout)
    style = out.setdefault('publication_style', {})
    met_font = float(style.get('metabolite_font_pt', 19))
    rxn_font = float(style.get('reaction_font_pt', 13))
    if not all(math.isfinite(v) and v > 0 for v in (met_font, rxn_font)):
        raise ValueError('publication font sizes must be finite and positive')
    met_scale, rxn_scale = 24 / met_font, 17 / rxn_font
    currency = set(out.get('hidden_metabolites', []))
    reactions = set(out['reactions'])

    def size(entry, *, cofactor=False, reaction=False):
        ratio = rxn_scale if reaction else met_scale
        for key in ('label_dx', 'label_dy'):
            if key in entry:
                entry[key] *= ratio
        if 'node_radius' in entry:
            entry['node_radius'] = 7 if cofactor else 12

    for mid, entry in out['metabolites'].items():
        size(entry, cofactor=mid in currency)
    for entry in out['reactions'].values():
        size(entry, reaction=True)
        for key in ('cofactors', 'local_instances'):
            for view in entry.get(key, {}).values():
                size(view, cofactor=key == 'cofactors')
    spec = out.get('circular_backbone', {})
    for key, entry in spec.get('view_overrides', {}).items():
        is_currency = key in currency or any(key == mid + '_' + rid for mid in currency for rid in reactions)
        size(entry, cofactor=is_currency, reaction=key in reactions)
    for collection in ('shared_metabolites', 'shared_instances'):
        for entry in spec.get(collection, {}).values():
            if 'label_offset' in entry:
                entry['label_offset'] = [v * met_scale for v in entry['label_offset']]
    for entry in spec.get('side_reactions', []):
        if 'label_offset' in entry:
            entry['label_offset'] = [v * rxn_scale for v in entry['label_offset']]
    style.update(metabolite_font_pt=24, reaction_font_pt=17, readability_profile='large_connected_v1',
                 auto_label_clearance=True)
    out['line_style']['width_range'] = [3.2, 7.6]
    return out
