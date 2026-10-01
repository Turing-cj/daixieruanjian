"""Reaction-source and split-use annotations; never infer atom-level tracing.

Dots express a parent/portion relationship. Numbered roots express independent
producing reactions of one chemical pool. Multi-source consumption deliberately
has no fabricated parent assignment. Geometry and chemistry are not inputs to
this presentation-only transformation beyond the signed incidence partition.
"""
from __future__ import annotations

import math


LEGEND = ('[A.1], [A.2] = portions of [A]; [B1], [B2] = independent sources; '
          '[B1+B2:pool] = pooled use from those possible sources; exact allocation unresolved.')


def _pool_marker(code, source_ids):
    """Expose known possible sources without inventing source-to-use tracing."""
    if len(source_ids) > 1:
        return f"{'+'.join(source_ids)}:pool"
    return f'{code}:pool'


def derive_source_split(views, pool_codes):
    pools, annotations = {}, {}
    for mid, code in sorted(pool_codes.items()):
        members = {vid: view for vid, view in sorted(views.items())
                   if view['metabolite'] == mid}
        production = sorted(
            (inc['reaction'], vid, inc['signed_rate'])
            for vid, view in members.items() for inc in view['incidences']
            if inc['signed_rate'] > 0)
        sources = [{'id': code if len(production) == 1 else f'{code}{index+1}',
                    'reaction': rid, 'view_id': vid, 'amount': amount}
                   for index, (rid, vid, amount) in enumerate(production)]
        ids = [source['id'] for source in sources]
        produced = math.fsum(view['produced'] for view in members.values())
        consumed = math.fsum(view['consumed'] for view in members.values())
        tolerance = 1e-10 * max(1, produced, consumed)
        # A net import is an additional, unresolved source, not a free split.
        partition = len(sources) == 1 and consumed <= produced + tolerance
        children = [vid for vid, view in members.items()
                    if view['produced'] == 0 and view['consumed'] > 0]
        child_ids = {vid: f'{ids[0]}.{index+1}' for index, vid in enumerate(children)} if partition else {}
        uses = []
        for vid, view in members.items():
            own = [source for source in sources if source['view_id'] == vid]
            if own and view['consumed'] > view['produced'] + tolerance:
                # max(P,C) is a pooled demand here, not the amount made by own.
                role, marker, parent_ids = 'source_and_pool_use', _pool_marker(code, ids), []
            elif own:
                role, marker = 'source', '+'.join(source['id'] for source in own)
                parent_ids = []
            elif vid in child_ids:
                role, marker, parent_ids = 'split_use', child_ids[vid], ids[:]
            elif view['consumed'] > 0:
                role, marker, parent_ids = 'pool_use', _pool_marker(code, ids), []
            else:
                role, marker, parent_ids = 'inactive', f'{code}:0', []
            annotations[vid] = {
                'metabolite': mid, 'pool_code': code, 'role': role, 'marker': marker,
                'source_ids': [source['id'] for source in own],
                'source_reactions': [source['reaction'] for source in own],
                'parent_source_ids': parent_ids,
                'possible_source_ids': ids[:],
            }
            for incidence in view['incidences']:
                if incidence['signed_rate'] < 0:
                    uses.append({'view_id': vid, 'reaction': incidence['reaction'],
                                 'amount': -incidence['signed_rate'],
                                 'parent_source_ids': ids[:] if partition else [],
                                 'possible_source_ids': ids[:]})
        pools[mid] = {'pool_code': code, 'sources': sources,
                      'produced': produced, 'consumed': consumed,
                      'net': produced-consumed,
                      'single_source_partition': partition,
                      'boundary_import_required': consumed > produced + tolerance,
                      'consumption_allocations': uses}
    return {'version': 1, 'pools': pools, 'views': annotations,
            'scope': 'signed reaction-source bookkeeping; not atom tracing or separate physical pools',
            'legend': LEGEND}
