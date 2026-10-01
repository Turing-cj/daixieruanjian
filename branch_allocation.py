"""Quantity labels for actual rendered views, separate from chemical pools.

Contract: one signed reaction incidence belongs to exactly one displayed view.
For each view, amount=max(sum(production),sum(consumption)); do not add the
two sides of a through-flow. Global pool totals remain separate audit data.
This module changes presentation metadata only, never chemistry or geometry.
"""
from __future__ import annotations

import math


def apply_local_flow_labels(model, flux, scale, counts, balances, relationship_mode='pool'):
    groups = {}
    assignments = []
    canonical = {m.id: m.notes.get('map_info', {}) for m in model.metabolites}

    def register(view_id, mid, kind, info):
        group = groups.setdefault(view_id, {'metabolite': mid, 'kind': kind,
                                           'incidences': [], 'entries': []})
        if group['metabolite'] != mid or group['kind'] != kind:
            raise ValueError(f'conflicting quantity view identity: {view_id}')
        group['entries'].append(info)
        return group

    for mid, info in canonical.items():
        register(mid, mid, 'canonical', info)
    for reaction in model.reactions:
        info = reaction.notes.get('map_info', {})
        local, cofactors = info.get('local_instances', {}), info.get('cofactors', {})
        if set(local) & set(cofactors):
            raise ValueError(f'ambiguous local/cofactor quantity view: {reaction.id}')
        stoich = {m.id: c for m, c in reaction.metabolites.items()}
        if (set(local) | set(cofactors)) - set(stoich):
            raise ValueError(f'quantity view without biochemical incidence: {reaction.id}')
        for mid, coefficient in stoich.items():
            entry = local.get(mid, cofactors.get(mid))
            if entry is None:
                group = groups[mid]
            else:
                kind = 'local_instance' if mid in local else 'cofactor'
                view_id = entry.get('view_instance_id') or mid + '_' + reaction.id
                group = register(view_id, mid, kind, entry)
            rate = float(coefficient) * float(flux[reaction.id]) / scale
            group['incidences'].append({'reaction': reaction.id, 'signed_rate': rate})

    views = {}
    for view_id, group in groups.items():
        mid = group['metabolite']
        incidences = sorted(group['incidences'], key=lambda row: row['reaction'])
        produced = math.fsum(max(0, x['signed_rate']) for x in incidences)
        consumed = math.fsum(max(0, -x['signed_rate']) for x in incidences)
        amount = max(produced, consumed)
        record = {'metabolite': mid, 'kind': group['kind'], 'produced': produced,
                  'consumed': consumed, 'net': produced-consumed, 'amount': amount,
                  'pool_total': counts[mid], 'incidences': incidences}
        views[view_id] = record
        fallback = canonical[mid].get('pool_label_base', canonical[mid].get('pathway_label_base', canonical[mid].get('display_name', mid)))
        bases = {entry.get('pool_label_base', entry.get('pathway_label_base', entry.get('display_name', fallback)))
                 for entry in group['entries']}
        if len(bases) != 1:
            raise ValueError(f'inconsistent shared-view quantity labels: {view_id}')
        base = bases.pop()
        for entry in group['entries']:
            assignments.append((entry, base, record, view_id))

    # Assert the partition on both sides, never by summing max(P_view,C_view).
    for mid, expected in balances.items():
        for side in ('produced', 'consumed'):
            actual = math.fsum(v[side] for v in views.values() if v['metabolite'] == mid)
            if not math.isclose(actual, expected[side], rel_tol=1e-10, abs_tol=1e-10):
                raise ValueError(f'nonconserving view allocation: {mid}/{side}')

    presentation = model.notes.get('circular_backbone', {})
    provenance = None
    if relationship_mode == 'source_split':
        from source_split import derive_source_split
        provenance = derive_source_split(views, presentation.get('pool_codes', {}))

    # Commit only after all identities, incidences, labels and numbers validate.
    for entry, base, record, view_id in assignments:
        entry['pool_label_base'] = base
        relation = provenance['views'].get(view_id) if provenance else None
        if relation:
            identity = canonical[record['metabolite']].get('instance_base_name', record['metabolite'])
            base = f"{identity} [{relation['marker']}]"
            entry['flow_relation'] = dict(relation)
        else:
            entry.pop('flow_relation', None)
        entry.update(pathway_label_base=base, pathway_multiplicity=record['amount'],
                     pool_total=record['pool_total'], amount_scope='local_flow',
                     amount_view_id=view_id,
                     display_name=f"{record['amount']:.12g} {base}")
    result = {'enabled': True, 'label_mode': 'local_flow', 'normalization_flux': scale,
              'convention': 'max(local gross production, local gross consumption) / normalization_flux',
              'counts': counts, 'balances': balances, 'view_flows': views,
              'partition_verified': True,
              'legend': 'Numbers = local flow, not repeated pool totals.'}
    if provenance is not None:
        result.update(relationship_mode='source_split', source_split=provenance)
        presentation['instance_legend'] = provenance['legend']
    elif presentation.get('pool_codes'):
        codes = presentation['pool_codes']
        presentation['instance_legend'] = (' / '.join(f'[{codes[mid]}] {mid}' for mid in sorted(codes))
                                            + ': matching labels share the same pool.')
    model.notes['pathway_multiplicity'] = result
    return result
