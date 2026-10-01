"""Unfold an ordered side pathway; never alter the biochemical model or main ring.

Repeated end pools retain their canonical IDs and explicit pool-code view labels.
A primary outlet feeding an open chain is a shared physical node, never a
detached display clone. No decorative return arc or invented reaction is added.
"""
from __future__ import annotations
import copy
import math
from collections import Counter


def unfold_side_chain(layout, reaction_ids, stoichiometry, protected_metabolites=()):
    """Return a layout copy with one continuous ordered trunk per reaction.

    ``stoichiometry`` uses the renderer's metabolite aliases and solved reaction
    orientation. It is read only, for positioning side incidences by their sign.
    The native compiler remains responsible for checking chemistry and coverage.
    """
    out = copy.deepcopy(layout)
    spec = out['circular_backbone']
    ids = list(reaction_ids)
    if not ids:
        return out
    if len(set(ids)) != len(ids):
        raise ValueError('side-chain reaction IDs must be unique')
    by_id = {s['reaction']: s for s in spec.get('side_reactions', [])}
    if any(rid not in by_id for rid in ids):
        raise ValueError('side chain references an unknown side reaction')
    chain = [by_id[rid] for rid in ids]
    for a, b in zip(chain, chain[1:]):
        if a['target'] != b['source']:
            raise ValueError('side chain must follow consecutive metabolite incidences')
    backbone = set(spec['items'][::2])
    protected = backbone | set(protected_metabolites)
    nodes = [chain[0]['source']] + [s['target'] for s in chain]
    counts = Counter(nodes)
    # A non-backbone primary outlet is free to move, but must not be cloned:
    # doing so disconnects the complete side pathway from its true producer.
    # Use topology rather than candidate IDs or particular metabolite names.
    shared_entry = nodes[0] if (nodes[0] in protected and nodes[0] not in backbone
                               and counts[nodes[0]] == 1) else None
    overrides = spec.setdefault('view_overrides', {})
    groups = spec.setdefault('shared_instances', {})
    for name, group in list(groups.items()):
        involved = set(group['reactions']) & set(ids)
        if involved:
            if involved != set(group['reactions']):
                raise ValueError('cannot unfold a shared view spanning unrelated reactions')
            groups.pop(name)
            overrides.pop(name, None)
    for name in list(overrides):
        if name in ids or any(name.endswith('_' + rid) for rid in ids):
            overrides.pop(name)

    cx, cy = spec['center']
    x = cx + 1120
    spacing = 240 if len(chain) <= 3 else 180 if len(chain) <= 6 else 160
    direction = 1
    y0 = max(180, cy - len(chain) * spacing / 2)
    if shared_entry:
        # Enter outside the lower edge of the primary ring, then go upward.
        # This leaves the real producing incidence clear of the ring and of
        # the side pathway's cofactor fan. It is not a second circular arc.
        direction = -1
        y0 = max(cy + spec['radius'] + 145, 180 + len(chain) * spacing)
    positions = [[x, y0 + direction * i * spacing] for i in range(len(nodes))]
    font = float(out['publication_style']['metabolite_font_pt']) * 4 / 3

    def text_offset(mid):
        name = out['metabolites'][mid].get('display_name', mid)
        return -((len(name) + 6) * font * .53 / 2 + 22)

    for i, side in enumerate(chain):
        side['start'], side['end'] = positions[i], positions[i+1]
        side.pop('control_points', None)
        side['label_offset'] = [-92, -4]
        for key, mid in (('local_source', side['source']), ('local_target', side['target'])):
            side[key] = mid != shared_entry and (mid in protected or counts[mid] > 1)
        # Discard old polar cofactor coordinates AND their final-view overrides.
        ym = (positions[i][1] + positions[i+1][1]) / 2
        reaction = out['reactions'][side['reaction']]
        cofactors = reaction.get('cofactors', {})
        equation = stoichiometry[side['reaction']]
        for sign in (-1, 1):
            members = [mid for mid in cofactors if equation[mid] * sign > 0]
            for rank, mid in enumerate(members):
                # A vertical fan, with input/output columns separated. Labels
                # face away from the trunk instead of lying across another edge.
                xx = x + (310 if sign < 0 else 610)
                yy = ym + direction * sign * (26 + 36 * rank)
                label = cofactors[mid].get('display_name', mid)
                dx = (len(label) + 2) * font * .53 / 2 + 18
                cofactors[mid].update(x=xx, y=yy, align='center', label_dx=dx, label_dy=-4)

    for i, (mid, pos) in enumerate(zip(nodes, positions)):
        incidences = []
        if i:
            incidences.append(chain[i-1])
        if i < len(chain):
            incidences.append(chain[i])
        local = mid != shared_entry and (mid in protected or counts[mid] > 1)
        off = [text_offset(mid), -4]
        if mid == shared_entry:
            off = [0, 42]
            spec.setdefault('shared_metabolites', {})[mid] = {
                'position': pos, 'label_offset': off,
            }
        if local:
            # Only neighbouring incidences share a view. In particular, the
            # first and last occurrences of a returning pool must stay apart.
            if len(incidences) == 2:
                groups[f'linear_{mid}_{i}'] = {
                    'metabolite': mid, 'reactions': [s['reaction'] for s in incidences],
                    'position': pos, 'label_offset': off,
                }
            for side in incidences:
                overrides[mid + '_' + side['reaction']] = {
                    'x': pos[0], 'y': pos[1], 'label_dx': off[0], 'label_dy': off[1], 'node_radius': 9,
                }
            if mid not in backbone:
                # Its canonical copy belongs to the primary cycle/other branch.
                # Do not carry a stale coordinate from the old second ring.
                prior = overrides.get(mid, {})
                prior.pop('x', None)
                prior.pop('y', None)
                prior.update(label_dx=80, label_dy=22)
                overrides[mid] = prior
        else:
            out['metabolites'][mid].update(x=pos[0], y=pos[1])
            overrides[mid] = {'x': pos[0], 'y': pos[1], 'label_dx': off[0], 'label_dy': off[1], 'node_radius': 9}

    # Additional carbon inputs/outputs stay local to their real reaction trunk.
    # This includes DHA -> DHAPT -> Pyruvate; it is not a second backbone.
    for i, side in enumerate(chain):
        ym = (positions[i][1] + positions[i+1][1]) / 2
        cofactors = out['reactions'][side['reaction']].get('cofactors', {})
        for mid, coefficient in stoichiometry[side['reaction']].items():
            if mid in cofactors or mid in (side['source'], side['target']) or mid in backbone:
                continue
            # Carbon branches use the left side; neighbouring labels use right.
            for trunk_mid in (side['source'], side['target']):
                dx = abs(text_offset(trunk_mid))
                if trunk_mid not in protected and counts[trunk_mid] == 1:
                    overrides[trunk_mid]['label_dx'] = dx
                for sibling in chain:
                    key = trunk_mid + '_' + sibling['reaction']
                    if key in overrides:
                        overrides[key]['label_dx'] = dx
                for group in groups.values():
                    if group['metabolite'] == trunk_mid and set(group['reactions']) <= set(ids):
                        group['label_offset'][0] = dx
            sign = -1 if coefficient < 0 else 1
            pos = [x - 300, ym + sign * 120]
            out['metabolites'][mid].update(x=pos[0], y=pos[1])
            overrides[mid] = {'x': pos[0], 'y': pos[1], 'label_dx': -15, 'label_dy': sign * 40, 'node_radius': 9}
            if mid in spec.get('shared_metabolites', {}):
                spec['shared_metabolites'][mid].update(position=pos, label_offset=[-15, sign*40])
            for other in spec['side_reactions']:
                if other['reaction'] in ids:
                    continue
                for endpoint, position_key, local_key in (('source', 'start', 'local_source'), ('target', 'end', 'local_target')):
                    if other[endpoint] == mid and not other.get(local_key):
                        other[position_key] = pos
                        if endpoint == 'target':
                            # A bridge connected to this moved carbon input must
                            # move its product satellites too, not leave them at
                            # the former circular-layout coordinates.
                            source_mid = other['source']
                            if source_mid not in backbone:
                                overrides.setdefault(source_mid, {}).update(label_dx=-85, label_dy=-10)
                            bridge_rid = other['reaction']
                            bridge_cofactors = out['reactions'][bridge_rid].get('cofactors', {})
                            products = [m for m in bridge_cofactors if stoichiometry[bridge_rid][m] > 0]
                            for rank, product in enumerate(products):
                                placement = {'x': pos[0]+80+rank*100, 'y': pos[1]+120,
                                             'label_dx': 0, 'label_dy': 32}
                                bridge_cofactors[product].update(placement)
                                overrides.setdefault(product+'_'+bridge_rid, {}).update(placement)

    out['canvas'] = [2800, max(1870, math.ceil(max(p[1] for p in positions) + 240))]
    spec['instance_legend_position'] = [out['canvas'][0]/2, out['canvas'][1]-135]
    out['publication_style']['legend_position'] = [out['canvas'][0]/2-210, out['canvas'][1]-100]
    spec['side_pathway_presentation'] = {'mode': 'open_chain', 'reaction_ids': ids,
        'direction': direction, 'shared_entry': shared_entry,
        'returning_pool': nodes[0] if nodes[0] == nodes[-1] else None,
        'same_pool_return_uses_labels': nodes[0] == nodes[-1]}
    return out
