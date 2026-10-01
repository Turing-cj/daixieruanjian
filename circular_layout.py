"""Compile a closed circular visual backbone without changing a biochemical model.

An entity may have several *view instances*, but still has one model ID. Extra
carbon substrates/products use local_instances, never the cofactor namespace.
The browser draws one main arc and one branch per additional endpoint, not
arbitrary substrate/product pairs. Legacy layout remains available unchanged.
"""
from __future__ import annotations

import math
from collections import Counter


def _info(obj):
    return obj.notes.setdefault('map_info', {})


def _positive(value, name):
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be finite and positive')
    return value


def _point(value, name):
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError(f'{name} must be [x, y]')
    result = [float(x) for x in value]
    if not all(math.isfinite(x) for x in result):
        raise ValueError(f'{name} must be finite')
    return result


def apply_pathway_multiplicity(model, flux, policy):
    """Annotate pool throughput, not concentration or per-incidence demand.

    For a canonical metabolite m, P=sum(max(S_mr*v_r, 0)) and
    C=sum(max(-S_mr*v_r, 0)). The displayed amount is max(P,C)/scale.
    P=C for balanced internal pools; max also covers boundary-only pools.
    Both sides and their residual are retained, so this is NOT a balance proof.
    Model IDs, stoichiometry, fluxes and visual geometry are never modified.
    """
    if not policy.get('enabled', False):
        return {'enabled': False}
    mode = policy.get('label_mode', 'pool_total')
    if mode not in ('pool_total', 'local_flow'):
        raise ValueError(f'unknown pathway label_mode: {mode}')
    relationship_mode = policy.get('relationship_mode', 'pool')
    if relationship_mode not in ('pool', 'source_split'):
        raise ValueError(f'unknown pathway relationship_mode: {relationship_mode}')
    if relationship_mode == 'source_split' and mode != 'local_flow':
        raise ValueError('source_split relationships require local_flow quantities')
    scale = _positive(policy.get('normalization_flux', 1), 'multiplicity normalization')
    sides = {m.id: {'produced': [], 'consumed': []} for m in model.metabolites}
    # Validate everything before mutating any presentation metadata.
    for reaction in model.reactions:
        try:
            value = float(flux[reaction.id])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f'missing or invalid multiplicity flux: {reaction.id}') from error
        if not math.isfinite(value):
            raise ValueError(f'nonfinite multiplicity flux: {reaction.id}')
        for metabolite, coefficient in reaction.metabolites.items():
            rate = float(coefficient) * value / scale
            if not math.isfinite(rate):
                raise ValueError(f'nonfinite multiplicity incidence: {reaction.id}/{metabolite.id}')
            sides[metabolite.id]['produced' if rate >= 0 else 'consumed'].append(abs(rate))
    counts, balances = {}, {}
    for mid, values in sides.items():
        produced, consumed = math.fsum(values['produced']), math.fsum(values['consumed'])
        count = max(produced, consumed)
        counts[mid] = count
        balances[mid] = {'produced': produced, 'consumed': consumed, 'net': produced-consumed}

    if mode == 'local_flow':
        from branch_allocation import apply_local_flow_labels
        return apply_local_flow_labels(model, flux, scale, counts, balances, relationship_mode)

    def annotate(info, mid, fallback):
        # Keep a separate base: a second compilation cannot produce "3 3 F6P".
        base = info.get('pool_label_base', info.get('pathway_label_base', info.get('display_name', fallback)))
        info['pathway_label_base'] = base
        info.pop('flow_relation', None)
        info['pathway_multiplicity'] = counts[mid]
        info['display_name'] = f'{counts[mid]:.12g} {base}'

    for metabolite in model.metabolites:
        annotate(_info(metabolite), metabolite.id, metabolite.id)
    for reaction in model.reactions:
        for key in ('local_instances', 'cofactors'):
            for mid, entry in _info(reaction).get(key, {}).items():
                canonical = _info(model.metabolites.get_by_id(mid))
                annotate(entry, mid, canonical['pathway_label_base'])
    result = {'enabled': True, 'normalization_flux': scale,
              'convention': 'max(gross production, gross consumption) / normalization_flux',
              'counts': counts, 'balances': balances,
              'legend': 'Numbers = whole-pathway pool totals, not branch amounts.'}
    presentation = model.notes.get('circular_backbone', {})
    codes = presentation.get('pool_codes', {})
    if codes:
        presentation['instance_legend'] = (' / '.join(f'[{codes[mid]}] {mid}' for mid in sorted(codes))
                                            + ': matching labels share the same pool.')
    model.notes['pathway_multiplicity'] = result
    return result


def compile_circular_layout(model, layout):
    spec = layout['circular_backbone']
    cx, cy = _point(spec.get('center'), 'circle center')
    radius = _positive(spec.get('radius'), 'circle radius')
    start = float(spec.get('start_angle_deg', -90))
    if not math.isfinite(start):
        raise ValueError('circle start angle must be finite')
    items = list(spec.get('items', []))
    if len(items) < 6 or len(items) % 2 or len(set(items)) != len(items):
        raise ValueError('backbone must alternate unique metabolites/reactions in a closed cycle')
    lookup = {o.id:o for o in list(model.metabolites) + list(model.reactions)}
    for index, item in enumerate(items):
        collection = model.metabolites if index % 2 == 0 else model.reactions
        if item not in collection:
            raise ValueError(f'backbone item {index}: unexpected or unknown ID {item}')
    # A circle must follow actual incidences; a decorative circle is not a cycle.
    for index in range(1, len(items), 2):
        reaction = model.reactions.get_by_id(items[index])
        stoich = {m.id:c for m,c in reaction.metabolites.items()}
        source, target = items[index-1], items[(index+1) % len(items)]
        if stoich.get(source, 0) >= 0 or stoich.get(target, 0) <= 0:
            raise ValueError(f'{reaction.id}: backbone source/target is not a substrate/product')

    visual = layout.get('scientific_layout', {})
    palette = visual.get('role_palette', {})
    for role, ids in visual.get('roles', {}).items():
        for item in ids:
            if item not in lookup:
                raise ValueError(f'unknown visual role ID: {item}')
            info = _info(lookup[item])
            info['scientific_role'] = role
            if role in palette:
                info['color'] = palette[role]
    for kind in ('metabolite', 'reaction'):
        for item, style in visual.get(kind + '_styles', {}).items():
            if item in lookup and 'display_name' in style:
                _info(lookup[item])['display_name'] = style['display_name']

    inset = _positive(spec.get('satellite_inset', 130), 'satellite inset')
    spacing = _positive(spec.get('satellite_spacing', 65), 'satellite spacing')
    local_carbon_outset = spec.get('local_carbon_outset')
    if local_carbon_outset is not None:
        local_carbon_outset = _positive(local_carbon_outset, 'local carbon outset')
    currency_outset = _positive(spec.get('currency_satellite_outset', 125), 'currency satellite outset')
    external_outset = _positive(spec.get('external_satellite_outset', 145), 'external satellite outset')
    if inset >= radius:
        raise ValueError('satellite inset must be smaller than circle radius')
    external = set(spec.get('external_metabolites', []))
    canonical_satellites = set(spec.get('canonical_satellites', []))
    step = 360.0 / len(items)
    clone_ids = set()
    represented = set()
    shared = spec.get('shared_metabolites', {})
    for mid, entry in shared.items():
        if mid not in model.metabolites or mid in items or _info(lookup[mid]).get('hidden'):
            raise ValueError(f'shared metabolite must be a visible, non-backbone metabolite: {mid}')
        _point(entry.get('position'), f'shared {mid} position')
        # A shared pool is a real node, never an unconnected display clone.
        for reaction in model.reactions:
            _info(reaction).get('local_instances', {}).pop(mid, None)

    # A canonical satellite is a non-backbone metabolite whose one backbone
    # incidence should use the real metabolite node. Other incidences may still
    # be local source/split views. This closes the representation gap between
    # "clone every incidence" and "merge every incidence into one shared hub".
    backbone_reactions = [model.reactions.get_by_id(items[index]) for index in range(1, len(items), 2)]
    for mid in canonical_satellites:
        if mid not in model.metabolites or mid in items or mid in shared or _info(lookup[mid]).get('hidden'):
            raise ValueError(f'canonical satellite must be a visible, non-backbone, non-shared metabolite: {mid}')
        owners = [reaction.id for reaction in backbone_reactions
                  if any(metabolite.id == mid for metabolite in reaction.metabolites)]
        if len(owners) != 1:
            raise ValueError(f'canonical satellite requires exactly one backbone incidence: {mid}/{owners}')

    def label(info, dx=0, dy=0):
        info.update(align='center', label_dx=dx, label_dy=dy)

    def instance(reaction, metabolite, x, y, cofactor=False):
        info = _info(reaction)
        key = 'cofactors' if cofactor else 'local_instances'
        original = _info(metabolite)
        entry = dict(info.get(key, {}).get(metabolite.id, {}))
        entry.update(x=x, y=y, display_name=entry.get('display_name', original.get('display_name', metabolite.id)))
        if 'node_radius' in original:
            entry.setdefault('node_radius', original['node_radius'])
        entry.pop('hidden', None)
        if not cofactor:
            entry.update(display_instance_of=metabolite.id, clone_marker=True,
                         scientific_role=original.get('scientific_role', 'regeneration'))
            clone_ids.add(metabolite.id)
        label(entry)
        info.setdefault(key, {})[metabolite.id] = entry
        return entry

    for index, item in enumerate(items):
        theta = math.radians(start + step * index)
        info = _info(lookup[item])
        info.update(x=cx+radius*math.cos(theta), y=cy+radius*math.sin(theta),
                    layout_structure='circular_backbone', layout_theta_deg=start+step*index,
                    layout_structure_center=[cx,cy], layout_structure_radius=[radius,radius])
        if info.get('hidden'):
            raise ValueError(f'backbone cannot contain a hidden item: {item}')
        offset = 43 if index % 2 == 0 else -48
        label(info, offset*math.cos(theta), offset*math.sin(theta)-4)

    for index in range(1, len(items), 2):
        reaction = model.reactions.get_by_id(items[index])
        info = _info(reaction)
        represented.add(reaction.id)
        source, target = items[index-1], items[(index+1) % len(items)]
        for key in list(info):
            if key.startswith('route_') or key in ('edge_routes', 'reaction_bus', 'edge_route_overrides'):
                del info[key]
        info['route_style'] = 'circle_hyperedge'
        info['circle_arc'] = {
            'center':[cx,cy], 'radius':radius, 'source':source, 'target':target,
            'start_angle_deg':start+step*(index-1),
            'end_angle_deg':start+step*(index+1),
            'reaction_angle_deg':start+step*index, 'join_half_length':26,
        }
        theta = math.radians(start + step * index)
        normal = (math.cos(theta), math.sin(theta))
        tangent = (-normal[1], normal[0])
        by_side = {-1:[], 1:[]}
        for metabolite, coefficient in reaction.metabolites.items():
            if metabolite.id not in (source, target):
                by_side[-1 if coefficient < 0 else 1].append(metabolite)
        for sign, metabolites in by_side.items():
            for rank, metabolite in enumerate(metabolites):
                mid = metabolite.id
                if mid in shared:
                    continue
                cofactor = mid in info.get('cofactors', {}) or _info(metabolite).get('hidden', False)
                outward = cofactor or mid in external or local_carbon_outset is not None
                if mid in external:
                    radial = radius + external_outset
                elif cofactor:
                    radial = radius + currency_outset
                elif local_carbon_outset is not None:
                    radial = radius + local_carbon_outset
                else:
                    radial = radius-inset
                along = 0 if mid in external else sign*(spacing+rank*spacing)
                # Multiple co-products fan out, rather than sharing a label position.
                x = cx + radial*normal[0] + along*tangent[0]
                y = cy + radial*normal[1] + along*tangent[1]
                if mid in external or mid in canonical_satellites:
                    entry = _info(metabolite)
                    entry.update(x=x, y=y)
                else:
                    entry = instance(reaction, metabolite, x, y, cofactor)
                label_offset = 30 if outward else -27
                label(entry, label_offset*normal[0], label_offset*normal[1]-3)

    for mid, entry in shared.items():
        x, y = _point(entry['position'], f'shared {mid} position')
        info = _info(lookup[mid])
        info.update(x=x, y=y)
        label(info, *_point(entry.get('label_offset', [0, -26]), f'shared {mid} label offset'))

    for side in spec.get('side_reactions', []):
        rid = side['reaction']
        if rid in represented or rid not in model.reactions:
            raise ValueError(f'duplicate or unknown side reaction: {rid}')
        reaction = model.reactions.get_by_id(rid)
        source, target = side['source'], side['target']
        stoich = {m.id:c for m,c in reaction.metabolites.items()}
        if stoich.get(source,0) >= 0 or stoich.get(target,0) <= 0:
            raise ValueError(f'{rid}: local side reaction requires a substrate-to-product trunk')
        def endpoint(mid, key, local_key):
            current = _info(lookup[mid])
            point = _point(side.get(key, [current.get('x'), current.get('y')]), f'side {key}')
            if not side.get(local_key) and mid in items and math.dist(point, [current['x'], current['y']]) > 1e-7:
                raise ValueError(f'{rid}: side reaction cannot displace backbone endpoint {mid}')
            return point
        a, b = endpoint(source, 'start', 'local_source'), endpoint(target, 'end', 'local_target')
        if math.dist(a,b) < 80:
            raise ValueError(f'{rid}: local side reaction is too short for readable labels')
        info = _info(reaction)
        for key in list(info):
            if key.startswith('route_') or key in ('edge_routes','reaction_bus','edge_route_overrides'):
                del info[key]
        info.update(x=(a[0]+b[0])/2, y=(a[1]+b[1])/2,
                    route_style='local_hyperedge', local_trunk={'source':source,'target':target})
        if 'control_points' in side:
            controls = side['control_points']
            if not isinstance(controls, (list, tuple)) or len(controls) != 2:
                raise ValueError(f'{rid}: control_points must contain two points')
            controls = [_point(p, f'{rid} control point') for p in controls]
            info['local_trunk']['control_points'] = controls
            info.update(x=(a[0]+3*controls[0][0]+3*controls[1][0]+b[0])/8,
                        y=(a[1]+3*controls[0][1]+3*controls[1][1]+b[1])/8)
        label(info, *_point(side.get('label_offset', [0, -35]), f'{rid} label offset'))
        for mid, point, key in ((source,a,'local_source'), (target,b,'local_target')):
            if side.get(key):
                entry = instance(reaction, model.metabolites.get_by_id(mid), *point)
            else:
                entry = _info(model.metabolites.get_by_id(mid))
                entry.update(x=point[0], y=point[1])
            if side.get(key) or mid not in items and mid not in shared:
                label(entry, 0, 28)
        local_satellites = side.get('local_satellites', {})
        if not isinstance(local_satellites, dict):
            raise ValueError(f'{rid}: local_satellites must be a mapping')
        for mid, satellite in local_satellites.items():
            if (mid in (source, target) or mid not in stoich or mid in shared
                    or mid not in model.metabolites or _info(model.metabolites.get_by_id(mid)).get('hidden')
                    or mid in info.get('cofactors', {})):
                raise ValueError(f'{rid}: invalid local satellite: {mid}')
            if not isinstance(satellite, dict) or 'position' not in satellite:
                raise ValueError(f'{rid}: local satellite requires a position: {mid}')
            x, y = _point(satellite['position'], f'{rid} local satellite {mid}')
            entry = instance(reaction, model.metabolites.get_by_id(mid), x, y)
            label(entry, *_point(satellite.get('label_offset', [0, 28]),
                                 f'{rid} local satellite {mid} label offset'))
        represented.add(rid)

    shared_instances = spec.get('shared_instances', {})
    assigned = set()
    for view_id, group in shared_instances.items():
        mid, reaction_ids = group.get('metabolite'), group.get('reactions', [])
        if not view_id or view_id in lookup or mid not in clone_ids or len(set(reaction_ids)) != len(reaction_ids) or len(reaction_ids) < 2:
            raise ValueError(f'invalid shared instance: {view_id}')
        x, y = _point(group.get('position'), f'shared instance {view_id} position')
        dx, dy = _point(group.get('label_offset', [0, -26]), f'shared instance {view_id} label offset')
        for rid in reaction_ids:
            if rid not in model.reactions or (rid, mid) in assigned:
                raise ValueError(f'invalid or duplicate shared instance incidence: {view_id}/{rid}')
            entry = _info(model.reactions.get_by_id(rid)).get('local_instances', {}).get(mid)
            if not entry or entry.get('display_instance_of') != mid:
                raise ValueError(f'shared instance must reuse an existing carbon incidence: {view_id}/{rid}')
            entry.update(x=x, y=y, view_instance_id=view_id)
            label(entry, dx, dy)
            assigned.add((rid, mid))

    expected = {r.id for r in model.reactions if not _info(r).get('hidden')}
    if represented != expected:
        raise ValueError('circular layout must represent every visible reaction: ' + str(sorted(expected-represented)))
    pool_codes = {mid: chr(65+i) if i < 26 else f'P{i+1}' for i, mid in enumerate(sorted(clone_ids))}
    use_codes = spec.get('instance_marker_style') == 'pool_code'
    for mid in clone_ids:
        info = _info(model.metabolites.get_by_id(mid))
        info['clone_marker'] = True
        base = info.setdefault('instance_base_name', info.get('display_name', mid).rstrip('*'))
        info['display_name'] = base + (f' [{pool_codes[mid]}]' if use_codes else '*')
    for reaction in model.reactions:
        for mid, entry in _info(reaction).get('local_instances', {}).items():
            entry['display_name'] = _info(model.metabolites.get_by_id(mid))['display_name']

    # Optional final view placement for batch layouts. Chemistry is never an
    # override target; primary-ring coordinates still face the radial QA gate.
    views = {obj.id: _info(obj) for obj in lookup.values()}
    for reaction in model.reactions:
        for kind in ('cofactors', 'local_instances'):
            for mid, entry in _info(reaction).get(kind, {}).items():
                views[mid + '_' + reaction.id] = entry
                if entry.get('view_instance_id'):
                    views[entry['view_instance_id']] = entry
    allowed = {'x', 'y', 'label_dx', 'label_dy', 'node_radius', 'align'}
    for view_id, overrides in spec.get('view_overrides', {}).items():
        if view_id not in views or not set(overrides) <= allowed:
            raise ValueError(f'invalid circular view override: {view_id}')
        for key, value in overrides.items():
            if key != 'align' and not math.isfinite(float(value)):
                raise ValueError(f'nonfinite view override: {view_id}/{key}')
        if 'node_radius' in overrides and float(overrides['node_radius']) <= 0:
            raise ValueError(f'nonpositive view radius: {view_id}')
        views[view_id].update(overrides)

    summary = {
        'enabled':True, 'center':[cx,cy], 'radius':radius, 'sweep_deg':360,
        'items':items, 'reaction_count':len(represented),
        'backbone_reaction_count':len(items)//2,
        'clone_ids':sorted(clone_ids),
        'local_instance_count':sum(len(_info(r).get('local_instances',{})) for r in model.reactions),
        'visible_local_instance_count':len({entry.get('view_instance_id', mid+'_'+r.id)
            for r in model.reactions for mid, entry in _info(r).get('local_instances', {}).items()}),
        'title':spec.get('title',''), 'subtitle':spec.get('subtitle',''),
        'title_position':_point(spec.get('title_position',[cx,cy-70]), 'title position'),
        'instance_legend_position':_point(spec.get('instance_legend_position', [layout['canvas'][0]/2, layout['canvas'][1]-65]), 'legend position'),
        'input_arrows':bool(spec.get('input_arrows', False)),
        'shared_metabolites':sorted(shared),
        'canonical_satellites':sorted(canonical_satellites),
        'pool_codes':pool_codes if use_codes else {},
        'shared_instances':shared_instances,
        'instance_legend':(' / '.join(f'[{pool_codes[mid]}] {mid}' for mid in sorted(clone_ids)) + ': matching labels share the same pool.') if use_codes else
                          '* Repeated labels denote the same metabolite pool; no extra reaction is implied.'
    }
    model.notes['circular_backbone'] = summary
    result = {'enabled':True, 'circular_backbone':summary,
              'roles':dict(Counter(_info(o).get('scientific_role','other') for o in lookup.values())),
              'routes':{'mode':'one_trunk_per_reaction', 'pairwise_routes':False}}
    model.notes['scientific_layout'] = result
    return result


def circular_publication_qa(model, layout):
    spec = layout['circular_backbone']
    cx, cy = _point(spec['center'], 'center')
    radius = _positive(spec['radius'], 'radius')
    width, height = layout['canvas']
    lookup = {o.id:o for o in list(model.metabolites)+list(model.reactions)}
    errors, invalid, outside = [], [], []
    for item in spec['items']:
        info = _info(lookup[item])
        errors.append(abs(math.hypot(info['x']-cx, info['y']-cy)-radius))
    for obj in lookup.values():
        info = _info(obj)
        entries = [] if info.get('hidden') else [(obj.id,info)]
        for key in ('cofactors','local_instances'):
            entries.extend((obj.id+':'+mid, entry) for mid,entry in info.get(key,{}).items())
        for name, entry in entries:
            x, y = entry.get('x',float('nan')), entry.get('y',float('nan'))
            if not (math.isfinite(x) and math.isfinite(y)):
                invalid.append(name)
            elif not (0 < x < width and 0 < y < height):
                outside.append(name)
    max_error = max(errors,default=float('inf'))
    return {'pass':max_error < 0.05 and not invalid and not outside,
            'mode':'circular_backbone', 'axis_ratio':1.0,
            'maximum_node_radial_error_px':max_error,
            'missing_coordinates':invalid, 'out_of_bounds':outside,
            'browser_geometry_verification_required':True}
