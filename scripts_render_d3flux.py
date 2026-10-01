#!/usr/bin/env python3
"""Render solved FreeFlux fluxes with the native d3flux example renderer."""
import argparse
import copy
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

from cobra import Model, Reaction, Metabolite
from d3flux import flux_map
from d3flux.core.flux_layouts import create_model_json

from scripts_standalone_d3flux import write_standalone, export_native
from scientific_layout import apply_scientific_layout, build_publication_qa
from circular_layout import apply_pathway_multiplicity

TOKEN = re.compile(r'(?:(\d+(?:\.\d+)?))?([^\W\d][\w-]*)', re.UNICODE)

PATHWAY_COLORS = {
    'CBB': '#93c8ed', 'PPP': '#f7c1c1', 'EMP': '#f9c780',
    'ED': '#aedbd2', 'TCA': '#bc9cc8', 'Anaplerotic': '#d5e9c6',
    'Exchange': '#e4f0f7', 'Other': '#999999',
}
PATHWAY_COLORS_STANDARD2 = {
    'CBB': '#aea7d1', 'PPP': '#dbe5f5', 'EMP': '#f2bcd1',
    'ED': '#fadad4', 'TCA': '#50b9b8', 'Anaplerotic': '#474a9d',
    'Exchange': '#8ec7d7', 'Other': '#e3c0da',
}
METABOLITE_COLORS = {
    'ATP': '#ffbb78', 'CO2': '#aec7e8', 'NADH': '#b5cf6b',
    'Q8': '#b5cf6b', 'Q8H2': '#b5cf6b', 'NADPH': '#8ca252',
}
PATHWAY_REACTIONS = {
    'CBB': {'v2_Rubisco', 'v12_prk'},
    'PPP': {'v3_sba', 'v4_sbp', 'v5_tkt1', 'v6_tkt2', 'v7_tal',
            'v10_rpi', 'v11_rpe', 'v14_g6pdh'},
    'EMP': {'v13_pgi', 'v15_pfk', 'v16_fba', 'v17_pgk', 'v18_tpi',
            'v19_gapdh', 'v20_gpm', 'v21_eno', 'v22_pk'},
    'TCA': {'v23_pdh', 'v24_cs', 'v25_can', 'v26_icd', 'v27_KGDH',
            'v27_AKGDH', 'v28_sdh', 'v29_fum', 'v30_mdh'},
    'Anaplerotic': {'v31_me', 'v32_ppc'},
    'Exchange': {'v1_fdh', 'v33_CO2_out'},
}

METABOLITE_CATEGORY_COLORS = {
    'Carbon source / exchange': '#4c78a8',
    'CBB intermediates': '#2f8f5b',
    'EMP / sugar phosphate': '#f28e2b',
    'PPP intermediates': '#b279a2',
    'TCA / organic acid': '#e15759',
    'Biomass / product': '#59a14f',
    'Dilution / dummy': '#9c9c9c',
    'Other': '#7f7f7f',
}
METABOLITE_CATEGORIES = {
    'FOR': 'Carbon source / exchange', 'CO2': 'Carbon source / exchange',
    'CO2_out': 'Carbon source / exchange',
    'RuBP': 'CBB intermediates', 'SBP': 'CBB intermediates',
    'BPG': 'EMP / sugar phosphate', 'PG3': 'EMP / sugar phosphate', 'GAP': 'EMP / sugar phosphate',
    'DHAP': 'EMP / sugar phosphate', 'F6P': 'EMP / sugar phosphate',
    'FBP': 'EMP / sugar phosphate', 'G6P': 'EMP / sugar phosphate',
    'PG2': 'EMP / sugar phosphate', 'PEP': 'EMP / sugar phosphate',
    'PYR': 'EMP / sugar phosphate',
    'E4P': 'PPP intermediates', 'S7P': 'PPP intermediates',
    'R5P': 'PPP intermediates', 'Ru5P': 'PPP intermediates',
    'Xu5P': 'PPP intermediates',
    'AcCoA': 'TCA / organic acid', 'OAA': 'TCA / organic acid',
    'Cit': 'TCA / organic acid', 'ICit': 'TCA / organic acid',
    'AKG': 'TCA / organic acid', 'Suc': 'TCA / organic acid',
    'Fum': 'TCA / organic acid', 'Mal': 'TCA / organic acid',
    'Biomass': 'Biomass / product',
    'Cits': 'Dilution / dummy', 'Citu': 'Dilution / dummy',
    'GAPs': 'Dilution / dummy', 'GAPu': 'Dilution / dummy',
    'sink': 'Dilution / dummy',
}


def pathway_for_reaction(reaction_id):
    return next((name for name, ids in PATHWAY_REACTIONS.items() if reaction_id in ids), 'Other')


def parse_side(text, sign):
    out = []
    for part in text.split('+'):
        part = re.sub(r'\([^)]*\)$', '', part.strip())
        if not part:
            continue
        match = TOKEN.fullmatch(part)
        if not match:
            raise ValueError(f'cannot parse metabolite token: {part}')
        out.append((match.group(2), sign * float(match.group(1) or 1.0)))
    return out


def read_model(path):
    model = Model('FreeFlux_V5_3')
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        rid, substrates, products, reversible = line.split('\t')[:4]
        if rid in model.reactions:
            raise ValueError(f'duplicate reaction: {rid}')
        stoich = defaultdict(float)
        for mid, coefficient in parse_side(substrates, -1) + parse_side(products, 1):
            stoich[mid] += coefficient
            if mid not in model.metabolites:
                model.add_metabolites([Metabolite(mid, compartment='c')])
        reaction = Reaction(rid, lower_bound=-1000 if int(reversible) else 0)
        reaction.add_metabolites({model.metabolites.get_by_id(mid): value
                                 for mid, value in stoich.items() if value})
        reaction.notes['freeflux_atom_equation'] = substrates + ' -> ' + products
        model.add_reactions([reaction])
    return model


def apply_layout(model, layout):
    template_path = Path(layout['template'])
    if not template_path.is_absolute():
        template_path = Path(__file__).resolve().parent / template_path
    template = json.loads(template_path.read_text(encoding='utf-8'))
    scale, dx, dy = layout['template_transform']
    # FreeFlux execution sanitizes the alpha character in this reaction ID.
    aliases = {'v27_KGDH': 'v27_AKGDH'}
    for kind in ('metabolites', 'reactions'):
        originals = {item['id']: item for item in template[kind]}
        for obj in getattr(model, kind):
            info = {}
            original = layout['template_' + kind].get(obj.id)
            if original is None:
                original = layout['template_' + kind].get(aliases.get(obj.id, ''))
            if original:
                source = originals[original]['notes']['map_info']
                info = {'x': source['x'] * scale + dx, 'y': source['y'] * scale + dy}
                if 'align' in source:
                    info['align'] = source['align']
            info.update(layout[kind].get(obj.id, layout[kind].get(aliases.get(obj.id, ''), {})))
            if obj.id in layout.get('hidden_' + kind, []):
                info['hidden'] = True
            if not info.get('hidden') and not {'x', 'y'} <= info.keys():
                raise ValueError(f'missing layout coordinates: {obj.id}')
            if kind == 'metabolites':
                info.setdefault('category', METABOLITE_CATEGORIES.get(obj.id, 'Other'))
                metabolite_id = re.sub(r'_[a-z]$', '', obj.id, flags=re.IGNORECASE).upper()
                if metabolite_id in METABOLITE_COLORS:
                    info['color'] = METABOLITE_COLORS[metabolite_id]
                else:
                    info.pop('color', None)
                info.setdefault('display_name', obj.id)
            else:
                info.setdefault('display_name', obj.id)
                info.setdefault('pathway', pathway_for_reaction(obj.id))
                info.setdefault('color', PATHWAY_COLORS[info['pathway']])
                for mid in info.get('cofactors', {}):
                    if mid not in {met.id for met in obj.metabolites}:
                        raise ValueError(f'{obj.id}: layout cofactor {mid} not in this reaction')
            obj.notes['map_info'] = info
    # Keep palette metadata outside map_info: d3flux forwards map_info as
    # render_model kwargs, so arbitrary keys there become invalid arguments.
    model.notes['metabolite_category_colors'] = dict(METABOLITE_CATEGORY_COLORS)
    model.notes['palette_schemes'] = {
        'metabolites': {
            'Standard': {'mode': 'standard'},
            'Categories': {'mode': 'category', 'colors': dict(METABOLITE_CATEGORY_COLORS)},
        },
        'reactions': {
            'Standard': dict(PATHWAY_COLORS),
            'Standard2': dict(PATHWAY_COLORS_STANDARD2),
        },
    }


def apply_line_style(model, flux, policy):
    mode = policy.get('mode', 'relative_dash')
    if mode not in {'relative_dash', 'flux_width_only'}:
        raise ValueError(f'unsupported line-style mode: {mode}')
    threshold = float(policy['threshold'])
    if not math.isfinite(threshold) or not 0 < threshold <= 1:
        raise ValueError('relative dash threshold must be in (0, 1]')
    visible = [r for r in model.reactions if not r.notes['map_info'].get('hidden')
               and any(not m.notes['map_info'].get('hidden') for m in r.reactants)
               and any(not m.notes['map_info'].get('hidden') for m in r.products)]
    reference = policy.get('reference_flux')
    automatic = reference is None
    if automatic:
        reference = max((abs(flux[r.id]) for r in visible), default=0.0)
    reference = float(reference)
    if not math.isfinite(reference) or reference < 0 or (not automatic and reference == 0):
        raise ValueError('fixed reference_flux must be finite and positive')
    dashed = []
    width_min, width_max = map(float, policy.get('width_range', [1.6, 4.0]))
    width_exponent = float(policy.get('width_exponent', 0.65))
    if not (0 < width_min <= width_max and width_exponent > 0):
        raise ValueError('line style width_range/exponent must be positive and ordered')
    sorted_abs = sorted((abs(flux[r.id]) for r in visible), reverse=True)
    for reaction in visible:
        info = reaction.notes['map_info']
        relative = abs(flux[reaction.id]) / reference if reference else 0.0
        info['relative_flux'] = relative
        info['stroke'] = width_min + (width_max - width_min) * (min(1.0, relative) ** width_exponent)
        # Publication mode encodes magnitude only with width. Dashes are
        # reserved for semantic uncertainty / knockouts, not low flux.
        info['dasharray'] = ('5,5' if mode == 'relative_dash' and relative < threshold else 'none')
        info['relative_rank'] = next(i for i, value in enumerate(sorted_abs, 1) if value == abs(flux[reaction.id])) / len(visible)
        if info.get('group') == 'ko':
            info['dasharray'] = '5,5,1,5'
        elif mode == 'relative_dash' and relative < threshold:
            dashed.append(reaction.id)
    zero_rule = ('minimum-width solid unless a semantic style explicitly requests dashes'
                 if mode == 'flux_width_only'
                 else 'all visible zero-flux reactions are dashed')
    result = {'mode': mode, 'threshold': threshold, 'reference_flux': reference,
              'reference_scope': 'visible_reactions' if automatic else 'fixed',
              'absolute_threshold': reference * threshold,
              'width_range': [width_min, width_max], 'width_exponent': width_exponent,
              'unit': policy.get('unit', 'mmol/gDCW/h'), 'low_relative_flux_reactions': dashed,
              'zero_reference_rule': zero_rule}
    model.notes['flux_map_line_style'] = result
    return result


def apply_concentration_sizes(model, concentrations):
    """Store log-normalized pool concentrations for radius rendering."""
    values = {m.id: float(concentrations[m.id]) for m in model.metabolites
              if m.id in concentrations and math.isfinite(float(concentrations[m.id]))
              and float(concentrations[m.id]) > 0}
    if not values:
        return {'mode': 'unavailable', 'unit': 'mmol/gDCW'}
    logs = {mid: math.log10(value) for mid, value in values.items()}
    lo, hi = min(logs.values()), max(logs.values())
    span = hi - lo
    for metabolite in model.metabolites:
        info = metabolite.notes['map_info']
        if metabolite.id not in values:
            continue
        info['concentration'] = values[metabolite.id]
        info['concentration_log10'] = logs[metabolite.id]
        info['concentration_normalized'] = (logs[metabolite.id] - lo) / span if span else 0.5
    result = {'mode': 'log10_minmax', 'unit': 'mmol/gDCW',
              'minimum': 10 ** lo, 'maximum': 10 ** hi,
              'node_radius_range': [4, 12], 'missing_value_rule': 'fallback_to_flux'}
    model.notes['flux_map_concentration_size'] = result
    return result


def apply_initial_map(model, path, inherit_coordinates=False):
    """Reuse saved visual style while keeping the current layout authoritative.

    Coordinates and label offsets are geometry, not style.  They are inherited
    only when explicitly requested; otherwise the current layout/scientific
    layout remains the single source of truth.
    """
    saved = json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(saved.get('metabolites'), list) or not isinstance(saved.get('reactions'), list):
        raise ValueError('initial map must be a saved D3Flux Map JSON')
    dynamic = {'flux', 'concentration', 'concentration_log10',
               'concentration_normalized', 'relative_flux', 'relative_rank',
               'dasharray', 'stroke'}
    geometry = {'x', 'y', 'label_dx', 'label_dy'}
    passive_style = {'align', 'display_name', 'color', 'category'}
    for kind in ('metabolites', 'reactions'):
        saved_by_id = {item.get('id'): item for item in saved[kind] if isinstance(item, dict)}
        for obj in getattr(model, kind):
            source = saved_by_id.get(obj.id, {}).get('notes', {}).get('map_info', {})
            current = obj.notes['map_info']
            if inherit_coordinates:
                carried = {key: value for key, value in source.items() if key not in dynamic}
                current.update(copy.deepcopy(carried))
            else:
                # The current layout owns geometry, visibility, cofactor
                # expansion, routing, pathway semantics, and explicit styles.
                # A saved map may only fill missing passive presentation keys.
                for key in passive_style:
                    if key in source and key not in current:
                        current[key] = copy.deepcopy(source[key])
    saved_notes = saved.get('notes', {})
    for key in ('metabolite_category_colors', 'palette_schemes'):
        if key in saved_notes:
            model.notes[key] = copy.deepcopy(saved_notes[key])
    saved_map_info = saved_notes.get('map_info', {})
    for key in ('pathway_colors', 'guide_circles'):
        if key in saved_map_info:
            model.notes.setdefault('map_info', {})[key] = copy.deepcopy(saved_map_info[key])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fit_result', type=Path)
    parser.add_argument('model', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--layout', type=Path, default=Path(__file__).with_name('v5_3_d3flux_layout.json'))
    parser.add_argument('--initial-map', type=Path, help='Saved Map JSON used as the initial visual template')
    parser.add_argument('--inherit-map-coordinates', action='store_true',
                        help='Explicitly let --initial-map override layout x/y and label offsets')
    args = parser.parse_args()
    fit = json.loads(args.fit_result.read_text(encoding='utf-8'))
    flux = fit.get('opt_net_fluxes') or {}
    sources = [args.fit_result, args.model, args.layout]
    if fit.get('flux_csv'):
        table = Path(fit['flux_csv'])
        if not table.is_absolute():
            candidates = (args.fit_result.parent / table, Path.cwd() / table)
            table = next((candidate for candidate in candidates if candidate.exists()), candidates[0])
        with table.open(encoding='utf-8-sig', newline='') as handle:
            flux = {row['reaction_id']: float(row['net_flux']) for row in csv.DictReader(handle)}
        sources.append(table)
    model = read_model(args.model)
    for reaction in model.reactions:
        if reaction.id not in flux or not math.isfinite(float(flux[reaction.id])):
            raise ValueError(f'missing or nonfinite solved flux: {reaction.id}')
    flux = {reaction.id: float(flux[reaction.id]) for reaction in model.reactions}
    layout = json.loads(args.layout.read_text(encoding='utf-8'))
    template_path = Path(layout['template'])
    if not template_path.is_absolute():
        template_path = Path(__file__).resolve().parent / template_path
    sources.append(template_path)
    apply_layout(model, layout)
    if args.initial_map:
        apply_initial_map(model, args.initial_map, inherit_coordinates=args.inherit_map_coordinates)
        sources.append(args.initial_map)
    scientific_layout = apply_scientific_layout(model, layout)
    multiplicity = apply_pathway_multiplicity(model, flux, layout.get('pathway_multiplicity', {}))
    concentration_size = apply_concentration_sizes(model, fit.get('opt_concentrations') or {})
    line_policy = layout.get('line_style', layout.get('relative_dash', {'threshold': 0.4}))
    line_style = apply_line_style(model, flux, line_policy)
    publication_qa = build_publication_qa(model, layout)
    html = flux_map(model, flux_dict=flux, figsize=layout['canvas'],
                    fontsize=14, hide_unused=False, figure_id='freeflux_map',
                    display_name_format=lambda met: met.id)
    if line_style['mode'] == 'flux_width_only':
        legend = ["Flux magnitude encoded by line width",
                  "Dashed lines reserved for semantic uncertainty / knockouts",
                  f"Reference: {line_style['reference_flux']:.6g} {line_style['unit']}"]
    else:
        legend = [f"Low relative flux (< {line_style['threshold'] * 100:g}% of reference)",
                  f"Relative flux >= {line_style['threshold'] * 100:g}% of reference",
                  f"Reference: {line_style['reference_flux']:.6g} {line_style['unit']}"]
    publication_style = layout.get('publication_style', {})
    default_met_font = int(publication_style.get('metabolite_font_pt', 15))
    default_rxn_font = int(publication_style.get('reaction_font_pt', 9))
    legend_x, legend_y = publication_style.get('legend_position', [1120, 910])
    legend_js = '''
    var controlRoot = document.getElementById('freeflux_dynamic_controls');
    function mountControl(element) {
      element.className = (element.className ? element.className + ' ' : '') + 'control-card';
      controlRoot.appendChild(element);
      return element;
    }
    var standardPathwayColors = PATHWAY_COLORS_PLACEHOLDER;
    var standardMetaboliteColors = METABOLITE_COLORS_PLACEHOLDER;
    var pathwayColors = Object.assign({}, standardPathwayColors);
    freeflux_mapmodel.notes = freeflux_mapmodel.notes || {};
    freeflux_mapmodel.notes.map_info = freeflux_mapmodel.notes.map_info || {};
    var savedPathwayColors = freeflux_mapmodel.notes.map_info.pathway_colors || {};
    Object.keys(savedPathwayColors).forEach(function(name) { pathwayColors[name] = savedPathwayColors[name]; });
    function persistPathwayColors() { freeflux_mapmodel.notes.map_info.pathway_colors = pathwayColors; }
    persistPathwayColors();
    var pathwayControls = document.createElement('div'); pathwayControls.id='pathway-color-controls';
    pathwayControls.innerHTML='<strong>Pathway colors</strong><p class="control-hint">Hover or focus a color to highlight its pathway.</p><div class="pathway-grid"></div>';
    Object.keys(pathwayColors).forEach(function(name) {
      var label=document.createElement('label'), title=document.createElement('span'), input=document.createElement('input');
      label.className='pathway-chip'; title.textContent=name;
      input.type='color'; input.dataset.pathway=name; input.value=pathwayColors[name];
      input.setAttribute('aria-label', name+' pathway color');
      label.appendChild(title); label.appendChild(input);
      function highlight(active) {
        svg.selectAll('path.link').classed('pathway-muted', function(d) { return active && d.rxn.notes.map_info.pathway!==name; })
          .classed('pathway-highlight', function(d) { return active && d.rxn.notes.map_info.pathway===name; });
      }
      label.onmouseenter=function() { highlight(true); };
      label.onmouseleave=function() { highlight(label.contains(document.activeElement)); };
      input.onfocus=function() { highlight(true); };
      input.onblur=function() { highlight(false); };
      pathwayControls.querySelector('.pathway-grid').appendChild(label);
    });
    var pathwaySection = document.createElement('section'); pathwaySection.id='pathway-section';
    pathwaySection.innerHTML='<strong>Pathways &amp; reactions</strong>';
    mountControl(pathwaySection);
    pathwaySection.appendChild(pathwayControls);
    function updatePathwayColors() {
      pathwayControls.querySelectorAll('input[data-pathway]').forEach(function(input) {
        pathwayColors[input.getAttribute('data-pathway')] = input.value;
      });
      persistPathwayColors();
      svg.selectAll('path.link').each(function(d) {
        var info=d.rxn.notes.map_info, color=pathwayColors[info.pathway] || info.color || '#999999';
        d3.select(this).attr('stroke', color);
      });
      svg.selectAll('marker').select('path').each(function() {
        var reaction=d3.select(this.parentNode).datum(), pathway=reaction && reaction.notes && reaction.notes.map_info.pathway;
        d3.select(this).attr('fill', pathwayColors[pathway] || '#999999');
      });
    }
    pathwayControls.querySelectorAll('input[data-pathway]').forEach(function(input) { input.oninput=updatePathwayColors; });

    function applyStandardPalette() {
      Object.keys(pathwayColors).forEach(function(name) { delete pathwayColors[name]; });
      Object.keys(standardPathwayColors).forEach(function(name) { pathwayColors[name]=standardPathwayColors[name]; });
      pathwayControls.querySelectorAll('input[data-pathway]').forEach(function(input) {
        input.value=pathwayColors[input.dataset.pathway] || '#999999';
      });
      freeflux_mapmodel.reactions.forEach(function(reaction) {
        var info=reaction.notes.map_info, color=pathwayColors[info.pathway] || pathwayColors.Other;
        info.color=color;
      });
      freeflux_mapmodel.metabolites.forEach(function(metabolite) {
        var info=metabolite.notes.map_info;
        var id=(info.source_metabolite_id || metabolite.id).toUpperCase().replace(/_[A-Z]$/, '');
        if (standardMetaboliteColors[id]) info.color=standardMetaboliteColors[id];
        else delete info.color;
      });
      updatePathwayColors();
      refreshStandardMetaboliteLegend();
      metaboliteCategoryLegend.style.display='grid';
      svg.selectAll('circle.metabolite').style('fill', function(d) {
        var source=categoryMetabolites[d.notes.orig_id || d.id];
        return source && source.notes.map_info.color || '#1f77b4';
      });
    }

    var metaboliteCategoryColors = Object.assign(Object.create(null),
      freeflux_mapmodel.notes.metabolite_category_colors || METABOLITE_CATEGORY_COLORS_PLACEHOLDER);
    metaboliteCategoryColors.Other = metaboliteCategoryColors.Other || '#7f7f7f';
    var categoryMetabolites = Object.create(null);
    freeflux_mapmodel.metabolites.forEach(function(m) {
      m.notes=m.notes||{}; m.notes.map_info=m.notes.map_info||{};
      var info=m.notes.map_info;
      info.category=info.category || 'Other';
      metaboliteCategoryColors[info.category]=metaboliteCategoryColors[info.category] || info.color || '#7f7f7f';
      categoryMetabolites[m.id]=m;
    });
    function persistMetaboliteCategoryColors() { freeflux_mapmodel.notes.metabolite_category_colors = metaboliteCategoryColors; }
    function repaintMetaboliteCategories() {
      freeflux_mapmodel.metabolites.forEach(function(m) {
        var info=m.notes.map_info; info.color=metaboliteCategoryColors[info.category] || metaboliteCategoryColors.Other;
      });
      svg.selectAll('circle.metabolite').style('fill', function(d) {
        // d3flux renders copies, including local copies of cofactor metabolites.
        var source=categoryMetabolites[d.notes.orig_id || d.id];
        var info=source ? source.notes.map_info : {category:'Other', color:metaboliteCategoryColors.Other};
        d.notes.map_info.category=info.category; d.notes.map_info.color=info.color;
        return info.color;
      }).attr('data-category', function(d) { return d.notes.map_info.category; });
      svg.selectAll('circle.metabolite').selectAll('title.category-title').data(function(d) {
        return [d.id + ': ' + d.notes.map_info.category];
      }).enter().append('title').attr('class', 'category-title');
      svg.selectAll('circle.metabolite title.category-title').text(function() {
        var d=d3.select(this.parentNode).datum(); return d.id + ': ' + d.notes.map_info.category;
      });
    }
    var metaboliteCategoryControls = document.createElement('div'); metaboliteCategoryControls.id='metabolite-category-controls';
    metaboliteCategoryControls.innerHTML='<strong>Metabolite categories</strong> '
      + '<select id="metCategory" aria-label="Category" style="max-width:100%"></select> '
      + '<input id="metCategoryName" aria-label="Category name" placeholder="Category name" maxlength="80"> '
      + '<input id="metCategoryColor" aria-label="Category color" title="Category color" type="color" value="#7f7f7f"> '
      + '<input id="metCategoryMembers" aria-label="Member metabolite IDs" title="Member metabolite IDs, comma-separated" placeholder="Metabolite IDs" style="flex:1 1 280px;min-width:0"> '
      + '<button id="metCategoryApply" type="button">Apply</button> '
      + '<button id="metCategoryAdd" type="button">Add</button> '
      + '<button id="metCategoryDelete" type="button">Delete</button> '
      + '<output id="metCategoryStatus" role="status" style="flex-basis:100%;overflow-wrap:anywhere"></output>';
    mountControl(metaboliteCategoryControls);
    var metaboliteCategoryLegend = document.createElement('div'); metaboliteCategoryLegend.id='metabolite-category-legend';
    metaboliteCategoryControls.appendChild(metaboliteCategoryLegend);
    metaboliteCategoryLegend.style.display='none';
    function refreshMetaboliteCategoryLegend() {
      metaboliteCategoryLegend.textContent='';
      Object.keys(metaboliteCategoryColors).sort().forEach(function(name) {
        var label=document.createElement('span'), swatch=document.createElement('i');
        label.className='metabolite-category-chip';
        label.style='display:flex;align-items:center;gap:4px;overflow-wrap:anywhere;min-width:0';
        swatch.style='width:11px;height:11px;flex-shrink:0;border-radius:50%';
        swatch.style.backgroundColor=metaboliteCategoryColors[name];
        label.appendChild(swatch); label.appendChild(document.createTextNode(name));
        label.title=freeflux_mapmodel.metabolites.filter(function(m){return m.notes.map_info.category===name;}).map(function(m){return m.id;}).join(', ');
        metaboliteCategoryLegend.appendChild(label);
      });
    }
    function refreshStandardMetaboliteLegend() {
      metaboliteCategoryLegend.textContent='';
      [['Default metabolites','#1f77b4'],['ATP',standardMetaboliteColors.ATP],['CO2',standardMetaboliteColors.CO2],['NADH / Q8 / Q8H2',standardMetaboliteColors.NADH],['NADPH',standardMetaboliteColors.NADPH]].forEach(function(entry) {
        var label=document.createElement('span'), swatch=document.createElement('i');
        label.className='metabolite-category-chip'; label.style='display:flex;align-items:center;gap:4px;overflow-wrap:anywhere;min-width:0';
        swatch.style='width:11px;height:11px;flex-shrink:0;border-radius:50%;background-color:'+entry[1];
        label.appendChild(swatch); label.appendChild(document.createTextNode(entry[0])); metaboliteCategoryLegend.appendChild(label);
      });
    }
    function refreshMetaboliteCategoryEditor(selected) {
      var select=document.getElementById('metCategory'); select.innerHTML='';
      Object.keys(metaboliteCategoryColors).sort().forEach(function(name) { var o=document.createElement('option'); o.value=name; o.textContent=name; select.appendChild(o); });
      if(selected && metaboliteCategoryColors[selected]) select.value=selected;
      loadMetaboliteCategoryEditor();
    }
    function loadMetaboliteCategoryEditor() {
      var name=document.getElementById('metCategory').value, color=metaboliteCategoryColors[name] || '#7f7f7f';
      document.getElementById('metCategoryName').value=name || ''; document.getElementById('metCategoryColor').value=color;
      document.getElementById('metCategoryMembers').value=freeflux_mapmodel.metabolites.filter(function(m){ return (m.notes.map_info||{}).category===name; }).map(function(m){return m.id;}).join(',');
      document.getElementById('metCategoryDelete').disabled=name==='Other';
    }
    document.getElementById('metCategory').onchange=loadMetaboliteCategoryEditor;
    function saveMetaboliteCategory(adding) {
      var old=adding ? null : document.getElementById('metCategory').value;
      var name=document.getElementById('metCategoryName').value.trim(), color=document.getElementById('metCategoryColor').value;
      var status=document.getElementById('metCategoryStatus'); status.textContent='';
      if(!name) { status.textContent='Category name is required.'; return; }
      if(old==='Other' && name!==old) { status.textContent='Other is the fallback category and cannot be renamed.'; return; }
      if(name!==old && metaboliteCategoryColors[name]) { status.textContent='Category already exists: '+name; return; }
      var members=document.getElementById('metCategoryMembers').value.split(',').map(function(x){return x.trim();}).filter(Boolean);
      var unknown=members.filter(function(id){return !categoryMetabolites[id];});
      if(unknown.length) { status.textContent='Unknown metabolite IDs: '+unknown.join(', '); return; }
      if(old && old!==name) delete metaboliteCategoryColors[old]; metaboliteCategoryColors[name]=color;
      freeflux_mapmodel.metabolites.forEach(function(m){
        if(members.indexOf(m.id)>=0) m.notes.map_info.category=name;
        else if(old && m.notes.map_info.category===old) m.notes.map_info.category='Other';
      });
      persistMetaboliteCategoryColors(); refreshMetaboliteCategoryEditor(name); refreshMetaboliteCategoryLegend(); repaintMetaboliteCategories();
    }
    document.getElementById('metCategoryApply').onclick=function() { saveMetaboliteCategory(false); };
    document.getElementById('metCategoryAdd').onclick=function() { saveMetaboliteCategory(true); };
    document.getElementById('metCategoryDelete').onclick=function() {
      var name=document.getElementById('metCategory').value; if(!name || name==='Other') return;
      delete metaboliteCategoryColors[name];
      freeflux_mapmodel.metabolites.forEach(function(m){if(m.notes.map_info.category===name)m.notes.map_info.category='Other';});
      document.getElementById('metCategoryStatus').textContent='';
      persistMetaboliteCategoryColors(); refreshMetaboliteCategoryEditor('Other'); refreshMetaboliteCategoryLegend(); repaintMetaboliteCategories();
    };
    persistMetaboliteCategoryColors(); refreshMetaboliteCategoryEditor(); refreshMetaboliteCategoryLegend(); refreshStandardMetaboliteLegend(); metaboliteCategoryLegend.style.display='grid';

    var savedPalettes = freeflux_mapmodel.notes.palette_schemes || {};
    var metabolitePalettes = Object.assign({'Standard': {mode:'standard'}, 'Categories': {mode:'category'}}, savedPalettes.metabolites || {});
    var reactionPalettes = Object.assign({'Standard': Object.assign({}, standardPathwayColors)}, savedPalettes.reactions || {});
    function applyMetabolitePalette(name) {
      var scheme=metabolitePalettes[name]; if(!scheme) return;
      if(name==='Categories' || scheme.mode==='category') {
        refreshMetaboliteCategoryLegend(); repaintMetaboliteCategories(); metaboliteCategoryLegend.style.display='grid';
        return;
      }
      applyStandardPalette();
    }
    function applyReactionPalette(name) {
      var colors=reactionPalettes[name]; if(!colors) return;
      Object.keys(pathwayColors).forEach(function(key) { delete pathwayColors[key]; });
      Object.keys(colors).forEach(function(key) { pathwayColors[key]=colors[key]; });
      pathwayControls.querySelectorAll('input[data-pathway]').forEach(function(input) { input.value=pathwayColors[input.dataset.pathway] || '#999999'; });
      freeflux_mapmodel.reactions.forEach(function(r) { var info=r.notes.map_info; info.color=pathwayColors[info.pathway] || pathwayColors.Other || '#999999'; });
      updatePathwayColors();
    }
    var paletteLibrary=document.createElement('div'); paletteLibrary.id='palette-library';
    paletteLibrary.innerHTML='<strong>Palette schemes</strong>'
      + '<label>Metabolite scheme<select id="metabolitePaletteScheme"></select></label>'
      + '<label>Reaction scheme<select id="reactionPaletteScheme"></select></label>';
    mountControl(paletteLibrary);
    function refreshPaletteSelects() {
      [['metabolitePaletteScheme', metabolitePalettes], ['reactionPaletteScheme', reactionPalettes]].forEach(function(item) {
        var select=document.getElementById(item[0]); select.innerHTML=''; Object.keys(item[1]).sort().forEach(function(name) { var option=document.createElement('option'); option.value=name; option.textContent=name; select.appendChild(option); });
      });
    }
    refreshPaletteSelects();
    document.getElementById('metabolitePaletteScheme').onchange=function() { applyMetabolitePalette(this.value); };
    document.getElementById('reactionPaletteScheme').onchange=function() { applyReactionPalette(this.value); };

    var splitControls=document.createElement('details'); splitControls.id='split-metabolite-controls';
    splitControls.className='control-card';
    splitControls.innerHTML='<summary>Split metabolite nodes</summary>'
      + '<label>Metabolite<div id="splitMetabolites" class="checkbox-list"></div></label>'
      + '<label>Move reactions to the new node<div id="splitReactions" class="checkbox-list"></div></label>'
      + '<button id="splitMetaboliteApply" type="button">Split node</button>'
      + '<output id="splitMetaboliteStatus" role="status"></output>';
    controlRoot.appendChild(splitControls);
    var mergeControls=document.createElement('details'); mergeControls.id='merge-metabolite-controls';
    mergeControls.className='control-card';
    mergeControls.innerHTML='<summary>Merge metabolite nodes</summary>'
      + '<label>Merge into node<select id="mergeTarget"></select></label>'
      + '<label>Copies to merge<div id="mergeMetabolites" class="checkbox-list"></div></label>'
      + '<button id="mergeMetaboliteApply" type="button">Merge selected nodes</button>'
      + '<output id="mergeMetaboliteStatus" role="status"></output>';
    controlRoot.appendChild(mergeControls);
    function connectedReactions(metaboliteId) {
      return freeflux_mapmodel.reactions.filter(function(reaction) {
        return Object.prototype.hasOwnProperty.call(reaction.metabolites, metaboliteId);
      });
    }
    function refreshSplitMetabolites() {
      var select=document.getElementById('splitMetabolites'), previous=[].map.call(select.querySelectorAll('input:checked'), function(o) { return o.value; });
      select.innerHTML='';
      freeflux_mapmodel.metabolites.filter(function(metabolite) {
        return !metabolite.notes.map_info.hidden && (metabolite.notes.map_info.source_metabolite_id || connectedReactions(metabolite.id).length > 1);
      }).sort(function(a,b) { return a.id.localeCompare(b.id); }).forEach(function(metabolite) {
        var label=document.createElement('label'); label.innerHTML='<input type="radio" name="split-metabolite" value="'+metabolite.id+'"> '+(metabolite.notes.map_info.display_name || metabolite.id)+' ('+metabolite.id+')'; label.querySelector('input').checked=previous.indexOf(metabolite.id)>=0; select.appendChild(label);
      });
      refreshSplitReactions();
      refreshMergeMetabolites();
    }
    function refreshSplitReactions() {
      var metaboliteIds=[].map.call(document.querySelectorAll('#splitMetabolites input:checked'), function(o) { return o.value; });
      var select=document.getElementById('splitReactions'); select.innerHTML='';
      var byId={};
      metaboliteIds.forEach(function(id) { connectedReactions(id).forEach(function(reaction) { byId[reaction.id]=reaction; }); });
      Object.keys(byId).sort().map(function(id) { return byId[id]; }).forEach(function(reaction) {
        var label=document.createElement('label'); label.innerHTML='<input type="checkbox" value="'+reaction.id+'"> '+(reaction.notes.map_info.display_name || reaction.id)+' | '+d3.format('.4g')(Math.abs(reaction.notes.map_info.flux || 0)); select.appendChild(label);
      });
      document.getElementById('splitMetaboliteStatus').textContent='';
    }
    function refreshMergeMetabolites() {
      var source=document.getElementById('mergeTarget'), old=source.value, sources={}; source.innerHTML='';
      freeflux_mapmodel.metabolites.filter(function(m) { return m.notes.map_info.source_metabolite_id; }).forEach(function(m) { sources[m.notes.map_info.source_metabolite_id]=true; });
      freeflux_mapmodel.reactions.forEach(function(r) { Object.keys((r.notes.map_info || {}).cofactors || {}).forEach(function(id) { sources[id]=true; }); });
      freeflux_mapmodel.metabolites.filter(function(m) { return sources[m.id]; }).forEach(function(m) { var o=document.createElement('option'); o.value=m.id; o.textContent=m.id; source.appendChild(o); });
      if ([].some.call(source.options, function(o) { return o.value===old; })) source.value=old;
      refreshMergeCopies();
    }
    function refreshMergeCopies() {
      var target=document.getElementById('mergeTarget').value, targetNode=freeflux_mapmodel.metabolites.filter(function(m) { return m.id===target; })[0];
      var source=targetNode && (targetNode.notes.map_info.source_metabolite_id || targetNode.id), box=document.getElementById('mergeMetabolites'); box.innerHTML='';
      freeflux_mapmodel.metabolites.filter(function(m) { return m.id!==target && m.notes.map_info.source_metabolite_id===source; }).forEach(function(m) {
        var label=document.createElement('label'); label.innerHTML='<input type="checkbox" value="copy:'+m.id+'"> '+m.id; box.appendChild(label);
      });
      freeflux_mapmodel.reactions.forEach(function(r) {
        var cofactors=(r.notes.map_info || {}).cofactors || {};
        if (!cofactors[source]) return;
        var label=document.createElement('label'); label.innerHTML='<input type="checkbox" value="cofactor:'+r.id+':'+source+'"> '+source+' | '+r.id+' (cofactor)'; box.appendChild(label);
      });
    }
    document.getElementById('splitMetabolites').onchange=refreshSplitReactions;
    document.getElementById('mergeTarget').onchange=refreshMergeCopies;
    document.getElementById('splitMetaboliteApply').onclick=function() {
      var metaboliteIds=[].map.call(document.querySelectorAll('#splitMetabolites input:checked'), function(o) { return o.value; });
      var selected=[].map.call(document.querySelectorAll('#splitReactions input:checked'), function(o) { return o.value; });
      var status=document.getElementById('splitMetaboliteStatus');
      if (!metaboliteIds.length || !selected.length) { status.textContent='Select metabolites and reactions.'; return; }
      metaboliteIds.forEach(function(metaboliteId, nodeIndex) {
        var source=freeflux_mapmodel.metabolites.filter(function(m) { return m.id===metaboliteId; })[0];
        var copy=JSON.parse(JSON.stringify(source)), index=2;
        while (freeflux_mapmodel.metabolites.some(function(m) { return m.id===source.id+'__copy'+index; })) index++;
        copy.id=source.id+'__copy'+index;
        copy.notes.map_info.source_metabolite_id=source.notes.map_info.source_metabolite_id || source.id;
        copy.notes.map_info.display_name=source.notes.map_info.display_name || source.id;
        copy.notes.map_info.x=Math.min(width, (+source.notes.map_info.x || width/2)+60+nodeIndex*30);
        copy.notes.map_info.y=Math.min(height, (+source.notes.map_info.y || height/2)+60+nodeIndex*30);
        delete copy.notes.map_info.hidden; delete copy.notes.map_info.label_hidden;
        freeflux_mapmodel.metabolites.push(copy);
        freeflux_mapmodel.reactions.forEach(function(reaction) {
          if (selected.indexOf(reaction.id)<0 || reaction.metabolites[metaboliteId]===undefined) return;
          if (reaction.notes && reaction.notes.map_info && reaction.notes.map_info.cofactors) delete reaction.notes.map_info.cofactors[metaboliteId];
          reaction.metabolites[copy.id]=reaction.metabolites[metaboliteId]; delete reaction.metabolites[metaboliteId];
        });
      });
      window.freefluxRebuild();
    };
    document.getElementById('mergeMetaboliteApply').onclick=function() {
      var copyIds=[].map.call(document.querySelectorAll('#mergeMetabolites input:checked'), function(o) { return o.value; });
      var status=document.getElementById('mergeMetaboliteStatus');
      if (!copyIds.length) { status.textContent='Select one or more split copies.'; return; }
      var targetId=document.getElementById('mergeTarget').value;
      if (!targetId || copyIds.indexOf(targetId)>=0) { status.textContent='Choose a target and different copies.'; return; }
      copyIds.forEach(function(copyId) {
        var parts=copyId.split(':'), isCofactor=parts[0]==='cofactor';
        var actualId=isCofactor ? parts.slice(2).join(':') : parts.slice(1).join(':'), sourceId=actualId;
        var copy=freeflux_mapmodel.metabolites.filter(function(m) { return m.id===actualId; })[0];
        freeflux_mapmodel.reactions.forEach(function(reaction) {
          if (reaction.notes && reaction.notes.map_info && reaction.notes.map_info.cofactors && (!isCofactor || reaction.id===parts[1])) delete reaction.notes.map_info.cofactors[sourceId];
          if (isCofactor || reaction.metabolites[actualId]===undefined) return;
          reaction.metabolites[targetId]=(reaction.metabolites[targetId] || 0)+reaction.metabolites[actualId]; delete reaction.metabolites[actualId];
        });
        if (!isCofactor) freeflux_mapmodel.metabolites=freeflux_mapmodel.metabolites.filter(function(m) { return m.id!==actualId; });
      });
      window.freefluxRebuild();
    };
    refreshSplitMetabolites();
    
    var publicationLineMode = LINE_MODE_PLACEHOLDER;
    var controls = document.createElement('div'); controls.id='flux-threshold-controls';
    if (publicationLineMode === 'relative_dash') {
      controls.innerHTML='<label>Bold top <output id="boldValue">20%</output> <input id="boldCut" type="range" min="0" max="100" value="20" step="1"></label> &nbsp; <label>Dashed bottom <output id="dashValue">20%</output> <input id="dashCut" type="range" min="0" max="100" value="20" step="1"></label>';
    } else {
      controls.innerHTML='<strong>Flux encoding</strong><p class="control-hint">Line width encodes flux magnitude. Dashes are reserved for semantic uncertainty / knockouts.</p>';
    }
    mountControl(controls);
    function updateThresholds() {
      var bold=+document.getElementById('boldCut').value/100, dash=+document.getElementById('dashCut').value/100;
      document.getElementById('boldValue').value=Math.round(bold*100)+'%'; document.getElementById('dashValue').value=Math.round(dash*100)+'%';
      svg.selectAll('path.link').each(function(d) {
        var info=d.rxn.notes.map_info, rank=info.relative_rank || 1;
        d3.select(this).attr('stroke-width', rank<=bold ? 6 : get_flux_width(d.rxn));
        d3.select(this).attr('stroke-dasharray', rank>1-dash ? '5,5' : (info.dasharray || 'none'));
      });
      svg.selectAll('.flux-legend text').filter(function(){return this.textContent.indexOf('top')>=0;}).text('Bold top '+Math.round(bold*100)+'%');
      svg.selectAll('.flux-legend text').filter(function(){return this.textContent.indexOf('bottom')>=0;}).text('Dashed bottom '+Math.round(dash*100)+'%');
    }
    if (publicationLineMode === 'relative_dash') {
      document.getElementById('boldCut').oninput=updateThresholds; document.getElementById('dashCut').oninput=updateThresholds;
    }
    var typography=document.createElement('div'); typography.id='typography-controls';
    typography.innerHTML='<strong>Label sizes</strong>'
      + '<label>Metabolites <output id="metFontValue">MET_FONT_PLACEHOLDER</output> pt <input id="metFontSize" type="range" min="10" max="28" value="MET_FONT_PLACEHOLDER" step="1"></label>'
      + '<label>Reactions <output id="rxnFontValue">RXN_FONT_PLACEHOLDER</output> pt <input id="rxnFontSize" type="range" min="8" max="20" value="RXN_FONT_PLACEHOLDER" step="1"></label>';
    mountControl(typography);
    function updateTypography() {
      var met=document.getElementById('metFontSize').value, rxn=document.getElementById('rxnFontSize').value;
      document.getElementById('metFontValue').value=met; document.getElementById('rxnFontValue').value=rxn;
      svg.selectAll('text.nodelabel:not(.rxn):not(.cofactor)').style('font-size', met+'pt');
      svg.selectAll('text.nodelabel.cofactor').style('font-size', Math.max(8, met*.7)+'pt');
      svg.selectAll('text.nodelabel.rxn').style('font-size', rxn+'pt');
    }
    document.getElementById('metFontSize').oninput=updateTypography;
    document.getElementById('rxnFontSize').oninput=updateTypography;
    updateTypography();

    var guideControls=document.createElement('div'); guideControls.id='circle-guide-controls';
    guideControls.innerHTML='<strong>Circle guides</strong>'
      + '<div class="circle-guide-actions"><button id="addCircleGuide" type="button">Add circle</button><button id="deleteCircleGuide" type="button" disabled>Delete selected</button></div>'
      + '<label>Radius <output id="circleRadiusValue">120</output><input id="circleRadius" type="range" min="20" max="500" value="120" step="5" disabled></label>'
      + '<label>Opacity <output id="circleOpacityValue">50%</output><input id="circleOpacity" type="range" min="5" max="100" value="50" step="5" disabled></label>'
      + '<label class="circle-snap-toggle"><input id="circleSnapEnabled" type="checkbox" checked> Snap metabolites and reactions</label>'
      + '<label>Snap distance <output id="circleSnapValue">30</output><input id="circleSnapDistance" type="range" min="5" max="100" value="30" step="5"></label>';
    mountControl(guideControls);
    var guideStore=freeflux_mapmodel.notes.map_info.guide_circles || [];
    freeflux_mapmodel.notes.map_info.guide_circles=guideStore;
    var selectedGuide=null, guideLayer=svg.insert('g', 'g').attr('class', 'circle-guide-layer');
    var guideDrag=d3.behavior.drag().on('dragstart', function(d) {
      d3.event.sourceEvent.stopPropagation(); selectGuide(d);
    }).on('drag', function(d) {
      d.x+=d3.event.dx; d.y+=d3.event.dy; d3.select(this).attr('cx', d.x).attr('cy', d.y);
    });
    function selectGuide(guide) {
      selectedGuide=guide;
      document.getElementById('deleteCircleGuide').disabled=!guide;
      document.getElementById('circleRadius').disabled=!guide;
      document.getElementById('circleOpacity').disabled=!guide;
      if(guide) {
        document.getElementById('circleRadius').value=guide.radius;
        document.getElementById('circleOpacity').value=Math.round(guide.opacity*100);
        document.getElementById('circleRadiusValue').value=guide.radius;
        document.getElementById('circleOpacityValue').value=Math.round(guide.opacity*100)+'%';
      }
      renderGuides();
    }
    function renderGuides() {
      var circles=guideLayer.selectAll('circle.circle-guide').data(guideStore, function(d){return d.id;});
      circles.enter().append('circle').attr('class','circle-guide').call(guideDrag).on('click', function(d){d3.event.stopPropagation();selectGuide(d);});
      circles.attr('cx',function(d){return d.x;}).attr('cy',function(d){return d.y;}).attr('r',function(d){return d.radius;})
        .style('opacity',function(d){return d.opacity;}).classed('selected',function(d){return d===selectedGuide;});
      circles.exit().remove();
    }
    document.getElementById('addCircleGuide').onclick=function() {
      var guide={id:'circle_'+Date.now(), x:width/2, y:height/2, radius:120, opacity:.5}; guideStore.push(guide); selectGuide(guide);
    };
    document.getElementById('deleteCircleGuide').onclick=function() {
      if(!selectedGuide)return; guideStore.splice(guideStore.indexOf(selectedGuide),1); selectGuide(null);
    };
    document.getElementById('circleRadius').oninput=function() {
      if(!selectedGuide)return; selectedGuide.radius=+this.value; document.getElementById('circleRadiusValue').value=this.value; renderGuides();
    };
    document.getElementById('circleOpacity').oninput=function() {
      if(!selectedGuide)return; selectedGuide.opacity=+this.value/100; document.getElementById('circleOpacityValue').value=this.value+'%'; renderGuides();
    };
    document.getElementById('circleSnapDistance').oninput=function() { document.getElementById('circleSnapValue').value=this.value; };
    function nearestGuide(x,y) {
      var nearest=null, gap=Infinity;
      guideStore.forEach(function(guide) {
        var dx=x-guide.x, dy=y-guide.y, distance=Math.sqrt(dx*dx+dy*dy);
        var candidate=Math.abs(distance-guide.radius); if(candidate<gap) { gap=candidate; nearest={guide:guide, dx:dx, dy:dy, distance:distance}; }
      });
      if(nearest) nearest.gap=gap;
      return nearest;
    }
    function placeNode(node,x,y) { node.x=node.px=x; node.y=node.py=y; }
    function placeNodeOnGuide(node,guide,x,y) {
      var dx=x-guide.x, dy=y-guide.y, distance=Math.sqrt(dx*dx+dy*dy), scale=guide.radius/(distance || 1);
      placeNode(node, guide.x+(distance ? dx*scale : guide.radius), guide.y+(distance ? dy*scale : 0));
    }
    window.freefluxMoveNodeWithGuide=function(node,x,y,attachedGuide) {
      if(!node) return null;
      if(!document.getElementById('circleSnapEnabled').checked) { placeNode(node,x,y); return null; }
      var snapDistance=+document.getElementById('circleSnapDistance').value;
      if(attachedGuide) {
        var attachedDistance=Math.sqrt(Math.pow(x-attachedGuide.x,2)+Math.pow(y-attachedGuide.y,2));
        if(Math.abs(attachedDistance-attachedGuide.radius)<=Math.max(snapDistance*2,snapDistance+20)) {
          placeNodeOnGuide(node,attachedGuide,x,y); return attachedGuide;
        }
        placeNode(node,x,y); return null;
      }
      var nearest=nearestGuide(x,y);
      if(nearest && nearest.gap<=snapDistance) { placeNodeOnGuide(node,nearest.guide,x,y); return nearest.guide; }
      placeNode(node,x,y); return null;
    };
    node_drag.on('dragstart.circleGuideSnap', function(d) {
      d._guidePointerX=d.x; d._guidePointerY=d.y;
      var nearest=nearestGuide(d.x,d.y); d._attachedGuide=nearest && nearest.gap<1 ? nearest.guide : null;
    }).on('drag.circleGuideSnap', function(d) {
      d._guidePointerX+=d3.event.dx; d._guidePointerY+=d3.event.dy;
      d._attachedGuide=window.freefluxMoveNodeWithGuide(d,d._guidePointerX,d._guidePointerY,d._attachedGuide);
      guideLayer.selectAll('circle.circle-guide').classed('snap-active',function(guide){return guide===d._attachedGuide;});
      d3.select(this).attr('transform','translate('+d.x+','+d.y+')');
    }).on('dragend.circleGuideSnap', function(d) {
      delete d._guidePointerX; delete d._guidePointerY; delete d._attachedGuide;
      guideLayer.selectAll('circle.circle-guide').classed('snap-active',false);
    });
    renderGuides();
    var legend = svg.append("g").attr("class", "flux-legend")
      .attr("transform", "translate(LEGEND_X_PLACEHOLDER,LEGEND_Y_PLACEHOLDER)").attr("font-family", "Arial").attr("font-size", 11).attr("fill", "#5f6b78");
    var legendLabels = LEGEND_LABELS_PLACEHOLDER;
    if (publicationLineMode === 'flux_width_only') {
      legend.append("line").attr("x1", 0).attr("x2", 38).attr("y1", 0).attr("y2", 0)
        .attr("stroke", "#657786").attr("stroke-width", 3.8).attr("stroke-dasharray", "none");
      legend.append("text").attr("x", 48).attr("y", 4).text(legendLabels[0]);
      legend.append("text").attr("x", 48).attr("y", 23).text(legendLabels[1]);
      legend.append("text").attr("x", 48).attr("y", 42).text(legendLabels[2]);
    } else {
      [0, 1].forEach(function(i) {
        legend.append("line").attr("x1", 0).attr("x2", 35).attr("y1", i*25).attr("y2", i*25)
          .attr("stroke", "#999").attr("stroke-width", 2).attr("stroke-dasharray", i ? "none" : "5,5");
        legend.append("text").attr("x", 45).attr("y", i*25+4).text(legendLabels[i]);
      });
      legend.append("text").attr("y", 60).text(legendLabels[2]);
    }
    '''.replace('LINE_MODE_PLACEHOLDER', json.dumps(line_style['mode']))\
       .replace('MET_FONT_PLACEHOLDER', str(default_met_font))\
       .replace('RXN_FONT_PLACEHOLDER', str(default_rxn_font))\
       .replace('LEGEND_X_PLACEHOLDER', str(float(legend_x)))\
       .replace('LEGEND_Y_PLACEHOLDER', str(float(legend_y)))\
       .replace('LEGEND_LABELS_PLACEHOLDER', json.dumps(legend))\
       .replace('PATHWAY_COLORS_PLACEHOLDER', json.dumps(PATHWAY_COLORS))\
       .replace('METABOLITE_COLORS_PLACEHOLDER', json.dumps(METABOLITE_COLORS))\
       .replace('METABOLITE_CATEGORY_COLORS_PLACEHOLDER', json.dumps(METABOLITE_CATEGORY_COLORS))
    interaction_js = r'''
    var deleteControls = document.createElement('div'); deleteControls.id='delete-controls';
    deleteControls.innerHTML='<button id="deleteMode" type="button" aria-pressed="false">Delete mode: off</button> <button id="undoDelete" type="button" title="Undo last deletion"><span class="button-icon" aria-hidden="true">↶</span> Undo last</button> <button id="restoreHidden" type="button" title="Restore hidden elements"><span class="button-icon" aria-hidden="true">↺</span> Restore hidden</button>';
    mountControl(deleteControls);
    controlRoot.prepend(deleteControls);
    deleteControls.appendChild(document.getElementById('freeflux_map_reactionbutton'));
    var loadControls = document.getElementById('load-json-controls') || document.createElement('div'); loadControls.id='load-json-controls';
    loadControls.innerHTML='<label class="load-json-button" for="load-json">↥ Load map JSON</label><input id="load-json" type="file" accept="application/json">';
    var options=document.getElementById('freeflux_map_options');
    options.insertBefore(loadControls, options.querySelector('.svgbutton'));
    document.getElementById('load-json').onchange=function(event) {
      var file=event.target.files[0]; if(!file) return;
      var reader=new FileReader();
      reader.onload=function() { try {
        var loaded=JSON.parse(reader.result);
        if(!loaded || !Array.isArray(loaded.metabolites) || !Array.isArray(loaded.reactions)) throw new Error('Choose a saved D3Flux model JSON');
        var ids=new Set();
        loaded.metabolites.forEach(function(m) {
          if(!m || typeof m.id!=='string' || ids.has(m.id)) throw new Error('Invalid or duplicate metabolite ID');
          ids.add(m.id); m.notes=m.notes||{}; m.notes.map_info=m.notes.map_info||{};
        });
        var reactionIds=new Set();
        loaded.reactions.forEach(function(r) {
          if(!r || typeof r.id!=='string' || reactionIds.has(r.id) || !r.metabolites || typeof r.metabolites!=='object') throw new Error('Invalid reaction');
          reactionIds.add(r.id);
          Object.keys(r.metabolites).forEach(function(id) { if(!ids.has(id) || !Number.isFinite(r.metabolites[id])) throw new Error('Invalid reaction metabolite: '+id); });
          r.notes=r.notes||{}; r.notes.map_info=r.notes.map_info||{};
        });
        freeflux_mapmodel=loaded; savedSnapshot=undefined; window.freefluxRebuild();
        try { localStorage.setItem('d3flux:' + location.pathname + ':freeflux_map', JSON.stringify(loaded)); }
        catch(error) { alert('Map loaded, but refresh recovery is unavailable: '+error.message); }
      } catch(error) { alert('Unable to load JSON: '+error.message); } };
      reader.readAsText(file);
    };
    var deleteMode = false;
    var deleteHistory = [];
    document.getElementById('deleteMode').onclick = function() { deleteMode=!deleteMode; this.textContent='Delete mode: '+(deleteMode?'on':'off'); this.setAttribute('aria-pressed', String(deleteMode)); svg.classed('delete-mode', deleteMode); };
    document.getElementById('undoDelete').onclick = function() {
      var action = deleteHistory.pop();
      if (!action) return;
      if (action.target && action.target.nodeType) d3.select(action.target).classed('user-hidden', false).style('display', action.displayBefore);
      if (action.data && action.data.notes && action.data.notes.map_info) action.data.notes.map_info.hidden = action.hiddenBefore;
      action.data.notes.map_info.label_hidden=action.labelHiddenBefore;
      action.links.forEach(function(item) { d3.select(item.element).classed('user-hidden', false).style('display', item.display); });
    };
    document.getElementById('restoreHidden').onclick = function() {
      nodes.forEach(function(d) { if (d.notes && d.notes.map_info) d.notes.map_info.hidden=false; });
      svg.selectAll('.node-group, .node, .nodelabel').classed('user-hidden', false).style('display', null);
      svg.selectAll('.link').classed('user-hidden', false).style('display', null);
      deleteHistory.length = 0;
    };
    function enableDelete(selection, hideLinks) { selection.classed('delete-target', true).on('click.delete', function(d) {
      if (!deleteMode) return; d3.event.stopPropagation();
      var links = [];
      if (hideLinks) svg.selectAll('.link').filter(function(link) { return link.rxn === d || link.source === d || link.target === d; }).each(function() {
        links.push({element:this, display:this.style.display || null});
      });
      var displayBefore = this.style.display || null;
      deleteHistory.push({target:this, data:d, displayBefore:displayBefore, hiddenBefore:d.notes && d.notes.map_info ? !!d.notes.map_info.hidden : false, labelHiddenBefore:!!d.notes.map_info.label_hidden, links:links});
      if (hideLinks && d.notes && d.notes.map_info) d.notes.map_info.hidden=true;
      if (!hideLinks) d.notes.map_info.label_hidden=true;
      d3.select(this).classed('user-hidden', true).style('display','none');
      if (hideLinks) svg.selectAll('.link').filter(function(link) { return link.rxn === d || link.source === d || link.target === d; }).style('display','none');
    }); }
    node.classed('node-group', true);
    enableDelete(node, true); enableDelete(text, false);
    text.style('display', function(d) { return d.notes.map_info.label_hidden ? 'none' : null; });

    // v3 pathway/reaction editor. Definitions are kept in the downloadable model JSON.
    var editor = document.createElement('div'); editor.id='pathway-editor';
    editor.innerHTML=''
      + '<select id="editKind"><option value="pathway">Pathway</option><option value="reaction">Reaction</option></select> '
      + '<select id="editItem"></select> <input id="editName" placeholder="name or id"> '
      + '<input id="editColor" type="color" value="#777777"> '
      + '<input id="editPathway" placeholder="pathway"> '
      + '<label>Reactants<select id="editReactants" multiple size="4" aria-label="Reactants"></select></label> '
      + '<label>Products<select id="editProducts" multiple size="4" aria-label="Products"></select></label> '
      + '<input id="editFlux" type="number" step="any" placeholder="flux"> '
      + '<button id="applyEdit" type="button">Apply</button> '
      + '<button id="addEdit" type="button">Add</button> '
      + '<button id="deleteEdit" type="button">Delete</button>';
    var editorDetails=document.createElement('details'); editorDetails.id='pathway-editor-details';
    editorDetails.innerHTML='<summary>Pathway / reaction editor</summary>';
    editorDetails.appendChild(editor); pathwaySection.appendChild(editorDetails);
    window.freefluxRebuild = function() {
      force.stop();
      document.getElementById('freeflux_map_options').appendChild(document.getElementById('freeflux_map_reactionbutton'));
      controlRoot.textContent = '';
      d3.select('#freeflux_map').selectAll('*').remove();
      main(freeflux_mapmodel);
    };
    function pathwayNames() {
      var names={}; reactions.forEach(function(r) { var p=(r.notes&&r.notes.map_info||{}).pathway; if (p) names[p]=1; });
      Object.keys(pathwayColors).forEach(function(p){names[p]=1;}); return Object.keys(names).sort();
    }
    function refreshEditorItems() {
      var kind=document.getElementById('editKind').value, select=document.getElementById('editItem'); select.innerHTML='';
      (kind==='pathway' ? pathwayNames() : reactions.map(function(r){return r.id;})).forEach(function(id) { var o=document.createElement('option'); o.value=id; o.textContent=id; select.appendChild(o); });
      loadEditorItem();
    }
    function loadEditorItem() {
      var kind=document.getElementById('editKind').value, id=document.getElementById('editItem').value;
      var item=kind==='pathway' ? {pathway:id,color:pathwayColors[id]||'#777777'} : reactions.filter(function(r){return r.id===id;})[0];
      var info=item && item.notes ? (item.notes.map_info||{}) : item||{};
      document.getElementById('editName').value=kind==='pathway'?id:(info.display_name||id);
      document.getElementById('editColor').value=info.color||pathwayColors[id]||'#777777';
      document.getElementById('editPathway').value=info.pathway||'';
      refreshReactionMetaboliteOptions();
      ['editReactants','editProducts'].forEach(function(field) {
        var select=document.getElementById(field);
        Array.prototype.forEach.call(select.options, function(option) {
          option.selected = !!(kind==='reaction' && item && item.metabolites &&
            ((field==='editReactants' && item.metabolites[option.value]<0) ||
             (field==='editProducts' && item.metabolites[option.value]>0)));
        });
      });
      document.getElementById('editFlux').value=kind==='reaction'&&item?(info.flux||0):'';
    }
    function refreshReactionMetaboliteOptions() {
      ['editReactants','editProducts'].forEach(function(field) {
        var select=document.getElementById(field), previous={};
        if (!select) return;
        Array.prototype.forEach.call(select.selectedOptions || [], function(option) { previous[option.value]=true; });
        select.innerHTML='';
        freeflux_mapmodel.metabolites.slice().sort(function(a,b) { return a.id.localeCompare(b.id); }).forEach(function(m) {
          var option=document.createElement('option'); option.value=m.id;
          option.textContent=(m.notes && m.notes.map_info && m.notes.map_info.display_name) || m.id;
          option.title=m.id; option.selected=!!previous[m.id]; select.appendChild(option);
        });
      });
    }
    function selectedMetaboliteIds(field) {
      return Array.prototype.map.call(document.getElementById(field).selectedOptions || [], function(option) { return option.value; });
    }
    function repaintPathways() {
      svg.selectAll('path.link').each(function(d) { var info=d.rxn.notes.map_info||{}; d3.select(this).attr('stroke', pathwayColors[info.pathway]||info.color||'#999'); });
      svg.selectAll('marker path').each(function() { var r=d3.select(this.parentNode).datum(), p=r&&r.notes&&r.notes.map_info.pathway; d3.select(this).attr('fill',pathwayColors[p]||'#999'); });
    }
    document.getElementById('editKind').onchange=refreshEditorItems;
    document.getElementById('editItem').onchange=loadEditorItem;
    document.getElementById('applyEdit').onclick=function() {
      var kind=document.getElementById('editKind').value, id=document.getElementById('editItem').value, name=document.getElementById('editName').value.trim(), color=document.getElementById('editColor').value, pathway=document.getElementById('editPathway').value.trim();
      if (!id) return;
      if (kind==='pathway') { pathwayColors[name||id]=color; if(name&&name!==id) delete pathwayColors[id]; reactions.forEach(function(r){var i=r.notes.map_info||{}; if(i.pathway===id){i.pathway=name||id; i.color=color;}}); freeflux_mapmodel.reactions=reactions; persistPathwayColors(); refreshEditorItems(); repaintPathways(); }
      else { var r=reactions.filter(function(x){return x.id===id;})[0]; if(!r)return; r.notes=r.notes||{}; r.notes.map_info=r.notes.map_info||{}; r.notes.map_info.display_name=name||id; r.notes.map_info.pathway=pathway||r.notes.map_info.pathway||'Other'; r.notes.map_info.color=color; r.notes.map_info.flux=+document.getElementById('editFlux').value||0; var mets={}; selectedMetaboliteIds('editReactants').forEach(function(m){mets[m]=-1;}); selectedMetaboliteIds('editProducts').forEach(function(m){mets[m]=1;}); if(Object.keys(mets).length) r.metabolites=mets; freeflux_mapmodel.reactions=reactions; if(window.freefluxRebuild)window.freefluxRebuild(); }
    };
    document.getElementById('addEdit').onclick=function() {
      var kind=document.getElementById('editKind').value, id=document.getElementById('editName').value.trim(), color=document.getElementById('editColor').value, pathway=document.getElementById('editPathway').value.trim()||id;
      if (!id) return;
      if(kind==='pathway') { pathwayColors[id]=color; persistPathwayColors(); refreshEditorItems(); return; }
      if(reactions.some(function(r){return r.id===id;})) return;
      var mets={}; selectedMetaboliteIds('editReactants').forEach(function(m){mets[m]=-1;}); selectedMetaboliteIds('editProducts').forEach(function(m){mets[m]=1;});
      reactions.push({id:id,name:id,metabolites:mets,lower_bound:0,upper_bound:1000,gene_reaction_rule:'',notes:{map_info:{display_name:id,pathway:pathway,color:color,flux:+document.getElementById('editFlux').value||0,x:width/2,y:height/2}}});
      pathwayColors[pathway]=pathwayColors[pathway]||color; if(window.freefluxRebuild) window.freefluxRebuild();
    };
    document.getElementById('deleteEdit').onclick=function() { var kind=document.getElementById('editKind').value,id=document.getElementById('editItem').value; if(!id)return; if(kind==='pathway'){delete pathwayColors[id]; reactions.forEach(function(r){if((r.notes.map_info||{}).pathway===id)r.notes.map_info.pathway='Other';}); freeflux_mapmodel.reactions=reactions; persistPathwayColors(); refreshEditorItems(); repaintPathways();} else {reactions = reactions.filter(function(r){return r.id!==id;}); freeflux_mapmodel.reactions=reactions; if(window.freefluxRebuild)window.freefluxRebuild();} };
    refreshEditorItems();
    // Independent label drag: node drag still moves circle and label together.
    text.each(function(d) {
      var info = d.notes.map_info || {}, dx = +info.label_dx || 0, dy = +info.label_dy || 0;
      d3.select(this).attr('data-label-id', d.id).attr('transform', 'translate('+dx+','+dy+')');
    });
    var label_drag = d3.behavior.drag().on('dragstart', function(d) {
      d3.event.sourceEvent.stopPropagation();
      d._label_dx = +d.notes.map_info.label_dx || 0;
      d._label_dy = +d.notes.map_info.label_dy || 0;
    }).on('drag', function(d) {
      d3.event.sourceEvent.stopPropagation();
      d._label_dx += d3.event.dx; d._label_dy += d3.event.dy;
      d.notes.map_info.label_dx = d._label_dx; d.notes.map_info.label_dy = d._label_dy;
      d3.select(this).attr('transform', 'translate('+d._label_dx+','+d._label_dy+')');
    });
    text.call(label_drag);
    var original_update_node = updateNode;
    updateNode = function() {
      original_update_node.call(this);
      this.select('text').each(function(d) {
        var info=d.notes.map_info||{};
        d3.select(this).attr('transform', 'translate('+(+info.label_dx||0)+','+(+info.label_dy||0)+')');
      });
    };
    '''
    page_source = html.data.replace('var updateNode = function()', 'var updateNode = function()')
    # Keep the renderer's data/geometry untouched; only wrap its existing DOM in a console shell.
    shell = re.search(r'(<div id="freeflux_map">.*?</div>)\s*(<div id="freeflux_map_options">.*?</div>)', page_source, re.S)
    if not shell:
        raise ValueError('d3flux HTML did not contain the expected map/options containers')
    css = '<style>' + Path(__file__).with_name('d3flux').joinpath('templates/console.css').read_text() + '</style>'
    page_source = page_source[:shell.start()] + (
        '<div class="freeflux-app"><aside id="freeflux_sidebar"><header class="console-brand"><span class="brand-mark" aria-hidden="true">F</span><div><div class="eyebrow">METABOLIC ATLAS</div><h1>FreeFlux Studio</h1></div></header><p class="console-intro">Explore the network. Shape the view.</p>'
        + shell.group(2) + '<div id="freeflux_dynamic_controls"></div></aside><main id="freeflux_canvas">'
        + shell.group(1) + '</main></div>'
    ) + page_source[shell.end():]
    page_source = css + page_source
    page_source = page_source.replace(
        'var nodewidth = metabolite_scale(Math.abs(d.notes.map_info.flux));',
        "var normalized = d.notes.map_info.concentration_normalized;\n"
        "        var nodewidth = normalized !== undefined ? 4 + 8 * normalized : metabolite_scale(Math.abs(d.notes.map_info.flux));")
    page_source = page_source.replace('    var updateAnchorLink = function()', legend_js + '\n' + interaction_js + '\n    var updateAnchorLink = function()', 1)
    write_standalone(page_source.replace('var node_drag = force.drag()', 'var node_drag = force.drag()'), args.output)
    map_data = json.loads(create_model_json(model, flux_dict=flux))
    args.output.with_suffix('.json').write_text(
        json.dumps(map_data, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    counts = export_native(args.output, layout['canvas'])
    report = {
        'renderer': 'native d3flux / WebCola / semantic hyperedge + cubic routing',
        'exports': {suffix: str(args.output.with_suffix('.' + suffix).resolve())
                    for suffix in ('html', 'svg', 'png', 'json')},
        'source_sha256': {str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
        'view_only_hidden_reactions': layout['hidden_reactions'],
        'view_only_hidden_metabolites': layout['hidden_metabolites'],
        'negative_irreversible_fluxes_preserved': {r.id: flux[r.id] for r in model.reactions
                                                  if not r.reversibility and flux[r.id] < 0},
        'rendered_elements': counts,
        'line_style': line_style,
        'scientific_layout': scientific_layout,
        'pathway_multiplicity': multiplicity,
        'publication_qa': publication_qa,
        'layout_authority': {
            'layout_coordinates_authoritative': not args.inherit_map_coordinates,
            'initial_map_coordinates_inherited': bool(args.initial_map and args.inherit_map_coordinates),
        },
        'concentration_size': concentration_size,
        'model_reaction_count': len(model.reactions),
        'note': 'Visualization only; biochemical stoichiometry/fluxes are unchanged and no fitting is performed.'
    }
    if layout.get('circular_backbone', {}).get('enabled'):
        browser_pass = counts.get('strict_visual_qa', {}).get('pass', False)
        report['acceptance'] = {
            'pass': bool(publication_qa['pass'] and browser_pass),
            'backend_geometry_pass': bool(publication_qa['pass']),
            'rendered_svg_pass': bool(browser_pass),
            'policy': 'equal-radius closed cycle + full incidences + strict rendered geometry'
        }
    args.output.with_suffix('.manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    if 'acceptance' in report and not report['acceptance']['pass']:
        raise RuntimeError('Circular render failed acceptance; inspect the manifest before using this figure.')


if __name__ == '__main__':
    main()
