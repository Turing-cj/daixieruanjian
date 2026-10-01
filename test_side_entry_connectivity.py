"""Side-chain layout must not split an existing primary carbon outlet."""
import copy
import json
import unittest
from pathlib import Path

from circular_layout import compile_circular_layout
from scripts_render_d3flux import read_model, apply_layout
from nine_paths_pdf_20260921.build_batch import make_layout, ALIASES, CURRENCY
from open_side_layout import unfold_side_chain

ROOT = Path(__file__).resolve().parent


class SideEntryConnectivityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pathways = [json.loads(line) for line in
            (ROOT / 'nine_paths_pdf_20260921/source/pathways.jsonl').read_text(encoding='utf-8').splitlines()
            if line.strip()]

    def test_primary_outlet_is_not_replaced_by_a_disconnected_source_clone(self):
        checked = 0
        for index, pathway in enumerate(self.pathways, 1):
            layout = make_layout(index, pathway)
            spec = layout['circular_backbone']
            ids = spec.get('side_pathway_presentation', {}).get('reaction_ids', [])
            if not ids:
                continue
            chain = [next(s for s in spec['side_reactions'] if s['reaction'] == rid) for rid in ids]
            source = chain[0]['source']
            primary = {ALIASES[mid] for r in pathway['reactions'] if r['reaction_id'] in spec['items'][1::2]
                       for mid in r['full_stoichiometry'] if ALIASES[mid] not in CURRENCY}
            occurrences = [source] + [s['target'] for s in chain]
            if source in primary and source not in spec['items'][::2] and occurrences.count(source) == 1:
                checked += 1
                with self.subTest(candidate=index, source=source):
                    self.assertFalse(chain[0].get('local_source'),
                                     'An existing primary outlet was split into an unconnected local view')
                    model = read_model(ROOT / f'nine_paths_pdf_20260921/candidate_{index:02d}/model.tsv')
                    apply_layout(model, copy.deepcopy(layout))
                    compile_circular_layout(model, layout)
                    self.assertFalse(any(source in r.notes['map_info'].get('local_instances', {})
                                         for r in model.reactions))
                    info = model.metabolites.get_by_id(source).notes['map_info']
                    self.assertEqual([info['x'], info['y']], chain[0]['start'])
        self.assertGreater(checked, 0)

    def test_entry_policy_does_not_depend_on_candidate_or_metabolite_names(self):
        layout = make_layout(2, self.pathways[1], side_pathway_mode='legacy_circle')
        # Rename all biochemical identifiers; geometry and incidences are unchanged.
        names = {mid: 'pool_' + str(i) for i, mid in enumerate(layout['metabolites'])}
        names.update({rid: 'step_' + str(i) for i, rid in enumerate(layout['reactions'])})
        def rename(value):
            if isinstance(value, str):
                if value in names:
                    return names[value]
                for mid in layout['metabolites']:
                    for rid in layout['reactions']:
                        if value == mid + '_' + rid:
                            return names[mid] + '_' + names[rid]
                return value
            if isinstance(value, list):
                return [rename(v) for v in value]
            if isinstance(value, dict):
                return {rename(k): rename(v) for k, v in value.items()}
            return value
        renamed = rename(layout)
        stoich = {names[r['reaction_id']]: {names[ALIASES[mid]]: c for mid, c in r['full_stoichiometry'].items()}
                  for r in self.pathways[1]['reactions']}
        spec = renamed['circular_backbone']
        ids = [names[rid] for rid in ('IML_MGSA', 'IML_GLYOX3', 'IML_LDH_D')]
        protected = {names[ALIASES[mid]] for r in self.pathways[1]['reactions']
                     if r['reaction_id'] in layout['circular_backbone']['items'][1::2]
                     for mid in r['full_stoichiometry'] if ALIASES[mid] not in CURRENCY}
        before = copy.deepcopy(renamed)
        result = unfold_side_chain(renamed, ids, stoich, protected)
        self.assertEqual(renamed, before)
        entry = next(s for s in result['circular_backbone']['side_reactions'] if s['reaction'] == ids[0])
        self.assertFalse(entry.get('local_source'), 'Connectivity repair must work without known candidate IDs')


if __name__ == '__main__':
    unittest.main()
