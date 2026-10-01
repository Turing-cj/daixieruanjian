"""Contracts for a real circular pathway, not merely an ellipse with zero crossings."""
import copy
import json
import math
import unittest
from pathlib import Path

from scripts_render_d3flux import read_model, apply_layout
from scientific_layout import apply_scientific_layout, build_publication_qa

ROOT = Path(__file__).resolve().parent
DEMO = ROOT / 'candidate01_demo'
CIRCLE = {
    'enabled': True, 'center': [760, 710], 'radius': 475,
    'start_angle_deg': -147,
    'items': ['Ru5P', 'RUMP_R05338', 'H6P', 'RUMP_R05339',
              'F6P', 'IML_PGI', 'G6P', 'IML_G6PDH2r',
              'SixPGL', 'IML_PGL', 'SixPG', 'IML_EDD',
              'KDPG', 'IML_EDA', 'G3P', 'IML_TKT2',
              'E4P', 'IML_TALA', 'S7P', 'IML_TKT1',
              'R5P', 'IML_RPI'],
    'external_metabolites': ['Formaldehyde', 'Pyruvate'],
    'satellite_inset': 130, 'satellite_spacing': 65,
    'side_reactions': [
        {'reaction': 'IML_RPE', 'source': 'Xu5P', 'target': 'Ru5P',
         'start': [650, 790], 'end': [885, 790], 'local_target': True}
    ],
    'title': 'RuMP-ED cycle', 'subtitle': 'Candidate 01',
    'title_position': [760, 615]
}


def circle_layout():
    layout = json.loads((DEMO / 'layout.json').read_text(encoding='utf-8'))
    layout['circular_backbone'] = copy.deepcopy(CIRCLE)
    layout['canvas'] = [1520, 1500]
    return layout


class CircularRendererTests(unittest.TestCase):
    def compile(self, layout=None):
        layout = layout or circle_layout()
        model = read_model(DEMO / 'model.tsv')
        apply_layout(model, layout)
        summary = apply_scientific_layout(model, layout)
        return model, summary, layout

    def test_backbone_is_a_complete_equal_radius_circle(self):
        model, summary, _ = self.compile()
        self.assertIn('circular_backbone', summary)
        for item in CIRCLE['items']:
            obj = (model.reactions if item in model.reactions else model.metabolites).get_by_id(item)
            info = obj.notes['map_info']
            self.assertAlmostEqual(475, math.hypot(info['x'] - 760, info['y'] - 710), places=7)
        self.assertEqual(360, summary['circular_backbone']['sweep_deg'])
        last = model.reactions.get_by_id('IML_RPI').notes['map_info']['circle_arc']
        self.assertEqual('Ru5P', last['target'])

    def test_two_to_two_has_one_trunk_and_local_instances(self):
        model, _, _ = self.compile()
        for rid in ('IML_TKT2', 'IML_TALA', 'IML_TKT1'):
            info = model.reactions.get_by_id(rid).notes['map_info']
            self.assertEqual('circle_hyperedge', info.get('route_style'))
            self.assertEqual(2, len(info['local_instances']))
            self.assertNotIn('edge_routes', info)
            for mid, instance in info['local_instances'].items():
                self.assertEqual(mid, instance['display_instance_of'])
                self.assertLess(math.hypot(instance['x']-info['x'], instance['y']-info['y']), 190)

    def test_no_biochemical_reactions_or_coefficients_are_changed(self):
        original = read_model(DEMO / 'model.tsv')
        model, _, _ = self.compile()
        stoich = lambda m: {r.id: {s.id:c for s,c in r.metabolites.items()} for r in m.reactions}
        self.assertEqual(stoich(original), stoich(model))
        self.assertEqual(12, len(model.reactions))

    def test_all_reactions_have_explicit_nonpaired_routing(self):
        model, summary, _ = self.compile()
        self.assertEqual(12, summary.get('circular_backbone', {}).get('reaction_count'))
        self.assertTrue(all(r.notes['map_info'].get('route_style') in
                            ('circle_hyperedge', 'local_hyperedge') for r in model.reactions))

    def test_invalid_topological_cycle_is_rejected(self):
        layout = circle_layout()
        layout['circular_backbone']['items'][2] = 'Pyruvate'
        with self.assertRaisesRegex(ValueError, 'backbone|substrate|product'):
            self.compile(layout)

    def test_nonfinite_or_nonpositive_radius_is_rejected(self):
        for radius in (0, -1, float('nan'), float('inf')):
            layout = circle_layout()
            layout['circular_backbone']['radius'] = radius
            with self.subTest(radius=radius), self.assertRaises(ValueError):
                self.compile(layout)

    def test_qa_rejects_one_node_displaced_from_circle(self):
        model, _, layout = self.compile()
        self.assertTrue(build_publication_qa(model, layout)['pass'])
        model.metabolites.get_by_id('F6P').notes['map_info']['x'] += 35
        self.assertFalse(build_publication_qa(model, layout)['pass'])

    def test_repeated_carbon_metabolites_are_not_called_cofactors(self):
        model, _, _ = self.compile()
        info = model.reactions.get_by_id('IML_TKT2').notes['map_info']
        self.assertIn('F6P', info.get('local_instances', {}))
        self.assertNotIn('F6P', info.get('cofactors', {}))
        self.assertTrue(info['local_instances']['F6P']['clone_marker'])

    def test_canonical_satellite_uses_real_node_for_single_backbone_incidence(self):
        layout = circle_layout()
        layout['circular_backbone']['external_metabolites'].remove('Pyruvate')
        layout['circular_backbone']['canonical_satellites'] = ['Pyruvate']
        layout['circular_backbone']['local_carbon_outset'] = 100
        model, summary, _ = self.compile(layout)
        reaction = model.reactions.get_by_id('IML_EDA').notes['map_info']
        self.assertNotIn('Pyruvate', reaction.get('local_instances', {}))
        info = model.metabolites.get_by_id('Pyruvate').notes['map_info']
        self.assertGreater(math.hypot(info['x'] - 760, info['y'] - 710), 475)
        self.assertEqual(['Pyruvate'], summary['circular_backbone']['canonical_satellites'])

    def test_canonical_satellite_rejects_multiple_backbone_incidences(self):
        layout = circle_layout()
        layout['circular_backbone']['canonical_satellites'] = ['Xu5P']
        with self.assertRaisesRegex(ValueError, 'exactly one backbone incidence'):
            self.compile(layout)


if __name__ == '__main__':
    unittest.main()
