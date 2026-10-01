"""Regression contracts for connected regeneration and explicit pool aliases."""
import math
import unittest

import test_circular_renderer as circle_tests


def connected_layout():
    layout = circle_tests.circle_layout()
    spec = layout['circular_backbone']
    spec.update(shared_metabolites={'Xu5P': {'position': [600, 770],
                                           'label_offset': [36, -20]}},
                instance_marker_style='pool_code', input_arrows=True)
    spec['shared_instances'] = {
        'G3P_regeneration': {'metabolite': 'G3P', 'reactions': ['IML_TALA', 'IML_TKT1'],
                             'position': [475, 900], 'label_offset': [45, -6]}
    }
    spec['side_reactions'] = [
        {'reaction': 'IML_RPE', 'source': 'Xu5P', 'target': 'Ru5P',
         'start': [600, 770],
         'end': [760 + 475 * math.cos(math.radians(-147)),
                 710 + 475 * math.sin(math.radians(-147))],
         'control_points': [[570, 640], [445, 485]], 'label_offset': [48, -8]}
    ]
    return layout


class RelationClarityTests(unittest.TestCase):
    def compile(self, layout=None):
        return circle_tests.CircularRendererTests().compile(layout or connected_layout())

    def test_adjacent_tal_tkt1_share_one_g3p_view_without_moving_ring_g3p(self):
        model, summary, _ = self.compile()
        for rid in ('IML_TALA', 'IML_TKT1'):
            entry = model.reactions.get_by_id(rid).notes['map_info']['local_instances']['G3P']
            self.assertEqual('G3P_regeneration', entry.get('view_instance_id'))
            self.assertEqual([475, 900], [entry['x'], entry['y']])
        self.assertEqual(3, summary['circular_backbone'].get('visible_local_instance_count'))

    def test_shared_instance_cannot_invent_an_incidence(self):
        layout = connected_layout()
        layout['circular_backbone']['shared_instances']['G3P_regeneration']['reactions'].append('IML_PGI')
        with self.assertRaisesRegex(ValueError, 'shared instance'):
            self.compile(layout)

    def test_tkt_products_share_the_actual_rpe_substrate(self):
        model, _, _ = self.compile()
        for rid in ('IML_TKT1', 'IML_TKT2', 'IML_RPE'):
            self.assertNotIn('Xu5P', model.reactions.get_by_id(rid).notes['map_info'].get('local_instances', {}))
        info = model.metabolites.get_by_id('Xu5P').notes['map_info']
        self.assertEqual([600, 770], [info['x'], info['y']])
        self.assertFalse(info.get('clone_marker', False))

    def test_rpe_returns_to_the_actual_outer_ring_receptor(self):
        model, _, _ = self.compile()
        info = model.reactions.get_by_id('IML_RPE').notes['map_info']
        self.assertNotIn('Ru5P', info.get('local_instances', {}))
        self.assertEqual([[570, 640], [445, 485]], info['local_trunk'].get('control_points'))
        self.assertEqual('Ru5P', info['local_trunk']['target'])
        receptor = model.metabolites.get_by_id('Ru5P').notes['map_info']
        self.assertFalse(receptor.get('clone_marker', False))

    def test_remaining_repeated_pools_have_distinct_explicit_codes(self):
        model, summary, _ = self.compile()
        self.assertEqual(['F6P', 'G3P'], summary['circular_backbone']['clone_ids'])
        self.assertEqual({'F6P': 'A', 'G3P': 'B'}, summary['circular_backbone'].get('pool_codes'))
        for mid, code in [('F6P', 'A'), ('G3P', 'B')]:
            name = model.metabolites.get_by_id(mid).notes['map_info']['display_name']
            self.assertEqual(f'{mid} [{code}]', name)
            for reaction in model.reactions:
                local = reaction.notes['map_info'].get('local_instances', {}).get(mid)
                if local:
                    self.assertEqual(name, local['display_name'])

    def test_shared_pool_cannot_move_a_backbone_node_or_hidden_currency(self):
        for mid in ('Ru5P', 'H2O', 'unknown'):
            layout = connected_layout()
            layout['circular_backbone']['shared_metabolites'][mid] = {'position': [600, 800]}
            with self.subTest(mid=mid), self.assertRaisesRegex(ValueError, 'shared'):
                self.compile(layout)

    def test_side_reaction_cannot_silently_displace_the_outer_circle(self):
        layout = connected_layout()
        layout['circular_backbone']['side_reactions'][0]['end'] = [885, 790]
        with self.assertRaisesRegex(ValueError, 'backbone'):
            self.compile(layout)

    def test_explicit_input_arrow_policy_survives_compilation(self):
        _, summary, _ = self.compile()
        self.assertTrue(summary['circular_backbone'].get('input_arrows'))

    def test_side_endpoints_can_reuse_existing_node_positions(self):
        layout = connected_layout()
        side = layout['circular_backbone']['side_reactions'][0]
        del side['start'], side['end']
        try:
            self.compile(layout)
        except KeyError as error:
            self.fail(f'Existing biochemical endpoints must not require duplicated coordinates: {error}')


if __name__ == '__main__':
    unittest.main()
