"""Regression contracts for rendering the nine archived RuMP pathways."""
import unittest
from cobra import Model, Reaction, Metabolite
from circular_layout import compile_circular_layout
from module47_batch_20260922 import (
    _bounded_circular_pull,
    extra_currency_fan,
    filter_candidate_rows,
    longest_side_path,
    main_side_currency_lane,
    should_merge_shared_instance,
)


class SideReactionContract(unittest.TestCase):
    def test_candidate_filter_preserves_source_order(self):
        rows = [
            {'candidate_id': 'A'},
            {'candidate_id': 'B'},
            {'candidate_id': 'C'},
        ]
        self.assertEqual(
            ['A', 'C'],
            [row['candidate_id'] for row in filter_candidate_rows(rows, {'C', 'A'})],
        )

    def test_candidate_filter_rejects_unknown_ids(self):
        with self.assertRaises(ValueError):
            filter_candidate_rows([{'candidate_id': 'A'}], {'A', 'missing'})

    def make_model(self):
        model = Model('side_contract')
        mets = {name: Metabolite(name, compartment='c')
                for name in ('A', 'B', 'C', 'D', 'ATP', 'ADP', 'H')}
        for rid, stoich in {
            'R1': {'A': -1, 'B': 1}, 'R2': {'B': -1, 'C': 1},
            'R3': {'C': -1, 'A': 1},
            'SIDE': {'C': -1, 'ATP': -1, 'D': 1, 'ADP': 1, 'H': 1},
        }.items():
            rxn = Reaction(rid)
            rxn.add_metabolites({mets[mid]: value for mid, value in stoich.items()})
            model.add_reactions([rxn])
        for obj in list(model.metabolites) + list(model.reactions):
            obj.notes['map_info'] = {'x': 700, 'y': 400}
        for mid in ('ATP', 'ADP', 'H'):
            mets[mid].notes['map_info']['hidden'] = True
        model.reactions.get_by_id('SIDE').notes['map_info']['cofactors'] = {
            mid: {'x': 650 + i*80, 'y': 250}
            for i, mid in enumerate(('ATP', 'ADP', 'H'))}
        layout = {'canvas': [1200, 1000], 'circular_backbone': {
            'enabled': True, 'center': [400, 400], 'radius': 220,
            'items': ['A', 'R1', 'B', 'R2', 'C', 'R3'],
            'side_reactions': [{'reaction': 'SIDE', 'source': 'C', 'target': 'D'}],
            'view_overrides': {'D': {'label_dx': 55, 'label_dy': 9, 'node_radius': 9}},
        }}
        return model, layout

    def test_multisubstrate_side_reaction_preserves_every_incidence(self):
        model, layout = self.make_model()
        rxn = model.reactions.get_by_id('SIDE')
        original = {m.id: c for m, c in rxn.metabolites.items()}
        result = compile_circular_layout(model, layout)
        self.assertEqual(result['circular_backbone']['reaction_count'], 4)
        self.assertEqual({m.id: c for m, c in rxn.metabolites.items()}, original)
        self.assertEqual(set(rxn.notes['map_info']['cofactors']), {'ATP', 'ADP', 'H'})
        self.assertEqual(rxn.notes['map_info']['local_trunk']['source'], 'C')

    def test_explicit_view_style_does_not_change_chemistry(self):
        model, layout = self.make_model()
        compile_circular_layout(model, layout)
        info = model.metabolites.get_by_id('D').notes['map_info']
        self.assertEqual((info['label_dx'], info['label_dy'], info['node_radius']), (55, 9, 9))

    def test_invalid_side_direction_still_rejected(self):
        model, layout = self.make_model()
        layout['circular_backbone']['side_reactions'][0]['source'] = 'D'
        with self.assertRaises(ValueError):
            compile_circular_layout(model, layout)

    def test_nontrunk_carbon_can_be_an_explicit_local_satellite(self):
        model, layout = self.make_model()
        extra = Metabolite('E', compartment='c')
        extra.notes['map_info'] = {'x': 700, 'y': 400}
        model.add_metabolites([extra])
        reaction = model.reactions.get_by_id('SIDE')
        reaction.add_metabolites({extra: 1})
        layout['circular_backbone']['side_reactions'][0]['local_satellites'] = {
            'E': {'position': [820, 330], 'label_offset': [42, -6]},
        }
        original = {m.id: c for m, c in reaction.metabolites.items()}
        result = compile_circular_layout(model, layout)
        info = reaction.notes['map_info']['local_instances']['E']
        self.assertEqual((info['x'], info['y']), (820.0, 330.0))
        self.assertEqual('E', info['display_instance_of'])
        self.assertIn('E', result['circular_backbone']['clone_ids'])
        self.assertEqual(original, {m.id: c for m, c in reaction.metabolites.items()})

    def test_extra_currency_fan_separates_same_sign_rays_and_moves_outward(self):
        placements = extra_currency_fan(
            [0.0, 0.0],
            [600.0, 0.0],
            [-1000.0, 0.0],
            {-1: ['NADP'], 1: ['CO2', 'NADPH']},
            canonical_at_source=True,
        )
        reaction_mid = (300.0, 0.0)
        nadp = placements['NADP']['position']
        co2 = placements['CO2']['position']
        nadph = placements['NADPH']['position']

        # All currency nodes are pushed toward the open end of the stub.
        self.assertGreater(min(nadp[0], co2[0], nadph[0]), reaction_mid[0])
        # Substrates/products occupy opposite lateral half-planes.
        self.assertLess(nadp[1], 0.0)
        self.assertGreater(co2[1], 0.0)
        self.assertGreater(nadph[1], 0.0)
        # Same-sign products must not lie on one ray through the reaction.
        v1 = (co2[0] - reaction_mid[0], co2[1] - reaction_mid[1])
        v2 = (nadph[0] - reaction_mid[0], nadph[1] - reaction_mid[1])
        cross = v1[0] * v2[1] - v1[1] * v2[0]
        self.assertNotAlmostEqual(cross, 0.0)

    def test_currency_fan_stays_outside_reserved_carbon_lane(self):
        placements = extra_currency_fan(
            [0.0, 0.0], [0.0, 600.0], [-1000.0, 300.0],
            {-1: ['H'], 1: []}, canonical_at_source=True,
            lane_floor_by_sign={-1: 150.0},
        )
        # The currency incidence is nearer to the trunk than the local carbon
        # satellite, so it must occupy a more lateral lane to preserve edge
        # ordering rather than crossing the farther carbon ray.
        h = placements['H']['position']
        self.assertGreater(abs(h[0]), 150.0)

    def test_main_side_currency_lane_stays_outside_same_sign_carbon(self):
        self.assertEqual(-32.0, main_side_currency_lane(-1, 0, {-1: 0.0, 1: 0.0}))
        self.assertEqual(32.0, main_side_currency_lane(1, 0, {-1: 0.0, 1: 0.0}))
        self.assertLess(main_side_currency_lane(-1, 0, {-1: 150.0}), -150.0)
        self.assertGreater(main_side_currency_lane(1, 1, {1: 95.0}), 95.0)

    def test_side_path_tie_prefers_real_module_sink(self):
        adj = {
            'm:RING': ['r:R'],
            'r:R': ['m:E4P', 'm:SINK'],
            'm:E4P': [],
            'm:SINK': [],
        }
        path = longest_side_path(adj, {'R'}, {'RING'}, {'SINK'})
        self.assertEqual(['m:RING', 'r:R', 'm:SINK'], path)

    def test_shared_instance_rejects_mixed_source_sink_fanout(self):
        stoich = {
            'PRODUCE': {'M': 1.0},
            'USE_A': {'M': -1.0},
            'USE_B': {'M': -2.0},
            'USE_C': {'M': -1.0},
        }
        self.assertFalse(should_merge_shared_instance(
            'M', ['PRODUCE', 'USE_A', 'USE_B', 'USE_C'], stoich))

    def test_shared_instance_accepts_only_a_homogeneous_pair(self):
        stoich = {
            'USE_A': {'M': -1.0},
            'USE_B': {'M': -2.0},
            'PRODUCE': {'M': 1.0},
        }
        self.assertTrue(should_merge_shared_instance('M', ['USE_A', 'USE_B'], stoich))
        self.assertFalse(should_merge_shared_instance('M', ['USE_A', 'PRODUCE'], stoich))
        self.assertFalse(should_merge_shared_instance('M', ['USE_A'], stoich))

    def test_shared_instance_requires_one_angular_sector(self):
        stoich = {'USE_A': {'M': -1.0}, 'USE_B': {'M': -1.0}}
        self.assertTrue(should_merge_shared_instance(
            'M', ['USE_A', 'USE_B'], stoich,
            {'USE_A': 20.0, 'USE_B': 50.0}, max_separation_deg=54.0))
        self.assertFalse(should_merge_shared_instance(
            'M', ['USE_A', 'USE_B'], stoich,
            {'USE_A': 126.0, 'USE_B': 234.0}, max_separation_deg=54.0))
        # Circular wrap-around is still one local sector.
        self.assertTrue(should_merge_shared_instance(
            'M', ['USE_A', 'USE_B'], stoich,
            {'USE_A': 350.0, 'USE_B': 10.0}, max_separation_deg=54.0))

    def test_shared_instance_rejects_bow_tie_with_second_shared_carbon_endpoint(self):
        stoich = {
            'R1': {'M': -1.0, 'X': 1.0, 'A': -1.0},
            'R2': {'M': -1.0, 'X': 1.0, 'B': -1.0},
        }
        self.assertFalse(should_merge_shared_instance(
            'M', ['R1', 'R2'], stoich,
            {'R1': 100.0, 'R2': 130.0}, max_separation_deg=54.0))

    def test_secondary_attachment_cannot_drag_shared_hub_out_of_primary_sector(self):
        # A residual side branch at 130 degrees used to drag a connector whose
        # backbone sector is near 50 degrees far enough to create a bow-tie.
        self.assertAlmostEqual(54.0, _bounded_circular_pull(50.0, 70.0, 4.0))
        # Circular wrap-around must use the short angular direction.
        self.assertAlmostEqual(354.0, _bounded_circular_pull(350.0, 20.0, 4.0))


if __name__ == '__main__':
    unittest.main()
