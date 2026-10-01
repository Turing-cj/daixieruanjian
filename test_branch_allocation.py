"""Local labels partition real incidences; they never duplicate pool totals."""
import copy
import json
import unittest

import circular_layout
from test_circular_renderer import DEMO
import test_relation_clarity as relation_tests


class BranchAllocationTests(unittest.TestCase):
    def setUp(self):
        self.model, _, _ = relation_tests.RelationClarityTests().compile()
        self.flux = json.loads((DEMO / 'fit_result.json').read_text())['opt_net_fluxes']
        self.policy = {'enabled': True, 'label_mode': 'local_flow', 'normalization_flux': 1}

    def apply(self):
        return circular_layout.apply_pathway_multiplicity(self.model, self.flux, self.policy)

    def info(self, mid, rid=None, kind='local_instances'):
        if rid is None:
            return self.model.metabolites.get_by_id(mid).notes['map_info']
        return self.model.reactions.get_by_id(rid).notes['map_info'][kind][mid]

    def test_f6p_source_is_three_but_each_remote_consumer_is_one(self):
        self.apply()
        self.assertEqual('3 F6P [A]', self.info('F6P')['display_name'])
        self.assertEqual('1 F6P [A]', self.info('F6P', 'IML_TALA')['display_name'])
        self.assertEqual('1 F6P [A]', self.info('F6P', 'IML_TKT2')['display_name'])

    def test_g3p_has_two_local_flows_of_one_not_two_copies_of_two(self):
        self.apply()
        self.assertEqual('1 G3P [B]', self.info('G3P')['display_name'])
        for rid in ('IML_TALA', 'IML_TKT1'):
            self.assertEqual('1 G3P [B]', self.info('G3P', rid)['display_name'])

    def test_shared_view_counts_inflow_or_outflow_not_their_sum(self):
        report = self.apply()
        self.assertIn('view_flows', report)
        local = report['view_flows']['G3P_regeneration']
        self.assertEqual(1, local['produced'])
        self.assertEqual(1, local['consumed'])
        self.assertEqual(1, local['amount'])
        self.assertEqual(2, local['pool_total'])
        self.assertEqual(['IML_TALA', 'IML_TKT1'], [x['reaction'] for x in local['incidences']])

    def test_currency_labels_are_local_reaction_amounts(self):
        self.apply()
        for rid in ('IML_G6PDH2r', 'IML_PGL'):
            self.assertEqual('1 H⁺', self.info('H', rid, 'cofactors')['display_name'])

    def test_pool_totals_remain_available_separately(self):
        report = self.apply()
        self.assertEqual(3, report['counts']['F6P'])
        self.assertEqual(2, report['counts']['G3P'])
        self.assertEqual(2, report['counts']['H'])
        self.assertEqual(3, self.info('F6P', 'IML_TALA').get('pool_total'))

    def test_every_signed_incidence_is_partitioned_exactly_once(self):
        report = self.apply()
        self.assertIn('view_flows', report)
        incidences = [(v['metabolite'], i['reaction'])
                      for v in report['view_flows'].values() for i in v['incidences']]
        expected = [(m.id, r.id) for r in self.model.reactions for m in r.metabolites]
        self.assertCountEqual(expected, incidences)
        for mid, total in report['balances'].items():
            views = [v for v in report['view_flows'].values() if v['metabolite'] == mid]
            for side in ('produced', 'consumed'):
                self.assertAlmostEqual(total[side], sum(v[side] for v in views))

    def test_reverse_flux_retains_amounts_but_reverses_direction(self):
        forward = self.apply()
        self.assertIn('view_flows', forward)
        self.flux = {rid: -v for rid, v in self.flux.items()}
        backward = self.apply()
        for key, old in forward['view_flows'].items():
            new = backward['view_flows'][key]
            self.assertEqual(old['amount'], new['amount'])
            self.assertEqual(old['produced'], new['consumed'])
            self.assertEqual(old['consumed'], new['produced'])

    def test_nonunit_stoichiometry_and_normalization_are_not_ignored(self):
        reaction = self.model.reactions.get_by_id('IML_TALA')
        reaction.add_metabolites({self.model.metabolites.get_by_id('F6P'): -1})
        self.policy['normalization_flux'] = 4
        self.apply()
        self.assertEqual('0.5 F6P [A]', self.info('F6P', 'IML_TALA')['display_name'])
        self.assertEqual('0.25 F6P [A]', self.info('F6P', 'IML_TKT2')['display_name'])

    def test_zero_and_fractional_flux_keep_real_numbers(self):
        self.flux = {rid: v / 8 for rid, v in self.flux.items()}
        self.apply()
        self.assertEqual('0.125 F6P [A]', self.info('F6P', 'IML_TALA')['display_name'])
        self.flux = {rid: 0 for rid in self.flux}
        self.apply()
        self.assertEqual('0 F6P [A]', self.info('F6P', 'IML_TALA')['display_name'])

    def test_invalid_mode_fails_before_mutation(self):
        before = copy.deepcopy(self.model.notes)
        names = [copy.deepcopy(m.notes) for m in self.model.metabolites]
        self.policy['label_mode'] = 'invented'
        with self.assertRaises(ValueError):
            self.apply()
        self.assertEqual(before, self.model.notes)
        self.assertEqual(names, [m.notes for m in self.model.metabolites])

    def test_invalid_flux_fails_before_any_metadata_change(self):
        before = copy.deepcopy([o.notes for o in list(self.model.metabolites)+list(self.model.reactions)])
        self.flux['IML_RPI'] = float('nan')
        with self.assertRaises(ValueError):
            self.apply()
        self.assertEqual(before, [o.notes for o in list(self.model.metabolites)+list(self.model.reactions)])

    def test_reapplication_and_switching_from_legacy_do_not_stack_prefixes(self):
        circular_layout.apply_pathway_multiplicity(self.model, self.flux, {'enabled': True})
        first = copy.deepcopy(self.apply())
        names = [copy.deepcopy(o.notes) for o in list(self.model.metabolites)+list(self.model.reactions)]
        self.assertEqual(first, self.apply())
        self.assertEqual(names, [o.notes for o in list(self.model.metabolites)+list(self.model.reactions)])
        self.assertEqual('1 F6P [A]', self.info('F6P', 'IML_TALA')['display_name'])

    def test_chemistry_flux_and_geometry_do_not_change(self):
        chemistry = {r.id: {m.id: c for m, c in r.metabolites.items()} for r in self.model.reactions}
        flux = copy.deepcopy(self.flux)
        def positions():
            infos = [o.notes['map_info'] for o in list(self.model.metabolites)+list(self.model.reactions)]
            infos += [v for r in self.model.reactions for k in ('local_instances', 'cofactors')
                      for v in r.notes['map_info'].get(k, {}).values()]
            return [{k: v[k] for k in ('x', 'y', 'label_dx', 'label_dy', 'node_radius') if k in v} for v in infos]
        before = positions()
        self.apply()
        self.assertEqual(chemistry, {r.id: {m.id: c for m, c in r.metabolites.items()} for r in self.model.reactions})
        self.assertEqual(flux, self.flux)
        self.assertEqual(before, positions())


if __name__ == '__main__':
    unittest.main()
