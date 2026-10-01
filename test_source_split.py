"""Source identity is not the same thing as a metabolite pool identity."""
import copy
import unittest

import circular_layout
import test_branch_allocation as branch_tests


class SourceSplitTests(unittest.TestCase):
    info = branch_tests.BranchAllocationTests.info
    apply = branch_tests.BranchAllocationTests.apply

    def setUp(self):
        branch_tests.BranchAllocationTests.setUp(self)
        self.policy['relationship_mode'] = 'source_split'

    def test_f6p_children_explicitly_reference_the_parent_not_plain_pool_aliases(self):
        self.apply()
        self.assertEqual('3 F6P [A]', self.info('F6P')['display_name'])
        self.assertEqual('1 F6P [A.1]', self.info('F6P', 'IML_TALA')['display_name'])
        self.assertEqual('1 F6P [A.2]', self.info('F6P', 'IML_TKT2')['display_name'])

    def test_g3p_independent_production_has_distinct_roots_not_child_ids(self):
        self.apply()
        self.assertEqual('1 G3P [B1]', self.info('G3P')['display_name'])
        self.assertEqual('1 G3P [B2]', self.info('G3P', 'IML_TALA')['display_name'])
        self.assertEqual('1 G3P [B2]', self.info('G3P', 'IML_TKT1')['display_name'])

    def test_audit_includes_main_consumption_and_both_remote_portions(self):
        report = self.apply()
        self.assertIn('source_split', report)
        pool = report['source_split']['pools']['F6P']
        self.assertEqual(['RUMP_R05339'], [s['reaction'] for s in pool['sources']])
        self.assertEqual(3, pool['sources'][0]['amount'])
        self.assertEqual(3, sum(x['amount'] for x in pool['consumption_allocations']))
        self.assertEqual({'IML_PGI', 'IML_TALA', 'IML_TKT2'},
                         {x['reaction'] for x in pool['consumption_allocations']})
        self.assertTrue(pool['single_source_partition'])

    def test_two_sources_never_receive_invented_source_to_consumer_assignments(self):
        report = self.apply()
        self.assertIn('source_split', report)
        pool = report['source_split']['pools']['G3P']
        self.assertFalse(pool['single_source_partition'])
        self.assertEqual(['IML_EDA', 'IML_TALA'], [s['reaction'] for s in pool['sources']])
        for use in pool['consumption_allocations']:
            self.assertEqual([], use['parent_source_ids'])
            self.assertEqual(['B1', 'B2'], use['possible_source_ids'])

    def test_multi_source_pool_use_exposes_possible_source_roots_in_label(self):
        from source_split import derive_source_split
        report = derive_source_split({
            'source_a': {'metabolite': 'M', 'produced': 1, 'consumed': 0,
                         'incidences': [{'reaction': 'PROD_A', 'signed_rate': 1}]},
            'source_b': {'metabolite': 'M', 'produced': 1, 'consumed': 0,
                         'incidences': [{'reaction': 'PROD_B', 'signed_rate': 1}]},
            'consumer': {'metabolite': 'M', 'produced': 0, 'consumed': 2,
                         'incidences': [{'reaction': 'USE', 'signed_rate': -2}]},
        }, {'M': 'Q'})
        self.assertEqual('Q1+Q2:pool', report['views']['consumer']['marker'])
        self.assertEqual([], report['views']['consumer']['parent_source_ids'])
        self.assertEqual(['Q1', 'Q2'], report['views']['consumer']['possible_source_ids'])

    def test_reversed_flux_changes_sources_instead_of_preserving_false_parent(self):
        self.flux = {rid: -value for rid, value in self.flux.items()}
        report = self.apply()
        self.assertIn('source_split', report)
        pool = report['source_split']['pools']['F6P']
        self.assertEqual(3, len(pool['sources']))
        self.assertFalse(pool['single_source_partition'])
        self.assertEqual('3 F6P [A1+A2+A3:pool]', self.info('F6P')['display_name'])

    def test_zero_flux_does_not_invent_a_source_or_parent(self):
        self.flux = {rid: 0 for rid in self.flux}
        report = self.apply()
        self.assertIn('source_split', report)
        self.assertEqual([], report['source_split']['pools']['F6P']['sources'])
        self.assertEqual('0 F6P [A:0]', self.info('F6P')['display_name'])

    def test_reapplication_and_round_trip_modes_do_not_stack_markers(self):
        first = copy.deepcopy(self.apply())
        labels = [copy.deepcopy(o.notes) for o in list(self.model.metabolites)+list(self.model.reactions)]
        self.assertEqual(first, self.apply())
        self.assertEqual(labels, [o.notes for o in list(self.model.metabolites)+list(self.model.reactions)])
        self.policy.pop('relationship_mode')
        self.apply()
        self.assertEqual('1 F6P [A]', self.info('F6P', 'IML_TALA')['display_name'])
        self.policy['relationship_mode'] = 'source_split'
        self.apply()
        self.assertEqual('1 F6P [A.1]', self.info('F6P', 'IML_TALA')['display_name'])

    def test_invalid_relationship_mode_fails_before_mutation(self):
        before = copy.deepcopy([o.notes for o in list(self.model.metabolites)+list(self.model.reactions)])
        self.policy['relationship_mode'] = 'invented'
        with self.assertRaises(ValueError):
            self.apply()
        self.assertEqual(before, [o.notes for o in list(self.model.metabolites)+list(self.model.reactions)])

    def test_switch_back_to_legacy_pool_totals_removes_source_markers(self):
        self.apply()
        circular_layout.apply_pathway_multiplicity(self.model, self.flux, {'enabled': True})
        self.assertEqual('3 F6P [A]', self.info('F6P', 'IML_TALA')['display_name'])
        self.assertNotIn('flow_relation', self.info('F6P', 'IML_TALA'))

    def test_missing_supply_does_not_claim_single_source_partition(self):
        self.flux['RUMP_R05339'] = 1
        report = self.apply()
        self.assertTrue(report['source_split']['pools']['F6P']['boundary_import_required'])
        self.assertFalse(report['source_split']['pools']['F6P']['single_source_partition'])
        self.assertEqual('1 F6P [A:pool]', self.info('F6P', 'IML_TALA')['display_name'])

    def test_generic_merged_sources_are_not_misrepresented_as_a_split(self):
        from source_split import derive_source_split
        views = {'generic_view': {'metabolite': 'M', 'produced': 2, 'consumed': 2,
                 'incidences': [{'reaction': 'produce_z', 'signed_rate': 1},
                                {'reaction': 'produce_a', 'signed_rate': 1},
                                {'reaction': 'consume', 'signed_rate': -2}]}}
        report = derive_source_split(views, {'M': 'Q'})
        self.assertEqual('Q1+Q2', report['views']['generic_view']['marker'])
        self.assertEqual(['produce_a', 'produce_z'], report['views']['generic_view']['source_reactions'])
        self.assertEqual([], report['pools']['M']['consumption_allocations'][0]['parent_source_ids'])


if __name__ == '__main__':
    unittest.main()
