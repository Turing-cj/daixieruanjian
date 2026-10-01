"""Pool totals come from signed stoichiometric fluxes, never view counts."""
import copy
import json
import unittest

import circular_layout
import test_relation_clarity as relation_tests
from test_circular_renderer import DEMO


class PathwayMultiplicityTests(unittest.TestCase):
    def setUp(self):
        self.model, _, _ = relation_tests.RelationClarityTests().compile()
        self.flux = json.loads((DEMO / 'fit_result.json').read_text())['opt_net_fluxes']
        self.policy = {'enabled': True, 'normalization_flux': 1.0}

    def apply(self):
        return circular_layout.apply_pathway_multiplicity(self.model, self.flux, self.policy)

    def test_gross_pool_totals_are_not_net_flux_or_number_of_views(self):
        result = self.apply()
        expected = {'Formaldehyde':3, 'Ru5P':3, 'H6P':3, 'F6P':3,
                    'G3P':2, 'Xu5P':2, 'H':2}
        for metabolite in self.model.metabolites:
            self.assertEqual(expected.get(metabolite.id, 1), result['counts'][metabolite.id])
        self.assertEqual({'produced':3.0, 'consumed':3.0, 'net':0.0}, result['balances']['F6P'])

    def test_all_carbon_views_inherit_canonical_count(self):
        self.apply()
        for reaction in self.model.reactions:
            for mid, view in reaction.notes['map_info'].get('local_instances', {}).items():
                canonical = self.model.metabolites.get_by_id(mid).notes['map_info']
                self.assertEqual(canonical['pathway_multiplicity'], view['pathway_multiplicity'])
                self.assertEqual(canonical['display_name'], view['display_name'])

    def test_currency_totals_preserve_local_chemical_names(self):
        self.apply()
        view = self.model.reactions.get_by_id('IML_G6PDH2r').notes['map_info']['cofactors']['NADP']
        self.assertEqual('1 NADP⁺', view['display_name'])

    def test_reapplication_is_idempotent_and_does_not_stack_prefixes(self):
        self.apply()
        before = copy.deepcopy(self.model.notes)
        names = [m.notes['map_info']['display_name'] for m in self.model.metabolites]
        self.apply()
        self.assertEqual(before, self.model.notes)
        self.assertEqual(names, [m.notes['map_info']['display_name'] for m in self.model.metabolites])

    def test_signed_reverse_flux_and_nonunit_stoichiometry(self):
        reaction = self.model.reactions.get_by_id('IML_RPE')
        xu = self.model.metabolites.get_by_id('Xu5P')
        reaction.add_metabolites({xu:-1})  # now consumes 2 Xu5P per turnover
        self.flux = {rid:-value for rid,value in self.flux.items()}
        result = self.apply()
        self.assertEqual(4, result['counts']['Xu5P'])
        self.assertEqual(4, result['balances']['Xu5P']['produced'])
        self.assertEqual(2, result['balances']['Xu5P']['consumed'])

    def test_fractional_and_zero_flux_are_not_rounded_to_integer_molecules(self):
        self.flux = {rid:value/8 for rid,value in self.flux.items()}
        self.apply()
        self.assertEqual('0.375 F6P [A]', self.model.metabolites.get_by_id('F6P').notes['map_info']['display_name'])
        self.flux = {rid:0 for rid in self.flux}
        self.apply()
        self.assertEqual('0 F6P [A]', self.model.metabolites.get_by_id('F6P').notes['map_info']['display_name'])

    def test_normalization_is_explicit_and_consistent(self):
        self.flux = {rid:value*10 for rid,value in self.flux.items()}
        self.policy['normalization_flux'] = 10
        self.assertEqual(3, self.apply()['counts']['F6P'])

    def test_invalid_inputs_fail_before_any_label_is_changed(self):
        before = [copy.deepcopy(m.notes) for m in self.model.metabolites]
        for value in (float('nan'), float('inf'), None):
            self.flux['IML_RPI'] = value
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.apply()
            self.assertEqual(before, [m.notes for m in self.model.metabolites])
        del self.flux['IML_RPI']
        with self.assertRaises(ValueError):
            self.apply()

    def test_invalid_normalization_is_rejected(self):
        for value in (0, -1, float('inf'), float('nan')):
            self.policy['normalization_flux'] = value
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.apply()

    def test_original_stoichiometry_flux_and_geometry_are_unchanged(self):
        chemistry = {r.id:{m.id:c for m,c in r.metabolites.items()} for r in self.model.reactions}
        flux = copy.deepcopy(self.flux)
        geometry = [(o.id, o.notes['map_info']['x'], o.notes['map_info']['y'])
                    for o in list(self.model.metabolites)+list(self.model.reactions)]
        self.apply()
        self.assertEqual(chemistry, {r.id:{m.id:c for m,c in r.metabolites.items()} for r in self.model.reactions})
        self.assertEqual(flux, self.flux)
        self.assertEqual(geometry, [(o.id, o.notes['map_info']['x'], o.notes['map_info']['y'])
                                  for o in list(self.model.metabolites)+list(self.model.reactions)])

    def test_disabled_policy_keeps_legacy_labels_untouched(self):
        self.policy['enabled'] = False
        self.assertFalse(self.apply()['enabled'])
        self.assertEqual('F6P [A]', self.model.metabolites.get_by_id('F6P').notes['map_info']['display_name'])


if __name__ == '__main__':
    unittest.main()
