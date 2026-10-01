"""The primary cycle is circular; an output sequence must not become a second ring."""
import copy
import json
import unittest
from pathlib import Path
from nine_paths_pdf_20260921.build_batch import make_layout


class OpenSidePathwayContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).parent / 'nine_paths_pdf_20260921/source/pathways.jsonl'
        cls.pathways = [json.loads(line) for line in path.read_text(encoding='utf-8').splitlines() if line.strip()]

    def side_chain(self, index):
        layout = make_layout(index, self.pathways[index-1])
        chain = [s for s in layout['circular_backbone']['side_reactions']
                 if s['reaction'] not in ('IML_RPE', 'IML_ALKP')]
        return layout, chain

    def test_open_product_chains_do_not_use_circular_arcs(self):
        for index in (2, 5):
            with self.subTest(candidate=index):
                _, chain = self.side_chain(index)
                self.assertTrue(chain)
                self.assertTrue(all('control_points' not in side for side in chain),
                                'An open output sequence is still being forced onto a large circular arc')
                deltas = [side['end'][1] - side['start'][1] for side in chain]
                self.assertTrue(all(d * deltas[0] > 0 for d in deltas))

    def test_returning_pool_is_unfolded_without_a_second_geometric_ring(self):
        for index in (3, 6, 7, 8, 9):
            with self.subTest(candidate=index):
                layout, chain = self.side_chain(index)
                self.assertEqual(chain[0]['source'], chain[-1]['target'])
                self.assertNotEqual(chain[0]['start'], chain[-1]['end'],
                                    'The same metabolite pool may have two labelled views; no second ring is necessary')
                self.assertNotIn('G3P_output', layout['circular_backbone']['shared_instances'])
                self.assertTrue(chain[0]['local_source'])
                self.assertTrue(chain[-1]['local_target'])

    def test_every_side_step_has_one_continuous_forward_trunk(self):
        for index in (2, 3, 5, 6, 7, 8, 9):
            with self.subTest(candidate=index):
                _, chain = self.side_chain(index)
                for side in chain:
                    self.assertEqual(side['start'][0], side['end'][0])
                    self.assertGreaterEqual(abs(side['end'][1] - side['start'][1]), 120)
                for left, right in zip(chain, chain[1:]):
                    self.assertEqual(left['target'], right['source'])
                    self.assertEqual(left['end'], right['start'])

    def test_source_chemistry_and_flux_are_not_mutated(self):
        for index, pathway in enumerate(self.pathways, 1):
            before = copy.deepcopy(pathway)
            layout = make_layout(index, pathway)
            self.assertEqual(pathway, before)
            self.assertEqual(layout['line_style']['width_range'], [3.2, 7.6])
            self.assertEqual(layout['publication_style']['metabolite_font_pt'], 24)
            self.assertEqual(layout['publication_style']['reaction_font_pt'], 17)

    def test_existing_single_circle_candidates_remain_single_circle(self):
        for index in (1, 4):
            layout, chain = self.side_chain(index)
            self.assertEqual(chain, [])
            self.assertEqual(layout['circular_backbone']['radius'], 500 if index == 1 else 600)

    def test_trunk_labels_face_away_from_extra_carbon_branches(self):
        for index in (3, 6, 7, 8, 9):
            layout, chain = self.side_chain(index)
            view = layout['circular_backbone']['view_overrides']
            self.assertGreater(view['PEP']['label_dx'], 0)
            dhap_view = 'DHAP_IML_DHAPT' if index == 6 else 'DHAP'
            self.assertGreater(view[dhap_view]['label_dx'], 0)

    def test_moved_bridge_keeps_currency_node_clear_of_carbon_output(self):
        import math
        layout, _ = self.side_chain(6)
        pi = layout['reactions']['IML_ALKP']['cofactors']['Pi']
        view = layout['circular_backbone']['view_overrides']
        pyr = view['Pyruvate']
        self.assertGreater(math.dist([pi['x'], pi['y']], [pyr['x'], pyr['y']]), 100)
        self.assertLess(view['DHAP']['label_dx'], 0)


if __name__ == '__main__':
    unittest.main()
