"""Readable publication defaults must reach every view without changing chemistry."""
import copy
import json
import unittest
from pathlib import Path
from nine_paths_pdf_20260921.build_batch import make_layout, CURRENCY


class PublicationReadabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pathways = [json.loads(line) for line in
            (Path(__file__).parent / 'nine_paths_pdf_20260921/source/pathways.jsonl').read_text(encoding='utf-8').splitlines()
            if line.strip()]

    def test_all_nine_layouts_use_larger_publication_defaults(self):
        for index, pathway in enumerate(self.pathways, 1):
            original = copy.deepcopy(pathway)
            layout = make_layout(index, pathway)
            with self.subTest(candidate=index):
                self.assertEqual(layout['publication_style']['metabolite_font_pt'], 24)
                self.assertEqual(layout['publication_style']['reaction_font_pt'], 17)
                self.assertEqual(layout['line_style']['width_range'], [3.2, 7.6])
                for mid, info in layout['metabolites'].items():
                    self.assertEqual(info['node_radius'], 7 if mid in CURRENCY else 12)
                self.assertEqual(pathway, original)

    def test_profile_is_idempotent_and_does_not_move_nodes_or_change_the_primary_ring(self):
        from publication_readability import enlarge_publication
        layout = make_layout(2, self.pathways[1])
        before = copy.deepcopy(layout)
        self.assertEqual(enlarge_publication(layout), layout)
        self.assertEqual(layout, before)


if __name__ == '__main__':
    unittest.main()
