"""Audit source/split semantics from the actual SVG endpoints, not caches."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from playwright.sync_api import sync_playwright
import test_circular_browser as legacy

ROOT = Path(__file__).resolve().parent


class SourceSplitBrowserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='source_split_test_', dir=legacy.DEMO)
        cls.directory = Path(cls.temporary.name)
        cls.html = cls.directory / 'test.html'
        source = ROOT / 'nine_paths_branch_amounts_20260922/candidate_01'
        layout = json.loads((source / 'layout.json').read_text(encoding='utf-8'))
        layout['pathway_multiplicity']['relationship_mode'] = 'source_split'
        path = cls.directory / 'layout.json'
        path.write_text(json.dumps(layout), encoding='utf-8')
        result = subprocess.run([sys.executable, '-X', 'utf8', str(ROOT / 'render_network.py'),
            str(source / 'fit_result.json'), str(source / 'model.tsv'), str(cls.html),
            '--layout', str(path)], cwd=ROOT, capture_output=True, text=True, encoding='utf-8', timeout=90)
        if not cls.html.exists():
            raise RuntimeError(result.stdout[-5000:] + result.stderr)
        cls.playwright = sync_playwright().start()
        cls.browser = cls.playwright.chromium.launch()
        cls.auditor = (ROOT / 'browser_circle_qa.js').read_text(encoding='utf-8')

    @classmethod
    def tearDownClass(cls):
        cls.browser.close()
        cls.playwright.stop()
        cls.temporary.cleanup()

    setUp = legacy.CircularBrowserTests.setUp
    tearDown = legacy.CircularBrowserTests.tearDown
    wait_for_scene = legacy.CircularBrowserTests.wait_for_scene
    inject_and_audit = legacy.CircularBrowserTests.inject_and_audit
    snapshot = legacy.CircularBrowserTests.snapshot
    geometry = legacy.CircularBrowserTests.geometry
    load_saved_map = legacy.CircularBrowserTests.load_saved_map

    def test_correct_source_split_svg_passes_independent_semantic_audit(self):
        qa = self.page.evaluate(self.auditor)
        self.assertFalse(qa['semantic_errors'], qa['semantic_errors'])
        self.assertEqual('source_split', qa.get('relationship_mode'))
        self.assertEqual(1, qa['visual_carbon_component_count'])

    def test_actual_labels_distinguish_children_and_independent_sources(self):
        labels = self.page.locator('g.node-group').evaluate_all('''gs=>Object.fromEntries(gs
          .filter(g=>g.__data__ && g.__data__.type!=='rxn')
          .map(g=>[g.__data__.id,[...g.querySelectorAll('tspan.text')].map(t=>t.textContent).join('')]))''')
        self.assertEqual('3 F6P [A]', labels['F6P'])
        self.assertEqual('1 F6P [A.1]', labels['F6P_IML_TALA'])
        self.assertEqual('1 F6P [A.2]', labels['F6P_IML_TKT2'])
        self.assertEqual('1 G3P [B1]', labels['G3P'])
        self.assertEqual('1 G3P [B2]', labels['G3P_regeneration'])

    def test_wrong_parent_is_rejected_even_when_label_and_cache_agree(self):
        qa = self.inject_and_audit("""
          const g=[...document.querySelectorAll('g.node-group')].find(g=>g.__data__.id==='F6P_IML_TALA');
          const info=g.__data__.notes.map_info;
          info.flow_relation.marker='B1.1'; info.flow_relation.parent_source_ids=['B1'];
          info.pathway_label_base='F6P [B1.1]';
          g.querySelector('tspan.text').textContent='1 F6P [B1.1]';
        """)
        self.assertFalse(qa['pass'])
        self.assertTrue(any('source/split' in str(error) for error in qa['semantic_errors']))

    def test_independent_source_cannot_be_disguised_as_child(self):
        qa = self.inject_and_audit("""
          const g=[...document.querySelectorAll('g.node-group')].find(g=>g.__data__.id==='G3P_regeneration');
          const info=g.__data__.notes.map_info;
          info.flow_relation.marker='B1.1'; info.flow_relation.role='split_use';
          info.flow_relation.parent_source_ids=['B1']; info.flow_relation.source_ids=[];
          info.pathway_label_base='G3P [B1.1]';
          g.querySelector('tspan.text').textContent='1 G3P [B1.1]';
        """)
        self.assertFalse(qa['pass'])
        self.assertTrue(any('source/split' in str(error) for error in qa['semantic_errors']))

    def test_save_load_preserves_relations_and_geometry(self):
        before, geometry = self.snapshot(), self.geometry()
        path = self.directory / 'saved_source_split.json'
        with self.page.expect_download() as download:
            self.page.locator('#freeflux_map_options .download').click()
        download.value.save_as(path)
        self.load_saved_map(path)
        self.assertEqual(before, self.snapshot())
        self.assertEqual(geometry, self.geometry())
        self.assertFalse(self.page.evaluate(self.auditor)['semantic_errors'])


if __name__ == '__main__':
    unittest.main()
