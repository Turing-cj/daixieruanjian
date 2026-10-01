"""Real SVG quantities are audited against actual incidence endpoints."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from playwright.sync_api import sync_playwright
import test_circular_browser as legacy

ROOT = Path(__file__).resolve().parent


class BranchBrowserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='branch_test_', dir=legacy.DEMO)
        cls.directory = Path(cls.temporary.name)
        cls.html = cls.directory / 'test.html'
        layout = json.loads((legacy.DEMO / 'layout.json').read_text(encoding='utf-8'))
        layout['pathway_multiplicity']['label_mode'] = 'local_flow'
        path = cls.directory / 'layout.json'
        path.write_text(json.dumps(layout), encoding='utf-8')
        result = subprocess.run([sys.executable, '-X', 'utf8', str(ROOT / 'render_network.py'),
            str(legacy.DEMO / 'fit_result.json'), str(legacy.DEMO / 'model.tsv'), str(cls.html),
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

    def test_correct_local_flow_scene_passes_independent_audit(self):
        qa = self.page.evaluate(self.auditor)
        self.assertTrue(qa['pass'], qa)
        self.assertEqual(1, qa['visual_carbon_component_count'])

    def test_actual_svg_contains_source_three_and_consumer_one(self):
        labels = self.page.locator('g.node-group').evaluate_all('''gs=>Object.fromEntries(gs
          .filter(g=>g.__data__ && g.__data__.type!=='rxn')
          .map(g=>[g.__data__.id,[...g.querySelectorAll('tspan.text')].map(t=>t.textContent).join('')]))''')
        self.assertEqual('3 F6P [A]', labels['F6P'])
        self.assertEqual('1 F6P [A]', labels['F6P_IML_TALA'])
        self.assertEqual('1 F6P [A]', labels['F6P_IML_TKT2'])
        self.assertEqual('1 G3P [B]', labels['G3P'])
        self.assertEqual('1 G3P [B]', labels['G3P_regeneration'])

    def test_repeating_pool_total_on_a_branch_is_rejected_even_with_matching_cache(self):
        qa = self.inject_and_audit("""
          const g=[...document.querySelectorAll('g.node-group')].find(g=>g.__data__.id==='F6P_IML_TALA');
          g.__data__.notes.map_info.pathway_multiplicity=3;
          g.querySelector('tspan.text').textContent='3 F6P [A]';
        """)
        self.assertFalse(qa['pass'])
        self.assertTrue(qa['semantic_errors'])

    def test_shared_flow_must_not_double_count_in_and_out(self):
        qa = self.inject_and_audit("""
          const g=[...document.querySelectorAll('g.node-group')].find(g=>g.__data__.id==='G3P_regeneration');
          g.__data__.notes.map_info.pathway_multiplicity=2;
          g.querySelector('tspan.text').textContent='2 G3P [B]';
        """)
        self.assertFalse(qa['pass'])
        self.assertTrue(qa['semantic_errors'])

    def test_wrong_separate_pool_total_metadata_is_rejected(self):
        qa = self.inject_and_audit("""
          const g=[...document.querySelectorAll('g.node-group')].find(g=>g.__data__.id==='F6P_IML_TALA');
          g.__data__.notes.map_info.pool_total=9;
        """)
        self.assertFalse(qa['pass'])
        self.assertTrue(qa['semantic_errors'])

    def test_save_load_preserves_local_amounts_and_all_geometry(self):
        before, geometry = self.snapshot(), self.geometry()
        path = self.directory / 'saved_local.json'
        with self.page.expect_download() as download:
            self.page.locator('#freeflux_map_options .download').click()
        download.value.save_as(path)
        self.load_saved_map(path)
        self.assertEqual(before, self.snapshot())
        self.assertEqual(geometry, self.geometry())
        self.assertTrue(self.page.evaluate(self.auditor)['pass'])

    def test_reversing_all_fluxes_preserves_local_amounts_and_audit(self):
        snapshot = self.snapshot()
        for r in snapshot['reactions']:
            r['notes']['map_info']['flux'] *= -1
        path = self.directory / 'reversed_local.json'
        path.write_text(json.dumps(snapshot), encoding='utf-8')
        self.load_saved_map(path)
        qa = self.page.evaluate(self.auditor)
        self.assertTrue(qa['pass'], qa)


if __name__ == '__main__':
    unittest.main()
