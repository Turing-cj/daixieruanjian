"""Real Chromium regression tests, including deliberate visual counterexamples."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from playwright.sync_api import sync_playwright
from scripts_render_d3flux import read_model

ROOT = Path(__file__).resolve().parent
DEMO = ROOT / 'candidate01_demo'


class CircularBrowserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='browser_test_', dir=DEMO)
        cls.directory = Path(cls.temporary.name)
        cls.html = cls.directory / 'test.html'
        command = [sys.executable, '-X', 'utf8', str(ROOT / 'render_network.py'),
                   str(DEMO / 'fit_result.json'), str(DEMO / 'model.tsv'), str(cls.html),
                   '--layout', str(DEMO / 'layout.json')]
        result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, encoding='utf-8', timeout=90)
        if result.returncode:
            raise RuntimeError(result.stdout[-7000:] + result.stderr)
        cls.playwright = sync_playwright().start()
        cls.browser = cls.playwright.chromium.launch()
        cls.auditor = (ROOT / 'browser_circle_qa.js').read_text(encoding='utf-8')

    @classmethod
    def tearDownClass(cls):
        cls.browser.close()
        cls.playwright.stop()
        cls.temporary.cleanup()

    def setUp(self):
        self.errors = []
        self.page = self.browser.new_page(viewport={'width':1520, 'height':1580})
        self.page.on('pageerror', lambda error:self.errors.append(str(error)))
        self.page.goto(self.html.as_uri())
        self.wait_for_scene()

    def tearDown(self):
        self.page.close()
        self.assertEqual([], self.errors)

    def wait_for_scene(self):
        self.page.wait_for_function('''() => {
          const paths=[...document.querySelectorAll('path.link')];
          if(!window.freefluxAuditSnapshot || paths.length!==26 || paths.some(p=>!p.getAttribute('d'))) return false;
          const geometry=paths.map(p=>p.getAttribute('d')).join('');
          if(window.qaGeometry!==geometry) {window.qaGeometry=geometry;window.qaStable=performance.now();}
          return performance.now()-window.qaStable>250;
        }''', timeout=15000)

    def inject_and_audit(self, mutation):
        return self.page.evaluate('() => {' + mutation + ';return (' + self.auditor + ')();}')

    def snapshot(self):
        return self.page.evaluate('() => window.freefluxAuditSnapshot()')

    def geometry(self):
        return self.page.locator('path.link').evaluate_all("ps=>ps.map(p=>p.getAttribute('d'))")

    def load_saved_map(self, path):
        # Wait for authoritative rebuilt DOM, not the still-valid old scene
        # while FileReader has not completed its asynchronous callback.
        self.page.evaluate("() => window.qaPreviousSvg=document.querySelector('#freeflux_map svg')")
        self.page.locator('#load-json').set_input_files(path)
        self.page.wait_for_function("() => document.querySelector('#freeflux_map svg')!==window.qaPreviousSvg")
        self.wait_for_scene()

    def test_rendered_scene_and_all_biochemical_incidences(self):
        qa = self.page.evaluate(self.auditor)
        self.assertTrue(qa['pass'], qa)
        self.assertEqual(11, qa['backbone_arc_count'])
        original = read_model(DEMO / 'model.tsv')
        actual = self.snapshot()
        self.assertEqual({r.id:{m.id:c for m,c in r.metabolites.items()} for r in original.reactions},
                         {r['id']:r['metabolites'] for r in actual['reactions']})
        flux = json.loads((DEMO / 'fit_result.json').read_text())['opt_net_fluxes']
        self.assertEqual(flux, {r['id']:r['notes']['map_info']['flux'] for r in actual['reactions']})

    def test_pathway_multiplicity_is_shown_on_every_visible_view(self):
        labels = self.page.locator('g.node-group').evaluate_all('''groups=>groups
          .filter(g=>g.__data__ && g.__data__.type!=='rxn')
          .map(g=>({id:g.__data__.notes.orig_id || g.__data__.id,
            label:[...g.querySelectorAll('text.nodelabel tspan.text')].map(t=>t.textContent).join('\\n')}))''')
        expected = {'Formaldehyde':'3 Formaldehyde', 'Ru5P':'3 Ru5P',
                    'H6P':'3 H6P', 'F6P':'3 F6P [A]', 'G6P':'1 G6P',
                    'SixPGL':'1 6PGL', 'SixPG':'1 6PG', 'KDPG':'1 KDPG',
                    'G3P':'2 G3P [B]', 'Pyruvate':'1 Pyruvate', 'E4P':'1 E4P',
                    'Xu5P':'2 Xu5P', 'S7P':'1 S7P', 'R5P':'1 R5P',
                    'NADP':'1 NADP⁺', 'NADPH':'1 NADPH', 'H':'2 H⁺', 'H2O':'1 H₂O'}
        self.assertEqual(set(expected), {entry['id'] for entry in labels})
        for entry in labels:
            with self.subTest(view=entry):
                self.assertEqual(expected[entry['id']], entry['label'])

    def test_wrong_numeric_prefix_is_rejected_even_when_pool_code_is_correct(self):
        qa = self.inject_and_audit("const g=[...document.querySelectorAll('g.node-group')].find(g=>g.__data__.id==='F6P');g.querySelector('tspan.text').textContent='9 F6P [A]'")
        self.assertFalse(qa['pass'])
        self.assertTrue(qa['semantic_errors'])

    def test_wrong_count_on_nonrepeated_node_is_rejected(self):
        qa = self.inject_and_audit("const g=[...document.querySelectorAll('g.node-group')].find(g=>g.__data__.id==='H6P');g.querySelector('tspan.text').textContent='9 H6P'")
        self.assertFalse(qa['pass'])
        self.assertTrue(qa['semantic_errors'])

    def test_consistently_wrong_counts_are_recomputed_from_chemistry(self):
        qa = self.inject_and_audit("[...document.querySelectorAll('g.node-group')].filter(g=>g.__data__ && (g.__data__.notes.orig_id || g.__data__.id)==='F6P').forEach(g=>{g.__data__.notes.map_info.pathway_multiplicity=9;g.querySelector('tspan.text').textContent='9 F6P [A]'})")
        self.assertFalse(qa['pass'])
        self.assertTrue(any('multiplicity' in str(error) for error in qa['semantic_errors']))

    def test_actual_download_and_load_preserve_instances_and_geometry(self):
        before = self.snapshot()
        geometry = self.geometry()
        with self.page.expect_download() as download:
            self.page.locator('#freeflux_map_options .download').click()
        saved = self.directory / 'saved_map.json'
        download.value.save_as(saved)
        self.assertEqual(before, json.loads(saved.read_text(encoding='utf-8')))
        self.load_saved_map(saved)
        self.assertEqual(before, self.snapshot())
        self.assertEqual(geometry, self.geometry())
        self.assertTrue(self.page.evaluate(self.auditor)['pass'])

    def test_reverse_flux_keeps_stoichiometry_and_reverses_every_incidence(self):
        original = self.snapshot()
        for reaction in original['reactions']:
            reaction['notes']['map_info']['flux'] *= -1
        reversed_map = self.directory / 'reversed.json'
        reversed_map.write_text(json.dumps(original), encoding='utf-8')
        self.load_saved_map(reversed_map)
        self.assertEqual(original, self.snapshot())
        qa = self.page.evaluate(self.auditor)
        self.assertTrue(qa['pass'], qa)

    def test_unequal_screen_axis_scaling_is_rejected(self):
        qa = self.inject_and_audit("document.querySelector('#freeflux_map svg').style.transform='scaleX(1.1)'")
        self.assertFalse(qa['pass'])
        self.assertGreater(qa['screen_axis_ratio'], 1.05)

    def test_standalone_svg_retains_caption_typography(self):
        # Exported SVG has no HTML body to inherit the editor's font from.
        self.page.goto(self.html.with_suffix('.svg').as_uri())
        fonts = self.page.locator('.circle-caption text').evaluate_all(
            'elements=>elements.map(e=>getComputedStyle(e).fontFamily)')
        self.assertEqual(4, len(fonts))
        self.assertTrue(all('Arial' in font for font in fonts), fonts)

    def test_elliptical_path_transform_is_rejected(self):
        qa = self.inject_and_audit("document.querySelector('path.link').setAttribute('transform','scale(1.08,1)')")
        self.assertFalse(qa['pass'])
        self.assertGreater(qa['maximum_rendered_radial_error_px'], 1)

    def test_missing_carbon_branch_is_rejected(self):
        qa = self.inject_and_audit("[...document.querySelectorAll('path.link')].find(p=>p.__data__.rxn.id==='IML_TKT2' && p.__data__.part==='input').remove()")
        self.assertFalse(qa['pass'])
        self.assertTrue(qa['semantic_errors'])

    def test_gap_inside_a_path_is_rejected_even_when_chemistry_is_complete(self):
        qa = self.inject_and_audit("""
          const p=[...document.querySelectorAll('path.link')].find(p=>p.__data__.part==='output');
          const length=p.getTotalLength(), points=Array.from({length:101},(_,i)=>p.getPointAtLength(length*i/100));
          p.setAttribute('d',points.map((q,i)=>i>40 && i<60 ? '' : (i===0 || i===60 ? 'M' : 'L')+q.x+','+q.y).join(' '));
        """)
        self.assertFalse(qa['pass'], 'A visually broken line passed acceptance')
        self.assertTrue(qa.get('path_continuity_errors'), qa)

    def test_unregistered_visual_endpoint_is_rejected_outside_declared_shared_pools(self):
        qa = self.inject_and_audit("""
          const p=[...document.querySelectorAll('path.link')].find(p=>p.__data__.part==='output' && p.__data__.target.id==='Pyruvate');
          const n=p.__data__.target;
          p.__data__.target=Object.assign({},n,{id:'detached_product',notes:Object.assign({},n.notes,{orig_id:'Pyruvate'})});
        """)
        self.assertFalse(qa['pass'], 'An endpoint with no matching visible node passed acceptance')
        self.assertTrue(qa.get('visual_connectivity_errors'), qa)

    def test_archived_disconnected_output_chains_are_rejected(self):
        for index in (2, 5):
            path=ROOT / f'nine_paths_open_chain_20260922/candidate_{index:02d}/candidate_{index:02d}.html'
            self.page.goto(path.as_uri())
            self.page.wait_for_function("() => window.freefluxAuditSnapshot && document.querySelectorAll('path.link').length>20 && [...document.querySelectorAll('path.link')].every(p=>p.getAttribute('d'))")
            self.page.wait_for_timeout(400)
            qa=self.page.evaluate(self.auditor)
            with self.subTest(candidate=index):
                self.assertFalse(qa['pass'])
                self.assertEqual(qa['biochemical_carbon_component_count'], 1)
                self.assertEqual(qa['visual_carbon_component_count'], 2)
                self.assertTrue(qa['visual_connectivity_errors'])

    def test_label_clearance_preserves_all_paths_and_node_positions(self):
        geometry = self.geometry()
        before = self.page.locator('g.node-group').evaluate_all('gs=>gs.map(g=>[g.__data__.id,g.__data__.x,g.__data__.y])')
        self.page.evaluate("() => {const e=document.querySelector('text.nodelabel[data-label-id=H6P]');e.__data__.notes.map_info.label_dx=0;e.__data__.notes.map_info.label_dy=0;e.setAttribute('transform','translate(0,0)');}")
        self.assertFalse(self.page.evaluate(self.auditor)['pass'])
        result = self.page.evaluate((ROOT / 'label_clearance.js').read_text(encoding='utf-8'))
        self.assertTrue(result['paths_unchanged'])
        self.assertEqual(result['unresolved'], [])
        self.assertEqual(geometry, self.geometry())
        self.assertEqual(before, self.page.locator('g.node-group').evaluate_all('gs=>gs.map(g=>[g.__data__.id,g.__data__.x,g.__data__.y])'))
        self.assertTrue(self.page.evaluate(self.auditor)['pass'])

    def test_duplicate_trunk_is_rejected(self):
        qa = self.inject_and_audit("const p=document.querySelector('path.link');const copy=p.cloneNode(true);copy.__data__=p.__data__;p.parentNode.appendChild(copy)")
        self.assertFalse(qa['pass'])
        self.assertTrue(qa['overlaps'])
        self.assertTrue(qa['semantic_errors'])

    def test_clipped_text_is_rejected(self):
        qa = self.inject_and_audit("document.querySelector('.circle-caption text').setAttribute('x','-100')")
        self.assertFalse(qa['pass'])
        self.assertTrue(qa['clipped_labels'])

    def test_overlapping_labels_are_rejected(self):
        qa = self.inject_and_audit("const t=document.querySelector('.circle-caption text');t.parentNode.appendChild(t.cloneNode(true))")
        self.assertFalse(qa['pass'])
        self.assertTrue(qa['label_collisions'])

    def test_label_overlapping_its_own_visible_node_is_rejected(self):
        qa = self.inject_and_audit("document.querySelector('circle.metabolite#H6P').setAttribute('r','50')")
        self.assertFalse(qa['pass'])
        self.assertTrue(qa['label_node_hits'])

    def test_consumption_branch_arrow_is_rejected(self):
        qa = self.inject_and_audit("const p=[...document.querySelectorAll('path.link')].find(p=>p.__data__.part==='input');p.setAttribute('marker-end','url(#unexpected)')")
        self.assertFalse(qa['pass'])
        self.assertTrue(qa['semantic_errors'])

    def test_missing_consumption_arrow_is_rejected_when_explicit(self):
        qa = self.inject_and_audit("const p=[...document.querySelectorAll('path.link')].find(p=>p.__data__.part==='input');p.removeAttribute('marker-end')")
        self.assertFalse(qa['pass'])
        self.assertTrue(qa['semantic_errors'])

    def test_shared_xu5p_cannot_be_replaced_with_an_unconnected_view(self):
        qa = self.inject_and_audit("const p=[...document.querySelectorAll('path.link')].find(p=>p.__data__.rxn.id==='IML_TKT2' && p.__data__.part==='output');const n=p.__data__.target;p.__data__.target=Object.assign({},n,{id:'Xu5P_disconnected',notes:Object.assign({},n.notes,{orig_id:'Xu5P'})})")
        self.assertFalse(qa['pass'])
        self.assertTrue(qa['semantic_errors'])

    def test_wrong_display_pool_code_is_rejected(self):
        qa = self.inject_and_audit("const g=[...document.querySelectorAll('g.node-group')].find(g=>g.__data__.id==='F6P');g.querySelector('tspan.text').textContent='F6P [B]'")
        self.assertFalse(qa['pass'])
        self.assertTrue(qa['semantic_errors'])

    def test_regeneration_uses_shared_visible_nodes(self):
        actual = self.page.evaluate('''() => {
          const ps=[...document.querySelectorAll('path.link')].map(p=>p.__data__);
          const tkt1=ps.find(d=>d.rxn.id==='IML_TKT1' && d.part==='output');
          const tkt2=ps.find(d=>d.rxn.id==='IML_TKT2' && d.part==='output');
          const rpe=ps.find(d=>d.rxn.id==='IML_RPE');
          const tal=ps.find(d=>d.rxn.id==='IML_TALA' && d.part==='output');
          const g3p=ps.find(d=>d.rxn.id==='IML_TKT1' && d.part==='input');
          return {xu:tkt1.target===tkt2.target && tkt1.target===rpe.source,
                  ru:rpe.target.id, g3p:tal.target===g3p.source,
                  sharedG3p:[...document.querySelectorAll('g.node-group')].filter(g=>g.__data__.id==='G3P_regeneration').length};
        }''')
        self.assertEqual({'xu':True, 'ru':'Ru5P', 'g3p':True, 'sharedG3p':1}, actual)

    def test_dragging_shared_g3p_updates_both_reactions_and_survives_reload(self):
        before = self.snapshot()
        node = self.page.locator('circle.metabolite#G3P_regeneration')
        box = node.bounding_box()
        self.assertIsNotNone(box)
        x, y = box['x'] + box['width']/2, box['y'] + box['height']/2
        self.page.mouse.move(x, y)
        self.page.mouse.down()
        self.page.mouse.move(x + 16, y - 8, steps=8)
        self.page.mouse.up()
        self.wait_for_scene()
        moved = self.snapshot()
        entries = [r['notes']['map_info']['local_instances']['G3P'] for r in moved['reactions']
                   if r['id'] in ('IML_TALA', 'IML_TKT1')]
        self.assertEqual(entries[0], entries[1])
        self.assertGreater(abs(entries[0]['x'] - 475), 5)
        self.assertEqual(next(m for m in before['metabolites'] if m['id']=='G3P'),
                         next(m for m in moved['metabolites'] if m['id']=='G3P'))
        with self.page.expect_download() as download:
            self.page.locator('#freeflux_map_options .download').click()
        saved = self.directory / 'moved_shared_map.json'
        download.value.save_as(saved)
        geometry = self.geometry()
        self.load_saved_map(saved)
        self.assertEqual(moved, self.snapshot())
        self.assertEqual(geometry, self.geometry())
        self.assertTrue(self.page.evaluate(self.auditor)['pass'])

    def test_old_spiderweb_is_not_misreported_as_zero_crossings(self):
        old = DEMO / 'candidate01_hyperedge_v4_1_0_final.html'
        self.assertTrue(old.exists())
        self.page.goto(old.as_uri())
        self.page.wait_for_function("() => [...document.querySelectorAll('path.link')].length===21 && [...document.querySelectorAll('path.link')].every(p=>p.getAttribute('d'))")
        qa = self.page.evaluate(self.auditor)
        self.assertGreater(len(qa['crossings']), 0, qa)


if __name__ == '__main__':
    unittest.main()
