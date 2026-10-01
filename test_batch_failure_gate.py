"""A failed render must not be reported successful using a stale manifest."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from nine_paths_pdf_20260921 import build_batch as batch


class BatchFailureGateTests(unittest.TestCase):
    def test_stale_success_manifest_cannot_mask_a_new_render_exception(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            folder = root / 'candidate_01'
            folder.mkdir()
            (folder / 'layout.json').write_text(json.dumps({'publication_style': {}}), encoding='utf-8')
            (folder / 'candidate_01.manifest.json').write_text(json.dumps({
                'acceptance': {'pass': True}, 'rendered_elements': {'strict_visual_qa': {'pass': True}}}), encoding='utf-8')
            with patch.object(batch, 'HERE', root), patch('scripts_render_d3flux.main', side_effect=RuntimeError('render failed')):
                with self.assertRaises(RuntimeError):
                    batch.render([1])


if __name__ == '__main__':
    unittest.main()
