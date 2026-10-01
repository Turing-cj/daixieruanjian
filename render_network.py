#!/usr/bin/env python3
"""Render a solved metabolic network as the interactive FreeFlux Studio HTML."""
import sys
from pathlib import Path

from scripts_render_d3flux import main


if __name__ == '__main__':
    if '--initial-map' not in sys.argv:
        sys.argv.extend(['--initial-map', str(Path(__file__).with_name('initial_template.json'))])
    main()
