#!/usr/bin/env python3
"""Wrap a d3flux/Jupyter HTML fragment with local browser dependencies."""
import argparse
import os
from pathlib import Path

FILES = ('require.min.js', 'd3.v3.min.js', 'jquery-2.2.4.min.js',
         'cola.v3.min.js', 'math.min.js', 'FileSaver.min.js')

VENDOR = Path(__file__).resolve().parent / 'vendor'


def write_standalone(text, output, vendor=VENDOR):
    replacements = {
        "'https://code.jquery.com/jquery-2.2.4.min'": "'vendor/jquery-2.2.4.min'",
        '"https://d3js.org/d3.v3.min"': '"vendor/d3.v3.min"',
        '"https://cdn.jsdelivr.net/gh/tgdwyer/WebCola@v3/WebCola/cola.v3.min"': '"vendor/cola.v3.min"',
        '"https://cdnjs.cloudflare.com/ajax/libs/mathjs/2.4.0/math.min"': '"vendor/math.min"',
        '"https://cdnjs.cloudflare.com/ajax/libs/FileSaver.js/2014-11-29/FileSaver.min"': '"vendor/FileSaver.min"',
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    text = text.replace('JSON.stringify(model);', 'JSON.stringify(model, null, 2) + "\\n";')
    # WebCola v3 expects the global d3 object before its shim executes.
    text = text.replace("exports: 'cola'", "deps: ['d3'], exports: 'cola'")
    # Packing disconnected components moves even fixed template coordinates.
    text = text.replace('.links(links)', '.links(links).handleDisconnected(!nodes.every(function(n) { return n.fixed; }))')
    page = '<!doctype html><html><head><meta charset="utf-8"><title>d3flux map</title>' \
        '<meta name="viewport" content="width=device-width,initial-scale=1">' \
        '<style>body{margin:0;background:white;font-family:Arial,sans-serif}button{margin:6px;padding:6px 12px;cursor:pointer}</style>' \
        '<script src="vendor/require.min.js"></script></head><body>' + text + '''
<script>
// Replace the example's obsolete external SVG Crowbar download, not its renderer.
document.addEventListener('click', function(event) {
  if (!event.target.closest('.svgbutton')) return;
  event.preventDefault(); event.stopImmediatePropagation();
  var svg = document.querySelector('svg').cloneNode(true);
  svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  var box = svg.viewBox.baseVal;
  svg.setAttribute('width', box.width); svg.setAttribute('height', box.height);
  svg.removeAttribute('style');
  var url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(svg)], {type:'image/svg+xml'}));
  var link = document.createElement('a'); link.href = url; link.download = 'flux-map.svg'; link.click();
  setTimeout(function() { URL.revokeObjectURL(url); }, 1000);
}, true);
</script></body></html>'''
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(page, encoding='utf-8')
    target = output.parent / 'vendor'
    target.mkdir(exist_ok=True)
    for name in FILES:
        if (target / name).resolve() != (vendor / name).resolve():
            (target / name).write_bytes((vendor / name).read_bytes())


def export_native(html, size, vendor=VENDOR):
    from playwright.sync_api import sync_playwright

    environment = dict(os.environ)
    local_libs = vendor / 'browser-libs/usr/lib/x86_64-linux-gnu'
    if local_libs.exists():
        environment['LD_LIBRARY_PATH'] = str(local_libs) + ':' + environment.get('LD_LIBRARY_PATH', '')
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(env=environment)
        page = browser.new_page(viewport={'width': size[0], 'height': size[1] + 80}, device_scale_factor=2)
        errors = []
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.goto(html.resolve().as_uri())
        page.wait_for_function('''() => {
          const paths = [...document.querySelectorAll('path.link')];
          if (!paths.length || paths.some(p => !p.getAttribute('d') || /NaN|Infinity/.test(p.getAttribute('d')))) return false;
          const geometry = paths.map(p => p.getAttribute('d')).join('');
          if (window.lastGeometry !== geometry) { window.lastGeometry = geometry; window.stableSince = performance.now(); }
          return performance.now() - window.stableSince > 500;
        }''', timeout=30000)
        if errors:
            raise RuntimeError('\n'.join(errors))
        svg = page.locator('svg').first
        markup = svg.evaluate('''svg => {
          const copy = svg.cloneNode(true), box = svg.viewBox.baseVal;
          copy.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
          copy.setAttribute('width', box.width); copy.setAttribute('height', box.height);
          copy.removeAttribute('style');
          return new XMLSerializer().serializeToString(copy);
        }''')
        html.with_suffix('.svg').write_text(markup, encoding='utf-8')
        # Export at the SVG's declared size instead of the sidebar-reduced size.
        raster = browser.new_page(viewport={'width':size[0], 'height':size[1]}, device_scale_factor=2)
        raster.set_content('<style>body{margin:0;background:white}svg{display:block}</style>' + markup)
        raster.locator('svg').screenshot(path=str(html.with_suffix('.png')))
        raster.close()
        geometry_qa = page.evaluate('''() => {
          const paths = [...document.querySelectorAll('path.link')];
          const sampleStep = 8;
          function meta(path) {
            const d = path.__data__ || {}, rxn = d.rxn || {};
            return {
              reaction: rxn.id || null,
              source: d.source && d.source.id || null,
              target: d.target && d.target.id || null
            };
          }
          function samples(path) {
            const length = path.getTotalLength();
            const count = Math.max(2, Math.ceil(length / sampleStep));
            const points = [];
            for (let i = 0; i <= count; i++) {
              const p = path.getPointAtLength(length * i / count);
              points.push({x:p.x, y:p.y});
            }
            return points;
          }
          function orient(a,b,c) {
            return (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x);
          }
          function properIntersect(a,b,c,d) {
            const o1=orient(a,b,c), o2=orient(a,b,d);
            const o3=orient(c,d,a), o4=orient(c,d,b);
            return o1*o2 < 0 && o3*o4 < 0;
          }
          function pointSegmentDistance(p,a,b) {
            const vx=b.x-a.x, vy=b.y-a.y, wx=p.x-a.x, wy=p.y-a.y;
            const l2=vx*vx+vy*vy;
            const t=l2 ? Math.max(0,Math.min(1,(wx*vx+wy*vy)/l2)) : 0;
            const qx=a.x+t*vx, qy=a.y+t*vy;
            return Math.hypot(p.x-qx,p.y-qy);
          }

          const sampled = paths.map(path => ({path, meta:meta(path), points:samples(path)}));
          const crossingPairs = [];
          for (let i=0; i<sampled.length; i++) {
            for (let j=i+1; j<sampled.length; j++) {
              const a=sampled[i], b=sampled[j];
              if (a.meta.reaction && a.meta.reaction===b.meta.reaction) continue;
              const endpointsA=new Set([a.meta.source,a.meta.target].filter(Boolean));
              if ([b.meta.source,b.meta.target].some(id => id && endpointsA.has(id))) continue;
              let hit=false;
              for (let x=0; x<a.points.length-1 && !hit; x++) {
                for (let y=0; y<b.points.length-1; y++) {
                  if (properIntersect(a.points[x],a.points[x+1],
                                      b.points[y],b.points[y+1])) {
                    hit=true; break;
                  }
                }
              }
              if (hit) crossingPairs.push({
                first: a.meta.reaction+':'+a.meta.source+'->'+a.meta.target,
                second: b.meta.reaction+':'+b.meta.source+'->'+b.meta.target
              });
            }
          }

          const nodes = [...document.querySelectorAll('g.node-group')].map(group => {
            const d=group.__data__ || {}, circle=group.querySelector('circle');
            return {
              id:d.id || null,
              x:+d.x, y:+d.y,
              radius:circle ? (+circle.getAttribute('r') || 5) : 5
            };
          }).filter(node => node.id && Number.isFinite(node.x) && Number.isFinite(node.y));
          const obstacleHits=[];
          sampled.forEach(item => {
            const excluded=new Set([item.meta.reaction,item.meta.source,item.meta.target].filter(Boolean));
            nodes.forEach(node => {
              if (excluded.has(node.id)) return;
              const limit=node.radius+12;
              let hit=false;
              for (let i=0;i<item.points.length-1;i++) {
                if (pointSegmentDistance(node,item.points[i],item.points[i+1]) < limit) {
                  hit=true; break;
                }
              }
              if (hit) obstacleHits.push({
                edge:item.meta.reaction+':'+item.meta.source+'->'+item.meta.target,
                node:node.id
              });
            });
          });
          return {
            sample_step_px: sampleStep,
            sampled_path_crossings: crossingPairs.length,
            crossing_pairs: crossingPairs,
            sampled_obstacle_hits: obstacleHits.length,
            obstacle_hits: obstacleHits
          };
        }''')
        counts = {
            'metabolite_circles': page.locator('circle.metabolite').count(),
            'reaction_paths': page.locator('path.link').count(),
            'cubic_paths': page.locator('path.link').evaluate_all(
                "ps => ps.filter(p => p.getAttribute('d').includes('C')).length"),
            'quadratic_paths': page.locator('path.link').evaluate_all(
                "ps => ps.filter(p => p.getAttribute('d').includes('Q')).length"),
            'linear_paths': page.locator('path.link').evaluate_all(
                "ps => ps.filter(p => p.getAttribute('d').includes('L')).length"),
            'geometry_qa': geometry_qa,
        }
        auditor = Path(__file__).with_name('browser_circle_qa.js')
        if auditor.exists():
            counts['strict_visual_qa'] = page.evaluate(auditor.read_text(encoding='utf-8'))
        browser.close()
    return counts


def main():
    p = argparse.ArgumentParser()
    p.add_argument('fragment', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--vendor', type=Path, default=VENDOR)
    a = p.parse_args()
    write_standalone(a.fragment.read_text(encoding='utf-8'), a.output, a.vendor)
    print(a.output.resolve())

if __name__ == '__main__': main()
