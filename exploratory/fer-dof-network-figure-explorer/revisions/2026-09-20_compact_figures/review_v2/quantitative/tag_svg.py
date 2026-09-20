#!/usr/bin/env python3
# Prototype; not yet carefully vetted.
"""Name generated SVG regions for editing; never change plotted geometry or text."""
from pathlib import Path
from lxml import etree
import json, re

out = Path(__file__).resolve().parent
ns = {'s': 'http://www.w3.org/2000/svg'}
checks = []
for path in sorted(out.glob('*review_v2.svg')):
    root = etree.parse(str(path))
    width, height = [float(v) for v in root.getroot().get('viewBox').split()[2:]]
    used = {}
    # svglite sometimes places all exterior labels in one clipping group.
    # Split its direct text elements into named subgroups without changing order.
    for g in root.xpath('/s:svg/s:g', namespaces=ns):
        for t in list(g):
            if etree.QName(t).localname != 'text':
                continue
            label = ''.join(t.itertext())
            semantic = None
            if label.startswith('Upper:'):
                semantic = 'wt-key'
            elif label.startswith('Lower:'):
                semantic = 'fer4-key'
            elif label.startswith('Nominated / scored zones'):
                semantic = 'zone-denominator-key'
            elif label.startswith('Grey:'):
                semantic = 'missingness-key'
            elif label.startswith('Filled:'):
                semantic = 'statistical-key'
            elif label.startswith('Odds ratio'):
                semantic = 'shared-axis-label'
            if semantic:
                sub = etree.Element('{%s}g' % ns['s'])
                sub.set('id', path.stem+'-'+semantic)
                g.insert(g.index(t), sub)
                g.remove(t)
                sub.append(t)
    for i, g in enumerate(root.xpath('/s:svg/s:g', namespaces=ns)):
        words = g.xpath('.//s:text/text()', namespaces=ns)
        text = ' '.join(words)
        raw = etree.tostring(g).decode()
        name = 'generated-region'
        if 'celltype' in path.name:
            if len(g.xpath('.//s:rect', namespaces=ns)) > 20:
                name = 'celltype-tiles-and-zone-counts'
            elif 'DOF4.6' in text:
                name = 'dof-labels'
            elif 'Atri.' in text:
                name = 'celltype-labels'
            elif 'Upper:' in text:
                name = 'wt-key'
            elif 'Lower:' in text:
                name = 'fer4-key'
            elif 'n/K' in text:
                name = 'zone-denominator-key'
            elif 'Grey:' in text:
                name = 'missingness-key'
        else:
            if 'Odds ratio' in text:
                name = 'shared-axis-label'
            elif 'Filled:' in text:
                name = 'statistical-key'
            elif 'DOF' in words:
                name = 'family-labels'
            elif 'TF family' in text:
                name = 'family-heading'
            elif 'WT-associated' in text:
                name = 'wt-axis-and-heading'
            elif 'fer-4-associated' in text:
                name = 'fer4-axis-and-heading'
            elif len(g.xpath('.//s:circle', namespaces=ns)):
                name = 'wt-estimates-and-confidence-intervals' if '#0072B2' in raw else 'fer4-estimates-and-confidence-intervals'
            elif any(re.fullmatch(r'\d+/\d+', w) for w in words):
                xs = [float(t.get('x')) for t in g.xpath('.//s:text', namespaces=ns)]
                name = 'wt-counts' if sum(xs)/len(xs) < width*.7 else 'fer4-counts'
            elif text.startswith('Nominated /'):
                x = float(g.xpath('.//s:text', namespaces=ns)[0].get('x'))
                name = 'wt-count-heading' if x < width*.7 else 'fer4-count-heading'
        used[name] = used.get(name, 0)+1
        g.set('id', f'{path.stem}-{name}-{used[name]:02d}')
    root.write(str(path), encoding='UTF-8', xml_declaration=True)
    sizes = [float(v) for t in root.xpath('//s:text', namespaces=ns)
             for v in re.findall(r'font-size:\s*([\d.]+)px', t.get('style', ''))]
    assert sizes and min(sizes) >= 7
    checks.append({'file': path.name, 'width_mm': width/72*25.4,
                   'height_mm': height/72*25.4, 'smallest_text_pt': min(sizes),
                   'editable_text_elements': len(sizes), 'named_groups': used})
(out/'svg_validation.json').write_text(json.dumps(checks, indent=2)+'\n')
print(json.dumps(checks, indent=2))
