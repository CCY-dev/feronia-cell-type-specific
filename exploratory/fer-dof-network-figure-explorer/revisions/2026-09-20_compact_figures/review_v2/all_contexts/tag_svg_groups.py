#!/usr/bin/env python3
# Prototype; not yet carefully vetted.
"""Add semantic group IDs to existing svglite output without changing geometry."""
from pathlib import Path
import csv
import xml.etree.ElementTree as ET
ns = 'http://www.w3.org/2000/svg'
ET.register_namespace('', ns)
ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')
base = Path(__file__).resolve().parent
for path in base.glob('*_compact.svg'):
    prefix = path.stem
    with (base / (prefix + '_tf_nodes.csv')).open() as h:
        tfs = {r['source_symbol'] for r in csv.DictReader(h)}
    with (base / (prefix + '_target_nodes.csv')).open() as h:
        targets = {r['target_name'] for r in csv.DictReader(h)}
    tree = ET.parse(path)
    root = tree.getroot()
    counter = {}
    for i, group in enumerate(list(root.findall(f'{{{ns}}}g')), 1):
        group.set('id', f'plot-region-{i}')
        children = list(group)
        for child in children:
            group.remove(child)
        last_key = None
        subgroup = None
        for child in children:
            tag = child.tag.rsplit('}', 1)[-1]
            text = ''.join(child.itertext())
            if tag == 'line': key = 'additional-pairs-from-excluded-fits' if 'stroke-dasharray' in child.get('style', '') else 'pairs-supported-in-seven-fits'
            elif tag == 'circle': key = 'dof-nodes' if '#253B49' in child.get('style', '') else 'target-family-nodes'
            elif tag == 'text':
                key = 'dof-labels' if text in tfs else 'target-labels' if text in targets else 'title' if text.startswith('DOF links') else 'family-legend'
            else: key = 'plot-layout'
            if key != last_key:
                counter[key] = counter.get(key, 0) + 1
                subgroup = ET.SubElement(group, f'{{{ns}}}g', {'id': f'{key}-{counter[key]}'})
                last_key = key
            subgroup.append(child)
    tree.write(path, encoding='utf-8', xml_declaration=True)
