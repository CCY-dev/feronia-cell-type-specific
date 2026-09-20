#!/usr/bin/env python3
# Prototype; not yet carefully vetted.
"""Name layers of the code-generated evidence bars without changing geometry."""
from pathlib import Path
import xml.etree.ElementTree as ET
NS='http://www.w3.org/2000/svg'
ET.register_namespace('',NS)
ET.register_namespace('xlink','http://www.w3.org/1999/xlink')
path=Path(__file__).with_name('dof_evidence_bars.svg')
tree=ET.parse(path);root=tree.getroot()
for i,g in enumerate(root.findall(f'{{{NS}}}g')):
    g.set('id',f'evidence-bars-region-{i}')
    run=[];last=None;group=None
    for e in list(g):
        tag=e.tag.rsplit('}',1)[-1]
        if tag in ('defs','clipPath'):continue
        txt=''.join(e.itertext())
        cat='counts-and-axis'
        if tag=='rect':cat='evidence-bars-and-key'
        if tag=='text':
            if '37 DOF' in txt:cat='title'
            elif 'call' in txt:cat='evidence-key'
            elif 'GRN' in txt or 'nominated' in txt:cat='nomination-row-labels'
        if cat!=last:
            group=ET.Element(f'{{{NS}}}g',{'id':f'{cat}-{i}-{len(run)}'})
            run.append(group);last=cat
        g.remove(e);group.append(e)
    for sub in run:g.append(sub)
tree.write(path,encoding='utf-8',xml_declaration=True)
