#!/usr/bin/env python3
# Prototype; not yet carefully vetted.
"""Group generated SVG for editing without changing plotted geometry."""
import sys
import xml.etree.ElementTree as ET
NS='http://www.w3.org/2000/svg'
ET.register_namespace('',NS)
ET.register_namespace('xlink','http://www.w3.org/1999/xlink')
path=sys.argv[1]
tree=ET.parse(path);root=tree.getroot()
panel=next(g for g in root.findall(f'{{{NS}}}g') if g.find(f'{{{NS}}}circle') is not None)
elements=list(panel)
header=next(e for e in elements if e.tag.endswith('text') and e.text=='Exclusive overlap')
bar_x=float(header.attrib['x'])
groups={name:ET.Element(f'{{{NS}}}g',{'id':name}) for name in ('row-backgrounds','membership-and-set-totals','overlap-bars-and-axis')}
for e in elements:
    panel.remove(e)
    x=float(e.attrib.get('x',e.attrib.get('x1',e.attrib.get('cx',0))))
    if '#F6F7F8' in e.attrib.get('style',''):
        label='row-backgrounds'
    elif x >= bar_x-0.1:
        label='overlap-bars-and-axis'
    else:
        label='membership-and-set-totals'
    groups[label].append(e)
for g in groups.values():panel.append(g)
tree.write(path,encoding='utf-8',xml_declaration=True)
