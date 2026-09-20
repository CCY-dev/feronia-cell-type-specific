#!/usr/bin/env python3
# PROTOTYPE: not yet carefully vetted. Checks cover display identity and geometry only.
from pathlib import Path
from collections import Counter
from lxml import etree
import csv, hashlib, json, re, subprocess

out=Path(__file__).resolve().parent
src=out.parent.parent/'review_v2'/'all_contexts'/'all_contexts_as_exported_union_edges.csv'
copy=out/'source_union_edges.csv'
assert src.read_bytes()==copy.read_bytes()
original=list(csv.DictReader(src.open()))
plotted=list(csv.DictReader((out/'plotted_edges.csv').open()))
core=list(original[0])
key=lambda r:(r['source'],r['target'])
assert {key(r):r for r in original}=={key(r):{k:r[k] for k in core} for r in plotted}
tfs=list(csv.DictReader((out/'plotted_tf_nodes.csv').open()))
targets=list(csv.DictReader((out/'plotted_target_nodes.csv').open()))
assert len(plotted)==180 and len(tfs)==10 and len(targets)==37
assert min(float(r['y']) for r in tfs)>max(float(r['y']) for r in targets)
counts=Counter(r['support_class'] for r in original)
assert counts=={'supported_in_current_seven':147,'additional_excluded_contexts_only':33}

svg=out/'all_contexts_panel_d.svg'
root=etree.parse(str(svg));ns={'s':'http://www.w3.org/2000/svg'}
tf_names={r['source_symbol'] for r in tfs};target_names={r['target_name'] for r in targets}
group_counts=Counter()
for j,g in enumerate(root.xpath('/s:svg/s:g',namespaces=ns)):
    g.set('id',f'panel-d-generated-region-{j+1}')
    last=None;sub=None
    for child in list(g):
        tag=etree.QName(child).localname
        txt=''.join(child.itertext());style=child.get('style','')
        if tag=='line':
            semantic=('network-' if 'stroke-opacity: 0.50' in style else 'support-key-')+('additional-pairs' if 'stroke-dasharray' in style else 'seven-fit-pairs')
        elif tag=='circle':
            semantic='dof-nodes' if '#253B49' in style else 'target-family-nodes'
        elif tag=='text':
            semantic='dof-labels' if txt in tf_names else 'target-labels' if txt in target_names else 'support-key-labels' if 'pairs:' in txt else 'family-key-labels'
        else:semantic='layout'
        if semantic!=last:
            group_counts[semantic]+=1
            sub=etree.Element('{%s}g'%ns['s'],id=f'panel-d-{semantic}-{group_counts[semantic]}')
            g.insert(g.index(child),sub);last=semantic
        g.remove(child);sub.append(child)
root.write(str(svg),encoding='UTF-8',xml_declaration=True)
texts=root.xpath('//s:text',namespaces=ns)
sizes=[float(v) for t in texts for v in re.findall(r'font-size:\s*([.\d]+)px',t.get('style',''))]
assert min(sizes)>=8.4
assert sum(''.join(t.itertext()) in tf_names for t in texts)==10
assert sum(''.join(t.itertext()) in target_names for t in texts)==37
lines=root.xpath('//s:line',namespaces=ns)
network=[l for l in lines if 'stroke-opacity: 0.50' in l.get('style','')]
assert len(network)==180 and sum('stroke-dasharray' in l.get('style','') for l in network)==33
assert all('stroke-width: 0.32;' in l.get('style','') for l in network)
pdf=out/'all_contexts_panel_d.pdf'
info=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
assert re.search(r'Pages:\s+1(?:\s|$)',info)
wh=[float(v)*25.4/72 for v in re.search(r'Page size:\s+([.\d]+) x ([.\d]+) pts',info).groups()]
assert abs(wh[0]-135)<.01 and abs(wh[1]-90)<.01
fonts=subprocess.check_output(['pdffonts',str(pdf)],text=True)
assert 'Type 3' not in fonts
(out/'font_check.txt').write_text(fonts)
for dpi,name in [(600,'check_pdf_600dpi'),(96,'check_native_size')]:
    subprocess.run(['pdftoppm','-png','-r',str(dpi),'-singlefile',str(pdf),str(out/name)],check=True)
report={'status':'PROTOTYPE: not yet carefully vetted','source_path':str(src),
        'source_sha256':hashlib.sha256(copy.read_bytes()).hexdigest(),
        'source_byte_identical':True,'all_original_edge_columns_identical':True,
        'n_pairs':180,'n_dofs':10,'n_targets':37,'support_classes':dict(counts),
        'pdf_width_mm':round(wh[0],3),'pdf_height_mm':round(wh[1],3),'pdf_pages':1,
        'smallest_text_pt':min(sizes),'network_stroke_pt':.32,'network_alpha':.5,
        'dof_point_diameter_mm':round(2*4.23*25.4/72,3),
        'target_point_diameter_mm':round(2*3.02*25.4/72,3),
        'all_dofs_above_all_targets':True,'editable_text_elements':len(texts),
        'named_groups':dict(group_counts)}
(out/'source_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
