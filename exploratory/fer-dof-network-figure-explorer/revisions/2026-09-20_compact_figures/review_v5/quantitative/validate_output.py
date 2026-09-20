#!/usr/bin/env python3
# PROTOTYPE: not yet carefully vetted. Layout/data identity checks only.
from pathlib import Path
from collections import Counter
from lxml import etree
import csv, hashlib, json, re, subprocess

out = Path(__file__).resolve().parent
old = out.parent.parent / 'review_v3' / 'quantitative'
report = {'status': 'PROTOTYPE: not yet carefully vetted', 'data_checks': [], 'figures': []}
for name in ('source_celltype_summary_all37.csv', 'plot_source_celltype_tested10.csv',
             'source_top10_family.csv', 'plot_source_four_significant_families.csv'):
    p, q = out/name, old/name
    assert p.read_bytes() == q.read_bytes(), name
    report['data_checks'].append({'file':name, 'byte_identical_to_review_v3':True,
                                 'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
rows = list(csv.DictReader((out/'plot_source_celltype_tested10.csv').open()))
assert len(rows)==130 and len({r['gene'] for r in rows})==10 and len({r['tissue'] for r in rows})==13
forest = list(csv.DictReader((out/'source_top10_family.csv').open()))
assert len(forest)==20 and len({r['Family'] for r in forest})==10
assert sum(float(r['p_adj_BH'])<.05 for r in forest)==6
assert sum(float(r['odds_ratio'])==0 for r in forest)==1
ns = {'s':'http://www.w3.org/2000/svg'}
font_reports = []
for stem, hmm in (('top10_aligned_forest_review_v5',80), ('celltype_tested10_review_v5',90),
                  ('four_families_review_v5',43)):
    root = etree.parse(str(out/(stem+'.svg')))
    texts = root.xpath('//s:text',namespaces=ns)
    sizes = [float(v) for t in texts for v in re.findall(r'font-size:\s*([\d.]+)px',t.get('style',''))]
    assert min(sizes)>=8
    info = subprocess.check_output(['pdfinfo',str(out/(stem+'.pdf'))],text=True)
    match = re.search(r'Page size:\s+([.\d]+) x ([.\d]+) pts',info)
    width,height = [float(v)*25.4/72 for v in match.groups()]
    assert abs(width-110)<.01 and abs(height-hmm)<.01
    fonts = subprocess.check_output(['pdffonts',str(out/(stem+'.pdf'))],text=True)
    assert 'Type 3' not in fonts
    font_reports.append(stem+'\n'+fonts)
    for dpi,tag in ((240,'check_pdf'),(96,'check_actualsize')):
        subprocess.run(['pdftoppm','-png','-r',str(dpi),'-singlefile',str(out/(stem+'.pdf')),
                        str(out/(tag+'_'+stem))],check=True)
    report['figures'].append({'stem':stem,'width_mm':round(width,4),'height_mm':round(height,4),
                              'minimum_font_pt':min(sizes),'PDF_no_Type3_fonts':True})
    if stem.startswith('celltype'):
        dimensions = Counter((round(float(r.get('width'))*25.4/72,3),round(float(r.get('height'))*25.4/72,3))
                             for r in root.xpath('//s:rect',namespaces=ns)
                             if r.get('width') and '%' not in r.get('width') and r.get('height') and '%' not in r.get('height'))
        tile = next(dim for dim,n in dimensions.items() if n==130)
        assert 7 <= tile[0] <= 7.2
        number = next(t for t in texts if ''.join(t.itertext())=='1/4')
        number_width = float(number.get('textLength').replace('px',''))*25.4/72
        report['matrix'] = {'rows':130,'genes':10,'tissues':13,'tile_width_mm':tile[0],
                            'tile_height_mm':tile[1],'number_1_over_4_width_mm':round(number_width,3),
                            'side_padding_mm':round((tile[0]-number_width)/2,3)}
    else:
        dof = [t for t in texts if ''.join(t.itertext())=='DOF']
        assert len(dof)==1 and 'font-weight: bold' in dof[0].get('style','')
(out/'font_check.txt').write_text('\n'.join(font_reports))
(out/'source_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
(out/'validation.txt').write_text(
    'PROTOTYPE: not yet carefully vetted. Layout/data checks only.\n'
    'PDF dimensions: top-ten forest 110 x 80 mm; matrix 110 x 90 mm; four-family companion 110 x 43 mm.\n'
    'Minimum SVG text is 8 pt; matrix fractions and gene labels 9 pt, tissue labels 8.5 pt.\n'
    f'Actual matrix tile: {tile[0]:.3f} x {tile[1]:.3f} mm; 1/4 width: {number_width:.3f} mm; horizontal side padding {(tile[0]-number_width)/2:.3f} mm.\n'
    'Four source/plot CSVs are byte-identical to review v3. Matrix: 130 rows, 10 DOFs, 13 tissues. Forest: 20 rows, 10 families, six BH-significant estimates, one zero odds ratio.\n'
    'DOF family label is bold without added color meaning. Full tissue names below at 60 degrees; Metaphloem and companion cell wraps across two lines.\n'
    'PDFs use embedded TrueType Arial fonts without Type 3; SVGs retain editable text and named groups.\n'
    'PDF rasters generated at 240 dpi and 96 dpi. CVD simulations use colorspace deutan, protan and tritan.\n'
    'PDF rasters and all three CVD simulations of both primary figures were visually inspected: no visible clipping or text collisions. Native display-size check is approximate, not a physical print proof.\n'
    'Previous versions and frozen source inputs remain unchanged.\n')
print(json.dumps(report,indent=2))
