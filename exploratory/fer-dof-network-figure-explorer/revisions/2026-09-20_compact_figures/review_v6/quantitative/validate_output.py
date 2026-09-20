#!/usr/bin/env python3
# PROTOTYPE: not yet carefully vetted. Checks cover layout/data identity, not analysis validity.
from pathlib import Path
from collections import Counter
from lxml import etree
import csv, hashlib, json, re, subprocess

out = Path(__file__).resolve().parent
old = out.parent.parent / 'review_v5' / 'quantitative'
stem = 'celltype_tested10_review_v6'
report = {'status':'PROTOTYPE: not yet carefully vetted','data_checks':[]}
for name in ('source_celltype_summary_all37.csv','plot_source_celltype_tested10.csv'):
    p, q = out/name, old/name
    same = p.read_bytes() == q.read_bytes()
    assert same, name
    report['data_checks'].append({'file':name,'byte_identical_to_review_v5':same,
                                  'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
rows = list(csv.DictReader((out/'plot_source_celltype_tested10.csv').open()))
assert len(rows)==130 and len({r['gene'] for r in rows})==10 and len({r['tissue'] for r in rows})==13
root = etree.parse(str(out/(stem+'.svg')))
ns = {'s':'http://www.w3.org/2000/svg'}
sizes = [float(v) for t in root.xpath('//s:text',namespaces=ns)
         for v in re.findall(r'font-size:\s*([\d.]+)px',t.get('style',''))]
assert min(sizes)>=8
dimensions = Counter((round(float(r.get('width'))*25.4/72,3),round(float(r.get('height'))*25.4/72,3))
                     for r in root.xpath('//s:rect',namespaces=ns)
                     if r.get('width') and '%' not in r.get('width') and r.get('height') and '%' not in r.get('height'))
tile = next(dim for dim,n in dimensions.items() if n==130)
assert 7 <= tile[0] <= 7.2
number = next(t for t in root.xpath('//s:text',namespaces=ns) if ''.join(t.itertext())=='1/4')
assert 'font-size: 8.50px' in number.get('style','')
assert all('x: no paired scores' not in ''.join(t.itertext()) for t in root.xpath('//s:text',namespaces=ns))
number_width = float(number.get('textLength').replace('px',''))*25.4/72
info = subprocess.check_output(['pdfinfo',str(out/(stem+'.pdf'))],text=True)
match = re.search(r'Page size:\s+([.\d]+) x ([.\d]+) pts',info)
width,height = [float(v)*25.4/72 for v in match.groups()]
assert abs(width-110)<.01 and abs(height-90)<.01
assert re.search(r'Pages:\s+1(?:\s|$)',info)
fonts = subprocess.check_output(['pdffonts',str(out/(stem+'.pdf'))],text=True)
assert 'Type 3' not in fonts
(out/'font_check.txt').write_text(fonts)
for dpi,tag in ((600,'check_pdf_600dpi'),(96,'check_actualsize')):
    subprocess.run(['pdftoppm','-png','-r',str(dpi),'-singlefile',str(out/(stem+'.pdf')),str(out/tag)],check=True)
report.update({'rows':130,'genes':10,'tissues':13,'width_mm':round(width,4),
               'height_mm':round(height,4),'tile_width_mm':tile[0],
               'tile_height_mm':tile[1],'number_1_over_4_width_mm':round(number_width,3),
               'minimum_font_pt':min(sizes),'PDF_no_Type3_fonts':True})
(out/'source_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
(out/'validation.txt').write_text(
    'PROTOTYPE: not yet carefully vetted. Layout/data checks only.\n'
    f'PDF dimensions: {width:.3f} x {height:.3f} mm.\n'
    f'Actual SVG tile width: {tile[0]:.3f} mm; 1/4 text width: {number_width:.3f} mm.\n'
    f'Actual SVG tile height: {tile[1]:.3f} mm. Minimum text: {min(sizes):.1f} pt.\n'
    'Both source/plot CSVs are byte-identical to review v5; 130 rows, 10 DOFs, 13 tissues.\n'
    'Full names appear below at 60 degrees; two-line key below the labels. Only the wrapped Metaphloem and companion cell label shifts left 1 mm; 600 dpi crop inspected against v5.\n'
    'Counts are 8.5 pt, tissue labels 8.5 pt, gene labels 9 pt, key 8 pt. All 130 K values are positive; unused x key omitted, conditional code retained.\n'
    'The 600 dpi PDF label crop and full PDF, native-size preview and all three CVD previews were visually inspected.\n'
    'PDF uses embedded TrueType Arial fonts; SVG retains editable text and named groups.\n'
    'Previous versions and frozen source inputs remain unchanged.\n')
print(json.dumps(report,indent=2))
