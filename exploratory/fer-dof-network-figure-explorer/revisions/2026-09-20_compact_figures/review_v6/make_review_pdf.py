#!/usr/bin/env python3
"""Prototype; not yet carefully vetted. Native-size figure review v6."""
# Layout checks do not establish scientific validation. Source reconciliation remains pending.
from pathlib import Path
from io import BytesIO
from reportlab.pdfgen import canvas
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph
from reportlab.lib.styles import ParagraphStyle
from pypdf import PdfReader, PdfWriter, Transformation

ROOT=Path(__file__).resolve().parent
FONT=Path('/System/Library/Fonts/Supplemental')
pdfmetrics.registerFont(TTFont('Arial',str(FONT/'Arial.ttf')))
pdfmetrics.registerFont(TTFont('Arial-Bold',str(FONT/'Arial Bold.ttf')))
STYLE=ParagraphStyle('body',fontName='Arial',fontSize=8.5,leading=11,spaceAfter=4)

def page(width,height,title,entries,notes=(),footer=None,show_dimensions=True):
    buf=BytesIO();c=canvas.Canvas(buf,pagesize=(width*mm,height*mm))
    c.setFont('Arial-Bold',12);c.drawString(12*mm,(height-11)*mm,title)
    c.setFont('Arial',8);c.drawString(12*mm,(height-17)*mm,'Prototype; not yet carefully vetted. Print at 100% / actual size.')
    for rel,x,y,label in entries:
        src=PdfReader(ROOT/rel).pages[0];w=float(src.mediabox.width);h=float(src.mediabox.height)
        assert x*mm+w<(width-7)*mm and y*mm+h<(height-21)*mm
        c.setFont('Arial-Bold',9);c.drawString(x*mm,y*mm+h+3*mm,label)
        if show_dimensions:
            c.setFont('Arial',7);c.drawString(x*mm,y*mm-3*mm,f'{w/mm:.0f} x {h/mm:.0f} mm; native vector panel')
    for text,x,ytop,w in notes:
        par=Paragraph(text,STYLE);_,h=par.wrap(w*mm,height*mm)
        par.drawOn(c,x*mm,ytop*mm-h)
    if footer:
        c.setFont('Arial',7.5);c.drawString(12*mm,7*mm,footer)
    c.showPage();c.save();buf.seek(0);dst=PdfReader(buf).pages[0]
    for rel,x,y,_ in entries:
        dst.merge_transformed_page(PdfReader(ROOT/rel).pages[0],Transformation().translate(x*mm,y*mm))
    return dst

w=PdfWriter()
w.add_page(page(297,210,'Figure candidates: family enrichment, DOF evidence, cellular context',[
    ('../review_v5/quantitative/top10_aligned_forest_review_v5.pdf',12,104,'a  TF family enrichment'),
    ('overlap/dof_evidence_bars.pdf',138,138,'b  Evidence within the DOF family'),
    ('quantitative/celltype_tested10_review_v6.pdf',12,8,'c  The ten assayed DOFs by cell type'),
],notes=[('Panel b uses all 37 DOFs. Panel a uses the 33 DOFs eligible for its paired-score TF background. These denominators answer different questions.',139,132,108)],
footer=None,show_dimensions=False))

net=ROOT/'../review_v2/all_contexts/all_contexts_overview_compact.pdf'
if net.exists():
    w.add_page(page(210,247,'All-context overview and a vascular example',[
        ('../review_v2/all_contexts/all_contexts_overview_compact.pdf',12,121,'d  All available cellular contexts considered'),
        ('../networks/vascular_overview_compact.pdf',12,20,'Supporting example: elongating pericycle and procambium'),
    ],footer='Basal-network unions summarize retained links across contexts; the lines do not represent intercellular transport.'))

w.add_page(page(210,170,'Retained UpSet option and data provenance',[
    ('../upset/compact_upset.pdf',12,83,'Alternative b  Exclusive intersections'),
],notes=[
    ('The inclusive DOF-DEG overlap is <b>12</b>. All 12 are also GRN-nominated.<br/><br/>The <b>0</b> row means DOF + DEG <b>without GRN nomination</b>. Empty dots denote exclusion, not an ignored set.<br/><br/>Singletons omitted from this display: 11 DOFs, 314 GRN-TFs and 5,282 scDEGs.',108,136,89),
    ('<b>Underlying data are available.</b> Live checks on 20 September 2026 verified the integrated scRNA-seq object, CellOracle inputs and saved links objects, and full DE tables on Resnick. The <a href="https://github.com/tmnolan/feronia-cell-type-specific/tree/5649a9485de0d50d03c8135c5ff254bdaa8d0a1a" color="#0072B2">paper GitHub repository</a> provides source code and exports; its HEAD matched the pinned revision at the recorded access check.',12,72,184),
    ('<b>The limitation is consistency, not missing data.</b> Three saved fer-4 endodermal graphs disagree with their stored degree scores. Separate phloem/xylem target selection has a substring-matching issue. The figure revision does not repair either finding. The independent code and visualization audit is complete and confirms both source issues; neither has been repaired. See the linked <a href="file:///Users/tmn23/projects/fer-sc-manuscript-2026/revisions/2026-09-20_compact_figures/review_v2/SOURCE_AVAILABILITY.md" color="#0072B2">source record</a> for exact files and checks.',12,43,184),
],footer='No model refit or source repair. Independent checks reproduce the saved-table calculations.'))
w.add_metadata({'/Title':'FER figure review v6: DOF emphasis and full cell-type labels; prototype','/Author':'Nolan Lab'})
out=ROOT/'review_at_print_size_v6.pdf'
with out.open('wb') as f:w.write(f)
print(out)
