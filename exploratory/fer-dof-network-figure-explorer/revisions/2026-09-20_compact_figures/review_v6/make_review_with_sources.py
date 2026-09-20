#!/usr/bin/env python3
# Prototype; not yet carefully vetted.
"""Append linked input/output/code provenance to the current native-size proof."""
from pathlib import Path
from html import escape
import csv
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.lib.colors import HexColor
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader, PdfWriter

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
V3 = HERE.parent / 'review_v3'
V5 = HERE.parent / 'review_v5'
FONT = Path('/System/Library/Fonts/Supplemental')
pdfmetrics.registerFont(TTFont('Arial', str(FONT/'Arial.ttf')))
pdfmetrics.registerFont(TTFont('Arial-Bold', str(FONT/'Arial Bold.ttf')))
styles = {
    'title': ParagraphStyle('title', fontName='Arial-Bold', fontSize=14, leading=17, spaceAfter=8),
    'head': ParagraphStyle('head', fontName='Arial-Bold', fontSize=10, leading=13, spaceBefore=6, spaceAfter=5),
    'body': ParagraphStyle('body', fontName='Arial', fontSize=10, leading=13, spaceAfter=6),
    'small': ParagraphStyle('small', fontName='Arial', fontSize=9.5, leading=12.4, spaceAfter=5),
}

def para(text, style='body'):
    return Paragraph(text, styles[style])

def linked(path, label=None):
    path = Path(path)
    if not path.is_absolute():
        path = REPO/path
    assert path.exists(), path
    return f'<a href="{escape(path.resolve().as_uri(), quote=True)}" color="#0072B2">{escape(label or path.name)}</a>'

def paths(value):
    return [p.strip() for p in value.split(';') if p.strip()]

def listlinks(value, code=False):
    result = []
    for raw in paths(value):
        original = REPO/raw
        copy = V3/'code_reference'/raw
        target = copy if code and copy.exists() else original
        result.append(linked(target, original.name))
    return '<br/>'.join(result)

audit = REPO/'.living/tasks/fer-fig6-independent-audit-2026-09-20/outputs/claude-science/AUDIT_REPORT.md'
story = [para('Validation status and source guide', 'title'),
    para('<b>Prototype; not yet carefully vetted.</b> The figures, derived tables and reference code can support discussion and follow-up. They are not final publication results.'),
    para('The latest presentation revision aligns panels a and c at 110 mm, makes DOF bold, and sets each matrix cell to about 7.1 mm with 8.5 pt numbers and about 1.46 mm side padding. Its source and plotted tables are byte-identical to the preceding version. Rendering and arithmetic checks do not establish scientific validity.'),
    para('Specific source issues and repair gates', 'head'),
    para('<b>Endodermis:</b> three fer-4 contexts have filtered graphs that disagree with their stored outgoing-degree scores. First compare the filtered graphs with the existing raw links and declared filter. Recompute scores from the confirmed graph in isolated working directories if needed. A fresh inference fit is not yet shown to be necessary. If scores change, rerun nominations, enrichment and GRN-dependent overlaps.'),
    para('<b>Phloem and xylem:</b> substring-based target-feature selection admitted genes selected from pole-pericycle DE contexts into the separate tissue fits. This concerns gene features, not accidental inclusion of pericycle cells. Correcting these models requires an explicit map and controlled refitting; deleting edges afterward is not equivalent.'),
    para('<b>Interpretation:</b> nominations reflect relative network ranks, not measured TF activity. Context unions do not demonstrate co-occurrence within one cell or intercellular transport. YUC5 was excluded before fitting in several tissues, including endodermis. Its missing edges there are not negative biological evidence.'),
    para('Audit status', 'head'),
    para(('The requested independent audit report is available: '+linked(audit)+' and '+linked(HERE/'AUDIT_RECEIPT_VERIFICATION.md', 'verification and wording correction')+'. Its scope and unresolved findings must be read before reuse. Audit completion is not publication approval.') if audit.exists() else 'The requested independent Claude Science audit is still in progress. No audit completion or source repair is claimed in this version.'),
    para('Detailed evidence: '+linked(HERE/'VALIDATION_STATUS.md', 'validation status')+'; '+linked(HERE.parent/'review_v2/SOURCE_AVAILABILITY.md', 'source availability')+'; '+linked(HERE.parent/'review_v2/audit_response/ENDODERMIS_TRIAGE.md', 'endodermal triage')+'.'),
    para('Available data and code', 'head'),
    para('The pinned <a href="https://github.com/tmnolan/feronia-cell-type-specific/tree/5649a9485de0d50d03c8135c5ff254bdaa8d0a1a" color="#0072B2">paper GitHub repository</a> contains analysis code and exported tables. Integrated scRNA-seq, CellOracle input/raw-link/filtered-link objects and full DE tables are present on Resnick. Exact remote paths, check dates and limits are in the linked source record; large objects were not all reloaded for this layout revision.'),
    para('The '+linked(HERE/'provenance/CURRENT_INPUT_OUTPUT_CODE.csv', 'complete input/output/code inventory')+' records exact repository-relative paths, calculation and plotting scripts, outputs, commands and dependencies for all six panel types. '+linked(HERE/'README.md', 'Current presentation notes')+' identify the newer forest and matrix.'),
    para('The '+linked(V3/'code_reference/README.md', 'code reference README')+' and '+linked(V3/'code_reference/source_copy_manifest.csv', 'source/copy hash manifest')+' accompany 30 reference scripts, each labelled prototype. They are source snapshots with original project dependencies, not a self-contained portable pipeline. Current plotting code is linked on the following pages. No GitLab release has been made.'),
    para('Additive reproducibility check: '+linked(V5/'provenance/reconstruct_memberships.R', 'membership reconstruction R script')+' and its '+linked(V5/'provenance/README.md', 'comparison report')+'. This reproduces memberships from the recorded nomination and DE tables; it does not repair upstream network scores.'),
    PageBreak()]

rows = list(csv.DictReader((HERE/'provenance/CURRENT_INPUT_OUTPUT_CODE.csv').open()))
for block, subset in enumerate((rows[:3], rows[3:])):
    story.append(para('Panel inputs, outputs and code' + (' (continued)' if block else ''), 'title'))
    story.append(para('Links open local files in this checkout. Original source paths and hashes are retained in the CSV inventory. Code snapshots are explicitly marked prototype.', 'small'))
    for row in subset:
        title = row['panel']
        plotcode = row['plotting_code']
        outputs = row['outputs']
        size = row['final_size_mm']
        outputs = ' ; '.join(p for p in paths(outputs) if p.endswith('.pdf'))
        # SVG helpers are recorded in the full inventory; the concise PDF links the main R entry points.
        plotcode = ' ; '.join(p for p in paths(plotcode) if p.endswith('.R'))
        story.extend([
            para(escape(title)+' | '+escape(size)+' mm', 'head'),
            para('<b>Inputs:</b> '+listlinks(row['input_tables_read_by_plot_script']), 'small'),
            para('<b>Calculation:</b> '+escape(row['calculation']), 'small'),
            para('<b>Plotting code:</b> '+listlinks(plotcode, code=True)+'<br/><b>Output:</b> '+listlinks(outputs)+' (PDF; PNG and editable SVG are saved alongside)', 'small'),
            para('<b>Scope:</b> '+escape(row['prototype_caveats']), 'small'),
            Spacer(1, 3*mm),
        ])
    if block == 0:
        story.append(PageBreak())

def footer(c, doc):
    c.saveState()
    c.setFont('Arial',7)
    c.setFillColor(HexColor('#555555'))
    c.drawString(14*mm, 10*mm, 'Prototype; not yet carefully vetted. Local source links require the project checkout.')
    c.drawRightString(196*mm,10*mm,f'Appendix {doc.page}')
    c.restoreState()

appendix = HERE/'sources_and_validation.pdf'
SimpleDocTemplate(str(appendix), pagesize=(210*mm,297*mm),
    leftMargin=14*mm,rightMargin=14*mm,topMargin=15*mm,bottomMargin=19*mm,
    title='FER Figure 6 prototype: data, code and validation status',author='Nolan Lab').build(story,onFirstPage=footer,onLaterPages=footer)
writer=PdfWriter()
writer.append(str(HERE/'review_at_print_size_v6.pdf'))
writer.append(str(appendix))
writer.add_metadata({'/Title':'FER Figure 6 review with source guide; prototype','/Author':'Nolan Lab'})
output=HERE/'review_with_sources_v6.pdf'
with output.open('wb') as stream:
    writer.write(stream)
print(output)
print(f'Appendix pages: {len(PdfReader(appendix).pages)}; total pages: {len(writer.pages)}')
