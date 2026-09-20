#!/usr/bin/env python3
# Prototype; not yet carefully vetted.
"""Prototype packaging helper. Add a cover and replace local proof links with GitLab links.

This does not recalculate any scientific result or modify the retained source PDF.
Run from any directory; paths resolve relative to this package by default.
Requires Python 3, reportlab and pypdf.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import quote, unquote, urlsplit

from pypdf import PdfReader, PdfWriter
from pypdf.generic import NameObject, TextStringObject
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph

NAME = "FER-DOF Network Figure Explorer"
DEFAULT_BASE = (
    "https://gitlab.caltech.edu/nolan_lab/feronia-cell-type-specific/-/blob/"
    "codex/fer-dof-figure-explorer/exploratory/fer-dof-network-figure-explorer/"
)
PROOF = Path("revisions/2026-09-20_compact_figures/review_v7/review_with_sources_v7.pdf")
PREVIEW = PROOF.parent / "review-1.png"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def package_path(uri, root):
    # Legacy file URIs store the source project name. Do not depend on its host path.
    raw = unquote(urlsplit(uri).path)
    marker = "/fer-sc-manuscript-2026/"
    if marker not in raw:
        raise ValueError(f"Unrecognized local link: {uri}")
    rel = Path(raw.split(marker, 1)[1])
    candidates = [rel]
    copied = "review_v3/code_reference/"
    if copied in str(rel):
        candidates.append(Path(str(rel).split(copied, 1)[1]))
    for candidate in candidates:
        resolved = (root / candidate).resolve()
        if resolved.is_relative_to(root) and resolved.is_file():
            return candidate.as_posix()
    raise FileNotFoundError(f"No packaged target for {uri}")


def make_cover(root, base, dimensions):
    width, height = dimensions
    buf = io.BytesIO()
    c = canvas.Canvas(buf, pagesize=dimensions, invariant=1)
    ink = colors.HexColor("#20313B")
    muted = colors.HexColor("#53616B")
    accent = colors.HexColor("#12669A")
    c.setFillColor(ink)
    c.setFont("Helvetica-Bold", 25)
    c.drawString(36, height - 51, NAME)
    c.setFillColor(muted)
    c.setFont("Helvetica", 12)
    c.drawString(36, height - 75, "Figure 6 exploration and source guide | Nolan Lab | 20 September 2026")
    c.setFillColor(colors.HexColor("#EEF2F4"))
    c.roundRect(36, height - 125, width - 72, 31, 4, stroke=0, fill=1)
    c.setFillColor(ink)
    c.setFont("Helvetica-Bold", 13)
    c.drawString(48, height - 113, "Prototype; not yet carefully vetted")

    style = ParagraphStyle("body", fontName="Helvetica", fontSize=11, leading=16, textColor=ink)
    def para(text, x, top, w):
        p = Paragraph(text, style)
        _, h = p.wrap(w, height)
        p.drawOn(c, x, top-h)
        return top-h

    x, w, top = 36, 362, height-147
    top = para("<b>Start here for Chia-Yu</b>", x, top, w) - 9
    top = para("The following six pages retain the current v7 figure proof and its source appendix, including the DOF-above, auxin-below panel D. Original proof page numbers are preserved.", x, top, w) - 12
    top = para("The package contains plotting code, seven immediate saved-table inputs, editable SVGs, panel PDFs and PNGs, an input/output/code map, and selected audit records.", x, top, w) - 12
    top = para("<b>What is checked:</b> saved-table arithmetic and memberships, recorded source-file identities, panel-D edge identity and support classes, and the rendered proof. This packaging step changes navigation only.", x, top, w) - 12
    top = para("<b>What remains unresolved:</b> three fer-4 endodermal graph/score mismatches and phloem/xylem gene-feature selection. Source repair and downstream recalculation remain pending. These figures do not establish TF activity or intercellular regulation.", x, top, w)

    preview = root / PREVIEW
    if preview.is_file():
        c.drawImage(ImageReader(str(preview)), 427, height-374, width=377, height=238,
                    preserveAspectRatio=True, anchor="c", mask="auto")
    c.setFillColor(muted)
    c.setFont("Helvetica", 9)
    c.drawString(437, height-387, "Current four-panel proof; full vector version follows.")

    links = [("Start-here README and reproduction notes", "README.md"),
             ("Panel inputs, outputs and code", "revisions/2026-09-20_compact_figures/review_v7/provenance/CURRENT_INPUT_OUTPUT_CODE.csv"),
             ("Validation status and next analytical checks", "VALIDATION_STATUS.md")]
    y = 84
    for label, rel in links:
        if not (root/rel).is_file():
            raise FileNotFoundError(rel)
        c.setFont("Helvetica-Bold", 10)
        c.setFillColor(accent)
        c.drawString(36, y, label)
        c.linkURL(base + quote(rel, safe="/"), (36, y-3, 406, y+12), relative=0)
        y -= 18
    c.setFillColor(muted)
    c.setFont("Helvetica", 8)
    c.drawRightString(width-36, 32, "Private lab GitLab links require repository access.")
    c.showPage()
    c.save()
    buf.seek(0)
    return PdfReader(buf)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--gitlab-base", default=DEFAULT_BASE)
    args = parser.parse_args()
    root = args.root.resolve()
    base = args.gitlab_base.rstrip("/") + "/"
    if urlsplit(base).scheme != "https":
        raise ValueError("GitLab base must use HTTPS")
    source = root/PROOF
    before = digest(source)
    reader = PdfReader(source)
    dimensions = tuple(float(v) for v in (reader.pages[0].mediabox.width, reader.pages[0].mediabox.height))
    writer = PdfWriter()
    writer.append(make_cover(root, base, dimensions))
    writer.append(reader)
    mapping = []
    for page_num, page in enumerate(writer.pages, 1):
        for ref in page.get("/Annots", []):
            annotation = ref.get_object()
            action = annotation.get("/A", {})
            uri = str(action.get("/URI", ""))
            if uri.startswith("file:"):
                rel = package_path(uri, root)
                url = base + quote(rel, safe="/")
                action[NameObject("/URI")] = TextStringObject(url)
                mapping.append({"page": page_num, "package_path": rel, "url": url})
    writer.add_metadata({"/Title": NAME, "/Author": "Nolan Lab", "/Subject": "Prototype; not yet carefully vetted. Figure 6 exploration and source guide."})
    writer.add_outline_item("Start here and validation status", 0)
    writer.add_outline_item("Current figure proof: panels a-d", 1)
    writer.add_outline_item("Network alternatives", 2)
    writer.add_outline_item("Retained UpSet alternative", 3)
    writer.add_outline_item("Sources and audit status", 4)
    output = root / "FER_DOF_Network_Figure_Explorer.pdf"
    with output.open("wb") as handle:
        writer.write(handle)
    final = PdfReader(output)
    for original, packaged in zip(reader.pages, final.pages[1:]):
        assert original.get_contents().get_data() == packaged.get_contents().get_data(), "Proof drawing content changed"
        assert original.extract_text() == packaged.extract_text(), "Proof text changed"
    links = [str(ref.get_object().get("/A", {}).get("/URI", ""))
             for page in final.pages for ref in page.get("/Annots", [])]
    assert not any(uri.startswith("file:") for uri in links)
    assert digest(source) == before, "Source PDF unexpectedly changed"
    report = {"status": "Prototype; not yet carefully vetted", "source_pdf": PROOF.as_posix(),
              "source_pdf_sha256": before, "output_pdf": output.name, "output_pdf_sha256": digest(output),
              "pages": len(writer.pages), "converted_local_links": len(mapping),
              "total_uri_links": len([u for u in links if u]), "local_uri_links_remaining": 0,
              "six_original_page_content_streams_and_text_unchanged": True,
              "linked_package_files_exist": True, "mapping": mapping}
    (root/"provenance").mkdir(exist_ok=True)
    (root/"provenance/SHARE_PDF_LINK_CHECK.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k != "mapping"}, indent=2))


if __name__ == "__main__":
    main()
