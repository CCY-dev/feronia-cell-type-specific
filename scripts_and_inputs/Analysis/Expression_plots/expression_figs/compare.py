"""Before/after review pages: the current deck page (vector, unchanged) next to the two draft options.

One A3-landscape page per figure: BEFORE | option A | option B, each scaled to fit its column (scale printed).
"""
from __future__ import annotations

import hashlib
import pathlib

MM = 72 / 25.4


def build_compare(spec, out):
    import fitz
    deck_path = pathlib.Path(spec['deck'])
    if not deck_path.exists():
        print('compare: deck not found, skipped:', deck_path)
        return []
    deck = fitz.open(deck_path)
    sha = hashlib.sha256(deck_path.read_bytes()).hexdigest()
    doc = fitz.open()
    W, H = 420 * MM, 297 * MM
    margin, top = 8 * MM, 16 * MM
    col_w = (W - 2 * margin - 2 * 6 * MM) / 3
    titles = {'dot': 'AFTER, option A: dot plot (percent expressing + mean)',
              'violin': 'AFTER, option B: split violin plot'}
    for pg in spec['pages']:
        page = doc.new_page(width=W, height=H)
        page.insert_text((margin, 9 * MM), 'FER single-cell expression panels: before / after (draft for review)',
                         fontsize=11, fontname='helv')
        cols = [('before', None)] + [('after', a) for a in pg['after']]
        for k, (kind, fid) in enumerate(cols):
            x0 = margin + k * (col_w + 6 * MM)
            if kind == 'before':
                src, pno = deck, pg['before_page'] - 1
                label = f"BEFORE: {pg['label']}, deck page {pg['before_page']}"
            else:
                src, pno = fitz.open(out / 'figures' / f'{fid}.pdf'), 0
                label = titles['dot' if 'dot' in fid else 'violin'] + f'  [{fid}]'
            r = src[pno].rect
            s = min(col_w / r.width, (H - top - 10 * MM) / r.height)
            dest = fitz.Rect(x0, top, x0 + r.width * s, top + r.height * s)
            page.draw_rect(fitz.Rect(dest.x0 - 1, dest.y0 - 1, dest.x1 + 1, dest.y1 + 1), color=(0.75, 0.75, 0.75),
                           width=0.4)
            page.show_pdf_page(dest, src, pno)
            page.insert_textbox(fitz.Rect(x0, top - 6 * MM, x0 + col_w, top - 0.5 * MM),
                                f'{label}  (shown at {s * 100:.0f}%)', fontsize=7.5, fontname='helv')
        page.insert_text((margin, H - 4 * MM), f'Deck: {deck_path.name} (sha256 {sha[:16]}...). Drafts are vector '
                         'PDFs; UMAP points are rasterised at 600 ppi inside them.', fontsize=6.5, fontname='helv')
    cdir = out / 'compare'
    cdir.mkdir(exist_ok=True)
    p = cdir / 'before_after.pdf'
    doc.save(p, garbage=3, deflate=True, no_new_id=True)
    paths = [p]
    for i, page in enumerate(doc):
        q = cdir / f'before_after_p{i + 1}.png'
        page.get_pixmap(dpi=110).save(q)
        paths.append(q)
    return [str(x) for x in paths]
