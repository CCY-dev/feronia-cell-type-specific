"""Check built figures against the journal profile (ferfigs/profiles/<profile>.yaml).

Reads the PDF actually produced (PyMuPDF) plus the build record (<stem>.build.json) and reports, per figure:
page size, every text size and font, text colours, stroke widths, embedded font types, effective resolution of
every image at printed size, panel-letter sequence and style, and gaps between juxtaposed images.
Errors break a journal rule; warnings are limits that a figure should normally meet (e.g. the 170-mm height that
leaves room for the legend) or items an author must confirm.
"""
from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path

from . import style as S

PT_MM = 25.4 / 72


def _rgb(c):
    return ((c >> 16 & 255) / 255, (c >> 8 & 255) / 255, (c & 255) / 255)


def check_pdf(pdf_path, build_record=None, profile=None):
    import fitz
    prof = profile or S.load_profile()
    pdf_path = Path(pdf_path)
    rec = build_record
    if rec is None:
        bj = pdf_path.with_suffix('.build.json')
        rec = json.loads(bj.read_text()) if bj.exists() else {}
    doc = fitz.open(pdf_path)
    pg = doc[0]
    err, warn, info, codes = [], [], {}, {}

    def E(code, msg):
        err.append(msg)
        codes[msg] = code

    def W(code, msg):
        warn.append(msg)
        codes[msg] = code
    w_mm, h_mm = pg.rect.width * PT_MM, pg.rect.height * PT_MM
    info['size_mm'] = [round(w_mm, 1), round(h_mm, 1)]
    if w_mm > prof['width_mm']['max'] + 0.05:
        E('width', f'width {w_mm:.1f} mm exceeds {prof["width_mm"]["max"]} mm')
    if h_mm > prof['height_mm']['hard_max'] + 0.05:
        E('page_depth', f'height {h_mm:.1f} mm exceeds the page depth {prof["height_mm"]["hard_max"]} mm')
    elif h_mm > prof['height_mm']['target_max'] + 0.05:
        W('height_over_target', f'height {h_mm:.1f} mm is above {prof["height_mm"]["target_max"]} mm (no room for the legend '
                    f'under the figure; production may shrink it)')
    # ---- text
    letters = set((rec.get('letters') or {}).keys())
    sizes, fonts, bad_col, letter_spans = Counter(), Counter(), Counter(), []
    tmin, tmax = prof['text_pt']['min'], prof['text_pt']['max']
    small, big = [], []
    for b in pg.get_text('dict')['blocks']:
        for l in b.get('lines', []):
            for s in l['spans']:
                t = s['text'].strip()
                if not t:
                    continue
                sz = round(s['size'], 2)
                is_letter = (t in letters and 'Bold' in s['font'] and abs(sz - prof['panel_label']['size_pt']) < 0.05)
                if is_letter:
                    letter_spans.append(t)
                    continue
                sizes[sz] += len(t)
                fonts[s['font']] += 1
                if sz < tmin - 0.01:
                    small.append((t[:20], sz))
                if sz > tmax + 0.01:
                    big.append((t[:20], sz))
                c = _rgb(s['color'])
                tol = prof['text_colour']['tolerance']
                if not any(all(abs(a - b_) <= tol for a, b_ in zip(c, ok)) for ok in prof['text_colour']['allowed']):
                    bad_col[tuple(round(v, 2) for v in c)] += 1
    info['text_sizes_pt'] = dict(sorted(sizes.items()))
    info['fonts'] = dict(fonts)
    if small:
        E('text_small', f'{len(small)} text run(s) below {tmin} pt, e.g. {small[:3]}')
    if big:
        E('text_big', f'{len(big)} text run(s) above {tmax} pt, e.g. {big[:3]}')
    allowed = prof['fonts']['allowed'] + prof['fonts'].get('allowed_monospace', [])
    badf = [f for f in fonts if not any(a.replace(' ', '').lower() in f.replace(' ', '').lower() for a in allowed)]
    if badf:
        E('font', f'font(s) outside {prof["fonts"]["allowed"]}: {badf}')
    if bad_col:
        W('coloured_text', f'coloured text (journal: avoid coloured text): {dict(bad_col)}')
    # embedded font types
    ftypes = Counter(f[2] for f in doc.get_page_fonts(0))
    info['embedded_font_types'] = dict(ftypes)
    if any(t in ('Type3',) for t in ftypes):
        E('type3', 'Type 3 fonts present (text not editable)')
    # ---- panel letters
    want = sorted(letters)
    if want:
        seq = [chr(ord('a') + i) for i in range(len(want))]
        if want != seq:
            E('letters_sequence', f'panel letters {want} are not a continuous sequence a, b, c ...')
        missing = [l for l in want if l not in letter_spans]
        if missing:
            E('letters_missing', f'panel letters {missing} not found as {prof["panel_label"]["size_pt"]} pt bold text')
    # ---- line widths
    widths = Counter()
    for d in pg.get_drawings():
        if d.get('type') in ('s', 'fs') and d.get('width'):
            widths[round(d['width'], 2)] += 1
    info['stroke_widths_pt'] = dict(sorted(widths.items()))
    lo, hi = prof['line_width_pt']['min'], prof['line_width_pt']['max']
    out_w = {k: v for k, v in widths.items() if k < lo - 0.005 or k > hi + 0.005}
    if out_w:
        W('stroke_width', f'stroke widths outside {lo}-{hi} pt: {out_w}')
    # ---- images
    imgs = []
    for im in pg.get_image_info():
        r = fitz.Rect(im['bbox'])
        if r.width <= 0:
            continue
        ppi_h = im['width'] / (r.width / 72)                  # pixels per inch across ...
        ppi_v = im['height'] / (r.height / 72) if r.height > 0 else ppi_h   # ... and down; the lower one counts
        imgs.append(round(min(ppi_h, ppi_v)))
    info['image_ppi'] = imgs
    lowres = [p for p in imgs if p < prof['images']['min_ppi'] - 0.5]
    info['image_ppi_low'] = lowres
    if lowres:
        E('image_ppi', f'{len(lowres)} image(s) below {prof["images"]["min_ppi"]} ppi at printed size: {lowres}')
    # gaps between juxtaposed images (build record)
    gap_min = prof['integrity']['juxtaposed_images_gap_mm']
    ims = rec.get('images', [])
    close = []
    for i, a in enumerate(ims):
        for b in ims[i + 1:]:
            ax0, ay0, ax1, ay1 = a['x_mm'], a['y_mm'], a['x_mm'] + a['w_mm'], a['y_mm'] + a['h_mm']
            bx0, by0, bx1, by1 = b['x_mm'], b['y_mm'], b['x_mm'] + b['w_mm'], b['y_mm'] + b['h_mm']
            yo = min(ay1, by1) - max(ay0, by0); xo = min(ax1, bx1) - max(ax0, bx0)
            if yo > 0.5:
                g = max(bx0 - ax1, ax0 - bx1)
                if -0.01 < g < gap_min - 0.01 or (g <= -0.01 and xo > 0):
                    close.append((a['name'], b['name'], round(g, 2)))
            elif xo > 0.5:
                g = max(by0 - ay1, ay0 - by1)
                if -0.01 < g < gap_min - 0.01:
                    close.append((a['name'], b['name'], round(g, 2)))
    if close:
        W('image_gap', f'{len(close)} image pair(s) closer than {gap_min} mm (demarcate juxtaposed images): {close[:4]}')
    return dict(figure=pdf_path.stem, errors=err, warnings=warn, info=info, codes=codes)


# Which check messages a figure config may accept on purpose ('accepted:' rules). Journal errors about page width,
# page depth, text size, fonts, Type 3 text and panel letters can never be accepted. Height above the legend guide and
# image resolution can be accepted only within stated bounds, so a later change cannot hide behind an old decision.
ACCEPTABLE = {'height_over_target': ('max_mm',), 'image_ppi': ('min_ppi', 'max_count'),
              'coloured_text': (), 'stroke_width': (), 'image_gap': ()}


def _within(code, rule, info):
    """(True, '') when the measured value is inside the bounds the rule accepts, else (False, why)."""
    if code == 'height_over_target':
        h = info.get('size_mm', [0, 0])[1]
        if h > float(rule['max_mm']) + 0.05:
            return False, f"the figure is now {h} mm high; the rule accepts up to {rule['max_mm']} mm"
    if code == 'image_ppi':
        low = info.get('image_ppi_low', [])
        if len(low) > int(rule['max_count']) or (low and min(low) < float(rule['min_ppi'])):
            return False, (f"{len(low)} image(s) down to {min(low)} ppi; the rule accepts at most {rule['max_count']} "
                           f"image(s) at {rule['min_ppi']} ppi or more")
    return True, ''


def apply_accepted(r, rules):
    """Move check messages that the figure config accepts on purpose (config 'accepted:', each rule with a 'match'
    substring, a 'reason' and, for height and resolution, the bounds it accepts) out of errors/warnings into
    r['accepted'], so they stay visible with their reason. These are author decisions, not journal waivers.
    A rule without a reason, a rule that matches a message that cannot be accepted, a rule without the bounds its
    message needs, or a measured value outside the bounds leaves the message in place and adds an error or warning
    saying why. Rules that match nothing are listed in r['accepted_unused']."""
    rules = list(rules or [])
    codes = r.get('codes') or {}
    good, notes_e, notes_w = [], [], []
    for rule in rules:
        if not str(rule.get('reason', '')).strip():
            notes_e.append(f"accepted rule {rule.get('match')!r} has no reason; give one in the config")
        elif not str(rule.get('match', '')).strip():
            notes_e.append(f"accepted rule with reason {str(rule['reason'])[:40]!r}... has no 'match' text")
        else:
            good.append(rule)
    acc, used = [], set()
    for key in ('errors', 'warnings'):
        keep = []
        for m in r.get(key, []):
            i = next((k for k, x in enumerate(good) if str(x['match']) in m), None)
            if i is None:
                keep.append(m)
                continue
            rule, code = good[i], codes.get(m)
            used.add(i)
            if code not in ACCEPTABLE:
                keep.append(m)
                notes_e.append(f"accepted rule {rule['match']!r} matches a journal error that cannot be accepted "
                               f"('{m}'); fix the figure or narrow the rule")
                continue
            missing = [b for b in ACCEPTABLE[code] if rule.get(b) is None]
            if missing:
                keep.append(m)
                notes_e.append(f"accepted rule {rule['match']!r} must state the bounds it accepts: "
                               f"{', '.join(missing)}")
                continue
            ok, why = _within(code, rule, r.get('info', {}))
            if not ok:
                keep.append(m)
                (notes_e if key == 'errors' else notes_w).append(
                    f"accepted rule {rule['match']!r} no longer covers this: {why}")
                continue
            acc.append(dict(level=key[:-1], message=m, reason=rule['reason'], decided_by=rule.get('decided_by', ''),
                            bounds={b: rule[b] for b in ACCEPTABLE[code]}))
        r[key] = keep
    r['errors'] = r.get('errors', []) + notes_e
    r['warnings'] = r.get('warnings', []) + notes_w
    r['accepted'] = acc
    r['accepted_unused'] = [str(x['match']) for k, x in enumerate(good) if k not in used]
    return r


def report_markdown(results):
    lines = ['| Figure | Size (mm) | Errors | Warnings | Data warnings | Accepted on purpose: author decision, not a '
             'journal waiver (reason) |', '|---|---|---|---|---|---|']
    for r in results:
        acc = '<br>'.join(f"{a['message']} ({a['reason']})" for a in r.get('accepted', [])) or 'none'
        dw = '<br>'.join(r.get('data_warnings', [])) or 'none'
        lines.append(f"| {r['figure']} | {r['info']['size_mm'][0]} x {r['info']['size_mm'][1]} | "
                     f"{'<br>'.join(r['errors']) or 'none'} | {'<br>'.join(r['warnings']) or 'none'} | {dw} | {acc} |")
    return '\n'.join(lines)
