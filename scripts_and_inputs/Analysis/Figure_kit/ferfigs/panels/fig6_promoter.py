"""Fig. 6c: pYUC5 promoter with the predicted DOF sites, drawn from the coordinate table
data/fig6/fig6c_promoter_elements.csv and the probe sequences data/fig6/fig6c_probe_sequences.csv
(both written by tools/fig6_promoter/; TAIR10 sequence, Araport11 model). Port of draw_promoter_panel.py.

Coordinates are TSS-relative (+1 = TSS, no position 0); base p occupies [u(p), u(p) + 1) on a continuous axis.
Sense-strand sites are triangles above the line pointing right (5'->3'); antisense-strand sites below, pointing
left. The EMSA probe is a red segment on the line, opened onto the probe and mutant sequences (monospace, the only
non-Arial text in the kit). Highlighted bases (the two TAAAG sites and their substitutions) are marked by a light
red box behind black letters (journal: no coloured text); `highlight: text` restores the candidate's red letters.
"""
from __future__ import annotations

import pandas as pd

from ..schemas import ValidationError
from .fig6_common import MM, cap_mm, mono_advance, mono_family, mono_text

MINUS = '\u2212'


def u(p):
    p = int(p)
    return p if p < 0 else p - 1


def _blocks(s):
    out = set()
    for part in str(s or '').split(';'):
        if part.strip():
            a, b = (int(v) for v in part.split('-'))
            out |= set(range(a - 1, b))
    return out


def draw(page, mm: MM, box, el, seqs, cfg, seq_fonts, gid='fig6-c'):
    st = page.st
    fs = cfg.get('text_pt', st['font']['size_pt'])
    x, y, w, h = box
    C = cfg['colours']
    get = {r['element']: r for _, r in el[el['element'] != 'dof_site'].iterrows()}
    for need in ('reference_interval', 'five_prime_utr', 'atg', 'emsa_probe'):
        if need not in get:
            raise ValidationError([f'fig6c_promoter_elements.csv: no {need!r} row'])
    sites = el[el['element'] == 'dof_site']
    ref, utr, atg, probe = get['reference_interval'], get['five_prime_utr'], get['atg'], get['emsa_probe']
    stub = cfg['cds_stub_bp']
    u_left, u_right = u(ref['tss_rel_start']), u(atg['tss_rel_end']) + 1 + stub
    x0, x1 = x + cfg['track_x'][0], x + w - cfg['track_x'][1]
    S = (x1 - x0) / (u_right - u_left)

    def X(uu):
        return x0 + (uu - u_left) * S
    cap = cap_mm(fs)
    yl, bh, g = y + cfg['y_line'], cfg['box_h'], cfg['gap']
    th, tw = cfg['tri_h'], cfg['tri_w']
    sense_top, anti_top = yl - bh / 2 - g - th, yl + bh / 2 + g
    lw_tri = cfg['tri_edge_pt']

    def tri(xc, top, right, gid_, w_=tw, h_=th, z=3.0):
        pts = ([(xc - w_ / 2, top), (xc + w_ / 2, top + h_ / 2), (xc - w_ / 2, top + h_)] if right else
               [(xc + w_ / 2, top), (xc - w_ / 2, top + h_ / 2), (xc + w_ / 2, top + h_)])
        mm.poly(pts, face='#ffffff', edge='#ffffff', lw=cfg['halo_pt'], join='round', gid=f'{gid_}-halo', z=z)
        mm.poly(pts, face=C['motif'], edge=C['motif'], lw=lw_tri, join='round', gid=gid_, z=z + 0.0005)

    # header: construct name, then the key (right-aligned) saying what the marks are
    yh = y + cfg['y_head']
    head_end = x0 + page.text(x0, yh, cfg['title'], size=fs, weight='bold', gid=f'{gid}-title')
    parts = [('text', cfg['key_text'][0]), ('tri', True), ('text', cfg['key_text'][1]), ('tri', False),
             ('text', cfg['key_text'][2])]
    widths = [page.text_width(v, fs) if k == 'text' else cfg['key_tri_w'] for k, v in parts]
    gaps = cfg['key_gaps']
    kx = x1 - sum(widths) - sum(gaps)
    if kx < head_end + 3:
        raise ValidationError(['c: the key collides with the title; widen panel c or shorten c.key_text.'])
    for n, ((k, v), wd) in enumerate(zip(parts, widths)):
        if k == 'text':
            page.text(kx, yh, v, size=fs, gid=f'{gid}-key')
        else:
            tri(kx + wd / 2, yh - cap / 2 - 0.8, v, f'{gid}-key', w_=cfg['key_tri_w'], h_=1.6, z=5)
        kx += wd + (gaps[n] if n < len(gaps) else 0)
    # coordinate axis: TSS = +1, ticks from config
    ya = y + cfg['y_axis']
    mm.line(X(u_left), ya, X(u(atg['tss_rel_end']) + 1), ya, lw=cfg['axis_pt'], colour=C['axis'], gid=f'{gid}-axis')
    for p in cfg['ticks']:
        xx = X(u(p) + 0.5)
        mm.line(xx, ya - 0.55, xx, ya, lw=cfg.get('tick_pt', cfg['axis_pt']), colour=C['axis'], gid=f'{gid}-axis')
        lab = f'+{p}' if p > 0 else f'{MINUS}{abs(p)}'
        page.text(xx, y + cfg['y_axlab'], lab, size=fs, ha='center', gid=f'{gid}-axis')
    xa = X(u(atg['tss_rel_start']) + 1.5)
    wa = page.text(xa, y + cfg['y_axlab'], cfg['atg_label'], size=fs, ha='center', gid=f'{gid}-atg')
    if xa + wa / 2 > x + w + 0.01:
        raise ValidationError(['c: the ATG label is clipped; reduce c.track_x right margin.'])
    # promoter line, 5' UTR, start of the CDS, TSS bent arrow
    mm.line(X(u_left), yl, X(u(atg['tss_rel_start'])), yl, lw=cfg['dna_pt'], colour=C['dna'], gid=f'{gid}-dna')
    xu0, xu1 = X(u(utr['tss_rel_start'])), X(u(utr['tss_rel_end']) + 1)
    mm.rect(xu0, yl - bh / 2, xu1 - xu0, bh, face=C['utr_fill'], edge=C['utr_edge'], lw=cfg['utr_edge_pt'],
            gid=f'{gid}-utr', z=3)
    wu = page.text((xu0 + xu1) / 2, yl + cap / 2, cfg['utr_label'], size=fs, ha='center', gid=f'{gid}-utr')
    if wu > xu1 - xu0 - 0.4:
        raise ValidationError(["c: the 5' UTR label does not fit in its box."])
    mm.rect(X(u(atg['tss_rel_start'])), yl - bh / 2, X(u_right) - X(u(atg['tss_rel_start'])), bh, face=C['cds'],
            gid=f'{gid}-cds', z=3)
    xt, yar = X(u(1)), y + cfg['y_arrow']
    mm.line(xt, yl - bh / 2, xt, yar, lw=cfg['arrow_pt'], colour='k', gid=f'{gid}-tss', z=4)
    mm.line(xt, yar, xt + 1.3, yar, lw=cfg['arrow_pt'], colour='k', gid=f'{gid}-tss', z=4)
    mm.poly([(xt + 2.0, yar), (xt + 1.2, yar - 0.4), (xt + 1.2, yar + 0.4)], face='k', gid=f'{gid}-tss', z=4)
    # predicted sites, 5' -> 3' so the white halos separate close neighbours
    centres = {'sense': [], 'antisense': []}
    for i, (_, s) in enumerate(sites.assign(_o=sites['tss_rel_start'].astype(int)).sort_values('_o').iterrows()):
        a, b = sorted((u(s['tss_rel_start']), u(s['tss_rel_end'])))
        xc = (X(a) + X(b + 1)) / 2
        right = s['strand'] == 'sense'
        centres[s['strand']].append(xc)
        tri(xc, sense_top if right else anti_top, right, f'{gid}-site', z=3 + i * 0.001)
    for strand, xs in centres.items():
        xs.sort()
        if len(xs) > 1 and min(b - a for a, b in zip(xs, xs[1:])) < 0.7 * tw:
            raise ValidationError([f'c: {strand} site glyphs closer than 0.7 x their width; widen panel c.'])
    # EMSA probe on the line, opened onto the aligned sequences
    pb0, pb1 = yl - cfg['probe_h'] / 2, yl + cfg['probe_h'] / 2
    xp0, xp1 = X(u(probe['tss_rel_start'])), X(u(probe['tss_rel_end']) + 1)
    mm.rect(xp0, pb0, xp1 - xp0, pb1 - pb0, face=C['probe'], gid=f'{gid}-probe', z=4)
    fam = mono_family(seq_fonts)
    fsq = cfg.get('seq_pt', fs)
    cw = mono_advance(page, fam, fsq)
    rows = {r['name']: r for _, r in seqs.iterrows()}
    x_lab = x0
    x_pre = x0 + cfg['seq_label_w']
    x_first = x_pre + 3 * cw
    n = max(len(r['sequence_5to3']) for r in rows.values())
    x_last = x_first + n * cw
    if x_last + 3 * cw > x + w + 0.01:
        raise ValidationError(['c: the probe sequences do not fit; widen panel c or reduce c.seq_pt.'])
    y_rows = [y + v for v in cfg['y_seq']]
    y_zoom = y_rows[0] - cap_mm(fsq) - 0.75
    mm.poly([(xp0, pb1), (xp1, pb1), (x_last + 0.4, y_zoom), (x_first - 0.4, y_zoom)], face=C['zoom'],
            gid=f'{gid}-zoom', z=0)
    hl_text = cfg.get('highlight', 'box') == 'text'
    for (name, lab), yy in zip(cfg['seq_rows'], y_rows):
        if name not in rows:
            raise ValidationError([f'fig6c_probe_sequences.csv: no row named {name!r}'])
        seq = str(rows[name]['sequence_5to3']).upper()
        hl = _blocks(rows[name].get('highlight_1based', ''))
        page.text(x_lab, yy, lab, size=fs, gid=f'{gid}-seqlabel')
        mono_text(page, x_pre, yy, '5\u2032-', fam, fsq, gid=f'{gid}-seq')
        i = 0
        while i < len(seq):
            j, on = i, i in hl
            while j < len(seq) and (j in hl) == on:
                j += 1
            if on and not hl_text:
                mm.rect(x_first + i * cw, yy - cap_mm(fsq) - 0.35, (j - i) * cw, cap_mm(fsq) + 0.7,
                        face=C['highlight_box'], gid=f'{gid}-highlight', z=1)
            col_t = C['probe'] if (on and hl_text) else (0, 0, 0)
            if cfg.get('seq_per_char'):        # proportional font (e.g. Arial): one character centred per cell
                for k in range(i, j):
                    page.text(x_first + (k + 0.5) * cw, yy, seq[k], size=fsq, ha='center', colour=col_t,
                              gid=f'{gid}-seq')
            else:
                mono_text(page, x_first + i * cw, yy, seq[i:j], fam, fsq, gid=f'{gid}-seq', colour=col_t)
            i = j
        mono_text(page, x_last, yy, '-3\u2032', fam, fsq, gid=f'{gid}-seq')
    return dict(type='promoter', mm_per_bp=round(S, 6), displayed_tss_rel=[int(ref['tss_rel_start']),
                int(atg['tss_rel_end']) + stub], n_sites=int(len(sites)),
                n_sense=int((sites['strand'] == 'sense').sum()), n_antisense=int((sites['strand'] == 'antisense').sum()),
                probe_tss_rel=[int(probe['tss_rel_start']), int(probe['tss_rel_end'])], sequence_font=fam,
                seq_char_mm=round(cw, 4))
