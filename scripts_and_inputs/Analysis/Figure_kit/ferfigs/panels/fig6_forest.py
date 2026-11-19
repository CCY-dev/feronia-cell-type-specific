"""Fig. 6a: TF-family enrichment forest plot, drawn from Chia-Yu Chen's export TF_family_enrichment_all_tests.csv
(columns genotype, family, nominated, eligible, OR, lo, hi, p, q), used unchanged.

Display rule (as coded in compact_panels.R and draw_enrichment_a.py; parameters in config fig6.yaml a.rule):
families with at least `min_eligible` eligible TFs; the `n_families` with the largest odds ratio in either genotype;
rows ordered by the odds ratio of `order_by`; WT above fer-4 in each row; filled point = q < `q_cut`; log axis.
An OR of 0 has no position on a log axis: its interval is drawn from the upper 95% limit towards zero, as an arrow.
Nothing is recomputed: OR, CI, q and counts are printed or plotted exactly as in the table.
"""
from __future__ import annotations

import math

import pandas as pd

from .. import style as S
from ..schemas import ValidationError
from .fig6_common import MM, PT, arrow_head, cap_mm

REQUIRED = ['genotype', 'family', 'nominated', 'eligible', 'OR', 'lo', 'hi', 'q']


def read(path):
    df = pd.read_csv(path)
    miss = [c for c in REQUIRED if c not in df.columns]
    if miss:
        raise ValidationError([f'{path}: missing column(s) {miss}; expected the TF_family_enrichment export columns '
                               f'{REQUIRED}.'])
    return df


def select(df, rule):
    """Same rule and tie handling as draw_enrichment_a.py (stable sort in table order)."""
    elig = df[df['eligible'] >= rule['min_eligible']]
    best = {}
    for fam, orr in zip(elig['family'], elig['OR']):
        best[fam] = max(best.get(fam, -1.0), float(orr))
    top = [f for f, _ in sorted(best.items(), key=lambda kv: -kv[1])[:rule['n_families']]]
    sel = elig[elig['family'].isin(top)].copy()
    ref = sel[sel['genotype'] == rule['order_by']]
    order = [f for f in ref.sort_values('OR', ascending=False, kind='stable')['family']]
    if len(order) != len(top):
        raise ValidationError([f"enrichment table: not every selected family has a '{rule['order_by']}' row."])
    return sel, order


def draw(page, mm: MM, box, df, cfg, gid='fig6-a'):
    st = page.st
    fs = cfg.get('text_pt', st['font']['size_pt'])            # names, headers, axis title
    fs_s = cfg.get('small_pt', st['font']['small_pt'])        # counts, tick labels, key
    x, y, w, h = box
    rule = cfg['rule']
    sel, order = select(df, rule)
    genos = cfg['genotypes']                          # [{name, label, colour}, ...] top to bottom within a row
    lo_or, hi_or = cfg['or_range']
    bad = sel[((sel['OR'] > 0) & (sel['lo'] < lo_or)) | (sel['hi'] > hi_or)]
    if len(bad):
        raise ValidationError([f"enrichment rows outside the axis range {cfg['or_range']}: "
                               f"{bad[['genotype', 'family']].values.tolist()}; widen a.or_range in config/fig6.yaml."])
    X0, X1 = x + cfg['x_axis'][0], x + cfg['x_axis'][1]
    lmin, lmax = math.log10(lo_or), math.log10(hi_or)

    def xp(v):
        return X0 + (math.log10(v) - lmin) / (lmax - lmin) * (X1 - X0)
    pitch, off = cfg['pitch'], cfg['offset']
    y0 = y + cfg['y_row0']
    yc = {f: y0 + i * pitch for i, f in enumerate(order)}
    x_name = x + cfg['x_name']
    labels = cfg.get('family_labels', {})
    bold = set(cfg.get('bold_families', []))
    cols = {g['name']: x + c for g, c in zip(genos, cfg['x_counts'])}
    colour = {g['name']: g['colour'] for g in genos}
    sign = {genos[0]['name']: -1, genos[1]['name']: 1}
    coloured_text = bool(cfg.get('coloured_text', False))
    cap, cap_s = cap_mm(fs), cap_mm(fs_s)
    lw_ci, lw_pt = cfg['ci_lw_pt'], cfg['point_edge_pt']
    d_pt = cfg['point_d']
    # row bands and family names
    for i, f in enumerate(order):
        if i % 2 == 0:
            mm.rect(x_name + cfg['band_x'][0], yc[f] - pitch / 2, x + cfg['band_x'][1] - x_name - cfg['band_x'][0],
                    pitch, face=cfg['band_colour'], gid=f'{gid}-band', z=0)
        page.text(x_name, yc[f] + cap / 2, labels.get(f, f.replace('_', ' ')), size=fs, ha='right',
                  weight='bold' if f in bold else None, gid=f'{gid}-family')
    # OR = 1 reference
    mm.line(xp(1), y0 - pitch / 2 - 0.3, xp(1), y0 + (len(order) - 1) * pitch + pitch / 2 + 0.3,
            lw=cfg['ref_lw_pt'], colour=cfg['ref_colour'], dashes=cfg['ref_dashes'], gid=f'{gid}-or1', z=1)
    rows = []
    for _, r in sel.iterrows():
        g = r['genotype']
        yy = yc[r['family']] + sign[g] * off
        col = colour[g]
        if float(r['OR']) == 0:
            xa, xb = xp(float(r['hi'])), xp(cfg['zero_arrow_to'])
            mm.line(xa, yy, xb + cfg.get('arrow_shaft_gap', 0.8), yy, lw=lw_ci, colour=col, gid=f'{gid}-ci', z=3)
            arrow_head(mm, xb, yy, math.pi, 1.0, 0.9, col, gid=f'{gid}-arrow', z=3)
        else:
            mm.line(xp(float(r['lo'])), yy, xp(float(r['hi'])), yy, lw=lw_ci, colour=col, gid=f'{gid}-ci', z=3)
            filled = float(r['q']) < rule['q_cut']
            mm.circle(xp(float(r['OR'])), yy, d_pt, face=col if filled else '#ffffff', edge=col, lw=lw_pt,
                      gid=f'{gid}-point', z=4)
        page.text(cols[g], yc[r['family']] + cap_s / 2, f"{int(r['nominated'])}/{int(r['eligible'])}", size=fs_s,
                  ha='center', colour=col if coloured_text else None, gid=f'{gid}-counts')
        rows.append(dict(row=order.index(r['family']) + 1, family=r['family'], genotype=g,
                         nominated=int(r['nominated']), eligible=int(r['eligible']), OR=float(r['OR']),
                         lo=float(r['lo']), hi=float(r['hi']), q=float(r['q']),
                         filled=bool(float(r['q']) < rule['q_cut']), zero_or_arrow=float(r['OR']) == 0))
    # axis
    ya = y + cfg['y_axis']
    lw_ax = cfg.get('axis_pt', st['lines']['axis_pt'])
    ax_col = cfg.get('axis_colour', 'k')
    mm.line(X0, ya, X1, ya, lw=lw_ax, colour=ax_col, gid=f'{gid}-axis')
    tick = cfg.get('tick_len', st['lines']['tick_len_pt'] * PT)
    for v in cfg['ticks']:
        mm.line(xp(v), ya, xp(v), ya + tick, lw=cfg.get('tick_pt', st['lines']['tick_pt']), colour=ax_col,
                gid=f'{gid}-axis')
        if 'y_tick' in cfg:                                   # tick-label baseline (mm below the panel top)
            page.text(xp(v), y + cfg['y_tick'], f'{v:g}', size=fs_s, ha='center', gid=f'{gid}-ticklabels')
        else:
            page.text(xp(v), ya + tick + 0.6, f'{v:g}', size=fs_s, ha='center', va='top', gid=f'{gid}-ticklabels')
    y_title = y + cfg['y_title']
    page.text((X0 + X1) / 2, y_title, cfg['axis_title'], size=fs, ha='center', gid=f'{gid}-axistitle')
    # header row and column labels (a marker in the series colour keys the colours; text stays black)
    yh, ys = y + cfg['y_head'], y + cfg['y_sub']
    page.text(x_name, yh, cfg['header_names'], size=fs, ha='right', weight='bold', gid=f'{gid}-header')
    cx = sum(cols.values()) / 2
    page.text(cx, yh, cfg['header_counts'], size=fs, ha='center', weight='bold', gid=f'{gid}-header')
    for g in genos:
        tw = page.text(cols[g['name']], ys, g['label'], size=fs, ha='center',
                       colour=g['colour'] if coloured_text else None, gid=f'{gid}-header')
        if not coloured_text:
            mm.circle(cols[g['name']] - tw / 2 - 1.2, ys - cap / 2, d_pt, face=g['colour'], edge=g['colour'],
                      lw=lw_pt, gid=f'{gid}-key', z=4)
    # significance key, right of the axis title on its baseline (key_x: one x per entry, or the first x and a gap)
    ky = y_title - cap_s / 2
    tgap = cfg.get('key_text_gap', 1.2)
    kxs = cfg['key_x'] if isinstance(cfg['key_x'], list) else None
    kx = x + (kxs[0] if kxs else cfg['key_x'])
    for i, (lab, filled) in enumerate(cfg['key']):
        if kxs:
            kx = x + kxs[i]
        mm.circle(kx, ky, d_pt, face=cfg['key_colour'] if filled else '#ffffff', edge=cfg['key_colour'], lw=lw_pt,
                  gid=f'{gid}-key', z=4)
        kx += tgap + page.text(kx + tgap, y_title, lab, size=fs_s, gid=f'{gid}-key') + cfg.get('key_gap', 2.4)
    return dict(type='forest', families=order, rows=sorted(rows, key=lambda d: (d['row'], d['genotype'] != 'WT')),
                rule=rule, n_rows_in_table=int(len(df)))
