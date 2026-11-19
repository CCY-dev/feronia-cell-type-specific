"""Fig. 6e: EMSA gel with its lane-condition grid. The gel is the native image (crop from the manifest, no
pixel change, framed); lane centres (gel pixels), lane groups and the condition of every lane (-, +, G, M, 50x,
100x) come from data/fig6/fig6e_emsa_lanes.csv; row and group labels come from the config (e.rows, e.groups).
A row may map table values to printed text (e.rows[].values, e.g. 50x -> 50 under a 'Competitor (x)' label).
Signs are true minus / plus.
Wide cells that would touch over neighbouring lanes are moved apart as a block; the -/+ signs stay on their lanes.
"""
from __future__ import annotations

from .. import style as S
from ..schemas import ValidationError
from . import images as I
from .fig6_common import MM, cap_mm, spread

MINUS = '\u2212'


def draw(page, mm: MM, fig_dir, row, lanes, cfg, x_gel, y_gel, w_gel, gid='fig6-e'):
    """Vertical positions are baselines in mm relative to the gel image top (negative = above the gel):
    row_baselines (one per condition row, top to bottom), group_baseline (Probe / Mutated Probe) with a plain
    line group_line_gap below it spanning the group's lanes +- group_line_pad; lane numbers lane_number_gap below
    the gel bottom; row_labels place a label beside the gel at a gel-pixel row (uncropped image)."""
    st = page.st
    fs = cfg.get('text_pt', st['font']['size_pt'])
    pm_pt = cfg.get('sign_pt', fs)
    crop = I.parse_crop(row.get('crop_px'))
    pl = I.place(page, fig_dir, row, x_gel, y_gel, w=w_gel, gid=f'{gid}-gel')
    o = cfg.get('frame_outset', 0.0)                   # frame drawn just outside the pixels
    mm.rect(pl.x0 - o, pl.y0 - o, pl.w_mm + 2 * o, pl.h_mm + 2 * o, face='none', edge='k', lw=cfg['frame_pt'],
            gid=f'{gid}-frame', z=3)
    cx0 = crop[0] if crop else 0.0
    cy0 = crop[1] if crop else 0.0
    lx = {int(r['lane']): pl.x0 + (float(r['x_px']) - cx0) * pl.mm_per_px for _, r in lanes.iterrows()}
    for r_ in cfg['rows']:
        if r_['column'] not in lanes.columns:
            raise ValidationError([f"fig6e_emsa_lanes.csv: no column {r_['column']!r} (config e.rows)"])
    cell = {r_['column']: dict(zip(lanes['lane'].astype(int), lanes[r_['column']].astype(str))) for r_ in cfg['rows']}
    x_lab = pl.x0 - cfg['label_gap']
    rows = list(cfg['rows'])
    base = cfg['row_baselines']
    if len(base) != len(rows):
        raise ValidationError([f'e: {len(rows)} condition rows but {len(base)} row_baselines in config/fig6.yaml'])
    for r_, b in zip(rows, base):
        yy = pl.y0 + b
        page.text(x_lab, yy, r_['label'], size=fs, ha='right', gid=f'{gid}-rowlabel')
        cells = []
        printed = r_.get('values') or {}             # optional table value -> printed text (e.g. 50x -> 50)
        for lane in sorted(lx):
            v = cell[r_['column']][lane].strip()
            v = str(printed.get(v, v))
            v = MINUS if v in ('-', MINUS) else v
            size = pm_pt if v in (MINUS, '+') else (cfg.get('word_pt', fs) if len(v) == 1 else cfg.get('amount_pt', fs))
            cells.append((lx[lane], page.text_width(v, size), v, size))
        pinned = [i for i, c in enumerate(cells) if c[2] in (MINUS, '+') or len(c[2]) == 1]
        xs = spread([(c[0], c[1]) for c in cells], min_gap=cfg['min_label_gap'], pinned=pinned)
        for (ax_, wd, v, size), xc in zip(cells, xs):
            page.text(xc, yy, v, size=size, ha='center', gid=f'{gid}-cells')
    # group labels over the lanes (Probe / Mutated Probe), each with a line under it across its lanes
    y_grp = pl.y0 + cfg['group_baseline']
    pad = cfg['group_line_pad']
    for gname, gl in cfg['groups'].items():
        ln = sorted(int(v) for v in lanes[lanes['group'] == gname]['lane'])
        if not ln:
            raise ValidationError([f'fig6e_emsa_lanes.csv: no lanes in group {gname!r} (config e.groups)'])
        xa, xb = lx[ln[0]] - pad, lx[ln[-1]] + pad
        page.text((xa + xb) / 2, y_grp, gl, size=fs, ha='center', gid=f'{gid}-groups')
        mm.line(xa, y_grp + cfg['group_line_gap'], xb, y_grp + cfg['group_line_gap'], lw=cfg['group_line_pt'],
                colour='k', gid=f'{gid}-groups', z=3)
    # labels beside gel rows (e.g. the free probe), lane numbers below
    for lab, ypx in cfg['row_labels'].items():
        yy = pl.y0 + (float(ypx) - cy0) * pl.mm_per_px
        page.text(x_lab, yy + cap_mm(fs) / 2, lab, size=fs, ha='right', gid=f'{gid}-bandlabel')
    yl = pl.y0 + pl.h_mm + cfg['lane_number_gap']
    page.text(x_lab, yl, cfg['lane_label'], size=fs, ha='right', gid=f'{gid}-lanes')
    for lane, xx in lx.items():
        page.text(xx, yl, str(lane), size=fs, ha='center', gid=f'{gid}-lanes')
    top = y_grp - cap_mm(fs)
    return dict(type='gel', lanes={k: round(v, 3) for k, v in lx.items()}, gel_mm=[round(pl.w_mm, 2), round(pl.h_mm, 2)],
                gel_ppi=pl.ppi, top_mm=round(top, 3), bottom_mm=round(yl, 3)), pl
