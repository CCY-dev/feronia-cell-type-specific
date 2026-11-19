"""Fig. 6f: RALF23 seedling tiles (one genotype per tile: mock seedlings left, RALF23 right), relabelled as in
draw_seedlings_f.py. Tiles are crops of the native image (manifest panel f; pixels unchanged) placed at one
scale with the journal's demarcation gap. Genotype names above the top row and below the bottom row, the -/+ signs
over the two seedling groups of each top-row tile (x from data/fig6/fig6f_seedlings.csv), RALF23 on the -/+
baseline left of the tiles, and the printed root-length ratio (same table) in white at the lower right of each tile.
The 1-cm bar burned into one tile is covered exactly by a vector bar from the manifest (scale_bar_px).
"""
from __future__ import annotations

from .. import style as S
from ..schemas import ValidationError
from . import images as I
from .fig6_common import cap_mm, display

MINUS = '\u2212'


def draw(page, fig_dir, rows, table, cfg, x, y_tiles, tile_w, gap, gid='fig6-f'):
    st = page.st
    fs = cfg.get('text_pt', st['font']['size_pt'])
    pm_pt = cfg.get('sign_pt', 7.0)
    rows = rows.sort_values(['row', 'col'])
    tab = {r['tile']: r for _, r in table.iterrows()}
    miss = [t for t in rows['tile'] if t not in tab]
    if miss:
        raise ValidationError([f'fig6f_seedlings.csv: no row for tile(s) {miss}'])
    placed = {}
    ymax_row = {}
    for _, r in rows.iterrows():
        c = I.parse_crop(r.get('crop_px'))
        wpx, hpx = c[2] - c[0], c[3] - c[1]
        th = tile_w * hpx / wpx
        ymax_row.setdefault(int(r['row']), th)
    y_rows, acc = {}, y_tiles
    for ri in sorted(ymax_row):
        y_rows[ri] = acc
        acc += ymax_row[ri] + gap
    for _, r in rows.iterrows():
        ri, ci = int(r['row']), int(r['col'])
        xx = x + ci * (tile_w + gap)
        pl = I.place(page, fig_dir, r, xx, y_rows[ri], w=tile_w, gid=f'{gid}-tile')
        placed[(ri, ci)] = (pl, r)
        I.draw_scalebar(page, pl, r, colour=(1, 1, 1), lw=cfg['scale_bar_pt'], gid=f'{gid}-scalebar')
        t = tab[r['tile']]
        page.text(pl.x0 + pl.w_mm - cfg['ratio_inset'][0], pl.y0 + pl.h_mm - cfg['ratio_inset'][1], str(t['ratio']),
                  size=fs, ha='right', colour=(1, 1, 1), gid=f'{gid}-ratio')
    top_row, bottom_row = min(y_rows), max(y_rows)
    y_pm = y_tiles - cfg['pm_above']
    y_name_top = y_pm - cfg['name_above_pm']
    for (ri, ci), (pl, r) in placed.items():
        lab = display(r.get('genotype') or r['tile'], cfg.get('display_labels'))
        if ri == top_row:
            page.text(pl.x0 + pl.w_mm / 2, y_name_top, lab, size=fs, ha='center', gid=f'{gid}-genotype')
            t = tab[r['tile']]
            for key, sgn in (('minus_x_px', MINUS), ('plus_x_px', '+')):
                page.text(pl.x0 + float(t[key]) * pl.mm_per_px, y_pm, sgn, size=pm_pt, ha='center', gid=f'{gid}-pm')
        if ri == bottom_row:
            page.text(pl.x0 + pl.w_mm / 2, pl.y0 + pl.h_mm + cfg['name_below'] + cap_mm(fs), lab, size=fs,
                      ha='center', gid=f'{gid}-genotype')
    page.text(x - cfg.get('pm_label_gap', 1.2), y_pm, cfg['pm_label'], size=fs, ha='right', gid=f'{gid}-pm')
    pls = [p for p, _ in placed.values()]
    x1 = max(p.x0 + p.w_mm for p in pls)
    y1 = max(p.y0 + p.h_mm for p in pls)
    return dict(type='images', tiles=[r['tile'] for _, r in rows.iterrows()],
                ratios={r['tile']: str(tab[r['tile']]['ratio']) for _, r in rows.iterrows()},
                tile_mm=[round(tile_w, 3), round(pls[0].h_mm, 3)], ppi=pls[0].ppi,
                top_mm=round(y_name_top - cap_mm(fs), 3), bottom_mm=round(y1 + cfg['name_below'] + cap_mm(fs), 3),
                right_mm=round(x1, 3))
