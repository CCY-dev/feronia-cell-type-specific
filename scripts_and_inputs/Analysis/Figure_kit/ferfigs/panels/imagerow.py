"""A row of image tiles at one magnification, with the header structure used across the FER figures:

    reporter header ----------------------------------   (optional, bracket over the whole row)
    group header ----------   group header ----------    (optional, bracket over each group of tiles)
    RALF23      -        +        -        +             (optional treatment row: true minus and plus)
              [tile]  [tile]   [tile]  [tile]  [key]

Tiles come from the image manifest (panel rows ordered by col). One mm_per_px for the whole row keeps every tile
at the same magnification when the tiles share a pixel size (um_per_px); the vector scale bar is drawn from the
manifest. Gaps: `within` between tiles of a group, `between` between groups (both >= the journal's 0.5-mm
demarcation for images of different samples).
"""
from __future__ import annotations

from .. import style as S
from . import images as I

MINUS = '\u2212'


def tile_px(fig_dir, r):
    from PIL import Image
    c = I.parse_crop(r.get('crop_px'))
    if c:
        return c[2] - c[0], c[3] - c[1]
    with Image.open(fig_dir / r['file']) as im:
        w, h = im.size
    if (r.get('orient') or 'identity') in ('rot90ccw', 'rot270ccw', 'mirrorLR+rot90ccw', 'mirrorLR+rot270ccw'):
        w, h = h, w
    return w, h


def row_width(fig_dir, rows, mm_per_px=None, tile_w=None, groups=None, within=0.6, between=1.6):
    xs, x = [], 0.0
    cols = list(rows['col'])
    gid = {c: gi for gi, (_, cs) in enumerate(groups or [('', cols)]) for c in cs}
    prev = None
    for _, r in rows.iterrows():
        if prev is not None:
            x += within if gid.get(r['col']) == gid.get(prev) else between
        w = tile_w if tile_w else tile_px(fig_dir, r)[0] * mm_per_px
        xs.append((x, w))
        x += w
        prev = r['col']
    return x, xs


def draw(page, fig_dir, rows, x, y_tiles, mm_per_px=None, tile_w=None, groups=None, within=0.6, between=1.6,
         reporter=None, pm=None, pm_label='RALF23', key=None, labels_above=None, overlay_colour=(1, 1, 1),
         gid='row', y_reporter=None, y_groups=None, y_pm=None):
    """rows: manifest rows of one panel (sorted by col). groups: [(markup, [cols])]. pm: {col: '-'|'+'}.
    key: dict(file=..., w=2.2, h=12, labels=('Max', 'Min')) drawn right of the last tile, bottom-aligned.
    labels_above: {col: markup} single labels centred over tiles (instead of group headers).
    y_* : baselines (mm) for the header lines; defaults stack upwards from the tile tops.
    Returns dict(placed={col: Placed}, x0, x1, y_top, y_bottom)."""
    fig_dir = fig_dir
    rows = rows.sort_values('col')
    st = page.st
    size = st['font']['size_pt']
    width, xs = row_width(fig_dir, rows, mm_per_px, tile_w, groups, within, between)
    placed = {}
    for (dx, w), (_, r) in zip(xs, rows.iterrows()):
        pl = I.place(page, fig_dir, r, x + dx, y_tiles, w=w, gid=f'{gid}-tile')
        placed[r['col']] = pl
        if r.get('overlays'):
            I.draw_overlays(page, pl, r['overlays'], colour=overlay_colour, gid=f'{gid}-overlay', origin=I.origin_of(r))
        I.draw_scalebar(page, pl, r, gid=f'{gid}-scalebar')
    x1 = x + width
    bottom = max(p.y0 + p.h_mm for p in placed.values())
    # header lines, stacked upwards: treatment row, group headers, reporter header
    yp = y_pm if y_pm is not None else y_tiles - 1.0
    if pm:
        for c, sgn in pm.items():
            p = placed[c]
            page.text(p.x0 + p.w_mm / 2, yp, MINUS if sgn in ('-', MINUS) else '+', size=7.0, ha='center',
                      gid=f'{gid}-pm')
        page.text(x - 1.2, yp, pm_label, ha='right', gid=f'{gid}-pm')
    yg = y_groups if y_groups is not None else (yp - 3.3 if pm else y_tiles - 1.2)
    if groups:
        for lab, cs in groups:
            a = placed[cs[0]]; b = placed[cs[-1]]
            page.header(a.x0, b.x0 + b.w_mm, yg, lab, gid=f'{gid}-header')
    if labels_above:
        for c, lab in labels_above.items():
            p = placed[c]
            page.text(p.x0 + p.w_mm / 2, yg, lab, ha='center', gid=f'{gid}-header')
    yr = y_reporter if y_reporter is not None else yg - 3.3
    if reporter:
        page.header(x, x1, yr, reporter, gid=f'{gid}-header')
    if key:
        kx = x1 + key.get('gap', 0.8)
        kh, kw = key['h'], key['w']
        kp = page.image(fig_dir / key['file'], kx, bottom - kh, h=kh, gid=f'{gid}-key', name=key['file'])
        page.rect(kx, bottom - kh, kp.w_mm, kh, lw=0.25, edge='k', gid=f'{gid}-key')
        hi, lo = key.get('labels', ('Max', 'Min'))
        page.text(kx + kp.w_mm + 0.5, bottom - kh + 0.36 * size / S.PT_PER_MM * 2, hi, gid=f'{gid}-key')
        page.text(kx + kp.w_mm + 0.5, bottom, lo, gid=f'{gid}-key')
        x1 = kx + kp.w_mm + 0.5 + page.text_width(hi)
    base = [v for v in (yr if reporter else None, yg if (groups or labels_above) else None, yp if pm else None)
            if v is not None]
    top = min(base) - 0.716 * size / S.PT_PER_MM if base else y_tiles
    return dict(placed=placed, x0=x, x1=x1, y_top=top, y_bottom=bottom)
