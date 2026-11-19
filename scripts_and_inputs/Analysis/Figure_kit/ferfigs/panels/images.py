"""Image tiles from an image manifest (one row per tile; see schemas.SCHEMAS['image_manifest']).

The kit never changes pixel intensities. It only orients (dihedral transforms, no resampling), crops, places at
one scale, and draws vector overlays (boxes, arrows) and vector scale bars on top. Swapping in a new image =
replacing the file and, if needed, its crop_px / um_per_px in the manifest.
"""
from __future__ import annotations

import math
from pathlib import Path

from .. import style as S


def parse_crop(s):
    if s is None or str(s).strip() == '':
        return None
    v = [float(t) for t in str(s).replace(',', ' ').split()]
    if len(v) != 4:
        raise ValueError(f"crop_px '{s}' must be four numbers: x0 y0 x1 y1")
    return tuple(v)


def parse_overlays(s):
    """'box 10 20 60 80; arrow 30 40 45 40; text 5 5 label' -> [(kind, [numbers], text)] in image pixels."""
    out = []
    if s is None or str(s).strip() == '':
        return out
    for item in str(s).split(';'):
        tok = item.split()
        if not tok:
            continue
        kind = tok[0]
        nums, rest = [], []
        for t in tok[1:]:
            try:
                nums.append(float(t))
            except ValueError:
                rest.append(t)
        out.append((kind, nums, ' '.join(rest)))
    return out


def place(page, fig_dir, row, x, y, w=None, h=None, gid=None):
    p = Path(fig_dir) / row['file']
    return page.image(p, x, y, w=w, h=h, crop=parse_crop(row.get('crop_px')), orient=row.get('orient') or 'identity',
                      gid=gid or f"{page.name}-{row['panel']}-{row['tile']}", name=row['file'])


def origin_of(row):
    """Pixel origin to subtract from overlay / scale-bar coordinates (crop corner when px_origin = file)."""
    if (row.get('px_origin') or 'placed') == 'file':
        c = parse_crop(row.get('crop_px'))
        if c:
            return c[0], c[1]
    return 0.0, 0.0


def draw_overlays(page, placed, spec, colour=(1, 1, 1), lw=0.5, gid=None, origin=(0.0, 0.0)):
    ox, oy = origin
    for kind, n, txt in parse_overlays(spec):
        if kind in ('box', 'poly', 'arrow', 'text'):
            n = [v - (ox if i % 2 == 0 else oy) for i, v in enumerate(n)] if kind != 'text' else \
                [n[0] - ox, n[1] - oy] + n[2:]
        elif kind == 'arrowhead':
            n = [n[0] - ox, n[1] - oy] + n[2:]
        if kind == 'box':
            x0, y0 = placed.mm(n[0], n[1]); x1, y1 = placed.mm(n[2], n[3])
            page.rect(x0, y0, x1 - x0, y1 - y0, lw=lw, edge=colour, gid=gid, zorder=6)
        elif kind == 'poly':           # closed outline through the given pixel corners (e.g. a rotated box)
            pts = [placed.mm(n[i], n[i + 1]) for i in range(0, len(n) - 1, 2)]
            page.polyline(pts, lw=lw, colour=colour, closed=True, gid=gid, zorder=6)
        elif kind == 'arrow':          # arrowhead at (x1, y1), tail at (x0, y0)
            x0, y0 = placed.mm(n[0], n[1]); x1, y1 = placed.mm(n[2], n[3])
            page.line(x0, y0, x1, y1, lw=lw, colour=colour, gid=gid, zorder=6)
            ang = math.atan2(y1 - y0, x1 - x0); L = 0.9; a = math.radians(25)
            pts = [(x1, y1), (x1 - L * math.cos(ang - a), y1 - L * math.sin(ang - a)),
                   (x1 - L * math.cos(ang + a), y1 - L * math.sin(ang + a))]
            page.polyline(pts, lw=0, colour=None, fill=colour, closed=True, gid=gid, zorder=6)
        elif kind == 'arrowhead':      # filled triangle pointing in direction (dx, dy) at (x, y)
            cx, cy = placed.mm(n[0], n[1]); dx, dy = n[2], n[3]
            ang = math.atan2(dy, dx); L = n[4] if len(n) > 4 else 1.0; a = math.radians(28)
            pts = [(cx, cy), (cx - L * math.cos(ang - a), cy - L * math.sin(ang - a)),
                   (cx - L * math.cos(ang + a), cy - L * math.sin(ang + a))]
            page.polyline(pts, lw=0, colour=None, fill=colour, closed=True, gid=gid, zorder=6)
        elif kind == 'text':
            cx, cy = placed.mm(n[0], n[1])
            page.text(cx, cy, txt, colour=colour, gid=gid)
        else:
            raise ValueError(f"unknown overlay '{kind}' (use box, poly, arrow, arrowhead, text)")


def draw_scalebar(page, placed, row, colour=(1, 1, 1), lw=1.0, gid=None):
    ox, oy = origin_of(row)
    """Vector scale bar from the manifest: scale_bar_um long at um_per_px; at scale_bar_px ('x0 y', pixels of the
    placed image) or bottom right. Returns the bar length in mm, or None when the row has no scale bar."""
    um, upp = row.get('scale_bar_um'), row.get('um_per_px')
    if not um or not upp:
        return None
    um, upp = float(um), float(upp)
    Lpx = um / upp
    pos = row.get('scale_bar_px')
    if pos:
        x0, yc = [float(v) for v in str(pos).split()]
        x0, yc = x0 - ox, yc - oy
        xa, ya = placed.mm(x0, yc); xb, _ = placed.mm(x0 + Lpx, yc)
        page.line(xa, ya, xb, ya, lw=lw, colour=colour, gid=gid or f'{page.name}-scalebar', zorder=7)
        return xb - xa
    return page.scalebar(placed, um, upp, colour=colour, lw=lw, gid=gid)


def scale_for_width(rows, fig_dir, tile_w):
    """mm per micrometre when the widest tile (in um) is tile_w mm: keeps every tile of a panel at one
    magnification when um_per_px is given for each tile."""
    from PIL import Image
    widest = 0.0
    for _, r in rows.iterrows():
        c = parse_crop(r.get('crop_px'))
        if c:
            wpx = c[2] - c[0]
        else:
            with Image.open(Path(fig_dir) / r['file']) as im:
                wpx = im.size[0]
        widest = max(widest, wpx * float(r['um_per_px']))
    return tile_w / widest


def grid(page, fig_dir, rows, x, y, tile_w=None, tile_h=None, gap=None, col_x=None, row_y=None, gid=None,
         um_per_mm=None):
    """Place tiles by their row/col numbers. Every tile gets the same width (tile_w) or height (tile_h), or, with
    um_per_mm, one magnification for all tiles (width = width_px * um_per_px / um_per_mm; needs um_per_px per tile).
    col_x / row_y: optional explicit offsets (mm, relative to x / y) per column / row index.
    Returns {(row, col): Placed}."""
    gap = page.st['gaps_mm']['tile'] if gap is None else gap
    out = {}
    from PIL import Image
    # tile heights for row placement
    sizes = {}
    for _, r in rows.iterrows():
        c = parse_crop(r.get('crop_px'))
        if c:
            wpx, hpx = c[2] - c[0], c[3] - c[1]
        else:
            with Image.open(Path(fig_dir) / r['file']) as im:
                wpx, hpx = im.size
            if (r.get('orient') or 'identity') in ('rot90ccw', 'rot270ccw', 'mirrorLR+rot90ccw', 'mirrorLR+rot270ccw'):
                wpx, hpx = hpx, wpx
        sizes[(r['row'], r['col'])] = (wpx, hpx)
    upp = {(r['row'], r['col']): float(r['um_per_px']) for _, r in rows.iterrows() if r.get('um_per_px')}
    def wh(key):
        wpx, hpx = sizes[key]
        if um_per_mm is not None:
            if key not in upp:
                raise ValueError(f'tile {key}: um_per_px is needed for one-magnification placement (um_per_mm)')
            w = wpx * upp[key] / um_per_mm
            return w, w * hpx / wpx
        if tile_w is not None:
            return tile_w, tile_w * hpx / wpx
        return tile_h * wpx / hpx, tile_h
    rows_idx = sorted({k[0] for k in sizes}); cols_idx = sorted({k[1] for k in sizes})
    row_h = {ri: max(wh(k)[1] for k in sizes if k[0] == ri) for ri in rows_idx}
    col_w = {ci: max(wh(k)[0] for k in sizes if k[1] == ci) for ci in cols_idx}
    xs, acc = {}, 0.0
    for ci in cols_idx:
        xs[ci] = col_x[ci] if col_x and ci in col_x else acc
        acc = xs[ci] + col_w[ci] + gap
    ys, acc = {}, 0.0
    for ri in rows_idx:
        ys[ri] = row_y[ri] if row_y and ri in row_y else acc
        acc = ys[ri] + row_h[ri] + gap
    for _, r in rows.iterrows():
        k = (r['row'], r['col'])
        tw, th = wh(k)
        pl = place(page, fig_dir, r, x + xs[k[1]], y + ys[k[0]], w=tw, gid=gid)
        out[k] = pl
        if r.get('overlays'):
            draw_overlays(page, pl, r['overlays'], gid=gid, origin=origin_of(r))
    return out
