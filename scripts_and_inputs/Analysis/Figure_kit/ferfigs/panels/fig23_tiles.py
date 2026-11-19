"""Image tiles at one magnification, tile labels, overlay marks, chart text size and chart keys (Figs 2, 3).

Shared by figures/fig2.py and figures/fig3.py; built only on the core primitives (Page, images.place,
images.draw_overlays, images.draw_scalebar, imagerow.tile_px). Candidate for promotion into the core
(docs/notes/fig2_NOTES.md and fig3_NOTES.md, core change requests).

Magnification: every tile of a panel is drawn at one scale in micrometres per millimetre (um_per_mm), so a tile's
printed width is  width_px * um_per_px / um_per_mm  and one scale bar is valid for the whole panel even when the
tiles have different pixel sizes (Fig. 2g). A tile listed in fit_row_height is instead scaled to the height of the
other tiles of its row (as the published layout shows the 965-px Fig. 2b tile); the build record lists
it. Pixels are never changed (orientation and crop only).

Extra manifest overlay kinds handled here (anything else goes to the core images.draw_overlays, which also accepts
an optional colour word at the end of an item, e.g. 'box 10 20 60 80 white'):
  openhead x y dx dy <colour>     open triangle with its tip at pixel (x, y) pointing along (dx, dy); size and
                                  colour from the figure config (marks: {<colour>: {length_mm, width_mm, rgb}})
  mask x0 y0 x1 y1 <colour>       filled rectangle (covers a burned-in label; listed in notes)
  vbar x y_top                    vertical scale bar from pixel (x, y_top) downwards, scale_bar_um long
"""
from __future__ import annotations

import copy
import math
from contextlib import contextmanager

import matplotlib as mpl
from matplotlib.lines import Line2D

from .. import style as S
from ..layout import ASCENDER
from . import images as I
from .imagerow import tile_px

DESC = 0.21          # Arial descender (fraction of the font size), for text below images


def cap_mm(size):
    return ASCENDER * size / S.PT_PER_MM


def line_pitch(size):
    return 1.18 * size / S.PT_PER_MM


# ------------------------------------------------------------------ placement
def scale_mm_per_px(row, um_per_mm):
    upp = row.get('um_per_px')
    if not upp:
        raise ValueError(f"tile {row['panel']}/{row['tile']}: um_per_px is empty; one magnification needs it")
    return float(upp) / float(um_per_mm)


def grid_geometry(fig_dir, rows, um_per_mm, col_gap, row_gap, fit_row_height=()):
    """Tile sizes, column widths and row heights (mm) of a manifest grid at one magnification. Tiles named in
    fit_row_height are scaled to the tallest other tile of their row (aspect kept)."""
    size = {}
    for _, r in rows.iterrows():
        wpx, hpx = tile_px(fig_dir, r)
        s = scale_mm_per_px(r, um_per_mm)
        size[(r['row'], r['col'])] = [wpx * s, hpx * s, r['tile']]
    for k, (w, h, t) in list(size.items()):
        if t in fit_row_height:
            others = [v[1] for kk, v in size.items() if kk[0] == k[0] and v[2] not in fit_row_height]
            if not others:
                raise ValueError(f'fit_row_height: tile {t} has no other tile in its row to take the height from')
            H = max(others)
            size[k] = [w * H / h, H, t]
    ri = sorted({k[0] for k in size}); ci = sorted({k[1] for k in size})
    col_w = {c: max(size[k][0] for k in size if k[1] == c) for c in ci}
    row_h = {r: max(size[k][1] for k in size if k[0] == r) for r in ri}
    width = sum(col_w.values()) + col_gap * (len(ci) - 1)
    height = sum(row_h.values()) + row_gap * (len(ri) - 1)
    return dict(size=size, col_w=col_w, row_h=row_h, width=width, height=height, rows=ri, cols=ci)


def place_grid(page, fig_dir, rows, x, y, um_per_mm, col_gap=0.6, row_gap=0.6, gid='grid', valign='top',
               fit_row_height=()):
    """Place tiles by their manifest row/col at one magnification, left-aligned in their column and top- (or
    bottom-) aligned in their row. Returns (placed {(row, col): Placed}, geometry)."""
    g = grid_geometry(fig_dir, rows, um_per_mm, col_gap, row_gap, fit_row_height)
    xs, acc = {}, x
    for c in g['cols']:
        xs[c] = acc; acc += g['col_w'][c] + col_gap
    ys, acc = {}, y
    for r in g['rows']:
        ys[r] = acc; acc += g['row_h'][r] + row_gap
    placed = {}
    for _, r in rows.sort_values(['row', 'col']).iterrows():
        k = (r['row'], r['col'])
        w, h, _ = g['size'][k]
        yy = ys[k[0]] + (g['row_h'][k[0]] - h if valign == 'bottom' else 0.0)
        placed[k] = I.place(page, fig_dir, r, xs[k[1]], yy, w=w, gid=f'{gid}-tile')
    return placed, g


def boxes(rows, placed):
    """{tile: [x, y, w, h]} and the union box (mm) of placed tiles, for the build record and layout checks."""
    tb = {}
    for _, r in rows.sort_values(['row', 'col']).iterrows():
        p = placed[(r['row'], r['col'])]
        tb[r['tile']] = [round(p.x0, 3), round(p.y0, 3), round(p.w_mm, 3), round(p.h_mm, 3)]
    x0 = min(b[0] for b in tb.values()); y0 = min(b[1] for b in tb.values())
    x1 = max(b[0] + b[2] for b in tb.values()); y1 = max(b[1] + b[3] for b in tb.values())
    return tb, [round(x0, 3), round(y0, 3), round(x1 - x0, 3), round(y1 - y0, 3)]


# ------------------------------------------------------------------ marks on images
def _colour(name, marks, default):
    if name in (marks or {}):
        return S.rgb(marks[name]['rgb'])
    return {'white': (1.0, 1.0, 1.0), 'black': (0.0, 0.0, 0.0)}.get(name, default)


def open_head(page, x, y, dx, dy, length, width, colour, lw, gid=None):
    """Open (outlined) triangle with its tip at (x, y) mm, pointing along (dx, dy)."""
    n = math.hypot(dx, dy) or 1.0
    ux, uy = dx / n, dy / n
    bx, by = x - ux * length, y - uy * length          # base midpoint
    px, py = -uy * width / 2, ux * width / 2
    page.polyline([(x, y), (bx + px, by + py), (bx - px, by - py)], lw=lw, colour=colour, closed=True, gid=gid,
                  zorder=6)


def draw_marks(page, placed, row, marks=None, bar_colour=(1, 1, 1), bar_lw=1.0, mark_lw=0.6, gid='marks'):
    """Overlays of one manifest row (custom kinds here, the rest through the core) and its scale bar."""
    upp = float(row['um_per_px']) if row.get('um_per_px') else None
    for item in str(row.get('overlays') or '').split(';'):
        tok = item.split()
        if not tok:
            continue
        kind = tok[0]
        nums = [float(t) for t in tok[1:] if _isnum(t)]
        words = [t for t in tok[1:] if not _isnum(t)]
        col = _colour(words[0], marks, None) if words else None
        if kind == 'openhead':
            spec = (marks or {}).get(words[0], {}) if words else {}
            x, y = placed.mm(nums[0], nums[1])
            open_head(page, x, y, nums[2], nums[3], spec.get('length_mm', 1.0), spec.get('width_mm', 1.15),
                      col or (0, 0, 0), spec.get('lw', mark_lw), gid=f'{gid}-arrowhead')
        elif kind == 'mask':
            x0, y0 = placed.mm(nums[0], nums[1]); x1, y1 = placed.mm(nums[2], nums[3])
            page.rect(x0, y0, x1 - x0, y1 - y0, lw=0, edge='none', face=col or (1, 1, 1), gid=f'{gid}-mask',
                      zorder=5)
        elif kind == 'vbar':
            if not (row.get('scale_bar_um') and upp):
                raise ValueError(f"tile {row['panel']}/{row['tile']}: vbar needs scale_bar_um and um_per_px")
            L = float(row['scale_bar_um']) / upp
            xa, ya = placed.mm(nums[0], nums[1]); _, yb = placed.mm(nums[0], nums[1] + L)
            page.line(xa, ya, xa, yb, lw=bar_lw, colour=bar_colour, gid=f'{gid}-scalebar', zorder=7)
        else:
            I.draw_overlays(page, placed, kind + ' ' + ' '.join(t for t in tok[1:] if _isnum(t) or kind == 'text'),
                            colour=col or (1, 1, 1), lw=mark_lw, gid=f'{gid}-overlay')
    if row.get('scale_bar_px'):
        I.draw_scalebar(page, placed, row, colour=bar_colour, lw=bar_lw, gid=f'{gid}-scalebar')


def _isnum(t):
    try:
        float(t)
        return True
    except ValueError:
        return False


def frames(page, placed, lw=0.5, colour=(0, 0, 0), gid='frame'):
    """Thin outline on the edge of each tile (white-background images need it to read as separate images)."""
    for p in placed.values():
        page.rect(p.x0, p.y0, p.w_mm, p.h_mm, lw=lw, edge=colour, gid=gid, zorder=6)


# ------------------------------------------------------------------ labels
def labels_above(page, centres, lines_by_key, y_base, size=None, pitch=None, gid='labels'):
    """Labels centred at centres[k]; lines_by_key[k] = ['line 1', 'line 2']; the LAST line sits on y_base
    (bottom-aligned stack, line pitch in mm). Returns the y of the top of the tallest stack."""
    size = size or page.st['font']['size_pt']
    pitch = pitch or line_pitch(size)
    top = y_base
    for k, lines in lines_by_key.items():
        n = len(lines)
        for i, s in enumerate(lines):
            yb = y_base - (n - 1 - i) * pitch
            page.text(centres[k], yb, s, size=size, ha='center', gid=gid)
        top = min(top, y_base - (n - 1) * pitch - cap_mm(size))
    return top


def labels_below(page, centres, lines_by_key, y_base, size=None, pitch=None, gid='labels'):
    """Labels centred at centres[k], the FIRST line on baseline y_base. Returns the lowest descender (mm)."""
    size = size or page.st['font']['size_pt']
    pitch = pitch or line_pitch(size)
    low = y_base
    for k, lines in lines_by_key.items():
        for i, s in enumerate(lines):
            yb = y_base + i * pitch
            page.text(centres[k], yb, s, size=size, ha='center', gid=gid)
            low = max(low, yb + DESC * size / S.PT_PER_MM)
    return low


def max_label_width(page, labels, size=None):
    return max(page.text_width(s, size) for lines in labels for s in lines)


def check_label_spacing(page, centres, lines_by_key, min_gap=0.6, size=None, where=''):
    """Neighbouring labels (same line index) must be at least min_gap mm apart; returns a list of problems."""
    keys = sorted(centres, key=lambda k: centres[k])
    probs = []
    for a, b in zip(keys, keys[1:]):
        for i in range(max(len(lines_by_key[a]), len(lines_by_key[b]))):
            la = lines_by_key[a][i] if i < len(lines_by_key[a]) else ''
            lb = lines_by_key[b][i] if i < len(lines_by_key[b]) else ''
            if not la or not lb:
                continue
            gap = (centres[b] - page.text_width(lb, size) / 2) - (centres[a] + page.text_width(la, size) / 2)
            if gap < min_gap:
                probs.append(f'{where}: labels {a!r} and {b!r} are {gap:.2f} mm apart (< {min_gap} mm); enlarge the '
                             f'tiles or split the label differently in the config')
    return probs


# ------------------------------------------------------------------ charts
def xlabel_depth(page, spec):
    """Distance (mm) from the bottom of a dotbar frame to the lowest point of its x labels (as dotbar draws them:
    tick, 0.6-mm gap, labels rotated by x_rotation and ending at their tick, or horizontal under it)."""
    size = page.st['font']['size_pt']
    tick = page.st['lines']['tick_len_pt'] / S.PT_PER_MM
    rot = math.radians(spec.get('x_rotation', 45))
    asc, desc = cap_mm(size), DESC * size / S.PT_PER_MM
    w = max(page.text_width(S.display(g, spec.get('label_kind', 'genotypes')), size) for g in spec['groups'])
    if rot:
        return tick + 0.6 + w * math.sin(rot) + (asc + desc) * math.cos(rot)
    return tick + 0.5 + asc + desc


@contextmanager
def chart_text(page, size_pt):
    """Draw a chart (dotbar) with a text size other than the house default, e.g. the approved layout's chart text
    capped at 7 pt: the page uses a copy of its style with font.size_pt = size_pt for the duration, and matplotlib
    creates the tick labels at that size. Nothing else of the style changes."""
    old = page.st
    st = copy.deepcopy(old)
    st['font']['size_pt'] = float(size_pt)
    page.st = st
    try:
        with mpl.rc_context({'xtick.labelsize': float(size_pt), 'ytick.labelsize': float(size_pt)}):
            yield
    finally:
        page.st = old


# ------------------------------------------------------------------ chart key (series markers)
def key(page, x_centre, y_base, items, title=None, gid='key', marker_gap=0.9, item_gap=2.4, title_gap=2.4,
        size=None):
    """Series key in the style of dotbar's legend (open circle + label), optionally preceded by a title naming
    the treatment (e.g. 'RALF23'). items: [(markup, rgb)]. Centred on x_centre; text baseline y_base.
    Returns (x0, x1) of the key."""
    st = page.st
    size = size or st['font']['size_pt']
    msize = st['markers']['diameter_pt']
    mk = msize / S.PT_PER_MM
    widths = [page.text_width(t, size) for t, _ in items]
    tw = page.text_width(title, size) if title else 0.0
    total = sum(widths) + len(items) * (mk + marker_gap) + item_gap * (len(items) - 1) + (tw + title_gap if title else 0)
    cx = x_centre - total / 2
    x0 = cx
    if title:
        page.text(cx, y_base, title, size=size, gid=gid)
        cx += tw + title_gap
    for (t, c), w in zip(items, widths):
        page.fig.add_artist(Line2D([page.fx(cx + mk / 2)], [page.fy(y_base - 0.36 * size / S.PT_PER_MM)],
                                   transform=page.fig.transFigure, marker='o', markersize=msize * 1.25,
                                   markerfacecolor='none', markeredgecolor=c,
                                   markeredgewidth=st['lines']['marker_edge_pt'] * 1.2, linestyle='none', gid=gid))
        page.text(cx + mk + marker_gap, y_base, t, size=size, gid=gid)
        cx += mk + marker_gap + w + item_gap
    return x0, cx - item_gap


def key_at(page, y_base, items, title=None, title_x=None, marker_gap=0.9, size=None, gid='key'):
    """Series key with every element at a given position (e.g. measured on an approved page): items =
    [(markup, rgb, x_marker_left_mm)]; the label starts marker_gap mm after the marker; title (e.g. 'RALF23') starts
    at title_x. Returns (x0, x1) of the key."""
    st = page.st
    size = size or st['font']['size_pt']
    msize = st['markers']['diameter_pt']
    mk = msize / S.PT_PER_MM
    x0 = title_x if title else min(x for _, _, x in items)
    if title:
        page.text(title_x, y_base, title, size=size, gid=gid)
    x1 = x0
    for t, c, xm in items:
        page.fig.add_artist(Line2D([page.fx(xm + mk / 2)], [page.fy(y_base - 0.36 * size / S.PT_PER_MM)],
                                   transform=page.fig.transFigure, marker='o', markersize=msize * 1.25,
                                   markerfacecolor='none', markeredgecolor=c,
                                   markeredgewidth=st['lines']['marker_edge_pt'] * 1.2, linestyle='none', gid=gid))
        w = page.text(xm + mk + marker_gap, y_base, t, size=size, gid=gid)
        x1 = max(x1, xm + mk + marker_gap + w)
    return x0, x1
