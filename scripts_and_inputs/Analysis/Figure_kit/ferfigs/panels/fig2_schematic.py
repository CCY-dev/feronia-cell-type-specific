"""Fig. 2a: how the cell type-specific FER knockout (TSKO) lines were made, drawn as vector line art from words in
the config (config/fig2.yaml, schematic:). Box colours come from the genotype palette in config/style.yaml
(colours.genotype_fill / genotype_stroke), so the boxes match the chart colours in d, f and h.

    [ fer-4 mutant ]                      top box
          |  Complementation
    [ pFER:FER-GFP/fer-4 ]                middle box
          |  Cell type-specific CRISPR
          |  pAtU6-FERgRNA
    ------+------------------------       bar
    |       |       |       |       |     one arrow per line, promoter label on the arrow
  [KO]    [KO]    [KO]    [KO]    [KO]    tissue boxes (name + target tissue)

Geometry (mm, relative to the panel top) and every word come from the config; nothing here is data. The defaults
reproduce the published schematic: 0.7-pt dark grey lines, 1.5 x 1.26-mm arrowheads, 0.6-pt box
outlines, corner radius 0.18 x the shorter box side, box text centred as on the approved page.
"""
from __future__ import annotations

from matplotlib.path import Path
from matplotlib.patches import PathPatch

from .. import style as S
from .fig23_tiles import line_pitch

KAPPA = 0.5523


def tint(rgb, f):
    """Colour at fraction f over white (f = 1: the colour itself)."""
    r = S.rgb(rgb)
    return tuple(c * f + (1 - f) for c in r)


def rounded_box(page, x, y, w, h, r, lw, edge, face, gid=None, zorder=3):
    r = min(r, w / 2, h / 2)
    k = KAPPA * r
    P = [(x + r, y), (x + w - r, y), (x + w - r + k, y), (x + w, y + r - k), (x + w, y + r),
         (x + w, y + h - r), (x + w, y + h - r + k), (x + w - r + k, y + h), (x + w - r, y + h),
         (x + r, y + h), (x + r - k, y + h), (x, y + h - r + k), (x, y + h - r),
         (x, y + r), (x, y + r - k), (x + r - k, y), (x + r, y), (x + r, y)]
    C = [Path.MOVETO, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4,
         Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4,
         Path.CURVE4, Path.CLOSEPOLY]
    path = Path([(page.fx(a), page.fy(b)) for a, b in P], C)
    patch = PathPatch(path, transform=page.fig.transFigure, lw=lw, edgecolor=edge, facecolor=face,
                      joinstyle='round', zorder=zorder)
    if gid:
        patch.set_gid(gid)
    page.fig.add_artist(patch)
    return patch


def text_block(page, cx, cy, lines, size, gid=None):
    """Lines centred on (cx, cy) as in the approved schematic: first baseline = cy - n * pitch / 2 + 0.82 * size
    (pitch 1.18 x size)."""
    n = len(lines)
    pitch = line_pitch(size)
    y = cy - n * pitch / 2 + 0.82 * size / S.PT_PER_MM
    for i, s in enumerate(lines):
        page.text(cx, y + i * pitch, s, size=size, ha='center', gid=gid)


def arrow_down(page, x, y0, y1, lw, colour, head, gid=None, gap=None):
    """Vertical arrow from y0 down to y1 (tip); head = (length, width) mm. gap=(ya, yb): the line is left open
    between ya and yb (a label sits there, as the white box behind the label on the approved page)."""
    L, W = head
    if gap:
        page.line(x, y0, x, gap[0], lw=lw, colour=colour, gid=gid)
        page.line(x, gap[1], x, y1 - L * 0.9, lw=lw, colour=colour, gid=gid)
    else:
        page.line(x, y0, x, y1 - L * 0.9, lw=lw, colour=colour, gid=gid)
    page.polyline([(x, y1), (x - W / 2, y1 - L), (x + W / 2, y1 - L)], lw=0, colour=None, fill=colour, closed=True,
                  gid=gid)


def draw(page, x, y, w, spec, gid='schematic'):
    """Schematic with its top-left at (x, y) mm and width w; every vertical position comes from spec['geometry']
    (mm below y). Returns the boxes drawn (mm) for the build record."""
    size = spec.get('text_pt', page.st['font']['size_pt'])
    psize = spec.get('promoter_pt', size)
    bsize = spec.get('box_pt', size)
    lw = spec.get('line_pt', 0.7)
    lc = S.rgb(spec.get('line_colour', [64, 64, 64]))
    head = spec.get('head_mm', [1.5, 1.26])
    blw = spec.get('box_lw', 0.6)
    frac = spec.get('corner_frac', 0.18)
    pal = page.st['colours']
    g = spec['geometry']
    cx = x + w / 2
    drawn = []

    def node(bx, by, bw, bh, lines, geno, t, tsize):
        face = tint(pal['genotype_fill'][geno], t)
        edge = S.rgb(pal['genotype_stroke'][geno])
        rounded_box(page, bx, by, bw, bh, frac * min(bw, bh), blw, edge, face, gid=f'{gid}-box')
        text_block(page, bx + bw / 2, by + bh / 2, lines, tsize, gid=f'{gid}-text')
        drawn.append([round(bx, 3), round(by, 3), round(bw, 3), round(bh, 3)])

    top, mid = spec['top'], spec['middle']
    t0, t1 = y + g['top_box'][0], y + g['top_box'][1]
    m0, m1 = y + g['middle_box'][0], y + g['middle_box'][1]
    node(cx - top['w'] / 2, t0, top['w'], t1 - t0, top['text'], top['genotype'], top.get('tint', 0.3), size)
    node(cx - mid['w'] / 2, m0, mid['w'], m1 - m0, mid['text'], mid['genotype'], mid.get('tint', 0.5), size)
    arrow_down(page, cx, t1, m0, lw, lc, head, gid=f'{gid}-arrow')
    lab_dx = spec.get('label_dx_mm', 1.6)
    for key in ('edge1', 'edge2'):
        e = spec[key]
        for t, yb in zip(e['text'], e['baselines']):
            page.text(cx + lab_dx, y + yb, t, size=size, gid=f'{gid}-text')
    ybar = y + g['bar']
    page.line(cx, m1, cx, ybar, lw=lw, colour=lc, gid=f'{gid}-arrow')
    br = spec['branches']
    n = len(br)
    gap = spec.get('box_gap_mm', 1.4)
    bw = (w - gap * (n - 1)) / n
    xs = [x + i * (bw + gap) + bw / 2 for i in range(n)]
    page.line(xs[0], ybar, xs[-1], ybar, lw=lw, colour=lc, gid=f'{gid}-arrow', cap='projecting')
    k0, k1 = y + g['ko_top'], y + g['ko_bottom']
    pb = y + g['promoter_baseline']
    above, below = spec.get('promoter_gap_mm', [1.76, 0.56])
    for xb, b in zip(xs, br):
        arrow_down(page, xb, ybar, k0, lw, lc, head, gid=f'{gid}-arrow', gap=(pb - above, pb + below))
        page.text(xb, pb, b['promoter'], size=psize, ha='center', gid=f'{gid}-text')
        node(xb - bw / 2, k0, bw, k1 - k0, b['text'], b['genotype'], b.get('tint', spec.get('tint', 0.3)), bsize)
    x0 = min(d[0] for d in drawn); y0 = min(d[1] for d in drawn)
    x1 = max(d[0] + d[2] for d in drawn); y1 = max(d[1] + d[3] for d in drawn)
    return dict(box_w=round(bw, 3), x_branch=[round(v, 3) for v in xs], boxes_mm=drawn,
                panel_box_mm=[round(x0, 3), round(y0, 3), round(x1 - x0, 3), round(y1 - y0, 3)])
