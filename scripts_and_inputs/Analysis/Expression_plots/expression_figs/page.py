"""A page in millimetres (origin top-left, y downwards) on one full-page matplotlib axes.

All text is vector (TrueType, pdf.fonttype 42, so it stays editable); markup '*...*' sets italic runs and '**...**'
bold runs. Widths are measured with matplotlib's TextToPath (unhinted font metrics), so centred and right-aligned
mixed-style strings land on their anchors. Images (UMAP maps, colour bars) are placed at their native pixel size
with interpolation 'none', so the PDF holds each array unresampled.
"""
from __future__ import annotations

import re

import matplotlib
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextToPath

MM_PER_PT = 25.4 / 72.0
MINUS = '\u2212'
_T2P = TextToPath()


def setup_rc(family):
    matplotlib.rcParams.update({
        'pdf.fonttype': 42, 'svg.fonttype': 'none', 'font.family': 'sans-serif', 'font.sans-serif': list(family),
        'path.simplify': False, 'image.interpolation': 'none', 'savefig.dpi': 600,
        'pdf.compression': 6,
    })
    matplotlib.rcParams['svg.hashsalt'] = 'fer-si-expression'      # deterministic SVG ids


def runs(s):
    """'a *b* **c**' -> [(text, style, weight)]."""
    out = []
    for part in re.split(r'(\*\*[^*]+\*\*|\*[^*]+\*)', s):
        if not part:
            continue
        if part.startswith('**'):
            out.append((part[2:-2], 'normal', 'bold'))
        elif part.startswith('*'):
            out.append((part[1:-1], 'italic', 'normal'))
        else:
            out.append((part, 'normal', 'normal'))
    return out


class Page:
    def __init__(self, w_mm, h_mm, family=('Arial', 'Helvetica', 'DejaVu Sans'), size_pt=6.0):
        setup_rc(family)
        self.w, self.h, self.family, self.size = float(w_mm), float(h_mm), list(family), size_pt
        self.fig = Figure(figsize=(w_mm / 25.4, h_mm / 25.4))
        self.ax = self.fig.add_axes((0, 0, 1, 1))
        self.ax.set_xlim(0, w_mm)
        self.ax.set_ylim(h_mm, 0)
        self.ax.set_axis_off()
        self.ax.set_aspect('auto')
        self.items = []          # (kind, gid, bbox mm) for checks
        self.bottom = 0.0        # lowest drawn y (mm), for fit_height

    def _b(self, *ys):
        self.bottom = max(self.bottom, *[float(v) for v in ys])

    def fit_height(self, margin=0.5):
        """Shrink/grow the page to the drawn content; the axes spans the figure, so positions keep their mm."""
        h = round(self.bottom + margin, 2)
        self.h = h
        self.fig.set_size_inches(self.w / 25.4, h / 25.4)
        self.ax.set_ylim(h, 0)
        return h

    # ---- text ---------------------------------------------------------------------------------------------------
    def _prop(self, size, style, weight):
        return FontProperties(family=self.family, size=size, style=style, weight=weight)

    def text_width(self, s, size=None, weight='normal'):
        size = size or self.size
        w = 0.0
        for t, st, wt in runs(s):
            wt = 'bold' if weight == 'bold' else wt
            w += _T2P.get_text_width_height_descent(t, self._prop(size, st, wt), ismath=False)[0]
        return w * MM_PER_PT

    def text(self, x, y, s, size=None, ha='left', va='baseline', rotation=0.0, colour='k', weight='normal', gid=None):
        """Markup text at (x, y) mm. rotation 0 or 90 (reads bottom to top); va: baseline | center | top."""
        size = size or self.size
        total = self.text_width(s, size, weight)
        cap = 0.716 * size * MM_PER_PT                    # Arial cap height
        off = {'left': 0.0, 'center': -total / 2, 'right': -total}[ha]
        voff = {'baseline': 0.0, 'center': cap / 2, 'top': cap}[va]
        pos = off
        for t, st, wt in runs(s):
            wt = 'bold' if weight == 'bold' else wt
            if rotation == 0:
                xx, yy = x + pos, y + voff
            else:                                          # 90 deg: advance runs upward; voff moves right
                xx, yy = x + voff, y - pos
            self.ax.text(xx, yy, t, fontproperties=self._prop(size, st, wt), ha='left', va='baseline',
                         rotation=rotation, rotation_mode='anchor', color=colour, gid=gid)
            pos += self.text_width(t, size, wt)
        if rotation == 0:
            box = (x + off, y + voff - cap, x + off + total, y + voff + 0.21 * size * MM_PER_PT)
        else:
            box = (x + voff - cap, y - off - total, x + voff + 0.21 * size * MM_PER_PT, y - off)
        self.items.append(('text', gid, box, s))
        self._b(box[3])
        return box

    # ---- shapes ------------------------------------------------------------------------------------------------
    def line(self, x0, y0, x1, y1, lw=0.5, colour='k', gid=None, zorder=3):
        self._b(y0, y1)
        self.ax.plot([x0, x1], [y0, y1], lw=lw, color=colour, solid_capstyle='butt', gid=gid, zorder=zorder)

    def polyline(self, pts, lw=0.5, colour='k', gid=None, zorder=3):
        xs, ys = zip(*pts)
        self._b(*ys)
        self.ax.plot(xs, ys, lw=lw, color=colour, solid_capstyle='butt', solid_joinstyle='miter', gid=gid,
                     zorder=zorder)

    def polygon(self, xs, ys, fill, edge=None, lw=0.0, gid=None, zorder=2):
        self._b(max(ys))
        self.ax.fill(xs, ys, facecolor=fill, edgecolor=edge if edge is not None else 'none', linewidth=lw,
                     joinstyle='round', gid=gid, zorder=zorder)

    def rect(self, x, y, w, h, fill='none', edge='none', lw=0.0, gid=None, zorder=1):
        from matplotlib.patches import Rectangle
        self._b(y + h)
        self.ax.add_patch(Rectangle((x, y), w, h, facecolor=fill, edgecolor=edge, linewidth=lw, gid=gid,
                                    zorder=zorder))

    def dots(self, xs, ys, d_mm, colours, edge=(0.45, 0.45, 0.45), lw=0.2, gid=None, zorder=4):
        """Filled circles with diameter d_mm (array) at (xs, ys) mm."""
        import numpy as np
        s = (np.asarray(d_mm, float) / MM_PER_PT) ** 2
        self._b(max(np.asarray(ys, float) + np.asarray(d_mm, float) / 2))
        self.ax.scatter(xs, ys, s=s, c=colours, edgecolors=[edge], linewidths=lw, marker='o', gid=gid,
                        zorder=zorder)

    def image(self, arr, x, y, w, h, gid=None, zorder=1):
        self.ax.imshow(arr, extent=(x, x + w, y + h, y), origin='upper', interpolation='none', aspect='auto',
                       gid=gid, zorder=zorder)
        self.items.append(('image', gid, (x, y, x + w, y + h), arr.shape))
        self._b(y + h)

    def save(self, stem, formats=('pdf', 'png'), png_dpi=600):
        paths = []
        for f in formats:
            p = f'{stem}.{f}'
            kw = dict(dpi=png_dpi) if f == 'png' else {}
            if f == 'pdf':
                kw['metadata'] = {'CreationDate': None, 'ModDate': None, 'Producer': None, 'Creator': None}
            if f == 'svg':
                kw['metadata'] = {'Date': None, 'Creator': None}
            self.fig.savefig(p, format=f, **kw)
            paths.append(p)
        return paths
