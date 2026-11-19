"""Shared drawing helpers for the Fig. 6 panels (new primitives; the core library is unchanged).

Everything is in millimetres from the top-left of the page, like ferfigs.layout.Page. Patches that need true
geometry (circles, rounded boxes) go on one full-page axes whose data units are millimetres (aspect 1:1), so a
circle stays a circle. All text still goes through page.text (Arial); only DNA sequences use mono_text.
"""
from __future__ import annotations

import math

import numpy as np
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyBboxPatch, Polygon, Rectangle
from matplotlib.text import Text

from .. import style as S
from ..schemas import ValidationError

PT = 25.4 / 72                 # mm per point
CAP = 0.716                    # Arial cap/ascender height, fraction of the font size


def display(name, local=None, kind='genotypes'):
    """Display markup: config/labels.yaml first, then the figure config's display_labels (Fig. 6 names that are
    not yet in labels.yaml), else the name itself."""
    lab = (S.labels().get(kind) or {}).get(name)
    if lab is not None:
        return lab
    return ((local or {}).get(kind) or {}).get(name, name)


def hexrgb(c):
    """'#rrggbb' or [r, g, b] (0-255 or 0-1) -> (r, g, b) in 0-1."""
    return S.rgb(c)


def cap_mm(size_pt):
    return CAP * size_pt * PT


class MM:
    """One full-page axes in millimetres (y down) for geometric patches. zorder 2: above image axes (0),
    below figure-level lines (4) and text (5) drawn by Page."""

    def __init__(self, page, gid='fig6-mm', zorder=2):
        self.page = page
        ax = page.axes(0, 0, page.W, page.H, gid=gid, frame=False)
        ax.set_xlim(0, page.W)
        ax.set_ylim(page.H, 0)
        ax.set_zorder(zorder)
        ax.patch.set_visible(False)
        self.ax = ax

    def line(self, x0, y0, x1, y1, lw, colour='k', alpha=None, dashes=None, gid=None, z=2, cap='butt'):
        ln = Line2D([x0, x1], [y0, y1], lw=lw, color=_c(colour), alpha=alpha, solid_capstyle=cap, zorder=z)
        if dashes:
            ln.set_dashes(dashes)
        if gid:
            ln.set_gid(gid)
        self.ax.add_line(ln)
        return ln

    def circle(self, x, y, d, face, edge=None, lw=0.0, gid=None, z=3):
        c = Circle((x, y), d / 2, facecolor=_c(face), edgecolor=_c(edge) if edge is not None else 'none',
                   lw=lw, zorder=z)
        if gid:
            c.set_gid(gid)
        self.ax.add_patch(c)
        return c

    def poly(self, pts, face, edge=None, lw=0.0, gid=None, z=3, join='miter', alpha=None):
        p = Polygon(pts, closed=True, facecolor=_c(face), edgecolor=_c(edge) if edge is not None else 'none',
                    lw=lw, joinstyle=join, zorder=z, alpha=alpha)
        if gid:
            p.set_gid(gid)
        self.ax.add_patch(p)
        return p

    def rect(self, x, y, w, h, face, edge=None, lw=0.0, gid=None, z=1, alpha=None):
        r = Rectangle((x, y), w, h, facecolor=_c(face), edgecolor=_c(edge) if edge is not None else 'none',
                      lw=lw, zorder=z, joinstyle='miter', alpha=alpha)
        if gid:
            r.set_gid(gid)
        self.ax.add_patch(r)
        return r

    def rbox(self, x, y, w, h, radius, face, edge=None, lw=0.0, gid=None, z=1, alpha=None):
        b = FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0,rounding_size={radius}', facecolor=_c(face),
                           edgecolor=_c(edge) if edge is not None else 'none', lw=lw, zorder=z, alpha=alpha)
        if gid:
            b.set_gid(gid)
        self.ax.add_patch(b)
        return b


def _c(c):
    if c is None or (isinstance(c, str) and c == 'none'):
        return 'none'
    if isinstance(c, str) and not c.startswith('#'):
        return c                                  # matplotlib colour name, e.g. 'k', 'white'
    return hexrgb(c)


# ------------------------------------------------------------------ monospace text (DNA sequences only)
def mono_family(candidates):
    """First installed family from the config list; the journal profile allows Courier, Courier New,
    Andale Mono and DejaVu Sans Mono for sequences."""
    installed = {f.name for f in fm.fontManager.ttflist}
    for fam in candidates:
        if fam in installed:
            return fam
    raise ValidationError([f'none of the sequence fonts {candidates} is installed (config fig6.yaml: '
                           f'fonts.sequence). DejaVu Sans Mono ships with matplotlib.'])


def mono_advance(page, family, size):
    """Advance width of one character (mm), measured with the page renderer on 40 characters."""
    t = Text(0, 0, 'ACGT' * 10, fontsize=size, fontfamily=family)
    t.set_figure(page.fig)
    w, _, _ = page._r.get_text_width_height_descent('ACGT' * 10, t.get_fontproperties(), ismath=False)
    return w * 72.0 / page.fig.dpi / 40 * PT


def mono_text(page, x, y, s, family, size, gid=None, colour=(0, 0, 0), zorder=5):
    """Left-aligned, baseline-anchored monospace text at (x, y) mm."""
    a = page.fig.text(x / page.W, 1 - y / page.H, s, fontsize=size, fontfamily=family, color=hexrgb(colour),
                      ha='left', va='baseline', zorder=zorder)
    if gid:
        a.set_gid(gid)
    return a


# ------------------------------------------------------------------ colour helpers
_M = np.array([[0.4124564, 0.3575761, 0.1804375], [0.2126729, 0.7151522, 0.0721750],
               [0.0193339, 0.1191920, 0.9503041]])
_WHITE = np.array([0.95047, 1.0, 1.08883])


def _lab(rgb255):
    c = np.asarray(rgb255, float) / 255
    lin = np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
    xyz = _M @ lin / _WHITE
    f = np.where(xyz > (6 / 29) ** 3, np.cbrt(xyz), xyz / (3 * (6 / 29) ** 2) + 4 / 29)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])


def _rgb(lab):
    fy = (lab[0] + 16) / 116
    fx, fz = fy + lab[1] / 500, fy - lab[2] / 200
    inv = [t ** 3 if t > 6 / 29 else 3 * (6 / 29) ** 2 * (t - 4 / 29) for t in (fx, fy, fz)]
    lin = np.clip(np.linalg.solve(_M, np.array(inv) * _WHITE), 0, 1)
    c = np.where(lin <= 0.0031308, 12.92 * lin, 1.055 * lin ** (1 / 2.4) - 0.055)
    return tuple(int(v) for v in np.clip(np.round(c * 255), 0, 255))


def lab_ramp(lo_hex, hi_hex):
    """f(t in 0-1) -> (r, g, b) 0-255, linear in CIELAB (ggplot2 scale_fill_gradient)."""
    lo, hi = _lab([int(lo_hex[i:i + 2], 16) for i in (1, 3, 5)]), _lab([int(hi_hex[i:i + 2], 16) for i in (1, 3, 5)])

    def f(t):
        t = min(max(float(t), 0.0), 1.0)
        return _rgb(lo + t * (hi - lo))
    return f


def luminance(rgb255):
    c = np.asarray(rgb255, float) / 255
    lin = np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
    return float(0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2])


def contrast_text(rgb255):
    """(white or black, contrast ratio): whichever contrasts more with the fill (WCAG 2)."""
    lc = luminance(rgb255)
    cw, cb = 1.05 / (lc + 0.05), (lc + 0.05) / 0.05
    return ((1, 1, 1), cw) if cw >= cb else ((0, 0, 0), cb)


def tint(hexc, t, base=(255, 255, 255)):
    c = [int(hexc[i:i + 2], 16) for i in (1, 3, 5)]
    return '#%02x%02x%02x' % tuple(round(b + (x - b) * t) for x, b in zip(c, base))


def gradient_image(func, w_mm, h_mm, ppi=600, top_value=1.0, bottom_value=0.0):
    """A vertical colour ramp as an RGB array at `ppi` for its printed size (continuous, no visible steps)."""
    nw, nh = max(2, int(round(w_mm / 25.4 * ppi))), max(2, int(round(h_mm / 25.4 * ppi)))
    vals = np.linspace(top_value, bottom_value, nh)
    col = np.array([func(v) for v in vals], dtype=np.uint8)
    return np.repeat(col[:, None, :], nw, axis=1)


# ------------------------------------------------------------------ 1-D label spreading
def spread(items, min_gap=0.25, pinned=None):
    """items: [(anchor_x, width)] in mm, sorted by anchor. Returns centre x per item so that neighbouring labels
    are at least min_gap apart. Overlapping free (not pinned) neighbours move apart as one block centred on their
    anchors; pinned items (e.g. -/+ signs over gel lanes) never move. Raises ValidationError when a block still
    hits a pinned neighbour (make the panel wider)."""
    pinned = set(pinned or [])
    clusters = [[i] for i in range(len(items))]

    def layout(cl):
        tot = sum(items[i][1] for i in cl) + min_gap * (len(cl) - 1)
        mid = sum(items[i][0] for i in cl) / len(cl)
        x, out = mid - tot / 2, {}
        for i in cl:
            out[i] = x + items[i][1] / 2
            x += items[i][1] + min_gap
        return out
    while True:
        pos = {}
        for cl in clusters:
            pos.update(layout(cl))
        merged = False
        for k in range(len(clusters) - 1):
            a, b = clusters[k][-1], clusters[k + 1][0]
            if a in pinned or b in pinned:
                continue
            if pos[a] + items[a][1] / 2 + min_gap > pos[b] - items[b][1] / 2 + 1e-9:
                clusters[k] = clusters[k] + clusters.pop(k + 1)
                merged = True
                break
        if not merged:
            break
    out = [pos[i] for i in range(len(items))]
    for i in range(len(items) - 1):
        if out[i] + items[i][1] / 2 > out[i + 1] - items[i + 1][1] / 2 + 1e-9:
            raise ValidationError([f'labels at x = {items[i][0]:.2f} and {items[i + 1][0]:.2f} mm overlap; widen the '
                                   f'panel (config) or shorten the labels.'])
    return out


def arrow_head(mm, x, y, ang, length, width, colour, gid=None, z=3):
    """Filled triangular head with its tip at (x, y), pointing along angle `ang` (radians, page coordinates)."""
    ux, uy = math.cos(ang), math.sin(ang)
    bx, by = x - ux * length, y - uy * length
    px, py = -uy * width / 2, ux * width / 2
    return mm.poly([(x, y), (bx + px, by + py), (bx - px, by - py)], face=colour, gid=gid, z=z)
