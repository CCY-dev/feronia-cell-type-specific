"""A figure page in millimetres.

All positions are millimetres from the top-left corner of the page. Every element is a matplotlib artist, so the
PDF keeps editable text (TrueType) and vector line art, and the SVG keeps <text> elements with named groups
(gid) for Illustrator. Images are embedded at their native pixels (no resampling) unless they exceed
style.output.image_max_ppi at their printed size, in which case they are downsampled once (Lanczos) and the
placement record says so.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Polygon
from matplotlib.text import Text

from . import style as S

MM = S.MM
PT_PER_MM = S.PT_PER_MM
# Arial metrics (fractions of the font size): ascender of b/d/h, cap height, x-height
ASCENDER = 0.716
XHEIGHT = 0.519


@dataclass
class Placed:
    """An image placed on the page: pixel (column, row) -> page mm via x0 + col * mm_per_px."""
    name: str
    x0: float
    y0: float
    mm_per_px: float
    w_px: int
    h_px: int
    src_px: tuple
    ppi: float
    resampled: bool = False
    sha256: str = ''
    crop: tuple = None
    orient: str = 'identity'

    @property
    def w_mm(self):
        return self.w_px * self.mm_per_px

    @property
    def h_mm(self):
        return self.h_px * self.mm_per_px

    def mm(self, col, row):
        return self.x0 + col * self.mm_per_px, self.y0 + row * self.mm_per_px


DIHEDRAL = {'identity': 0, 'rot90ccw': 1, 'rot180': 2, 'rot270ccw': 3, 'mirrorLR': 4, 'mirrorLR+rot90ccw': 5,
            'mirrorTB': 6, 'mirrorLR+rot270ccw': 7}


def oriented(arr, orient):
    k = DIHEDRAL[orient] if isinstance(orient, str) else int(orient or 0)
    a = np.fliplr(arr) if k >= 4 else arr
    return np.rot90(a, k % 4)


class Page:
    def __init__(self, width_mm: float, height_mm: float, name: str):
        self.fam = S.apply_rc()
        self.st = S.load_style()
        self.W, self.H, self.name = float(width_mm), float(height_mm), name
        self.fig = plt.figure(figsize=(self.W * MM, self.H * MM), dpi=300)
        self._r = self.fig.canvas.get_renderer()
        self.images: list[Placed] = []
        self.letters: dict[str, tuple] = {}
        self.records: dict = {}            # per-panel records for checks and source data

    # ------------------------------------------------------------ coordinates
    def fx(self, x):
        return x / self.W

    def fy(self, y):
        return 1 - y / self.H

    def axes(self, x, y, w, h, gid=None, frame=True):
        ax = self.fig.add_axes((self.fx(x), self.fy(y + h), w / self.W, h / self.H))
        if gid:
            ax.set_gid(gid)
        if not frame:
            ax.set_axis_off()
        return ax

    # ------------------------------------------------------------ text
    def _width_pt(self, s, size, fontstyle='normal', weight='normal', family=None):
        if not s:
            return 0.0
        t = Text(0, 0, s, fontsize=size, fontstyle=fontstyle, fontweight=weight, fontfamily=family or self.fam)
        t.set_figure(self.fig)
        w, _, _ = self._r.get_text_width_height_descent(s, t.get_fontproperties(), ismath=False)
        return w * 72.0 / self.fig.dpi

    def text_width(self, s, size=None):
        """Width in mm of markup text."""
        size = size or self.st['font']['size_pt']
        return sum(self._width_pt(t, S.script_size(size, sc), fs, fw) for t, fs, fw, sc in S.parse_runs(s)) / PT_PER_MM

    def text(self, x, y, s, size=None, ha='left', va='baseline', rotation=0.0, colour=None, gid=None,
             weight=None, fontstyle=None, zorder=5, family=None):
        """Markup text ('*italic*', '**bold**') with its anchor at (x, y) mm. ha: left|center|right along the
        baseline; va: baseline|center|top|bottom (center = middle of the x-height... cap height for digits).
        Returns the text width in mm."""
        size = size or self.st['font']['size_pt']
        fam = S.resolve_mono() if family == 'mono' else self.fam
        parts = S.parse_runs(s) if family != 'mono' else [(s, 'normal', 'normal', None)]
        if weight or fontstyle:
            parts = [(t, fontstyle or fs, weight or fw, sc) for t, fs, fw, sc in parts]
        widths = [self._width_pt(t, S.script_size(size, sc), fs, fw, fam) for t, fs, fw, sc in parts]
        total = sum(widths)
        th = math.radians(rotation)
        ux, uy = math.cos(th), math.sin(th)            # baseline direction (points, y up)
        nx, ny = -uy, ux                               # up direction
        start = {'left': 0.0, 'center': -total / 2, 'right': -total}[ha]
        lift = {'baseline': 0.0, 'center': -0.36 * size, 'top': -ASCENDER * size, 'bottom': 0.21 * size}[va]
        col = S.rgb(colour) if colour is not None else S.rgb(self.st['colours']['text'])
        x0_pt, y0_pt = x * PT_PER_MM, (self.H - y) * PT_PER_MM
        off = start
        for (t, fs, fw, sc), w in zip(parts, widths):
            up = lift + S.script_shift(size, sc)
            px = x0_pt + off * ux + up * nx
            py = y0_pt + off * uy + up * ny
            a = self.fig.text(px / (self.W * PT_PER_MM), py / (self.H * PT_PER_MM), t, fontsize=S.script_size(size, sc),
                              fontstyle=fs, fontweight=fw, fontfamily=fam, color=col, rotation=rotation,
                              rotation_mode='anchor', ha='left', va='baseline', zorder=zorder)
            if gid:
                a.set_gid(gid)
            off += w
        return total / PT_PER_MM

    def panel_letter(self, letter, x, y_top):
        """Nature style: 8 pt bold upright lower case. y_top = top of the tallest lower-case letter (b, d, f, h)."""
        size = self.st['font']['panel_label_pt']
        base = y_top + ASCENDER * size / PT_PER_MM
        self.text(x, base, letter, size=size, weight='bold', gid=f'{self.name}-letter-{letter}')
        self.letters[letter] = (round(x, 3), round(y_top, 3))

    # ------------------------------------------------------------ line art
    def line(self, x0, y0, x1, y1, lw=None, colour='k', gid=None, zorder=4, cap='butt'):
        lw = lw if lw is not None else self.st['lines']['bracket_pt']
        ln = Line2D([self.fx(x0), self.fx(x1)], [self.fy(y0), self.fy(y1)], transform=self.fig.transFigure,
                    lw=lw, color=S.rgb(colour) if not isinstance(colour, str) else colour, solid_capstyle=cap,
                    zorder=zorder)
        if gid:
            ln.set_gid(gid)
        self.fig.add_artist(ln)
        return ln

    def polyline(self, pts, lw=None, colour='k', gid=None, closed=False, fill=None, zorder=4):
        lw = lw if lw is not None else self.st['lines']['bracket_pt']
        xy = [(self.fx(x), self.fy(y)) for x, y in pts]
        p = Polygon(xy, closed=closed, transform=self.fig.transFigure, lw=lw,
                    edgecolor=(S.rgb(colour) if colour is not None and not isinstance(colour, str) else colour),
                    facecolor=(S.rgb(fill) if fill is not None and not isinstance(fill, str) else (fill or 'none')),
                    joinstyle='miter', zorder=zorder)
        if gid:
            p.set_gid(gid)
        self.fig.add_artist(p)
        return p

    def rect(self, x, y, w, h, lw=None, edge='k', face='none', gid=None, zorder=4):
        lw = lw if lw is not None else self.st['lines']['bracket_pt']
        r = Rectangle((self.fx(x), self.fy(y + h)), w / self.W, h / self.H, transform=self.fig.transFigure, lw=lw,
                      edgecolor=(S.rgb(edge) if not isinstance(edge, str) else edge),
                      facecolor=(S.rgb(face) if not isinstance(face, str) else face), zorder=zorder)
        if gid:
            r.set_gid(gid)
        self.fig.add_artist(r)
        return r

    def bracket(self, x0, x1, y, tick=0.7, lw=None, gid=None, down=True):
        """Horizontal line at y from x0 to x1 with end ticks (downward by default)."""
        d = tick if down else -tick
        return self.polyline([(x0, y + d), (x0, y), (x1, y), (x1, y + d)], lw=lw, gid=gid)

    def header(self, x0, x1, y_base, s, size=None, bracket=True, gap=0.9, gid=None):
        """Centred label with an optional bracket `gap` mm below its baseline spanning x0..x1."""
        self.text((x0 + x1) / 2, y_base, s, size=size, ha='center', gid=gid)
        if bracket:
            self.bracket(x0 + 0.15, x1 - 0.15, y_base + gap, gid=gid)

    # ------------------------------------------------------------ images
    def image(self, src, x, y, w=None, h=None, crop=None, orient='identity', gid=None, name=None, zorder=1):
        """Place an image file (or array) with its top-left at (x, y) mm. Give w or h (mm); aspect is kept.
        crop = (x0, y0, x1, y1) in pixels of the oriented image."""
        if isinstance(src, (str, Path)):
            raw = Path(src).read_bytes()
            with Image.open(io.BytesIO(raw)) as im0:
                im0.load()
                im = im0.copy()
            name = name or Path(src).name
        else:
            im = Image.fromarray(np.asarray(src))
            raw = np.asarray(src).tobytes()
        if im.mode not in ('RGB', 'L'):
            im = im.convert('RGBA')
            bg = Image.new('RGBA', im.size, (255, 255, 255, 255))
            im = Image.alpha_composite(bg, im).convert('RGB')
        arr = oriented(np.asarray(im), orient)
        if crop is not None:
            x0c, y0c, x1c, y1c = [int(round(v)) for v in crop]
            arr = arr[y0c:y1c, x0c:x1c]
        hp, wp = arr.shape[:2]
        if w is None and h is None:
            raise ValueError('give w or h')
        if w is None:
            w = h * wp / hp
        mm_per_px = w / wp
        ppi = wp / (w / 25.4)
        resampled = False
        maxppi = self.st['output'].get('image_max_ppi')
        if maxppi and ppi > maxppi * 1.02:
            nw = int(round(w / 25.4 * maxppi)); nh = int(round(hp * nw / wp))
            arr = np.asarray(Image.fromarray(np.ascontiguousarray(arr)).resize((nw, nh), Image.LANCZOS))
            resampled = True
        hh = w * hp / wp
        ax = self.axes(x, y, w, hh, gid=gid, frame=False)
        cmap = 'gray' if arr.ndim == 2 else None
        ax.imshow(arr, extent=(0, wp, hp, 0), interpolation='none', aspect='auto', cmap=cmap,
                  vmin=0 if cmap else None, vmax=255 if cmap else None, zorder=zorder)
        ax.set_xlim(0, wp); ax.set_ylim(hp, 0)
        pl = Placed(name or 'array', x, y, mm_per_px, wp, hp, tuple(arr.shape[:2]), round(ppi, 1), resampled,
                    hashlib.sha256(raw).hexdigest()[:16], tuple(crop) if crop is not None else None,
                    orient if isinstance(orient, str) else str(orient))
        self.images.append(pl)
        return pl

    def scalebar(self, placed: Placed, length_um, um_per_px, corner='br', pad=0.8, lw=1.0, colour='white',
                 label=None, gid=None):
        """Vector scale bar on top of an image (never burned into pixels). Bar length = length_um / um_per_px
        pixels, drawn in the image's own mm scale. Optional label above the bar."""
        L = length_um / um_per_px * placed.mm_per_px
        x1 = placed.x0 + placed.w_mm - pad if corner.endswith('r') else placed.x0 + pad + L
        y = placed.y0 + placed.h_mm - pad if corner.startswith('b') else placed.y0 + pad
        self.line(x1 - L, y, x1, y, lw=lw, colour=colour, gid=gid or f'{self.name}-scalebar', zorder=6)
        if label:
            self.text(x1 - L / 2, y - 0.6, label, ha='center', colour=colour, size=self.st['font']['size_pt'])
        return L

    # ------------------------------------------------------------ output
    def save(self, out_dir, stem=None, formats=None, check=True):
        out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
        stem = stem or self.name
        formats = formats or self.st['output']['formats']
        paths = {}
        for fmt in formats:
            p = out / f'{stem}.{fmt}'
            if fmt == 'pdf':
                self.fig.savefig(p, format='pdf', metadata={'Creator': 'ferfigs', 'CreationDate': None,
                                                            'ModDate': None})
            elif fmt == 'svg':
                self.fig.savefig(p, format='svg', metadata={'Date': None, 'Creator': 'ferfigs'})
            elif fmt == 'png':
                self.fig.savefig(p, format='png', dpi=self.st['output']['png_dpi'], metadata={'Software': 'ferfigs'})
            paths[fmt] = p
        rec = dict(name=self.name, size_mm=[self.W, self.H], font=self.fam, letters=self.letters,
                   images=[dict(name=i.name, x_mm=round(i.x0, 2), y_mm=round(i.y0, 2), w_mm=round(i.w_mm, 2),
                                h_mm=round(i.h_mm, 2), px=[i.w_px, i.h_px], ppi=i.ppi, resampled=i.resampled,
                                sha256_16=i.sha256, crop_px=list(i.crop) if i.crop else None, orient=i.orient)
                           for i in self.images],
                   panels=self.records)
        (out / f'{stem}.build.json').write_text(json.dumps(rec, indent=1, default=str))
        plt.close(self.fig)
        return paths, rec
