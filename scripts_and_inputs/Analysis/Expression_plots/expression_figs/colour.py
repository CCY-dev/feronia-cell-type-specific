# Vendored unchanged from the FER figure kit v3.2 (ferfigs/panels/fig1_colour.py), so the expression colours are
# identical to Fig. 1a (Seurat FeaturePlot defaults lightgrey -> blue, interpolated in CIE-LAB).
"""Colour scales used by Fig. 1, implemented to give the same colours as the R code that defines them.

ColorRamp2 reproduces circlize::colorRamp2(breaks, colors) with its default space = "LAB" (her notebook's
sessionInfo lists circlize 0.4.16; the verification re-run in provenance/fig1/r_session.txt used circlize 0.4.18 with
colorspace 2.1.3, and its output is what this module is checked against): values are clamped to the break range,
binned with .bincode(right = TRUE,
include.lowest = TRUE), interpolated linearly in CIE-LAB between the two neighbouring break colours, converted back
to sRGB, clipped to [0, 1] and written as hex with colorspace's rounding (int)(255 * v + 0.5). The sRGB <-> LAB
conversion uses the constants of the R package colorspace (D65 white 95.047 / 100 / 108.883, the colorspace RGB
matrices and gamma 2.4 with the 0.03928 / 0.00304 thresholds). tools/fig1_tables.py --check-r compares this
module with R output on a 211-value grid and on all 272 cells of Fig. 1d (identical hex codes required).

Gradient2 is the two-colour LAB gradient used for the FER expression maps (Seurat FeaturePlot defaults,
lightgrey to blue, minimum to maximum of the data).
"""
from __future__ import annotations

import math

import numpy as np

# R colour names used by the Fig. 1 code (grDevices::col2rgb values)
R_COLOURS = {
    'white': (255, 255, 255), 'beige': (245, 245, 220), 'red': (255, 0, 0), 'lightgrey': (211, 211, 211),
    'grey': (190, 190, 190), 'blue': (0, 0, 255), 'black': (0, 0, 0),
}

_WHITE = (95.047, 100.000, 108.883)
_M = np.array([[0.412453, 0.357580, 0.180423],
               [0.212671, 0.715160, 0.072169],
               [0.019334, 0.119193, 0.950227]])
_MI = np.array([[3.240479, -1.537150, -0.498535],
                [-0.969256, 1.875992, 0.041556],
                [0.055648, -0.204043, 1.057311]])
_EPS = 216.0 / 24389.0
_KAP = 24389.0 / 27.0
_GAMMA = 2.4


def col2rgb(c):
    """R colour name or '#rrggbb' -> (r, g, b) in 0-255 integers."""
    if isinstance(c, str):
        if c.startswith('#'):
            return tuple(int(c[i:i + 2], 16) for i in (1, 3, 5))
        if c not in R_COLOURS:
            raise KeyError(f"unknown R colour name '{c}' (add it to R_COLOURS with its grDevices::col2rgb value)")
        return R_COLOURS[c]
    return tuple(int(v) for v in c)


def _ftrans(u):          # sRGB -> linear (colorspace)
    return np.where(u > 0.03928, ((u + 0.055) / 1.055) ** _GAMMA, u / 12.92)


def _gtrans(u):          # linear -> sRGB (colorspace)
    u = np.asarray(u, float)
    return np.where(u > 0.00304, 1.055 * np.power(np.abs(u), 1 / _GAMMA) * np.sign(u) - 0.055, 12.92 * u)


def _f(t):
    return np.where(t > _EPS, np.cbrt(t), (_KAP * t + 16.0) / 116.0)


def srgb_to_lab(rgb01):
    rgb = _ftrans(np.asarray(rgb01, float))
    xyz = (rgb @ _M.T) * _WHITE[1]
    xr, yr, zr = xyz[..., 0] / _WHITE[0], xyz[..., 1] / _WHITE[1], xyz[..., 2] / _WHITE[2]
    L = np.where(yr > _EPS, 116.0 * np.cbrt(yr) - 16.0, _KAP * yr)
    return np.stack([L, 500.0 * (_f(xr) - _f(yr)), 200.0 * (_f(yr) - _f(zr))], -1)


def lab_to_srgb(lab):
    lab = np.asarray(lab, float)
    L, A, B = lab[..., 0], lab[..., 1], lab[..., 2]
    XN, YN, ZN = _WHITE
    Y = np.where(L <= 0, 0.0, np.where(L <= 8.0, L * YN / _KAP, np.where(L <= 100, YN * ((L + 16) / 116) ** 3, YN)))
    fy = np.where(Y <= _EPS * YN, (_KAP * Y / YN + 16) / 116, np.cbrt(Y / YN))
    fx = fy + A / 500.0
    fz = fy - B / 200.0
    X = np.where(fx ** 3 <= _EPS, XN * (fx * 116 - 16) / _KAP, XN * fx ** 3)
    Z = np.where(fz ** 3 <= _EPS, ZN * (fz * 116 - 16) / _KAP, ZN * fz ** 3)
    rgb = (np.stack([X, Y, Z], -1) @ _MI.T) / YN
    return _gtrans(rgb)


def to_hex(rgb01):
    """colorspace::hex rounding: (int)(255 * v + 0.5), values already clipped to [0, 1]."""
    v = np.floor(255 * np.clip(np.asarray(rgb01, float), 0, 1) + 0.5).astype(int)
    if v.ndim == 1:
        return '#%02X%02X%02X' % tuple(v)
    return ['#%02X%02X%02X' % tuple(r) for r in v]


def hex_to_rgb01(h):
    return tuple(int(h[i:i + 2], 16) / 255 for i in (1, 3, 5))


class ColorRamp2:
    """circlize::colorRamp2(breaks, colors) with space = 'LAB'. Call with values -> list of '#RRGGBB'."""

    def __init__(self, breaks, colours):
        od = np.argsort(np.asarray(breaks, float), kind='stable')
        self.breaks = np.asarray(breaks, float)[od]
        self.colours = [colours[i] for i in od]
        self.rgb01 = np.array([col2rgb(c) for c in self.colours], float) / 255
        self.lab = srgb_to_lab(self.rgb01)

    def __call__(self, x):
        x = np.atleast_1d(np.asarray(x, float))
        b = self.breaks
        x2 = np.clip(x, b[0], b[-1])
        # .bincode(x2, breaks, right = TRUE, include.lowest = TRUE): bin i holds (b[i], b[i+1]], first bin closed
        ib = np.searchsorted(b, x2, side='left') - 1
        ib = np.clip(ib, 0, len(b) - 2)
        out = []
        for v, i in zip(x2, ib):
            c1, c2 = self.lab[i], self.lab[i + 1]
            lab = (v - b[i + 1]) * (c2 - c1) / (b[i + 1] - b[i]) + c2
            out.append(to_hex(np.clip(lab_to_srgb(lab), 0, 1)))
        return out


class Gradient2:
    """Two-colour gradient in LAB from vmin (low colour) to vmax (high colour); values clamped."""

    def __init__(self, low, high, vmin, vmax):
        self.vmin, self.vmax = float(vmin), float(vmax)
        self.lab = srgb_to_lab(np.array([col2rgb(low), col2rgb(high)], float) / 255)

    def rgb01(self, x):
        t = np.clip((np.asarray(x, float) - self.vmin) / (self.vmax - self.vmin), 0, 1)[..., None]
        lab = self.lab[0] * (1 - t) + self.lab[1] * t
        return np.floor(255 * np.clip(lab_to_srgb(lab), 0, 1) + 0.5) / 255


def luminance_contrast_black(rgb01):
    """WCAG contrast ratio of black text on this fill (>= 4.5 keeps black text legible)."""
    def lin(c):
        return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (lin(c) for c in rgb01)
    L = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return (L + 0.05) / 0.05


def _selftest():
    r = ColorRamp2([0, -math.log10(0.05), 10], ['white', 'beige', 'red'])
    return r([0, -math.log10(0.05), 10])


if __name__ == '__main__':
    print(_selftest())
