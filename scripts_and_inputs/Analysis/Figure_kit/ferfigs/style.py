"""House style and journal profile.

Everything visual that is not data comes from two YAML files:
  config/style.yaml                       fonts, sizes, line widths, palettes (house style)
  ferfigs/profiles/<profile>.yaml          journal limits used by `check` (Nature / Nature Plants by default)
"""
from __future__ import annotations

import re
import warnings
from functools import lru_cache
from pathlib import Path

import matplotlib
import yaml

matplotlib.use('Agg')
import matplotlib.font_manager as fm  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

KIT = Path(__file__).resolve().parents[1]
MM = 1 / 25.4            # inches per mm
PT_PER_MM = 72 / 25.4


@lru_cache(maxsize=None)
def load_style(path: str | None = None) -> dict:
    p = Path(path) if path else KIT / 'config' / 'style.yaml'
    with open(p) as f:
        return yaml.safe_load(f)


@lru_cache(maxsize=None)
def load_profile(name: str | None = None) -> dict:
    name = name or load_style()['profile']
    with open(Path(__file__).resolve().parent / 'profiles' / f'{name}.yaml') as f:
        return yaml.safe_load(f)


def rgb(c):
    """[0-255]*3 or [0-1]*3 or '#rrggbb' -> (r, g, b) floats in 0-1."""
    if isinstance(c, str):
        c = c.lstrip('#')
        return tuple(int(c[i:i + 2], 16) / 255 for i in (0, 2, 4))
    c = list(c)
    if any(v > 1 for v in c):
        return tuple(v / 255 for v in c[:3])
    return tuple(float(v) for v in c[:3])


def series_colour(name: str):
    cols = load_style()['colours']['series']
    if name not in cols:
        raise KeyError(f"No colour for series '{name}'. Add it under colours: series: in config/style.yaml "
                       f"(known: {', '.join(cols)}).")
    return rgb(cols[name])


def genotype_colours(name: str):
    """(fill, stroke) for a genotype; '_L1'-style line suffixes fall back to the base genotype."""
    st = load_style()['colours']
    base = re.sub(r'_L\d+$', '', name)
    if base not in st['genotype_fill']:
        raise KeyError(f"No colour for genotype '{name}'. Add it under colours: genotype_fill/genotype_stroke in "
                       f"config/style.yaml.")
    return rgb(st['genotype_fill'][base]), rgb(st['genotype_stroke'][base])


@lru_cache(maxsize=None)
def resolve_font() -> str:
    """First installed family from style.font.family. Warns when it is not Arial/Helvetica (journal requirement)."""
    fams = load_style()['font']['family']
    installed = {f.name for f in fm.fontManager.ttflist}
    for fam in fams:
        if fam in installed:
            if fam not in load_profile()['fonts']['allowed']:
                warnings.warn(f"Arial/Helvetica not installed; using '{fam}'. Text widths will differ slightly and "
                              f"`check` will flag the font. Install Arial before making submission files.")
            return fam
    warnings.warn('None of the configured fonts is installed; matplotlib default font used.')
    return matplotlib.rcParams['font.family'][0]


@lru_cache(maxsize=None)
def resolve_mono() -> str:
    """First installed monospace family from style.font.mono (sequences only; the profile allows these)."""
    fams = load_style()['font'].get('mono', ['Courier New', 'Courier', 'Andale Mono', 'DejaVu Sans Mono'])
    installed = {f.name for f in fm.fontManager.ttflist}
    for fam in fams:
        if fam in installed:
            return fam
    return 'monospace'


def apply_rc():
    """matplotlib settings for editable, deterministic, journal-style vector output."""
    st = load_style()
    fam = resolve_font()
    ln = st['lines']
    matplotlib.rcParams.update({
        'font.family': [fam],
        'font.size': st['font']['size_pt'],
        'axes.labelsize': st['font']['size_pt'],
        'xtick.labelsize': st['font']['size_pt'],
        'ytick.labelsize': st['font']['size_pt'],
        'legend.fontsize': st['font']['size_pt'],
        'axes.titlesize': st['font']['size_pt'],
        'axes.linewidth': ln['axis_pt'],
        'xtick.major.width': ln['tick_pt'],
        'ytick.major.width': ln['tick_pt'],
        'xtick.major.size': ln['tick_len_pt'],
        'ytick.major.size': ln['tick_len_pt'],
        'xtick.major.pad': 1.5,
        'ytick.major.pad': 1.5,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.unicode_minus': True,        # true minus signs in tick labels
        'axes.labelpad': 2.0,
        'axes.facecolor': 'none',
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
        'lines.solid_capstyle': 'butt',
        'pdf.fonttype': 42,                # TrueType: text stays editable, never outlined
        'ps.fonttype': 42,
        'svg.fonttype': 'none',            # SVG keeps <text>
        'svg.hashsalt': 'ferfigs',         # stable SVG ids (deterministic output)
        'pdf.compression': 6,
        'mathtext.fontset': 'custom',
        'mathtext.rm': fam,
        'mathtext.it': f'{fam}:italic',
        'mathtext.bf': f'{fam}:bold',
        'mathtext.default': 'regular',
        'image.interpolation': 'none',
        'image.resample': False,
        'legend.frameon': False,
        'path.simplify': False,            # keep every data point exactly where it is
    })
    plt.close('all')
    return fam


# ------------------------------------------------------------------ inline markup
_MARK = re.compile(r'(\*\*[^*]+\*\*|\*[^*]+\*|_\{[^{}]+\}|\^\{[^{}]+\})')
SCRIPT_SCALE = 0.85            # sub/superscript size relative to the text; never below 5 pt (journal minimum)


def parse_runs(s: str):
    """Markup -> runs (text, fontstyle, fontweight, script). '*italic*', '**bold**', 'log_{10}' (subscript) and
    'x^{2}' (superscript); script is None, 'sub' or 'sup'. Underscores without braces stay literal (pWER-KO_L1)."""
    out = []
    for part in _MARK.split(str(s)):
        if not part:
            continue
        if re.fullmatch(r'\*\*[^*]+\*\*', part):
            out.append((part[2:-2], 'normal', 'bold', None))
        elif re.fullmatch(r'\*[^*]+\*', part):
            out.append((part[1:-1], 'italic', 'normal', None))
        elif re.fullmatch(r'_\{[^{}]+\}', part):
            out.append((part[2:-1], 'normal', 'normal', 'sub'))
        elif re.fullmatch(r'\^\{[^{}]+\}', part):
            out.append((part[2:-1], 'normal', 'normal', 'sup'))
        else:                                   # literal text, including asterisks such as '****'
            out.append((part, 'normal', 'normal', None))
    return out


def script_size(size: float, script) -> float:
    return size if script is None else max(SCRIPT_SCALE * size, min(size, 5.0))


def script_shift(size: float, script) -> float:
    """Baseline shift in points (up positive) for a run."""
    return {None: 0.0, 'sub': -0.15 * size, 'sup': 0.38 * size}[script]


def parse_markup(s: str):
    """'FER/*fer-4*' -> [('FER/', 'normal', 'normal'), ('fer-4', 'italic', 'normal')]; **bold**; '-' kept as hyphen.
    (Sub/superscript runs are returned as plain runs here; Page.text draws them with parse_runs.)"""
    return [(t, fs, fw) for t, fs, fw, _ in parse_runs(s)]


def plain(s: str) -> str:
    return ''.join(t for t, _, _ in parse_markup(s))


@lru_cache(maxsize=None)
def labels() -> dict:
    p = KIT / 'config' / 'labels.yaml'
    if not p.exists():
        return {}
    with open(p) as f:
        return yaml.safe_load(f) or {}


def display(name: str, kind: str = 'genotypes') -> str:
    """Display markup for a data name (config/labels.yaml); the name itself when no entry exists."""
    return (labels().get(kind) or {}).get(name, name)
