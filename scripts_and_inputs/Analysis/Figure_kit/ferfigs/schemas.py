"""Tidy table schemas and plain-language validation.

Every data table the kit draws is one of these kinds. `python build.py validate` checks all of them and prints
each problem with the file, row and column, and how to fix it. The collaborator templates (templates/) use the
same column names.
"""
from __future__ import annotations

import difflib
import math
import re
from pathlib import Path

import pandas as pd

from .style import KIT

SCHEMAS = {
    # one row per individual measurement (one root, seedling, plate or image)
    'measurements': dict(
        required={'group': 'str', 'series': 'str', 'replicate': 'str', 'value': 'float'},
        optional={'sample_id': 'str', 'experiment': 'str', 'unit': 'str', 'notes': 'str'},
        doc='One row per individual measurement. group = x-axis category (genotype or time point); series = '
            'treatment or condition (colour); replicate = an identifier unique within group and series; value = the '
            'number plotted.'),
    # letters / asterisks / ns printed above bars, or on a bracket between two bars
    'annotations': dict(
        required={'group': 'str', 'series': 'str', 'label': 'str', 'kind': 'str'},
        optional={'group2': 'str', 'series2': 'str', 'test': 'str', 'p_value': 'float', 'source': 'str', 'notes': 'str',
                  'experimental_unit': 'str'},
        allowed={'kind': ['letters', 'stars', 'ns', 'text'],
                 'source': ['computed_by_authors', 'transcribed_from_deck_chart', 'proposed_by_helper',
                            'placeholder', '']},
        doc='Statistics labels exactly as they should be printed. The kit draws these and computes nothing. Leave '
            'group2/series2 empty for a label over a single bar; fill them for a bracket between two bars. test = the '
            'test that produced the label (for the legend and the Source Data file). experimental_unit = what n counts '
            'in that test (e.g. roots, seedlings, plates, independent experiments).'),
    # one row per image tile placed in a panel
    'image_manifest': dict(
        required={'panel': 'str', 'tile': 'str', 'file': 'str', 'row': 'int', 'col': 'int'},
        optional={'genotype': 'str', 'treatment': 'str', 'time': 'str', 'channel': 'str', 'orient': 'str',
                  'crop_px': 'str', 'um_per_px': 'float', 'scale_bar_um': 'float', 'scale_bar_px': 'str',
                  'overlays': 'str', 'px_origin': 'str',
                  'acquisition': 'str', 'raw_file': 'str', 'source': 'str', 'notes': 'str'},
        allowed={'orient': ['identity', 'mirrorLR', 'mirrorTB', 'rot90ccw', 'rot180', 'rot270ccw',
                            'mirrorLR+rot90ccw', 'mirrorLR+rot270ccw', ''],
                 'px_origin': ['placed', 'file', '']},
        doc='One row per image tile. file is relative to the figure data folder. crop_px = "x0 y0 x1 y1" in pixels '
            '(after orient). um_per_px = microscope calibration (micrometres per pixel) used for the scale bar; '
            'scale_bar_um = bar length to draw on this tile; scale_bar_px = "x0 y" left end of the bar in pixels '
            '(default: bottom right). Overlay and scale-bar pixels are counted in the placed (cropped) image by default; set '
            'px_origin = file to count them in the full oriented file, so editing crop_px keeps them on the tissue. '
            'overlays = boxes/arrows in image pixels, e.g. "box 10 20 60 80; arrow 30 40 45 40". No intensity '
            'changes are made by the kit.'),
}


_UNIT_SYNONYMS = {
    'um': 'µm', 'μm': 'µm', 'µm': 'µm', 'micrometre': 'µm', 'micrometer': 'µm', 'micrometres': 'µm',
    'micrometers': 'µm', 'micron': 'µm', 'microns': 'µm',
    'mm': 'mm', 'millimetre': 'mm', 'millimeter': 'mm', 'millimetres': 'mm', 'millimeters': 'mm',
    'cm': 'cm', 'centimetre': 'cm', 'centimeter': 'cm',
    'a.u.': 'a.u.', 'a.u': 'a.u.', 'au': 'a.u.', 'arbitrary units': 'a.u.', 'arbitrary unit': 'a.u.',
}


def canonical_unit(u) -> str:
    """Unit text compared without case and spacing; common spellings of the same unit are equal (um = µm)."""
    s = ' '.join(str(u).strip().split()).lower()
    return _UNIT_SYNONYMS.get(s, s)


def unit_problems(df, unit, name):
    """Every filled unit cell must be the figure's unit; the kit never converts. Empty cells mean the figure's unit."""
    if not unit or 'unit' not in df.columns:
        return []
    want = canonical_unit(unit)
    bad = {}
    for i, v in df['unit'].items():
        if str(v).strip() and canonical_unit(v) != want:
            bad.setdefault(str(v).strip(), []).append(str(i + 2))
    return [f"{name} line(s) {', '.join(rows[:8])}: unit '{u}' is not the figure's unit '{unit}'. Convert these "
            f"values to {unit} before importing (the kit does not convert units), or change the unit in the figure "
            f"config if the whole chart should be in {u}." for u, rows in sorted(bad.items())]


class ValidationError(Exception):
    def __init__(self, problems):
        self.problems = problems
        super().__init__('\n'.join(problems))


def _suggest(value, options):
    m = difflib.get_close_matches(str(value), [str(o) for o in options], n=1, cutoff=0.6)
    return f" Did you mean '{m[0]}'?" if m else ''


def read_table(path, kind, levels=None, unit=None, allow_proposed=False):
    """Read a CSV (or the first sheet of an .xlsx) and check it against a schema. Returns a DataFrame.
    levels: optional {column: allowed names} (the order lists of the figure config), checked in the same pass.
    unit: the figure's unit for a measurements table; every filled unit cell must equal it.
    allow_proposed: statistics proposed by the helper (source proposed_by_helper, *_PROPOSED files) are refused
    unless the figure config opts in explicitly.
    Raises ValidationError listing every problem found."""
    p = Path(path)
    if not p.is_absolute():
        p = KIT / p
    if not p.exists():
        raise ValidationError([f'{p}: file not found. Check the path in the figure config (config/*.yaml).'])
    if p.suffix.lower() in ('.xlsx', '.xls'):
        df = pd.read_excel(p, dtype=str)
    else:
        df = pd.read_csv(p, dtype=str, keep_default_na=False, encoding='utf-8-sig')
    df.columns = [c.strip() for c in df.columns]
    sch = SCHEMAS[kind]
    probs = []
    missing = [c for c in sch['required'] if c not in df.columns]
    if missing:
        probs.append(f"{p.name}: missing column(s) {missing}. Required columns for a {kind} table: "
                     f"{list(sch['required'])}.")
        raise ValidationError(probs)
    unknown = [c for c in df.columns if c not in sch['required'] and c not in sch['optional']]
    if unknown:
        probs.append(f"{p.name}: unknown column(s) {unknown} (ignored). Allowed optional columns: "
                     f"{list(sch['optional'])}." + ''.join(_suggest(u, list(sch['optional'])) for u in unknown))
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip().replace({'nan': ''})
    types = {**sch['required'], **sch['optional']}
    for c, t in types.items():
        if c not in df.columns:
            continue
        for i, v in df[c].items():
            row = i + 2       # header is line 1
            if v == '':
                if c in sch['required']:
                    probs.append(f"{p.name} line {row}: column '{c}' is empty.")
                continue
            if t == 'float':
                try:
                    x = float(v.replace('\u2212', '-'))
                    if not math.isfinite(x):
                        raise ValueError
                except ValueError:
                    probs.append(f"{p.name} line {row}: '{c}' = '{v}' is not a number.")
            elif t == 'int':
                if not re.fullmatch(r'-?\d+', v):
                    probs.append(f"{p.name} line {row}: '{c}' = '{v}' is not a whole number.")
    for c, opts in sch.get('allowed', {}).items():
        if c in df.columns:
            for i, v in df[c].items():
                if v not in opts:
                    probs.append(f"{p.name} line {i + 2}: '{c}' = '{v}' is not one of {[o for o in opts if o]}."
                                 + _suggest(v, opts))
    for col, allowed in (levels or {}).items():
        if col in df.columns:
            bad = [(i, v) for i, v in df[col].items() if v != '' and v not in allowed]
            for v in sorted({v for _, v in bad}):
                rows = [str(i + 2) for i, w in bad if w == v][:8]
                probs.append(f"{p.name} line(s) {', '.join(rows)}: {col} '{v}' is not in the figure config list "
                             f"{list(allowed)}." + _suggest(v, allowed))
    if kind == 'measurements':
        dup = df.duplicated(['group', 'series', 'replicate'], keep=False)
        if dup.any():
            rows = [str(i + 2) for i in df.index[dup]][:10]
            probs.append(f"{p.name}: the same group/series/replicate appears more than once (lines {', '.join(rows)}). "
                         f"Give each measurement its own replicate id.")
    if kind == 'measurements':
        probs += unit_problems(df, unit, p.name)
    if kind == 'annotations' and not allow_proposed:
        if p.stem.upper().endswith('_PROPOSED') or ('source' in df.columns and (df['source'] == 'proposed_by_helper').any()):
            probs.append(f"{p.name}: contains statistics proposed by the helper (source proposed_by_helper). The kit "
                         f"draws only labels the authors computed: review them, set source to computed_by_authors with "
                         f"the test, or opt in explicitly with allow_proposed_stats: true in the panel config.")
    if kind == 'annotations':
        for i, r in df.iterrows():
            lab, kd = r['label'], r['kind']
            ok = {'letters': re.fullmatch(r'[a-z]+', lab), 'stars': re.fullmatch(r'\*{1,4}', lab),
                  'ns': lab == 'ns', 'text': True}.get(kd, True)
            if not ok:
                probs.append(f"{p.name} line {i + 2}: label '{lab}' does not look like kind '{kd}'.")
    hard = [q for q in probs if 'ignored' not in q]
    if hard:
        raise ValidationError(probs)
    if kind == 'measurements':
        df['value'] = df['value'].str.replace('\u2212', '-').astype(float)
    if kind == 'image_manifest':
        df['row'] = df['row'].astype(int)
        df['col'] = df['col'].astype(int)
    try:                                           # messages name kit files (data/...), not this copy's location
        df.attrs['path'] = p.resolve().relative_to(KIT.resolve()).as_posix()
    except ValueError:
        df.attrs['path'] = str(p)
    df.attrs['warnings'] = probs
    return df


def check_levels(df, column, allowed, path=''):
    """Every value of df[column] must be one of `allowed` (the order list in the figure config)."""
    bad = sorted(set(df[column]) - set(allowed))
    if bad:
        raise ValidationError([f"{Path(path).name}: {column} value(s) {bad} are not in the figure config order list "
                               f"{list(allowed)}." + ''.join(_suggest(b, allowed) for b in bad)
                               + ' Fix the table or add the name to the config.'])
