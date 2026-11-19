"""Fig. 6d: pYUC5-LUC activity heatmap, sideways (effector columns x replicate rows), port of draw_luc_sideways.py.

Values: data/fig6/luc.csv when present (measurements; see from_measurements), otherwise the displayed values
data/fig6/fig6d_luc_displayed_values.csv (the numbers printed in the original figure, checked against the native
PowerPoint table by tools/seed_fig6_from_sources.py). Colour = log10 fold change, white -> blue, linear in CIELAB
over [0, vmax] (ggplot2 scale_fill_gradient); the colour bar is one continuous gradient image. Cell text is white or
black, whichever contrasts more with the tile (WCAG 2, >= 4.5:1 checked).
"""
from __future__ import annotations

import math
import re

import numpy as np
import pandas as pd

from ..schemas import ValidationError
from .fig6_common import MM, PT, cap_mm, contrast_text, gradient_image, lab_ramp


def _rep(v):
    """'Rep2' or '2' -> 'Rep2'; anything else (e.g. '2.5', 'R2', blank) -> None, so it is reported, never coerced."""
    v = str(v).strip()
    m = re.fullmatch(r'(?:Rep)?([1-9][0-9]*)', v)
    return f'Rep{m.group(1)}' if m else None


def from_measurements(path, order, reps):
    """luc.csv: effector, replicate and EITHER luc + ren (one row per well) OR fold_change.
    With luc/ren: fold change = mean(luc/ren) of the effector in a replicate / mean Control (luc/ren) of the same
    replicate (the first effector in the config order is the control)."""
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    probs = []
    if not {'effector', 'replicate'} <= set(df.columns):
        probs.append(f'{path}: need columns effector, replicate and either luc + ren or fold_change')
    raw = {'luc', 'ren'} <= set(df.columns) and (df.get('luc', pd.Series(dtype=str)).str.strip() != '').any()
    if not raw and 'fold_change' not in df.columns:
        probs.append(f'{path}: need luc + ren (raw readings) or fold_change')
    if probs:
        raise ValidationError(probs)
    acc = {}
    for i, r in df.iterrows():
        eff, rep = r['effector'].strip(), _rep(r['replicate'])
        if eff not in order:
            probs.append(f'{path} line {i + 2}: effector {eff!r} not in the config order {order}')
            continue
        if rep is None or rep not in reps:     # never drop or merge a supplied replicate silently
            probs.append(f"{path} line {i + 2}: replicate {r['replicate']!r} is not one of {list(reps)}. Add it to "
                         f"panels.d.replicates in config/fig6.yaml (one heatmap row per replicate) or correct the id.")
            continue
        try:
            v = float(r['luc']) / float(r['ren']) if raw else float(r['fold_change'])
        except (ValueError, ZeroDivisionError):
            probs.append(f'{path} line {i + 2}: not a number')
            continue
        if not v > 0:
            probs.append(f'{path} line {i + 2}: values must be positive')
        acc.setdefault((eff, rep), []).append(v)
    for e in order:
        for rp in reps:
            if (e, rp) not in acc:
                probs.append(f'{path}: no value for {e} {rp}')
    if probs:
        raise ValidationError(probs)
    out = {}
    for e in order:
        for rp in reps:
            m = float(np.mean(acc[(e, rp)]))
            fc = m / float(np.mean(acc[(order[0], rp)])) if raw else m
            out[(e, rp)] = (fc, f'{fc:.1f}')
    how = ('luc/ren per well, fold change = effector mean / control mean within replicate' if raw
           else 'fold_change as supplied')
    return out, how


def from_displayed(path, order, reps):
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    out = {}
    for i, r in df.iterrows():
        rep = _rep(r['replicate'])
        if rep is None or rep not in reps:
            raise ValidationError([f"{path} line {i + 2}: replicate {r['replicate']!r} is not one of {list(reps)}"])
        out[(r['effector'], rep)] = (float(r['fold_change']), r['label'])
    miss = [f'{e} {rp}' for e in order for rp in reps if (e, rp) not in out]
    if miss:
        raise ValidationError([f'{path}: missing {miss}'])
    for (e, rp), (v, lab) in out.items():
        if f'{v:.1f}' != lab:
            raise ValidationError([f'{path}: label {lab!r} is not the value {v} at one decimal ({e} {rp})'])
    return out, 'displayed values (printed in the original figure)'


def draw(page, mm: MM, box, vals, cfg, gid='fig6-d'):
    st = page.st
    fs = cfg.get('text_pt', st['font']['size_pt'])
    fs_s = cfg.get('small_pt', st['font']['small_pt'])        # colour-bar tick labels
    x, y, w, h = box
    order, reps = cfg['effectors'], cfg['replicates']
    vmax = float(cfg['vmax'])
    ramp = lab_ramp(cfg['colour_low'], cfg['colour_high'])
    lfc = {k: math.log10(v[0]) for k, v in vals.items()}
    if max(lfc.values()) > vmax + 1e-9:
        raise ValidationError([f'd: log10 fold change {max(lfc.values()):.2f} is above the colour-scale maximum '
                               f'{vmax}; raise d.vmax in config/fig6.yaml.'])
    hx0 = x + cfg['x_tiles']
    ch, g = cfg['cell_h'], cfg['gutter']
    if 'right_pad' in cfg:                       # tiles fill the box up to right_pad (colour bar and labels there)
        cw = (w - cfg['right_pad'] - cfg['x_tiles']) / len(order)
    else:
        cw = cfg['cell_w']
    hy0 = y + cfg['y_tiles']
    cap = cap_mm(fs)
    # title (bold; *pYUC5* italic)
    page.text(hx0 + g / 2, y + cfg['y_title'], cfg['title'], size=fs, weight='bold', gid=f'{gid}-title')
    contrasts, widest = [], 0.0
    for j, e in enumerate(order):
        for i, rp in enumerate(reps):
            v, lab = vals[(e, rp)]
            fill = ramp(lfc[(e, rp)] / vmax)
            tc, cr = contrast_text(fill)
            contrasts.append(cr)
            x0, y0 = hx0 + j * cw + g / 2, hy0 + i * ch + g / 2
            white = fill == (255, 255, 255)
            mm.rbox(x0, y0, cw - g, ch - g, cfg['radius'], face=[c / 255 for c in fill],
                    edge=cfg['white_tile_edge'] if white else None, lw=cfg['white_tile_edge_pt'] if white else 0.0,
                    gid=f'{gid}-tile', z=1)
            page.text(x0 + (cw - g) / 2, y0 + (ch - g) / 2 + cap / 2, lab, size=fs, ha='center', colour=tc,
                      gid=f'{gid}-value')
            widest = max(widest, page.text_width(lab, fs))
    if min(contrasts) < 4.5:
        raise ValidationError([f'd: text contrast {min(contrasts):.2f} below 4.5:1'])
    if widest > cw - g - 2 * 0.5:
        raise ValidationError([f'd: value labels ({widest:.2f} mm) do not fit the tiles; widen d.cell_w.'])
    for i, rp in enumerate(reps):
        page.text(hx0 - 0.9, hy0 + i * ch + ch / 2 + cap / 2, cfg['replicate_labels'].get(rp, rp), size=fs,
                  ha='right', gid=f'{gid}-rows')
    hy1 = hy0 + len(reps) * ch
    for j, e in enumerate(order):
        page.text(hx0 + (j + 0.5) * cw + 0.35, hy1 + cfg['name_gap'], cfg.get('effector_labels', {}).get(e, e),
                  size=fs, rotation=cfg['name_angle'], ha='right', va='top', gid=f'{gid}-names')
    # continuous colour bar
    bx = hx0 + len(order) * cw + cfg['bar_gap']
    by0, by1 = hy0 + g / 2, hy1 - g / 2
    img = gradient_image(lambda t: ramp(t / vmax), cfg['bar_w'], by1 - by0, ppi=st['output']['image_max_ppi'],
                         top_value=vmax, bottom_value=0.0)
    page.image(img, bx, by0, w=cfg['bar_w'], gid=f'{gid}-colourbar', name='fig6d_colour_bar')
    mm.rect(bx, by0, cfg['bar_w'], by1 - by0, face='none', edge=cfg['white_tile_edge'], lw=cfg['bar_frame_pt'],
            gid=f'{gid}-colourbar', z=3)
    for v in cfg['bar_ticks']:
        yy = by1 - v / vmax * (by1 - by0)
        mm.line(bx + cfg['bar_w'], yy, bx + cfg['bar_w'] + 0.5, yy, lw=cfg['bar_frame_pt'],
                colour=cfg.get('bar_tick_colour', 'k'), gid=f'{gid}-colourbar', z=3)
        page.text(bx + cfg['bar_w'] + 0.8, yy + cap_mm(fs_s) / 2, f'{v:g}', size=fs_s, gid=f'{gid}-colourbar')
    # bar title on the title baseline, right-aligned: markup string ('log_{10}FC'; the kit sets the subscript at
    # 0.85 x, never below 5 pt) or the older [log, 10, FC] runs
    sub = cfg['bar_title']
    yt = y + cfg['y_title']
    if isinstance(sub, str):
        page.text(x + w - cfg.get('bar_title_right_pad', 0.0), yt, sub, size=fs, ha='right', gid=f'{gid}-bartitle')
    else:
        w_all = page.text_width(sub[0], fs) + page.text_width(sub[1], 5.0) + page.text_width(sub[2], fs)
        xr = min(x + w, bx + cfg['bar_w'] + 4.0) - w_all
        xr += page.text(xr, yt, sub[0], size=fs, gid=f'{gid}-bartitle')
        xr += page.text(xr, yt + cfg['subscript_drop'], sub[1], size=5.0, gid=f'{gid}-bartitle')
        page.text(xr, yt, sub[2], size=fs, gid=f'{gid}-bartitle')
    return dict(type='heatmap', effectors=order, replicates=reps, vmax=vmax,
                values={f'{e}|{rp}': vals[(e, rp)][1] for e in order for rp in reps},
                min_text_contrast=round(min(contrasts), 2), colour_bar_px=list(img.shape[:2]))
