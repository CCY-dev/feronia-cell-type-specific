"""Charts for Fig. 1, drawn with Page primitives (vector rectangles, lines and text; no axes text).

  count_bars        b: vertical bars with the value printed over every bar, highlighted bins, centred summary text
  legend_key        a key: coloured dots + names in columns under a bold title
  rows_bars_heatmap c + d: one row band per subpopulation shared by a stacked horizontal bar chart (c) and a heatmap
                    of counts coloured by a statistic (d), with a continuous colour key
Numbers drawn are exactly the table values; nothing is computed here except positions.
"""
from __future__ import annotations

import numpy as np

from .. import style as S
from .fig1_colour import ColorRamp2, hex_to_rgb01, luminance_contrast_black

MINUS = '\u2212'
ASC = 0.716     # Arial ascender (fraction of size)
DESC = 0.21


def _mm(pt):
    return pt / S.PT_PER_MM


def fmt_int(v):
    return f'{int(v):,}'


def axis_ticks(vmax, step):
    return np.arange(0, vmax + 1e-9, step)


# ------------------------------------------------------------------ b
def count_bars(page, box, labels, values, spec, geo, lines, gid='fig1-b'):
    """box = (x0, y_top, x1, y_axis) of the plot frame in mm (y_top = level of y_max). labels: x tick labels (markup);
    values: bar heights (ints). geo: text offsets (mm) and sizes; lines: line widths and tick length."""
    st = page.st
    x0, yt, x1, ya = box
    size, small, title_pt = st['font']['size_pt'], st['font']['small_pt'], geo['title_pt']
    n = len(values)
    ymax, ystep = float(spec['y_max']), float(spec['y_step'])
    if max(values) > ymax:
        from ..schemas import ValidationError
        raise ValidationError([f'fig1b: largest value {max(values)} is above y_max {ymax:g}; raise y_max in the config.'])
    pad = geo['x_pad_slots']
    slot = (x1 - x0) / (n + 2 * pad)
    xc = [x0 + slot * (pad + i + 0.5) for i in range(n)]

    def Y(v):
        return ya - v / ymax * (ya - yt)
    hi = set(spec['highlight_bins'])
    bw = slot * spec.get('bar_width', 0.8)
    tops = []
    for i, (lab_key, v) in enumerate(zip(spec['_bins'], values)):
        face = spec['colours']['highlight'] if lab_key in hi else spec['colours']['other']
        page.rect(xc[i] - bw / 2, Y(v), bw, ya - Y(v), lw=0, edge='none', face=S.rgb(face), gid=f'{gid}-bar', zorder=3)
        page.text(xc[i], Y(v) - geo['value_gap'], fmt_int(v), size=small, ha='center', gid=f'{gid}-values')
        tops.append(Y(v))
    # axes: left spine with ticks, bottom spine with a tick under every bar
    tl = _mm(lines['tick_len_pt'])
    page.line(x0, yt, x0, ya, lw=lines['axis_pt'], gid=f'{gid}-axis', cap='projecting')
    page.line(x0, ya, x1, ya, lw=lines['axis_pt'], gid=f'{gid}-axis', cap='projecting')
    for t in axis_ticks(ymax, ystep):
        page.line(x0 - tl, Y(t), x0, Y(t), lw=lines['tick_pt'], gid=f'{gid}-axis')
        page.text(x0 - tl - geo['ytick_pad'], Y(t) + geo['ytick_dy'], fmt_int(t), ha='right', gid=f'{gid}-axis')
    for i, lab in enumerate(labels):
        page.line(xc[i], ya, xc[i], ya + tl, lw=lines['tick_pt'], gid=f'{gid}-axis')
        page.text(xc[i], ya + geo['xtick_dy'], lab, ha='center', gid=f'{gid}-axis')
    tw = max(page.text_width(fmt_int(t)) for t in axis_ticks(ymax, ystep))
    y_label_x = x0 - tl - geo['ytick_pad'] - tw - geo['ylabel_gap']
    page.text(y_label_x, (yt + ya) / 2, spec['y_label'], size=title_pt, rotation=90, ha='center', gid=f'{gid}-axis')
    x_title_base = ya + geo['xtitle_dy']
    page.text((x0 + x1) / 2, x_title_base, spec['x_label'], size=title_pt, ha='center', gid=f'{gid}-axis')
    return dict(xc=xc, tops=tops, slot=slot, bar_w=bw, Y=Y, bottom=x_title_base + _mm(DESC * title_pt),
                y_label_x=y_label_x - _mm(ASC * title_pt))


# ------------------------------------------------------------------ key for a
def legend_key(page, x, y_title, title, entries, geo, col_dx, rows_per_col, edge=None, gid='fig1-key'):
    """entries: [(markup, '#rrggbb')]. Title baseline at y_title; entry i in column i // rows_per_col, its baseline
    (r + 1) pitches below the title. geo: pitch, marker_d, marker_cx / marker_cy (dot centre right of x / above the
    baseline), text_dx, title_pt, entry_pt. edge: None (no outline) or [grey, width_pt]."""
    from matplotlib.patches import Ellipse
    page.text(x, y_title, title, size=geo['title_pt'], weight='bold', gid=gid)
    d = geo['marker_d']
    for i, (lab, hexc) in enumerate(entries):
        c, r = divmod(i, rows_per_col)
        xb = x + col_dx[c]
        yb = y_title + (r + 1) * geo['pitch']
        # a circle in page mm: an ellipse in figure fractions (the page is not square)
        circ = Ellipse((page.fx(xb + geo['marker_cx']), page.fy(yb - geo['marker_cy'])), d / page.W, d / page.H,
                       transform=page.fig.transFigure, facecolor=hex_to_rgb01(hexc),
                       edgecolor=(edge[0],) * 3 if edge else 'none', lw=edge[1] if edge else 0, zorder=5)
        circ.set_gid(gid)
        page.fig.add_artist(circ)
        page.text(xb + geo['text_dx'], yb, lab, size=geo['entry_pt'], gid=gid)
    n_rows = min(len(entries), rows_per_col)
    return y_title + n_rows * geo['pitch'] + _mm(DESC * geo['entry_pt'])


# ------------------------------------------------------------------ c + d
def rows_bars_heatmap(page, geo, c_rows, d_counts, d_values, scale, spec_c, spec_d, lines, gid='fig1-cd'):
    """c_rows: [(label markup, up, down, total)] in display order; d_counts / d_values: arrays (rows x hormones);
    scale: (breaks, colours) of the colour ramp; geo: layout.cd of the config (mm); lines: line widths.
    Vertical positions are offsets from the top of the first row (the c axis line). Returns geometry and what was
    drawn."""
    st = page.st
    size, small, title_pt = st['font']['size_pt'], st['font']['small_pt'], geo['title_pt']
    n = len(c_rows)
    p = geo['row_pitch']
    y_rows = geo['y_top'] + geo['rows_top']
    em = geo['row_text_em']                  # row texts: baseline = row centre + em x font size
    tl = _mm(lines['tick_len_pt'])

    def RY(i):              # centre of row i
        return y_rows + (i + 0.5) * p
    # alternating bands across the label column, c and d
    x_band0, x_band1 = geo['band_x0'], geo['d_x0'] + geo['tile_w'] * len(spec_d['hormones']) + geo['band_pad_right']
    for i in range(0, n, 2):
        page.rect(x_band0, y_rows + i * p, x_band1 - x_band0, p, lw=0, edge='none',
                  face=S.rgb(spec_c['band_colour']), gid=f'{gid}-band', zorder=1)
    # row labels (right-aligned)
    for i, (lab, *_r) in enumerate(c_rows):
        page.text(geo['label_right'], RY(i) + _mm(em * size), lab, ha='right', gid=f'{gid}-labels')
    # ---- c: stacked bars (down first, then up), totals at the bar ends, axis on top
    cx0, cw, cmax = geo['c_x0'], geo['c_w'], float(spec_c['x_max'])

    def X(v):
        return cx0 + v / cmax * cw
    if max(r[3] for r in c_rows) > cmax:
        from ..schemas import ValidationError
        raise ValidationError([f'fig1c: largest total above x_max {cmax:g}; raise x_max in the config.'])
    bh = p * spec_c['bar_height']
    col_dn, col_up = S.rgb(spec_c['colours']['down']), S.rgb(spec_c['colours']['up'])
    for i, (_l, up, dn, tot) in enumerate(c_rows):
        yb = RY(i) - bh / 2
        if dn:
            page.rect(X(0), yb, X(dn) - X(0), bh, lw=0, edge='none', face=col_dn, gid=f'{gid}-c-down', zorder=3)
        if up:
            page.rect(X(dn), yb, X(dn + up) - X(dn), bh, lw=0, edge='none', face=col_up, gid=f'{gid}-c-up', zorder=3)
        page.text(X(tot) + geo['total_gap'], RY(i) + _mm(em * small), fmt_int(tot), size=small, gid=f'{gid}-c-totals')
    page.line(X(0), y_rows, X(cmax), y_rows, lw=lines['axis_pt'], gid=f'{gid}-c-axis', cap='projecting')
    for t in axis_ticks(cmax, float(spec_c['x_step'])):
        page.line(X(t), y_rows - tl, X(t), y_rows, lw=lines['tick_pt'], gid=f'{gid}-c-axis')
        page.text(X(t), y_rows + geo['c_tick_dy'], fmt_int(t), ha='center', gid=f'{gid}-c-axis')
    page.text(X(0) + cw / 2, y_rows + geo['c_title_dy'], spec_c['title'], size=title_pt, ha='center',
              gid=f'{gid}-c-axis')
    # up/down key over the label column: swatches at a fixed pitch, names on one baseline
    k = geo['c_key']
    for j, (lab, col) in enumerate(((spec_c['key_labels']['down'], col_dn), (spec_c['key_labels']['up'], col_up))):
        kx = k['x'] + j * k['pitch']
        page.rect(kx, y_rows + k['swatch_dy'], k['swatch'][0], k['swatch'][1], lw=0, edge='none', face=col,
                  gid=f'{gid}-c-key')
        page.text(kx + k['swatch'][0] + k['text_dx'], y_rows + k['text_dy'], lab, gid=f'{gid}-c-key')
    # ---- d: heatmap tiles with counts
    ramp = ColorRamp2(*scale)
    hs = spec_d['hormones']
    dx0, tw = geo['d_x0'], geo['tile_w']
    th = p * spec_d['tile_height']
    fills = np.array(ramp(d_values.ravel())).reshape(d_values.shape)
    edge = S.rgb(spec_d['tile_edge'])
    for i in range(n):
        for j in range(len(hs)):
            fx = dx0 + j * tw + tw * (1 - spec_d['tile_width']) / 2
            rgb = hex_to_rgb01(fills[i, j])
            page.rect(fx, RY(i) - th / 2, tw * spec_d['tile_width'], th, lw=spec_d['tile_edge_pt'], edge=edge,
                      face=rgb, gid=f'{gid}-d-tile', zorder=3)
            if luminance_contrast_black(rgb) < 4.5:
                raise ValueError(f'black count text would be illegible on {fills[i, j]} (contrast below 4.5)')
            page.text(dx0 + (j + 0.5) * tw, RY(i) + _mm(em * small), str(int(d_counts[i, j])), size=small,
                      ha='center', gid=f'{gid}-d-counts')
    for j, h in enumerate(hs):
        page.text(dx0 + (j + 0.5) * tw, y_rows + geo['d_label_dy'], spec_d['hormone_labels'].get(h, h), ha='center',
                  gid=f'{gid}-d-axis')
    page.text(dx0 + tw * len(hs) / 2, y_rows + geo['d_title_dy'], spec_d['title'], size=title_pt, ha='center',
              gid=f'{gid}-d-axis')
    # ---- colour key: continuous bar (raster at ppi), vector ticks, title (markup, e.g. '-log_{10} *P*')
    k = geo['d_key']
    kx = dx0 + tw * len(hs) + k['gap']
    kmid = y_rows + n * p * k['center_frac']
    kpx_w, kpx_h = int(round(k['w'] / 25.4 * k['ppi'])), int(round(k['h'] / 25.4 * k['ppi']))
    kw, kh = kpx_w * 25.4 / k['ppi'], kpx_h * 25.4 / k['ppi']
    ky0 = kmid - kh / 2
    lo, hi_ = scale[0][0], scale[0][-1]
    vals = hi_ - (np.arange(kpx_h) + 0.5) / kpx_h * (hi_ - lo)
    grad = np.array([[int(h_[i:i + 2], 16) for i in (1, 3, 5)] for h_ in ramp(vals)], np.uint8)
    page.image(np.repeat(grad[:, None, :], kpx_w, axis=1), kx, ky0, w=kw, gid=f'{gid}-d-key', name='fig1d_key')
    grey = (k['frame_grey'],) * 3
    page.rect(kx, ky0, kw, kh, lw=k['frame_pt'], edge=grey, gid=f'{gid}-d-key')
    ticks = sorted(set(list(scale[0]) + list(k.get('extra_ticks', []))))
    labs = []
    for t in ticks:
        yy = ky0 + kh - (t - lo) / (hi_ - lo) * kh
        page.line(kx + kw, yy, kx + kw + _mm(k['tick_len_pt']), yy, lw=lines['tick_pt'], colour=grey,
                  gid=f'{gid}-d-key')
        lab = f'{t:.1f}'.rstrip('0').rstrip('.') if abs(t - round(t)) > 1e-9 else f'{t:g}'
        if abs(t - hi_) < 1e-9:
            lab = k['cap_prefix'] + lab
        page.text(kx + kw + k['label_dx'], yy + k['label_dy'], lab, gid=f'{gid}-d-key')
        labs.append(lab)
    # title ends level with the right end of the widest tick label (approved page), on a baseline above the bar
    x_right = kx + kw + k['label_dx'] + max(page.text_width(l_) for l_ in labs)
    page.text(x_right, ky0 + k['title_dy'], spec_d['key_title'], ha='right', gid=f'{gid}-d-key')
    drawn = dict(rows=n, hormones=list(hs), fills=[[fills[i, j] for j in range(len(hs))] for i in range(n)])
    return dict(y_rows=y_rows, bottom=y_rows + n * p, d_x0=dx0, key_x=kx, key_box=(kx, ky0, kw, kh),
                key_right=x_right, drawn=drawn, RY=RY)
