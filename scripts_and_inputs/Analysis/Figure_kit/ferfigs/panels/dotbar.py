"""Bar (mean) + error bar (SD or SEM) + every individual value, grouped by x category and coloured by series or
by group. Statistics labels come from an annotations table and are drawn exactly as given; nothing is computed
here except mean, SD/SEM and n for the bars.

spec (from the figure config):
  groups: [..]            x categories in order (names as in the data)
  series: [..]            series in order within each group (use ['all'] for one series)
  colour_by: series|group
  y_label: 'text' or ['line 1', 'line 2'] (markup: *italic*)
  y_max: 60  y_step: 20   (optional: automatic when missing)
  error: sd|sem|none      (default sd)
  legend: true|false      (series key above the plot; default true when colour_by is series)
  x_rotation: 45|0
  group_gap: 0.7          (empty slots between groups)
  star_rotation: 90|0     (asterisks above narrow bars read bottom to top, as in Prism)
  letter_rotation: 90|0   (letter strings longer than letter_rotation_min_len read bottom to top; default 0)
  require_labels: true    (every bar must have a letters label)
  bar_width: 0.8          (bar width as a fraction of one slot; default style bars.width)
  stats_pt: 7.0           (size of letters / asterisks; default style font.stats_pt, else font.size_pt)
"""
from __future__ import annotations

import math

import numpy as np

from .. import style as S
from ..schemas import ValidationError

PHI = 0.6180339887498949


def nice_axis(vmax, n=5):
    raw = vmax / n
    mag = 10 ** math.floor(math.log10(raw))
    step = min((s * mag for s in (1, 2, 2.5, 5, 10) if s * mag >= raw), default=raw)
    return step * math.ceil(vmax / step - 1e-9), step


def jitter(n, width):
    """Deterministic, value-independent horizontal spread (golden-ratio sequence), centred on 0."""
    if n == 1:
        return np.zeros(1)
    k = np.arange(n)
    return ((k * PHI) % 1 - 0.5) * width


def draw(page, box, df, spec, ann=None, gid='panel'):
    st = page.st
    x, y, w, h = box
    groups, series = list(spec['groups']), list(spec['series'])
    S_ = len(series)
    gap = spec.get('group_gap', 0.7)
    colour_by = spec.get('colour_by', 'series')
    err_kind = spec.get('error', 'sd')
    bw = spec.get('bar_width', st['bars']['width'])   # per-chart override
    # ---- check the data against the config
    probs = []
    for col, allowed in (('group', groups), ('series', series)):
        bad = sorted(set(df[col]) - set(allowed))
        if bad:
            probs.append(f"{df.attrs.get('path', 'data')}: {col} {bad} not in the config list {allowed}.")
    if probs:
        raise ValidationError(probs)
    # ---- bars
    bars = []
    for gi, g in enumerate(groups):
        for si, s in enumerate(series):
            v = df[(df['group'] == g) & (df['series'] == s)]['value'].to_numpy(float)
            if len(v) == 0:
                continue
            xc = gi * (S_ + gap) + si
            mean = float(v.mean())
            sd = float(v.std(ddof=1)) if len(v) > 1 else 0.0
            e = {'sd': sd, 'sem': sd / math.sqrt(len(v)) if len(v) > 1 else 0.0, 'none': 0.0}[err_kind]
            if colour_by == 'group':
                fill, stroke = S.genotype_colours(g)
            else:
                fill = stroke = S.series_colour(s)
            bars.append(dict(group=g, series=s, x=xc, values=v, n=len(v), mean=mean, sd=sd, err=e,
                             fill=fill, stroke=stroke))
    if not bars:
        raise ValidationError([f'{df.attrs.get("path", "data")}: no rows match the config groups/series.'])
    vmax = max(max(b['values'].max(), b['mean'] + b['err']) for b in bars)
    if spec.get('y_max') is None:
        ymax, ystep = nice_axis(vmax * 1.05)
    else:
        ymax, ystep = float(spec['y_max']), float(spec.get('y_step') or nice_axis(float(spec['y_max']))[1])
    ymin = float(spec.get('y_min', 0.0))
    if vmax > ymax + 1e-9:
        raise ValidationError([f"{df.attrs.get('path', 'data')}: largest value {vmax:.3g} is above y_max {ymax:g}; "
                               f"raise y_max in the config."])
    xlo = -0.5 - 0.35
    xhi = (len(groups) - 1) * (S_ + gap) + S_ - 1 + 0.5 + 0.35
    ax = page.axes(x, y, w, h, gid=f'{gid}-axes')
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ymin, ymax)
    ticks = np.arange(ymin, ymax + ystep * 0.001, ystep)
    ax.set_yticks(ticks)
    fmt = spec.get('y_tick_format', '{:g}')
    ax.set_yticklabels([fmt.format(t).replace('-', '\u2212') for t in ticks])
    centres = [gi * (S_ + gap) + (S_ - 1) / 2 for gi in range(len(groups))]
    ax.set_xticks(centres)
    ax.set_xticklabels([])
    ax.spines['bottom'].set_bounds(xlo, xhi)
    ln = st['lines']
    msize = st['markers']['diameter_pt']
    for b in bars:
        ax.bar(b['x'], b['mean'], width=bw, color=b['fill'] + (st['bars']['alpha'],), edgecolor='black',
               linewidth=ln['bar_edge_pt'], zorder=1, gid=f'{gid}-bar')
        xs = b['x'] + jitter(b['n'], bw * 0.62)
        ax.scatter(xs, b['values'], s=msize ** 2,
                   facecolors='none', edgecolors=[b['stroke']], linewidths=ln['marker_edge_pt'], zorder=2,
                   clip_on=False, gid=f'{gid}-points')
        if err_kind != 'none' and b['n'] > 1:
            ax.errorbar(b['x'], b['mean'], yerr=b['err'], fmt='none', ecolor='black', elinewidth=ln['error_pt'],
                        capsize=1.6, capthick=ln['error_pt'], zorder=3)
    # ---- mm helpers
    def X(xd):
        return x + (xd - xlo) / (xhi - xlo) * w

    def Y(yd):
        return y + h - (yd - ymin) / (ymax - ymin) * h
    size = st['font']['size_pt']
    tick_mm = ln['tick_len_pt'] / S.PT_PER_MM
    # ---- x labels (markup, optionally rotated so that the label ends at its tick)
    rot = spec.get('x_rotation', 45)
    kind = spec.get('label_kind', 'genotypes')
    for gi, g in enumerate(groups):
        lab = S.display(g, kind)
        if rot:
            page.text(X(centres[gi]) + 0.2, y + h + tick_mm + 0.6, lab, rotation=rot, ha='right', va='top',
                      gid=f'{gid}-xlabels')
        else:
            page.text(X(centres[gi]), y + h + tick_mm + 0.5, lab, ha='center', va='top', gid=f'{gid}-xlabels')
    # ---- y label left of the tick labels
    tl_w = max(page.text_width(fmt.format(t).replace('-', '\u2212'), size) for t in ticks)
    ylab = spec.get('y_label', '')
    lines = ylab if isinstance(ylab, list) else [ylab]
    xl = x - tick_mm - 0.5 - tl_w - 0.9
    for i, line in enumerate(lines[::-1]):
        page.text(xl - i * size * 1.18 / S.PT_PER_MM, y + h / 2, line, rotation=90, ha='center', va='baseline',
                  gid=f'{gid}-ylabel')
    # ---- series key above the plot
    if spec.get('legend', colour_by == 'series') and series != ['all']:
        short = spec.get('legend_labels') or {}
        items = [(short.get(s, S.display(s, 'series')), S.series_colour(s)) for s in series]
        widths = [page.text_width(t, size) for t, _ in items]
        mk = msize / S.PT_PER_MM
        sp = 2.4
        title = spec.get('legend_title')
        tw0 = page.text_width(title, size) + 1.6 if title else 0.0
        total = tw0 + sum(widths) + len(items) * (mk + 0.9) + sp * (len(items) - 1)
        cx = x + w / 2 - total / 2
        yb = y - spec.get('legend_gap', 1.2)
        if title:
            page.text(cx, yb, title, gid=f'{gid}-legend')
            cx += tw0
        for (t, c), tw in zip(items, widths):
            ax2 = page.fig
            from matplotlib.lines import Line2D
            ax2.add_artist(Line2D([page.fx(cx + mk / 2)], [page.fy(yb - 0.36 * size / S.PT_PER_MM)],
                                  transform=page.fig.transFigure, marker='o', markersize=msize * 1.25,
                                  markerfacecolor='none', markeredgecolor=c, markeredgewidth=ln['marker_edge_pt'] * 1.2,
                                  linestyle='none', gid=f'{gid}-legend'))
            page.text(cx + mk + 0.9, yb, t, gid=f'{gid}-legend')
            cx += mk + 0.9 + tw + sp
    # ---- statistics labels, exactly as given
    drawn = []
    if ann is not None and len(ann):
        key = {(b['group'], b['series']): b for b in bars}
        probs = []
        for i, r in ann.iterrows():
            if (r['group'], r['series']) not in key:
                probs.append(f"{ann.attrs.get('path', 'stats')} line {i + 2}: no bar for group '{r['group']}', "
                             f"series '{r['series']}'.")
            g2 = r.get('group2', '') or ''
            if g2 and (g2, r.get('series2', '')) not in key:
                probs.append(f"{ann.attrs.get('path', 'stats')} line {i + 2}: no bar for group2 '{g2}', series2 "
                             f"'{r.get('series2', '')}'.")
        if spec.get('require_labels'):
            have = {(r['group'], r['series']) for _, r in ann.iterrows() if not (r.get('group2', '') or '')}
            miss = [f'{g}/{s}' for (g, s) in key if (g, s) not in have]
            if miss:
                probs.append(f"{ann.attrs.get('path', 'stats')}: no label for bar(s) {miss}. Every bar needs a label "
                             f"(or set require_labels: false in the config).")
        if probs:
            raise ValidationError(probs)

        def top_mm(b):
            return Y(max(b['values'].max(), b['mean'] + b['err']))
        off = spec.get('label_offset_mm', 0.8)
        # statistics labels (letters / asterisks): own size, default style font.stats_pt (7.0 pt = the journal maximum)
        ssz = spec.get('stats_pt', page.st['font'].get('stats_pt', page.st['font']['size_pt']))
        srot = spec.get('star_rotation', 90)
        lrot = spec.get('letter_rotation', 0)          # e.g. 90: letter strings read bottom to top over narrow bars
        lrot_min = spec.get('letter_rotation_min_len', 1)   # only strings longer than this are turned
        if 'group2' not in ann:
            ann = ann.assign(group2='', series2='')
        single = ann[ann['group2'] == '']
        for _, r in single.iterrows():
            b = key[(r['group'], r['series'])]
            lab = r['label']
            if r['kind'] == 'stars' and srot:
                page.text(X(b['x']), top_mm(b) - off, lab, rotation=srot, ha='left', va='center', gid=f'{gid}-stats', size=ssz)
            elif r['kind'] == 'letters' and lrot and len(lab) > lrot_min:
                page.text(X(b['x']), top_mm(b) - off, lab, rotation=lrot, ha='left', va='center', gid=f'{gid}-stats', size=ssz)
            else:
                page.text(X(b['x']), top_mm(b) - off, lab, ha='center', gid=f'{gid}-stats', size=ssz)
            drawn.append(dict(group=r['group'], series=r['series'], label=lab))
        # brackets between two bars (statistics rows with group2/series2)
        br = ann[ann['group2'] != '']
        spans = []
        for _, r in br.iterrows():
            b1, b2 = key[(r['group'], r['series'])], key[(r['group2'], r['series2'])]
            xa, xb = sorted((b1['x'], b2['x']))
            spans.append((xa, xb, r))
        # brackets that touch or overlap share one height (Prism draws touching ones as one line)
        spans.sort(key=lambda t: (t[0], t[1]))
        clusters = []
        for sp_ in spans:
            if clusters and sp_[0] <= clusters[-1]['xb'] + 1e-9:
                clusters[-1]['items'].append(sp_); clusters[-1]['xb'] = max(clusters[-1]['xb'], sp_[1])
            else:
                clusters.append(dict(xa=sp_[0], xb=sp_[1], items=[sp_]))
        for c in clusters:
            inside = [b for b in bars if c['xa'] - 1e-9 <= b['x'] <= c['xb'] + 1e-9]
            yline = min(top_mm(b) for b in inside) - off - 0.6
            levels, extents = [], []
            for xa, xb, r in c['items']:
                lab_w = page.text_width(r['label'], ssz)
                rot_lab = r['kind'] == 'stars' and srot and lab_w > X(xb) - X(xa)
                cxm = X((xa + xb) / 2)
                ext = (cxm - 0.4, cxm + 0.4) if rot_lab else (cxm - lab_w / 2 - 0.25, cxm + lab_w / 2 + 0.25)
                lvl = 0
                while (any(not (xb <= a or xa >= bb) and L == lvl for a, bb, L in levels) or
                       any(L == lvl and not (ext[1] <= e0 or ext[0] >= e1) for e0, e1, L in extents)):
                    lvl += 1
                levels.append((xa, xb, lvl))
                extents.append((ext[0], ext[1], lvl))
                yy = yline - lvl * 2.6
                page.bracket(X(xa), X(xb), yy, tick=0.6, gid=f'{gid}-stats')
                lab = r['label']
                if r['kind'] == 'stars' and srot and page.text_width(lab, ssz) > X(xb) - X(xa):
                    page.text(X((xa + xb) / 2), yy - 0.4, lab, rotation=srot, ha='left', va='center',
                              gid=f'{gid}-stats', size=ssz)
                else:
                    page.text(X((xa + xb) / 2), yy - 0.45, lab, ha='center', gid=f'{gid}-stats', size=ssz)
                drawn.append(dict(group=r['group'], series=r['series'], group2=r['group2'],
                                  series2=r['series2'], label=lab))
    if rot:
        import math as _m
        longest = max(page.text_width(S.display(g, kind)) for g in groups)
        xl_bottom = y + h + tick_mm + 0.6 + longest * _m.sin(_m.radians(rot)) + size / S.PT_PER_MM * _m.cos(_m.radians(rot))
    else:
        xl_bottom = y + h + tick_mm + 0.5 + size * 0.95 / S.PT_PER_MM
    rec = dict(type='dotbar', box_mm=[round(v, 3) for v in box], y_range=[ymin, ymax], y_step=ystep,
               x_label_bottom_mm=round(xl_bottom, 2),
               error=err_kind, labels_drawn=drawn,
               bars=[dict(group=b['group'], series=b['series'], n=b['n'], mean=round(b['mean'], 6),
                          sd=round(b['sd'], 6), err=round(b['err'], 6), x_mm=round(X(b['x']), 3),
                          top_mm=round(Y(b['mean']), 3)) for b in bars])
    return ax, rec
