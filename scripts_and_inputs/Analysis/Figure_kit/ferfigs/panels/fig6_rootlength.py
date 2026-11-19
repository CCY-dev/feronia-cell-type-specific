"""Fig. 6g until per-seedling data arrive: a clearly marked placeholder in the exact frame the dot-bar chart will
use (same axes, ticks, labels and key as ferfigs.panels.dotbar draws for this config, but no bars, points or
statistics), inside a dashed outline with a note naming the file that switches the panel to the real chart.
When data/fig6/fig6g_root_length.csv exists, figures/fig6.py draws the chart with dotbar.draw instead.

mockup(): optional layout check before the data arrive (config g.mockup). Bars, SD and the
printed labels come from a table traced from the original chart (tools/fig6g_trace_mockup.py); the points are
simulated at normal quantiles so each bar's mean and SD equal the traced values exactly. The plot carries a visible
mock-up note, Source Data shows only the traced table under the status 'mockup', and the build warns.
"""
from __future__ import annotations

import numpy as np

from .. import style as S
from .dotbar import nice_axis
from .fig6_common import MM, display


def placeholder(page, mm: MM, box, outline, spec, note, gid='fig6-g'):
    st = page.st
    fs = st['font']['size_pt']
    x, y, w, h = box
    groups, series = list(spec['groups']), list(spec['series'])
    n_s, gap = len(series), spec.get('group_gap', 0.7)
    xlo = -0.5 - 0.35
    xhi = (len(groups) - 1) * (n_s + gap) + n_s - 1 + 0.5 + 0.35
    ymax, ystep = float(spec['y_max']), float(spec.get('y_step') or nice_axis(float(spec['y_max']))[1])
    ax = page.axes(x, y, w, h, gid=f'{gid}-axes')
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(0, ymax)
    ticks = np.arange(0, ymax + ystep * 0.001, ystep)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{t:g}' for t in ticks])
    centres = [gi * (n_s + gap) + (n_s - 1) / 2 for gi in range(len(groups))]
    ax.set_xticks(centres)
    ax.set_xticklabels([])
    tick_mm = st['lines']['tick_len_pt'] / S.PT_PER_MM

    def X(xd):
        return x + (xd - xlo) / (xhi - xlo) * w
    for gi, g in enumerate(groups):
        page.text(X(centres[gi]), y + h + tick_mm + 0.5, display(g, spec.get('display_labels'),
                                                               spec.get('label_kind', 'genotypes')),
                  ha='center', va='top', gid=f'{gid}-xlabels')
    tl_w = max(page.text_width(f'{t:g}', fs) for t in ticks)
    lines = spec['y_label'] if isinstance(spec['y_label'], list) else [spec['y_label']]
    xl = x - tick_mm - 0.5 - tl_w - 0.9
    for i, line in enumerate(lines[::-1]):
        page.text(xl - i * fs * 1.18 / S.PT_PER_MM, y + h / 2, line, rotation=90, ha='center', gid=f'{gid}-ylabel')
    # series key above the frame, as dotbar draws it
    msize = st['markers']['diameter_pt']
    items = [(S.display(s, 'series'), S.series_colour(s)) for s in series]
    widths = [page.text_width(t, fs) for t, _ in items]
    mk = msize / S.PT_PER_MM
    total = sum(widths) + len(items) * (mk + 0.9) + 2.4 * (len(items) - 1)
    cx = x + w / 2 - total / 2
    yb = y - spec.get('legend_gap', 1.2)
    for (t, c), tw in zip(items, widths):
        mm.circle(cx + mk / 2, yb - 0.36 * fs / S.PT_PER_MM, mk * 1.25, face='#ffffff', edge=c,
                  lw=st['lines']['marker_edge_pt'] * 1.2, gid=f'{gid}-legend', z=4)
        page.text(cx + mk + 0.9, yb, t, gid=f'{gid}-legend')
        cx += mk + 0.9 + tw + 2.4
    # the placeholder marking: dashed outline and a centred note in the empty plot area
    ox, oy, ow, oh = outline
    mm.rect(ox, oy, ow, oh, face='none', edge=note['outline_colour'], lw=note['outline_pt'], gid=f'{gid}-placeholder',
            z=0).set_linestyle((0, tuple(note['dashes'])))
    yy = y + h / 2 - (len(note['lines']) - 1) * note['line_mm'] / 2
    widest = max(page.text_width(t, note.get('size_pt', fs)) for t in note['lines'])
    if widest > w - 1.0:
        raise ValueError(f'g placeholder note ({widest:.1f} mm) is wider than the plot frame ({w:.1f} mm)')
    for i, t in enumerate(note['lines']):
        page.text(x + w / 2, yy + i * note['line_mm'], t, ha='center', size=note.get('size_pt', fs),
                  weight='bold' if i == 0 else None, gid=f'{gid}-placeholder')
    return dict(type='placeholder', box_mm=[round(v, 3) for v in box], groups=groups, series=series,
                note=note['lines'])


def mockup(page, ctx, spec, gid='fig6-g'):
    """Mock-up chart in the final frame (see module docstring). Returns the panel record."""
    from statistics import NormalDist

    import pandas as pd

    from .dotbar import draw
    m = spec['mockup']
    t = pd.read_csv(ctx.path(m['table']), dtype=str, keep_default_na=False)
    rows = []
    for _, r in t.iterrows():
        n, mu, sd = int(r['n']), float(r['mean']), float(r['sd'])
        z = np.array([NormalDist().inv_cdf((i + 0.5) / n) for i in range(n)])
        z = (z - z.mean()) / z.std(ddof=1)                     # exact mean 0, SD 1
        for i, v in enumerate(mu + sd * z):
            rows.append(dict(group=r['group'], series=r['series'], replicate=str(i + 1), value=float(v),
                             unit=spec.get('unit', 'mm')))
    df = pd.DataFrame(rows)
    ann = t[t['label'] != ''][['group', 'series', 'label', 'kind']].copy()
    ann['group2'], ann['series2'] = '', ''
    ann = ann.reset_index(drop=True)
    ann.attrs['path'] = m['table']
    ax, rec = draw(page, tuple(spec['plot_box']), df, dict(spec, require_labels=False), ann if len(ann) else None,
                   gid=gid)
    x, y, w, h = spec['plot_box']
    fs = m.get('size_pt', 5.5)
    for i, line in enumerate(m['note']):
        page.text(x + w - 0.6, y + m.get('first_line_mm', 2.4) + i * m.get('line_mm', 2.2), line, size=fs, ha='right',
                  colour=m.get('colour'), fontstyle='italic', gid=f'{gid}-mockup')   # black: no coloured text
    rec.pop('bars', None)                                      # the simulated points are not source data
    rec.update(type='mockup', note=list(m['note']),
               source_tables={'MOCK-UP values traced from the original chart (not data)': m['table']})
    return rec
