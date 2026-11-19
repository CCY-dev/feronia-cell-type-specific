"""Figure 6: RALF23 induces YUC5 expression through DOF TFs.

  a  root length, dof mutants        dotbar
  b  seedlings + printed ratios      ferfigs.panels.fig6_seedlings
  c  pYUC5:GFP-GUS in dof mutants    ferfigs.panels.imagerow + LUT key
  d  quantification of c             dotbar
  e  DOF qPCR, one chart per gene    dotbar x 4 with one shared series key
  f  pDAG1 / pDOF2.2 reporters       imagerow + LUT key
  g  quantification of f             dotbar x 2 (bracket + stars) with one shared series key
  h  root image + working model      ferfigs.panels.fig6_model
Data: data/fig6/ (measurement tables from the Source Data by tools/fig6_tables.py). Every number on a chart is
computed from a measurement table at build time; labels are drawn exactly as given in the *_stats.csv tables.
"""
import copy

import pandas as pd
from matplotlib.lines import Line2D

from ferfigs import style as S
from ferfigs.panels import dotbar, fig6_model, fig6_seedlings, imagerow
from ferfigs.panels.fig6_common import MM


def _runs(rows, column):
    out = []
    for _, r in rows.sort_values('col').iterrows():
        lab = S.display(str(r[column]))
        if out and out[-1][0] == lab:
            out[-1][1].append(r['col'])
        else:
            out.append((lab, [r['col']]))
    return out


def _key(page, fig_dir, k, x, y_tiles_bottom, gid):
    """Colour key as in kit Fig. 4: LUT strip in a thin frame, 'Max' above, 'Min' below level with the tile bottoms."""
    y1 = y_tiles_bottom - k['bottom_offset']
    y0 = y1 - k['h']
    kp = page.image(fig_dir / k['file'], x, y0, h=k['h'], gid=gid, name=k['file'])
    page.rect(x, y0, kp.w_mm, k['h'], lw=0.25, edge='k', gid=gid)
    hi, lo = k.get('labels', ('Max', 'Min'))
    cx = x + kp.w_mm / 2
    page.text(cx, y0 - k['max_gap'], hi, ha='center', gid=gid)
    page.text(cx, y_tiles_bottom, lo, ha='center', gid=gid)
    return x + kp.w_mm


def _series_key(page, cx_mid, yb, series, gid):
    """The dotbar series key (same marker, size and spacing as ferfigs.panels.dotbar), centred at cx_mid."""
    st = page.st
    size = st['font']['size_pt']
    msize = st['markers']['diameter_pt']
    mk = msize / S.PT_PER_MM
    items = [(S.display(s, 'series'), S.series_colour(s)) for s in series]
    widths = [page.text_width(t, size) for t, _ in items]
    sp = 2.4
    total = sum(widths) + len(items) * (mk + 0.9) + sp * (len(items) - 1)
    cx = cx_mid - total / 2
    for (t, c), tw in zip(items, widths):
        page.fig.add_artist(Line2D([page.fx(cx + mk / 2)], [page.fy(yb - 0.36 * size / S.PT_PER_MM)],
                                   transform=page.fig.transFigure, marker='o', markersize=msize * 1.25,
                                   markerfacecolor='none', markeredgecolor=c,
                                   markeredgewidth=st['lines']['marker_edge_pt'] * 1.2, linestyle='none', gid=gid))
        page.text(cx + mk + 0.9, yb, t, gid=gid)
        cx += mk + 0.9 + tw + sp


def _chart(ctx, page, name, box, fs):
    s = ctx.cfg['panels'][name]
    df = ctx.table(s['data'], 'measurements')
    st = ctx.table(s['stats'], 'annotations') if s.get('stats') else None
    ax, rec = dotbar.draw(page, tuple(box), df, s, st, gid=f'fig6-{name}')
    ax.tick_params(axis='y', labelsize=fs)
    if s.get('title'):
        page.text(box[0] + box[2] / 2, box[1] - ctx.cfg['panels'][name.split('_')[0]].get('title_gap', 1.6),
                  s['title'], ha='center', gid=f'fig6-{name}-title')
    ctx.record(name, rec)
    return rec


def _image_row(ctx, page, fig_dir, man, p, mm_per_px=None, reporter=None):
    rc, im = ctx.cfg['panels'][p], ctx.cfg['images']
    rows = man[man['panel'] == p].sort_values('col')
    first = rows.iloc[0]
    if mm_per_px is None:
        mm_per_px = rc['tile_h'] / imagerow.tile_px(fig_dir, first)[1]
    out = imagerow.draw(page, fig_dir, rows, rc['x0'], rc['y'], mm_per_px=mm_per_px,
                        groups=_runs(rows, 'genotype'), within=rc['within'], between=rc['between'],
                        reporter=reporter, pm={r['col']: im['signs'][r['treatment']] for _, r in rows.iterrows()},
                        pm_label=rc['treatment_label'], y_reporter=rc.get('y_reporter'), y_groups=rc['y_groups'],
                        y_pm=rc['y_pm'], gid=f'fig6-{p}')
    xk = _key(page, fig_dir, im['key'], out['x1'] + im['key']['gap'], out['y_bottom'], gid=f'fig6-{p}-key')
    upp = float(first['um_per_px'])
    ctx.record(p, dict(type='images', tiles=list(rows['file']), mm_per_px=round(mm_per_px, 6),
                       um_per_mm=round(upp / mm_per_px, 3), right_mm=round(xk, 2),
                       bottom_mm=round(out['y_bottom'], 2)))
    return mm_per_px, out


def build(ctx):
    cfg = ctx.cfg
    P = cfg['panels']
    page = ctx.page()
    page.st = copy.deepcopy(page.st)
    page.st['font']['size_pt'] = float(cfg.get('font_pt', page.st['font']['size_pt']))
    fs = page.st['font']['size_pt']
    mm = MM(page)
    fig_dir = ctx.path('data/fig6')
    man = ctx.table(cfg['manifest'], 'image_manifest')

    # ---- a: root length (authors' per-seedling data)
    rec_a = _chart(ctx, page, 'a', P['a']['box'], fs)

    # ---- b: seedlings with printed ratios; printed vs computed ratio recorded (source rows)
    pb = P['b']
    tab_b = pd.read_csv(ctx.path(pb['table']), dtype=str, keep_default_na=False)
    rec = fig6_seedlings.draw(page, fig_dir, man[man['panel'] == 'b'], tab_b, pb, pb['x_tiles'], pb['y_tiles'],
                              pb['tile_w'], pb['gap'], gid='fig6-b')
    bars = {(b['group'], b['series']): b['mean'] for b in rec_a['bars']}
    rows = []
    for _, r in tab_b.iterrows():
        if not r['ratio']:          # public copy without the printed ratios (data/LABELS_NOT_INCLUDED)
            continue
        comp = bars[(r['genotype'], 'RALF23')] / bars[(r['genotype'], 'Mock')]
        rows.append([r['genotype'], r['ratio'], round(comp, 4), f'{comp:.2f}' == r['ratio']])
        if f'{comp:.2f}' != r['ratio']:
            ctx.warnings.append(f"fig6 b: printed ratio {r['ratio']} for {r['genotype']} differs from "
                                f"mean(RALF23)/mean(Mock) = {comp:.4f} of panel a (printed value drawn)")
    rec['source_tables'] = {'printed ratios and -/+ positions': pb['table']}
    rec['source_rows'] = {'printed vs computed ratio': {
        'columns': ['genotype', 'printed', 'mean(RALF23)/mean(Mock) of panel a', 'printed = computed (2 dp)'],
        'rows': rows}}
    ctx.record('b', rec)

    # ---- c: pYUC5:GFP-GUS in dof mutants; d: quantification
    mpp, out_c = _image_row(ctx, page, fig_dir, man, 'c', reporter=P['c']['reporter'])
    _chart(ctx, page, 'd', P['d']['box'], fs)

    # ---- e: qPCR, one chart per gene, one shared key
    pe = P['e']
    for name in pe['charts']:
        _chart(ctx, page, name, (pe['x'][name], pe['y'], pe['w'], pe['h']), fs)
    xs = [pe['x'][n] for n in pe['charts']]
    _series_key(page, (xs[0] + xs[-1] + pe['w']) / 2, pe['key_y'], P[pe['charts'][0]]['series'], gid='fig6-e-key')

    # ---- f: DOF reporters at the magnification of c (same 1.887 um/px); g: quantification
    _image_row(ctx, page, fig_dir, man, 'f', mm_per_px=mpp)
    pg = P['g']
    for name in pg['charts']:
        _chart(ctx, page, name, (pg['x'][name], pg['y'], pg['w'], pg['h']), fs)
    xs = [pg['x'][n] for n in pg['charts']]
    _series_key(page, (xs[0] + xs[-1] + pg['w']) / 2, pg['key_y'], P[pg['charts'][0]]['series'], gid='fig6-g-key')

    # ---- h: root image + model
    ph = dict(P['h'])
    row_h = man[man['panel'] == 'h'].iloc[0]
    ctx.record('h', fig6_model.draw(page, mm, fig_dir, row_h, ph, ph['box']))

    for letter, (x, y_top) in cfg['letters'].items():
        page.panel_letter(letter, x, y_top)
    return page


def validate(ctx):
    ctx.custom_table(ctx.cfg['panels']['b']['table'], ['tile', 'genotype', 'row', 'col', 'ratio', 'minus_x_px',
                                                       'plus_x_px'], ['row', 'col', 'minus_x_px', 'plus_x_px'],
                     hint='Regenerate with ROOT/scripts/26_fig6_tables.py.')
    return 1
