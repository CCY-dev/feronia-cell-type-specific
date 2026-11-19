"""Figure 1: UMAPs of WT and fer-4 roots (a), genes by number of subpopulations in which they are DEGs (b), DEGs per
subpopulation (c) and DEGs in hormone-response gene sets (d).

Reads config/fig1.yaml and the tables in data/fig1/ written by tools/fig1_tables.py from Chia-Yu Chen's pinned
repository (commit d9d516e; upstream/chia-yu_d9d516e/). Nothing is computed from raw data here: counts, P values and
colours are drawn exactly as tabled; the colour of d is her colorRamp2 applied to the tabled statistic.
Panel d shows her as-coded P values unless the config sets panels.d.statistic: q_proposed (opt-in proposal).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ferfigs import style as S
from ferfigs.layout import ASCENDER
from ferfigs.schemas import ValidationError
from ferfigs.panels import fig1_umap, fig1_charts

SCHEMA = {
    'umap': ['umap_1', 'umap_2', 'genotype', 'celltype.anno.Li', 'time.anno.Li', 'FER_SCT'],
    'palette': ['map', 'name', 'hex', 'level', 'draw_order'],
    'b': ['n_subpopulations', 'n_genes'],
    'c': ['rank', 'subpopulation', 'up', 'down', 'total'],
    'd': ['subpopulation', 'hormone', 'overlap', 'n_deg', 'n_set', 'genome_size', 'p_value', 'mlog10p'],
    'd_q_proposed': ['subpopulation', 'hormone', 'overlap', 'n_universe', 'n_set_in_universe', 'p_proposed',
                     'q_proposed', 'mlog10q'],
    'd_q_authors': ['subpopulation', 'hormone', 'overlap', 'n_universe', 'n_set_in_universe', 'p', 'q'],
    'd_colour_scale': ['value', 'colour', 'hex'],
    'summary': ['key', 'value'],
}
NUMERIC = {'umap': ['umap_1', 'umap_2', 'FER_SCT'], 'b': ['n_genes'], 'c': ['rank', 'up', 'down', 'total'],
           'd': ['overlap', 'n_deg', 'n_set', 'genome_size', 'p_value', 'mlog10p'],
           'd_q_proposed': ['overlap', 'p_proposed', 'q_proposed', 'mlog10q'], 'd_colour_scale': ['value'],
           'd_q_authors': ['overlap', 'n_universe', 'n_set_in_universe', 'p', 'q']}
OPTIONAL = {'d_q_authors'}          # read only when panels.d.statistic asks for it


def read_table(ctx, key):
    return ctx.custom_table(ctx.cfg['tables'][key], SCHEMA[key], NUMERIC.get(key, []),
                            hint='Run  python tools/fig1_tables.py  (after python tools/fetch_upstream.py) or fix the '
                                 'path under tables: in config/fig1.yaml.')


def validate(ctx):
    """Called by  python build.py validate : checks every table listed under tables: in config/fig1.yaml."""
    keys = [k for k in ctx.cfg['tables'] if k in SCHEMA and not (k in OPTIONAL and not ctx.path(ctx.cfg['tables'][k]).exists())]
    for key in keys:
        read_table(ctx, key)
    return len(keys)


def build(ctx):
    cfg = ctx.cfg
    P, G, N = cfg['panels'], cfg['layout'], cfg.get('names', {})
    page = ctx.page()
    st = page.st

    def name(kind, x):
        """Printed name: config/labels.yaml (kit-wide) first, then names: in config/fig1.yaml, then the data name."""
        lab, own = (S.labels().get(kind) or {}).get(x), (N.get(kind) or {}).get(x)
        if lab is not None and own is not None and lab != own:
            ctx.warnings.append(f"Fig. 1: '{x}' is printed as '{lab}' (config/labels.yaml); config/fig1.yaml says "
                                f"'{own}'. Remove one of the two.")
        return lab if lab is not None else (own if own is not None else x)

    umap = read_table(ctx, 'umap')
    pal = read_table(ctx, 'palette')
    # opt-in colour changes (empty by default = her notebook palettes): panels.a.palette_override.{celltype,stage}
    over = (P['a'].get('palette_override') or {})
    for kind, m in over.items():
        for nm, hexc in (m or {}).items():
            sel = (pal['map'] == kind) & (pal['name'] == nm)
            if not sel.any():
                raise ValidationError([f"config/fig1.yaml palette_override.{kind}: '{nm}' is not a {kind} name in "
                                       f"fig1a_palette.csv."])
            pal.loc[sel, 'hex'] = hexc
            ctx.warnings.append(f'Fig. 1a: {kind} colour of {nm} overridden in the config ({hexc}); the upstream '
                                f'palette is in data/fig1/fig1a_palette.csv.')
    tb, tc, td = read_table(ctx, 'b'), read_table(ctx, 'c'), read_table(ctx, 'd')
    scale_t = read_table(ctx, 'd_colour_scale')
    summary = read_table(ctx, 'summary').set_index('key')['value']
    T = cfg['tables']

    # ------------------------------------------------------------ a: six maps, header, row labels, FER key, glyph
    lines = {**st['lines'], **(G.get('lines') or {})}          # figure overrides (approved page) of the house widths

    def letter(lt, x, baseline):
        """8-pt letter whose baseline sits on the approved letter's baseline (layout gives x and baseline)."""
        page.panel_letter(lt, x, baseline - ASCENDER * st['font']['panel_label_pt'] / S.PT_PER_MM)

    pa = P['a']
    ga = dict(G['a'])
    # checked before drawing, so a mismatched or too-small cells table is refused with a message, not a drawing error
    if int(summary.get('n_cells', len(umap))) != len(umap):
        raise ValidationError([f"fig1_summary.csv n_cells {summary['n_cells']} differs from the UMAP table "
                               f"({len(umap):,} cell{'' if len(umap) == 1 else 's'}). Supply the matching summary with "
                               f"the cells (template sheet Fig1a_summary) or re-run tools/fig1_tables.py."])
    for c in ('umap_1', 'umap_2'):
        if umap[c].nunique() < 2:
            raise ValidationError([f"{T['umap']}: {c} has fewer than two distinct values; the maps need cells spread "
                                   f"in both UMAP coordinates."])
    out_a = fig1_umap.draw_grid(page, umap, pal, pa, ga, dict(genotype=lambda g: name('genotypes', g)), lines)
    letter('a', G['letters_x'], ga['letter_baseline'])
    ctx.record('a', dict(out_a['rec'], source_tables={'per-cell UMAP coordinates, annotations, FER': T['umap'],
                                                      'colours (her notebook)': T['palette']}))
    # ------------------------------------------------------------ key for a (cell type, stage), above b
    K = G['key']
    bottoms = []
    for kind, dx, cols in (('celltype', 0.0, K['celltype_col_dx']), ('stage', K['stage_dx'], [0.0])):
        sub = pal[pal['map'] == kind].copy()
        present = set(umap[pa['columns'][kind]])
        sub = sub[sub['name'].isin(present)].sort_values('level', key=lambda s: s.astype(int))
        entries = [(name('celltypes' if kind == 'celltype' else 'stages', r['name']), r['hex']) for _, r in sub.iterrows()]
        bottoms.append(fig1_charts.legend_key(page, K['x'] + dx, K['y'], pa['key'][kind]['title'], entries, K, cols,
                                              pa['key'][kind]['rows_per_column'], edge=K['edge'].get(kind),
                                              gid=f'fig1-a-key-{kind}'))
    # ------------------------------------------------------------ b
    pb, gb = dict(P['b']), G['b']
    bins = list(tb['n_subpopulations'])
    pb['_bins'] = bins
    labels = [pb['bin_labels'].get(x, x) for x in bins]
    vals = [int(v) for v in tb['n_genes']]
    total = int(sum(vals))
    if 'n_deg_unique' in summary and int(summary['n_deg_unique']) != total:
        raise ValidationError([f"fig1b counts sum to {total:,} but fig1_summary.csv n_deg_unique is "
                               f"{int(summary['n_deg_unique']):,}; re-run tools/fig1_tables.py."])
    # plot frame: x axis a fixed distance above the bottom of panel a (approved page), top at y_max
    ya = out_a['bottom'] - gb['axis_above_a_bottom']
    out_b = fig1_charts.count_bars(page, (gb['x0'], gb['y_top'], gb['x1'], ya), labels, vals, pb, gb, lines)
    n12 = int(sum(v for b_, v in zip(bins, vals) if b_ in pb['highlight_bins']))
    pct = int(round(100 * n12 / total))
    # summary centred in the empty region: between the last highlighted bar and the end of the axis, and between
    # the top of the first other bar and y_max
    i3 = len(pb['highlight_bins'])
    cx = (out_b['xc'][i3 - 1] + out_b['bar_w'] / 2 + gb['x1']) / 2
    cy = (out_b['tops'][i3] + gb['y_top']) / 2
    (d1, d2), (s1, s2) = gb['summary_dy'], gb['summary_pt']
    page.text(cx, cy + d1, pb['total_text'].format(total=f'{total:,}'), size=s1, ha='center', weight='bold',
              gid='fig1-b-summary')
    page.text(cx, cy + d2, pb['share_text'].format(n12=f'{n12:,}', pct=pct), size=s2, ha='center',
              gid='fig1-b-summary')
    letter('b', *gb['letter'])
    ctx.record('b', dict(type='count_bars', bins=bins, values=vals, total=total, one_or_two=n12, percent=pct,
                         y_range=[0, pb['y_max']], source_tables={'genes by number of subpopulations': T['b']}))
    # ------------------------------------------------------------ c + d (shared rows)
    pc, pd_ = P['c'], dict(P['d'])
    gcd = dict(G['cd'])
    gcd['y_top'] = out_a['bottom'] + gcd['gap_after_a']
    tc = tc.sort_values('rank').reset_index(drop=True)
    rows = list(tc['subpopulation'])
    unknown = sorted(set(td['subpopulation']) - set(rows))
    if unknown:
        raise ValidationError([f'fig1d: subpopulations {unknown} are not in fig1c; both tables must come from the same '
                               f'run of tools/fig1_tables.py.'])
    hs = list(pd_['hormones'])
    stat = pd_.get('statistic', 'p_as_coded')
    dd = td.set_index(['subpopulation', 'hormone'])
    missing = [(s_, h) for s_ in rows for h in hs if (s_, h) not in dd.index]
    if missing:
        raise ValidationError([f'fig1d_hormone_overlap.csv: no row for {missing[:3]} (subpopulation, hormone).'])
    counts = np.array([[int(dd.loc[(s_, h), 'overlap']) for h in hs] for s_ in rows])
    if stat == 'p_as_coded':
        values = np.array([[float(dd.loc[(s_, h), 'mlog10p']) for h in hs] for s_ in rows])
        stat_table = T['d']
    elif stat == 'q_proposed':
        tq = read_table(ctx, 'd_q_proposed').set_index(['subpopulation', 'hormone'])
        if not all(int(tq.loc[(s_, h), 'overlap']) == counts[i, j] for i, s_ in enumerate(rows) for j, h in enumerate(hs)):
            raise ValidationError(['fig1d q table: overlap counts differ from the as-coded table.'])
        values = np.array([[float(tq.loc[(s_, h), 'mlog10q']) for h in hs] for s_ in rows])
        stat_table = T['d_q_proposed']
        ctx.warnings.append('Fig. 1d is coloured by the PROPOSED q statistic (panels.d.statistic: q_proposed), not by '
                            'the as-coded P of the upstream notebook. Agree with Chia-Yu before using it.')
    elif stat == 'q_authors':
        if not ctx.path(T['d_q_authors']).exists():
            raise ValidationError([f"panels.d.statistic is q_authors but {T['d_q_authors']} does not exist: import the "
                                   f"Fig1d_overlap_q template sheet first."])
        tq = read_table(ctx, 'd_q_authors')
        dup = tq.duplicated(['subpopulation', 'hormone'])
        if dup.any():
            raise ValidationError([f"{T['d_q_authors']}: (subpopulation, hormone) repeated on line(s) "
                                   f"{[int(i) + 2 for i in tq.index[dup]][:5]}."])
        tq = tq.set_index(['subpopulation', 'hormone'])
        miss = [(s_, h) for s_ in rows for h in hs if (s_, h) not in tq.index]
        if miss:
            raise ValidationError([f"{T['d_q_authors']}: no row for {miss[:3]} (all {len(rows) * len(hs)} cells needed)."])
        if not all(int(tq.loc[(s_, h), 'overlap']) == counts[i, j] for i, s_ in enumerate(rows) for j, h in enumerate(hs)):
            raise ValidationError([f"{T['d_q_authors']}: overlap counts differ from the as-coded table; the colour "
                                   f"statistic may change, the printed counts may not."])
        qv = np.array([[float(tq.loc[(s_, h), 'q']) for h in hs] for s_ in rows])
        if not ((qv >= 0) & (qv <= 1)).all():
            raise ValidationError([f"{T['d_q_authors']}: q must be between 0 and 1."])
        values = -np.log10(qv + 1e-200)            # same transform as her notebook applies to P
        stat_table = T['d_q_authors']
    else:
        raise ValidationError([f"config/fig1.yaml panels.d.statistic = '{stat}': use p_as_coded, q_proposed or "
                               f"q_authors."])
    pd_['key_title'] = pd_['key_titles'][stat]
    pd_['hormone_labels'] = pd_.get('hormone_labels') or {}
    sc = scale_t.sort_values('value')
    scale = (list(sc['value'].astype(float)), list(sc['colour']))
    c_rows = [(name('subpopulations', s_), int(r.up), int(r.down), int(r.total)) for s_, r in zip(rows, tc.itertuples())]
    out_cd = fig1_charts.rows_bars_heatmap(page, gcd, c_rows, counts, values, scale, pc, pd_, lines)
    letter('c', G['letters_x'], gcd['y_top'] + gcd['letter_baseline'])
    letter('d', gcd['d_letter_x'], gcd['y_top'] + gcd['letter_baseline'])
    ctx.record('c', dict(type='stacked_rows', rows=rows, up=[int(v) for v in tc['up']], down=[int(v) for v in tc['down']],
                         total=[int(v) for v in tc['total']], x_range=[0, pc['x_max']],
                         source_tables={'DEGs per subpopulation': T['c']}))
    ctx.record('d', dict(type='overlap_heatmap', statistic=stat, rows=rows, hormones=hs,
                         counts=counts.tolist(), fills=out_cd['drawn']['fills'],
                         source_tables={'overlap counts and P (as coded)': T['d'], 'colour scale': T['d_colour_scale'],
                                        **({'proposed q (opt-in)': T['d_q_proposed']} if stat == 'q_proposed' else {}),
                                        **({"authors' adjusted q": T['d_q_authors']} if stat == 'q_authors' else {})},
                         colour_statistic_table=stat_table))
    bottom = out_cd['bottom']
    if bottom > cfg['size_mm'][1] + 1e-6:
        raise ValidationError([f'Fig. 1 content ends at {bottom:.1f} mm, below the page height {cfg["size_mm"][1]} mm: '
                               f'reduce layout.cd.row_pitch or layout.a.map_w.'])
    ctx.note('layout_mm', dict(a_bottom=round(out_a['bottom'], 3), b_axis=round(ya, 3), cd_top=round(gcd['y_top'], 3),
                               rows_top=round(out_cd['y_rows'], 3), rows_bottom=round(bottom, 3)))
    return page
