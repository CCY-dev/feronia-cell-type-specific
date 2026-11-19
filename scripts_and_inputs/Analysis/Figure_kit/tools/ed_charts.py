"""Extended Data charts (ED Figs 2, 4, 5, 6, 7, 9, 10) re-plotted from the paper's Source Data in the kit chart style
(ferfigs.panels.dotbar + config/style.yaml + config/labels.yaml), with the statistical test named in each legend
recomputed from the source data.

usage (from the kit folder):
  python tools/source_data_to_long.py "<Source data.xlsx>"      # -> source_data/long/
  python tools/ed_charts.py

Reads   source_data/long/ALL_long.csv (+ summary.csv: n and mean per group x series are asserted)
        data/ed/ED_chart_sizes.csv     (page size of each chart, mm)
        data/ed/ED_author_labels.csv   (the statistics labels printed on the published charts; optional)
Writes  output/ed_charts/<ID>.{pdf,png,svg,build.json}   (one small page per chart)
        output/ed_charts/data/<ID>.csv                   (every plotted value, with its source cell)
        output/stats_checks/STATS_CHECK_ED.csv            (label printed on the chart vs label computed with the test
                                                           named in the legend)
        output/stats_checks/NOTES_ED.txt                  (family / experimental-unit checks)

The printed statistics labels are the authors' (from their analysis software). When data/ed/ED_author_labels.csv is
not present, every chart is drawn without statistics labels (same layout otherwise) and STATS_CHECK_ED.csv still
lists the computed labels, with agree = 'labels not in this repository'.

Nothing in the kit is edited: series colours / display names that the kit does not know (RALF1, _L2 lines,
yuc mutants, ...) are added to the IN-MEMORY copies of style.yaml / labels.yaml only.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as SS

KIT = Path(__file__).resolve().parents[1]
LONG_DIR = KIT / 'source_data/long'
OUT = KIT / 'output/ed_charts'
OUT_STATS = KIT / 'output/stats_checks'
LABELS_CSV = KIT / 'data/ed/ED_author_labels.csv'
SIZES_CSV = KIT / 'data/ed/ED_chart_sizes.csv'
sys.path.insert(0, str(KIT))

from ferfigs import style as S  # noqa: E402
from ferfigs.layout import Page  # noqa: E402
from ferfigs.panels import dotbar  # noqa: E402
from ferfigs.stats_helper import cld, stars  # noqa: E402

# ---------------------------------------------------------------- in-memory additions to the kit style / labels
ST = S.load_style()
ST['colours']['series'].update({
    'RALF1': [0, 158, 115], 'RALF31': [204, 121, 167], 'RALF33': [139, 121, 243],
    # genotypes used as a SERIES (qPCR, ED10c): one colour per genotype (stroke of the kit genotype palette;
    # Col-0 lightened so that it differs from fer-4 at 50 % fill)
    'Col-0': [150, 150, 150], 'fer-4': [40, 40, 40], 'FER/fer-4': [120, 120, 120],
    'pWER-KO_L1': [213, 94, 0], 'pCO2-KO_L1': [40, 140, 200], 'pSCR-KO_L1': [200, 135, 0],
    'pSHR-KO_L1': [0, 158, 115], 'pSMB-KO_L1': [175, 160, 0],
})
for name, rgbv in {'Mock': [0, 134, 192], 'RALF1': [0, 158, 115], 'RALF23': [223, 115, 0],
                   'RALF31': [204, 121, 167], 'RALF33': [139, 121, 243]}.items():
    # treatments on the x axis (ED7d, ED7f): colour_by group reads the genotype palette
    ST['colours']['genotype_fill'][name] = rgbv
    ST['colours']['genotype_stroke'][name] = rgbv
for g in ('yuc5', 'yuc2589', 'yuc35789 (yucQ)', 'DOF4.6-OE #17', 'DOF4.6-OE #24'):
    ST['colours']['genotype_fill'].setdefault(g, [255, 255, 255])
    ST['colours']['genotype_stroke'].setdefault(g, [90, 90, 90])
import contextlib


@contextlib.contextmanager
def ed2_colours():
    """ED 2b, 2c, 2e: Col-0 blue / fer-4 orange (two greys are hard to tell apart in these panels).
    Scoped to ED 2 so that the genotype colours of the other charts (e.g. ED 10c) do not change."""
    keep = {k: {g: ST['colours'][k].get(g) for g in ('Col-0', 'fer-4')} for k in ('series', 'genotype_fill', 'genotype_stroke')}
    for k in keep:
        ST['colours'][k]['Col-0'], ST['colours'][k]['fer-4'] = [0, 134, 192], [223, 115, 0]
    try:
        yield
    finally:
        for k, v in keep.items():
            for g, c in v.items():
                if c is None:
                    ST['colours'][k].pop(g, None)
                else:
                    ST['colours'][k][g] = c


LB = S.labels()
LB.setdefault('genotypes', {})
for p in ('pWER-KO', 'pCO2-KO', 'pSCR-KO', 'pSHR-KO', 'pSMB-KO'):
    for L in ('L1', 'L2', 'L3'):
        LB['genotypes'].setdefault(f'{p}_{L}', f'*{p}_{L}*')
LB['genotypes'].update({'yuc5': '*yuc5*', 'yuc2589': '*yuc2589*', 'yuc35789 (yucQ)': '*yucQ*',
                        'DOF4.6-OE #17': '*DOF4.6-OE* #17', 'DOF4.6-OE #24': '*DOF4.6-OE* #24',
                        'Epidermis': 'Epidermis', 'Cortex': 'Cortex', 'Endodermis': 'Endodermis', 'Stele': 'Stele'})
LB.setdefault('series', {})
LB['series'].update({k: k for k in ('RALF1', 'RALF31', 'RALF33')})
LB['series'].update({'Col-0': 'Col-0', 'fer-4': '*fer-4*', 'FER/fer-4': '*FER/fer-4*',
                     **{f'{p}_L1': f'*{p}_L1*' for p in ('pWER-KO', 'pCO2-KO', 'pSCR-KO', 'pSHR-KO', 'pSMB-KO')}})

# ---------------------------------------------------------------- data
LONG = pd.read_csv(LONG_DIR / 'ALL_long.csv', dtype={'replicate': str, 'sample_id': str, 'notes': str},
                   keep_default_na=False)
LONG['value'] = LONG['value'].astype(float)
SUMMARY = pd.read_csv(LONG_DIR / 'summary.csv', keep_default_na=False)


def take(panel, series=None):
    d = LONG[LONG['panel'] == panel].copy()
    if series is not None:
        d = d[d['series'].isin(series if isinstance(series, (list, tuple)) else [series])]
    assert len(d), panel
    return d


def check_summary(d, panel):
    """n and mean per group x series equal source_data/long/summary.csv (all listed values)."""
    sm = SUMMARY[SUMMARY['panel'] == panel]
    for (g, s), sub in d.groupby(['group', 's_orig' if 's_orig' in d else 'series']):
        r = sm[(sm['group'] == g) & (sm['series'] == s)]
        assert len(r) == 1, (panel, g, s)
        assert int(r['n'].iloc[0]) == len(sub), (panel, g, s, len(sub), r['n'].iloc[0])
        assert abs(float(r['mean'].iloc[0]) - sub['value'].mean()) < 1e-4 * max(1, abs(sub['value'].mean())), \
            (panel, g, s)


# ---------------------------------------------------------------- statistics
CHECK = []          # rows of STATS_CHECK_ED.csv


def chk(chart, panel, comparison, label, test, p='', p_alt='', note='', printed='', agree=None, drawn=None):
    """label = computed with the legend's test; drawn = what the chart shows (the authors' printed label;
    '' = nothing printed for this bar). Without the labels file nothing is drawn and agree is not assessed."""
    drawn = label if drawn is None else drawn
    if not HAVE_LABELS:
        drawn, printed, agree = '', '', NO_LABELS
    if agree is None:
        a, b = drawn.replace(' ', ''), label.replace(' ', '')
        if a == '' and b == 'ns':          # ED10c: 'no asterisk' means not significant
            agree = 'yes'
        elif a in ('*', '**', '***', '****') and b in ('*', '**', '***', '****') and chart.startswith('ED2'):
            agree = 'yes (significant; author prints one * per discovery)' if a == '*' else ('yes' if a == b else 'significance yes, level no')
        else:
            agree = 'yes' if a == b else ('significance yes, level no' if a != 'ns' and b != 'ns' and a and b and
                                          not a[0].isalpha() else 'no')
    CHECK.append(dict(chart=chart, panel=panel, comparison=comparison, label_drawn=drawn, label_computed=label,
                      author_printed=printed or drawn, agree=agree, test=test,
                      p_value=(f'{p:.3g}' if isinstance(p, float) else p),
                      p_alternative=(f'{p_alt:.3g}' if isinstance(p_alt, float) else p_alt),
                      source='authors\' printed label (drawn) vs computed from the source data', note=note))


# Statistics labels printed on the published charts, one row per chart x group (x series):
# data/ed/ED_author_labels.csv (label_type 'stats'; label_type 'printed_ratio' = numbers printed beside ED 9d).
# AUTH[chart] is keyed as the drawing code needs it: (group, series) for two-way charts, 'group | series' for the
# Tukey charts over group x series cells, group otherwise.
NO_LABELS = 'labels not in this repository'
HAVE_LABELS = LABELS_CSV.exists()
_TUPLE_KEY = {'ED4f', 'ED5g', 'ED6c', 'ED7b', 'ED10c'}
_CELL_KEY = ('ED9d', 'ED9e_')
AUTH, PRINTED = {}, {}
if HAVE_LABELS:
    for r in csv.DictReader(open(LABELS_CSV, newline='', encoding='utf-8')):
        c, g, s_ = r['chart'], r['group'], r['series']
        if r['label_type'] == 'printed_ratio':
            PRINTED.setdefault(c, {})[g] = r['label']
            continue
        key = (g, s_) if c in _TUPLE_KEY else (f'{g} | {s_}' if c.startswith(_CELL_KEY) else g)
        AUTH.setdefault(c, {})[key] = r['label']


def auth(chart):
    """Printed labels of a chart ({} when the labels file is absent)."""
    return AUTH.get(chart, {}) if HAVE_LABELS else {}


def tukey_letters(d, cellcol, order, alpha=0.05):
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    tk = pairwise_tukeyhsd(d['value'].to_numpy(), d[cellcol].to_numpy(), alpha=alpha)
    res = pd.DataFrame(tk.summary().data[1:], columns=tk.summary().data[0])
    sig = {frozenset((r['group1'], r['group2'])) for _, r in res.iterrows() if bool(r['reject'])}
    means = [d[d[cellcol] == c]['value'].mean() for c in order]
    F, p = SS.f_oneway(*[d[d[cellcol] == c]['value'] for c in order])
    return cld(order, means, sig), res, (F, p)


def pooled(d, cells):
    """Pooled within-cell MSE and its df (= residual of the full-interaction two-way ANOVA)."""
    ss, df_ = 0.0, 0
    for c in cells:
        v = d[d['cell'] == c]['value'].to_numpy()
        ss += ((v - v.mean()) ** 2).sum()
        df_ += len(v) - 1
    return ss / df_, df_


def control_comparisons(d, pairs, dfam_rows, method):
    """pairs: [(row, ctrl_cell, trt_cell)], d has a 'cell' column, all cells of the ANOVA in d.
    method 'dunnett' (multivariate t, one family) or 'sidak'. Returns per pair: t, raw p, p (family = all
    comparisons of the panel, Prism's default), p (family = each row)."""
    cells = list(dict.fromkeys(d['cell']))
    mse, dfe = pooled(d, cells)
    info = []
    for row, c0, c1 in pairs:
        a, b = d[d['cell'] == c0]['value'], d[d['cell'] == c1]['value']
        se = math.sqrt(mse * (1 / len(a) + 1 / len(b)))
        t = (b.mean() - a.mean()) / se
        info.append(dict(row=row, c0=c0, c1=c1, n0=len(a), n1=len(b), t=t, p_raw=2 * SS.t.sf(abs(t), dfe)))

    def adjust(idx):
        k = len(idx)
        out = {}
        if method == 'sidak':
            for i in idx:
                out[i] = 1 - (1 - info[i]['p_raw']) ** k
            return out
        R = np.eye(k)
        for a_, i in enumerate(idx):
            for b_, j in enumerate(idx):
                if a_ != b_ and info[i]['c0'] == info[j]['c0']:
                    n0, ni, nj = info[i]['n0'], info[i]['n1'], info[j]['n1']
                    R[a_, b_] = math.sqrt(ni * nj / ((ni + n0) * (nj + n0)))
        if k > 1:
            # Dunnett: P(max_j |T_j| >= |t|), T ~ multivariate t (df = pooled ANOVA df, correlation R), by seeded
            # simulation (400 000 draws); clipped to its exact bounds [raw P, Sidak P].
            rng = np.random.default_rng(1)
            L = np.linalg.cholesky(R)
            mx = []
            for _ in range(4):
                z = rng.standard_normal((100000, k)) @ L.T
                s = np.sqrt(rng.chisquare(dfe, 100000) / dfe)
                mx.append(np.abs(z / s[:, None]).max(axis=1))
            mx = np.sort(np.concatenate(mx))
        for a_, i in enumerate(idx):
            tt = abs(info[i]['t'])
            if k == 1:
                p = info[i]['p_raw']
            else:
                p = 1 - np.searchsorted(mx, tt, side='left') / len(mx)
                p = min(max(p, info[i]['p_raw']), 1 - (1 - info[i]['p_raw']) ** k)
            out[i] = p
        return out
    p_all = adjust(list(range(len(info))))
    p_row = {}
    for row in dict.fromkeys(x['row'] for x in info):
        p_row.update(adjust([i for i, x in enumerate(info) if x['row'] == row]))
    for i, x in enumerate(info):
        x['p'] = p_all[i]
        x['p_row'] = p_row[i]
    return info, mse, dfe


def bky(pvals, Q=0.01):
    from statsmodels.stats.multitest import multipletests
    rej, q, _, _ = multipletests(pvals, alpha=Q, method='fdr_tsbky')
    return rej, q


# ---------------------------------------------------------------- chart page sizes
SIZES = {r['chart']: (float(r['w_mm']), float(r['h_mm'])) for r in csv.DictReader(open(SIZES_CSV, newline=''))}
DRAWN = []


# ---------------------------------------------------------------- drawing
FS = ST['font']['size_pt']


def bbox_mm(page):
    r = page.fig.canvas.get_renderer()
    bb = page.fig.get_tightbbox(r)          # inches, origin bottom-left
    x0, x1 = bb.x0 * 25.4, bb.x1 * 25.4
    y_top, y_bot = page.H - bb.y1 * 25.4, page.H - bb.y0 * 25.4
    return x0, y_top, x1, y_bot


def render(chart, W, H, df, spec, ann=None, title=None, top=None, extra=None, left=12.0, right=0.8, title_pt=None):
    """Two passes: draw with a guessed frame, measure the ink, then fit the frame so that the ink fills the page
    (0.3 mm margin) without clipping. Returns (page, record)."""
    PAD = 0.3
    has_key = spec.get('legend', spec.get('colour_by', 'series') == 'series') and list(spec['series']) != ['all']
    top = top if top is not None else (2.0 + (3.2 if has_key else 0) + (3.0 if title else 0))
    box = [left, top, W - left - right, H - top - 14.0]
    for it in range(6):
        page = Page(W, H, chart)
        ax, rec = dotbar.draw(page, tuple(box), df, spec, ann, gid=chart)
        if title:
            page.text(box[0] + box[2] / 2, 0.6, title, ha='center', va='top',
                      gid=f'{chart}-title', **({'size': title_pt} if title_pt else {}))
        if extra:
            extra(page, rec, box)
        x0, y0, x1, y1 = bbox_mm(page)
        dl, dr, db, dt = x0 - PAD, (W - PAD) - x1, (H - PAD) - y1, y0 - PAD
        if (max(abs(dl), abs(dr), abs(db)) < 0.15 and dt >= -0.05) or it == 5:
            break
        box[0] -= dl
        box[2] += dl + dr
        box[3] += db
        if dt < 0:              # something above the page top: push the frame down
            box[1] -= dt
            box[3] += dt
        import matplotlib.pyplot as plt
        plt.close(page.fig)
    x0, y0, x1, y1 = bbox_mm(page)
    assert x0 >= 0 and y0 >= 0 and x1 <= W + 0.01 and y1 <= H + 0.01, (chart, x0, y0, x1, y1, W, H)
    page.save(OUT, chart)
    rec['frame_mm'] = [round(v, 2) for v in box]
    rec['ink_mm'] = [round(v, 2) for v in (x0, y0, x1, y1)]
    return page, rec


def axis(vals_max, headroom):
    return dotbar.nice_axis(vals_max * headroom)


def ann_df(rows):
    a = pd.DataFrame(rows, columns=['group', 'series', 'label', 'kind', 'group2', 'series2'])
    return a.fillna('')


SOURCE_ROWS = {}


def source(chart, d, cols_extra=()):
    keep = ['sheet', 'panel', 'unit', 'group', 'series', 'replicate', 'value', 'sample_id', 'cell', 'notes']
    out = d[keep + [c for c in cols_extra if c not in keep]].copy()
    out.insert(0, 'chart', chart)
    (OUT / 'data').mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT / 'data' / f'{chart}.csv', index=False, lineterminator='\n')
    SOURCE_ROWS[chart] = len(out)


def register(chart):
    """Page size (mm) of a chart."""
    DRAWN.append(chart)
    return SIZES[chart]


NOTES = []


# ================================================================ ED2b / ED2c  qPCR (multiple unpaired t, BKY 1 %)
def ed2_qpcr(chart, panel):
    W, H = register(chart)
    d = take(panel)
    check_summary(d, panel)
    genes = list(dict.fromkeys(d['group']))
    pv, pw = [], []
    for g in genes:
        a = d[(d['group'] == g) & (d['series'] == 'Col-0')]['value']
        b = d[(d['group'] == g) & (d['series'] == 'fer-4')]['value']
        pv.append(SS.ttest_ind(a, b, equal_var=True).pvalue)
        pw.append(SS.ttest_ind(a, b, equal_var=False).pvalue)
    rej, q = bky(pv, 0.01)
    rows = []
    for g, p, p2, r_, qq in zip(genes, pv, pw, rej, q):
        lab = stars(p) if r_ else 'ns'
        au = auth(chart).get(g, '')
        if HAVE_LABELS:
            rows.append([g, 'fer-4', au, 'stars' if au != 'ns' else 'ns', '', ''])
        chk(chart, panel, f'{g}: fer-4 vs Col-0', lab,
            "multiple unpaired t tests (Student, per gene) + two-stage step-up BKY FDR, Q = 1 %",
            p=float(p), p_alt=float(p2), drawn=au,
            note=f'q(BKY) = {qq:.3g}; discovery = {bool(r_)}; stars from the raw t-test P of discoveries '
                 f'(Prism prints discoveries); p_alternative = Welch; n = '
                 + '/'.join(str(int(((d['group'] == g) & (d['series'] == s)).sum())) for s in ('Col-0', 'fer-4')))
    vmax = d.groupby(['group', 'series'])['value'].agg(lambda v: max(v.max(), v.mean() + v.std(ddof=1))).max()
    ymax, ystep = axis(vmax, 1.1)
    spec = dict(groups=genes, series=['Col-0', 'fer-4'], colour_by='series', legend=True,
                y_label='Relative expression', y_max=ymax, y_step=ystep, x_rotation=45, group_gap=0.7,
                star_rotation=0, label_kind='genes', y_tick_format='{:g}')
    title = {'ED2b': 'Upregulated', 'ED2c': 'Downregulated'}[chart]
    with ed2_colours():
        render(chart, W, H, d, spec, ann_df(rows), title=title, title_pt=8.0)
    source(chart, d)


# ================================================================ ED2e  reporter quantification (Student's t)
def ed2e():
    d0 = take('ED2e')
    check_summary(d0, 'ED2e')
    for chart, rep in (('ED2e_pWOL', 'pWOL:erGFP'), ('ED2e_AUX1', 'pAUX1:AUX1-Venus')):
        # both charts are drawn at ONE page size (data/ed/ED_chart_sizes.csv) so that bars, text and axes match
        W, H = register(chart)
        d = d0[d0['series'] == rep].copy()
        d['s_orig'] = d['series']
        d['series'] = 'all'
        a, b = d[d['group'] == 'Col-0']['value'], d[d['group'] == 'fer-4']['value']
        p = SS.ttest_ind(a, b, equal_var=True).pvalue
        lab = stars(p)
        au = auth(chart).get('fer-4', '')
        chk(chart, 'ED2e', f'{rep}: fer-4 vs Col-0', lab, "Student's t test (two-sided, equal variance)", p=float(p),
            drawn=au,
            p_alt=float(SS.ttest_ind(a, b, equal_var=False).pvalue),
            note=f'n = {len(a)}/{len(b)} (Col-0/fer-4); legend n = 14-15'
                 + ('; AUX1-Venus n = 16 per genotype is OUTSIDE the legend range' if len(a) > 15 or len(b) > 15 else ''))
        vmax = max(a.max(), b.max(), a.mean() + a.std(), b.mean() + b.std())
        ymax, ystep = axis(vmax, 1.15)
        spec = dict(groups=['Col-0', 'fer-4'], series=['all'], colour_by='group', legend=False,
                    y_label=['Fluorescence intensity', 'per unit area (a.u.)'], y_max=ymax, y_step=ystep,
                    x_rotation=0, group_gap=0.4, star_rotation=0)
        title = {'pWOL:erGFP': '*pWOL:erGFP*', 'pAUX1:AUX1-Venus': '*pAUX1:AUX1-Venus*'}[rep]
        rows = [['Col-0', 'all', au, 'stars', 'fer-4', 'all']] if HAVE_LABELS else []
        with ed2_colours():
            render(chart, W, H, d, spec, ann_df(rows), title=title)
        source(chart, d.assign(series=d['s_orig']))


# ================================================================ ED4a / ED4b / ED4d  one-way ANOVA + Tukey letters
L1 = [f'{p}_L1' for p in ('pWER-KO', 'pCO2-KO', 'pSCR-KO', 'pSHR-KO', 'pSMB-KO')]
L2 = [f'{p}_L2' for p in ('pWER-KO', 'pCO2-KO', 'pSCR-KO', 'pSHR-KO', 'pSMB-KO')]
L3 = [f'{p}_L3' for p in ('pWER-KO', 'pCO2-KO', 'pSCR-KO', 'pSHR-KO', 'pSMB-KO')]
CTRL = ['Col-0', 'fer-4', 'FER/fer-4']
FIG2D_LABELS = KIT / 'data/fig2/fig2d_meristem_cell_number_stats.csv'
FIG2F_LABELS = KIT / 'data/fig2/fig2f_root_hair_length_stats.csv'


def calibrate_letters(panel, groups_shown, groups_all, author_csv, what):
    """Do the printed letters of the main-figure panel come from an ANOVA over the shown bars or over every line in
    the sheet? Writes the comparison to NOTES."""
    auth_ = pd.read_csv(author_csv)
    auth_ = dict(zip(auth_['group'], auth_['label']))
    d = take(panel)
    res = {}
    for fam, gs in (('shown', groups_shown), ('all', groups_all)):
        dd = d[d['group'].isin(gs)].copy()
        dd['cell'] = dd['group']
        lt, _, _ = tukey_letters(dd, 'cell', gs)
        match = sum(_same_partition(lt, auth_, groups_shown, a, b) for a in groups_shown for b in groups_shown if a < b)
        tot = len(groups_shown) * (len(groups_shown) - 1) // 2
        res[fam] = (match, tot, {g: lt[g] for g in groups_shown})
    NOTES.append(f'calibration {what}: author letters {auth_}; Tukey over shown bars -> {res["shown"][2]} '
                 f'({res["shown"][0]}/{res["shown"][1]} pairwise sig. calls agree); over all lines -> '
                 f'{res["all"][2]} ({res["all"][0]}/{res["all"][1]})')
    return res


def _same_partition(lt, auth_, gs, a, b):
    """1 when the 'share a letter' call for a pair agrees between two letter displays."""
    s1 = bool(set(lt[a]) & set(lt[b]))
    s2 = bool(set(auth_[a]) & set(auth_[b]))
    return int(s1 == s2)


def pair_agree(chart, comp, auth_, cells):
    """Letters are not unique, so compare the pairwise 'differ / do not differ' calls of the two displays."""
    bad = [f'{a} vs {b}' for i, a in enumerate(cells) for b in cells[i + 1:]
           if bool(set(comp[a]) & set(comp[b])) != bool(set(auth_[a]) & set(auth_[b]))]
    tot = len(cells) * (len(cells) - 1) // 2
    NOTES.append(f'{chart}: author vs computed letters agree on {tot - len(bad)}/{tot} pairwise calls'
                 + (f'; differ: {bad}' if bad else ''))
    return tot - len(bad), tot, bad


def oneway_chart(chart, panel, series_val, groups, y_label, legend_n, unit_note='',
                 cell_values=None, y_fmt='{:g}', headroom=1.18, letter_rot=None, family_groups=None, label_offset=None):
    W, H = register(chart)
    d = take(panel, series_val)
    check_summary(d, panel)
    d = d[d['group'].isin(groups)].copy()
    fam = family_groups or groups
    dt = (cell_values if cell_values is not None else take(panel, series_val))
    dt = dt[dt['group'].isin(fam)].copy()
    dt['cell'] = dt['group']
    lt, res, (F, p) = tukey_letters(dt, 'cell', fam)
    rows = []
    au = auth(chart)
    if HAVE_LABELS:
        ok, tot, bad = pair_agree(chart, lt, au, groups)
    for g in groups:
        n = int((d['group'] == g).sum())
        if HAVE_LABELS:
            rows.append([g, d['series'].iloc[0], au[g], 'letters', '', ''])
            agree = ('yes' if au[g] == lt[g] else
                     ('equivalent' if not [b for b in bad if g in b.split(' vs ')] else 'no'))
            pre = f'pairwise calls agree {ok}/{tot}; '
        else:
            agree, pre = None, ''
        chk(chart, panel, g, lt[g], "one-way ANOVA + Tukey's multiple comparisons (alpha 0.05), compact letters",
            p=float(p), drawn=au.get(g, ''), agree=agree,
            note=pre + f'ANOVA F = {F:.4g} over {len(fam)} groups; n = {n} values' + unit_note
                 + ('' if n == legend_n else f'; legend n = {legend_n}'))
    vmax = d.groupby('group')['value'].agg(lambda v: max(v.max(), v.mean() + v.std(ddof=1))).max()
    ymax, ystep = axis(vmax, headroom)
    longest = max(len(v) for v in (au if HAVE_LABELS else {g: lt[g] for g in groups}).values())
    spec = dict(groups=groups, series=[d['series'].iloc[0]], colour_by='group', legend=False, y_label=y_label,
                y_max=ymax, y_step=ystep, x_rotation=45, group_gap=0.275, require_labels=HAVE_LABELS,
                y_tick_format=y_fmt,
                letter_rotation=(letter_rot if letter_rot is not None else (90 if longest > 2 else 0)),
                letter_rotation_min_len=1)
    if label_offset is not None:                 # gap between the highest dot and its letters (kit default 0.8 mm)
        spec['label_offset_mm'] = label_offset
    render(chart, W, H, d, spec, ann_df(rows))
    source(chart, d)
    return lt


# ================================================================ two-way ANOVA charts with control comparisons
def control_chart(chart, panel, groups, series, pairs_fn, method, test_name, legend_n, y_label,
                  bracket=False, d=None, headroom=1.12, legend_title=None, x_rot=45, group_gap=0.7, extra=None):
    W, H = register(chart)
    if d is None:
        d = take(panel)
        check_summary(d, panel)
    d = d[d['group'].isin(groups) & d['series'].isin(series)].copy()
    d['cell'] = d['group'] + ' | ' + d['series']
    pairs = pairs_fn(groups, series)
    info, mse, dfe = control_comparisons(d, pairs, None, method)
    rows = []
    # family: Dunnett one per row (genotype / cell type), Sidak all comparisons of the panel -- the choices that
    # reproduce every printed label (checked against the other choice in p_alternative)
    if method == 'sidak':
        for x in info:
            x['p'], x['p_row'] = x['p_row'], x['p']
    fam_txt = 'all comparisons in the panel' if method == 'sidak' else (
        'one per ' + ('genotype' if not panel.startswith('ED10') else 'cell type'))
    alt_txt = ('one per genotype' if method == 'sidak' else 'all comparisons in the panel')
    for x in info:
        g0, s0 = x['c0'].split(' | ')
        g1, s1 = x['c1'].split(' | ')
        lab = stars(x['p_row'])     # one family per row (genotype / cell type): reproduces the printed labels
        au = auth(chart).get((g1, s1), '')
        kind = 'ns' if au == 'ns' else 'stars'
        if au:
            if bracket:
                rows.append([g0, s0, au, kind, g1, s1])
            else:
                rows.append([g1, s1, au, kind, '', ''])
        ns = f"n = {x['n0']}/{x['n1']}"
        flag = '' if stars(x['p']) == stars(x['p_row']) else f"; LABEL DIFFERS with family = {alt_txt}: {stars(x['p'])}"
        chk(chart, panel, f"{x['c1']} vs {x['c0']}", lab, test_name, p=float(x['p_row']), p_alt=float(x['p']),
            drawn=au, note=f"t = {x['t']:.3g}, raw P = {x['p_raw']:.3g}, pooled MSE df = {dfe}; family = {fam_txt} "
                 f"(p_alternative = family {alt_txt}){flag}; "
                 f"{ns}; legend n = {legend_n}")
    vmax = d.groupby('cell')['value'].agg(lambda v: max(v.max(), v.mean() + v.std(ddof=1))).max()
    ymax, ystep = axis(vmax, headroom)
    spec = dict(groups=groups, series=series, colour_by='series', legend=True, y_label=y_label, y_max=ymax,
                y_step=ystep, x_rotation=x_rot, group_gap=group_gap, star_rotation=90,
                legend_title=legend_title)
    if HAVE_LABELS:
        assert len(rows) == len(AUTH[chart]), chart
    render(chart, W, H, d, spec, ann_df(rows), extra=extra)
    source(chart, d)
    return d


def vs_first_series(groups, series):
    return [(g, f'{g} | {series[0]}', f'{g} | {s}') for g in groups for s in series[1:]]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    OUT_STATS.mkdir(parents=True, exist_ok=True)
    if not HAVE_LABELS:
        print(f'note: {LABELS_CSV.relative_to(KIT)} not found -> charts are drawn without statistics labels; '
              f'STATS_CHECK_ED.csv lists the computed labels only')
    # ---- ED2
    ed2_qpcr('ED2b', 'ED2b')
    ed2_qpcr('ED2c', 'ED2c')
    ed2e()

    # ---- ED4a / b: does Fig 2d (same sheet) use an ANOVA over the 8 shown bars or all 13 lines?
    fam_a = CTRL + L2                      # default (and the result of the check below): the 8 shown bars
    if HAVE_LABELS:
        if FIG2D_LABELS.exists():
            calibrate_letters('Fig2d;ED4a', CTRL + L1, CTRL + L1 + L2, FIG2D_LABELS, 'Fig 2d meristem cell number')
        da = take('Fig2d;ED4a')
        sc = {}
        for nm, gs in (('8 shown bars', CTRL + L2), ('all 13 lines', CTRL + L2 + L1)):
            x = da[da['group'].isin(gs)].assign(cell=lambda t: t['group'])
            lt_, _, _ = tukey_letters(x, 'cell', gs)
            sc[nm] = (sum(_same_partition(lt_, AUTH['ED4a'], None, a, b) for a in CTRL + L2 for b in CTRL + L2 if a < b), gs)
        NOTES.append('ED4a family check vs the ED4a author letters (pairwise calls of 28): '
                     + '; '.join(f'{k}: {v[0]}' for k, v in sc.items()))
        fam_a = max(sc.values(), key=lambda v: v[0])[1]
    NOTES.append(f'ED4a letters computed over {"all 13 lines (better match to the Fig 2d author letters)" if len(fam_a) == 13 else "the 8 shown bars"}')
    oneway_chart('ED4a', 'Fig2d;ED4a', 'all', CTRL + L2, 'Meristem cell number', 30, family_groups=fam_a)
    order_b = CTRL + [x for p in zip(L1, L2) for x in p]
    oneway_chart('ED4b', 'ED4b', 'all', order_b, 'Meristem length (µm)', 30)

    # ---- ED4d root hair L3 set. Unit: 100 hairs (as plotted) or 10 root means (legend n = 10)?
    hair = take('Fig2f;ED4d', 'L3 set')
    check_summary(hair, 'Fig2f;ED4d')
    use_means = False                      # default (and the result of the check below): Tukey on the 100 hairs
    if HAVE_LABELS and FIG2F_LABELS.exists():
        l1 = take('Fig2f;ED4d', 'L1 set')
        auth_ = pd.read_csv(FIG2F_LABELS)
        auth_ = dict(zip(auth_['group'], auth_['label']))
        gl1 = CTRL + L1
        cal_rh = {}
        g3 = CTRL + L3
        for unit, dd, dd3 in (('hairs', l1, hair),
                              ('root means', l1.groupby(['group', 'sample_id'], as_index=False)['value'].mean(),
                               hair.groupby(['group', 'sample_id'], as_index=False)['value'].mean())):
            dd = dd.copy(); dd['cell'] = dd['group']
            lt, _, _ = tukey_letters(dd, 'cell', gl1)
            m = sum(_same_partition(lt, auth_, gl1, a, b) for a in gl1 for b in gl1 if a < b)
            dd3 = dd3.copy(); dd3['cell'] = dd3['group']
            lt3, _, _ = tukey_letters(dd3, 'cell', g3)
            m3 = sum(_same_partition(lt3, AUTH['ED4d'], g3, a, b) for a in g3 for b in g3 if a < b)
            cal_rh[unit] = (m + m3, lt, m)
            NOTES.append(f'ED4d unit check: Tukey on {unit}: Fig 2f L1 {m}/28, ED4d L3 {m3}/28 pairwise calls agree '
                         f'with the author letters')
        NOTES.append(f'calibration Fig 2f root hair (L1 set) author letters {auth_}: Tukey on 100 hairs -> {cal_rh["hairs"][1]} '
                     f'({cal_rh["hairs"][2]}/28 pairs agree); on 10 root means -> {cal_rh["root means"][1]} '
                     f'({cal_rh["root means"][2]}/28)')
        use_means = cal_rh['root means'][0] > cal_rh['hairs'][0]   # ties -> 100 hairs (what is plotted)
    tvals = (hair.groupby(['group', 'sample_id'], as_index=False)['value'].mean() if use_means else hair)
    oneway_chart('ED4d', 'Fig2f;ED4d', 'L3 set', CTRL + L3, 'Root hair length (mm)', 10,
                 unit_note=('; Tukey run on the 10 per-root means (legend n = 10 roots); dots = all 100 hairs'
                            if use_means else '; Tukey run on all 100 hairs (10 roots x 10 hairs); legend n = 10 roots'),
                 cell_values=tvals.assign(series='L3 set'), y_fmt='{:.1f}', label_offset=2.0)

    # ---- ED4f RALF23 TSKO L3; two-way ANOVA + Dunnett, RALF23 vs Mock within genotype
    g4f = CTRL + L3
    control_chart('ED4f', 'ED4f', g4f, ['Mock', 'RALF23'], vs_first_series, 'dunnett',
                  "two-way ANOVA (genotype x treatment) + Dunnett's test, RALF23 vs Mock within each genotype",
                  '14-16', 'Root length (mm)')

    # ---- ED5g yuc mutants: two-way ANOVA + Sidak (legend also says 't test')
    control_chart('ED5g', 'ED5g', ['Col-0', 'yuc5', 'yuc2589', 'yuc35789 (yucQ)'], ['Mock', 'RALF23'],
                  vs_first_series, 'sidak',
                  "two-way ANOVA (genotype x treatment) + Sidak's multiple comparisons, RALF23 vs Mock within genotype "
                  "(legend ALSO says 't test': conflict, ANOVA used)", '13-15', 'Root length (mm)', x_rot=0)

    # ---- ED6c PEO-IAA: two-way ANOVA + Dunnett vs RALF23 within genotype
    def vs_ralf(groups, series):
        return [(g, f'{g} | RALF23', f'{g} | {s}') for g in groups for s in ('Mock', 'RALF23 & PEO-IAA')]
    control_chart('ED6c', 'ED6c', ['FER/fer-4', 'fer-4'], ['Mock', 'RALF23', 'RALF23 & PEO-IAA'],
                  vs_ralf, 'dunnett',
                  "two-way ANOVA (genotype x treatment) + Dunnett's test, Mock and RALF23 & PEO-IAA vs RALF23 within "
                  "genotype", '13-15', 'Root length (mm)', bracket=True, x_rot=0, headroom=1.15)

    # ---- ED7b multi-RALF root length
    control_chart('ED7b', 'ED7b', ['Col-0', 'fer-4'], ['Mock', 'RALF1', 'RALF23', 'RALF31', 'RALF33'],
                  vs_first_series, 'dunnett',
                  "two-way ANOVA (genotype x treatment) + Dunnett's test, each RALF vs Mock within genotype",
                  '13-15', 'Root length (mm)', x_rot=0)

    # ---- ED7d / ED7f: one-way ANOVA + Tukey, treatments on x
    for chart, panel, ylab, legn in (
            ('ED7d', 'ED7d', ['*pYUC5:GFP-GUS* intensity', 'per unit area (a.u.)'], 15),
            ('ED7f', 'ED7f', ['*DR5rev::GFP* intensity', 'per unit area (a.u.)'], 10)):
        W, H = register(chart)
        d = take(panel)
        check_summary(d, panel)
        tr = list(dict.fromkeys(d['series']))
        d = d.assign(s_orig=d['series'], group=d['series'], series='all')
        dd = d.assign(cell=d['group'])
        lt, _, (F, p) = tukey_letters(dd, 'cell', tr)
        au = auth(chart)
        if HAVE_LABELS:
            pair_agree(chart, lt, au, tr)
        for t in tr:
            n = int((d['group'] == t).sum())
            chk(chart, panel, t, lt[t], "one-way ANOVA + Tukey's multiple comparisons (alpha 0.05), compact letters",
                p=float(p), drawn=au.get(t, ''),
                note=f'ANOVA F = {F:.4g}; n = {n}' + ('' if n == legn else f'; legend n = {legn}'))
        vmax = d.groupby('group')['value'].agg(lambda v: max(v.max(), v.mean() + v.std(ddof=1))).max()
        ymax, ystep = axis(vmax, 1.15)
        spec = dict(groups=tr, series=['all'], colour_by='group', legend=False, y_label=ylab, y_max=ymax,
                    y_step=ystep, x_rotation=0, group_gap=0.3, label_kind='series', require_labels=HAVE_LABELS)
        render(chart, W, H, d, spec, ann_df([[t, 'all', au[t], 'letters', '', ''] for t in tr] if HAVE_LABELS else []))
        source(chart, d.assign(group='Col-0', series=d['s_orig']))

    # ---- ED9d DOF4.6-OE: one-way ANOVA + Tukey on all 6 bars; ratio labels in the RALF23 columns
    W, H = register('ED9d')
    d = take('ED9d')
    check_summary(d, 'ED9d')
    gs = ['Col-0', 'DOF4.6-OE #17', 'DOF4.6-OE #24']
    d['cell'] = d['group'] + ' | ' + d['series']
    cells = [f'{g} | {s}' for g in gs for s in ('Mock', 'RALF23')]
    lt, _, (F, p) = tukey_letters(d, 'cell', cells)
    au = auth('ED9d')
    if HAVE_LABELS:
        pair_agree('ED9d', lt, au, cells)
    printed = PRINTED.get('ED9d', {})      # ratios printed beside the published chart
    ratios = {}
    for g in gs:
        ratios[g] = d[d['cell'] == f'{g} | RALF23']['value'].mean() / d[d['cell'] == f'{g} | Mock']['value'].mean()
    for c in cells:
        n = int((d['cell'] == c).sum())
        chk('ED9d', 'ED9d', c, lt[c], "one-way ANOVA + Tukey's multiple comparisons over the 6 bars (alpha 0.05)",
            p=float(p), drawn=au.get(c, ''), note=f'ANOVA F = {F:.4g}; n = {n}; legend n = 15-17'
                             + (' -> OUTSIDE the legend range' if not 15 <= n <= 17 else ''))
    for g in gs:
        pr = printed.get(g, '')
        chk('ED9d', 'ED9d', f'{g}: ratio RALF23/Mock', f'{ratios[g]:.2f}', 'mean(RALF23) / mean(Mock), source data',
            printed=pr, agree='yes' if f'{ratios[g]:.2f}' == pr else 'no',
            note='printed value = ratio printed on the published chart (panel d); drawn = computed')
    vmax = d.groupby('cell')['value'].agg(lambda v: max(v.max(), v.mean() + v.std(ddof=1))).max()
    ymax, ystep = axis(vmax, 1.2)

    def ratio_text(page, rec, box):
        for b in rec['bars']:
            if b['series'] == 'RALF23':
                ymm = box[1] + box[3] - 0.35 * (box[1] + box[3] - b['top_mm'])
                page.text(b['x_mm'], ymm, f"{ratios[b['group']]:.2f}", ha='center', va='center',
                          size=FS * 0.92, gid='ED9d-ratio')
    spec = dict(groups=gs, series=['Mock', 'RALF23'], colour_by='series', legend=True, y_label='Root length (mm)',
                y_max=ymax, y_step=ystep, x_rotation=45, group_gap=0.6, require_labels=HAVE_LABELS)
    render('ED9d', W, H, d, spec, ann_df([[c.split(' | ')[0], c.split(' | ')[1], au[c], 'letters', '', '']
                                          for c in cells] if HAVE_LABELS else []), extra=ratio_text)
    source('ED9d', d)

    # ---- ED9e DOF qPCR: per gene two-way ANOVA (genotype x treatment) + Tukey over the 4 bars
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm
    d0 = take('ED9e')
    check_summary(d0, 'ED9e')
    genes = list(dict.fromkeys(d0['sample_id']))      # DOF4.1, COG1, OBP3, DOF6, OBP1, OBP2 (left to right)
    for gi, gene in enumerate(genes):
        chart = f'ED9e_{gene}'
        W, H = register(chart)
        d = d0[d0['sample_id'] == gene].copy()
        d['cell'] = d['group'] + ' | ' + d['series']
        cells = [f'{g} | {s}' for g in ('Col-0', 'fer-4') for s in ('Mock', 'RALF23')]
        lt, _, _ = tukey_letters(d, 'cell', cells)
        au = auth(chart)
        if HAVE_LABELS:
            pair_agree(chart, lt, au, cells)
        aov = anova_lm(smf.ols('value ~ C(group) * C(series)', data=d).fit(), typ=2)
        pint = float(aov.loc['C(group):C(series)', 'PR(>F)'])
        for c in cells:
            n = int((d['cell'] == c).sum())
            chk(chart, 'ED9e', f'{gene}: {c}', lt[c],
                "two-way ANOVA (genotype x treatment, per gene) + Tukey's multiple comparisons over the 4 bars",
                p=pint, drawn=au.get(c, ''), note=f'p_value = interaction P; n = {n}' + ('' if 2 <= n <= 3 else '; legend n = 2-3'))
        vmax = d.groupby('cell')['value'].agg(lambda v: max(v.max(), v.mean() + v.std(ddof=1))).max()
        ymax, ystep = axis(vmax, 1.2)
        spec = dict(groups=['Col-0', 'fer-4'], series=['Mock', 'RALF23'], colour_by='series', legend=False,
                    y_label=('Relative expression' if gi == 0 else ''), y_max=ymax, y_step=ystep, x_rotation=0,
                    group_gap=0.5, require_labels=HAVE_LABELS, y_tick_format=('{:.1f}' if ystep < 1 else '{:g}'))
        render(chart, W, H, d, spec, ann_df([[c.split(' | ')[0], c.split(' | ')[1], au[c], 'letters', '', '']
                                             for c in cells] if HAVE_LABELS else []), title=f'*{gene}*', top=5.4,
               title_pt=10.0)          # gene names at 10 pt so that they stay legible at the printed size
        source(chart, d)
    # ED9e key: a separate key-only page (Mock / RALF23) placed under the six charts
    W, H = register('ED9e_key')
    page = Page(W, H, 'ED9e_key')
    from matplotlib.lines import Line2D
    msize = ST['markers']['diameter_pt']
    items = [('Mock', S.series_colour('Mock')), ('RALF23', S.series_colour('RALF23'))]
    widths = [page.text_width(t) for t, _ in items]
    mk = msize / S.PT_PER_MM
    total = sum(widths) + len(items) * (mk + 0.9) + 2.4
    cx, yb = W / 2 - total / 2, H / 2 + 0.36 * FS / S.PT_PER_MM
    for (t, c), tw in zip(items, widths):
        page.fig.add_artist(Line2D([page.fx(cx + mk / 2)], [page.fy(yb - 0.36 * FS / S.PT_PER_MM)],
                                   transform=page.fig.transFigure, marker='o', markersize=msize * 1.25,
                                   markerfacecolor='none', markeredgecolor=c,
                                   markeredgewidth=ST['lines']['marker_edge_pt'] * 1.2, linestyle='none'))
        page.text(cx + mk + 0.9, yb, t)
        cx += mk + 0.9 + tw + 2.4
    page.save(OUT, 'ED9e_key')

    # ---- ED10c PI puncta, RALF23-treated roots, by cell type; Dunnett vs FER/fer-4 within cell type.
    # Plotted: RALF23-treated roots only, the 7 genotypes of the image row; the sheet also holds Mock and Col-0
    # columns, which are not drawn.
    d = take('ED10c')
    check_summary(d, 'ED10c')
    d = d.copy()
    d['tissue'] = d['notes'].str.replace('tissue=', '', regex=False)
    d = d[d['series'] == 'RALF23']
    gens = ['FER/fer-4', 'fer-4'] + L1
    tissues = ['Epidermis', 'Cortex', 'Endodermis', 'Stele']
    dd = d.assign(treatment=d['series'], genotype=d['group'], group=d['tissue'], series=d['group'])
    for t in tissues:
        for g in gens:
            assert int(((dd['group'] == t) & (dd['series'] == g)).sum()) == 10, (t, g)

    def vs_fer(groups, series):
        return [(t, f'{t} | FER/fer-4', f'{t} | {g}') for t in groups for g in series[1:]]
    control_chart('ED10c', 'ED10c', tissues, gens, vs_fer, 'dunnett',
                  "two-way ANOVA (genotype x cell type, RALF23-treated roots) + Dunnett's test vs FER/fer-4 within "
                  "each cell type", 10, 'Number of PI puncta', d=dd, x_rot=0, group_gap=1.0, headroom=1.1)

    # ---- tables
    pd.DataFrame(CHECK).to_csv(OUT_STATS / 'STATS_CHECK_ED.csv', index=False, lineterminator='\n')
    (OUT_STATS / 'NOTES_ED.txt').write_text('\n'.join(NOTES) + '\n')
    print('\n'.join(NOTES))
    missing = sorted(set(SIZES) - set(DRAWN))
    assert not missing, f'charts in {SIZES_CSV.name} that were not drawn: {missing}'
    print(f'{len(DRAWN)} charts -> {OUT.relative_to(KIT)}; source rows: {SOURCE_ROWS}')


if __name__ == '__main__':
    main()
