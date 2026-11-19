#!/usr/bin/env python
"""Fig. 2 measurement tables and statistics check.

Reads the tidy source data in source_data/long/ (written by tools/source_data_to_long.py from the Source data
workbook) and writes data/fig2/:

  fig2d_meristem_cell_number.csv   sheet 'Fig_2d&Ext_Fig_4a Meristem#', the 8 *_L1 / control groups
  fig2f_root_hair_length.csv       sheet 'Fig_2f&Ext_Fig_4d Root hair(mm)', 'L1 set', one row per root hair
  fig2g_root_length_ralf23.csv     sheet 'Fig_2g RALF_CRISPR-L1(mm_', Mock + 1 uM RALF23

Every n and mean per group x series is asserted equal to source_data/long/summary.csv.

Statistics: the labels drawn in the figure are the authors' labels, stored in data/fig2/*_stats.csv and
data/fig2/fig2h_ratios.csv (read only, never rewritten here). This script recomputes them with the test named in the
legend and writes the comparison to output/stats_checks/STATS_CHECK_fig2.csv:
  d  one-way ANOVA + Tukey HSD on the 8 groups drawn (compact letter display)
  f  one-way ANOVA + Tukey HSD per root hair (the unit drawn), plus per-root means as an alternative unit
  g  two-way ANOVA (genotype x treatment), 1 uM RALF23 vs Mock within each genotype on the pooled residual MS,
     Sidak over the 6 genotypes
  h  ratios printed on the seedling photograph vs mean(1 uM RALF23) / mean(Mock) of the g table
If the label files are not present, the computed labels are written with label_drawn empty.

Usage (from the kit folder, after `python tools/source_data_to_long.py "<Source data.xlsx>"`):
  python tools/fig2_tables.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

KIT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(KIT))
from ferfigs.stats_helper import cld, stars  # noqa: E402

SRC = KIT / 'source_data' / 'long'
OUT = KIT / 'data' / 'fig2'
CHK = KIT / 'output' / 'stats_checks'
NO_LABELS = 'labels not in this repository'

L1 = ['Col-0', 'fer-4', 'FER/fer-4', 'pWER-KO_L1', 'pCO2-KO_L1', 'pSCR-KO_L1', 'pSHR-KO_L1', 'pSMB-KO_L1']
RALF_G = ['Col-0', 'pWER-KO_L1', 'pCO2-KO_L1', 'pSCR-KO_L1', 'pSHR-KO_L1', 'pSMB-KO_L1']
T1 = '1 µM RALF23'
CHECK = []
MISSING = []


def write(df, name):
    df.to_csv(OUT / name, index=False, lineterminator='\n')


def assert_summary(df, panel, series_map=None):
    """n and mean per group x series of the written table == summary.csv (source series names via series_map)."""
    s = pd.read_csv(SRC / 'summary.csv')
    s = s[s['panel'] == panel]
    for (g, ser), sub in df.groupby(['group', 'series']):
        src_ser = (series_map or {}).get(ser, ser)
        r = s[(s['group'] == g) & (s['series'] == src_ser)]
        assert len(r) == 1, (panel, g, ser)
        r = r.iloc[0]
        assert len(sub) == int(r['n']), (panel, g, ser, len(sub), r['n'])
        assert abs(sub['value'].mean() - r['mean']) < 1e-4, (panel, g, ser, sub['value'].mean(), r['mean'])


def drawn_labels(name, groups, series=None):
    """Authors' labels {group: label} from data/fig2/<name>; None when the file is not in this copy of the kit."""
    p = OUT / name
    if not p.exists():
        MISSING.append(name)
        return None
    st = pd.read_csv(p, dtype=str, keep_default_na=False)
    if series is not None:
        st = st[st['series'] == series]
    lab = dict(zip(st['group'], st['label']))
    assert set(groups) <= set(lab), (name, set(groups) - set(lab))
    return lab


def agree(drawn, computed):
    if drawn is None:
        return NO_LABELS
    return 'yes' if drawn == computed else 'no'


def letters_tukey(df, alpha=0.05):
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    from scipy.stats import f_oneway
    order = list(dict.fromkeys(df['group']))
    F, p = f_oneway(*[df[df['group'] == g]['value'] for g in order])
    tk = pairwise_tukeyhsd(df['value'], df['group'], alpha=alpha)
    res = pd.DataFrame(tk.summary().data[1:], columns=tk.summary().data[0])
    sig = {frozenset((r['group1'], r['group2'])) for _, r in res.iterrows() if bool(r['reject'])}
    lt = cld(order, [df[df['group'] == g]['value'].mean() for g in order], sig)
    return lt, F, p


def same_letters(a, b):
    """Two compact letter displays are equivalent when they imply the same 'shares a letter' relation."""
    ks = list(a)
    return all((bool(set(a[i]) & set(a[j])) == bool(set(b[i]) & set(b[j]))) for i in ks for j in ks)


def whole(drawn, computed):
    if drawn is None:
        return NO_LABELS
    return 'yes' if same_letters(drawn, computed) else 'no'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    CHK.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------- d: meristem cell number (L1 lines only)
    m = pd.read_csv(SRC / 'Fig_2d_Ext_Fig_4a_Meristem.csv')
    d = m[m['group'].isin(L1)].copy()
    d['group'] = pd.Categorical(d['group'], L1, ordered=True)
    d = d.sort_values(['group', 'replicate']).astype({'group': str})
    d_out = pd.DataFrame(dict(group=d['group'], series='all', replicate=d['replicate'].astype(int), value=d['value'],
                              unit='cells', sample_id='', experiment='',
                              notes='mean of left and right cell files; source cell ' + d['cell'] + ' (' + d['notes'] + ')'))
    assert sorted(set(d_out['group'])) == sorted(L1)
    assert_summary(d_out, 'Fig2d;ED4a')
    write(d_out, 'fig2d_meristem_cell_number.csv')
    auth = drawn_labels('fig2d_meristem_cell_number_stats.csv', L1)
    lt8, F8, p8 = letters_tukey(d_out)
    lt13, _, _ = letters_tukey(m[['group', 'value']])
    for g in L1:
        a = auth[g] if auth else ''
        CHECK.append(dict(panel='Fig2d', comparison=f'{g} (letters across the 8 groups drawn)', label_drawn=a,
                          label_computed=lt8[g], agree=agree(auth and a, lt8[g]),
                          note='legend names no test; computed = one-way ANOVA + Tukey HSD on the 8 groups drawn '
                               f'(F = {F8:.1f}, P = {p8:.2g}); with the 5 _L2 lines included (13 groups, as in the sheet) '
                               f'{g} gets "{lt13[g]}"'))
    CHECK.append(dict(panel='Fig2d', comparison='whole letter display', label_drawn='author' if auth else '',
                      label_computed='Tukey 8 groups', agree=whole(auth, lt8),
                      note='equivalent = same pairs share a letter (letter names may differ)'))

    # ---------------------------------------------------------------- f: root hair length, L1 set, one row per hair
    h = pd.read_csv(SRC / 'Fig_2f_Ext_Fig_4d_Root_hair_mm.csv')
    f = h[h['series'] == 'L1 set'].copy()
    assert sorted(set(f['group'])) == sorted(L1)
    f['group'] = pd.Categorical(f['group'], L1, ordered=True)
    f = f.sort_values(['group']).astype({'group': str})
    f_out = pd.DataFrame(dict(group=f['group'], series='all', replicate=f['replicate'], value=f['value'], unit='mm',
                              sample_id=f['sample_id'], experiment='',
                              notes='one of the 10 longest hairs of the root in sample_id; source cell ' + f['cell']))
    assert f_out.groupby('group').size().eq(100).all()
    assert f_out.groupby('group')['sample_id'].nunique().eq(10).all()
    assert_summary(f_out, 'Fig2f;ED4d', series_map={'all': 'L1 set'})
    write(f_out, 'fig2f_root_hair_length.csv')
    auth = drawn_labels('fig2f_root_hair_length_stats.csv', L1)
    ltf, Ff, pf = letters_tukey(f_out)
    roots = f_out.groupby(['group', 'sample_id'], sort=False)['value'].mean().reset_index()
    ltr, Fr, pr = letters_tukey(roots)
    for g in L1:
        a = auth[g] if auth else ''
        CHECK.append(dict(panel='Fig2f', comparison=f'{g} (letters across 8 genotypes)', label_drawn=a,
                          label_computed=ltf[g], agree=agree(auth and a, ltf[g]),
                          note=f'computed per hair (n = 100: F = {Ff:.1f}, P = {pf:.2g}); per-root means (n = 10) give "{ltr[g]}"'))
    CHECK.append(dict(panel='Fig2f', comparison='whole letter display, per hair (unit drawn)',
                      label_drawn='author' if auth else '', label_computed='Tukey per hair', agree=whole(auth, ltf),
                      note='the chart plots every hair (100 points per genotype); legend says n = 10 (roots)'))
    CHECK.append(dict(panel='Fig2f', comparison='whole letter display, per-root means',
                      label_drawn='author' if auth else '', label_computed='Tukey per-root mean', agree=whole(auth, ltr),
                      note='alternative unit matching the legend n = 10'))

    # ---------------------------------------------------------------- g: root length, Mock vs 1 uM RALF23
    r = pd.read_csv(SRC / 'Fig_2g_RALF_CRISPR_L1_mm.csv')
    assert set(r['series']) == {'Mock', 'RALF23'}
    r['series'] = r['series'].replace({'RALF23': T1})
    r['group'] = pd.Categorical(r['group'], RALF_G, ordered=True)
    r['series'] = pd.Categorical(r['series'], ['Mock', T1], ordered=True)
    r = r.sort_values(['group', 'series']).astype({'group': str, 'series': str})
    g_out = pd.DataFrame(dict(group=r['group'], series=r['series'], replicate=r['replicate'].astype(int),
                              value=r['value'], unit='mm', sample_id='', experiment='',
                              notes='seedling root length; source cell ' + r['cell'] + " (sheet header '1μM RALF23')"))
    assert_summary(g_out, 'Fig2g', series_map={T1: 'RALF23'})
    write(g_out, 'fig2g_root_length_ralf23.csv')
    auth = drawn_labels('fig2g_root_length_ralf23_stats.csv', RALF_G, series=T1)
    import statsmodels.formula.api as smf
    from scipy import stats as SS
    mdl = smf.ols('value ~ C(group) * C(series)', data=g_out).fit()
    mse, dfe = mdl.mse_resid, mdl.df_resid
    for gname in RALF_G:
        a = g_out[(g_out['group'] == gname) & (g_out['series'] == 'Mock')]['value']
        b = g_out[(g_out['group'] == gname) & (g_out['series'] == T1)]['value']
        t = (b.mean() - a.mean()) / np.sqrt(mse * (1 / len(a) + 1 / len(b)))
        p1 = 2 * SS.t.sf(abs(t), dfe)                       # uncorrected, pooled MS
        p6 = 1 - (1 - p1) ** len(RALF_G)                    # Sidak over the 6 genotypes
        lab = auth[gname] if auth else ''
        CHECK.append(dict(panel='Fig2g', comparison=f'{gname}: {T1} vs Mock', label_drawn=lab,
                          label_computed=stars(p6), agree=agree(auth and lab, stars(p6)),
                          note=f'two-way ANOVA residual MS (df {dfe:.0f}); Sidak over the 6 genotypes (legend test) '
                               f'P = {p6:.2g}; uncorrected pooled-MS t test P = {p1:.2g} ({stars(p1)})'))

    # ---------------------------------------------------------------- h: ratios printed on the seedling photograph
    printed = drawn_labels('fig2h_ratios.csv', RALF_G, series=T1)
    mean = g_out.groupby(['group', 'series'])['value'].mean()
    for gname in RALF_G:
        ratio = mean[(gname, T1)] / mean[(gname, 'Mock')]
        lab = printed[gname] if printed else ''
        CHECK.append(dict(panel='Fig2h', comparison=f'{gname}: printed ratio 1 µM / Mock', label_drawn=lab,
                          label_computed=f'{ratio:.2f}', agree=agree(printed and lab, f'{ratio:.2f}'),
                          note=f'mean(1 µM RALF23)/mean(Mock) from the Fig. 2g source data = {ratio:.4f}; printed value kept'))

    pd.DataFrame(CHECK, columns=['panel', 'comparison', 'label_drawn', 'label_computed', 'agree', 'note']).to_csv(
        CHK / 'STATS_CHECK_fig2.csv', index=False, lineterminator='\n')
    print('wrote data/fig2/: fig2d_meristem_cell_number.csv, fig2f_root_hair_length.csv, fig2g_root_length_ralf23.csv; '
          'output/stats_checks/STATS_CHECK_fig2.csv')
    if MISSING:
        print(f'note: label files not found ({", ".join(sorted(set(MISSING)))}); computed labels written with label_drawn empty')
    for c in CHECK:
        if c['agree'] == 'no':
            print('MISMATCH', c['panel'], c['comparison'], c['label_drawn'], '->', c['label_computed'])


if __name__ == '__main__':
    main()
