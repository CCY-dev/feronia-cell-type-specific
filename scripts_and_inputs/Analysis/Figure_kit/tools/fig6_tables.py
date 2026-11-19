#!/usr/bin/env python
"""Fig. 6 measurement tables and statistics check.

Reads the tidy source data in source_data/long/ (written by tools/source_data_to_long.py from the Source data
workbook) and writes the kit's measurement schema into data/fig6/:
  a  fig6a_root_length.csv     sheet 'Fig_6a RALF_dof'       (6 genotypes x Mock / RALF23)
  d  fig6d_pyuc5_dof.csv       sheet 'Fig_6d pYUC5_dof 6h'   (Col-0, cdf4, dag1dof4.6 x Mock / RALF23)
  e  fig6e_qpcr_<gene>.csv     sheet 'Fig_6e DOFs qPCR'      (one table per gene: DOF4.6, CDF4, DOF2.2, DAG1)
  g  fig6g_<reporter>.csv      sheet 'Fig_6g DOFs-GFP'       (pDAG1:erGFP, pDOF2.2:erGFP)

Asserts: n and mean per group x series equal source_data/long/summary.csv (Fig6e: summary.csv pools the four genes,
so the pooled n / mean are asserted).

Statistics: the labels drawn are stored in data/fig6/*_stats.csv and the printed ratios of panel b in
data/fig6/fig6b_seedlings.csv (read only). The tests named in the legends are recomputed from the source data and
compared in output/stats_checks/STATS_CHECK_fig6.csv:
  a  two-way ANOVA (genotype x treatment), RALF23 vs Mock within each genotype on the pooled MSE, Sidak over the
     6 genotypes (Student's t test reported as well)
  b  printed ratio vs mean(RALF23) / mean(Mock) of the panel a table
  d, e  two-way ANOVA + Tukey HSD on all cell means (compact letters, 'a' on the largest mean); 'same grouping'
     when the letter names differ but the letters encode the same pairs
  g  Student's t test (unpaired, two-tailed; Welch reported as well)
If the label files are not present, the computed labels are written with label_drawn empty.

Usage (from the kit folder, after `python tools/source_data_to_long.py "<Source data.xlsx>"`):
  python tools/fig6_tables.py
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
OUT = KIT / 'data' / 'fig6'
CHK = KIT / 'output' / 'stats_checks'
NO_LABELS = 'labels not in this repository'
COLS = ['group', 'series', 'replicate', 'value', 'sample_id', 'experiment', 'unit', 'notes']
MISSING = []

summary = pd.read_csv(SRC / 'summary.csv')
summary.columns = [c.strip() for c in summary.columns]
checks = []          # rows of STATS_CHECK_fig6.csv


def load(sheet_csv):
    d = pd.read_csv(SRC / sheet_csv, dtype={'replicate': str, 'sample_id': str, 'notes': str}, keep_default_na=False)
    d['value'] = d['value'].astype(float)
    return d


def assert_summary(df, panel):
    """n and mean per group x series (pooled over genes for Fig6e) equal summary.csv."""
    s = summary[summary['panel'] == panel]
    assert len(s), panel
    got = df.groupby(['group', 'series'])['value'].agg(['count', 'mean']).reset_index()
    m = got.merge(s, on=['group', 'series'], how='outer', indicator=True)
    assert (m['_merge'] == 'both').all(), m
    assert (m['count'] == m['n']).all(), m
    assert np.allclose(m['mean_x'], m['mean_y'], atol=1e-4), m


def tidy(d, sheet, unit, group_col='group'):
    out = pd.DataFrame(dict(group=d[group_col], series=d['series'], replicate=d['replicate'].astype(str),
                            value=d['value'], sample_id=sheet + '!' + d['cell'], experiment='',
                            unit=unit, notes=d['notes']))
    assert not out.duplicated(['group', 'series', 'replicate']).any()
    return out[COLS]


def write(df, name):
    df.to_csv(OUT / name, index=False, lineterminator='\n')


def labels(name):
    """Authors' labels {(group, series): label} from data/fig6/<name>; None when not in this copy of the kit."""
    p = OUT / name
    if not p.exists():
        MISSING.append(name)
        return None
    st = pd.read_csv(p, dtype=str, keep_default_na=False)
    return {(r['group'], r['series']): r['label'] for _, r in st.iterrows()}


def tukey_letters(df):
    """Two-way ANOVA (group x series with interaction) + Tukey on all cell means (= pooled MSE of the full model,
    Tukey-Kramer), compact letters with 'a' on the largest mean."""
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    cell = df['group'] + ' | ' + df['series']
    tk = pairwise_tukeyhsd(df['value'], cell, alpha=0.05)
    res = pd.DataFrame(tk.summary().data[1:], columns=tk.summary().data[0])
    order = list(dict.fromkeys(cell))
    sig = {frozenset((r['group1'], r['group2'])) for _, r in res.iterrows() if bool(r['reject'])}
    lt = cld(order, [df[cell == c]['value'].mean() for c in order], sig)
    return {tuple(c.split(' | ')): lt[c] for c in order}


def same_partition(author, computed):
    """Letters agree when they encode the same 'shares a letter' relation (letter names may differ)."""
    keys = list(author)
    for i in keys:
        for j in keys:
            if bool(set(author[i]) & set(author[j])) != bool(set(computed[i]) & set(computed[j])):
                return False
    return True


def letter_checks(panel, df, drawn, note):
    comp = tukey_letters(df)
    keys = list(drawn) if drawn is not None else list(comp)
    same = drawn is not None and same_partition({k: drawn[k] for k in keys}, comp)
    ns_ = df.groupby(['group', 'series']).size().to_dict()
    for k in keys:
        lab = drawn[k] if drawn is not None else ''
        checks.append(dict(panel=panel, comparison=f'{k[0]} {k[1]}', label_drawn=lab, label_computed=comp[k],
                           agree=NO_LABELS if drawn is None else
                           ('yes' if comp[k] == lab else ('same grouping' if same else 'no')),
                           note=note.format(n=ns_[k])))


def main():
    import statsmodels.formula.api as smf
    from scipy import stats as ss
    OUT.mkdir(parents=True, exist_ok=True)
    CHK.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------- a: root length (two-way ANOVA + Sidak)
    sheet = 'Fig_6a RALF_dof'
    d = load('Fig_6a_RALF_dof.csv')
    assert_summary(d, 'Fig6a')
    GEN_A = ['Col-0', 'fer-4', 'dof4.6', 'cdf4', 'cog1', 'obp3']
    a = tidy(d, sheet, 'mm')
    write(a, 'fig6a_root_length.csv')
    drawn = labels('fig6a_root_length_stats.csv')
    mdl = smf.ols('value ~ C(group) * C(series)', data=a).fit()
    mse, dfr = mdl.mse_resid, mdl.df_resid
    for g in GEN_A:
        m = a[(a.group == g) & (a.series == 'Mock')]['value']
        r = a[(a.group == g) & (a.series == 'RALF23')]['value']
        t = (r.mean() - m.mean()) / np.sqrt(mse * (1 / len(m) + 1 / len(r)))
        p_raw = 2 * ss.t.sf(abs(t), dfr)
        p_sid = 1 - (1 - p_raw) ** len(GEN_A)
        p_t = ss.ttest_ind(m, r, equal_var=True).pvalue
        lab_c = stars(p_sid)
        lab = drawn[(g, 'RALF23')] if drawn is not None else ''
        checks.append(dict(panel='a', comparison=f'{g}: Mock vs RALF23', label_drawn=lab, label_computed=lab_c,
                           agree=NO_LABELS if drawn is None else ('yes' if lab_c == lab else 'no'),
                           note=f"two-way ANOVA (pooled MSE, df {dfr:.0f}) + Sidak over 6 comparisons: P = "
                                f"{(f'{p_sid:.3g}' if p_sid > 0 else '< 1e-15')}; "
                                f"Student's t test (legend also names a t test): P = {p_t:.3g} ({stars(p_t)}); "
                                f"n Mock/RALF23 = {len(m)}/{len(r)} (legend n = 14-17)"))

    # ---------------------------------------------------------------- b: printed ratios vs data a
    pb = OUT / 'fig6b_seedlings.csv'
    printed = None
    if pb.exists():
        printed = pd.read_csv(pb, dtype=str, keep_default_na=False).set_index('tile')['ratio'].to_dict()
    else:
        MISSING.append(pb.name)
    for g in GEN_A:
        mm_ = a[(a.group == g) & (a.series == 'Mock')]['value'].mean()
        rr = a[(a.group == g) & (a.series == 'RALF23')]['value'].mean()
        ratio = rr / mm_
        comp = f'{ratio:.2f}'
        lab = printed[g] if printed is not None else ''
        checks.append(dict(panel='b', comparison=f'{g}: printed ratio RALF23/Mock', label_drawn=lab,
                           label_computed=comp, agree=NO_LABELS if printed is None else ('yes' if comp == lab else 'no'),
                           note=f'mean(RALF23) {rr:.4f} / mean(Mock) {mm_:.4f} = {ratio:.4f} (sheet {sheet}); printed '
                                f'value kept on the image'))

    # ---------------------------------------------------------------- d: pYUC5 in dof mutants (two-way + Tukey)
    sheet = 'Fig_6d pYUC5_dof 6h'
    d = load('Fig_6d_pYUC5_dof_6h.csv')
    assert_summary(d, 'Fig6d')
    dd = tidy(d, sheet, 'a.u.')
    write(dd, 'fig6d_pyuc5_dof.csv')
    letter_checks('d', dd, labels('fig6d_pyuc5_dof_stats.csv'),
                  'two-way ANOVA + Tukey on 6 cell means, alpha 0.05; n = {n} '
                  '(legend n = 5-15; source Col-0 n = 20 per treatment)')

    # ---------------------------------------------------------------- e: DOF qPCR, one chart per gene
    sheet = 'Fig_6e DOFs qPCR'
    d = load('Fig_6e_DOFs_qPCR.csv')
    assert_summary(d, 'Fig6e')                     # summary.csv pools the four genes
    GENES = ['DOF4.6', 'CDF4', 'DOF2.2', 'DAG1']
    assert set(d['sample_id']) == set(GENES)
    for gene in GENES:
        sub = d[d['sample_id'] == gene].copy()
        e = tidy(sub, sheet, 'relative expression')
        e['experiment'] = gene
        fn = f"fig6e_qpcr_{gene.lower().replace('.', '')}.csv"
        write(e, fn)
        letter_checks(f'e {gene}', e, labels(fn.replace('.csv', '_stats.csv')),
                      'two-way ANOVA + Tukey on 4 cell means, alpha 0.05; n = {n} '
                      '(legend n = 7-9 from 3 biological replicates)')

    # ---------------------------------------------------------------- g: DOF reporters (Student's t test)
    sheet = 'Fig_6g DOFs-GFP'
    d = load('Fig_6g_DOFs_GFP.csv')
    assert_summary(d, 'Fig6g')
    for rep in ('pDAG1:erGFP', 'pDOF2.2:erGFP'):
        sub = d[d['group'] == rep]
        g_ = tidy(sub, sheet, 'a.u.')
        fn = f"fig6g_{rep.split(':')[0].lower().replace('.', '')}.csv"
        write(g_, fn)
        drawn = labels(fn.replace('.csv', '_stats.csv'))
        m = g_[g_.series == 'Mock']['value']
        r = g_[g_.series == 'RALF23']['value']
        p = ss.ttest_ind(m, r, equal_var=True).pvalue
        pw = ss.ttest_ind(m, r, equal_var=False).pvalue
        lab = drawn[(rep, 'Mock')] if drawn is not None else ''
        checks.append(dict(panel='g', comparison=f'{rep}: Mock vs RALF23', label_drawn=lab, label_computed=stars(p),
                           agree=NO_LABELS if drawn is None else ('yes' if stars(p) == lab else 'no'),
                           note=f"Student's t test (equal variances, two-tailed): P = {p:.3g}; Welch: P = {pw:.3g} "
                                f"({stars(pw)}); n = {len(m)}/{len(r)} (legend n = 6)"))

    ck = pd.DataFrame(checks)[['panel', 'comparison', 'label_drawn', 'label_computed', 'agree', 'note']]
    ck.to_csv(CHK / 'STATS_CHECK_fig6.csv', index=False, lineterminator='\n')
    print(ck[['panel', 'comparison', 'label_drawn', 'label_computed', 'agree']].to_string())
    if MISSING:
        print(f'note: label files not found ({", ".join(MISSING)}); computed labels written with label_drawn empty')


if __name__ == '__main__':
    main()
