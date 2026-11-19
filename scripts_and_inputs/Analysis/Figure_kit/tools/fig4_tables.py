#!/usr/bin/env python
"""Fig. 4 (DR5rev::GFP, PEO-IAA, afb1 / tir1-1/afb2-3) measurement tables and statistics check.

Reads the tidy source data in source_data/long/ (ALL_long.csv, summary.csv; written by tools/source_data_to_long.py
from the Source data workbook) and writes data/fig4/:

  fig4b_dr5_col0_fer4.csv              panel Fig4b (DR5rev::GFP, Col-0 / fer-4 x Mock / RALF23)
  fig4d_dr5_tsko.csv                   panel Fig4d (DR5rev::GFP, tissue-specific knockouts)
  fig4e_root_length_peo_iaa_L3.csv     panel Fig4e (root length, Mock / RALF23 / RALF23 & PEO-IAA, _L3 lines)
  fig4g_root_length_afb1_tir1.csv      panel Fig4g (root length, Col-0 / afb1 / tir1-1/afb2-3 x Mock / RALF23)

Every n and mean per group x series is asserted equal to summary.csv.

Statistics: the labels drawn are the authors' labels, stored in data/fig4/*_stats.csv, fig4e_ratios.csv and
fig4f_root_length_ratios.csv (read only). They are recomputed with the tests named in the legend:
  b, d  two-way ANOVA + Tukey HSD on all bars (letters, compared per bar and as a significance pattern)
  e     two-way ANOVA + Dunnett within each genotype (Mock and RALF23 & PEO-IAA vs RALF23, pooled MSE);
        Sidak over all 12 comparisons as a sensitivity check
  g     two-way ANOVA + Sidak (Mock vs RALF23 within genotype, pooled MSE); Student and Welch t tests reported
  e, f  printed ratios vs mean(treated) / mean(Mock) (f uses the panel g table)
Output: output/stats_checks/STATS_CHECK_fig4.csv and STATS_CHECK_fig4_anova.txt.
If the label files are not present, the computed labels are written with label_drawn empty.

Usage (from the kit folder, after `python tools/source_data_to_long.py "<Source data.xlsx>"`):
  python tools/fig4_tables.py
"""
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

KIT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(KIT))
from ferfigs.stats_helper import cld, stars  # noqa: E402

SRC = KIT / 'source_data' / 'long'
OUT = KIT / 'data' / 'fig4'
CHK = KIT / 'output' / 'stats_checks'
NO_LABELS = 'labels not in this repository'
MISSING = []

TABLES = {  # Fig 4 panel: (file stem, source panel, unit)
    'b': ('fig4b_dr5_col0_fer4', 'Fig4b', 'a.u.'),
    'd': ('fig4d_dr5_tsko', 'Fig4d', 'a.u.'),
    'e': ('fig4e_root_length_peo_iaa_L3', 'Fig4e', 'mm'),
    'g': ('fig4g_root_length_afb1_tir1', 'Fig4g', 'mm'),
}


def write_csv(df, path):
    tmp = path.with_suffix('.tmp')
    df.to_csv(tmp, index=False, lineterminator='\n')
    tmp.replace(path)


def labels(name):
    """Authors' labels {(group, series): label} from data/fig4/<name>; None when not in this copy of the kit."""
    p = OUT / name
    if not p.exists():
        MISSING.append(name)
        return None
    st = pd.read_csv(p, dtype=str, keep_default_na=False)
    return {(r['group'], r['series']): r['label'] for _, r in st.iterrows()}


def agree(drawn, key, computed):
    if drawn is None:
        return '', NO_LABELS
    return drawn[key], ('yes' if drawn[key] == computed else 'no')


def measurements():
    src = pd.read_csv(SRC / 'ALL_long.csv')
    src = src[src['value'].notna()]
    summ = pd.read_csv(SRC / 'summary.csv')
    out = {}
    for p, (stem, panel, unit) in TABLES.items():
        s = src[src['panel'] == panel].copy()
        assert len(s), panel
        assert set(s['unit']) == {unit}, (panel, set(s['unit']))
        df = pd.DataFrame(dict(group=s['group'], series=s['series'], replicate=s['replicate'].astype(int),
                               value=s['value'].astype(float), unit=unit,
                               sample_id=s['sheet'] + '!' + s['cell'], experiment='', notes=''))
        assert not df.duplicated(['group', 'series', 'replicate']).any(), panel
        sm = summ[summ['panel'] == panel]
        for _, r in sm.iterrows():
            v = df[(df['group'] == r['group']) & (df['series'] == r['series'])]['value']
            assert len(v) == int(r['n']), (panel, r['group'], len(v), r['n'])
            assert abs(v.mean() - float(r['mean'])) < 1e-4, (panel, r['group'], v.mean(), r['mean'])
        assert len(sm) == df.groupby(['group', 'series']).ngroups
        write_csv(df, OUT / f'{stem}.csv')
        out[p] = df
        print(f'{stem}: ' + ', '.join(f'{g}/{s} n={n}' for (g, s), n in
                                      df.groupby(['group', 'series'], sort=False).size().items()))
    return out


def anova_fit(df):
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm
    mdl = smf.ols('value ~ C(group) * C(series)', data=df).fit()
    return mdl, anova_lm(mdl, typ=2)


def dunnett_p(t_abs, ns, n_ctl, dfe):
    """Two-sided Dunnett-adjusted P for comparisons sharing one control (pooled MSE, df = dfe)."""
    from scipy.stats import multivariate_t
    k = len(ns)
    R = np.eye(k)
    for i in range(k):
        for j in range(k):
            if i != j:
                R[i, j] = 1 / math.sqrt((1 + n_ctl / ns[i]) * (1 + n_ctl / ns[j]))
    mvt = multivariate_t(loc=np.zeros(k), shape=R, df=dfe, seed=1)
    inside = mvt.cdf(np.full(k, t_abs), lower_limit=np.full(k, -t_abs), maxpts=200000)
    return float(min(1.0, max(0.0, 1 - inside)))


def letters_check(p, df, drawn, rows, txt):
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    df = df.copy()
    df['cell'] = df['group'] + ' | ' + df['series']
    mdl, aov = anova_fit(df)
    txt += [f'== Fig. 4{p}: two-way ANOVA (type II), value ~ group * series', aov.to_string()]
    tk = pairwise_tukeyhsd(df['value'], df['cell'], alpha=0.05)
    res = pd.DataFrame(tk.summary().data[1:], columns=tk.summary().data[0])
    txt += ['Tukey HSD on all bars', res.to_string(), '']
    order = list(dict.fromkeys(df['cell']))
    sig = {frozenset((r['group1'], r['group2'])) for _, r in res.iterrows() if bool(r['reject'])}
    lt = cld(order, [df[df['cell'] == c]['value'].mean() for c in order], sig)
    for c in order:
        g, s = c.split(' | ')
        lab, ag = agree(drawn, (g, s), lt[c])
        rows.append(dict(panel=f'4{p}', comparison=f'{g} {s} (letter vs all bars)', label_drawn=lab,
                         label_computed=lt[c], agree=ag,
                         note="authors' letter kept; computed = two-way ANOVA + Tukey HSD on all bars from the source data"))
    if drawn is None:
        rows.append(dict(panel=f'4{p}', comparison='all pairs (significance pattern)', label_drawn='',
                         label_computed='', agree=NO_LABELS, note=''))
        return
    dis = []
    for i, c1 in enumerate(order):
        for c2 in order[i + 1:]:
            k1, k2 = tuple(c1.split(' | ')), tuple(c2.split(' | '))
            d_diff = not (set(drawn[k1]) & set(drawn[k2]))
            c_diff = frozenset((c1, c2)) in sig
            if d_diff != c_diff:
                pr = res[((res['group1'] == c1) & (res['group2'] == c2)) | ((res['group1'] == c2) & (res['group2'] == c1))]
                dis.append(f"{c1} vs {c2}: drawn {'differ' if d_diff else 'same'}, computed "
                           f"{'differ' if c_diff else 'same'} (Tukey p-adj {float(pr['p-adj'].iloc[0]):.3g})")
    rows.append(dict(panel=f'4{p}', comparison='all pairs (significance pattern)', label_drawn='', label_computed='',
                     agree='yes' if not dis else 'no',
                     note='; '.join(dis) if dis else 'drawn and computed letters imply the same significant pairs'))


def stats(data):
    from scipy.stats import t as tdist, ttest_ind
    rows, txt = [], []
    # ---- b, d: letters
    for p in ('b', 'd'):
        letters_check(p, data[p], labels(f'{TABLES[p][0]}_stats.csv'), rows, txt)
    # ---- e: two-way ANOVA + Dunnett within genotype, control = RALF23
    drawn = labels(f'{TABLES["e"][0]}_stats.csv')
    df = data['e']
    mdl, aov = anova_fit(df)
    mse, dfe = mdl.mse_resid, mdl.df_resid
    txt += ['== Fig. 4e: two-way ANOVA (type II), root length ~ genotype * treatment', aov.to_string(),
            f'pooled MSE = {mse:.5g}, df = {dfe:g}',
            "Dunnett within each genotype (family = the 2 comparisons of one genotype), control = RALF23; P below 1e-6 "
            "reported as the Bonferroni bound 2 x P(unadj) (the multivariate-t integral has no precision there); "
            "sensitivity: one family of all 12 comparisons"]
    allp = []
    for g in list(dict.fromkeys(df['group'])):
        sub = df[df['group'] == g]
        c = sub[sub['series'] == 'RALF23']['value']
        ns, res = [], []
        for s in ('Mock', 'RALF23 & PEO-IAA'):
            v = sub[sub['series'] == s]['value']
            t = (v.mean() - c.mean()) / math.sqrt(mse * (1 / len(v) + 1 / len(c)))
            ns.append(len(v)); res.append((s, v, t))
        for s, v, t in res:
            pu = 2 * tdist.sf(abs(t), dfe)
            pd_ = dunnett_p(abs(t), ns, len(c), dfe)
            if pd_ < 1e-6:            # numerical cdf floor: report the Bonferroni upper bound instead
                pd_ = min(1.0, 2 * pu)
            allp.append((g, s, abs(t)))
            lab = stars(pd_)
            au, ag = agree(drawn, (g, s), lab)
            txt.append(f'{g}: {s} vs RALF23: diff = {v.mean() - c.mean():.3f}, t = {t:.3f}, P(unadj) = {pu:.3g}, '
                       f'P(Dunnett, 2 per family) = {pd_:.3g} -> {lab}')
            rows.append(dict(panel='4e', comparison=f'{g}: {s} vs RALF23', label_drawn=au, label_computed=lab,
                             agree=ag, note=f"computed = two-way ANOVA + Dunnett (family = genotype, pooled MSE) "
                             f"P = {pd_:.3g}; Sidak over all 12 = {stars(-math.expm1(12 * math.log1p(-pu)))}"))
    txt.append('sensitivity (all 12 comparisons as one family, Sidak on unadjusted P):')
    for g, s, t in allp:
        pu = 2 * tdist.sf(t, dfe)
        ps12 = -math.expm1(12 * math.log1p(-pu))
        txt.append(f'  {g}: {s} vs RALF23: Sidak12 P = {ps12:.3g} -> {stars(ps12)}')
    txt.append('')
    # ---- g: two-way ANOVA + Sidak (Mock vs RALF23 within genotype); Student / Welch t tests reported
    drawn = labels(f'{TABLES["g"][0]}_stats.csv')
    df = data['g']
    mdl, aov = anova_fit(df)
    mse, dfe = mdl.mse_resid, mdl.df_resid
    txt += ['== Fig. 4g: two-way ANOVA (type II), root length ~ genotype * treatment', aov.to_string(),
            f'pooled MSE = {mse:.5g}, df = {dfe:g}; Sidak over the 3 mock-vs-RALF23 comparisons']
    gg = list(dict.fromkeys(df['group']))
    for g in gg:
        a = df[(df['group'] == g) & (df['series'] == 'Mock')]['value']
        b = df[(df['group'] == g) & (df['series'] == 'RALF23')]['value']
        t = (a.mean() - b.mean()) / math.sqrt(mse * (1 / len(a) + 1 / len(b)))
        pu = 2 * tdist.sf(abs(t), dfe)
        ps = -math.expm1(len(gg) * math.log1p(-pu))
        tw = ttest_ind(a, b, equal_var=False).pvalue
        ts = ttest_ind(a, b, equal_var=True).pvalue
        lab = stars(ps)
        au, ag = agree(drawn, (g, 'Mock'), lab)
        txt.append(f'{g}: Mock vs RALF23: t = {t:.3f}, P(Sidak, 3) = {ps:.3g} -> {lab}; Student t P = {ts:.3g}; '
                   f'Welch t P = {tw:.3g} -> {stars(ts)} / {stars(tw)}')
        rows.append(dict(panel='4g', comparison=f'{g}: Mock vs RALF23', label_drawn=au, label_computed=lab,
                         agree=ag, note=f"two-way ANOVA + Sidak P = {ps:.3g}; Student t P = {ts:.3g} "
                         f"({stars(ts)}); Welch t P = {tw:.3g} ({stars(tw)}); legend names both tests (conflict)"))
    return rows, txt


def ratios(data, rows):
    for p, fname in (('e', 'fig4e_ratios.csv'), ('f', 'fig4f_root_length_ratios.csv')):
        df = data['e' if p == 'e' else 'g']
        m = df.groupby(['group', 'series'])['value'].mean()
        printed = labels(fname)
        if printed is None:     # no printed values: report every treated/Mock ratio at 2 decimals
            keys = [(g, s) for (g, s) in m.index if s != 'Mock' and (p == 'e' or s == 'RALF23')]
            printed_ = {k: None for k in keys}
        else:
            printed_ = printed
        for (g, s), lab in printed_.items():
            v = m[(g, s)] / m[(g, 'Mock')]
            dec = len(lab.split('.')[1]) if lab else 2
            rows.append(dict(panel=f'4{p}', comparison=f'ratio {g} {s}/Mock', label_drawn=lab or '',
                             label_computed=f'{v:.{dec}f} ({v:.4f})',
                             agree=NO_LABELS if lab is None else ('yes' if f'{v:.{dec}f}' == lab else 'no'),
                             note='printed value kept; computed = mean(treated)/mean(mock) of the source data'
                                  + ('' if p == 'e' else ' (panel g table)')))


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    CHK.mkdir(parents=True, exist_ok=True)
    data = measurements()
    rows, txt = stats(data)
    ratios(data, rows)
    write_csv(pd.DataFrame(rows), CHK / 'STATS_CHECK_fig4.csv')
    (CHK / 'STATS_CHECK_fig4_anova.txt').write_text('\n'.join(txt) + '\n')
    if MISSING:
        print(f'note: label files not found ({", ".join(MISSING)}); computed labels written with label_drawn empty')
    for r in rows:
        if r['agree'] == 'no':
            print('MISMATCH', r)
