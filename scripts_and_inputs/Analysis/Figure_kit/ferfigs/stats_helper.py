"""OPTIONAL helper that PROPOSES a statistics table for author review. It is never drawn automatically.

    python build.py stats-propose data/fig4/fig4g_yuc5_reporter_tsko.csv --test twoway_tukey

writes data/fig4/fig4g_yuc5_reporter_tsko_stats_PROPOSED.csv (+ _PROPOSED.txt with the ANOVA table). The build
only ever reads the file named in the figure config (..._stats.csv). To adopt a proposal, an author checks it
against their own analysis, copies the labels into ..._stats.csv and sets source = computed_by_authors.

Tests offered (pick the one that matches the experimental design and the Methods):
  tukey_all       one-way ANOVA across every bar, Tukey HSD on all pairs, compact letter display
  twoway_tukey    two-way ANOVA (group x series, with interaction), Tukey HSD on all bars, compact letters
  tukey_within    within each group, one-way ANOVA and Tukey HSD among the series; letters restart for every
                  group (written to ..._stats_within_PROPOSED.csv)
  vs_control      within each group, each series against the first (control) series: Dunnett's test
                  (scipy), asterisks (* P<0.05, ** <0.01, *** <0.001, **** <0.0001, ns)
Values recovered from figure drawings are not original measurements: a proposal computed from them is a
plausibility check only.
"""
from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd

from .schemas import read_table


def cld(names, means, sig_pairs):
    """Compact letter display (insert-and-absorb, Piepho 2004). names in any order; letters assigned so that 'a'
    goes to the group containing the largest mean. sig_pairs: set of frozenset({i, j}) that differ."""
    cols = [set(names)]
    for pair in sig_pairs:
        i, j = tuple(pair)
        new = []
        for c in cols:
            if i in c and j in c:
                new.append(c - {i}); new.append(c - {j})
            else:
                new.append(c)
        # absorb
        new = [c for c in new if c]
        keep = []
        for c in sorted(new, key=len, reverse=True):
            if not any(c <= k for k in keep):
                keep.append(c)
        cols = keep
    m = dict(zip(names, means))
    cols.sort(key=lambda c: -max(m[n] for n in c))
    letters = {n: '' for n in names}
    for k, c in enumerate(cols):
        for n in c:
            letters[n] += chr(ord('a') + k)
    return letters


def stars(p):
    return '****' if p < 1e-4 else '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 0.05 else 'ns'


def propose(csv_path, test='twoway_tukey', control=None, alpha=0.05):
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    df = read_table(csv_path, 'measurements')
    df['cell'] = df['group'] + ' | ' + df['series']
    order = list(dict.fromkeys(df['cell']))
    lines = [f'# PROPOSED statistics for {Path(csv_path).name} (test: {test}); not for drawing without author review']
    rows = []
    if test in ('tukey_all', 'twoway_tukey'):
        if test == 'twoway_tukey':
            import statsmodels.formula.api as smf
            from statsmodels.stats.anova import anova_lm
            mdl = smf.ols('value ~ C(group) * C(series)', data=df).fit()
            lines.append(anova_lm(mdl, typ=2).to_string())
        else:
            from scipy.stats import f_oneway
            F, p = f_oneway(*[df[df['cell'] == c]['value'] for c in order])
            lines.append(f'one-way ANOVA F = {F:.4g}, P = {p:.3g}')
        tk = pairwise_tukeyhsd(df['value'], df['cell'], alpha=alpha)
        res = pd.DataFrame(tk.summary().data[1:], columns=tk.summary().data[0])
        lines.append(res.to_string())
        sig = {frozenset((r['group1'], r['group2'])) for _, r in res.iterrows() if bool(r['reject'])}
        means = [df[df['cell'] == c]['value'].mean() for c in order]
        lt = cld(order, means, sig)
        label = ('two-way ANOVA with Tukey\'s multiple comparisons' if test == 'twoway_tukey'
                 else 'one-way ANOVA with Tukey\'s multiple comparisons')
        for c in order:
            g, s = c.split(' | ')
            rows.append(dict(group=g, series=s, label=lt[c], kind='letters', group2='', series2='', test=label,
                             p_value='', source='proposed_by_helper', notes=f'alpha = {alpha}'))
    elif test == 'tukey_within':
        from scipy.stats import f_oneway
        label = "one-way ANOVA with Tukey's multiple comparisons within each genotype"
        for g in dict.fromkeys(df['group']):
            sub = df[df['group'] == g]
            ser = list(dict.fromkeys(sub['series']))
            F, p = f_oneway(*[sub[sub['series'] == s]['value'] for s in ser])
            lines.append(f'{g}: one-way ANOVA F = {F:.4g}, P = {p:.3g}')
            tk = pairwise_tukeyhsd(sub['value'], sub['series'], alpha=alpha)
            res = pd.DataFrame(tk.summary().data[1:], columns=tk.summary().data[0])
            lines.append(res.to_string())
            sig = {frozenset((r['group1'], r['group2'])) for _, r in res.iterrows() if bool(r['reject'])}
            lt = cld(ser, [sub[sub['series'] == s]['value'].mean() for s in ser], sig)
            for s in ser:
                rows.append(dict(group=g, series=s, label=lt[s], kind='letters', group2='', series2='', test=label,
                                 p_value='', source='proposed_by_helper',
                                 notes=f'alpha = {alpha}; letters compare treatments within {g} only'))
    elif test == 'vs_control':
        from scipy.stats import dunnett
        for g in dict.fromkeys(df['group']):
            sub = df[df['group'] == g]
            ser = list(dict.fromkeys(sub['series']))
            ctl = control or ser[0]
            others = [s for s in ser if s != ctl]
            if not others:
                continue
            r = dunnett(*[sub[sub['series'] == s]['value'].to_numpy() for s in others],
                        control=sub[sub['series'] == ctl]['value'].to_numpy())
            for s, p in zip(others, r.pvalue):
                rows.append(dict(group=g, series=ctl, label=stars(p), kind='ns' if p >= 0.05 else 'stars', group2=g,
                                 series2=s, test=f"Dunnett's test against {ctl} within each genotype", p_value=f'{p:.3g}',
                                 source='proposed_by_helper', notes=''))
                lines.append(f'{g}: {s} vs {ctl}: P = {p:.3g}')
    else:
        raise ValueError(f'unknown test {test}')
    out = pd.DataFrame(rows)
    stem = Path(csv_path).with_suffix('')
    tag = '_within' if test == 'tukey_within' else ''          # both letter schemes can sit side by side
    out_csv = Path(f'{stem}_stats{tag}_PROPOSED.csv')
    out.to_csv(out_csv, index=False, lineterminator='\n')
    Path(f'{stem}_stats{tag}_PROPOSED.txt').write_text('\n'.join(lines) + '\n')
    return out_csv, out
