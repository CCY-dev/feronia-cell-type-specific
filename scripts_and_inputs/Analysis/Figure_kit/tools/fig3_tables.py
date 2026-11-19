#!/usr/bin/env python
"""Fig. 3 (pYUC5:GFP-GUS reporter) measurement tables and statistics check.

Reads the tidy source data in source_data/long/ (ALL_long.csv, summary.csv; written by tools/source_data_to_long.py
from the Source data workbook) and writes data/fig3/:

  fig3b_yuc5_reporter_timecourse.csv   panel Fig3b (time course; groups '1h' '3h' '6h' renamed '1 h' '3 h' '6 h')
  fig3d_yuc5_reporter_col0_fer4.csv    panel Fig3d (Col-0 / fer-4 x Mock / RALF23)
  fig3g_yuc5_reporter_tsko.csv         panel Fig3g (tissue-specific knockouts x Mock / RALF23)

Every n and mean per group x series is asserted equal to summary.csv.

Statistics: the letters drawn are the authors' letters, stored in data/fig3/*_stats.csv (read only). They are checked
against a two-way ANOVA (group x series, with interaction) + Tukey HSD on all bars, as named in the legend; per bar
and as a whole significance pattern (letters are compared as a partition: does each pair share a letter?).
Output: output/stats_checks/STATS_CHECK_fig3.csv and STATS_CHECK_fig3_anova.txt (ANOVA tables + Tukey results).
If the label files are not present, the computed letters are written with label_drawn empty.

Usage (from the kit folder, after `python tools/source_data_to_long.py "<Source data.xlsx>"`):
  python tools/fig3_tables.py
"""
import sys
from pathlib import Path

import pandas as pd

KIT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(KIT))
from ferfigs.stats_helper import cld  # noqa: E402

SRC = KIT / 'source_data' / 'long'
OUT = KIT / 'data' / 'fig3'
CHK = KIT / 'output' / 'stats_checks'
NO_LABELS = 'labels not in this repository'

TABLES = {  # panel letter: (file stem, source panel, group renames)
    'b': ('fig3b_yuc5_reporter_timecourse', 'Fig3b', {'1h': '1 h', '3h': '3 h', '6h': '6 h'}),
    'd': ('fig3d_yuc5_reporter_col0_fer4', 'Fig3d', {}),
    'g': ('fig3g_yuc5_reporter_tsko', 'Fig3g', {}),
}


def write_csv(df, path):
    tmp = path.with_suffix('.tmp')
    df.to_csv(tmp, index=False, lineterminator='\n')
    tmp.replace(path)


def measurements():
    src = pd.read_csv(SRC / 'ALL_long.csv')
    src = src[src['value'].notna()]
    summ = pd.read_csv(SRC / 'summary.csv')
    out = {}
    for p, (stem, panel, ren) in TABLES.items():
        s = src[src['panel'] == panel].copy()
        assert len(s), panel
        s['group'] = s['group'].map(lambda g: ren.get(g, g))
        df = pd.DataFrame(dict(group=s['group'], series=s['series'], replicate=s['replicate'].astype(int),
                               value=s['value'].astype(float), unit='a.u.',
                               sample_id=s['sheet'] + '!' + s['cell'], experiment='', notes=''))
        assert not df.duplicated(['group', 'series', 'replicate']).any(), panel
        sm = summ[summ.iloc[:, 0] == panel]
        for _, r in sm.iterrows():
            g = ren.get(r.iloc[1], r.iloc[1])
            v = df[(df['group'] == g) & (df['series'] == r.iloc[2])]['value']
            assert len(v) == int(r.iloc[3]), (panel, g, len(v), r.iloc[3])
            assert abs(v.mean() - float(r.iloc[4])) < 1e-4, (panel, g, v.mean(), r.iloc[4])
        assert len(sm) == df.groupby(['group', 'series']).ngroups
        write_csv(df, OUT / f'{stem}.csv')
        out[p] = df
        print(f'{stem}: ' + ', '.join(f'{g}/{s} n={n}' for (g, s), n in
                                      df.groupby(['group', 'series'], sort=False).size().items()))
    return out


def stats(data):
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    rows, txt, missing = [], [], []
    for p, (stem, panel, _) in TABLES.items():
        lab_file = OUT / f'{stem}_stats.csv'
        if lab_file.exists():
            st = pd.read_csv(lab_file, dtype=str, keep_default_na=False)
            drawn = {(r['group'], r['series']): r['label'] for _, r in st.iterrows()}
        else:
            missing.append(lab_file.name)
            drawn = None
        df = data[p].copy()
        df['cell'] = df['group'] + ' | ' + df['series']
        mdl = smf.ols('value ~ C(group) * C(series)', data=df).fit()
        txt += [f'== Fig. 3{p} ({panel}): two-way ANOVA (type II), value ~ group * series', anova_lm(mdl, typ=2).to_string()]
        tk = pairwise_tukeyhsd(df['value'], df['cell'], alpha=0.05)
        res = pd.DataFrame(tk.summary().data[1:], columns=tk.summary().data[0])
        txt += ['Tukey HSD on all bars (MSE of the full factorial model)', res.to_string(), '']
        order = list(dict.fromkeys(df['cell']))
        sig = {frozenset((r['group1'], r['group2'])) for _, r in res.iterrows() if bool(r['reject'])}
        lt = cld(order, [df[df['cell'] == c]['value'].mean() for c in order], sig)
        for c in order:
            g, s = c.split(' | ')
            d_ = drawn[(g, s)] if drawn is not None else ''
            rows.append(dict(panel=f'3{p}', comparison=f'{g} {s} (letter vs all bars)', label_drawn=d_,
                             label_computed=lt[c],
                             agree=NO_LABELS if drawn is None else ('yes' if d_ == lt[c] else 'no'), note=''))
        if drawn is None:
            rows.append(dict(panel=f'3{p}', comparison='all pairs (significance pattern)', label_drawn='',
                             label_computed='', agree=NO_LABELS, note=''))
            continue
        dis = []
        for c1, c2 in [(a, b) for i, a in enumerate(order) for b in order[i + 1:]]:
            k1, k2 = tuple(c1.split(' | ')), tuple(c2.split(' | '))
            d_diff = not (set(drawn[k1]) & set(drawn[k2]))
            c_diff = frozenset((c1, c2)) in sig
            if d_diff != c_diff:
                pr = res[((res['group1'] == c1) & (res['group2'] == c2)) | ((res['group1'] == c2) & (res['group2'] == c1))]
                dis.append(f"{c1} vs {c2}: drawn {'differ' if d_diff else 'same'}, computed "
                           f"{'differ' if c_diff else 'same'} (Tukey p-adj {float(pr['p-adj'].iloc[0]):.3g})")
        rows.append(dict(panel=f'3{p}', comparison='all pairs (significance pattern)', label_drawn='',
                         label_computed='', agree='yes' if not dis else 'no',
                         note='; '.join(dis) if dis else 'drawn and computed letters imply the same significant pairs'))
    write_csv(pd.DataFrame(rows), CHK / 'STATS_CHECK_fig3.csv')
    (CHK / 'STATS_CHECK_fig3_anova.txt').write_text('\n'.join(txt) + '\n')
    if missing:
        print(f'note: label files not found ({", ".join(missing)}); computed letters written with label_drawn empty')
    for r in rows:
        if r['agree'] == 'no':
            print('STATS', r)


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    CHK.mkdir(parents=True, exist_ok=True)
    stats(measurements())
