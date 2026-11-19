#!/usr/bin/env python
"""Fig. 5d (pYUC5:LUC transactivation) measurement tables and statistics check.

Fig. 5a, b tables are written by tools/fig5ab_tables.py (CellOracle outputs); this script writes only panel d.

Reads source_data/long/Fig_5d_LUC_assay.csv and summary.csv (written by tools/source_data_to_long.py from the Source
data workbook, sheet 'Fig_5d LUC assay'; one row per well = Firefly/Renilla, plus one '<date>_summary' row per
construct x experiment carrying the workbook's 'Normalized' fold change) and the effector order of panel d in
config/fig5.yaml. Writes data/fig5/:

  fig5d_luc_wells.csv     every LUC well (Firefly/Renilla) in the kit measurement schema
  fig5d_luc.csv           fold change per DOF x replicate = mean F/R of the DOF wells / mean F/R of the GUS-GFP
                          control wells of the same experiment (Control = 1)
  fig5d_luc_rep_map.csv   which experiment (date) is Rep1 / Rep2 for each DOF, well counts, the computed fold change,
                          the workbook 'Normalized' value where the workbook prints one, and the value shown in the
                          heatmap (1 decimal)

Asserts: n and mean per construct x experiment == summary.csv; every DOF has exactly two experiments; the computed
fold change equals the workbook 'Normalized' value (1 decimal) wherever the workbook prints one.

Statistics: the heatmap shows fold change only and the legend names no test. The workbook prints a P category per
DOF x experiment; it is compared with a two-sided t test on per-well F/R (DOF vs GUS-GFP of the same experiment,
Student and Welch) in output/stats_checks/STATS_CHECK_fig5.csv. The workbook P categories are not carried by the
long tables, so they are read from the Source data workbook only when its path is given:
  python tools/fig5_tables.py ["<Source data.xlsx>"]
Without it, label_workbook is empty and agree = 'labels not in this repository'.

Usage (from the kit folder, after `python tools/source_data_to_long.py "<Source data.xlsx>"`):
  python tools/fig5_tables.py                          # tables + computed tests
  python tools/fig5_tables.py "<Source data.xlsx>"     # also compare with the workbook P categories
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

KIT = Path(__file__).resolve().parents[1]
SRC = KIT / 'source_data' / 'long'
OUT = KIT / 'data' / 'fig5'
CHK = KIT / 'output' / 'stats_checks'
NO_LABELS = 'labels not in this repository'

# workbook construct names -> names used in the paper (DOF6 = AtDof3.2 = AT3G45610)
NAME = {'GUS-GFP': 'Control', 'Dof4.6': 'DOF4.6', 'Dof3.2': 'DOF6', 'Dof2.2': 'DOF2.2', 'Dof4.1': 'DOF4.1'}
# replicate shown as Rep1 / Rep2 in the heatmap = experiment date (each DOF was assayed in two experiments)
REP = {'DOF4.6': ('2026-03-04', '2026-02-27'), 'OBP1': ('2025-10-24', '2025-10-25'),
       'OBP2': ('2026-03-04', '2026-02-27'), 'OBP3': ('2025-10-21', '2025-10-25'),
       'DAG1': ('2026-03-04', '2026-02-27'), 'COG1': ('2025-10-25', '2026-02-27'),
       'CDF4': ('2026-03-12', '2026-04-08'), 'DOF2.2': ('2026-04-08', '2026-03-04'),
       'DOF6': ('2026-03-04', '2026-02-27'), 'DOF4.1': ('2026-04-08', '2025-11-27')}


def name(x):
    return NAME.get(x, x)


def cat(p):
    return 'ns' if p >= 0.05 else '*' if p >= 0.01 else '**' if p >= 0.001 else '***' if p >= 1e-4 else '****'


def workbook_p_labels(xlsx):
    """{(experiment date, effector): 'P < 0.05, *'} from the row under 'Normalized' of each experiment block."""
    import openpyxl
    ws = openpyxl.load_workbook(xlsx, data_only=True)['Fig_5d LUC assay']
    plab, names, date = {}, {}, None
    for row in ws.iter_rows():
        for c in row:
            if c.column == 1 and c.value is not None and hasattr(c.value, 'strftime'):
                date = c.value.strftime('%Y-%m-%d')
            elif c.value is not None and isinstance(c.value, str) and c.row > 1 and date:
                v = c.value.strip()
                if v.startswith('P '):
                    plab[(date, c.column)] = v
                elif v not in ('Firefly LUC', 'Renilla LUC', 'Firefly/Renilla', 'Normalized'):
                    names.setdefault(date, {})[c.column + 2] = name(v.replace('GUSGFP', 'GUS-GFP'))
    return {(d, e): plab[(d, col)] for d, nm in names.items() for col, e in nm.items() if (d, col) in plab}


def main(xlsx=None):
    OUT.mkdir(parents=True, exist_ok=True)
    CHK.mkdir(parents=True, exist_ok=True)
    effs = yaml.safe_load(open(KIT / 'config' / 'fig5.yaml'))['panels']['d']['effectors']
    assert effs[0] == 'Control' and set(effs[1:]) == set(REP), effs

    luc = pd.read_csv(SRC / 'Fig_5d_LUC_assay.csv')
    summ = luc[luc.replicate.str.endswith('_summary')].copy()
    wells = luc[~luc.replicate.str.endswith('_summary')].copy()
    s = pd.read_csv(SRC / 'summary.csv')
    s = s[s.panel == 'Fig5d'].set_index(['group', 'series'])
    g = wells.groupby(['group', 'series']).value
    chk = pd.DataFrame({'n': g.size(), 'mean': g.mean()})
    assert len(chk) == len(s)
    for k, r in chk.iterrows():
        assert r.n == s.loc[k, 'n'] and abs(r['mean'] - s.loc[k, 'mean']) < 1e-4, k

    # every well is used (column kept for the schema: a well excluded from the fold change would be False)
    wells['used_for_fold'] = True
    w = wells.assign(effector=wells.group.map(name), experiment=wells.series)
    w_out = w.rename(columns={'group': 'group_workbook'})
    w_out = w_out.assign(group=w_out.effector, series=w_out.experiment, sample_id=w_out.experiment,
                         notes=w_out.notes.fillna(''))
    w_out[['group', 'series', 'replicate', 'value', 'sample_id', 'experiment', 'unit', 'notes', 'used_for_fold',
           'group_workbook', 'cell']].to_csv(OUT / 'fig5d_luc_wells.csv', index=False)

    u = w[w.used_for_fold]
    mean = u.groupby(['effector', 'experiment']).value.mean()
    norm = summ.assign(effector=summ.group.map(name),
                       v=summ.notes.str.replace('normalized=', '').astype(float)).set_index(['effector', 'series']).v
    rows, rep_map, sc = [], [], []
    for e in effs:
        if e == 'Control':
            for rp in ('Rep1', 'Rep2'):
                rows.append(dict(effector='Control', replicate=rp, fold_change=1.0,
                                 notes='GUS-GFP effector; every fold change is relative to the GUS-GFP control of the '
                                       'same experiment'))
            continue
        exps = sorted(u[u.effector == e].experiment.unique())
        assert len(exps) == 2 and set(exps) == set(REP[e]), (e, exps)
        for rp, x in zip(('Rep1', 'Rep2'), REP[e]):
            fc = mean[(e, x)] / mean[('Control', x)]
            wb = norm.get((e, x), np.nan)
            if not np.isnan(wb):
                assert f'{fc:.1f}' == f'{wb:.1f}', (e, x, fc, wb)
            nd = int(((u.effector == e) & (u.experiment == x)).sum())
            nc = int(((u.effector == 'Control') & (u.experiment == x)).sum())
            rows.append(dict(effector=e, replicate=rp, fold_change=round(fc, 6),
                             notes=f'experiment {x}; mean F/R {nd} DOF wells / mean F/R {nc} GUS-GFP wells'))
            rep_map.append(dict(effector=e, replicate=rp, experiment=x, n_dof_wells=nd, n_control_wells=nc,
                                fold_computed=round(fc, 3), fold_workbook_normalized=wb,
                                fold_kit_displayed=float(f'{fc:.1f}')))
            a = u[(u.effector == e) & (u.experiment == x)].value
            b = u[(u.effector == 'Control') & (u.experiment == x)].value
            sc.append(dict(effector=e, replicate=rp, experiment=x, p_student=stats.ttest_ind(a, b).pvalue,
                           p_welch=stats.ttest_ind(a, b, equal_var=False).pvalue))
    pd.DataFrame(rows).to_csv(OUT / 'fig5d_luc.csv', index=False)
    pd.DataFrame(rep_map).to_csv(OUT / 'fig5d_luc_rep_map.csv', index=False)

    plab = workbook_p_labels(xlsx) if xlsx else None
    out = []
    for r in sc:
        wb = plab.get((r['experiment'], r['effector']), '') if plab is not None else ''
        wbcat = wb.split(',')[-1].strip() if wb else ''
        out.append(dict(panel='Fig5d', comparison=f"{r['effector']} vs GUS-GFP control, {r['replicate']} "
                                                   f"(experiment {r['experiment']})",
                        label_drawn='none (heatmap shows fold change only; legend names no test)',
                        label_workbook=wb,
                        label_computed=f"Student t P = {r['p_student']:.2g} ({cat(r['p_student'])}); "
                                       f"Welch P = {r['p_welch']:.2g} ({cat(r['p_welch'])})",
                        agree=NO_LABELS if plab is None else
                        ('yes' if wbcat in (cat(r['p_student']), cat(r['p_welch'])) else 'no'),
                        note='two-sided t test on per-well Firefly/Renilla, DOF vs GUS-GFP of the same experiment; '
                             'not drawn'))
    pd.DataFrame(out).to_csv(CHK / 'STATS_CHECK_fig5.csv', index=False)
    print(pd.DataFrame(rep_map).to_string())
    print(pd.DataFrame(out)[['comparison', 'label_workbook', 'label_computed', 'agree']].to_string())
    if plab is None:
        print('note: no Source data workbook given; workbook P categories not compared (label_workbook empty)')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
