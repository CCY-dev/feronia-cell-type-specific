#!/usr/bin/env python
"""Derive the Fig. 5a, b tables in data/fig5/ from this repository's CellOracle outputs.

Inputs (all in scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs/):
  TF_family_enrichment/TF_family_enrichment_all_tests.csv    Stage 5 script 01 (110 Fisher tests, BH)
  TF_family_enrichment/DOF_auxin_pair_states.csv            Stage 5 script 04 (one row per DOF-target pair)
  TF_family_enrichment/DOF_auxin_edges_all_contexts.csv     Stage 5 script 04 (target symbols and classes)
  20261021_..._all_metrics_rank_change_wide_TF_GE.txt       Stage 4 step 2 (WT- / fer-4-specific labels)

Fig. 5b rules (= Stage 5 script 05 and the Methods):
  DOFs    = the ten DOFs tested in the pYUC5:LUC assay that are genotype-specific in at least one genotype (9; CDF4 is
            nominated in neither genotype and is not shown); group = WT only / fer-4 only / Both
  targets = counted over those DOFs: every YUC/TAA1 and GH3 gene reached, plus the 5 Aux/IAA and 3 SAUR genes reached
            by the most DOFs, both INCLUDING ties at the cut-off
  outdeg  = displayed targets reached by the DOF; indeg = DOFs reaching the target; one edge per DOF-target pair
            (union over cell types and genotypes)
  layout  = node y positions from the kit's own automatic layout (ferfigs.panels.fig6_network.auto_layout)

Usage (from the kit folder):
  python tools/fig5ab_tables.py           # write data/fig5/fig5a_*, fig5b_*
  python tools/fig5ab_tables.py --check   # recompute and compare with data/fig5/ without writing (exit 1 if not equal)
"""
from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

import pandas as pd
import yaml

KIT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(KIT))
from ferfigs.panels.fig6_network import auto_layout  # noqa: E402

CO = KIT.parent / 'CellOracle_GRN_analysis' / 'outputs'
ENR = CO / 'TF_family_enrichment'
WIDE = CO / '20261021_feronia_salinity_WT_untreated_fer4_untreated_CellOracle_all_metrics_rank_change_wide_TF_GE.txt'
OUT = KIT / 'data' / 'fig5'

# the ten DOFs tested in the pYUC5:LUC assay (AGI -> name used in the paper)
DOF10 = {'AT2G34140': 'CDF4', 'AT3G45610': 'DOF6', 'AT3G50410': 'OBP1', 'AT1G29160': 'COG1', 'AT3G61850': 'DAG1',
         'AT4G00940': 'DOF4.1', 'AT4G24060': 'DOF4.6', 'AT1G07640': 'OBP2', 'AT3G55370': 'OBP3', 'AT2G28810': 'DOF2.2'}


def derive():
    texts = {}
    # ---- a: the enrichment table, unchanged
    enr = pd.read_csv(ENR / 'TF_family_enrichment_all_tests.csv')
    assert len(enr) == 110 and enr.groupby('genotype').size().to_dict() == {'WT': 55, 'fer-4': 55}
    texts['fig5a_TF_family_enrichment_all_tests.csv'] = enr.to_csv(index=False)

    # ---- b: nomination groups from the rank-change table (out-degree centrality, TFs)
    w = pd.read_csv(WIDE, sep='\t')
    o = w[(w.metric == 'degree_centrality_out') & (w.Is_TF == 'Yes')]
    g_wt = set(o.gene[o.Label == 'WT_specific'].str.upper())
    g_fe = set(o.gene[o.Label == 'fer4_specific'].str.upper())

    def group(g):
        return ('Both' if g in g_wt and g in g_fe else 'fer-4 only' if g in g_fe else 'WT only' if g in g_wt
                else 'Not nominated')

    pr = pd.read_csv(ENR / 'DOF_auxin_pair_states.csv')
    ed = pd.read_csv(ENR / 'DOF_auxin_edges_all_contexts.csv')
    sy = ed[['target', 'target_sym', 'target_class']].drop_duplicates()
    assert not sy.target.duplicated().any()
    show = {g: n for g, n in DOF10.items() if group(g) != 'Not nominated'}
    k = pr[pr.source.isin(show)].merge(sy, on='target')
    assert set(k.source) == set(show), 'a shown DOF has no auxin edge'
    ind = k.groupby(['target', 'target_sym', 'target_class']).size().rename('indeg').reset_index()

    def top(cls, n):                                   # the n most-connected, including ties at the cut-off
        x = ind[ind.target_class == cls].sort_values(['indeg', 'target_sym'], ascending=[False, True])
        return x[x.indeg >= x.indeg.iloc[n - 1]] if len(x) >= n else x

    tn = pd.concat([ind[ind.target_class.isin(['YUC/TAA', 'GH3'])], top('Aux/IAA', 5), top('SAUR', 3)])
    e = k[k.target.isin(tn.target)].copy()
    dn = e.groupby('source').size().rename('outdeg').reset_index()
    dn['label'] = dn.source.map(DOF10)
    dn['group'] = dn.source.map(group)
    assert set(dn.source) == set(show), 'a shown DOF reaches none of the displayed targets'
    tn = tn.reset_index(drop=True)
    # drawing order (overlapping anti-aliased lines, so it changes pixels): edges grouped by target, targets in order of
    # first appearance in the pair table, pair-table order within a target (= what older pandas inner merges returned,
    # which is how the published panel's table was written)
    first = {t: i for i, t in reversed(list(enumerate(e.target)))}
    e = e.assign(_o=e.target.map(first)).sort_values('_o', kind='stable').drop(columns='_o').reset_index(drop=True)
    texts['fig5b_network_dof_nodes.csv'] = dn[['source', 'outdeg', 'label', 'group']].to_csv(index=False)
    texts['fig5b_network_target_nodes.csv'] = tn[['target', 'target_sym', 'target_class', 'indeg']].to_csv(index=False)
    texts['fig5b_network_edges.csv'] = e[['source', 'target', 'target_sym', 'target_class']].to_csv(index=False)

    # ---- b: node layout (the kit's automatic layout, groups in this order)
    cfg_b = yaml.safe_load(open(KIT / 'config' / 'fig5ab.yaml'))['panels']['b']
    cfg_b['group_order'] = ['WT only', 'fer-4 only', 'Both']
    lay = auto_layout(dn, tn, cfg_b)
    lab = dict(zip(dn.source, dn.label)) | dict(zip(tn.target, tn.target_sym))
    lay.insert(2, 'label', lay.node.map(lab))
    lay['y_mm'] = lay.y_mm.round(3)
    texts['fig5b_network_layout.csv'] = lay.to_csv(index=False)
    print(f"Fig. 5b: {len(dn)} DOFs ({', '.join(sorted(dn.label))}), {len(tn)} targets, {len(e)} edges; "
          f"not shown (nominated in neither genotype): {', '.join(DOF10[g] for g in sorted(set(DOF10) - set(show))) or 'none'}")
    return texts


def main(argv=None):
    ap = argparse.ArgumentParser(description='Derive and verify the Fig. 5a, b tables')
    ap.add_argument('--check', action='store_true', help='compare with data/fig5/ without writing')
    a = ap.parse_args(argv)
    texts = derive()
    diff = []
    for name, txt in texts.items():
        p = OUT / name
        if a.check:
            if not p.exists() or not pd.read_csv(io.StringIO(txt)).equals(pd.read_csv(p)):
                diff.append(name)
        else:
            OUT.mkdir(parents=True, exist_ok=True)
            p.write_text(txt)
    if a.check:
        print('data/fig5 equals a fresh derivation' if not diff else 'data/fig5 differs: ' + ', '.join(diff))
        return 1 if diff else 0
    print(f'wrote {len(texts)} tables to {OUT}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
