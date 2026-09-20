#!/usr/bin/env python3
"""Prototype; not yet carefully vetted.

READ-ONLY remote check, NOT executed by claude-science (no Resnick credentials in that session).
Purpose: decide which saved component of endodermis_rawRNA_p01_20k.links.celloracle.links is stale in the three
inconsistent fer-4 clusters, by re-filtering the UNFILTERED raw links object with the run script's own parameters
(p=0.01, weight=coef_abs, threshold_number=20000) and comparing per-gene out-degree to (a) saved filtered_links and
(b) saved merged_score.degree_out. Loads ~0.5-1 GB of pickles; run on a compute node or with enough login memory.
Writes nothing remotely; prints JSON to stdout. Requires h5py, pandas. Does not import celloracle.
Expected outcome classes: raw-refilter == filtered_links (scores stale) | raw-refilter == merged_score (edges stale) |
matches neither (both from different fits / raw object itself re-saved).
"""
import json, pickle, zlib, sys
import h5py, pandas as pd
D = '/resnick/groups/NolanLab/cchen/Feronia/analysis/CellOracle/CellOracle_output/20260317_feronia_salinity_timecourse_WT_untreated_fer4_untreated_newIntegratedSeurat3_NoSuspiciousCluster/'
RAW = D + 'endodermis_rawRNA.links.celloracle.links'; FIL = D + 'endodermis_rawRNA_p01_20k.links.celloracle.links'
def read(f, key): return pickle.loads(zlib.decompress(f['&'+key][()].tobytes()))
def refilter(df, p=0.01, n=20000):
    d = df[df.p < p].copy(); d['coef_abs'] = d.coef_mean.abs()
    return d.sort_values('coef_abs', ascending=False).head(n)
out = {'raw': RAW, 'filtered': FIL, 'clusters': []}
with h5py.File(RAW, 'r') as fr, h5py.File(FIL, 'r') as ff:
    raw_links = read(fr, 'links_dict'); fil_links = read(ff, 'filtered_links'); score = read(ff, 'merged_score')
    out['raw_keys'] = sorted(k.lstrip('&') for k in fr.keys()); out['filtered_keys'] = sorted(k.lstrip('&') for k in ff.keys())
    for cl in fil_links:
        rf = refilter(raw_links[cl]); a = rf.groupby('source').size(); b = fil_links[cl].groupby('source').size()
        s = score[score.cluster == cl].degree_out
        idx = a.index.union(b.index).union(s.index)
        A, B, S = a.reindex(idx).fillna(0), b.reindex(idx).fillna(0), s.reindex(idx).fillna(0)
        out['clusters'].append({'cluster': cl, 'raw_rows': len(raw_links[cl]), 'raw_p_lt_0.01': int((raw_links[cl].p < 0.01).sum()), 'refiltered_edges': len(rf),
            'refilter_vs_filtered_links_genes_differ': int((A != B).sum()), 'refilter_vs_merged_score_genes_differ': int((A != S).sum()),
            'filtered_links_vs_merged_score_genes_differ': int((B != S).sum()),
            'refiltered_edge_set_equals_filtered_links': set(zip(rf.source, rf.target)) == set(zip(fil_links[cl].source, fil_links[cl].target))})
print(json.dumps(out, indent=1, default=str))
