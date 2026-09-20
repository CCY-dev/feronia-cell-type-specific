#!/usr/bin/env python3
"""Prototype; not yet carefully vetted.

Independent numeric checks for the FER Figure 6 audit (task 4704c686-0f68-4d94-b928-d8e5540222eb,
claude-science, frame 5763b37b-3775-44c4-ba93-96f258753129). Read-only: reads frozen review_v2 tables,
q04 result tables, the ten-DOF inventory and the pinned public clone; writes only into ./results/.
Usage: python3 audit_checks.py [repo_root]   (pandas, numpy, scipy)
"""
import sys, os, re, json
import pandas as pd, numpy as np
from scipy.stats import fisher_exact, false_discovery_control

ROOT = sys.argv[1] if len(sys.argv) > 1 else '/Users/tmn23/projects/fer-sc-manuscript-2026'
R = f'{ROOT}/revisions/2026-09-20_compact_figures/review_v2'
Q = f'{ROOT}/revisions/2026-09-08_routeA_execution/q04_network/results/'
PUB = f'{ROOT}/revisions/2026-09-08_fresh_audit/work/feronia-cell-type-specific/scripts_and_inputs/Analysis/'
CO = PUB + 'CellOracle_GRN_analysis/outputs/'
INV = f'{ROOT}/revisions/2026-09-20_dof_subnetworks/inventory/'
ANN = f'{ROOT}/revisions/2026-09-20_dof_subnetworks/annotations/'
BS = f'{ROOT}/revisions/2026-09-20_fig6_brainstorm/'
CS = f'{ROOT}/.living/tasks/fer-fig6-brainstorm-2026-09-20/outputs/claude-science/data/'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results'); os.makedirs(OUT, exist_ok=True)
TISSUES = ['atrichoblast','columella','cortex','endodermis','lrc','pericycle','phloem','procambium','trichoblast','xylem']
BAD_CTX = ['Maturation_Endodermis','Proliferation Domain_Endodermis','Transition Domain_Endodermis']
res = {}
J = lambda o: o.item() if hasattr(o, 'item') else (bool(o) if isinstance(o, np.bool_) else str(o))

# 1. Sets and exclusive intersections
rt = pd.read_csv(Q+'network_rank_table_DE_repaired.csv', low_memory=False)
out = rt[rt.metric=='degree_centrality_out']; on = out[out.Label.notna()]; onTF = on[on.Is_TF=='Yes']
wt = set(onTF[onTF.Label=='WT_specific'].gene); f4 = set(onTF[onTF.Label=='fer4_specific'].gene); GRN = wt|f4
pub_wt = set(pd.read_csv(CO+'20260317_feronia_salinity_WT_untreated_fer4_untreated_CellOracle_WT_specific_TF.txt', sep='\t').gene)
pub_f4 = set(pd.read_csv(CO+'20260317_feronia_salinity_WT_untreated_fer4_untreated_CellOracle_fer4_specific_TF.txt', sep='\t').gene)
sde = pd.read_csv(PUB+'Differential_expression_analysis/outputs/sig_de_fil.csv', low_memory=False); DEG = set(sde.gene)
m = pd.read_csv(BS+'fig6_dof_membership.csv'); DOF = set(m.gene)
inter = pd.Series([(g in DOF, g in GRN, g in DEG) for g in DOF|GRN|DEG]).value_counts().sort_index()
ref = pd.read_csv(BS+'fig6_full_intersections.csv')
ref_map = {(bool(r.DOF), bool(r.GRN_nominated), bool(r.called_DEG)): int(r.genes) for r in ref.itertuples()}
res['sets'] = dict(WT_TF=len(wt), fer4_TF=len(f4), GRN_union=len(GRN), nonTF_nominee_excluded=int(on[on.Is_TF=='No'].gene.nunique()),
    WT_matches_public_export=wt==pub_wt, fer4_matches_public_export=f4==pub_f4, DEG=len(DEG), DOF=len(DOF), union=len(DOF|GRN|DEG),
    intersections_recomputed={str(k): int(v) for k, v in inter.items()},
    intersections_match_frozen=all(ref_map.get(k, 0)==v for k, v in inter.items()) and all(inter.get(k, 0)==v for k, v in ref_map.items()),
    membership_flags_match=bool((m.gene.isin(wt)==(m.nominated_WT.astype(str)=='True')).all() and (m.gene.isin(f4)==(m.nominated_fer4.astype(str)=='True')).all() and (m.gene.isin(DEG)==(m.called_DE_any.astype(str)=='True')).all()))

# 2. Family enrichment
fam = pd.read_csv(Q+'TF_family_background_membership.csv')
fp = fam[(fam.universe=='paired_score_KayTF')&(fam.definition=='PlantTFDB')]; U = set(fp.gene)
en = pd.read_csv(Q+'TF_family_enrichment_universe_sensitivity.csv'); e = en[(en.universe=='paired_score_KayTF')&(en.definition=='PlantTFDB')]
rows = []
for geno, fg in [('WT_specific', wt), ('fer4_specific', f4)]:
    F = fg & U
    for famname, grp in fp.groupby('Family'):
        K = set(grp.gene); a = len(K&F); b = len(F-K); c = len(K-F); d = len(U-K-F)
        rows.append(dict(Family=famname, genotype=geno, universe_n=len(K), foreground_n=a, p_greater=fisher_exact([[a,b],[c,d]], alternative='greater')[1]))
r = pd.DataFrame(rows)
for geno in r.genotype.unique():
    idx = r.genotype==geno; r.loc[idx, 'BH'] = false_discovery_control(r.loc[idx, 'p_greater'], method='bh')
cmp = r.merge(e[['Family','genotype','universe_n','foreground_n','odds_ratio','p_greater','p_adj_BH']], on=['Family','genotype'], suffixes=('','_src'))
cmp.to_csv(f'{OUT}/family_enrichment_recomputed_vs_source.csv', index=False)
mx = e.groupby('Family').agg(universe_n=('universe_n','first'), maxOR=('odds_ratio','max'), minBH=('p_adj_BH','min')).reset_index()
el = mx[mx.universe_n>=10].sort_values(['maxOR','Family'], ascending=[False, True]).reset_index(drop=True); el['display_rank'] = range(1, len(el)+1)
z = pd.read_csv(f'{R}/quantitative/source_top10_family.csv')
res['family_enrichment'] = dict(universe_genes=len(U), families=int(fp.Family.nunique()), multi_family_genes=int((fp.groupby('gene').Family.nunique()>1).sum()),
    WT_classifiable=len(wt&U), WT_unclassifiable=len(wt-U), fer4_classifiable=len(f4&U), fer4_unclassifiable=len(f4-U),
    counts_match=bool((cmp.universe_n==cmp.universe_n_src).all() and (cmp.foreground_n==cmp.foreground_n_src).all()),
    max_abs_p_diff=float((cmp.p_greater-cmp.p_greater_src).abs().max()), max_abs_BH_diff=float((cmp.BH-cmp.p_adj_BH).abs().max()),
    eligible_families_ge10=len(el), top10_recomputed=el.Family.head(10).tolist(), top10_frozen=z.sort_values('display_rank').Family.drop_duplicates().tolist(),
    rank10_vs_rank11_OR=[float(el.maxOR.iloc[9]), float(el.maxOR.iloc[10])], significant_families=mx[mx.minBH<.05].Family.tolist(),
    CI_clipping=dict(min_nonzero_lower=float(z.odds_lower[z.odds_lower>0].min()), max_upper=float(z.odds_upper.max()), axis=[0.125, 24]))
el.to_csv(f'{OUT}/eligible_family_ranking_recomputed.csv', index=False)

# 3. Edge-vs-score consistency, all 87 genotype-contexts (public exports)
audit = []
for t in TISSUES:
    E = pd.read_csv(CO+f'{t}_rawRNA_celloracle_network_iGRN.csv', low_memory=False, usecols=['source','target','cluster','coef_mean','p'])
    S = pd.read_csv(CO+f'{t}_rawRNA_celloracle_gene_score_iGRN.csv', low_memory=False, usecols=['Unnamed: 0','degree_out','degree_in','degree_all','cluster']).rename(columns={'Unnamed: 0':'gene'})
    for cl, Ec in E.groupby('cluster'):
        Sc = S[S.cluster==cl].set_index('gene')
        oc = Ec.groupby('source').size(); ic = Ec.groupby('target').size()
        so = Sc.degree_out[Sc.degree_out>0]; si = Sc.degree_in[Sc.degree_in>0]
        common = list(set(oc.index)&set(so.index)); ci = list(set(ic.index)&set(si.index))
        nodes_e = set(Ec.source)|set(Ec.target); nodes_s = set(Sc.index)
        cmpd = pd.DataFrame({'edges': oc, 'score': Sc.degree_out}).fillna(0)
        audit.append(dict(tissue=t, cluster=cl, n_edges=len(Ec), dup_edges=int(Ec.duplicated(['source','target']).sum()), max_p=float(Ec.p.max()),
            edge_sources=len(set(oc.index)), scored_out_sources=len(set(so.index)), src_in_edges_not_scored_out=len(set(oc.index)-set(so.index)),
            src_scored_out_not_in_edges=len(set(so.index)-set(oc.index)), out_degree_mismatch_common=int((oc.reindex(common)!=so.reindex(common)).sum()),
            in_degree_mismatch_common=int((ic.reindex(ci)!=si.reindex(ci)).sum()), score_rows=len(Sc), sum_degree_out=int(Sc.degree_out.sum()), sum_degree_in=int(Sc.degree_in.sum()),
            edge_nodes=len(nodes_e), score_nodes=len(nodes_s), score_nodes_not_in_edges=len(nodes_s-nodes_e), edge_nodes_not_in_scores=len(nodes_e-nodes_s),
            outdegree_corr_edges_vs_scores=float(cmpd.corr().iloc[0,1]) if len(cmpd)>1 else None,
            json_def_edge_sources_absent_from_scores=len(set(oc.index)-nodes_s),
            json_def_scored_sources_count_mismatch=int(((oc.reindex(Sc.index).fillna(0))!=Sc.degree_out).sum())))
A = pd.DataFrame(audit); A.to_csv(f'{OUT}/public_export_degree_audit_recomputed.csv', index=False)
bad = A[(A.out_degree_mismatch_common>0)|(A.in_degree_mismatch_common>0)|(A.src_in_edges_not_scored_out>0)|(A.src_scored_out_not_in_edges>0)|(A.dup_edges>0)]
res['edge_score_consistency'] = dict(genotype_contexts=len(A), all_20000_edges=bool((A.n_edges==20000).all()), any_duplicate_edges=bool((A.dup_edges>0).any()), max_p=float(A.max_p.max()),
    inconsistent_clusters=bad.cluster.tolist(), inconsistent_degree_out_sums=bad.sum_degree_out.tolist(),
    inconsistent_outdegree_corr=bad.outdegree_corr_edges_vs_scores.round(3).tolist(),
    json_definition_counts=bad[['cluster','json_def_edge_sources_absent_from_scores','json_def_scored_sources_count_mismatch']].to_dict('records'))

# 4. Exposure of nominations to the inconsistent contexts (exposure, NOT a prediction of corrected results)
exp_rows = onTF[onTF.cell_type.isin(BAD_CTX)]
only_bad = set(exp_rows.gene) - set(onTF[~onTF.cell_type.isin(BAD_CTX)].gene)
pc = onTF.groupby('cell_type').size().sort_values(ascending=False)
cs = pd.read_csv(CS+'dof_context_status.csv'); nom = cs[cs.status.isin(['WT-associated','fer-4-associated'])]
endo = nom[(nom.tissue=='Endodermis')&nom.tested]
sub = lambda L: set(onTF[(onTF.Label==L)&~onTF.cell_type.isin(BAD_CTX)].gene)
res['exposure_to_inconsistent_contexts'] = dict(TF_nomination_rows_total=len(onTF), rows_in_bad_contexts=len(exp_rows), TFs_nominated_only_in_bad_contexts=len(only_bad),
    WT_TFs_all=len(wt), WT_TFs_excluding_bad=len(sub('WT_specific')), fer4_TFs_all=len(f4), fer4_TFs_excluding_bad=len(sub('fer4_specific')),
    WT_DOF_all=len(wt&DOF), WT_DOF_excluding_bad=len(sub('WT_specific')&DOF), fer4_DOF_all=len(f4&DOF), fer4_DOF_excluding_bad=len(sub('fer4_specific')&DOF),
    DOF_context_nominations_total=len(nom), DOF_context_nominations_in_bad=int(nom.context.isin(BAD_CTX).sum()),
    tested10_endodermal_nominations=len(endo), tested10_endodermal_in_bad=int(endo.context.isin(BAD_CTX).sum()),
    tested10_endodermis_tiles=int(endo.symbol.nunique()), tested10_endodermis_tiles_entirely_bad=int((endo.groupby('symbol').context.apply(lambda s: s.isin(BAD_CTX).all())).sum()),
    top5_contexts_by_nominations=pc.head(5).to_dict(), median_context_nominations=float(pc.median()))
pc.rename('TF_nominations').reset_index().to_csv(f'{OUT}/nominations_per_context.csv', index=False)
exp_rows[exp_rows.gene.isin(DOF)].merge(m[['gene','symbol']], on='gene')[['symbol','gene','cell_type','Label']].to_csv(f'{OUT}/DOF_nominations_in_inconsistent_contexts.csv', index=False)

# 5. Matrix source vs public score exports
allS = pd.concat([pd.read_csv(CO+f'{t}_rawRNA_celloracle_gene_score_iGRN.csv', low_memory=False, usecols=['Unnamed: 0','cluster']).rename(columns={'Unnamed: 0':'gene'}) for t in TISSUES])
allS['context'] = allS.cluster.str.replace('_WT_untreated','',regex=False).str.replace('_fer4_untreated','',regex=False)
allS['geno'] = np.where(allS.cluster.str.endswith('WT_untreated'), 'WT', 'fer4')
ng = allS[allS.gene.isin(DOF)].groupby(['gene','context']).geno.nunique().rename('n_geno').reset_index()
chk = cs.merge(ng, on=['gene','context'], how='left').fillna({'n_geno':0}); isnom = chk.status.isin(['WT-associated','fer-4-associated'])
pred = np.where(isnom, 'nom', np.where(chk.n_geno==2, 'scored, not nominated', np.where(chk.n_geno==1, 'scored in one genotype only', 'not in network')))
lab = out[out.Label.notna()&out.gene.isin(DOF)][['gene','cell_type','Label']].rename(columns={'cell_type':'context'})
chk2 = cs.merge(lab, on=['gene','context'], how='left')
exp = np.where(chk2.Label=='WT_specific','WT-associated', np.where(chk2.Label=='fer4_specific','fer-4-associated','other'))
q = pd.read_csv(f'{R}/quantitative/plot_source_celltype_tested10.csv')
res['matrix'] = dict(score_contexts=int(allS.context.nunique()), rank_table_contexts=int(rt.cell_type.nunique()), context_only_scored_one_genotype=sorted(set(allS.context)-set(rt.cell_type)),
    paired_status_matches_public_scores=bool((pred==np.where(isnom,'nom',chk.status)).all()),
    nominations_match_rank_table=bool(((exp=='other')==~chk2.status.isin(['WT-associated','fer-4-associated'])).all() and (chk2.status[exp!='other']==exp[exp!='other']).all()),
    nominated_always_paired=bool((chk[isnom].n_geno==2).all()), tiles=len(q), tiles_K0=int((q.K==0).sum()), both_direction_tiles=int(((q.n_WT>0)&(q.n_fer>0)).sum()),
    any_gene_context_metric_with_both_labels=int((rt[rt.Label.notna()].groupby(['gene','cell_type','metric']).Label.nunique()>1).sum()),
    identities=q[['gene','symbol']].drop_duplicates().set_index('gene').symbol.to_dict())

# 6. Phloem/xylem feature selection (reconstructs the public run script's rule)
u = pd.read_csv(Q+'h5ad_gene_universes.csv'); KAY = set(pd.read_csv(PUB+'CellOracle_GRN_analysis/inputs/Kay_TF_thalemine_annotations.csv').GeneID)
clusters = sorted(sde.cluster_id.unique())
def de_genes(pat, excl=None):
    cl = [c for c in clusters if pat in c.lower() and (excl is None or excl not in c.lower())]; return set(sde[sde.cluster_id.isin(cl)].gene), cl
fx = {}
for t in ['phloem','xylem']:
    Ut = set(u[u.tissue==t].gene); sub_all, cl_all = de_genes(t); sub_clean, _ = de_genes(t, 'pericycle')
    fx[t] = dict(universe=len(Ut), universe_subset_of_TF_union_substringDE=Ut<=(sub_all|KAY), substring_DE_clusters=cl_all,
        admitted_only_via_pole_pericycle=len(((sub_all-sub_clean)-KAY)&Ut), clean_universe=len((sub_clean|KAY)&Ut), YUC5_in_universe='AT5G43890' in Ut, YUC5_in_clean='AT5G43890' in ((sub_clean|KAY)&Ut))
res['phloem_xylem_feature_selection'] = fx

# 7. Ten-DOF inventory, all-context union, regex coverage
te = pd.read_csv(INV+'tested10_retained_edges.csv', low_memory=False); tested10 = set(q.gene)
allE = pd.concat([pd.read_csv(CO+f'{t}_rawRNA_celloracle_network_iGRN.csv', low_memory=False).query('source in @tested10') for t in TISSUES])
ue = pd.read_csv(f'{R}/all_contexts/all_contexts_as_exported_union_edges.csv'); se = pd.read_csv(f'{R}/all_contexts/all_contexts_as_exported_source_context_edges.csv')
add = ue[ue.support_class=='additional_excluded_contexts_only']
res['inventory_and_union'] = dict(inventory_rows=len(te), rebuilt_rows=len(allE), edge_keys_identical=set(zip(allE.source,allE.target,allE.cluster))==set(zip(te.source,te.target,te.cluster)),
    mismatch_flag_clusters=sorted(te[te.edge_score_integrity=='edge_export_score_mismatch'].cluster.unique()),
    union_pairs=len(ue), union_targets=int(ue.target.nunique()), source_rows=len(se), contexts=int(se.cell_type.nunique()), genotype_contexts=int(se.cluster.nunique()),
    supported_seven=int((ue.support_class=='supported_in_current_seven').sum()), additions=len(add), additions_with_phloem_xylem=int((add.n_phloem_xylem_map_issue>0).sum()),
    additions_other_endodermis_only=int(((add.n_phloem_xylem_map_issue==0)&(add.n_other_endodermis>0)).sum()), additions_solely_mismatch=int(((add.n_phloem_xylem_map_issue==0)&(add.n_other_endodermis==0)&(add.n_endodermis_score_mismatch>0)).sum()),
    pairs_with_any_mismatch_support=int((ue.n_endodermis_score_mismatch>0).sum()), pairs_with_any_phloem_xylem_support=int((ue.n_phloem_xylem_map_issue>0).sum()))
rx = re.compile(r'^(IAA[0-9]|SAUR[0-9]|GH3[.]|TAA[0-9]|TAR[0-9]|YUC[0-9])')
gaf = pd.read_csv(ANN+'gaf_gene_symbols.tsv', sep='\t').set_index('gene'); names = te.drop_duplicates('target').set_index('target').target_name; sel = set(se.target)
miss = []
for agi, nm in names.items():
    if agi in sel or agi not in gaf.index: continue
    sym = str(gaf.loc[agi,'symbol']); syns = str(gaf.loc[agi,'synonym']).split('|'); full = str(gaf.loc[agi,'name'])
    if rx.match(sym) or any(rx.match(s) for s in syns) or re.search(r'auxin-responsive|Aux/IAA|SAUR|small auxin|YUCCA|indole-3-pyruvate|tryptophan aminotransferase|GH3', full, re.I):
        miss.append(dict(target=agi, export_name=nm, gaf_symbol=sym, gaf_name=full))
res['inventory_and_union']['targets_AGI_only_names'] = int((names.index==names.values).sum()); res['inventory_and_union']['regex_possible_misses'] = miss

json.dump(res, open(f'{OUT}/audit_checks_results.json','w'), indent=1, default=J)
print(json.dumps(res, indent=1, default=J))
