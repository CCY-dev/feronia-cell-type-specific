# Prototype; not yet carefully vetted.
"""Prepare tidy inputs for fig6a_candidates.R from existing q04_network result tables.

Claude Science frame cbb4e1ad-18c1-4945-a745-291b60aef9bb, 2026-09-20.
Reads only: revisions/2026-09-08_routeA_execution/q04_network/results/*.csv and the shared
36-gene DOF list (dof_family_list.csv, recovered from Science frame 53a87dd3 and copied to ../data/).
Writes: ../data/*.csv and yuc5_summary.json. No primary data, no network refit.

Usage: python3 prepare_inputs.py [output_root]
"""
import json, os, sys
import pandas as pd
from scipy.stats import fisher_exact

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Q = "/Users/tmn23/projects/fer-sc-manuscript-2026/revisions/2026-09-08_routeA_execution/q04_network/results/"
D = os.path.join(ROOT, "data"); os.makedirs(D, exist_ok=True)
YUC5 = "AT5G43890"

rt = pd.read_csv(Q + "network_rank_table_DE_repaired.csv", low_memory=False)
fam = pd.read_csv(Q + "TF_family_background_membership.csv")
sc = pd.read_csv(Q + "DOF_YUC5_all_network_scores.csv")
de = pd.read_csv(Q + "DOF_YUC5_current_DE_expression.csv")
en = pd.read_csv(Q + "TF_family_enrichment_universe_sensitivity.csv")
ye = pd.read_csv(Q + "all_retained_YUC5_edges.csv")
pr = pd.read_csv(Q + "supplied_prior_YUC5_regulators.csv")
u = pd.read_csv(Q + "h5ad_gene_universes.csv")
ev = pd.read_csv(Q + "tested_DOF_integrated_evidence_matrix.csv")
dl = pd.read_csv(os.path.join(D, "dof_family_list.csv"))

dofs_tfdb = set(fam[(fam.Family == "Dof") & (fam.definition == "PlantTFDB")].gene)
# AT1G65935 is the 37th coauthor-listed DOF (revisions/2026-09-20_fig6_brainstorm/fig6_dof_membership.csv);
# it is absent from the historical 36-gene list and from every network universe, DE table and score table.
all_dof = sorted(dofs_tfdb | set(dl.gene) | {"AT1G65935"})
assert len(all_dof) == 37, len(all_dof)
sym = dict(zip(dl.gene, dl.label)); sym["AT1G65935"] = "AT1G65935"
for g, s in zip(ev.gene, ev.symbol):
    sym.setdefault(g, s)
tested = set(ev.gene)

zone_order = ["Proliferation Domain", "Transition Domain", "Elongation", "Maturation", "Distal", "Proximal"]
def split_ctx(c):
    if c.startswith("Distal") or c.startswith("Proximal"):
        return c.split(" ", 1)
    return c.split("_", 1)

# nominations use the out-degree metric (Methods); Label is the published WT_specific / fer4_specific call
out = rt[rt.metric == "degree_centrality_out"]
nomd = out[out.gene.isin(all_dof) & out.Label.notna()][["gene", "cell_type", "Label"]].rename(columns={"cell_type": "context"})
assert nomd.groupby("Label").gene.nunique().to_dict() == {"WT_specific": 19, "fer4_specific": 23}

sc = sc[sc.gene.isin(all_dof)].copy()
sc["context"] = sc.cluster.str.replace("_WT_untreated", "").str.replace("_fer4_untreated", "")
scored = sc.groupby(["gene", "context"]).genotype.nunique().rename("n_geno").reset_index()
ctxs = sorted(set(scored.context) | set(out.cell_type.unique()),
              key=lambda c: (zone_order.index(split_ctx(c)[0]), split_ctx(c)[1]))
assert len(ctxs) == 44

rows = []
for g in all_dof:
    for c in ctxs:
        ng = scored[(scored.gene == g) & (scored.context == c)].n_geno
        ng = int(ng.iloc[0]) if len(ng) else 0
        lab = nomd[(nomd.gene == g) & (nomd.context == c)].Label
        if len(lab):
            status = "WT-associated" if lab.iloc[0] == "WT_specific" else "fer-4-associated"
        elif ng == 2:
            status = "scored, not nominated"
        elif ng == 1:
            status = "scored in one genotype only"
        else:
            status = "not in network"
        z, t = split_ctx(c)
        rows.append(dict(gene=g, symbol=sym.get(g, g), context=c, zone=z, tissue=t, status=status, tested=g in tested))
M = pd.DataFrame(rows)

de = de[de.gene.isin(all_dof)].copy()
norm = lambda c: c.replace(" ", "").replace("_", "").lower()
ctx_norm = {norm(c): c for c in ctxs}
de["context"] = [ctx_norm.get(norm(c)) for c in de.cell_type]
assert de.context.notna().all()
M = M.merge(de[["gene", "context", "current_DE_status", "logFC", "WT_mean_CPM"]], on=["gene", "context"], how="left")
M["de_called"] = M.current_DE_status.isin(["called_up", "called_down"])
M.to_csv(os.path.join(D, "dof_context_status.csv"), index=False)

per = M.groupby(["symbol", "gene", "tested"]).agg(
    n_WT=("status", lambda s: (s == "WT-associated").sum()),
    n_fer4=("status", lambda s: (s == "fer-4-associated").sum()),
    n_scored=("status", lambda s: s.str.startswith("scored").sum() + (s == "WT-associated").sum() + (s == "fer-4-associated").sum()),
    n_DE=("de_called", "sum")).reset_index()
per.to_csv(os.path.join(D, "dof_summary.csv"), index=False)
assert (per.n_WT + per.n_fer4 > 0).sum() == 26

en[(en.universe == "paired_score_KayTF") & (en.definition == "PlantTFDB")].to_csv(
    os.path.join(D, "family_enrichment_paired_universe.csv"), index=False)

# direction of nominations by tissue, DOF vs other TFs (all nominated TFs on the out-degree metric)
outn = out[out.Label.notna()].copy()
outn["tissue"] = [split_ctx(c)[1] for c in outn.cell_type]
outn["is_dof"] = outn.gene.isin(all_dof)
tf = outn[outn.Is_TF.astype(str) == "Yes"]
t2 = tf.groupby(["Label", "tissue", "is_dof"]).size().unstack(fill_value=0)
t2.columns = ["otherTF", "DOF"]; t2["DOF_frac"] = t2.DOF / (t2.DOF + t2.otherTF)
t2.reset_index().to_csv(os.path.join(D, "nomination_by_tissue_DOF_vs_otherTF.csv"), index=False)
tf.groupby(["tissue", "is_dof", "Label"]).size().unstack(fill_value=0).reset_index().to_csv(
    os.path.join(D, "direction_by_tissue.csv"), index=False)
endo = tf[tf.tissue == "Endodermis"]
ct = pd.crosstab(endo.is_dof, endo.Label)
odds, p = fisher_exact(ct.values)

# YUC5 regulators
ye["is_dof"] = ye.source.isin(all_dof)
yt = ye.groupby("tissue").agg(edges=("source", "size"), uniq_src=("source", "nunique"),
                              uniq_dof_src=("source", lambda s: s[s.isin(all_dof)].nunique())).reset_index()
yt.to_csv(os.path.join(D, "yuc5_regulators_by_tissue.csv"), index=False)
ye.to_csv(os.path.join(D, "yuc5_edges_with_dof_flag.csv"), index=False)
tissues_all = sorted(u.tissue.unique()); yuc_in = set(u[u.gene == YUC5].tissue)
pd.DataFrame({"tissue": tissues_all, "YUC5_in_target_universe": [t in yuc_in for t in tissues_all]}).to_csv(
    os.path.join(D, "yuc5_universe_by_tissue.csv"), index=False)
pr["is_dof"] = pr.source.isin(all_dof)
a, n = int(ye[ye.is_dof].source.nunique()), int(ye.source.nunique())
b, N = int(pr.is_dof.sum()), len(pr)
# retained regulators are nested within the prior set, so no independence test is computed; percentages only
json.dump({"prior_YUC5_regulators": N, "prior_DOF": b, "retained_unique_sources": n, "retained_DOF_sources": a,
           "retained_DOF_pct": round(100 * a / n, 1), "prior_DOF_pct": round(100 * b / N, 1),
           "endodermis_direction_table": {"DOF": ct.loc[True].to_dict(), "otherTF": ct.loc[False].to_dict()},
           "endodermis_direction_fisher_OR_exploratory_nonindependent_rows": odds,
           "endodermis_direction_fisher_p_exploratory_nonindependent_rows": p},
          open(os.path.join(D, "yuc5_summary.json"), "w"), indent=1, default=float)
print("ok", ct.values.tolist(), round(odds, 3), round(p, 4), a, n, b, N)
