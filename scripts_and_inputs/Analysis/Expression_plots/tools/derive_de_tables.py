"""Derive the two small tables behind Supplementary Fig. 2 from files in this repository.

  data/de_dof_yuc5.csv     muscat rows (every tested gene x cluster) for the ten assayed DOFs and YUC5,
                           with DE = the notebook's sig_DE_fil rule
  data/dof_nominations.csv CellOracle out-degree nominations of the ten DOFs (Methods: delta-rank < -123 -> more
                           central in fer-4, > +132 -> more central in WT, among genes ranked top-100 in either network)

Run from this folder:  python tools/derive_de_tables.py  [--repo ../../..]
"""
import argparse, pathlib
import numpy as np, pandas as pd

GENES = {"DOF4.6": "AT4G24060", "CDF4": "AT2G34140", "COG1": "AT1G29160", "OBP3": "AT3G55370",
         "DOF6": "AT3G45610", "DOF4.1": "AT4G00940", "DOF2.2": "AT2G28810", "DAG1": "AT3G61850",
         "OBP1": "AT3G50410", "OBP2": "AT1G07640", "YUC5": "AT5G43890"}
DE_P, DE_FC, DE_FRQ = 0.05, 1.5, 0.1          # muscat notebook: p_val_thres, FC_thres, filter_freq
NOM_LO, NOM_HI = -123, 132                     # Methods: 10th / 90th percentile delta-rank thresholds
MIN_CELLS = 5                                  # muscat pbDS(min_cells = 5)

ap = argparse.ArgumentParser(); ap.add_argument("--repo", default=str(pathlib.Path(__file__).resolve().parents[4]))
a = ap.parse_args(); repo = pathlib.Path(a.repo); here = pathlib.Path(__file__).resolve().parents[1]
inv = {v: k for k, v in GENES.items()}

sd = pd.read_csv(repo / "scripts_and_inputs/Analysis/Differential_expression_analysis/outputs/sig_de.csv.gz")
t = sd[sd.gene.isin(GENES.values())][["gene", "cluster_id", "logFC", "p_adj.loc", "WT.frq", "fer_4.frq"]].copy()
t["symbol"] = t.gene.map(inv)
t["DE"] = (t["p_adj.loc"] <= DE_P) & (t.logFC.abs() > np.log2(DE_FC)) & ((t["fer_4.frq"] >= DE_FRQ) | (t["WT.frq"] >= DE_FRQ))
t = t.rename(columns={"gene": "agi", "cluster_id": "cluster"})[["symbol", "agi", "cluster", "logFC", "p_adj.loc", "WT.frq", "fer_4.frq", "DE"]]
t.to_csv(here / "data/de_dof_yuc5.csv", index=False, float_format="%.6g")

rc = pd.read_csv(repo / "scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs/"
                 "20260317_feronia_salinity_WT_untreated_fer4_untreated_CellOracle_all_metrics_rank_change_wide_TF_GE.txt", sep="\t")
od = rc[(rc.metric == "degree_centrality_out") & rc.gene.isin(set(GENES.values()) - {GENES["YUC5"]})].rename(columns={
    "fer4_untreated-WT_untreated: rank_diff": "drank", "fer4_untreated-WT_untreated: gene_rank_1": "r_fer4",
    "fer4_untreated-WT_untreated: gene_rank_2": "r_WT"})
od = od[(od.drank < NOM_LO) | (od.drank > NOM_HI)].copy()
od["more_central_in"] = np.where(od.drank < NOM_LO, "fer-4", "WT")
cells = pd.read_csv(here / "data/si_expression_cells.csv.gz", usecols=["orig.ident", "time.celltype.anno.Li"])
ct = pd.crosstab(cells["time.celltype.anno.Li"], cells["orig.ident"])
elig = (ct >= MIN_CELLS).all(axis=1)
key = t.set_index(["agi", "cluster"])
def status(r):
    if not elig.get(r.cell_type, False): return "cluster not DE-eligible"
    if (r.gene, r.cell_type) not in key.index: return "gene not tested here"
    return "DE" if bool(key.loc[(r.gene, r.cell_type), "DE"]) else "tested, not DE"
od["status"] = od.apply(status, axis=1)
od["cells_WT"] = od.cell_type.map(ct.get("sc_291", 0) + ct.get("sc_293", 0))
od["cells_fer4"] = od.cell_type.map(ct.get("sc_292", 0) + ct.get("sc_294", 0))
od["symbol"] = od.gene.map(inv)
out = od.rename(columns={"gene": "agi", "cell_type": "cluster", "p_adj": "p_adj"})[
    ["symbol", "agi", "cluster", "more_central_in", "r_WT", "r_fer4", "drank", "log2FC", "p_adj", "status", "cells_WT", "cells_fer4"]]
out.sort_values(["symbol", "more_central_in", "cluster"]).to_csv(here / "data/dof_nominations.csv", index=False)
print(f"de_dof_yuc5.csv: {len(t)} rows, {int(t.DE.sum())} DE | dof_nominations.csv: {len(out)} nominations, "
      f"{int((out.status == 'DE').sum())} also DE, {int((out.status == 'cluster not DE-eligible').sum())} in non-eligible clusters")
