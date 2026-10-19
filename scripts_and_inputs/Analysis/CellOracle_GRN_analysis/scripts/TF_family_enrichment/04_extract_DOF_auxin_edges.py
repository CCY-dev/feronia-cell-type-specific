#!/usr/bin/env python
# 04 — extract every DOF -> auxin-gene link from the CellOracle links objects (input for Fig. 6b). Run in the CellOracle environment.
"""
20260924 - Rebuild the panel-d edge set from OUR CellOracle links (Trevor's
`all_contexts_as_exported_union_edges.csv` was not shipped in from_Trevor/).

For every cell type x genotype context, keep the filtered (p<0.01, top-20k) edges whose SOURCE is a
DOF (PlantTFDB Family=="Dof") and whose TARGET is an auxin-pathway gene (YUC/TAA, GH3, Aux/IAA, SAUR).
Writes a tidy edge table so the panel can be redesigned without reloading 3.9 GB of .links objects.

The key design number this answers: what fraction of DOF->auxin edges are SHARED between genotypes
vs WT-only vs fer-4-only. That decides whether colouring edges by genotype declutters or worsens.
"""
import os, re, sys, glob, gzip
import pandas as pd

# Paths (override with environment variables; run from this folder). Mirrors 00_paths.R.
OUT = os.environ.get("CELLORACLE_DIR", os.path.join("..", "..", "outputs")) + "/"
FEAT = os.environ.get("FEATURES_TSV_GZ", os.path.join("..", "..", "inputs", "features.tsv.gz"))
PLANTTFDB = os.environ.get("PLANTTFDB_LIST", os.path.join("..", "..", "inputs", "Ath_TF_list.txt"))
RES = os.environ.get("OUT_DIR", os.path.join("..", "..", "outputs", "TF_family_enrichment"))
os.makedirs(RES, exist_ok=True)

# ---- gene sets -------------------------------------------------------------
fam = pd.read_csv(PLANTTFDB, sep="\t")
DOF = set(fam.loc[fam.Family == "Dof", "Gene_ID"].str.upper())

sym = {}
with gzip.open(FEAT, "rt") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            sym[p[0].upper()] = p[1]

def auxin_class(name):
    n = (name or "").upper()
    if re.fullmatch(r"YUC\d+", n) or n == "TAA1":            return "YUC/TAA"
    if re.fullmatch(r"GH3\.?\d+", n):                        return "GH3"
    if re.fullmatch(r"IAA\d+", n):                           return "Aux/IAA"
    if re.fullmatch(r"SAUR\d+", n):                          return "SAUR"
    return None

AUX = {g: auxin_class(s) for g, s in sym.items() if auxin_class(s)}
print(f"DOF sources: {len(DOF)} | auxin targets: {len(AUX)}", flush=True)
for k in ["YUC/TAA", "GH3", "Aux/IAA", "SAUR"]:
    print(f"   {k:8s} {sum(1 for v in AUX.values() if v==k)}", flush=True)

# ---- walk the links objects one at a time ----------------------------------
import celloracle as co   # noqa: E402

rows = []
files = sorted(glob.glob(OUT + "*_p01_20k.links.celloracle.links"))
print(f"\n{len(files)} links files", flush=True)
for f in files:
    base = os.path.basename(f).split("_rawRNA")[0]
    try:
        li = co.load_hdf5(f)
    except Exception as e:                                   # noqa: BLE001
        print(f"  !! {base}: {e}", flush=True); continue
    fl = getattr(li, "filtered_links", None) or {}
    for ctx, df in fl.items():
        if df is None or not len(df):
            continue
        d = df[df["source"].str.upper().isin(DOF) & df["target"].str.upper().isin(AUX)]
        if not len(d):
            continue
        d = d.copy()
        d["context_raw"] = ctx
        d["file"] = base
        rows.append(d)
    print(f"  {base}: {sum(len(r) for r in rows)} cumulative DOF->auxin edges", flush=True)
    del li

if not rows:
    sys.exit("no DOF->auxin edges found - check the source/target orientation")

E = pd.concat(rows, ignore_index=True)
E["source"] = E["source"].str.upper(); E["target"] = E["target"].str.upper()
E["genotype"] = E["context_raw"].str.extract(r"_(WT|fer4)_untreated$", expand=False)
E["cell_type"] = E["context_raw"].str.replace(r"_(WT|fer4)_untreated$", "", regex=True)
E["target_class"] = E["target"].map(AUX)
E["source_sym"] = E["source"].map(lambda g: sym.get(g, g))
E["target_sym"] = E["target"].map(lambda g: sym.get(g, g))
E = E.dropna(subset=["genotype"])
E.to_csv(os.path.join(RES, "DOF_auxin_edges_all_contexts.csv"), index=False)
print(f"\nwrote {len(E)} edge rows", flush=True)

# ---- the design number -----------------------------------------------------
pair = (E.groupby(["source", "target"])["genotype"].agg(lambda s: frozenset(s)).reset_index())
pair["state"] = pair["genotype"].map(lambda s: "shared" if len(s) == 2 else ("WT_only" if "WT" in s else "fer4_only"))
print("\n=== unique DOF->auxin PAIRS, pooled over cell types ===")
print(pair["state"].value_counts().to_string())
print(f"total unique pairs: {len(pair)}")

ctx = E.groupby(["source", "target", "cell_type"])["genotype"].agg(lambda s: frozenset(s)).reset_index()
ctx["state"] = ctx["genotype"].map(lambda s: "shared" if len(s) == 2 else ("WT_only" if "WT" in s else "fer4_only"))
print("\n=== per (pair x cell type) CONTEXT edges ===")
print(ctx["state"].value_counts().to_string())
print(f"total context-edges: {len(ctx)}")

print("\n=== out-degree per DOF (distinct auxin targets) ===")
od = E.groupby(["source_sym", "genotype"])["target"].nunique().unstack(fill_value=0)
od["total"] = od.sum(axis=1)
print(od.sort_values("total", ascending=False).head(15).to_string())
# write the genotype set SORTED: a raw frozenset prints in hash-randomised order, so the file would differ
# byte-for-byte between runs although its content is identical
pair_out = pair.assign(genotype=pair["genotype"].map(lambda s: ";".join(sorted(s))))
pair_out.to_csv(os.path.join(RES, "DOF_auxin_pair_states.csv"), index=False)
print("\nwrote pair-state table")
