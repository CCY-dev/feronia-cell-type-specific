# FERONIA-cell-type-specific

Code and processed data for the single-cell analyses in:

**Cell type-specific gene manipulation uncovers FERONIA receptor kinase signaling specificity in Arabidopsis roots**

The analyses compare 5-day-old wild-type (WT, Col-0) and *fer-4* mutant Arabidopsis roots by single-cell RNA-seq (scRNA-seq). They cover cell-type annotation and integration, sample quality checks, cell-type-resolved differential expression, and CellOracle gene regulatory networks (GRNs) leading to the DOF transcription factor (TF) analysis in Fig. 6.

- **Raw data:** NCBI GEO [GSE325470](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325470). This holds FASTQ files and the per-sample Seurat objects `sc_291`–`sc_294`, which are too large for this repository.
- **Samples:** `sc_291` WT rep 1, `sc_292` *fer-4* rep 1, `sc_293` WT rep 2, `sc_294` *fer-4* rep 2 ([sample sheet](scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/Nolan_single_cell_lab_samples_sc291-294.csv)).

---

## Find the code for a figure

| Figure panel | What it shows | Code | Key output / input |
|---|---|---|---|
| Fig. 1a | UMAP of 44,664 cells by cell type and developmental stage | [Stage 1 – integration notebook](scripts_and_inputs/Preprocessing/scripts/202501226_Ref_integration_feronia_salinity_time_course_seu3.ipynb) | integrated Seurat object (GEO) |
| Fig. 1b, c | DEGs per subpopulation and per cell type, up and down | [Stage 3 – muscat notebook](scripts_and_inputs/Analysis/Differential_expression_analysis/scripts/20260317_feronia_salinity_timecourse_muscat_newIntegratedSeurat3_NoSuspiciousCluster.ipynb) | [`sig_de_fil.csv`](scripts_and_inputs/Analysis/Differential_expression_analysis/outputs/sig_de_fil.csv) |
| Fig. 1d | Hormone-regulated genes among DEGs | [Stage 3 – hormone overlap notebook](scripts_and_inputs/Analysis/Differential_expression_analysis/scripts/20260317_feronia_salinity_DE_hormone_geneoverlap_newIntegratedSeurat3_NoSuspiciousCluster.ipynb) | [hormone gene lists](scripts_and_inputs/Analysis/Differential_expression_analysis/inputs) |
| Extended Data Fig. 7d, e; Supplementary Figs 1–2 | *YUC5* and the ten assayed DOFs: UMAPs and violin plots, WT vs *fer-4* (SCT-normalised) | code being added | integrated Seurat object |
| Extended Data Fig. 1a | Pseudobulk sample-to-sample correlation | [Stage 2 – `20260612_pseudobulk_sample_correlation.R`](scripts_and_inputs/SC_data_check/20260612_pseudobulk_sample_correlation.R) | integrated Seurat object (GEO) |
| Extended Data Fig. 1b | Curated cell-type markers by genotype | [Stage 2 – `20260615_Hsu_FG_markers_dotplot_byGeno.R`](scripts_and_inputs/SC_data_check/20260615_Hsu_FG_markers_dotplot_byGeno.R) | integrated Seurat object (GEO) |
| Extended Data Fig. 2a | log₂FC heatmap of the 13 validated genes | [Stage 3 – `20260615_feronia_13genes_muscat_logFC_heatmap.R`](scripts_and_inputs/Analysis/Differential_expression_analysis/scripts/20260615_feronia_13genes_muscat_logFC_heatmap.R) | full muscat result (see *Not in this repository*) |
| Extended Data Fig. 3 | GO enrichment of all DEGs per cell-type/stage cluster (top 15 GO:BP terms per cluster) | [Stage 3 – muscat notebook](scripts_and_inputs/Analysis/Differential_expression_analysis/scripts/20260317_feronia_salinity_timecourse_muscat_newIntegratedSeurat3_NoSuspiciousCluster.ipynb) (`compareDE()`, GO heatmap step) | [`sig_de_fil.csv`](scripts_and_inputs/Analysis/Differential_expression_analysis/outputs/sig_de_fil.csv); g:Profiler release `e113_eg59_p19` |
| (proposed ED Fig. 3 version) | Up- and down-regulated DEGs, 15 selected GO terms each | [`20260317_Feronia_GO_terms_plot_new_integrated_seu3.R`](scripts_and_inputs/Analysis/Differential_expression_analysis/scripts/20260317_Feronia_GO_terms_plot_new_integrated_seu3.R) | [GO tables](scripts_and_inputs/Analysis/Differential_expression_analysis/outputs/GO), [term selections](scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/GO_term_selection) |
| Fig. 6a | TF-family enrichment among WT- and *fer-4*-nominated TFs | [Stage 5 – `02_Fig6a_TF_family_forest.R`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/TF_family_enrichment/02_Fig6a_TF_family_forest.R) | [`TF_family_enrichment_all_tests.csv`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs/TF_family_enrichment/TF_family_enrichment_all_tests.csv) |
| Fig. 6b | DOF → auxin-gene network for the ten assayed DOFs | [Stage 5 – `05_Fig6b_DOF_auxin_network.R`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/TF_family_enrichment/05_Fig6b_DOF_auxin_network.R) | [`Fig6b_DOF_auxin_network_stats.csv`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs/TF_family_enrichment/Fig6b_DOF_auxin_network_stats.csv) |
| Supplementary Fig. SX | Scheme of the TF-family enrichment method | [Stage 5 – `03_SupplFigSX_method_scheme.R`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/TF_family_enrichment/03_SupplFigSX_method_scheme.R) | [`SupplFigSX_method_scheme_stats.csv`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs/TF_family_enrichment/SupplFigSX_method_scheme_stats.csv) |
| Supplementary tables | WT- and *fer-4*-specific TFs; all 110 family tests | [Stage 4 – TF selection notebook](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/20260317_feronia_salinity_time_course_WT_untreated_fer4_untreated_CellOracle_newIntegratedSeurat3_NoSuspiciousCluster_select_TF_genes.ipynb); [Stage 5 – `01_TF_family_enrichment.R`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/TF_family_enrichment/01_TF_family_enrichment.R) | [`…_WT_specific_TF.txt`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs/20260317_feronia_salinity_WT_untreated_fer4_untreated_CellOracle_WT_specific_TF.txt), [`…_fer4_specific_TF.txt`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs/20260317_feronia_salinity_WT_untreated_fer4_untreated_CellOracle_fer4_specific_TF.txt), [`TF_family_enrichment_all_tests.csv`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs/TF_family_enrichment/TF_family_enrichment_all_tests.csv) |

---

## Workflow at a glance

```
GEO GSE325470 (FASTQ → scKB-COPILOT → per-sample Seurat objects sc_291–sc_294)
        │
        ▼
Stage 1  Preprocessing ............ annotation transfer + integration with the root atlas  → Fig. 1a
        │
        ├──► Stage 2  SC_data_check ........ pseudobulk correlation, marker dot plot      → ED Fig. 1
        │
        ├──► Stage 3  Differential expression (muscat) → DEGs, GO, hormone overlap        → Fig. 1b–d, ED Figs 2a, 3
        │
        └──► Stage 4  CellOracle GRNs per cell type → network scores → genotype-specific TFs
                    │
                    ▼
             Stage 5  TF-family enrichment + DOF → auxin network                         → Fig. 6a,b, Suppl. Fig. SX
```

Each stage below lists its inputs, the scripts in run order, and its outputs.

---

## Stage 1 – Preprocessing: annotation and integration

Folder: [`scripts_and_inputs/Preprocessing/`](scripts_and_inputs/Preprocessing)

Reads were demultiplexed with Cell Ranger `mkfastq` (v3.1.0) and quality-filtered with the scKB-COPILOT pipeline. Both steps run upstream of this repository, and their outputs are on GEO. Cell type and developmental stage were assigned by label transfer from the published Arabidopsis root atlas (GEO [GSE152766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152766)).

1. [`202501226_feronia_salinity_time_course_update_metadata.ipynb`](scripts_and_inputs/Preprocessing/scripts/202501226_feronia_salinity_time_course_update_metadata.ipynb) (conda env `seu4_label_transfer`). Exports the label-transfer annotations from the Seurat 4 objects so they can be added to the Seurat 3 objects.
2. [`202501226_Ref_integration_feronia_salinity_time_course_seu3.ipynb`](scripts_and_inputs/Preprocessing/scripts/202501226_Ref_integration_feronia_salinity_time_course_seu3.ipynb) (conda env `scKB_copilot`, Seurat v3.1.5). Updates the metadata, integrates the four samples with the reference atlas, and runs PCA, UMAP and clustering. It then removes the reference cells and draws the UMAPs (**Fig. 1a**). One cluster found only in `sc_291` was flagged as suspicious and removed. All downstream analyses use this object (`…_seu3_sc291suspiciousClusterRemoved.rds`, "NoSuspiciousCluster" in file names).

Input: [`Protoplasting_DEgene_FC2_list.txt`](scripts_and_inputs/Preprocessing/inputs/Protoplasting_DEgene_FC2_list.txt), the protoplasting-responsive genes (fold change > 2) used by both notebooks.

## Stage 2 – Sample quality checks

Folder: [`scripts_and_inputs/SC_data_check/`](scripts_and_inputs/SC_data_check)

Both scripts read the integrated object from Stage 1 and run in conda env `muscat` (Seurat 5).

- [`20260612_pseudobulk_sample_correlation.R`](scripts_and_inputs/SC_data_check/20260612_pseudobulk_sample_correlation.R). Sums raw counts per sample, converts to log₂(CPM + 1) and keeps genes with ≥ 10 counts. It then computes Pearson and Spearman sample-by-sample correlations and draws clustered heatmaps (**Extended Data Fig. 1a**).
- [`20260615_Hsu_FG_markers_dotplot_byGeno.R`](scripts_and_inputs/SC_data_check/20260615_Hsu_FG_markers_dotplot_byGeno.R). Dot plot of literature-curated cell-type markers (SEMITONES ∪ novoSpaRc sets) across cell types, faceted by genotype (**Extended Data Fig. 1b**). Marker-to-cell-type mapping: [`inputs/20260615_Hsu_FG_union_celltype_mapping.csv`](scripts_and_inputs/SC_data_check/inputs/20260615_Hsu_FG_union_celltype_mapping.csv).

## Stage 3 – Differential expression, GO and hormone overlap

Folder: [`scripts_and_inputs/Analysis/Differential_expression_analysis/`](scripts_and_inputs/Analysis/Differential_expression_analysis)

1. [`20260317_feronia_salinity_timecourse_muscat_newIntegratedSeurat3_NoSuspiciousCluster.ipynb`](scripts_and_inputs/Analysis/Differential_expression_analysis/scripts/20260317_feronia_salinity_timecourse_muscat_newIntegratedSeurat3_NoSuspiciousCluster.ipynb) (conda env `muscat`). This step:
   - runs pseudobulk differential expression with muscat (v1.16.0, edgeR), comparing *fer-4* with WT within each cell-type × developmental-stage cluster;
   - draws the DEG summary panels (**Fig. 1b, c**);
   - runs GO enrichment with gprofiler2 (v0.2.3). g:Profiler is pinned to the archived release `e113_eg59_p19` (Ensembl 113 / Ensembl Genomes 59) with `set_base_url()`; this release reproduces the published GO tables exactly (all terms and clusters, |Δ(−log10 P)| < 1e-13) and the 68 terms of Extended Data Fig. 3. The live g:Profiler server gives slightly different results.
   - Outputs: [`outputs/sig_de_fil.csv`](scripts_and_inputs/Analysis/Differential_expression_analysis/outputs/sig_de_fil.csv) (significant DEGs), [`outputs/sig_de.csv.gz`](scripts_and_inputs/Analysis/Differential_expression_analysis/outputs/sig_de.csv.gz) (full muscat result, every tested gene × cluster), and the GO tables in [`outputs/GO/`](scripts_and_inputs/Analysis/Differential_expression_analysis/outputs/GO) (the published up/down tables plus the complete g:Profiler `e113_eg59_p19` results). The GO background is [`inputs/all_bg.csv`](scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/all_bg.csv), the expressed genes per cluster.
2. [`20260317_Feronia_GO_terms_plot_new_integrated_seu3.R`](scripts_and_inputs/Analysis/Differential_expression_analysis/scripts/20260317_Feronia_GO_terms_plot_new_integrated_seu3.R). Heatmaps of 15 selected GO terms for up- and down-regulated genes (a proposed version of Extended Data Fig. 3; selections in [`inputs/GO_term_selection/`](scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/GO_term_selection)). Note: the selections were made from an earlier GO table, and four selected terms are absent from the current tables.
3. [`20260317_feronia_salinity_DE_hormone_geneoverlap_newIntegratedSeurat3_NoSuspiciousCluster.ipynb`](scripts_and_inputs/Analysis/Differential_expression_analysis/scripts/20260317_feronia_salinity_DE_hormone_geneoverlap_newIntegratedSeurat3_NoSuspiciousCluster.ipynb) (conda env `muscat`). Tests the overlap of DEGs in each cluster with hormone-regulated gene sets using GeneOverlap (v1.38.0), with Fisher's exact test (**Fig. 1d**). Gene sets: [`IAA`](scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/IAA.csv), [`BR`](scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/BR.csv), [`JA`](scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/Hickman_2017_JA.csv), [`ABA`](scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/ABA.csv), [`CK`](scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/CK.csv), [`ACC`](scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/ACC.csv), [`SA`](scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/SA.csv), [`GA`](scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/GA.csv).
4. [`20260615_feronia_13genes_muscat_logFC_heatmap.R`](scripts_and_inputs/Analysis/Differential_expression_analysis/scripts/20260615_feronia_13genes_muscat_logFC_heatmap.R). Cluster-level log₂FC heatmap of the 13 genes validated by qRT-PCR, with asterisks from `p_adj.loc` (**Extended Data Fig. 2a**).
   - The companion [`20260615_feronia_13genes_umap_featureplot_q95.R`](scripts_and_inputs/Analysis/Differential_expression_analysis/scripts/20260615_feronia_13genes_umap_featureplot_q95.R) draws UMAP expression plots of the same genes.

Gene identifiers are mapped to symbols with [`features.tsv`](scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/features.tsv).

## Stage 4 – CellOracle gene regulatory networks

Folder: [`scripts_and_inputs/Analysis/CellOracle_GRN_analysis/`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis)

1. [`…_slurmArray_run_celloracle_grn.py`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/20260317_feronia_salinity_time_course_WT_untreated_fer4_untreated_CellOracle_newIntegratedSeurat3_NoSuspiciousCluster_slurmArray_run_celloracle_grn.py) (CellOracle v0.20.0). This step:
   - runs once per cell type as a SLURM array job, using the base GRN described in Methods;
   - fits networks for each genotype and developmental stage;
   - keeps links with p < 0.01, the top 20,000 by |coefficient|;
   - computes network scores, including out-degree centrality.
   - Outputs per cell type in [`outputs/`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs): `<cell type>_rawRNA_celloracle_network_iGRN.csv` (filtered links) and `<cell type>_rawRNA_celloracle_gene_score_iGRN.csv` (network scores). Example: [cortex links](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs/cortex_rawRNA_celloracle_network_iGRN.csv) and [scores](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs/cortex_rawRNA_celloracle_gene_score_iGRN.csv).
2. [`…_select_TF_genes.ipynb`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/20260317_feronia_salinity_time_course_WT_untreated_fer4_untreated_CellOracle_newIntegratedSeurat3_NoSuspiciousCluster_select_TF_genes.ipynb) (conda env `mast`). This step:
   - ranks genes by out-degree centrality within each network;
   - computes the rank change Δrank = rank(*fer-4*) − rank(WT) in each matched context;
   - joins TF annotations and expression log₂FC;
   - nominates genotype-specific TFs.
   - Outputs: [`…_all_metrics_rank_change_wide_TF_GE.txt`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs/20260317_feronia_salinity_WT_untreated_fer4_untreated_CellOracle_all_metrics_rank_change_wide_TF_GE.txt), [`…_WT_specific_TF.txt`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs/20260317_feronia_salinity_WT_untreated_fer4_untreated_CellOracle_WT_specific_TF.txt) and [`…_fer4_specific_TF.txt`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs/20260317_feronia_salinity_WT_untreated_fer4_untreated_CellOracle_fer4_specific_TF.txt) (supplementary tables).

**Reproducibility note.** Link coefficients from CellOracle depend on Python's string-hash seed (see [`envs/README.md`](envs/README.md#celloracle-set-pythonhashseed)); run the script with `export PYTHONHASHSEED=0`. The published networks were made without a fixed seed, so they cannot be regenerated exactly. In a rerun of columella from the same input, code and environment (20 Oct 2026), 87.9–89.1% of the top-20,000 links per network were the same, out-degree centrality agreed with Spearman ρ = 0.997 in all four networks, and 44–47 of the 50 genes with the highest out-degree were shared.

Inputs in [`inputs/`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/inputs): [`Kay_TF_thalemine_annotations.csv`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/inputs/Kay_TF_thalemine_annotations.csv) (the TF list used to build the networks) and [`features.tsv.gz`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/inputs/features.tsv.gz) (AGI → gene symbol).

## Stage 5 – TF-family enrichment and the DOF → auxin network (Fig. 6a, b)

Folder: [`scripts/TF_family_enrichment/`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/TF_family_enrichment). Its own **[README](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/TF_family_enrichment/README.md)** gives the run order, thresholds, background definition and package versions. In brief:

| Run | Script | Produces |
|---|---|---|
| 0 | [`00_paths.R`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/TF_family_enrichment/00_paths.R) | paths, which can be overridden by environment variables |
| 1 | [`01_TF_family_enrichment.R`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/TF_family_enrichment/01_TF_family_enrichment.R) | family × genotype Fisher tests over 1,292 background TFs, BH over 110 tests |
| 2 | [`02_Fig6a_TF_family_forest.R`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/TF_family_enrichment/02_Fig6a_TF_family_forest.R) | **Fig. 6a** |
| 3 | [`03_SupplFigSX_method_scheme.R`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/TF_family_enrichment/03_SupplFigSX_method_scheme.R) | **Supplementary Fig. SX** |
| 4 | [`04_extract_DOF_auxin_edges.py`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/TF_family_enrichment/04_extract_DOF_auxin_edges.py) | [`DOF_auxin_edges_all_contexts.csv`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs/TF_family_enrichment/DOF_auxin_edges_all_contexts.csv), [`DOF_auxin_pair_states.csv`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs/TF_family_enrichment/DOF_auxin_pair_states.csv) |
| 5 | [`05_Fig6b_DOF_auxin_network.R`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/TF_family_enrichment/05_Fig6b_DOF_auxin_network.R) | **Fig. 6b** |

PlantTFDB family assignments come from [`Ath_TF_list.txt`](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/inputs/Ath_TF_list.txt).

---

## Running the code

- **Environments.** The exact software environments are recorded in [`envs/`](envs) (see [`envs/README.md`](envs/README.md)). They were captured on 20 Oct 2026 from the conda environments that produced the results, and they match the package versions printed by `sessionInfo()` in the notebooks:

  | Environment | Used for | R / Python |
  |---|---|---|
  | `scKB_copilot` | Stage 1 integration (Seurat 3.1.5) | R 4.1.3 |
  | `seu4_label_transfer` | Stage 1 label-transfer metadata (Seurat 4.4.0) | R 4.4.3 |
  | `muscat` | Stages 2–3 (Seurat 5.2.1, muscat 1.16.0, edgeR 4.0.16, GeneOverlap 1.38.0, gprofiler2 0.2.3) | R 4.3.3 |
  | `mast` | Stage 4 TF selection notebook (Seurat 5.0.0) | R 4.4.3 |
  | `cell_oracle` | Stage 4 CellOracle run and Stage 5 script 04 (CellOracle 0.20.0) | Python 3.10.14 |

  For each environment, `envs/` holds the conda specification (`.yml`), the exact package list for linux-64 (`.linux-64.explicit.txt`), every installed R package with its version and source (`.R_packages.csv`), `pip freeze`, and a script that reinstalls the R packages that were installed by R rather than conda (`install_R_packages_<env>.R`). Exact archives of the environments (conda-pack) are deposited with the data.
- **Pinned external service.** GO enrichment uses g:Profiler release `e113_eg59_p19` (set in the muscat notebook).
- **Determinism.** Seurat's PCA, UMAP and clustering use fixed default seeds, and CellOracle 0.20.0 fixes `random_state = 123` in its bagging ridge models, so reruns in the same environments give the same results.
- **Inputs and outputs.** Scripts read their inputs from this repository: set the environment variable `FERONIA_REPO` to your clone (default: the current folder). Scripts that write files use `OUT_DIR` (default: the current folder). The sample sheet is [`Nolan_single_cell_lab_samples_sc291-294.csv`](scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/Nolan_single_cell_lab_samples_sc291-294.csv) (columns `sample`, `genotype`, `rep`).
- **Large objects.** Several notebooks and scripts read the Seurat objects and CellOracle inputs from the lab's cluster paths (`/resnick/groups/NolanLab/...`); point these at your downloaded copies (see below). The Stage 5 scripts use relative paths (see their [README](scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/TF_family_enrichment/README.md)).
- **Heavy steps.** Stage 1, the Seurat-loading scripts in Stages 2–3, and Stage 4 step 1 need the integrated object (about 20 GB in memory) or a cluster. Stage 3 plotting and all of Stage 5 run from the tables in this repository.

## Large files deposited outside this repository

| Item | Size | Needed by | Where |
|---|---|---|---|
| Per-sample Seurat objects (`sc_291`–`sc_294`) and FASTQ | | Stage 1 | GEO GSE325470 |
| Integrated, annotated Seurat object (`Feronia_salinity_timecourse_4S_20251226_seu3_sc291suspiciousClusterRemoved.rds`) | 19.5 GB | Stages 2–4 | data deposit (DOI to be added) |
| CellOracle input (`20250509_4S_feronia_salinity_time_course_untreated_celloracle_input.h5ad`) | 5.4 GB | Stage 4 step 1 | data deposit (DOI to be added) |
| Base GRN (`base_GRN_TG_to_TF_dictionary.pkl`; built as in Methods from scATAC-seq GSE155304 and DAP-seq) | 1.3 GB | Stage 4 step 1 | data deposit (DOI to be added) |
| CellOracle links objects (`*_rawRNA.links.celloracle.links`, `*_p01_20k.links.celloracle.links`) | 11 GB with the per-cell-type outputs | Stage 5 script 04 | data deposit (DOI to be added) |
| Exact environment archives (conda-pack, linux-64) | | all stages | data deposit (DOI to be added) |

Sizes and SHA-256 checksums of the 30 files to deposit (39.3 GB in total) are in [`DATA_DEPOSIT_MANIFEST.tsv`](DATA_DEPOSIT_MANIFEST.tsv).
