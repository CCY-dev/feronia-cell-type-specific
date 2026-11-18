# Single-cell expression figures: Supplementary Figs 1–2 (and a *YUC5* dot-plot alternative)

This folder draws the four single-cell expression figures for *YUC5* and the ten DOF transcription factors tested in the manuscript (DOF4.6, CDF4, COG1, OBP3, DOF6, DOF4.1, DOF2.2, DAG1, OBP1, OBP2). One command rebuilds all four from small tables kept here, so no Seurat object is needed to redraw them.

| Output (`output/figures/`) | Manuscript panel | Size (mm) |
|---|---|---|
| `ED7d_YUC5_UMAP.pdf` | *YUC5* expression on the UMAP, WT and *fer-4*. **Not used in the paper:** Extended Data Fig. 5d,e are drawn by [`../YUC5_expression_ED5/`](../YUC5_expression_ED5) | 92.5 × 47.9 |
| `ED7e_YUC5_dotplot.pdf` | *YUC5* dot plot by cell type and stage, WT and *fer-4* (not used, see above) | 180 × 34.5 |
| `SupplFig1_DOF_expression.pdf` | Supplementary Fig. 1: expression of the ten DOFs (a UMAPs, b dot plot) | 180 × 137.4 |
| `SupplFig2_DOF_DE_and_network.pdf` | Supplementary Fig. 2: where the DOFs and *YUC5* are differentially expressed (a) and where CellOracle nominated the DOFs (b) | 180 × 116.3 |

## Rebuild the figures

The figure files in `output/figures/` are not tracked in git; this command regenerates them (checked 18 Nov 2026: Supplementary Figs 1 and 2 are pixel-identical to the manuscript versions).

Needs Python ≥ 3.10 with the packages in `requirements.txt` (PyMuPDF only for the optional before/after pages).

```bash
cd scripts_and_inputs/Analysis/Expression_plots
python build.py                      # all four figures, source data and checks, written to output/
python build.py --only SupplFig2_DOF_DE_and_network --formats pdf png
python build.py --out /tmp/check     # write somewhere else
```

Each figure prints `[ok]` when its content fits its page. `output/build.json` records the SHA-256 of every input, the config and the code, the package versions and the per-figure checks. Two builds give byte-identical files (checked 20 Oct 2026; PDFs carry no creation date). Fonts: Arial, embedded as TrueType.

## Where the inputs come from

| File in `data/` | What it is | Made by |
|---|---|---|
| `si_expression_cells.csv.gz` | 44,664 cells (sc_291–sc_294): sample, genotype, cell type, stage, cluster (`time.celltype.anno.Li`), UMAP coordinates (the `umap` reduction of Fig. 1a) and SCT `data`-slot values for the 11 genes. SHA-256 in `.sha256`. | `tools/hpc/extract_si_expression.R` (below) |
| `si_expression_object_summary.txt`, `si_expression_sessionInfo.txt` | Object summary (assays, reductions, cells per cluster and genotype, gene ranges) and R `sessionInfo()` of the export job | same job |
| `de_dof_yuc5.csv` (+ README) | muscat result rows for the 11 genes in every tested cluster, with `DE` = the notebook's `sig_DE_fil` rule (`p_adj.loc` ≤ 0.05, \|fold change\| > 1.5, expressed in ≥ 10 % of cells in either genotype) | `tools/derive_de_tables.py` from `../Differential_expression_analysis/outputs/sig_de.csv.gz` |
| `dof_nominations.csv` (+ README) | CellOracle out-degree nominations of the ten DOFs from the seeded rerun (Δrank < −98 → more central in *fer-4*; > +105 → more central in WT; the 10th / 90th percentiles, computed from the table; genes in the top 100 of either network), with the DE status of the same cluster | `tools/derive_de_tables.py` from `../CellOracle_GRN_analysis/outputs/…_all_metrics_rank_change_wide_TF_GE.txt` |
| `de_tested_clusters.txt` | The 37 clusters muscat tested. The build checks that the cell-count rule below selects exactly these. | from `sig_de.csv.gz` |
| `fig1a_palette.csv` | Fig. 1a cell-type colours (not used by the four figures; kept for the optional reference map) | figure kit v3.2 |

**Re-export the per-cell table** (only if the Seurat object changes). On Resnick, from a copy of `tools/hpc/`:

```bash
sbatch run_extract.sh      # conda env scKB_copilot (R 4.1.3, Seurat 3.1.5); reads the integrated object, never writes it
```

Set the object path at the top of `extract_si_expression.R`. The run on 20 Oct 2026 (SLURM job 5098874) took 2.7 min with 128 GB; its log and stdout are in `tools/hpc/`.

**Re-derive the DE and nomination tables** (after any change to the muscat or CellOracle outputs):

```bash
python tools/derive_de_tables.py --repo ../../..      # rewrites data/de_dof_yuc5.csv and data/dof_nominations.csv
python build.py
```

The Δrank thresholds are computed from the rank-change table (10th / 90th percentile of out-degree Δrank) and checked against its `Label` column, so nothing needs editing after a CellOracle rerun.

## Rules the code follows

- **Expression** is the SCT assay `data` slot, as in the repository's other expression plots.
- **Clusters** are the `time.celltype.anno.Li` cell-type × stage clusters, ordered as in the Fig. 1a key and then by stage. Only clusters eligible for differential expression are drawn: every one of the four samples has at least 5 cells (muscat `pbDS(min_cells = 5)`), which gives 37 of 46 clusters. All 46 clusters remain in the source data with a flag.
- **UMAPs** use the coordinates of Fig. 1a. Every map has the same frame and equal x/y scaling. Colours run from light grey to blue (the Seurat FeaturePlot defaults), from the gene's minimum to its maximum over all cells. Expressing cells are drawn on top. Points are rasterised at 600 ppi inside the vector PDF.
- **Dot plots:** dot area ∝ % of cells with SCT > 0, and colour = mean SCT value of all cells in the cluster (arithmetic mean). One colour scale covers all genes and figures. WT and *fer-4* are grouped in separate rows. `–` marks a genotype with no cells in that cluster.
- **Supplementary Fig. 2a:** colour = log₂ fold change (*fer-4* vs WT), shown only where the gene is DE, with the fold change printed (↑ up, ↓ down in *fer-4*). Grey = tested, not DE. White = not tested by muscat (too few expressing cells).
- **Supplementary Fig. 2b:** nominations use the Fig. 5a colours (WT blue circles, *fer-4* orange diamonds). A box marks a nomination in a cluster where the same DOF is DE. Nine nominations fall in 4 of the 9 clusters that are not drawn; they are listed in `data/dof_nominations.csv`.
- No statistics are computed in this folder. All tests come from the muscat and CellOracle stages.

## Source data

`output/source_data/si_expression_cluster_summary.csv` holds, per gene × cluster × genotype: n cells, n and % expressing, arithmetic mean, Seurat DotPlot mean, median and maximum. `si_cluster_cells_per_sample.csv` holds the cells per cluster and sample behind the eligibility rule. For Supplementary Fig. 2, the source data are `data/de_dof_yuc5.csv` and `data/dof_nominations.csv`.

## Resolved (21 Oct 2026): fer-4 endodermis network scores

The March CellOracle run had three *fer-4* endodermis score tables overwritten by concurrent array tasks sharing CellOracle's temporary `./network_analysis/` folder. All networks were rerun with per-task working folders and a fixed hash seed (Stage 4, `20261020_celloracle_grn_seeded_isolated.py`), and `dof_nominations.csv` is derived from the rerun. Supplementary Fig. 2b shows the rerun nominations.

## Draft legends

**Supplementary Fig. 1.** Expression of the ten DOF transcription factors tested in this study. **a**, UMAPs coloured by expression (SCT-normalised) in WT and *fer-4*. **b**, Dot plot across cell types and developmental stages, grouped by genotype. Dot size, percentage of cells expressing; colour, mean expression. Only clusters eligible for differential expression (≥ 5 cells in each sample) are shown.

**Supplementary Fig. 2.** Differential expression and network nomination of the ten DOFs. **a**, Differential expression in *fer-4* vs WT for each cell type and developmental stage (muscat pseudobulk; adjusted P ≤ 0.05, fold change > 1.5, expressed in ≥ 10 % of cells in either genotype). Colour and numbers, fold change where differentially expressed. Grey, tested and not differentially expressed. White, not tested. *YUC5* is shown for reference. **b**, Clusters in which each DOF was nominated by the CellOracle out-degree rank change (Fig. 5a; blue circles, more central in WT; orange diamonds, more central in *fer-4*). Boxes mark nominations in clusters where the same DOF is differentially expressed (**a**). Nine nominations in four clusters with fewer than 5 cells in a sample are not shown (Source Data).
