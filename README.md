# FERONIA-cell-type-specific

Scripts and data for the paper: **Cell type-specific gene manipulation uncovers FERONIA receptor kinase signaling specificity in Arabidopsis roots**

## Overview

This repository contains scripts and input/output files for the analysis of cell type-specific transcriptional responses in *fer-4* and wild-type (WT) Arabidopsis roots, using single-cell RNA-seq data.

## Repository Structure

All scripts, inputs, and outputs are organized under `scripts_and_inputs/`, which is divided into two main sections:

- **Preprocessing**: notebooks and input files for sample integration and metadata assignment.
- **Analysis**: scripts, inputs, and outputs for downstream analyses, including differential expression and CellOracle gene regulatory network (GRN) inference. Each analysis subdirectory contains `scripts/`, `inputs/`, and `outputs/` folders.

## Preprocessing

Located in `scripts_and_inputs/Preprocessing/`.

Notebooks perform:
- **Metadata update**: Assigns cell type labels and sample metadata to the integrated Seurat object (`202501226_feronia_salinity_time_course_update_metadata.ipynb`)
- **Reference-based integration**: Integrates sc291–sc294 samples using Seurat v3 reference integration (`202501226_Ref_integration_feronia_salinity_time_course_seu3.ipynb`)

**Key input:**
- `Protoplasting_DEgene_FC2_list.txt` — list of genes differentially expressed due to protoplasting stress (used to filter out protoplasting artifacts)

## Analysis

### Differential Expression Analysis

Located in `scripts_and_inputs/Analysis/Differential_expression_analysis/`.

Notebooks perform:
- **Pseudobulk DE analysis**: Cell type-specific differential expression between *fer-4* and WT under untreated conditions using `muscat` (`20260317_feronia_salinity_timecourse_muscat_newIntegratedSeurat3_NoSuspiciousCluster.ipynb`)
- **Hormone gene overlap**: Overlap of DE genes with curated hormone-responsive gene sets (`20260317_feronia_salinity_DE_hormone_geneoverlap_newIntegratedSeurat3_NoSuspiciousCluster.ipynb`)

**Key inputs:**
- `Feronia_salinity_timecourse_4S_20251226_seu3_sc291suspiciousClusterRemoved.rds` — integrated Seurat object (gitignored)
- `features.tsv` — gene feature list
- `Nolan_single_cell_lab_samples_sc291-294.csv` — sample metadata
- `all_bg.csv` — background gene list for enrichment
- Hormone gene sets: `ABA.csv`, `ACC.csv`, `BR.csv`, `CK.csv`, `GA.csv`, `IAA.csv`, `SA.csv`, `Hickman_2017_JA.csv`

**Key output:**
- `sig_de_fil.csv` — filtered table of significant DE genes across cell types

### CellOracle GRN Analysis

Located in `scripts_and_inputs/Analysis/CellOracle_GRN_analysis/`.

Scripts perform:
- **GRN inference**: Runs CellOracle gene regulatory network (GRN) inference per cell type in a SLURM array job (`*_slurmArray_run_celloracle_grn.py`)
- **TF selection**: Identifies WT-specific and *fer-4*-specific transcription factors from inferred GRNs (`*_select_TF_genes.ipynb`)

**Key inputs:**
- `features.tsv.gz` — gene feature list
- `Kay_TF_thalemine_annotations.csv` — Arabidopsis TF annotations from ThaleMine (Kay et al.)

**Key outputs** (per cell type: atrichoblast, columella, cortex, endodermis, lrc, pericycle, phloem, procambium, trichoblast, xylem):
- `{cell_type}_rawRNA_celloracle.h5ad` — CellOracle AnnData object (gitignored)
- `{cell_type}_rawRNA_celloracle_network_iGRN.csv` — inferred GRN edge table
- `{cell_type}_rawRNA_celloracle_gene_score_iGRN.csv` — gene-level GRN scores
- `{cell_type}_rawRNA.links.celloracle.links` — raw CellOracle links object
- `{cell_type}_rawRNA_p01_20k.links.celloracle.links` — filtered links (p < 0.01, top 20k)

Summary tables across all cell types:
- `*_CellOracle_all_metrics_rank_change_wide_TF_GE.txt` — TF rank change metrics (wide format)
- `*_CellOracle_WT_specific_TF.txt` — WT-specific TFs
- `*_CellOracle_fer4_specific_TF.txt` — *fer-4*-specific TFs
