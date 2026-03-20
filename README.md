# FERONIA-cell-type-specific

Scripts and data for the paper: **Cell type-specific gene manipulation uncovers FERONIA receptor kinase signaling specificity in Arabidopsis roots**

## Overview

This repository contains scripts and input/output files for the analysis of cell type-specific transcriptional responses in *fer-4* and wild-type (WT) Arabidopsis roots, using single-cell RNA-seq data.

## Repository Structure

All scripts, inputs, and outputs are organized under `scripts_and_inputs/`, which is divided into two main sections:

- **Preprocessing**: notebooks and input files for sample annotation and integration. Note: the raw Seurat objects (`sc_291`–`sc_294` COPILOT seu3 `.rds` files) are too large to include here and are deposited in the GEO repository instead.
- **Analysis**: scripts, inputs, and outputs for downstream analyses, including differential expression and CellOracle gene regulatory network (GRN) inference.

