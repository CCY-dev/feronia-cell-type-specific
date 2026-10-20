#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 20260612_pseudobulk_sample_correlation.R
#
# Sample-wise correlation matrix from pseudobulk of the single-cell data.
#   - pseudobulk = sum of RNA raw counts across all cells of each sample (orig.ident)
#   - normalize to log2(CPM+1); keep genes with total count >= 10
#   - Pearson and Spearman correlation across genes -> sample x sample matrices
#   - clustered heatmaps annotated by genotype / rep
# Run in conda env `muscat` (Seurat 5) via SLURM.
# ---------------------------------------------------------------------------

rm(list = ls())
setwd("/resnick/groups/NolanLab/cchen/Feronia/analysis")
suppressMessages({
  library(tidyverse); library(Seurat); library(Matrix)
  library(ComplexHeatmap); library(circlize)
})
out_dir <- Sys.getenv("OUT_DIR", ".")   # outputs (was the analyst's cluster folder)

read_seu <- function(dir) {
  seu <- readRDS(dir)
  seu@assays$spliced_RNA <- NULL; seu@assays$spliced_SCT <- NULL
  seu@assays$unspliced_RNA <- NULL; seu@assays$unspliced_SCT <- NULL
  seu <- suppressWarnings(suppressMessages(UpdateSlots(seu)))
  seu <- suppressWarnings(suppressMessages(UpdateSeuratObject(seu)))
  seu
}

cat("Reading object...\n")
seu <- read_seu("/resnick/groups/NolanLab/shared_resources/annotated_objects/rds/integrated/Feronia_salinity_timecourse_4S_20251226_seu3_sc291suspiciousClusterRemoved.rds")
cat("Samples (orig.ident):\n"); print(table(seu$orig.ident))

# --- pseudobulk: sum RNA counts per sample (sparse matrix multiply) ---------
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")   # genes x cells (sparse)
samp   <- factor(seu$orig.ident)
mm     <- sparse.model.matrix(~ 0 + samp); colnames(mm) <- levels(samp)
pb     <- as.matrix(counts %*% mm)                              # genes x samples (summed counts)
cat("Pseudobulk matrix:", nrow(pb), "genes x", ncol(pb), "samples\n")

# --- normalize to log2 CPM, filter low genes --------------------------------
keep   <- rowSums(pb) >= 10
cpm    <- sweep(pb, 2, colSums(pb), "/") * 1e6
logcpm <- log2(cpm[keep, ] + 1)
cat("Genes kept (total count >= 10):", sum(keep), "\n")

# --- correlation matrices ---------------------------------------------------
pear  <- cor(logcpm, method = "pearson")
spear <- cor(logcpm, method = "spearman")
write.csv(round(pear, 4),  file.path(out_dir, "20260612_pseudobulk_sample_corr_pearson.csv"))
write.csv(round(spear, 4), file.path(out_dir, "20260612_pseudobulk_sample_corr_spearman.csv"))
cat("\nPearson correlation matrix:\n"); print(round(pear, 3))

# --- sample annotation (genotype / rep) from the sample sheet ---------------
# Sample sheet for sc_291-sc_294, included in this repository; set FERONIA_REPO to your clone (default: current directory).
samples <- suppressMessages(read_csv(
  file.path(Sys.getenv("FERONIA_REPO", "."), "scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/Nolan_single_cell_lab_samples_sc291-294.csv"),
  col_names = TRUE, show_col_types = FALSE))
sdict <- samples[match(colnames(pear), samples$sample), c(1, 6, 12)]
colnames(sdict) <- c("sample", "genotype", "rep")
ann_df <- data.frame(genotype = sdict$genotype, rep = as.character(sdict$rep))
rownames(ann_df) <- colnames(pear)

# --- annotated, clustered heatmaps ------------------------------------------
mk_ha <- function() {
  cols_g <- structure(scales::hue_pal()(length(unique(na.omit(ann_df$genotype)))),
                      names = unique(na.omit(ann_df$genotype)))
  cols_r <- structure(scales::hue_pal()(length(unique(na.omit(ann_df$rep)))),
                      names = unique(na.omit(ann_df$rep)))
  HeatmapAnnotation(genotype = ann_df$genotype, rep = ann_df$rep,
                    col = list(genotype = cols_g, rep = cols_r))
}
draw_hm <- function(m, title) {
  cf <- colorRamp2(c(min(m), (min(m) + 1) / 2, 1), c("#2166AC", "#F7F7F7", "#B2182B"))
  Heatmap(m, name = "r", top_annotation = mk_ha(), col = cf,
          column_title = title,
          cell_fun = function(j, i, x, y, w, h, fill)
            grid::grid.text(sprintf("%.3f", m[i, j]), x, y, gp = grid::gpar(fontsize = 8)),
          clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean")
}
pdf(file.path(out_dir, "20260612_pseudobulk_sample_corr_heatmap.pdf"), width = 8, height = 7)
draw(draw_hm(pear,  "Pseudobulk sample correlation (Pearson, log2 CPM)"))
draw(draw_hm(spear, "Pseudobulk sample correlation (Spearman, log2 CPM)"))
dev.off()

cat("\nSaved CSVs + heatmap PDF to", out_dir, "\nDONE\n")
