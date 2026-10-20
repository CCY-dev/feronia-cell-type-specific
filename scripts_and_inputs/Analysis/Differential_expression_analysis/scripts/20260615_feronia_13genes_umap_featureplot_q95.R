#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 20260615_feronia_13genes_umap_featureplot_q95.R
#
# REVISION of 20260610_feronia_13genes_umap_featureplot.R per user request:
#   1. max.cutoff = "q95"  -> cap color scale at the 95th percentile so a few
#      outlier high-expressing cells stop washing out the rest of the UMAP.
#   2. pt.size = 1.25x Seurat's auto point size -> larger dots.
#   3. order = TRUE -> plot high-expressing cells LAST so they sit on top of
#      the low/zero cells instead of being hidden underneath.
#
# Outputs go to NEW filenames (..._q95.pdf) so the originals are preserved
# for side-by-side comparison.
#
# Everything else (data prep, gene list, titles, layout) is unchanged from
# the 20260610 script.
#
# Run inside conda env `muscat` (Seurat 5) on a SLURM compute node.
# ---------------------------------------------------------------------------

rm(list = ls())
setwd("/resnick/groups/NolanLab/cchen/Feronia/analysis")

suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
})

cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")

# --- read Seurat 3/4 object inside a Seurat 5 env (from the notebook) --------
read_seu <- function(dir) {
  seu <- readRDS(dir)
  seu@assays$spliced_RNA   <- NULL
  seu@assays$spliced_SCT   <- NULL
  seu@assays$unspliced_RNA <- NULL
  seu@assays$unspliced_SCT <- NULL
  seu <- suppressWarnings(suppressMessages(UpdateSlots(seu)))
  seu <- suppressWarnings(suppressMessages(UpdateSeuratObject(seu)))
  return(seu)
}

cat("Reading integrated Seurat object ...\n")
seu <- read_seu("/resnick/groups/NolanLab/shared_resources/annotated_objects/rds/integrated/Feronia_salinity_timecourse_4S_20251226_seu3_sc291suspiciousClusterRemoved.rds")

# --- subset to the 4 timecourse samples -------------------------------------
seu_sub <- subset(seu, subset = orig.ident %in% c("sc_291", "sc_292", "sc_293", "sc_294"))
rm(seu); gc()

# --- attach genotype (geno) from the sample sheet ---------------------------
# Sample sheet for sc_291-sc_294, included in this repository; set FERONIA_REPO to your clone (default: current directory).
samples <- read_csv(file.path(Sys.getenv("FERONIA_REPO", "."), "scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/Nolan_single_cell_lab_samples_sc291-294.csv"),
                    col_names = TRUE, show_col_types = FALSE)
sdict <- samples[samples$sample %in% c("sc_291", "sc_292", "sc_293", "sc_294"), c(1, 6, 12)]

oi <- as_tibble(seu_sub$orig.ident)
colnames(oi) <- "sample"
joined_data <- left_join(oi, sdict, by = "sample")

seu_sub$geno <- joined_data$genotype
seu_sub$date <- joined_data$rep

# order panels WT (left) then fer-4 (right)
seu_sub$geno <- factor(seu_sub$geno, levels = c("WT", "fer-4"))

DefaultAssay(seu_sub) <- "SCT"

cat("Reductions available:", paste(Reductions(seu_sub), collapse = ", "), "\n")
cat("geno counts:\n"); print(table(seu_sub$geno, useNA = "ifany"))

# --- point size: 1.25x Seurat's AutoPointSize -------------------------------
# Seurat sizes each split panel by AutoPointSize() = min(1583 / n_cells, 1).
# With split.by the relevant n is the per-genotype cell count; we use the
# denser (larger) panel so neither panel gets over-crowded, then x1.25.
# A single explicit pt.size keeps WT and fer-4 dots identical for comparison.
geno_n   <- table(seu_sub$geno)
n_panel  <- max(geno_n)
auto_pt  <- min(1583 / n_panel, 1)
pt_size  <- 1.25 * auto_pt
cat(sprintf("Per-genotype cell counts: %s\n",
            paste(names(geno_n), geno_n, sep = "=", collapse = ", ")))
cat(sprintf("AutoPointSize(n=%d) = %.4f ; pt.size (x1.25) = %.4f\n",
            n_panel, auto_pt, pt_size))

# --- 13 genes from Slide3 Figure A, grouped by direction --------------------
gene_groups <- list(
  upregulated   = c("AT3G05170", "AT2G48090", "AT3G10930", "AT5G22330", "AT3G55970", "AT2G29510", "AT3G30720"),
  downregulated = c("AT5G46295", "AT3G10910", "AT3G01420", "AT3G59370", "AT3G50640", "AT2G37750")
)
selected_genes <- unlist(gene_groups, use.names = FALSE)
cat("Number of genes:", length(selected_genes),
    "(", length(gene_groups$upregulated), "up,", length(gene_groups$downregulated), "down )\n")

# --- common names from features.tsv (falls back to gene ID if no symbol) ----
feature_names <- read.delim("/resnick/groups/NolanLab/shared_resources/source_files/features.tsv",
                            sep = "\t", header = FALSE)
feature_names <- feature_names[, 1:2]
colnames(feature_names) <- c("gene", "common_name")

gene_map_df <- feature_names[match(selected_genes, feature_names$gene), ]
miss_name <- is.na(gene_map_df$gene)
gene_map_df$gene[miss_name]        <- selected_genes[miss_name]
gene_map_df$common_name[is.na(gene_map_df$common_name)] <-
  gene_map_df$gene[is.na(gene_map_df$common_name)]

# --- which genes are actually present in the SCT assay -----------------------
present <- selected_genes %in% rownames(seu_sub)
cat("Genes present in SCT assay:", sum(present), "/", length(selected_genes), "\n")
if (any(!present)) cat("MISSING from SCT assay:", paste(selected_genes[!present], collapse = ", "), "\n")

# --- helper: render one gene's UMAP page (or a placeholder) -----------------
plot_gene_page <- function(current_gene) {

  # common name as the page title (gene ID if no symbol)
  display_name <- gene_map_df$common_name[match(current_gene, gene_map_df$gene)]
  full_title   <- display_name

  if (current_gene %in% rownames(seu_sub)) {

    # ---- FeaturePlot block: q95 cutoff, larger dots, high expr on top ----
    p <- FeaturePlot(
      seu_sub, features = current_gene,
      reduction  = "umap",
      split.by   = "geno",
      ncol       = 2,
      max.cutoff = "q95",     # cap color scale at 95th percentile
      order      = TRUE,      # high-expressing cells plotted on top
      pt.size    = pt_size    # 1.25x Seurat auto size
    ) + theme(legend.position = c(0.9, 0.2))

    p <- p + plot_annotation(
      title = full_title,
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
    print(p)
    cat("  Plotted:", current_gene, "(", display_name, ")\n")

  } else {
    # placeholder page so absent genes are still represented
    pmsg <- ggplot() + theme_void() +
      annotate("text", x = 0, y = 0.15,
               label = paste0(full_title, " (", current_gene, ")"),
               size = 6, fontface = "bold") +
      annotate("text", x = 0, y = -0.1,
               label = "Not detected in SCT assay - no UMAP available",
               size = 5, color = "red")
    print(pmsg)
    cat("  Skipped (absent from SCT):", current_gene, "\n")
  }
}

# --- one multi-page PDF per direction (up / down) ---------------------------
out_dir <- Sys.getenv("OUT_DIR", ".")   # outputs (was the analyst's cluster folder)
out_files <- c(
  upregulated   = file.path(out_dir, "20260615_feronia_umap_upregulated_q95.pdf"),
  downregulated = file.path(out_dir, "20260615_feronia_umap_downregulated_q95.pdf")
)

for (grp in names(gene_groups)) {
  out_pdf <- out_files[[grp]]
  cat("== ", grp, " (", length(gene_groups[[grp]]), " genes) -> ", out_pdf, " ==\n", sep = "")
  # width/height match options(repr.plot.width = 10, repr.plot.height = 5)
  pdf(out_pdf, width = 10, height = 5)
  for (current_gene in gene_groups[[grp]]) plot_gene_page(current_gene)
  dev.off()
  cat("Saved PDF:", out_pdf, "\n")
}
cat("DONE\n")
