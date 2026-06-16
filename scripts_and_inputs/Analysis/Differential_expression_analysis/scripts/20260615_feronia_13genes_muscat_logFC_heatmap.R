#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 20260615_feronia_13genes_muscat_logFC_heatmap.R
#
# Heatmap of the 13 Slide3 Fig A DEGs (genes x cell-type cluster) using the
# muscat per-cluster logFC (contrast fer_4-WT). Companion to the UMAP
# FeaturePlots. Design choices (discussed with user):
#   * value = muscat logFC (NOT recomputed expression) -> consistent with DE.
#     Sign convention: logFC>0 = up in fer-4 (red), <0 = down in fer-4 (blue).
#   * source = sig_de.csv = the FULL muscat result (all tested gene x cluster,
#     718502 rows, sig + non-sig) so we can tell "tested but not significant"
#     (colored, no asterisk) apart from "not tested in that cluster" (grey/NA).
#     NB: all_bg.csv is the expression-filtered GO *background* (a subset), NOT
#     the tested universe -- using it wrongly greyed out tested cells (e.g.
#     DOX1 x Distal Columella, which is significant). Always use sig_de.csv here.
#   * significance overlay: asterisks from p_adj.loc (*<0.05 **<0.01 ***<0.001).
#   * diverging scale centered at 0, clipped symmetrically at the 95th pct of
#     |logFC| so a few large-logFC restricted genes don't wash out the rest.
#   * rows split (faceted) into Up- vs Down-in-fer-4 blocks, Slide3 order.
#
# No Seurat / object load needed -- pure table -> plot. Still run via SLURM.
# ---------------------------------------------------------------------------

rm(list = ls())

suppressMessages({
  library(tidyverse)
  library(scales)   # squish for out-of-bounds color clipping
})

# --- inputs -----------------------------------------------------------------
de_path  <- "/resnick/groups/NolanLab/cchen/Feronia/analysis/muscat_result/20260317_salinity_timecourse_useDateT_noSuspiciousClusters/group_timecelltypeannoLi/fer_4-WT/sig_de.csv"
feat_tsv <- "/resnick/groups/NolanLab/shared_resources/source_files/features.tsv"
out_dir  <- "/resnick/groups/NolanLab/cchen/Feronia/claude_result"

pdf_path    <- file.path(out_dir, "20260615_feronia_13genes_muscat_logFC_heatmap.pdf")
matrix_csv  <- file.path(out_dir, "20260615_feronia_13genes_muscat_logFC_matrix.csv")
long_csv    <- file.path(out_dir, "20260615_feronia_13genes_muscat_logFC_long.csv")

# --- the 13 genes, grouped by Slide3 direction (top->bottom order) ----------
up_genes   <- c("AT3G05170", "AT2G48090", "AT3G10930", "AT5G22330", "AT3G55970", "AT2G29510", "AT3G30720")
down_genes <- c("AT5G46295", "AT3G10910", "AT3G01420", "AT3G59370", "AT3G50640", "AT2G37750")
sel_genes  <- c(up_genes, down_genes)

gene_dir <- tibble(
  gene = sel_genes,
  direction = c(rep("Up in fer-4", length(up_genes)),
                rep("Down in fer-4", length(down_genes)))
)

# --- gene-symbol labels (features.tsv: AGI id, common name) -----------------
feat <- read.delim(feat_tsv, sep = "\t", header = FALSE)[, 1:2]
colnames(feat) <- c("gene", "symbol")
lab_map <- gene_dir %>%
  left_join(feat, by = "gene") %>%
  mutate(label = ifelse(!is.na(symbol) & symbol != "" & symbol != gene,
                        paste0(symbol, " (", gene, ")"), gene)) %>%
  select(gene, label, direction)

# y-axis order: first gene (up_genes[1]) at TOP -> reverse factor levels
label_top_to_bottom <- lab_map$label[match(sel_genes, lab_map$gene)]

# --- read muscat DE table (only needed cols) --------------------------------
needed <- c("gene", "cluster_id", "logFC", "p_adj.loc")
cat("Reading muscat sig_de.csv (full tested set) ...\n")
if (requireNamespace("data.table", quietly = TRUE)) {
  de <- as.data.frame(data.table::fread(de_path, select = needed, showProgress = FALSE))
  cat("  (read via data.table::fread)\n")
} else {
  de <- read.csv(de_path, header = TRUE, check.names = FALSE)[, needed]
  cat("  (read via read.csv)\n")
}
cat("  sig_de rows:", nrow(de), "\n")

de <- de %>% filter(gene %in% sel_genes)
cat("  rows for the 13 genes:", nrow(de), "\n")

clusters <- sort(unique(de$cluster_id))   # clusters where >=1 of the 13 tested
cat("  cell-type clusters with >=1 of the 13 genes tested:", length(clusters), "\n")
print(clusters)

# --- complete the gene x cluster grid (missing = not tested = NA -> grey) ----
df <- expand_grid(gene = sel_genes, cluster_id = clusters) %>%
  left_join(de, by = c("gene", "cluster_id")) %>%
  left_join(lab_map, by = "gene") %>%
  mutate(
    label      = factor(label, levels = rev(label_top_to_bottom)),
    direction  = factor(direction, levels = c("Up in fer-4", "Down in fer-4")),
    cluster_id = factor(cluster_id, levels = clusters),
    sig        = cut(p_adj.loc, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                     labels = c("***", "**", "*", "")),
    sig        = ifelse(is.na(p_adj.loc), "", as.character(sig))
  )

n_tested <- sum(!is.na(df$logFC))
n_sig    <- sum(df$sig != "" , na.rm = TRUE)
cat(sprintf("Grid cells: %d (genes %d x clusters %d); tested(colored)=%d; significant=%d\n",
            nrow(df), length(sel_genes), length(clusters), n_tested, n_sig))

# --- symmetric color cap (95th pct of |logFC|) ------------------------------
cap <- ceiling(quantile(abs(df$logFC), 0.95, na.rm = TRUE))
cap <- max(cap, 1)
cat(sprintf("Color scale clipped symmetrically at +/-%g (95th pct |logFC|)\n", cap))

# --- plot -------------------------------------------------------------------
p <- ggplot(df, aes(cluster_id, label, fill = logFC)) +
  geom_tile(color = "grey92", linewidth = 0.2) +
  geom_text(aes(label = sig), size = 2.6, color = "black", na.rm = TRUE) +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-cap, cap),
                       oob = scales::squish, na.value = "grey85",
                       name = "logFC\n(fer-4 vs WT)") +
  labs(
    title = "fer-4 vs WT logFC - 13 Slide3 DEGs (muscat, per cell-type cluster)",
    x = "Cell-type cluster (time.celltype.anno.Li)", y = NULL,
    caption = sprintf(paste0(
      "Color = muscat logFC (red = up in fer-4, blue = down in fer-4). ",
      "Grey = gene not tested in that cluster. ",
      "Asterisks = p_adj.loc (* <0.05, ** <0.01, *** <0.001). Scale clipped at +/-%g."), cap)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y  = element_text(size = 8),
    panel.grid   = element_blank(),
    strip.text.y = element_blank(),
    plot.title   = element_text(face = "bold", size = 12),
    plot.caption = element_text(hjust = 0, size = 7),
    legend.position = "right"
  )

w <- max(7, 2.6 + 0.30 * length(clusters))
h <- max(6, 2.6 + 0.34 * length(sel_genes))
ggsave(pdf_path, p, width = w, height = h, limitsize = FALSE)
cat(sprintf("Saved heatmap (%.1f x %.1f in): %s\n", w, h, pdf_path))

# --- companion CSVs (the numbers behind the figure) -------------------------
wide <- df %>%
  mutate(label = factor(label, levels = label_top_to_bottom)) %>%  # natural order in CSV
  arrange(label) %>%
  select(label, cluster_id, logFC) %>%
  pivot_wider(names_from = cluster_id, values_from = logFC)
write.csv(wide, matrix_csv, row.names = FALSE)
cat("Saved logFC matrix:", matrix_csv, "\n")

long_out <- df %>%
  select(gene, label, direction, cluster_id, logFC, p_adj.loc, sig) %>%
  arrange(direction, factor(label, levels = label_top_to_bottom), cluster_id)
write.csv(long_out, long_csv, row.names = FALSE)
cat("Saved long table:", long_csv, "\n")

cat("DONE\n")
