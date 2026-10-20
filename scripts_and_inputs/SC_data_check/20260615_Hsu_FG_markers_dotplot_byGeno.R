#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 20260615_Hsu_FG_markers_dotplot_byGeno.R
# DotPlot of the 55 curated F+G markers across celltype.anno.Li, FACETED by
# genotype (WT top / fer-4 bottom). Revisions per user vs the pooled version:
#   * expression from the SCT assay (SCT @data), not RNA log-norm
#   * diverging color bar #bcd0f5 -> white(0) -> #9c1727 on SCALED avg expr
#     (z-score per gene computed JOINTLY across all celltype x genotype groups
#      so WT vs fer-4 colors are comparable between facets)
#   * no grid lines
#   * facet rows = genotype (WT, fer-4)
# dot SIZE = % cells expressing ; dot COLOR = scaled avg expression.
# Run on SLURM (loads the 19.5 GB integrated object).
# ---------------------------------------------------------------------------

rm(list = ls())
setwd("/resnick/groups/NolanLab/cchen/Feronia/analysis")
suppressMessages({ library(Seurat); library(tidyverse); library(scales) })
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")

read_seu <- function(dir) {
  seu <- readRDS(dir)
  seu@assays$spliced_RNA   <- NULL
  seu@assays$spliced_SCT   <- NULL
  seu@assays$unspliced_RNA <- NULL
  seu@assays$unspliced_SCT <- NULL
  seu <- suppressWarnings(suppressMessages(UpdateSlots(seu)))
  seu <- suppressWarnings(suppressMessages(UpdateSeuratObject(seu)))
  seu
}

cat("Reading integrated object ...\n")
seu <- read_seu("/resnick/groups/NolanLab/shared_resources/annotated_objects/rds/integrated/Feronia_salinity_timecourse_4S_20251226_seu3_sc291suspiciousClusterRemoved.rds")
seu_sub <- subset(seu, subset = orig.ident %in% c("sc_291", "sc_292", "sc_293", "sc_294"))
rm(seu); gc()

# --- attach genotype (geno) from the sample sheet ---------------------------
# Sample sheet for sc_291-sc_294, included in this repository; set FERONIA_REPO to your clone (default: current directory).
samples <- read.csv(file.path(Sys.getenv("FERONIA_REPO", "."), "scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/Nolan_single_cell_lab_samples_sc291-294.csv"))
sdict <- samples[samples[[1]] %in% c("sc_291","sc_292","sc_293","sc_294"), c(1, 6)]
colnames(sdict) <- c("sample", "genotype")
seu_sub$geno <- factor(sdict$genotype[match(as.character(seu_sub$orig.ident), sdict$sample)],
                       levels = c("WT", "fer-4"))
cat("geno x sample:\n"); print(table(seu_sub$orig.ident, seu_sub$geno, useNA = "ifany"))

# --- cell-type order --------------------------------------------------------
ct_col <- "celltype.anno.Li"
stopifnot(ct_col %in% colnames(seu_sub@meta.data))
ct_order <- c("Quiescent Center", "Columella", "Lateral Root Cap",
              "Atrichoblast", "Trichoblast", "Cortex", "Endodermis",
              "Xylem Pole Pericycle", "Phloem Pole Pericycle", "Procambium",
              "Protophloem", "Metaphloem & Companion Cell",
              "Protoxylem", "Metaxylem")
seu_sub$ct <- factor(as.character(seu_sub@meta.data[[ct_col]]), levels = ct_order)

# combined celltype x genotype identity (one DotPlot -> shared per-gene scaling)
seu_sub$ct_geno <- paste(as.character(seu_sub$ct), as.character(seu_sub$geno), sep = "@@")
seu_sub <- subset(seu_sub, subset = ct_geno %in%
                  apply(expand.grid(ct_order, c("WT","fer-4")), 1, paste, collapse = "@@"))
Idents(seu_sub) <- "ct_geno"
cat("n cells per celltype x genotype:\n"); print(table(seu_sub$ct, seu_sub$geno))

# --- markers (finalized F+G union, 55) --------------------------------------
m <- read.csv(file.path(Sys.getenv("FERONIA_REPO", "."), "scripts_and_inputs/SC_data_check/inputs/20260615_Hsu_FG_union_celltype_mapping.csv"),
              check.names = FALSE)
m$primary <- trimws(sub(";.*", "", m$proposed_celltype.anno.Li))
m$label   <- ifelse(m$Gene != m$Locus, paste0(m$Locus, " (", m$Gene, ")"), m$Locus)
m <- m %>% mutate(primary = factor(primary, levels = ct_order)) %>% arrange(primary, Gene)

DefaultAssay(seu_sub) <- "SCT"
present <- m$Locus %in% rownames(seu_sub[["SCT"]])
cat("Markers present in SCT assay:", sum(present), "/", nrow(m), "\n")
if (any(!present)) cat("MISSING from SCT (dropped):", paste(m$Locus[!present], collapse = ", "), "\n")
m_use  <- m[present, ]
genes  <- m_use$Locus
labvec <- setNames(m_use$label, m_use$Locus)

# --- DotPlot over the 28 ct x geno groups (SCT, scaled) ---------------------
dp <- DotPlot(seu_sub, features = genes, assay = "SCT", scale = TRUE)
d  <- dp$data
# split the combined id back into celltype + genotype
d$celltype <- factor(sub("@@.*", "", as.character(d$id)), levels = rev(ct_order)) # QC at top
d$geno     <- factor(sub(".*@@", "", as.character(d$id)), levels = c("WT", "fer-4"))
d$gene     <- factor(as.character(d$features.plot), levels = genes)

# --- plot -------------------------------------------------------------------
p <- ggplot(d, aes(x = gene, y = celltype)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  facet_wrap(vars(geno), ncol = 1, scales = "free_x") +   # repeat x axis on both panels
  scale_x_discrete(labels = labvec[genes]) +
  scale_color_gradientn(colours = c("#0d409e", "#bcd0f5", "white", "#9c1727"),
                        values = scales::rescale(c(-2.5, -1.0, 0, 2.5)),  # white at 0
                        limits = c(-2.5, 2.5), oob = scales::squish,
                        na.value = "grey90", name = "Avg expression\n(SCT, scaled)") +
  scale_size(range = c(0, 6), name = "% expressed") +
  labs(title = "Hsu 2022 F+G curated markers - SCT expression by cell type and genotype",
       x = "Marker gene  [AT ID (symbol)]", y = "Cell type (celltype.anno.Li)") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.y = element_text(size = 8.5),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey95", color = NA),
        plot.title = element_text(face = "bold", size = 12))

out_pdf <- file.path(Sys.getenv("OUT_DIR", "."), "20260615_Hsu_FG_markers_dotplot_byGeno.pdf")
w <- max(10, 1.8 + 0.21 * length(genes))
h <- max(9, 3.2 + 0.34 * length(ct_order) * 2)        # two stacked facets + repeated x axis
ggsave(out_pdf, p, width = w, height = h, limitsize = FALSE)
cat(sprintf("Saved dot plot (%.1f x %.1f in): %s\n", w, h, out_pdf))

d$gene_label <- labvec[as.character(d$features.plot)]
out_csv <- file.path(Sys.getenv("OUT_DIR", "."), "20260615_Hsu_FG_markers_dotplot_byGeno_data.csv")
write.csv(d, out_csv, row.names = FALSE)
cat("Saved dot-plot data:", out_csv, "\n")
cat("DONE\n")
