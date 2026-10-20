# FER Extended Data Fig. 7d,e and Supplementary Figs 1-2: per-cell export (read-only on the shared object).
# Pattern reused from the figure kit's tools/hpc/extract_fig1a_umap.R: Seurat v3 slots are read directly (no
# UpdateSeuratObject), so every value is exactly what is stored in the object.
# Output: si_expression_cells.csv.gz (one row per cell of sc_291..sc_294) with
#   cell, orig.ident, genotype, celltype.anno.Li, time.anno.Li, time.celltype.anno.Li, umap_1, umap_2,
#   then SCT data-slot values for the 11 genes (columns named <SYMBOL>_<AGI>),
# plus si_expression_object_summary.txt and si_expression_sessionInfo.txt.
suppressMessages(library(Seurat))
# NOTE: lab cluster path below -- change it to your own path
f <- "/resnick/groups/NolanLab/shared_resources/annotated_objects/rds/integrated/Feronia_salinity_timecourse_4S_20251226_seu3_sc291suspiciousClusterRemoved.rds"
samples <- c("sc_291", "sc_292", "sc_293", "sc_294")
genotype_of <- c(sc_291 = "WT", sc_293 = "WT", sc_292 = "fer4", sc_294 = "fer4")   # task brief; checked below
genes <- c(YUC5 = "AT5G43890", DOF4.6 = "AT4G24060", CDF4 = "AT2G34140", COG1 = "AT1G29160", OBP3 = "AT3G55370",
           DOF6 = "AT3G45610", DOF4.1 = "AT4G00940", DOF2.2 = "AT2G28810", DAG1 = "AT3G61850", OBP1 = "AT3G50410",
           OBP2 = "AT1G07640")
t0 <- Sys.time()
seu <- readRDS(f)
sink("si_expression_object_summary.txt")
cat("object:", f, "\n"); cat("file size (bytes):", file.size(f), "\n"); cat("mtime:", format(file.mtime(f)), "\n")
cat("object class:", class(seu), " version slot:", as.character(seu@version), "\n")
md_all <- seu@meta.data
cat("cells in object:", nrow(md_all), "\n"); print(table(md_all$orig.ident))
cat("assays:", paste(names(seu@assays), collapse = ", "), " active:", seu@active.assay, "\n")
for (a in names(seu@assays)) cat("  assay", a, "data dims:", dim(seu@assays[[a]]@data), "\n")
cat("reductions:", paste(names(seu@reductions), collapse = ", "), "\n")
for (r in names(seu@reductions)) { e <- seu@reductions[[r]]@cell.embeddings
  cat("  reduction", r, "dims", dim(e), "cols", colnames(e)[1:min(3, ncol(e))], "\n") }
cat("meta columns:\n"); print(colnames(md_all))
md <- md_all[md_all$orig.ident %in% samples, , drop = FALSE]
cat("\ncells kept (sc_291..sc_294):", nrow(md), "\n")
gt <- unname(genotype_of[as.character(md$orig.ident)])
if ("genotype" %in% colnames(md)) {
  cat("\nmeta 'genotype' vs orig.ident:\n"); print(table(md$genotype, md$orig.ident))
  stopifnot(all(as.character(md$genotype) == gt))
  cat("meta genotype agrees with sc_291/sc_293 = WT, sc_292/sc_294 = fer4 for every kept cell\n")
} else cat("no meta 'genotype' column; genotype assigned from orig.ident\n")
if (!"time.celltype.anno.Li" %in% colnames(md)) {
  cat("\nNOTE: no meta column time.celltype.anno.Li; built here as paste(time.anno.Li, celltype.anno.Li, sep = '_')\n")
  md$time.celltype.anno.Li <- paste(md$time.anno.Li, md$celltype.anno.Li, sep = "_")
}
for (cn in c("celltype.anno.Li", "time.anno.Li", "time.celltype.anno.Li")) {
  stopifnot(cn %in% colnames(md))
  cat("\n== ", cn, " (", length(unique(md[[cn]])), " levels)\n", sep = "")
  if (is.factor(md[[cn]])) cat("factor levels:", paste(levels(md[[cn]]), collapse = " | "), "\n")
  print(table(md[[cn]], gt))
}
sct <- seu@assays$SCT
cat("\nSCT data class:", class(sct@data), "\n")
miss <- setdiff(genes, rownames(sct@data))
cat("genes missing from SCT data rows:", if (length(miss)) paste(miss, collapse = ", ") else "none", "\n")
out <- data.frame(cell = rownames(md), orig.ident = as.character(md$orig.ident), genotype = gt,
                  celltype.anno.Li = as.character(md$celltype.anno.Li), time.anno.Li = as.character(md$time.anno.Li),
                  time.celltype.anno.Li = as.character(md$time.celltype.anno.Li), stringsAsFactors = FALSE)
e <- seu@reductions$umap@cell.embeddings[rownames(md), 1:2]
out$umap_1 <- e[, 1]; out$umap_2 <- e[, 2]
for (r in setdiff(names(seu@reductions), c("umap", "pca"))) {   # any other 2-D embedding, for the deck comparison
  em <- seu@reductions[[r]]@cell.embeddings
  if (ncol(em) == 2) { out[[paste0(r, "_1")]] <- em[rownames(md), 1]; out[[paste0(r, "_2")]] <- em[rownames(md), 2] }
}
cat("\nSCT data-slot summary per gene (kept cells):\n")
for (s in names(genes)) {
  id <- genes[[s]]; col <- paste0(s, "_", id)
  v <- if (id %in% rownames(sct@data)) as.numeric(sct@data[id, rownames(md)]) else rep(NA_real_, nrow(md))
  out[[col]] <- v
  cat(sprintf("  %-7s %s  min %.4f  max %.4f  mean %.4f  frac>0 %.4f | WT mean %.4f frac>0 %.4f | fer4 mean %.4f frac>0 %.4f\n",
              s, id, min(v), max(v), mean(v), mean(v > 0), mean(v[gt == "WT"]), mean(v[gt == "WT"] > 0),
              mean(v[gt == "fer4"]), mean(v[gt == "fer4"] > 0)))
}
fer <- as.numeric(sct@data["AT3G51550", rownames(md)])
cat("\nFER (AT3G51550) SCT data, for comparison with the kit's fig1a export: min", min(fer), "max", max(fer),
    "sum", sum(fer), "\n")
sink()
write.csv(out, gzfile("si_expression_cells.csv.gz"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), "si_expression_sessionInfo.txt")
cat("rows written:", nrow(out), " columns:", ncol(out), "\n")
cat("done in", format(Sys.time() - t0), "\n")
