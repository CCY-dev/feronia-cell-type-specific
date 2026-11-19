# Fig. 1a source extraction (read-only on the shared object).
# Object: the one read by the repository muscat notebook
#   scripts_and_inputs/Analysis/Differential_expression_analysis/scripts/
#   20260317_feronia_salinity_timecourse_muscat_newIntegratedSeurat3_NoSuspiciousCluster.ipynb (cell 7)
# Writes per-cell UMAP coordinates, sample, annotations and SCT-normalised FER (AT3G51550) to
# fig1a_umap_cells.csv.gz, plus fig1a_object_summary.txt. Slots are read directly (no object update),
# so values are exactly those stored in the object.
suppressMessages(library(Seurat))
# NOTE: lab cluster path below -- change it to your own path
f <- "/resnick/groups/NolanLab/shared_resources/annotated_objects/rds/integrated/Feronia_salinity_timecourse_4S_20251226_seu3_sc291suspiciousClusterRemoved.rds"
t0 <- Sys.time()
seu <- readRDS(f)
sink("fig1a_object_summary.txt")
cat("object:", f, "\n"); cat("file size:", file.size(f), "\n")
md <- seu@meta.data   # Seurat v3 object: read slots directly (generics such as ncol() need UpdateSlots)
cat("cells:", nrow(md), "\n")
cat("assays:", paste(names(seu@assays), collapse = ", "), " active:", seu@active.assay, "\n")
for (a in names(seu@assays)) cat("  assay", a, "data dims:", dim(seu@assays[[a]]@data), "\n")
cat("reductions:", paste(names(seu@reductions), collapse = ", "), "\n")
cat("meta columns:\n"); print(colnames(md))
for (cn in colnames(md)) {
  v <- md[[cn]]
  if ((is.character(v) || is.factor(v)) && length(unique(v)) <= 60) {
    cat("\n== ", cn, " (", length(unique(v)), " levels)\n", sep = "")
    if (is.factor(v)) cat("factor levels:", paste(levels(v), collapse = " | "), "\n")
    print(table(v, md$orig.ident))
  }
}
for (r in names(seu@reductions)) { e <- seu@reductions[[r]]@cell.embeddings; cat("\nreduction", r, "dims", dim(e), "first cols", colnames(e)[1:min(3, ncol(e))], "\n") }
sct <- seu@assays$SCT
fer <- as.numeric(sct@data["AT3G51550", rownames(md)])
cat("\nFER SCT data: min", min(fer), "max", max(fer), "q95", quantile(fer, 0.95), "q99", quantile(fer, 0.99), "\n")
print(tapply(fer, md$orig.ident, function(x) c(mean = mean(x), frac_pos = mean(x > 0))))
sink()
keep <- names(Filter(function(v) (is.character(v) || is.factor(v)) && length(unique(v)) <= 60, md))
out <- data.frame(cell = rownames(md), md[, keep, drop = FALSE], stringsAsFactors = FALSE)
for (r in intersect(c("umap", "umap_2D"), names(seu@reductions))) {
  e <- seu@reductions[[r]]@cell.embeddings[rownames(md), 1:2]
  out[[paste0(r, "_1")]] <- e[, 1]; out[[paste0(r, "_2")]] <- e[, 2]
}
out$FER_SCT <- fer
write.csv(out, gzfile("fig1a_umap_cells.csv.gz"), row.names = FALSE)
cat("done in", format(Sys.time() - t0), "\n")
