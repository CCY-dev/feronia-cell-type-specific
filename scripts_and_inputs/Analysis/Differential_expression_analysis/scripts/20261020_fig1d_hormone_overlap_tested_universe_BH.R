# Fig. 1d: overlap of each subpopulation's DEGs with eight hormone-response gene sets (20 Oct 2026).
#
# Statistic (decided by Trevor Nolan, 20 Oct 2026):
#   universe  = the 17,554 genes tested by muscat in any subpopulation (unique genes in inputs/all_bg.csv);
#   gene sets = each hormone list restricted to that universe (genes that could never be DEGs are removed);
#   test      = GeneOverlap one-sided Fisher's exact test (as in the hormone-overlap notebook);
#   q         = Benjamini-Hochberg over all 272 tests (34 subpopulations x 8 hormone sets).
# The published notebook (20260317_..._hormone_geneoverlap_...ipynb) used the same universe size but did not restrict
# the lists to it and reported unadjusted P; it is kept unchanged as the record of that version.
#
# Run from anywhere:  Rscript 20261020_fig1d_hormone_overlap_tested_universe_BH.R
#   FERONIA_REPO = repository root (default: three folders above this script's folder, or the current folder)
#   OUT_DIR      = output folder (default: <repo>/scripts_and_inputs/Analysis/Differential_expression_analysis/outputs/hormone_overlap)
suppressPackageStartupMessages(library(GeneOverlap))   # recorded version 1.38.0 (env muscat)

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1])))
repo <- Sys.getenv("FERONIA_REPO", unset = normalizePath(file.path(here, "..", "..", "..", "..")))
de   <- file.path(repo, "scripts_and_inputs", "Analysis", "Differential_expression_analysis")
out  <- Sys.getenv("OUT_DIR", unset = file.path(de, "outputs", "hormone_overlap"))
dir.create(out, showWarnings = FALSE, recursive = TRUE)

sig  <- read.delim(file.path(de, "outputs", "sig_de_fil.csv"), sep = ",", check.names = FALSE)
bg   <- read.delim(file.path(de, "inputs", "all_bg.csv"), sep = ",", check.names = FALSE)
universe <- unique(bg$gene)
stopifnot(length(universe) == 17554)

sets <- c(IAA = "IAA.csv", BR = "BR.csv", JA = "Hickman_2017_JA.csv", ABA = "ABA.csv",
          CK = "CK.csv", ACC = "ACC.csv", SA = "SA.csv", GA = "GA.csv")
read_set <- function(f) trimws(read.csv(file.path(de, "inputs", f), fileEncoding = "UTF-8-BOM")$GeneID)
hormone <- lapply(sets, read_set)
hormone_in <- lapply(hormone, function(g) intersect(unique(g), universe))

clusters <- sort(unique(sig$cluster_id))
degs <- lapply(setNames(clusters, clusters), function(cl) unique(sig$gene[sig$cluster_id == cl]))
stopifnot(all(unlist(degs) %in% universe), length(clusters) == 34)

rows <- list()
for (cl in clusters) for (h in names(sets)) {
  go <- testGeneOverlap(newGeneOverlap(degs[[cl]], hormone_in[[h]], genome.size = length(universe)))
  rows[[length(rows) + 1]] <- data.frame(subpopulation = cl, hormone = h, overlap = length(go@intersection),
    n_deg = length(degs[[cl]]), n_universe = length(universe), n_set = length(unique(hormone[[h]])),
    n_set_in_universe = length(hormone_in[[h]]), p = go@pval)
}
res <- do.call(rbind, rows)
res$q <- p.adjust(res$p, method = "BH")
write.csv(res, file.path(out, "Fig1d_overlap_q.csv"), row.names = FALSE)
cat(sprintf("%d tests; q < 0.05: %d; P < 0.05: %d\n", nrow(res), sum(res$q < 0.05), sum(res$p < 0.05)))
writeLines(capture.output(sessionInfo()), file.path(out, "Fig1d_overlap_q.sessionInfo.txt"))
