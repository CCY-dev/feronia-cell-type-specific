# Step 1 of the post-rerun rebuild: rank-change wide table + WT-/fer-4-specific TF lists from CellOracle gene scores.
# A script version of the 20260317 select_TF_genes notebook (cells 7, 29, 35-101), same rules:
#   - per cell type (time.celltype.anno.Li without genotype) and per metric, rank genes within each network
#     (higher score = rank 1, ties "min");
#   - keep genes in the top 100 (slice_max, ties kept) of at least one of that cell type's networks;
#   - pair fer-4 with WT (inner join), rank_diff = rank(fer-4) - rank(WT); cell types with one genotype are skipped;
#   - thresholds = 10th / 90th percentile of rank_diff per metric over all rows of all cell types;
#     rank_diff < low -> fer4_specific, > high -> WT_specific;
#   - Is_TF = in the Kay TF list; DE columns joined from a muscat sig_de.csv by gene + cell type;
#   - specific-TF lists = degree_centrality_out, Is_TF, label.
# Usage: Rscript <this> <CellOracle output dir> <layout: flat|subfolders> <muscat sig_de.csv> <out dir> <file prefix>
suppressMessages({library(dplyr); library(data.table)})
a <- commandArgs(TRUE)
CO <- a[1]; LAYOUT <- a[2]; MUSCAT <- a[3]; OUT <- a[4]; PFX <- a[5]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
# NOTE: lab cluster path below -- change it to your own path
KAY <- "/resnick/groups/NolanLab/shared_resources/source_files/celloracle_inputs/Kay_TF_thalemine_annotations.csv"
CRUDE <- c("cortex", "endodermis", "trichoblast", "atrichoblast", "pericycle", "phloem", "xylem", "procambium", "lrc", "columella")
METRICS <- c("degree_centrality_all", "degree_centrality_in", "degree_centrality_out", "betweenness_centrality")
TOPN <- 100

long <- list()
for (cr in CRUDE) {
  f <- if (LAYOUT == "flat") file.path(CO, paste0(cr, "_rawRNA_celloracle_gene_score_iGRN.csv")) else
    file.path(CO, cr, paste0(cr, "_rawRNA_celloracle_gene_score_iGRN.csv"))
  gs <- read.csv(f, header = TRUE)
  colnames(gs)[1] <- "gene"
  gs$name <- paste0(gs$gene, "_", gsub("[[:punct:][:space:][:cntrl:]]", "", as.character(gs$name)))
  gs <- gs %>% mutate(cell_type = sub("_(WT|fer4)_.*", "", cluster),
                      group = sub(".*((WT|fer4)_.*)", "\\1", cluster))
  for (ct in unique(gs$cell_type)) {
    s0 <- gs[gs$cell_type == ct, ]
    if (length(unique(s0$group)) < 2) { message("skip (one genotype): ", ct); next }
    for (m in METRICS) {
      s <- s0 %>% select(gene, name, cluster, group, value = all_of(m)) %>%
        group_by(cluster) %>% mutate(gene_rank = rank(-value, ties.method = "min")) %>% ungroup()
      top <- s %>% group_by(cluster) %>% slice_max(order_by = value, n = TOPN, with_ties = TRUE) %>% ungroup()
      s <- s %>% filter(gene %in% unique(top$gene))
      d1 <- s %>% filter(group == "fer4_untreated") %>% select(gene, name, gene_rank)
      d2 <- s %>% filter(group == "WT_untreated") %>% select(gene, gene_rank)
      j <- inner_join(d1, d2, by = "gene", suffix = c("_1", "_2")) %>%
        mutate(rank_diff = gene_rank_1 - gene_rank_2, rank_diff_rank = rank(rank_diff, ties.method = "min"))
      long[[length(long) + 1]] <- j %>% mutate(crude_cell_type = cr, cell_type = ct, metric = m)
    }
  }
}
tb <- bind_rows(long)

thr <- tb %>% group_by(metric) %>%
  summarise(low = quantile(rank_diff, 0.10), high = quantile(rank_diff, 0.90), n = n(), .groups = "drop")
print(as.data.frame(thr))
tb <- tb %>% left_join(thr %>% select(metric, low, high), by = "metric") %>%
  mutate(Label = case_when(rank_diff < low ~ "fer4_specific", rank_diff > high ~ "WT_specific", TRUE ~ "None"),
         NA_group_count = 0) %>% select(-low, -high)

kay <- read.csv(KAY); colnames(kay) <- c("gene", "TF_name", "TF_description")
stopifnot(!anyDuplicated(kay$gene))
tb <- tb %>% left_join(kay, by = "gene") %>% mutate(Is_TF = ifelse(!is.na(TF_name), "Yes", "No"))

mu <- fread(MUSCAT, select = c(2, 3, 18, 14), data.table = FALSE)
colnames(mu) <- c("gene", "cell_type", "p_adj", "log2FC")
stopifnot(!anyDuplicated(mu[, c("gene", "cell_type")]))
tb <- tb %>% left_join(mu, by = c("gene", "cell_type")) %>%
  mutate(DE_label = ifelse(p_adj < 0.05 & abs(log2FC) > log2(1.5), "DE", "None"))
cnt <- tb %>% group_by(gene, metric) %>%
  summarise(fer4_specific_count = sum(Label == "fer4_specific"), WT_specific_count = sum(Label == "WT_specific"), .groups = "drop")
tb <- tb %>% left_join(cnt, by = c("gene", "metric"))

C <- "fer4_untreated-WT_untreated: "
wide <- tb %>% transmute(gene, name, crude_cell_type, cell_type, metric,
                         !!paste0(C, "gene_rank_1") := gene_rank_1, !!paste0(C, "gene_rank_2") := gene_rank_2,
                         !!paste0(C, "rank_diff") := rank_diff, !!paste0(C, "rank_diff_rank") := rank_diff_rank,
                         Label, fer4_specific_count, WT_specific_count, log2FC, p_adj, DE_label, NA_group_count,
                         Is_TF, TF_name, TF_description) %>%
  arrange(factor(metric, levels = sort(METRICS)), crude_cell_type, cell_type, `fer4_untreated-WT_untreated: rank_diff_rank`, gene)
write.table(wide, file.path(OUT, paste0(PFX, "_CellOracle_all_metrics_rank_change_wide_TF_GE.txt")), sep = "\t", row.names = FALSE, col.names = TRUE)
for (k in c("fer4", "WT")) {
  x <- wide %>% filter(metric == "degree_centrality_out", Is_TF == "Yes", Label == paste0(k, "_specific")) %>% arrange(cell_type, desc(log2FC))
  write.table(x, file.path(OUT, paste0(PFX, "_CellOracle_", k, "_specific_TF.txt")), sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)
  cat(k, "-specific: ", nrow(x), " rows, ", n_distinct(x$gene), " TFs, ", n_distinct(x$cell_type), " cell types\n", sep = "")
}
o <- wide %>% filter(metric == "degree_centrality_out")
cat("out-degree rows:", nrow(o), "| cell types:", n_distinct(o$cell_type), "| genes:", n_distinct(o$gene),
    "| TFs:", n_distinct(o$gene[o$Is_TF == "Yes"]), "\n")
