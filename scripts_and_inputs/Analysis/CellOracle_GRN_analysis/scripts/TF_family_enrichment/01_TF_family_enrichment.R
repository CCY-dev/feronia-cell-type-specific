#!/usr/bin/env Rscript
## 01 — TF-family enrichment among genotype-specific CellOracle regulators (Fig X1, Supplementary Table).
##
## Genotype-specific TFs come from the CellOracle rank-change table: genes are ranked by out-degree centrality within
## each network; genes in the top 100 of at least one network (per cell type, either genotype) enter the ranking;
## Δrank = rank(fer-4) − rank(WT); below the 10th percentile → fer-4-specific, above the 90th → WT-specific; a TF is
## nominated if this holds in at least one cell type.
## Background: every TF that could have entered the ranking = TFs of the CellOracle TF list scored for out-degree in
## BOTH genotypes in at least one cell type AND assigned a family by PlantTFDB (1,292 TFs).
## Test: per family and genotype, one-sided Fisher's exact test (family × nominated), two-sided 95% CI, BH across all
## 110 tests (55 families × 2 genotypes). Families with >=10 eligible TFs are displayed (filter applied after BH).
suppressMessages({library(dplyr)})
source("00_paths.R")

fam  <- read.delim(PLANTTFDB, stringsAsFactors = FALSE)
famv <- setNames(fam$Family, toupper(fam$Gene_ID))          # first family per gene (isoforms collapsed)
kay  <- toupper(read.csv(KAY_TF, stringsAsFactors = FALSE)$GeneID)

## ---- genotype-specific TFs, and the Δrank thresholds ----
w <- read.delim(WIDE_TABLE, check.names = FALSE, stringsAsFactors = FALSE); w$gene <- toupper(w$gene)
RD  <- "fer4_untreated-WT_untreated: rank_diff"              # = rank(fer-4) − rank(WT)
od  <- w[w$metric == "degree_centrality_out", ]
q   <- quantile(od[[RD]][is.finite(od[[RD]])], c(.10, .90))  # over all out-degree rows, as in the pipeline
o   <- od[od$Is_TF == "Yes", ]
gWT <- unique(o$gene[o$Label == "WT_specific"]); gFER <- unique(o$gene[o$Label == "fer4_specific"])
stopifnot(identical(sort(gFER), sort(unique(o$gene[o[[RD]] < q[1]]))),    # the labels ARE the percentile calls
          identical(sort(gWT),  sort(unique(o$gene[o[[RD]] > q[2]]))))
cat(sprintf("Δrank thresholds: 10th %.0f, 90th %.0f | ranked TFs %d | WT-specific %d | fer-4-specific %d | both %d\n",
            q[1], q[2], length(unique(o$gene)), length(gWT), length(gFER), length(intersect(gWT, gFER))))

## ---- background: TFs scored in both genotypes in >=1 cell type, PlantTFDB-classifiable ----
pres <- bind_rows(lapply(Sys.glob(SCORE_GLOB), function(f) {
  x <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE); names(x)[1] <- "gene"
  x <- x[!is.na(x$degree_centrality_out) & x$cluster != "", c("gene", "cluster")]
  x$gene <- toupper(x$gene)
  x$geno <- ifelse(grepl("_WT_untreated$", x$cluster), "WT", ifelse(grepl("_fer4_untreated$", x$cluster), "fer4", NA))
  x$ctx  <- sub("_(WT|fer4)_untreated$", "", x$cluster); x[!is.na(x$geno), c("gene", "ctx", "geno")] }))
paired <- pres %>% distinct() %>% group_by(gene, ctx) %>% summarise(ok = all(c("WT", "fer4") %in% geno), .groups = "drop") %>% filter(ok)
scored <- intersect(unique(paired$gene), kay)
BG <- scored[!is.na(famv[scored])]
cat(sprintf("TF list %d | scored in both genotypes %d | + PlantTFDB family = background %d (%d DOFs, %d families)\n",
            length(unique(kay)), length(scored), length(BG), sum(famv[BG] == "Dof"), length(unique(famv[BG]))))
stopifnot(all(intersect(union(gWT, gFER), names(famv)[!is.na(famv)]) %in% BG))   # every testable nominee is in BG
write.csv(data.frame(gene = BG, family = unname(famv[BG])), file.path(OUT_DIR, "enrichment_background_TFs.csv"), row.names = FALSE)

## ---- Fisher per family × genotype, BH over all tests ----
fisher_fam <- function(nominated, geno) {
  fv <- famv[BG]; hit <- BG %in% nominated
  bind_rows(lapply(sort(unique(fv)), function(f) {
    inF <- fv == f
    m <- matrix(c(sum(inF & hit), sum(inF & !hit), sum(!inF & hit), sum(!inF & !hit)), 2)
    f1 <- fisher.test(m, alternative = "greater"); f2 <- fisher.test(m)
    data.frame(genotype = geno, family = f, nominated = m[1], eligible = m[1] + m[2],
               OR = unname(f1$estimate), lo = f2$conf.int[1], hi = f2$conf.int[2], p = f1$p.value) })) }
res <- bind_rows(fisher_fam(gWT, "WT"), fisher_fam(gFER, "fer-4")) %>% mutate(q = p.adjust(p, "BH"))
write.csv(res, file.path(OUT_DIR, "TF_family_enrichment_all_tests.csv"), row.names = FALSE)
d <- res %>% filter(family == "Dof")
cat(sprintf("%d tests | %d with q < 0.05\n", nrow(res), sum(res$q < 0.05)))
for (i in seq_len(nrow(d))) cat(sprintf("  DOF %-6s %d/%d  OR %.2f (%.2f-%.2f)  q = %.2g\n",
                                        d$genotype[i], d$nominated[i], d$eligible[i], d$OR[i], d$lo[i], d$hi[i], d$q[i]))
