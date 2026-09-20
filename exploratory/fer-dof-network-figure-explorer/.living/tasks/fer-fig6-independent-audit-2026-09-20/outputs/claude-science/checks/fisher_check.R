#!/usr/bin/env Rscript
# Prototype; not yet carefully vetted.
# Recompute odds ratios, 95% CIs and one-sided p for the frozen top-ten family table with fisher.test.
# Read-only. Usage: Rscript fisher_check.R [repo_root]
a <- commandArgs(TRUE); root <- if (length(a)) a[1] else '/Users/tmn23/projects/fer-sc-manuscript-2026'
z <- read.csv(file.path(root,'revisions/2026-09-20_compact_figures/review_v2/quantitative/source_top10_family.csv'))
res <- do.call(rbind, lapply(seq_len(nrow(z)), function(i){
  a<-z$foreground_n[i]; K<-z$universe_n[i]; F<-z$foreground_total[i]; N<-z$universe_total[i]
  m<-matrix(c(a, F-a, K-a, N-K-F+a),2); ft2<-fisher.test(m); ft1<-fisher.test(m,alternative='greater')
  data.frame(Family=z$Family[i],genotype=z$genotype[i],OR_r=unname(ft2$estimate),lo_r=ft2$conf.int[1],hi_r=ft2$conf.int[2],p1_r=ft1$p.value,
             OR_src=z$odds_ratio[i],lo_src=z$odds_lower[i],hi_src=z$odds_upper[i],p_src=z$p_greater[i],onesided_CI_upper=ft1$conf.int[2])}))
out <- file.path(dirname(normalizePath(sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))][1]))),'results')
dir.create(out, showWarnings=FALSE); write.csv(res, file.path(out,'fisher_check_results.csv'), row.names=FALSE)
cat(sprintf('max|OR diff| %.2e  max|lower diff| %.2e  max|upper diff| %.2e  max|p diff| %.2e\n', max(abs(res$OR_r-res$OR_src)), max(abs(res$lo_r-res$lo_src)), max(abs(res$hi_r-res$hi_src)), max(abs(res$p1_r-res$p_src))))
cat('Source OR/CI = two-sided conditional-MLE fisher.test; source p = one-sided (greater). One-sided CI upper all Inf:', all(is.infinite(res$onesided_CI_upper)), '\n')
