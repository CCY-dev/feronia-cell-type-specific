#!/usr/bin/env Rscript
# Descriptive as-exported union only. No model fitting or source repairs.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(FALSE)
out <- dirname(normalizePath(sub('^--file=', '', a[grepl('^--file=', a)][1])))
root <- normalizePath(file.path(out, '../../../..'))
inp <- file.path(root, 'revisions/2026-09-20_dof_subnetworks/inventory/tested10_retained_edges.csv')
ann <- file.path(root, 'revisions/2026-09-20_dof_subnetworks/annotations/clean_programme_target_edges.tsv')
e <- fread(inp)
s <- unique(fread(ann)[, .(target, symbol)])
stopifnot(!anyDuplicated(s$target), !anyDuplicated(e[, .(source, target, cluster)]))
e[, target_export_name := target_name]
e[target %in% s$target, target_name := s$symbol[match(target, s$target)]]
e <- e[grepl('^(IAA[0-9]|SAUR[0-9]|GH3[.]|TAA[0-9]|TAR[0-9]|YUC[0-9])', target_name)]
e[, support_category := fifelse(tissue %in% c('phloem', 'xylem'), 'phloem_xylem_target_map_issue',
 fifelse(edge_score_integrity == 'edge_export_score_mismatch', 'endodermis_edge_score_mismatch',
 fifelse(tissue == 'endodermis', 'other_endodermis_conservatively_excluded', 'current_seven_included_tissues')))]
pasteu <- function(x) paste(sort(unique(x)), collapse='; ')
u <- e[, .(n_contexts=uniqueN(cell_type), contexts=pasteu(cell_type),
 contexts_WT=pasteu(cell_type[genotype=='WT']), contexts_fer4=pasteu(cell_type[genotype=='fer4']),
 n_context_genotype_links=.N,
 n_current_seven=sum(support_category=='current_seven_included_tissues'),
 n_phloem_xylem_map_issue=sum(support_category=='phloem_xylem_target_map_issue'),
 n_endodermis_score_mismatch=sum(support_category=='endodermis_edge_score_mismatch'),
 n_other_endodermis=sum(support_category=='other_endodermis_conservatively_excluded'),
 support_categories=pasteu(support_category)), by=.(source,source_symbol,target,target_name)]
u[, support_class := fifelse(n_current_seven>0, 'supported_in_current_seven', 'additional_excluded_contexts_only')]
reference <- fread(file.path(root, 'revisions/2026-09-20_compact_figures/networks/all_tissues_overview_compact_union_edges.csv'))
key <- function(z) paste(z$source,z$target,sep=':')
stopifnot(setequal(key(u[n_current_seven>0]), key(reference)))
fwrite(e,file.path(out,'all_contexts_as_exported_source_context_edges.csv'))
fwrite(u,file.path(out,'all_contexts_as_exported_union_edges.csv'))
fwrite(e[,.(retained_links=.N,distinct_pairs=uniqueN(paste(source,target))),by=.(tissue,cell_type,genotype,support_category)],file.path(out,'context_support_counts.csv'))
fwrite(u[,.(pairs=.N,TFs=uniqueN(source),targets=uniqueN(target)),by=support_class],file.path(out,'union_support_counts.csv'))
summary <- data.table(TFs=uniqueN(e$source),targets=uniqueN(e$target),pairs=nrow(u),
 contexts=uniqueN(e$cell_type),genotype_contexts=uniqueN(e$cluster),tissue_fits=uniqueN(e$tissue),source_rows=nrow(e))
fwrite(summary,file.path(out,'validation.csv'))
writeLines(c(paste('source:',inp),paste('source MD5:',tools::md5sum(inp)),
 paste('symbols:',ann),paste('symbols MD5:',tools::md5sum(ann)),
 'Same auxin-family regex and AGI symbol normalization as the compact overview.',
 'All original tissue_map and edge_score flags retained. Nothing repaired or refitted.'),file.path(out,'input_provenance.txt'))
print(summary)
print(u[,.(pairs=.N,TFs=uniqueN(source),targets=uniqueN(target)),by=support_class])
print(e[,.(links=.N,contexts=uniqueN(cell_type),genotype_contexts=uniqueN(cluster)),by=support_category])
