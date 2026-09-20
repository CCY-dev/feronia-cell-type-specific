#!/usr/bin/env Rscript
# Prototype; not yet carefully vetted.
# Inventory existing filtered GRN edges. No network refitting or new hypothesis tests.
suppressPackageStartupMessages(library(data.table))
setDTthreads(2)
args <- commandArgs(TRUE)
root <- if(length(args)) normalizePath(args[1]) else normalizePath('.')
src <- file.path(root, 'revisions/2026-09-08_routeA_execution/q04_network/results')
out <- file.path(root, 'revisions/2026-09-20_dof_subnetworks/inventory')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
write_csv <- function(d, name) fwrite(d, file.path(out, name), na='NA')
inputs <- c('DOF_edges_and_all_YUC5_edges.csv', 'tested_DOF_integrated_evidence_matrix.csv',
            'all_target_gene_eligibility.csv', 'DOF_YUC5_all_network_scores.csv',
            'DOF_YUC5_relationship_context_evidence.csv')
e0 <- fread(file.path(src, inputs[1]))
d <- fread(file.path(src, inputs[2]))
stopifnot(nrow(d)==10, !anyDuplicated(d$gene))
write_csv(d[,.(gene, symbol, alias)], 'tested10_identities.csv')
e <- e0[source %in% d$gene]
e[, source_symbol := d$symbol[match(source, d$gene)]]
e[, `:=`(cell_type=sub('_(WT|fer4)_untreated$', '', cluster),
         genotype=sub('.*_(WT|fer4)_untreated$', '\\1', cluster),
         tissue_mapping_status=ifelse(tissue %in% c('phloem','xylem'),
           'affected_fit_target_mapping', 'unaffected_fit'))]
bad_contexts <- c('Maturation_Endodermis_fer4_untreated',
                  'Proliferation Domain_Endodermis_fer4_untreated',
                  'Transition Domain_Endodermis_fer4_untreated')
e[, edge_score_integrity:=ifelse(cluster %in% bad_contexts,
                                'edge_export_score_mismatch', 'verified_outdegree_match')]
stopifnot(all(e$genotype %in% c('WT','fer4')),
          !anyDuplicated(e[,.(tissue, source, target, cluster)]),
          all(is.finite(e$coef_mean)), all(is.finite(e$coef_abs)),
          max(abs(abs(e$coef_mean)-e$coef_abs)) < 1e-7,
          all(e$p < .01))
e[, coefficient_sign := ifelse(coef_mean > 0, 'positive', 'negative')]
elig <- fread(file.path(src, inputs[3]))
stopifnot(!anyDuplicated(elig[,.(tissue,gene)]))
ek <- paste(elig$tissue,elig$gene)
e[, `:=`(source_in_actual_input=elig$actual_input[match(paste(tissue,source),ek)],
         target_in_actual_input=elig$actual_input[match(paste(tissue,target),ek)],
         target_in_explicit_input=elig$explicit_eligible[match(paste(tissue,target),ek)])]
stopifnot(all(e$source_in_actual_input), all(e$target_in_actual_input),
          all(e$target_excluded_by_explicit_map == !e$target_in_explicit_input))
# These symbol-defined groups are annotation leads, not curated pathways or enrichment results.
patterns <- c(auxin_biosynthesis='^(YUC[0-9]+|TAA[0-9]+|TAR[0-9]+)$',
              auxin_transport='^(PIN[0-9]+|AUX1|LAX[0-9]+)$',
              auxin_response='^(ARF[0-9]+.*|IAA[0-9]+|GH3.*|SAUR[0-9]+)$',
              cell_wall='^(EXPA[0-9]+|EXPB[0-9]+|XTH[0-9]+|CESA[0-9]+|CSL[A-Z][0-9]+|PME[0-9]+|PMEI[0-9]+)$',
              RALF_FER_LRX='^(RALF.*|FER|LRX[0-9]+)$',
              RBOH='^RBOH[A-Z]$')
tm <- unique(e[,.(target,target_name)])
stopifnot(!anyDuplicated(tm$target))
leads <- rbindlist(lapply(names(patterns), function(p) {
  z <- tm[grepl(patterns[[p]],target_name)]
  z[, programme_lead:=p]; z
}))
write_csv(data.table(programme_lead=names(patterns),symbol_regex=patterns), 'programme_lead_definitions.csv')
write_csv(leads, 'programme_lead_gene_dictionary.csv')
pl <- leads[,.(programme_leads=paste(programme_lead,collapse=';')),by=target]
e <- merge(e,pl,by='target',all.x=TRUE,sort=FALSE)
e[is.na(programme_leads),programme_leads:='unclassified_by_symbol_screen']
setcolorder(e,c('source','source_symbol','source_name','target','target_name','tissue','cell_type','genotype','cluster',
 'coef_mean','coef_abs','coefficient_sign','p','programme_leads','tissue_mapping_status',
 'source_in_actual_input','target_in_actual_input','target_in_explicit_input','target_excluded_by_explicit_map',
 'edge_score_integrity'))
setorder(e,source_symbol,tissue,cell_type,genotype,-coef_abs,target)
write_csv(e,'tested10_retained_edges.csv')
write_csv(e[tissue_mapping_status=='unaffected_fit'],'tested10_retained_edges_eight_unaffected_tissues.csv')
write_csv(e[tissue_mapping_status=='unaffected_fit' & tissue!='endodermis'],
 'tested10_retained_edges_seven_unaffected_consistent_tissues.csv')
write_csv(e[programme_leads!='unclassified_by_symbol_screen'],'programme_lead_edges.csv')
write_csv(e[programme_leads=='auxin_biosynthesis'],'auxin_biosynthesis_exact_edges.csv')
summ <- function(z) z[,.(retained_edges=.N, distinct_targets=uniqueN(target),
  distinct_source_target_pairs=uniqueN(paste(source,target)), contexts=uniqueN(cluster),
  positive_coefficients=sum(coef_mean>0),negative_coefficients=sum(coef_mean<0),
  mapping_excluded_target_edges=sum(target_excluded_by_explicit_map)),by=.(source,source_symbol)]
sall <- summ(e); sall[,scope:='all_10_tissue_fits']
sclean <- summ(e[tissue_mapping_status=='unaffected_fit']); sclean[,scope:='eight_unaffected_tissue_fits']
write_csv(rbind(sall,sclean),'tested10_summary.csv')
ctx <- e[,.(retained_targets=.N, target_genes=paste(target,collapse=';'),
            target_symbols=paste(target_name,collapse=';'),
            positive_coefficients=sum(coef_mean>0),negative_coefficients=sum(coef_mean<0)),
         by=.(source,source_symbol,tissue,cell_type,genotype,cluster,tissue_mapping_status)]
scores <- fread(file.path(src,inputs[4]))[gene %in% d$gene]
scorejoin <- scores[,.(source=gene,tissue,cell_type,genotype,cluster,score_degree_out=degree_out,
                       score_degree_centrality_out=degree_centrality_out)]
ctx <- merge(scorejoin,ctx,by=c('source','tissue','cell_type','genotype','cluster'),all=TRUE,sort=FALSE)
ctx[is.na(source_symbol),source_symbol:=d$symbol[match(source,d$gene)]]
ctx[is.na(tissue_mapping_status),tissue_mapping_status:=ifelse(tissue %in% c('phloem','xylem'),
   'affected_fit_target_mapping','unaffected_fit')]
ctx[is.na(retained_targets),retained_targets:=0L]
stopifnot(!anyNA(ctx$score_degree_out))
ctx[,count_matches_score:=score_degree_out==retained_targets]
ctx[,edge_score_integrity:=ifelse(cluster %in% bad_contexts,
 'edge_export_score_mismatch','verified_outdegree_match')]
stopifnot(all(ctx[!cluster %in% bad_contexts,count_matches_score]))
write_csv(ctx[count_matches_score==FALSE], 'tested10_edge_score_count_mismatches.csv')
setorder(ctx,source_symbol,tissue,cell_type,genotype)
write_csv(ctx,'tested10_context_target_counts_and_lists.csv')
pair <- dcast(ctx,source+source_symbol+tissue+cell_type+tissue_mapping_status~genotype,
              value.var='retained_targets')
setnames(pair,c('WT','fer4'),c('WT_retained_targets','fer4_retained_targets'))
pair[,both_genotypes_scored:=!is.na(WT_retained_targets)&!is.na(fer4_retained_targets)]
pair[,fer4_minus_WT_retained_targets:=fer4_retained_targets-WT_retained_targets]
pair[,comparison_export_integrity:=ifelse(tissue=='endodermis' &
 cell_type %in% sub('_fer4_untreated$','',bad_contexts),'unresolved_edge_export_score_mismatch',
 'verified_outdegree_match')]
write_csv(pair,'tested10_paired_context_target_counts.csv')
ps <- function(z,scope) {
  a <- z[programme_leads!='unclassified_by_symbol_screen',
    .(retained_edges=.N, distinct_targets=uniqueN(target), distinct_DOFS=uniqueN(source),
      source_symbols=paste(sort(unique(source_symbol)),collapse=';'),
      target_symbols=paste(sort(unique(target_name)),collapse=';'),
      tissues=paste(sort(unique(tissue)),collapse=';')),by=programme_leads]
  a[,scope:=scope]; a
}
write_csv(rbind(ps(e,'all_10_tissue_fits'),ps(e[tissue_mapping_status=='unaffected_fit'],
 'eight_unaffected_tissue_fits')),'programme_lead_summary.csv')
tg <- e[programme_leads!='unclassified_by_symbol_screen',
 .(retained_edges=.N, distinct_DOFS=uniqueN(source),
   source_symbols=paste(sort(unique(source_symbol)),collapse=';'),
   unaffected_fit_edges=sum(tissue_mapping_status=='unaffected_fit'),
   tissues=paste(sort(unique(tissue)),collapse=';'),
   contexts=paste(sort(unique(cell_type)),collapse=';')),
 by=.(programme_leads,target,target_name)]
setorder(tg,programme_leads,target_name)
write_csv(tg,'programme_target_summary.csv')
pcs <- e[programme_leads!='unclassified_by_symbol_screen',
 .(retained_edges=.N, distinct_targets=uniqueN(target), distinct_DOFs=uniqueN(source),
   DOFs=paste(sort(unique(source_symbol)),collapse=';'),
   targets=paste(sort(unique(target_name)),collapse=';')),
 by=.(programme_leads,tissue,cell_type,genotype,tissue_mapping_status)]
write_csv(pcs[order(programme_leads,tissue,cell_type,genotype)],'programme_context_summary.csv')
yg <- unique(e[programme_leads=='auxin_biosynthesis',.(target,target_name)])
ye <- CJ(tissue=sort(unique(e$tissue)),target=yg$target)
ye[,target_name:=yg$target_name[match(target,yg$target)]]
ye[,`:=`(actual_input=elig$actual_input[match(paste(tissue,target),ek)],
         explicit_eligible=elig$explicit_eligible[match(paste(tissue,target),ek)])]
# Genes absent from the union of legacy/corrected input sets are in neither input.
ye[is.na(actual_input),`:=`(actual_input=FALSE,explicit_eligible=FALSE)]
write_csv(ye,'observed_auxin_biosynthesis_target_input_eligibility.csv')
audit <- data.table(check=c('tested_DOFs','retained_edges','distinct_target_genes','distinct_source_target_pairs',
 'duplicate_source_target_cluster_keys','affected_tissue_edges','edges_to_mapping_excluded_targets',
 'coefficient_abs_max_difference','maximum_retained_p','candidate_context_count_mismatches'),
 value=c(nrow(d),nrow(e),uniqueN(e$target),uniqueN(paste(e$source,e$target)),
 anyDuplicated(e[,.(tissue,source,target,cluster)]),sum(e$tissue %in% c('phloem','xylem')),
 sum(e$target_excluded_by_explicit_map),max(abs(abs(e$coef_mean)-e$coef_abs)),max(e$p),sum(!ctx$count_matches_score)))
write_csv(audit,'validation.csv')
# Independent cross-check against the full public edge and score exports, all TFs.
pub <- file.path(root,'revisions/2026-09-08_fresh_audit/work/feronia-cell-type-specific',
 'scripts_and_inputs/Analysis/CellOracle_GRN_analysis/outputs')
fa <- lapply(sort(unique(e$tissue)),function(t) {
 ef <- file.path(pub,paste0(t,'_rawRNA_celloracle_network_iGRN.csv'))
 sf <- file.path(pub,paste0(t,'_rawRNA_celloracle_gene_score_iGRN.csv'))
 if(!file.exists(ef)||!file.exists(sf)) return(NULL)
 ee <- fread(ef,select=c('source','target','cluster'))
 ss <- fread(sf);setnames(ss,names(ss)[1],'source')
 z <- merge(ss[,.(source,cluster,degree_out)],ee[,.(edge_count=.N),by=.(source,cluster)],
  by=c('source','cluster'),all=TRUE)
 z[is.na(edge_count),edge_count:=0L]
 zz <- z[,.(score_genes=sum(!is.na(degree_out)),edge_source_genes=sum(edge_count>0),
  sources_absent_from_scores=sum(is.na(degree_out)),
  scored_sources_with_count_mismatch=sum(degree_out!=edge_count,na.rm=TRUE),
  total_retained_edges=sum(edge_count),total_score_outdegree=sum(degree_out,na.rm=TRUE)),by=cluster]
 zz[,tissue:=t];zz
})
full_audit <- rbindlist(fa)
write_csv(full_audit,'full_public_export_degree_audit.csv')
de_file <- file.path(out,'matching_basal_DE.csv')
if(file.exists(de_file)) {
 de <- fread(de_file)
 stopifnot(!anyDuplicated(de[,.(gene,cluster_id)]))
 de[,detection_eligible:=fer_4.frq>=0.1|WT.frq>=0.1]
 de[,DE_called:=detection_eligible & p_adj.loc<=.05 & abs(logFC)>log2(1.5)]
 de[,DE_status:=fifelse(!detection_eligible,'below_detection_threshold',
  fifelse(DE_called,fifelse(logFC>0,'called_up','called_down'),'tested_not_called'))]
 write_csv(de,'matching_basal_DE_annotated.csv')
 write_csv(data.table(source='/resnick/groups/NolanLab/cchen/Feronia/analysis/muscat_result/20260317_salinity_timecourse_useDateT_noSuspiciousClusters/group_timecelltypeannoLi/fer_4-WT/sig_de.csv',
  extracted_rows=nrow(de),extracted_gene_ids=uniqueN(de$gene),cell_contexts=uniqueN(de$cluster_id),
  method='Read-only remote Python CSV stream; subset to all ten DOFs plus all their retained fitted target IDs'),
  'basal_DE_extraction_provenance.csv')
}
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
write_csv(data.table(file=inputs,md5=unname(tools::md5sum(file.path(src,inputs)))), 'input_md5.csv')
print(audit); print(sall); print(e[programme_leads=='auxin_biosynthesis',
 .(source_symbol,target_name,cell_type,genotype,coef_mean,tissue_mapping_status)])
