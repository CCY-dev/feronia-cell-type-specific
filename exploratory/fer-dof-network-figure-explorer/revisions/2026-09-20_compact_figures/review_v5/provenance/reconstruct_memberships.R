#!/usr/bin/env Rscript
# Prototype; not yet carefully vetted.
# Reproduce memberships from the exact sources documented in network_audit.md.
# This is a provenance check. No score repair, network fitting or source modification.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(FALSE)
script <- normalizePath(sub('^--file=','',args[grepl('^--file=',args)][1]))
out <- dirname(script)
root <- normalizePath(file.path(out,'../../../..'))
q04 <- file.path(root,'revisions/2026-09-08_routeA_execution/q04_network/results')
pub <- file.path(root,'revisions/2026-09-08_fresh_audit/work/feronia-cell-type-specific')
old <- file.path(root,'revisions/2026-09-20_fig6_brainstorm')
input <- data.table(
 role=c('calculation: WT nominations','calculation: fer4 nominations','calculation: coauthor DOF family',
        'calculation: called DEG set','annotation: ten assayed DOFs and YUC5 contexts',
        'annotation: nine-threshold nomination support','annotation: display symbols',
        'comparison only: original DOF membership','comparison only: original intersections'),
 path=c(file.path(q04,'WT_specific_TF_DE_repaired.csv'),
        file.path(q04,'fer4_specific_TF_DE_repaired.csv'),
        file.path(q04,'DOF_family_definition_comparison.csv'),
        file.path(pub,'scripts_and_inputs/Analysis/Differential_expression_analysis/outputs/sig_de_fil.csv'),
        file.path(q04,'tested_DOF_integrated_evidence_matrix.csv'),
        file.path(q04,'TF_nomination_threshold_sensitivity.csv'),
        file.path(q04,'DOF_YUC5_current_DE_expression.csv'),
        file.path(old,'fig6_dof_membership.csv'),file.path(old,'fig6_full_intersections.csv')))
stopifnot(all(file.exists(input$path)))
sha256 <- function(p) {
 x <- system2('shasum',c('-a','256',shQuote(p)),stdout=TRUE)
 sub(' .*','',x[1])
}
input[,`:=`(bytes=file.info(path)$size,sha256=vapply(path,sha256,character(1)))]
fwrite(input,file.path(out,'input_sha256.csv'))
write_csv <- function(x,name)fwrite(x,file.path(out,name),na='NA')
as_bool <- function(x) {
 z<-tolower(as.character(x))
 stopifnot(all(is.na(z)|z%in%c('true','false','1','0')))
 ifelse(is.na(z),NA,z%in%c('true','1'))
}

# Derive set membership independently, before reading the reference output tables.
wt_raw <- fread(input$path[1]);fer_raw <- fread(input$path[2])
family <- fread(input$path[3]);de <- fread(input$path[4])
stopifnot(!anyDuplicated(family$gene),all(c('gene','metric','Is_TF','Label')%in%names(wt_raw)),
          all(c('gene','metric','Is_TF','Label')%in%names(fer_raw)))
wt <- wt_raw[metric=='degree_centrality_out' & Is_TF=='Yes' & Label=='WT_specific']
fer <- fer_raw[metric=='degree_centrality_out' & Is_TF=='Yes' & Label=='fer4_specific']
selection <- data.table(input=c('WT','fer4'),rows_in_file=c(nrow(wt_raw),nrow(fer_raw)),
 rows_after_outdegree_TF_label_filter=c(nrow(wt),nrow(fer)),
 other_metric_rows=c(sum(wt_raw$metric!='degree_centrality_out'),sum(fer_raw$metric!='degree_centrality_out')))
write_csv(selection,'nomination_input_selection.csv')
sets <- list(DOF=sort(unique(family$gene[as_bool(family$coauthor_DOF)])),
             GRN_WT=sort(unique(wt$gene)),GRN_fer4=sort(unique(fer$gene)),
             called_DEG=sort(unique(de$gene)))
sets$GRN_nominated <- sort(union(sets$GRN_WT,sets$GRN_fer4))
stopifnot(all(vapply(sets,function(g)!anyNA(g)&&all(grepl('^AT[1-5CM]G[0-9]{5}$',g)),logical(1))))
for(n in names(sets))write_csv(data.table(gene=sets[[n]]),paste0('set_',n,'.csv'))
all_genes <- sort(Reduce(union,sets))
u <- data.table(gene=all_genes)
u[,`:=`(DOF=as.integer(gene%in%sets$DOF),GRN_nominated=as.integer(gene%in%sets$GRN_nominated),
         called_DEG=as.integer(gene%in%sets$called_DEG),
         nominated_WT=gene%in%sets$GRN_WT,nominated_fer4=gene%in%sets$GRN_fer4)]
write_csv(u,'reconstructed_full_union_membership.csv')
patterns <- CJ(DOF=0:1,GRN_nominated=0:1,called_DEG=0:1)
patterns <- patterns[DOF+GRN_nominated+called_DEG>0]
counts <- u[,.(genes=.N),by=.(DOF,GRN_nominated,called_DEG)]
counts <- merge(patterns,counts,by=c('DOF','GRN_nominated','called_DEG'),all.x=TRUE)
counts[is.na(genes),genes:=0L]
setorder(counts,DOF,GRN_nominated,called_DEG)
write_csv(counts,'reconstructed_fig6_full_intersections.csv')

# Optional metadata reproduces the remaining original display columns.
# It is not used to decide DOF, GRN or DEG set membership.
tested <- fread(input$path[5]);sens <- fread(input$path[6]);expr <- fread(input$path[7])
stopifnot(!anyDuplicated(tested$gene))
sens <- sens[metric=='degree_centrality_out' & Label%in%c('WT_specific','fer4_specific')]
settings <- unique(sens[,.(gene,topN,tail_fraction)])[,.(n_nomination_settings_of_9=.N),by=gene]
symbol_map <- unique(expr[!is.na(TF_name)&nzchar(TF_name),.(gene,symbol=TF_name)])
ambiguous <- symbol_map[,.(symbols=paste(sort(unique(symbol)),collapse=';'),n_symbols=uniqueN(symbol)),by=gene][n_symbols>1]
write_csv(ambiguous,'ambiguous_source_symbols.csv')
if(nrow(ambiguous))stop('Source symbol ambiguity: memberships are saved; display metadata needs review.')
m <- data.table(gene=sets$DOF)
m[,symbol:=symbol_map$symbol[match(gene,symbol_map$gene)]]
m[is.na(symbol)|!nzchar(symbol),symbol:=gene]
m[gene%in%tested$gene,symbol:=tested$symbol[match(gene,tested$gene)]]
m[,`:=`(nominated_any=gene%in%sets$GRN_nominated,nominated_WT=gene%in%sets$GRN_WT,
         nominated_fer4=gene%in%sets$GRN_fer4,called_DE_any=gene%in%sets$called_DEG,
         tested_LUC=gene%in%tested$gene,
         n_nomination_settings_of_9=settings$n_nomination_settings_of_9[match(gene,settings$gene)],
         retained_YUC5_contexts=tested$contexts_retained_YUC5_edge[match(gene,tested$gene)],
         retained_YUC5_tissues=tested$tissues_retained_YUC5_edge[match(gene,tested$gene)])]
m[is.na(n_nomination_settings_of_9),n_nomination_settings_of_9:=0L]
setcolorder(m,c('gene','symbol','nominated_any','nominated_WT','nominated_fer4','called_DE_any','tested_LUC',
 'n_nomination_settings_of_9','retained_YUC5_contexts','retained_YUC5_tissues'))
write_csv(m,'reconstructed_fig6_dof_membership.csv')

# Only now consult the historical products for comparisons, not construction.
old_m <- fread(input$path[8]);old_c <- fread(input$path[9])
stopifnot(!anyDuplicated(old_m$gene),!anyDuplicated(old_c[,.(DOF,GRN_nominated,called_DEG)]))
flags <- c('nominated_any','nominated_WT','nominated_fer4','called_DE_any','tested_LUC')
for(n in flags)old_m[,(n):=as_bool(get(n))]
cmp <- merge(m,old_m,by='gene',all=TRUE,suffixes=c('_reconstructed','_reference'))
cols <- setdiff(names(m),'gene')
comparisons <- rbindlist(lapply(cols,function(n) {
 a<-as.character(cmp[[paste0(n,'_reconstructed')]])
 b<-as.character(cmp[[paste0(n,'_reference')]])
 # Blank/NA optional annotation fields have the same missing-data meaning.
 if(n=='retained_YUC5_tissues'){a[is.na(a)]<-'';b[is.na(b)]<-''}
 data.table(gene=cmp$gene,field=n,reconstructed=a,reference=b,
            matches=(is.na(a)&is.na(b))|(!is.na(a)&!is.na(b)&a==b))
}))
write_csv(comparisons,'membership_field_comparison.csv')
write_csv(comparisons[matches==FALSE|is.na(matches)],'membership_differences.csv')
ccomp <- merge(counts,old_c,by=c('DOF','GRN_nominated','called_DEG'),all=TRUE,
               suffixes=c('_reconstructed','_reference'))
ccomp[,matches:=!is.na(genes_reconstructed)&!is.na(genes_reference)&genes_reconstructed==genes_reference]
write_csv(ccomp,'intersection_comparison.csv')
write_csv(ccomp[matches==FALSE|is.na(matches)],'intersection_differences.csv')
same_genes <- setequal(m$gene,old_m$gene)
core_match <- same_genes && all(comparisons[field%in%flags,matches])
intersection_match <- all(ccomp$matches)
metadata_match <- all(comparisons[!field%in%flags,matches])
validation <- data.table(check=c('DOF_unique','WT_nominated_unique','fer4_nominated_unique','GRN_union_unique',
 'called_DEG_unique','three_set_union','nominated_DOF_unique','DEG_DOF_unique','tested_DOF_unique',
 'tested_nominated_and_DEG','tested_nominated_without_DEG','same_DOF_identifiers',
 'core_membership_flags_match','seven_intersections_match','all_optional_metadata_match',
 'membership_field_differences','intersection_pattern_differences'),
 value=as.character(c(length(sets$DOF),length(sets$GRN_WT),length(sets$GRN_fer4),length(sets$GRN_nominated),
 length(sets$called_DEG),nrow(u),sum(m$nominated_any),sum(m$called_DE_any),sum(m$tested_LUC),
 sum(m$tested_LUC&m$nominated_any&m$called_DE_any),sum(m$tested_LUC&m$nominated_any&!m$called_DE_any),
 same_genes,core_match,intersection_match,metadata_match,sum(!comparisons$matches),sum(!ccomp$matches))))
write_csv(validation,'validation.csv')
commit <- system2('git',c('-C',shQuote(pub),'rev-parse','HEAD'),stdout=TRUE)
writeLines(c('Prototype; not yet carefully vetted.',paste('Pinned public clone:',commit),
 'Existing score-derived nomination sources retained; no score reconciliation or refit.',
 'Membership is reconstructed before historical output tables are read for comparison.',
 'Unique genes, not gene-context rows, define all sets. Nominations use degree_centrality_out only.',
 'DOFs use coauthor_DOF TRUE. DEGs use unique gene IDs in pinned sig_de_fil.csv.',
 'No new DE threshold is applied to the already-filtered DE source.',
 paste('Core membership comparison:',if(core_match)'PASS' else 'FAIL'),
 paste('Seven exclusive intersections:',if(intersection_match)'PASS' else 'FAIL'),
 paste('Optional metadata comparison:',if(metadata_match)'PASS' else 'DIFFERENCES RECORDED')),
 file.path(out,'RECONSTRUCTION_STATUS.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
print(validation)
if(!core_match || !intersection_match)stop('Reconstruction conflict. Original files unchanged; differences saved for review.')
if(!metadata_match)warning('Core memberships agree, but optional metadata differences are preserved.')
