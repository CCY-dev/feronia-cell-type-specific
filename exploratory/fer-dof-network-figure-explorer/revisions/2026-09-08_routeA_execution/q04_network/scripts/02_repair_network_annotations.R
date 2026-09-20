#!/usr/bin/env Rscript
# Prototype; not yet carefully vetted.
# Q04: repair derived DE annotations; audit target selection and candidate context.
# Sources are immutable. Existing network coefficients/scores are never refitted.
suppressPackageStartupMessages(library(data.table))
setDTthreads(2)
a <- commandArgs(TRUE); ROOT <- a[1]
FER <- '/resnick/groups/NolanLab/cchen/Feronia'
CO <- file.path(FER,'analysis/CellOracle/CellOracle_output/20260317_feronia_salinity_timecourse_WT_untreated_fer4_untreated_newIntegratedSeurat3_NoSuspiciousCluster')
TD <- file.path(CO,'analysis/tables/TF_gene')
tag <- '20260317_feronia_salinity_WT_untreated_fer4_untreated_CellOracle_'
SRC <- '/resnick/groups/NolanLab/shared_resources/source_files/celloracle_inputs'
DEFILE <- file.path(FER,'analysis/muscat_result/20260317_salinity_timecourse_useDateT_noSuspiciousClusters/group_timecelltypeannoLi/fer_4-WT/sig_de.csv')
OLDFILE <- file.path(FER,'analysis/muscat_result/20250228_salinity_timecourse_useDateT/group_timecelltypeannoLi/fer_4-WT/sig_de.csv')
FILFILE <- sub('sig_de.csv','sig_de_fil.csv',DEFILE,fixed=TRUE)
write <- function(d,n) fwrite(d,file.path(ROOT,'results',n),na='NA')
source_paths <- c(DEFILE,OLDFILE,FILFILE,file.path(TD,paste0(tag,c('all_metrics_rank_change_wide_TF_GE.txt','WT_specific_TF.txt','fer4_specific_TF.txt'))),file.path(FER,'Ath_TF_list.txt'),file.path(FER,'claude_script/coauthor_AtDOF_list_AGI_fixed.txt'),file.path(SRC,'Kay_TF_thalemine_annotations.csv'))
wide <- fread(source_paths[4]);stopifnot(!anyDuplicated(wide[,.(gene,cell_type,metric)]))
wide[,row_original:=.I]
DE <- fread(DEFILE,select=c('gene','cluster_id','sc_292.cpm','sc_294.cpm','sc_291.cpm','sc_293.cpm','fer_4.frq','WT.frq','logFC','logCPM','p_val','p_adj.loc','p_adj.glb'))
setnames(DE,'cluster_id','cell_type');stopifnot(!anyDuplicated(DE[,.(gene,cell_type)]))
DE[,`:=`(WT_mean_CPM=(sc_291.cpm+sc_293.cpm)/2,fer4_mean_CPM=(sc_292.cpm+sc_294.cpm)/2,
 eligible_detection=fer_4.frq>=0.1|WT.frq>=0.1)]
DE[,current_called:=eligible_detection&p_adj.loc<=0.05&abs(logFC)>log2(1.5)]
DE[,current_DE_status:=fifelse(!eligible_detection,'tested_below_detection_threshold',fifelse(current_called,fifelse(logFC>0,'called_up','called_down'),'tested_not_called'))]
old <- fread(OLDFILE,select=c('gene','cluster_id','logFC','p_adj.loc'));setnames(old,c('cluster_id','logFC','p_adj.loc'),c('cell_type','old_source_logFC','old_source_p_adj'))
stopifnot(!anyDuplicated(old[,.(gene,cell_type)]))
fixed <- merge(wide,DE,by=c('gene','cell_type'),all.x=TRUE,sort=FALSE)
fixed <- merge(fixed,old,by=c('gene','cell_type'),all.x=TRUE,sort=FALSE);setorder(fixed,row_original)
stopifnot(nrow(fixed)==nrow(wide))
setnames(fixed,c('log2FC','p_adj','DE_label'),c('legacy_log2FC','legacy_p_adj','legacy_DE_label'))
fixed[is.na(current_DE_status),current_DE_status:='not_in_current_DE_output']
fixed[,`:=`(log2FC=logFC,p_adj=p_adj.loc,DE_label=fifelse(is.na(current_called),NA_character_,fifelse(current_called,'DE','None')),
 current_label_same_legacy_rule=fifelse(is.na(logFC),NA_character_,fifelse(p_adj.loc<0.05&abs(logFC)>log2(1.5),'DE','None')))]
unchanged <- setdiff(names(wide),c('log2FC','p_adj','DE_label'))
stopifnot(identical(as.data.frame(wide[,..unchanged]),as.data.frame(fixed[,..unchanged])))
write(fixed,'network_rank_table_DE_repaired.csv')
pairs <- unique(fixed[,.(gene,cell_type,legacy_log2FC,legacy_p_adj,legacy_DE_label,log2FC,p_adj,DE_label,current_label_same_legacy_rule,current_DE_status,old_source_logFC,old_source_p_adj)])
write(pairs,'annotation_repair_gene_contexts.csv')
summary <- data.table(question=c('wide_rows','unique_gene_contexts','comparable_current_contexts','logFC_difference_gt_0.01','legacy_vs_current_same_rule_label_changes','legacy_vs_current_final_rule_label_changes','legacy_matches_old_logFC_1e_10','legacy_comparable_old_contexts'),value=c(nrow(wide),nrow(pairs),sum(!is.na(pairs$log2FC)&!is.na(pairs$legacy_log2FC)),sum(abs(pairs$log2FC-pairs$legacy_log2FC)>0.01,na.rm=TRUE),sum(pairs$legacy_DE_label!=pairs$current_label_same_legacy_rule,na.rm=TRUE),sum(pairs$legacy_DE_label!=pairs$DE_label,na.rm=TRUE),sum(abs(pairs$legacy_log2FC-pairs$old_source_logFC)<1e-10,na.rm=TRUE),sum(!is.na(pairs$legacy_log2FC)&!is.na(pairs$old_source_logFC))))
write(summary,'annotation_repair_summary.csv');print(summary)
fg <- list()
for (label in c('WT_specific','fer4_specific')) {
 d <- fixed[metric=='degree_centrality_out'&Is_TF=='Yes'&Label==label]
 orig <- fread(file.path(TD,paste0(tag,label,'_TF.txt')))
 stopifnot(setequal(d$gene,orig$gene),nrow(d)==nrow(orig))
 fg[[label]] <- unique(d$gene);setorder(d,cell_type,-log2FC)
 write(d,paste0(label,'_TF_DE_repaired.csv'))
}
# Explicit tissue mapping: all observed labels were independently enumerated.
map <- fread(file.path(ROOT,'inputs/explicit_subgroup_tissue_map.csv'))
obs <- fread(file.path(ROOT,'results/source_subgroup_counts.csv'))
stopifnot(setequal(map$cell_type,obs$cell_type),!anyDuplicated(map$cell_type),!any(map$tissue=='UNMAPPED'))
DE <- merge(DE,map,by='cell_type',all.x=TRUE,sort=FALSE);stopifnot(!anyNA(DE$tissue))
fil <- fread(FILFILE,select=c('gene','cluster_id'));setnames(fil,'cluster_id','cell_type')
stopifnot(setequal(paste(DE[current_called==TRUE]$gene,DE[current_called==TRUE]$cell_type),paste(fil$gene,fil$cell_type)))
fil <- merge(fil,map,by='cell_type',all.x=TRUE)
tf <- fread(file.path(SRC,'Kay_TF_thalemine_annotations.csv'));setnames(tf,names(tf)[1:3],c('gene','TF_name','TF_description'));stopifnot(!anyDuplicated(tf$gene))
vars <- fread(file.path(ROOT,'results/h5ad_gene_universes.csv'));allvar<-vars[tissue=='all',gene]
tissues <- c('cortex','endodermis','trichoblast','atrichoblast','pericycle','phloem','xylem','procambium','lrc','columella')
allscore <- list(); targetsummary<-list();membership<-list();eligibility<-list();edgesums<-list();candedge<-list();dof_targets<-list()
fam <- unique(fread(file.path(FER,'Ath_TF_list.txt'))[,.(gene=Gene_ID,Family)])
ambig <- fam[, .N,by=gene][N>1];write(fam[gene%in%ambig$gene],'PlantTFDB_multifamily_genes.csv')
fam_unambiguous<-fam[!gene%in%ambig$gene]
# Reproduce the original alphabetical-first transcript-to-family resolution,
# and separately exclude ambiguous genes in the sensitivity below.
setorder(fam,gene,Family);fam<-fam[,.(Family=Family[1]),by=gene]
coDof <- scan(file.path(FER,'claude_script/coauthor_AtDOF_list_AGI_fixed.txt'),what='character',quiet=TRUE)
candidate <- unique(c(coDof,fam[Family=='Dof',gene],'AT5G43890'))
for (t in tissues) {
 oldgroups <- map[grepl(if(t=='lrc')'lateral root cap' else t,cell_type,ignore.case=TRUE),cell_type]
 if (t=='trichoblast') oldgroups <- oldgroups[!grepl('atrichoblast',oldgroups,ignore.case=TRUE)]
 newgroups <- map[tissue==t,cell_type]
 membership[[t]]<-data.table(tissue=t,cell_type=map$cell_type,legacy_substring=map$cell_type%in%oldgroups,explicit_membership=map$cell_type%in%newgroups)
 oldgenes <- intersect(allvar,union(tf$gene,fil[cell_type%in%oldgroups,gene]))
 newgenes <- intersect(allvar,union(tf$gene,fil[cell_type%in%newgroups,gene]))
 actual <- vars[tissue==t,gene]
 targetsummary[[t]] <- data.table(tissue=t,actual_input_genes=length(actual),legacy_reconstructed=length(oldgenes),corrected_input_genes=length(newgenes),actual_only=length(setdiff(actual,oldgenes)),reconstructed_only=length(setdiff(oldgenes,actual)),legacy_only_genes=length(setdiff(oldgenes,newgenes)),corrected_only_genes=length(setdiff(newgenes,oldgenes)))
 eligibility[[t]] <- data.table(tissue=t,gene=union(oldgenes,newgenes),legacy_eligible=union(oldgenes,newgenes)%in%oldgenes,explicit_eligible=union(oldgenes,newgenes)%in%newgenes,actual_input=union(oldgenes,newgenes)%in%actual)
 scorefile<-file.path(CO,paste0(t,'_rawRNA_celloracle_gene_score_iGRN.csv'));source_paths<-c(source_paths,scorefile)
 s<-fread(scorefile);setnames(s,names(s)[1],'gene');s[,`:=`(tissue=t,cell_type=sub('_(WT|fer4)_.*','',cluster),genotype=sub('.*_(WT|fer4)_untreated$','\\1',cluster))]
 stopifnot(!anyDuplicated(s[,.(gene,cluster)]));allscore[[t]]<-s
 edgefile<-file.path(CO,paste0(t,'_rawRNA_celloracle_network_iGRN.csv'));source_paths<-c(source_paths,edgefile)
 e<-fread(edgefile);e[,tissue:=t]
 e[,target_excluded_by_explicit_map:=target%in%setdiff(oldgenes,newgenes)]
 edgesums[[t]]<-e[,.(edges=.N,edges_to_excluded_targets=sum(target_excluded_by_explicit_map),excluded_target_genes=uniqueN(target[target_excluded_by_explicit_map]),affected_source_TFs=uniqueN(source[target_excluded_by_explicit_map])),by=.(tissue,cluster)]
 candedge[[t]]<-e[target=='AT5G43890'|source%in%candidate][,.(tissue,source,source_name,target,target_name,cluster,coef_mean,coef_abs,p,target_excluded_by_explicit_map)]
 dof_targets[[t]]<-e[source%in%candidate,.(retained_targets=.N,excluded_targets=sum(target_excluded_by_explicit_map)),by=.(tissue,source,cluster)]
 rm(e);gc()
}
scores <- rbindlist(allscore);rm(allscore)
write(rbindlist(targetsummary),'tissue_target_universe_audit.csv')
write(rbindlist(membership),'explicit_vs_substring_membership.csv')
elig<-rbindlist(eligibility);write(elig,'all_target_gene_eligibility.csv')
write(elig[gene%in%candidate],'DOF_YUC5_target_eligibility.csv')
write(rbindlist(edgesums),'network_edges_affected_by_target_mapping.csv')
candedge<-rbindlist(candedge);write(candedge,'DOF_edges_and_all_YUC5_edges.csv')
write(candedge[target=='AT5G43890'],'all_retained_YUC5_edges.csv')
write(rbindlist(dof_targets),'DOF_outdegree_edge_deletion_sensitivity.csv')
write(scores[gene%in%candidate],'DOF_YUC5_all_network_scores.csv')
write(merge(DE[gene%in%candidate],tf,by='gene',all.x=TRUE),'DOF_YUC5_current_DE_expression.csv')
# Every TF that has scores in both genotypes in at least one cell-stage context is
# eligible for a paired rank comparison before the top100 filter.
pair_scores<-dcast(scores, gene+tissue+cell_type~genotype,value.var='degree_centrality_out')
stopifnot(all(c('WT','fer4')%in%names(pair_scores)))
pair_genes<-unique(pair_scores[is.finite(WT)&is.finite(fer4),gene])
backgrounds<-list(legacy_selected_wide=unique(wide$gene),all_score_PlantTFDB=unique(scores$gene),paired_score_KayTF=intersect(pair_genes,tf$gene),legacy_selected_KayTF=intersect(unique(wide$gene),tf$gene))
famco<-copy(fam);famco[Family=='Dof'&!gene%in%coDof,Family:='Dof_PlantTFDB_only'];famco[gene%in%coDof,Family:='Dof'];famco<-rbind(famco,data.table(gene=setdiff(coDof,famco$gene),Family='Dof'))
write(merge(data.table(gene=union(fam[Family=='Dof',gene],coDof)),fam,by='gene',all.x=TRUE)[,coauthor_DOF:=gene%in%coDof],'DOF_family_definition_comparison.csv')
enrichment <- list(); universegenes<-list()
for (definition in c('PlantTFDB','coauthor_DOF_override','PlantTFDB_unambiguous')) {
 fm<-switch(definition,PlantTFDB=fam,coauthor_DOF_override=famco,PlantTFDB_unambiguous=fam_unambiguous)
 for (universe in names(backgrounds)) {
  u<-fm[gene%in%backgrounds[[universe]]];N<-nrow(u)
  universegenes[[paste(definition,universe)]]<-cbind(u,data.table(definition=definition,universe=universe))
  for(label in names(fg)) {
   foreground<-intersect(fg[[label]],u$gene);n<-length(foreground)
   ans<-rbindlist(lapply(sort(unique(u$Family)),function(f) {
    K<-sum(u$Family==f);k<-sum(u$Family==f&u$gene%in%foreground)
    mat<-matrix(c(k,n-k,K-k,N-n-K+k),nrow=2,byrow=TRUE)
    ft<-fisher.test(mat,alternative='greater');two<-fisher.test(mat,alternative='two.sided')
    data.table(Family=f,universe_n=K,foreground_n=k,universe_total=N,foreground_total=n,foreground_unclassifiable_or_ineligible=length(fg[[label]])-n,fold_enrich=(k/n)/(K/N),odds_ratio=unname(ft$estimate),odds_lower=two$conf.int[1],odds_upper=two$conf.int[2],p_greater=ft$p.value)
   }))
   ans[,`:=`(p_adj_BH=p.adjust(p_greater,'BH'),definition=definition,universe=universe,genotype=label)]
   enrichment[[paste(definition,universe,label)]]<-ans
  }
 }
}
enrich<-rbindlist(enrichment);write(enrich,'TF_family_enrichment_universe_sensitivity.csv');write(rbindlist(universegenes),'TF_family_background_membership.csv');print(enrich[Family=='Dof'])
# Reproduce top100 rank selection from raw score files, then vary top-N/quantiles.
metrics<-c('degree_centrality_all','degree_centrality_in','degree_centrality_out','betweenness_centrality')
long<-melt(scores,id.vars=c('gene','tissue','cell_type','cluster','genotype'),measure.vars=metrics,variable.name='metric',value.name='value')
long[,gene_rank:=rank(-value,ties.method='min'),by=.(cluster,metric)]
pairrank<-dcast(long,gene+tissue+cell_type+metric~genotype,value.var='gene_rank')
pairrank<-pairrank[is.finite(WT)&is.finite(fer4)];pairrank[,rank_diff:=fer4-WT]
ranksets<-list();thresholds<-list();sensitivity<-list()
for (topN in c(50,100,200)) for (tail in c(.05,.1,.2)) {
 d<-copy(pairrank[pmin(WT,fer4)<=topN])
 thr<-d[,.(low=as.numeric(quantile(rank_diff,tail)),high=as.numeric(quantile(rank_diff,1-tail))),by=metric]
 d<-merge(d,thr,by='metric');d[,Label:=fifelse(rank_diff<low,'fer4_specific',fifelse(rank_diff>high,'WT_specific','None'))]
 if (topN==100&&tail==.1) {
  check<-merge(wide[,.(gene,tissue=crude_cell_type,cell_type,metric,legacy_Label=Label,legacy_diff=get('fer4_untreated-WT_untreated: rank_diff'))],d,by=c('gene','tissue','cell_type','metric'),all=TRUE)
  write(check,'top100_exact_reproduction_comparison.csv')
  stopifnot(nrow(check)==nrow(wide),!anyNA(check$legacy_Label),!anyNA(check$Label),all(check$legacy_Label==check$Label),all(check$legacy_diff==check$rank_diff))
 }
 selected<-d[metric=='degree_centrality_out'&gene%in%tf$gene&Label!='None'];selected[,`:=`(topN=topN,tail_fraction=tail)]
 ranksets[[paste(topN,tail)]]<-selected
 thresholds[[paste(topN,tail)]]<-cbind(thr,data.table(topN=topN,tail_fraction=tail))
 for (lab in names(fg)) {
  g<-unique(selected[Label==lab,gene]);basefg<-fg[[lab]]
  sensitivity[[paste(topN,tail,lab)]]<-data.table(topN=topN,tail_fraction=tail,genotype=lab,n_TFs=length(g),original_recovered=length(intersect(g,basefg)),original_n=length(basefg),Jaccard=length(intersect(g,basefg))/length(union(g,basefg)),n_Dof_PlantTFDB=sum(g%in%fam[Family=='Dof',gene]))
 }
}
write(rbindlist(ranksets),'TF_nomination_threshold_sensitivity.csv');write(rbindlist(thresholds),'rank_selection_thresholds.csv');write(rbindlist(sensitivity),'TF_nomination_sensitivity_summary.csv')
# Baseline degree/expression-adjusted DOF association, descriptive sensitivity.
# Gene is the unit; aggregate across contexts, do not treat contexts as replicates.
base_degree<-scores[genotype=='WT',.(WT_median_outdegree=median(degree_out),scored_contexts=.N),by=gene]
base_expr<-DE[,.(WT_median_CPM=median(WT_mean_CPM),WT_median_frequency=median(WT.frq)),by=gene]
u<-merge(merge(fam[gene%in%backgrounds$paired_score_KayTF],base_degree,by='gene'),base_expr,by='gene',all.x=TRUE)
write(u,'paired_TF_baseline_covariates.csv')
adjusted<-list()
for (lab in names(fg)) {
 dd<-copy(u);dd[,`:=`(selected=as.integer(gene%in%fg[[lab]]),Dof=as.integer(Family=='Dof'))]
 dd<-dd[complete.cases(dd)]
 fit<-glm(selected~Dof+log1p(WT_median_outdegree)+log1p(WT_median_CPM)+log1p(scored_contexts),data=dd,family=binomial())
 cf<-coef(summary(fit))['Dof',]
 adjusted[[lab]]<-data.table(genotype=lab,n_genes=nrow(dd),n_Dof=sum(dd$Dof),odds_ratio=exp(cf[1]),Wald_lower=exp(cf[1]-1.96*cf[2]),Wald_upper=exp(cf[1]+1.96*cf[2]),p_Wald=cf[4],interpretation='Exploratory conditional association; prior degree unavailable; no causal or confirmatory test')
 write(as.data.table(coef(summary(fit)),keep.rownames='term'),paste0('exploratory_covariate_model_',lab,'.csv'))
}
write(rbindlist(adjusted),'DOF_covariate_adjusted_sensitivity.csv')
# Save reusable compact objects for figures, not the large source matrices.
saveRDS(list(enrich=enrich,repair=summary,target=rbindlist(targetsummary),edge=rbindlist(edgesums),sensitivity=rbindlist(sensitivity)),file.path(ROOT,'results/q04_figure_data.rds'))
write(data.table(path=source_paths,size_bytes=file.info(source_paths)$size,mtime=as.character(file.info(source_paths)$mtime)), 'source_file_manifest.csv')
writeLines(capture.output(sessionInfo()),file.path(ROOT,'logs/R_sessionInfo.txt'))
cat('All q04 numerical checks passed.\n')
