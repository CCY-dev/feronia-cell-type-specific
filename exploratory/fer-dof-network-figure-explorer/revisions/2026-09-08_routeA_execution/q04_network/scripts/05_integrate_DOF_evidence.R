#!/usr/bin/env Rscript
# Prototype; not yet carefully vetted.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(cowplot)})
setDTthreads(2);ROOT<-commandArgs(TRUE)[1]
read<-function(f)fread(file.path(ROOT,'results',f))
write<-function(d,f)fwrite(d,file.path(ROOT,'results',f),na='NA')
prior<-read('supplied_prior_TF_target_counts.csv');yucprior<-read('supplied_prior_YUC5_regulators.csv')
assays<-fread(file.path(ROOT,'inputs/experimental_DOF_evidence.csv'))
base<-merge(read('paired_TF_baseline_covariates.csv'),prior,by='gene',all.x=TRUE);base[is.na(prior_target_count),prior_target_count:=0]
fitout<-list()
for (label in c('WT_specific','fer4_specific')) {
 fg<-unique(read(paste0(label,'_TF_DE_repaired.csv'))$gene)
 dd<-copy(base);dd[,`:=`(selected=as.integer(gene%in%fg),Dof=as.integer(Family=='Dof'))];dd<-dd[complete.cases(dd)]
 fit<-glm(selected~Dof+log1p(prior_target_count)+log1p(WT_median_outdegree)+log1p(WT_median_CPM)+log1p(scored_contexts),data=dd,family=binomial())
 cf<-coef(summary(fit))['Dof',]
 fitout[[label]]<-data.table(genotype=label,n_genes=nrow(dd),n_Dof=sum(dd$Dof),odds_ratio=exp(cf[1]),Wald_lower=exp(cf[1]-1.96*cf[2]),Wald_upper=exp(cf[1]+1.96*cf[2]),p_Wald=cf[4],interpretation='Exploratory gene-level conditional association; not a confirmatory test or causal estimate')
 write(as.data.table(coef(summary(fit)),keep.rownames='term'),paste0('exploratory_prior_adjusted_model_',label,'.csv'))
}
write(rbindlist(fitout),'DOF_prior_degree_adjusted_sensitivity.csv');print(rbindlist(fitout))
write(base,'paired_TF_covariates_including_prior.csv')
# Relationship eligibility is explicit at each existing network context.
scores<-read('DOF_YUC5_all_network_scores.csv');ctx<-unique(scores[,.(tissue,cell_type,genotype,cluster)])
expected<-as.data.table(merge(as.data.frame(assays[,.(gene,symbol,alias)]),as.data.frame(ctx),by=NULL))
yucelig<-read('all_target_gene_eligibility.csv')[gene=='AT5G43890',.(tissue,YUC5_legacy_eligible=legacy_eligible,YUC5_explicit_eligible=explicit_eligible,YUC5_actual_input=actual_input)]
# An absent eligibility row means excluded by both reconstructions.
expected<-merge(expected,yucelig,by='tissue',all.x=TRUE)
for (k in c('YUC5_legacy_eligible','YUC5_explicit_eligible','YUC5_actual_input'))set(expected,which(is.na(expected[[k]])),k,FALSE)
expected[,supplied_prior_contains_DOF_to_YUC5:=gene%in%yucprior$source]
scores_sel<-scores[,.(gene,cluster,DOF_degree_out=degree_out,DOF_degree_centrality_out=degree_centrality_out)]
expected<-merge(expected,scores_sel,by=c('gene','cluster'),all.x=TRUE)
edges<-read('all_retained_YUC5_edges.csv');setnames(edges,'source','gene')
expected<-merge(expected,edges[,.(gene,cluster,coef_mean,coef_abs,p,edge_target_excluded_by_explicit_map=target_excluded_by_explicit_map)],by=c('gene','cluster'),all.x=TRUE)
expected[,relationship_status:=fifelse(!YUC5_actual_input,'YUC5 excluded before fitting',fifelse(!supplied_prior_contains_DOF_to_YUC5,'Relationship absent from supplied prior',fifelse(!is.na(coef_mean),'Retained fitted edge','Eligible supplied relationship; no retained edge')))]
DE<-read('DOF_YUC5_current_DE_expression.csv')
expected<-merge(expected,DE[,.(gene,cell_type,DOF_logFC=logFC,DOF_local_FDR=p_adj.loc,DOF_current_DE_status=current_DE_status,DOF_WT_mean_CPM=WT_mean_CPM,DOF_fer4_mean_CPM=fer4_mean_CPM,DOF_WT_frequency=WT.frq,DOF_fer4_frequency=fer_4.frq)],by=c('gene','cell_type'),all.x=TRUE)
write(expected,'DOF_YUC5_relationship_context_evidence.csv')
# One row per experimentally tested DOF; stage-specific detail remains separately linked.
network_summary<-expected[,.(supplied_prior_YUC5=any(supplied_prior_contains_DOF_to_YUC5),network_contexts=nrow(.SD),contexts_YUC5_in_input=sum(YUC5_actual_input),contexts_retained_YUC5_edge=sum(!is.na(coef_mean)),tissues_retained_YUC5_edge=paste(sort(unique(tissue[!is.na(coef_mean)])),collapse=';'),contexts_excluded_YUC5=sum(!YUC5_actual_input)),by=gene]
DEsummary<-DE[gene%in%assays$gene,.(basal_called_DE_contexts=paste(cell_type[current_called==TRUE],collapse=';'),basal_detection_eligible_contexts=sum(eligible_detection),basal_tested_contexts=.N),by=gene]
nom<-read('TF_nomination_threshold_sensitivity.csv')[topN==100&tail_fraction==.1&gene%in%assays$gene]
nom<-nom[,.(nominated_WT_contexts=paste(cell_type[Label=='WT_specific'],collapse=';'),nominated_fer4_contexts=paste(cell_type[Label=='fer4_specific'],collapse=';')),by=gene]
integrated<-Reduce(function(a,b)merge(a,b,by='gene',all.x=TRUE),list(assays,network_summary,DEsummary,nom))
write(integrated,'tested_DOF_integrated_evidence_matrix.csv')
# Assay coverage graphic uses evidence presence and untested states, not significance stars.
mat<-rbindlist(lapply(seq_len(nrow(assays)),function(i){d<-assays[i];data.table(symbol=d$symbol,assay=c('YUC5 promoter\nactivation','YUC5 fragment\nbinding in vitro','WT RALF23\ntranscript induction','RALF23 growth\nmutant assay','YUC5 induction\nin DOF mutant'),status=c('Observed',if(d$EMSA=='Not tested')'Not tested'else'Observed',if(grepl('No significant',d$root_RALF_qPCR))'No induction shown'else'Observed',if(d$mutant_RALF_growth=='Not tested')'Not tested'else'Tested; contrast unresolved','Not tested'))}))
mat[,symbol:=factor(symbol,levels=rev(assays$symbol))];mat[,assay:=factor(assay,levels=unique(assay))]
p<-ggplot(mat,aes(assay,symbol,fill=status))+geom_tile(colour='white',linewidth=1)+scale_fill_manual(values=c('Observed'='#0072B2','No induction shown'='#E69F00','Not tested'='#EEEEEE','Tested; contrast unresolved'='#999999'))+labs(title='The DOF evidence is complementary across assays',subtitle='No tested DOF has a measured endogenous YUC5 induction defect.\nDOF3.2 and DOF6 are aliases for AT3G45610.',x=NULL,y=NULL,fill=NULL)+theme_cowplot(font_size=11)+theme(legend.position='bottom',axis.text.x=element_text(size=10),plot.title=element_text(size=14),plot.subtitle=element_text(size=10),plot.title.position='plot')
ggsave(file.path(ROOT,'figures/06_DOF_assay_evidence.pdf'),p,width=9,height=5.7,device=cairo_pdf)
ggsave(file.path(ROOT,'figures/06_DOF_assay_evidence.png'),p,width=9,height=5.7,dpi=220,bg='white')
writeLines(capture.output(sessionInfo()),file.path(ROOT,'logs/R_prior_integration_sessionInfo.txt'))
cat('DOF evidence integration finished.\n')
