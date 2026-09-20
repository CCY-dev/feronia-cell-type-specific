#!/usr/bin/env Rscript
# Prototype; not yet carefully vetted.
# Bounded follow-up after source repairs: identical-case adjustment and omission
# of both target-mapping-affected tissues, without any network refit.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(cowplot)})
setDTthreads(2);ROOT<-commandArgs(TRUE)[1]
read<-function(f)fread(file.path(ROOT,'results',f))
write<-function(d,n)fwrite(d,file.path(ROOT,'results',n),na='NA')
base<-read('paired_TF_covariates_including_prior.csv')
models<-list('Crude, identical genes'=selected~Dof,'Prior coverage only'=selected~Dof+log1p(prior_target_count),'Baseline degree only'=selected~Dof+log1p(WT_median_outdegree),'Expression and contexts'=selected~Dof+log1p(WT_median_CPM)+log1p(scored_contexts),'All four covariates'=selected~Dof+log1p(prior_target_count)+log1p(WT_median_outdegree)+log1p(WT_median_CPM)+log1p(scored_contexts))
fits<-list()
for (lab in c('WT_specific','fer4_specific')) {
 fg<-unique(read(paste0(lab,'_TF_DE_repaired.csv'))$gene)
 d<-copy(base);d[,`:=`(Dof=as.integer(Family=='Dof'),selected=as.integer(gene%in%fg))];d<-d[complete.cases(d)]
 for(m in names(models)) {
  f<-glm(models[[m]],data=d,family=binomial());co<-coef(summary(f))['Dof',]
  fits[[paste(lab,m)]]<-data.table(genotype=lab,model=m,n_genes=nrow(d),n_Dof=sum(d$Dof),Dof_foreground=sum(d$Dof&d$selected),foreground_n=sum(d$selected),OR=exp(co[1]),CI_low=exp(co[1]-1.96*co[2]),CI_high=exp(co[1]+1.96*co[2]),p_Wald=co[4])
 }
}
out<-rbindlist(fits);write(out,'DOF_identical_case_covariate_sensitivity.csv');print(out)
d<-copy(out);d[,model:=factor(model,levels=rev(names(models)))];d[,genotype:=factor(genotype,levels=c('WT_specific','fer4_specific'),labels=c('Higher rank in WT','Higher rank in fer-4'))]
p<-ggplot(d,aes(OR,model))+geom_vline(xintercept=1,linetype=2,colour='grey55')+geom_errorbarh(aes(xmin=CI_low,xmax=CI_high),height=.12,colour='#0072B2')+geom_point(size=2.5,colour='#0072B2')+facet_wrap(~genotype,ncol=1)+scale_x_log10()+labs(title='Prior coverage changes the apparent DOF enrichment',subtitle='Same 1,213 TFs and 32 DOFs in every model; exploratory logistic models.\nIntervals are Wald 95% intervals, with no claim of causal adjustment.',x='DOF odds ratio (log scale)',y=NULL)+theme_cowplot(font_size=11)+theme(plot.title.position='plot',plot.title=element_text(size=14),plot.subtitle=element_text(size=10))
ggsave(file.path(ROOT,'figures/07_DOF_covariate_sensitivity.pdf'),p,width=8.3,height=7,device=cairo_pdf,bg='white');ggsave(file.path(ROOT,'figures/07_DOF_covariate_sensitivity.png'),p,width=8.3,height=7,dpi=220,bg='white')
# Omit phloem and xylem entirely, rederive pooled metric thresholds among the
# eight unaffected tissue fits. This is not a corrected vascular-network refit.
w<-read('network_rank_table_DE_repaired.csv')[!crude_cell_type%in%c('phloem','xylem')]
w[,rank_diff:=get('fer4_untreated-WT_untreated: rank_diff')]
thr<-w[,.(low=as.numeric(quantile(rank_diff,.1)),high=as.numeric(quantile(rank_diff,.9))),by=metric]
w<-merge(w,thr,by='metric');w[,unaffected_label:=fifelse(rank_diff<low,'fer4_specific',fifelse(rank_diff>high,'WT_specific','None'))]
sel<-w[metric=='degree_centrality_out'&Is_TF=='Yes'&unaffected_label!='None'];write(sel,'TF_nominations_eight_unaffected_tissues.csv');write(thr,'eight_tissue_rank_thresholds.csv')
sources<-read('source_file_manifest.csv')$path
scorepaths<-sources[grepl('gene_score_iGRN.csv$',sources)&!grepl('/(phloem|xylem)_',sources)]
scores<-rbindlist(lapply(scorepaths,function(p){z<-fread(p);setnames(z,names(z)[1],'gene');z[,.(gene,cluster,genotype=sub('.*_(WT|fer4)_untreated$','\\1',cluster),cell_type=sub('_(WT|fer4)_.*','',cluster),degree_centrality_out)]}))
paired<-dcast(scores,gene+cell_type~genotype,value.var='degree_centrality_out');pairedgenes<-unique(paired[is.finite(WT)&is.finite(fer4),gene])
tf<-fread(sources[grepl('Kay_TF_thalemine_annotations.csv$',sources)]);wanted<-tf[[1]]
fam<-unique(fread(sources[grepl('Ath_TF_list.txt$',sources)])[,.(gene=Gene_ID,Family)]);setorder(fam,gene,Family);fam<-fam[,.(Family=Family[1]),by=gene]
u<-fam[gene%in%intersect(pairedgenes,wanted)];write(u,'eight_tissue_paired_TF_universe.csv')
ans<-list();counts<-list()
for(lab in c('WT_specific','fer4_specific')) {
 g<-unique(sel[unaffected_label==lab,gene]);orig<-unique(read(paste0(lab,'_TF_DE_repaired.csv'))$gene)
 counts[[lab]]<-data.table(genotype=lab,TFs=length(g),overlap_submitted=length(intersect(g,orig)),original_TFs=length(orig),Jaccard=length(intersect(g,orig))/length(union(g,orig)),Dof_count=sum(g%in%u[Family=='Dof',gene]))
 fgs<-intersect(g,u$gene);N<-nrow(u);n<-length(fgs)
 r<-rbindlist(lapply(sort(unique(u$Family)),function(f){K<-sum(u$Family==f);k<-sum(u$Family==f&u$gene%in%fgs);mat<-matrix(c(k,n-k,K-k,N-n-K+k),nrow=2,byrow=TRUE);ft<-fisher.test(mat,alternative='greater');data.table(Family=f,K=K,k=k,N=N,n=n,OR=unname(ft$estimate),p=ft$p.value)}));r[,`:=`(BH=p.adjust(p,'BH'),genotype=lab)];ans[[lab]]<-r
}
write(rbindlist(counts),'eight_tissue_nomination_sensitivity_summary.csv');write(rbindlist(ans),'eight_tissue_family_enrichment.csv');print(rbindlist(counts));print(rbindlist(ans)[Family=='Dof'])
# Final validation summaries avoid ambiguous denominators.
fix<-read('network_rank_table_DE_repaired.csv')
ctx<-unique(fix[,.(gene,cell_type,legacy_log2FC,log2FC,legacy_DE_label,current_label_same_legacy_rule,DE_label,current_DE_status,eligible_detection)])
labelsummary<-ctx[,.(n=.N,comparable=sum(is.finite(legacy_log2FC)&is.finite(log2FC)),FC_changed_gt_0.01=sum(abs(legacy_log2FC-log2FC)>.01,na.rm=TRUE),same_rule_label_changes=sum(legacy_DE_label!=current_label_same_legacy_rule,na.rm=TRUE),final_rule_label_changes=sum(legacy_DE_label!=DE_label,na.rm=TRUE)),by=current_DE_status];write(labelsummary,'annotation_changes_by_detection_status.csv')
e<-read('network_edges_affected_by_target_mapping.csv');write(e[,.(networks=.N,all_edges=sum(edges),affected_edges=sum(edges_to_excluded_targets),minimum_percent=100*min(edges_to_excluded_targets/edges),maximum_percent=100*max(edges_to_excluded_targets/edges)),by=tissue],'target_mapping_effect_summary.csv')
y<-read('all_retained_YUC5_edges.csv');write(y[,.(retained_edges=.N,sources=uniqueN(source),contexts=uniqueN(cluster)),by=tissue],'YUC5_retained_edge_summary.csv')
mat<-read('DOF_YUC5_relationship_context_evidence.csv');write(mat[,.(contexts=.N),by=.(symbol,relationship_status)],'DOF_relationship_status_summary.csv')
writeLines(capture.output(sessionInfo()),file.path(ROOT,'logs/R_validation_sessionInfo.txt'))
cat('Final bounded q04 sensitivity and validation complete.\n')
