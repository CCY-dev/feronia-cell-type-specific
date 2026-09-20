#!/usr/bin/env Rscript
# Prototype; not yet carefully vetted.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(cowplot)})
ROOT<-commandArgs(TRUE)[1];x<-readRDS(file.path(ROOT,'results/q04_figure_data.rds'))
read<-function(f)fread(file.path(ROOT,'results',f))
theme_set(theme_cowplot(font_size=11)+theme(plot.title=element_text(size=14,face='bold'),plot.subtitle=element_text(size=10),strip.background=element_blank(),strip.text=element_text(size=10),legend.position='bottom',plot.title.position='plot'))
saveplot<-function(p,n,w,h){p$labels$title<-paste(strwrap(p$labels$title,width=65),collapse='\n');p$labels$subtitle<-paste(strwrap(p$labels$subtitle,width=100),collapse='\n');ggsave(file.path(ROOT,'figures',paste0(n,'.pdf')),p,width=w,height=h,device=cairo_pdf,bg='white');ggsave(file.path(ROOT,'figures',paste0(n,'.png')),p,width=w,height=h,dpi=220,bg='white')}
blue<-'#0072B2';orange<-'#D55E00'
# Show unique gene-context pairs, avoiding repetition over four network metrics.
d<-read('annotation_repair_gene_contexts.csv')[is.finite(log2FC)&is.finite(legacy_log2FC)]
p<-ggplot(d,aes(legacy_log2FC,log2FC))+geom_abline(slope=1,intercept=0,colour='grey65')+geom_point(size=.7,alpha=.25,colour=blue)+coord_equal()+labs(title='DE annotation repair preserves network scores',subtitle='Each point is one gene-by-cell-stage context; network scores are unchanged',x='Legacy log2 fold change',y='Current log2 fold change')
saveplot(p,'01_DE_annotation_repair',6,5.4)
# Family enrichment interpreted conditionally on the specified gene universe.
d<-x$enrich[Family=='Dof'&definition=='PlantTFDB'&universe!='legacy_selected_KayTF']
d[,universe_label:=factor(universe,levels=c('legacy_selected_wide','all_score_PlantTFDB','paired_score_KayTF'),labels=c('Published selected-gene background','All scored, classifiable TFs','Paired-score TFs eligible for nomination'))]
d[,genotype:=factor(genotype,levels=c('WT_specific','fer4_specific'),labels=c('Higher rank in WT','Higher rank in fer-4'))]
p<-ggplot(d,aes(odds_ratio,universe_label))+geom_vline(xintercept=1,linetype=2,colour='grey65')+geom_errorbarh(aes(xmin=odds_lower,xmax=odds_upper),height=.15,colour=blue)+geom_point(size=2.6,colour=blue)+facet_wrap(~genotype,ncol=1)+scale_x_log10()+labs(title='DOF enrichment persists across gene backgrounds',subtitle='Gene-level Fisher odds ratios and two-sided 95% intervals; BH tests are tabulated',x='Odds ratio (log scale)',y=NULL)
saveplot(p,'02_DOF_enrichment_backgrounds',8.5,5.8)
# Actual consequences of target-subgroup mismatch; this is not network re-inference.
d<-x$edge[tissue%in%c('phloem','xylem')];d[,context:=sub('_(WT|fer4)_untreated$','',cluster)];d[,genotype:=sub('.*_(WT|fer4)_untreated$','\\1',cluster)]
p<-ggplot(d,aes(100*edges_to_excluded_targets/edges,context,colour=genotype,shape=genotype))+geom_point(size=2.5,position=position_dodge(width=.4))+facet_wrap(~tissue,ncol=1,scales='free_y')+scale_colour_manual(values=c(WT=blue,fer4=orange))+labs(title='Pericycle-derived target eligibility reaches vascular networks',subtitle='Retained edges whose target enters only through the broad substring match',x='Edges to targets excluded by the explicit tissue map (%)',y=NULL,colour=NULL,shape=NULL)
saveplot(p,'03_target_membership_consequences',8.6,7)
# Basal DE context of the experimentally selected DOFs. Sparse detection is explicit.
genes<-data.table(gene=c('AT4G24060','AT3G50410','AT1G07640','AT3G55370','AT3G61850','AT1G29160','AT2G34140','AT2G28810','AT3G45610','AT4G00940'),symbol=c('DOF4.6','OBP1','OBP2','OBP3','DAG1','COG1','CDF4','DOF2.2','DOF3.2','DOF4.1'))
raw<-read('DOF_YUC5_current_DE_expression.csv')
grid<-CJ(gene=genes$gene,cell_type=sort(unique(raw$cell_type)))
d<-merge(merge(grid,raw,by=c('gene','cell_type'),all.x=TRUE),genes,by='gene')
stopifnot(nrow(d)==nrow(genes)*uniqueN(raw$cell_type))
d[,symbol:=factor(symbol,levels=genes$symbol)]
ord<-d[order(tissue,cell_type),unique(cell_type)];d[,cell_type:=factor(cell_type,levels=rev(ord))]
d[,display_FC:=fifelse(!is.na(eligible_detection)&eligible_detection,pmax(-2,pmin(2,logFC)),NA_real_)]
p<-ggplot(d,aes(symbol,cell_type))+geom_tile(aes(fill=display_FC),colour='white',linewidth=.2)+geom_point(data=d[is.na(logFC)],shape=4,size=1.6,colour='black')+geom_point(data=d[current_called==TRUE],shape=1,size=1.6,colour='black')+scale_fill_gradient2(low=blue,mid='white',high=orange,midpoint=0,limits=c(-2,2),na.value='#D9D9D9',name='fer-4 / WT\nlog2 FC')+labs(title='Basal DOF responses vary with cell context',subtitle='Grey: below 10% detection; x: unavailable DE result. Circle: current local-FDR DEG call.\nValues clipped at ±2; these samples received no RALF23 treatment.',x=NULL,y=NULL)+theme(axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),legend.position='right')
saveplot(p,'04_tested_DOF_expression_context',9.2,11)
# Sensitivity of nominated TF sets to the two arbitrary selection cutoffs.
d<-x$sensitivity;d[,genotype:=factor(genotype,levels=c('WT_specific','fer4_specific'),labels=c('Higher rank in WT','Higher rank in fer-4'))];d[,tail_fraction:=factor(tail_fraction,levels=c(.05,.1,.2),labels=c('5% tails','10% tails','20% tails'))]
p<-ggplot(d,aes(topN,Jaccard,colour=tail_fraction,shape=tail_fraction))+geom_line()+geom_point(size=2.3)+facet_wrap(~genotype)+scale_x_continuous(breaks=c(50,100,200))+scale_colour_manual(values=c('#0072B2','#000000','#D55E00'))+coord_cartesian(ylim=c(0,1))+labs(title='TF nominations change with rank-selection thresholds',subtitle='Jaccard overlap with the submitted top100 / 10% tail set; same fitted networks',x='Top genes retained per context and metric',y='Overlap with submitted TF set',colour=NULL,shape=NULL)
saveplot(p,'05_TF_nomination_threshold_sensitivity',8,4.8)
# Color vision checks use the same palettes; shape encodes genotype/threshold too.
checks<-c('Color keys use Okabe-Ito blue/orange plus black; group comparisons also use shape.','PDF exports use cairo_pdf; PNG exports are 220 dpi.','Heatmap grey denotes low detection, not zero expression or zero effect.')
if(requireNamespace('colorspace',quietly=TRUE)){
 pal<-c(blue,'#FFFFFF',orange,'#D9D9D9');checks<-c(checks,paste('Deutan:',paste(colorspace::deutan(pal),collapse=',')),paste('Protan:',paste(colorspace::protan(pal),collapse=',')),paste('Tritan:',paste(colorspace::tritan(pal),collapse=',')))
}else checks<-c(checks,'colorspace unavailable in inherited environment; no package installed; simulation not performed.')
writeLines(checks,file.path(ROOT,'logs/figure_design_checks.txt'))
writeLines(capture.output(sessionInfo()),file.path(ROOT,'logs/R_figures_sessionInfo.txt'))
