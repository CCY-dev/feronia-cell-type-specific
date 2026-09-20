#!/usr/bin/env Rscript
# Review v2: plain-language count headers. All estimates and statistical definitions unchanged.
suppressPackageStartupMessages({library(ggplot2);library(cowplot);library(colorspace)})
options(device=function(...) grDevices::pdf(file=NULL,...))
pdfFonts(Arial=pdfFonts()$Helvetica)
argv <- commandArgs(FALSE)
out <- dirname(normalizePath(sub('^--file=','',argv[grepl('^--file=',argv)][1])))
repo <- normalizePath(file.path(out,'../../../..'))
family_src <- file.path(repo,'revisions/2026-09-20_compact_figures/quantitative/source_top10_family.csv')
file.copy(family_src,file.path(out,'source_top10_family.csv'),overwrite=TRUE)
z <- read.csv(family_src,stringsAsFactors=FALSE)
stopifnot(nrow(z)==20, length(unique(z$Family))==10, sum(z$p_adj_BH<.05)==6, sum(z$odds_ratio==0)==1)
cols <- c(WT='#0072B2',fer='#D55E00')
palette_for <- function(cvd=NULL) if(is.null(cvd))cols else setNames(get(cvd,asNamespace('colorspace'))(cols),names(cols))
pt <- 72.27/25.4
theme_set(theme_cowplot(font_family='Arial',font_size=7.5)+
 theme(plot.margin=margin(1,1,1,1),axis.text=element_text(size=7),
       plot.title=element_text(size=8,face='bold'),axis.title=element_text(size=7.5)))
families<-z$Family[order(z$display_rank)];families<-unique(families)
lab_family<-function(x)ifelse(x=='Dof','DOF',gsub('_',' ',x))
base_y<-function(n=10)scale_y_continuous(limits=c(.5,n+.5),expand=c(0,0),breaks=1:n,labels=NULL)
forest_panel<-function(genotype,cvd=NULL,dat=z){
  a<-dat[dat$genotype==genotype,]; n<-length(unique(dat$Family))
  col<-unname(palette_for(cvd)[if(genotype=='WT_specific')'WT' else 'fer'])
  title<-if(genotype=='WT_specific')'WT-associated' else 'fer-4-associated'
  p<-ggplot(a,aes(y=row))+
    geom_vline(xintercept=1,linetype='dashed',size=.3,colour='#999999')+
    geom_segment(aes(x=pmax(odds_lower,.125),xend=odds_upper,yend=row),colour=col,size=.45)+
    geom_point(data=a[a$odds_ratio>0,],aes(x=odds_ratio,shape=p_adj_BH<.05),colour=col,size=1.75,stroke=.55)+
    scale_shape_manual(values=c(`TRUE`=16,`FALSE`=1),guide='none')+
    scale_x_log10(limits=c(.125,24),breaks=c(.25,1,4,16),labels=c('0.25','1','4','16'),expand=c(0,0))+
    base_y(n)+labs(x=NULL,y=NULL,title=title)+
    theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),
          panel.grid.major.x=element_line(colour='#EEEEEE',size=.2),
          plot.title=element_text(size=7.5,colour=col),
          axis.ticks.length=grid::unit(1,'mm'),axis.text.x=element_text(margin=margin(t=1.5)))
  if(any(a$odds_ratio==0))p<-p+
    geom_segment(data=a[a$odds_ratio==0,],aes(x=.21,xend=.127,yend=row),colour=col,size=.45,
       arrow=grid::arrow(length=grid::unit(1.2,'mm'),type='closed'))+
    geom_text(data=a[a$odds_ratio==0,],aes(x=.16,label='0'),hjust=0,vjust=-.65,size=7/pt,family='Arial')
  p
}
count_panel<-function(genotype,dat=z){
  a<-dat[dat$genotype==genotype,];n<-length(unique(dat$Family))
  ggplot(a,aes(x=1,y=row,label=count))+geom_text(hjust=1,size=7/pt,family='Arial')+
    base_y(n)+scale_x_continuous(limits=c(0,1.05),expand=c(0,0))+
    labs(title='Nominated /\neligible',x=NULL,y=NULL)+theme_void(base_family='Arial')+
    theme(plot.title=element_text(size=7,hjust=1,lineheight=.95),plot.margin=margin(1,0,1,0))
}
family_panel<-function(dat=z){
  a<-dat[!duplicated(dat$Family),];n<-nrow(a)
  ggplot(a,aes(x=1,y=row,label=lab_family(Family)))+
    geom_text(hjust=1,size=7.5/pt,family='Arial')+base_y(n)+
    scale_x_continuous(limits=c(0,1.025),expand=c(0,0))+theme_void(base_family='Arial')+
    labs(title='TF family')+
    theme(plot.title=element_text(size=7.5,hjust=1),plot.margin=margin(1,1,1,0))
}
forest_plot<-function(cvd=NULL,dat=z){
  panels<-plot_grid(family_panel(dat),forest_panel('WT_specific',cvd,dat),count_panel('WT_specific',dat),
    forest_panel('fer4_specific',cvd,dat),count_panel('fer4_specific',dat),
    nrow=1,rel_widths=c(.16,.24,.13,.24,.13),align='h',axis='tb')
  ggdraw()+draw_plot(panels,x=.002,y=.13,width=.996,height=.87)+
    draw_label('Odds ratio (95% CI)',x=.585,y=.086,size=7.5,fontfamily='Arial')+
    draw_label('Filled: BH < 0.05   Open: BH >= 0.05',x=.5,y=.026,size=7,fontfamily='Arial')
}

save_all<-function(p,stem,wmm,hmm){
  w<-wmm/25.4;h<-hmm/25.4
  grDevices::cairo_pdf(file=file.path(out,paste0(stem,'.pdf')),width=w,height=h,bg='white')
  print(p);dev.off()
  ggsave(file.path(out,paste0(stem,'.png')),p,width=w,height=h,dpi=300,bg='white')
}
save_all(forest_plot(),'top10_aligned_forest_review_v2',110,76)
sigfam<-unique(z$Family[z$p_adj_BH<.05]);zz<-z[z$Family%in%sigfam,]
sigord<-families[families%in%sigfam];zz$row<-length(sigord)+1-match(zz$Family,sigord)
write.csv(zz,file.path(out,'plot_source_four_significant_families.csv'),row.names=FALSE)
save_all(forest_plot(dat=zz),'four_families_review_v2',110,43)
for(v in c('deutan','protan','tritan')) {
  ggsave(file.path(out,paste0('CVD_',v,'_top10.png')),forest_plot(v),width=110/25.4,height=76/25.4,dpi=180,bg='white')
}
if(requireNamespace('colorblindr',quietly=TRUE)) colorblindr::cvd_grid(forest_plot())
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
cat('Review v2 PNG/PDF/SVG ready. Count headers now Nominated / eligible; source data unchanged.\n')
