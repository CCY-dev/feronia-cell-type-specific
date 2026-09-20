#!/usr/bin/env Rscript
# Prototype; not yet carefully vetted.
# Design-only revision. All source values, orders and analytical definitions are retained.
suppressPackageStartupMessages({library(ggplot2);library(cowplot);library(colorspace)})
options(device=function(...) grDevices::pdf(file=NULL,...))
pdfFonts(Arial=pdfFonts()$Helvetica)
argv <- commandArgs(FALSE)
out <- dirname(normalizePath(sub('^--file=','',argv[grepl('^--file=',argv)][1])))
repo <- normalizePath(file.path(out,'../../..'))
src <- file.path(repo,'revisions/2026-09-20_fig6_brainstorm')
cell_src <- file.path(src,'celltype_prototypes/celltype_summary_all37.csv')
family_src <- file.path(src,'top10_family_prototypes/plot_input_top10.csv')
file.copy(cell_src,file.path(out,'source_celltype_summary_all37.csv'),overwrite=TRUE)
file.copy(family_src,file.path(out,'source_top10_family.csv'),overwrite=TRUE)
d <- read.csv(cell_src,stringsAsFactors=FALSE)
q <- d[d$tested,]
z <- read.csv(family_src,stringsAsFactors=FALSE)
stopifnot(nrow(q)==130, length(unique(q$gene))==10, nrow(z)==20,
          all(q$n_WT+q$n_fer<=q$K),
          all(q$K+q$n_one_genotype+q$n_absent==q$n_zone),
          sum(z$p_adj_BH<.05)==6, sum(z$odds_ratio==0)==1)
tissues <- c('Atrichoblast','Trichoblast','Cortex','Endodermis','Xylem Pole Pericycle',
 'Phloem Pole Pericycle','Procambium','Protophloem','Metaphloem & Companion Cell',
 'Protoxylem','Metaxylem','Lateral Root Cap','Columella')
tlab <- c('Atri.','Tricho.','Cortex','Endo.','XPP','PPP','Procam.','Proto-\nphloem',
          'Meta-\nphloem/CC','Proto-\nxylem','Meta-\nxylem','LRC','Colum.')
ord <- c('DOF4.6','OBP1','OBP2','OBP3','DAG1','COG1','CDF4','DOF2.2','DOF6','DOF4.1')
q$xx <- match(q$tissue,tissues);q$yy <- match(q$symbol,rev(ord))
stopifnot(!anyNA(q$xx),!anyNA(q$yy))
write.csv(q,file.path(out,'plot_source_celltype_tested10.csv'),row.names=FALSE)
cols <- c(WT='#0072B2',fer='#D55E00')
palette_for <- function(cvd=NULL) if(is.null(cvd))cols else setNames(get(cvd,asNamespace('colorspace'))(cols),names(cols))
pt <- 72.27/25.4

cell_plot <- function(cvd=NULL){
  pal <- palette_for(cvd)
  a<-q[q$n_WT>0,];b<-q[q$n_fer>0,]
  lab<-rev(ord);lab[lab=='DOF6']<-'DOF3.2'
  g<-ggplot(q,aes(xx,yy))+
    geom_tile(aes(fill=K>0),width=.97,height=.94,colour='white',size=.2)+
    scale_fill_manual(values=c(`TRUE`='#E9E9E9',`FALSE`='white'),guide='none')+
    geom_rect(data=a,aes(xmin=xx-.485,xmax=xx+.485,ymin=yy,ymax=yy+.47),inherit.aes=FALSE,fill=pal['WT'])+
    geom_rect(data=b,aes(xmin=xx-.485,xmax=xx+.485,ymin=yy-.47,ymax=yy),inherit.aes=FALSE,fill=pal['fer'])+
    geom_text(data=a,aes(y=yy+.235,label=paste0(n_WT,'/',K)),size=7/pt,colour='white',family='Arial')+
    geom_text(data=b,aes(y=yy-.235,label=paste0(n_fer,'/',K)),size=7/pt,colour='white',family='Arial')+
    scale_x_continuous(breaks=1:13,labels=tlab,position='top',expand=c(0,0),limits=c(.5,13.5))+
    scale_y_continuous(breaks=1:10,labels=lab,expand=c(0,0),limits=c(.5,10.5))+
    labs(x=NULL,y=NULL)+theme_map(font_family='Arial',font_size=7.5)+
    theme(axis.text.x=element_text(size=7.2,colour='black',lineheight=.92,margin=margin(b=2)),
          axis.text.y=element_text(size=8,colour='black',face='italic',margin=margin(r=3)),
          axis.ticks=element_blank(),plot.margin=margin(0,1,0,0))
  if(any(q$K==0))g<-g+geom_text(data=q[q$K==0,],label='x',size=7/pt,colour='#666666')
  # A short positional key preserves meaning independently of colour.
  ggdraw()+draw_plot(g,x=0,y=.095,width=1,height=.905)+
    draw_label('Upper: WT-associated',x=.005,y=.035,hjust=0,size=7,fontfamily='Arial',colour=pal['WT'])+
    draw_label('Lower: fer-4-associated',x=.235,y=.035,hjust=0,size=7,fontfamily='Arial',colour=pal['fer'])+
    draw_label('n/K zones',x=.49,y=.035,hjust=0,size=7,fontfamily='Arial')+
    draw_label('Grey: no nomination   x: no paired scores',x=.62,y=.035,hjust=0,size=7,fontfamily='Arial')
}

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
    labs(title='n/K',x=NULL,y=NULL)+theme_void(base_family='Arial')+
    theme(plot.title=element_text(size=7.5,hjust=1),plot.margin=margin(1,0,1,0))
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
    nrow=1,rel_widths=c(.19,.265,.085,.265,.085),align='h',axis='tb')
  ggdraw()+draw_plot(panels,x=.002,y=.13,width=.996,height=.87)+
    draw_label('Odds ratio (95% CI)',x=.585,y=.086,size=7.5,fontfamily='Arial')+
    draw_label('Filled: BH < 0.05   Open: BH >= 0.05',x=.5,y=.026,size=7,fontfamily='Arial')
}

save_all<-function(p,stem,wmm,hmm){
  w<-wmm/25.4;h<-hmm/25.4
  grDevices::quartz(type='pdf',file=file.path(out,paste0(stem,'.pdf')),width=w,height=h,family='Arial',bg='white')
  print(p);dev.off()
  ggsave(file.path(out,paste0(stem,'.svg')),p,device=svglite::svglite,width=w,height=h,bg='white')
  ggsave(file.path(out,paste0(stem,'.png')),p,width=w,height=h,dpi=300,bg='white')
}
save_all(cell_plot(),'celltype_tested10_compact',170,68)
save_all(forest_plot(),'top10_aligned_forest_compact',110,76)
# Same statistics, subset to families with BH < 0.05 in at least one genotype.
sigfam<-unique(z$Family[z$p_adj_BH<.05]);zz<-z[z$Family%in%sigfam,]
sigord<-families[families%in%sigfam];zz$row<-length(sigord)+1-match(zz$Family,sigord)
write.csv(zz,file.path(out,'plot_source_four_significant_families.csv'),row.names=FALSE)
save_all(forest_plot(dat=zz),'four_families_compact',110,43)
for(v in c('deutan','protan','tritan')){
  ggsave(file.path(out,paste0('CVD_',v,'_celltype.png')),cell_plot(v),width=170/25.4,height=68/25.4,dpi=160,bg='white')
  ggsave(file.path(out,paste0('CVD_',v,'_top10.png')),forest_plot(v),width=110/25.4,height=76/25.4,dpi=180,bg='white')
}
if(requireNamespace('colorblindr',quietly=TRUE))colorblindr::cvd_grid(cell_plot())
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
writeLines(c('Design-only revision: no fitting, selection or significance changes.',
 'celltype_tested10_compact: 170 x 68 mm; minimum label 7 pt.',
 'top10_aligned_forest_compact: 110 x 76 mm; minimum label 7 pt.',
 'four_families_compact: 110 x 43 mm; optional same-statistics subset.',
 'Source scores and nominations are inherited. Endodermis/phloem/xylem audit issues remain unresolved.'),
 file.path(out,'VALIDATION.txt'))
cat('Compact PNG/PDF/SVG ready. 130 cell-type tiles; 20 family estimates; source values unchanged.\n')
