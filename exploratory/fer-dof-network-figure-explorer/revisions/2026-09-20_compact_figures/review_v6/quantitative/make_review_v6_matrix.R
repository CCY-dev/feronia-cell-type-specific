#!/usr/bin/env Rscript
# PROTOTYPE: not yet carefully vetted. Design/data checks do not resolve inherited analysis limitations.
# Design-only revision. All source values, orders and analytical definitions are retained.
suppressPackageStartupMessages({library(ggplot2);library(cowplot);library(colorspace)})
options(device=function(...) grDevices::pdf(file=NULL,...))
pdfFonts(Arial=pdfFonts()$Helvetica)
argv <- commandArgs(FALSE)
out <- dirname(normalizePath(sub('^--file=','',argv[grepl('^--file=',argv)][1])))
repo <- normalizePath(file.path(out,'../../../..'))
src <- file.path(repo,'revisions/2026-09-20_fig6_brainstorm')
cell_src <- file.path(src,'celltype_prototypes/celltype_summary_all37.csv')
file.copy(cell_src,file.path(out,'source_celltype_summary_all37.csv'),overwrite=TRUE)
d <- read.csv(cell_src,stringsAsFactors=FALSE)
q <- d[d$tested,]
stopifnot(nrow(q)==130, length(unique(q$gene))==10, length(unique(q$tissue))==13,
          all(q$n_WT+q$n_fer<=q$K),
          all(q$K+q$n_one_genotype+q$n_absent==q$n_zone))
tissues <- c('Atrichoblast','Trichoblast','Cortex','Endodermis','Xylem Pole Pericycle',
 'Phloem Pole Pericycle','Procambium','Protophloem','Metaphloem & Companion Cell',
 'Protoxylem','Metaxylem','Lateral Root Cap','Columella')
tlab <- c('Atrichoblast','Trichoblast','Cortex','Endodermis','Xylem pole pericycle',
          'Phloem pole pericycle','Procambium','Protophloem',
          'Metaphloem and\ncompanion cell','Protoxylem','Metaxylem','Lateral root cap','Columella')
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
    geom_text(data=a,aes(y=yy+.235,label=paste0(n_WT,'/',K)),size=8.5/pt,colour='white',family='Arial')+
    geom_text(data=b,aes(y=yy-.235,label=paste0(n_fer,'/',K)),size=8.5/pt,colour='white',family='Arial')+
    scale_x_continuous(breaks=1:13,labels=tlab,position='bottom',expand=c(0,0),limits=c(.5,13.5))+
    scale_y_continuous(breaks=1:10,labels=lab,expand=c(0,0),limits=c(.5,10.5))+
    labs(x=NULL,y=NULL)+theme_map(font_family='Arial',font_size=8.5)+
    theme(axis.text.x=element_text(size=8.5,colour='black',angle=60,hjust=1,vjust=1,lineheight=.94,margin=margin(t=3)),
          axis.text.y=element_text(size=9,colour='black',face='italic',margin=margin(r=3)),
          axis.ticks=element_blank(),plot.margin=margin(0,1,0,0))
  if(any(q$K==0))g<-g+geom_text(data=q[q$K==0,],label='x',size=9/pt,colour='#666666')
  # Full bottom labels are angled at the user's explicit request.
  # Shift only the wrapped long label left by 1 mm to clear Protoxylem.
  nudge_long_label <- function(z) {
    if(inherits(z,'text') && !is.null(z$label)) {
      i <- which(z$label=='Metaphloem and\ncompanion cell')
      if(length(i))z$x[i] <- z$x[i]-grid::unit(1,'mm')
    }
    if(!is.null(z$children))z$children <- do.call(grid::gList,lapply(z$children,nudge_long_label))
    if(!is.null(z$grobs))z$grobs <- lapply(z$grobs,nudge_long_label)
    z
  }
  g <- nudge_long_label(ggplotGrob(g))
  grey_key <- if(any(q$K==0))'Grey: no nomination   x: no paired scores' else 'Grey: no nomination'
  # Canvas width matches panel A; count font and true tile width are checked.
  ggdraw()+draw_plot(g,x=.008,y=.12,width=.978,height=.873)+
    draw_label('Upper: WT-associated',x=.006,y=.075,hjust=0,size=8,fontfamily='Arial',colour=pal['WT'])+
    draw_label('Lower: fer-4-associated',x=.50,y=.075,hjust=0,size=8,fontfamily='Arial',colour=pal['fer'])+
    draw_label('Nominated / scored zones',x=.006,y=.025,hjust=0,size=8,fontfamily='Arial')+
    draw_label(grey_key,x=.50,y=.025,hjust=0,size=8,fontfamily='Arial')
}

save_all<-function(p,stem,wmm,hmm){
  # Evaluate grob construction before opening Quartz, avoiding a blank first page.
  force(p)
  w<-wmm/25.4;h<-hmm/25.4
  grDevices::quartz(type='pdf',file=file.path(out,paste0(stem,'.pdf')),width=w,height=h,family='Arial',bg='white')
  print(p);dev.off()
  ggsave(file.path(out,paste0(stem,'.svg')),p,device=svglite::svglite,width=w,height=h,bg='white')
  ggsave(file.path(out,paste0(stem,'.png')),p,width=w,height=h,dpi=300,bg='white')
}
save_all(cell_plot(),'celltype_tested10_review_v6',110,90)
for(v in c('deutan','protan','tritan')) {
  ggsave(file.path(out,paste0('CVD_',v,'_celltype.png')),cell_plot(v),width=110/25.4,height=90/25.4,dpi=160,bg='white')
}
if(requireNamespace('colorblindr',quietly=TRUE))colorblindr::cvd_grid(cell_plot())
cat('Matrix review v6 ready. Nominated / scored zones; all data unchanged.\n')
