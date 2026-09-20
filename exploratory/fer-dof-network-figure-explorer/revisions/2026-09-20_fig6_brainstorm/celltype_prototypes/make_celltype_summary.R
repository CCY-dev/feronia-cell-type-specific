#!/usr/bin/env Rscript
# Prototype; not yet carefully vetted.
# Descriptive collapse of existing zone-resolved nominations. No network refit or significance test.
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(colorspace)})
pdfFonts(Arial=pdfFonts()$Helvetica)
options(device=function(...) grDevices::pdf(file=NULL,...))
arg <- commandArgs(trailingOnly=FALSE)
script <- normalizePath(sub('^--file=','',arg[grepl('^--file=',arg)][1]))
out <- dirname(script); repo <- normalizePath(file.path(out,'../../..'))
src <- file.path(repo,'.living/tasks/fer-fig6-brainstorm-2026-09-20/outputs/claude-science/data/dof_context_status.csv')
x <- read.csv(src,stringsAsFactors=FALSE)
x$tested <- tolower(x$tested)=='true';x$de_called <- tolower(x$de_called)=='true'
stopifnot(!anyDuplicated(x[c('gene','context')]), nrow(x)==37*44)
paired <- c('WT-associated','fer-4-associated','scored, not nominated')
d <- x %>% group_by(gene,symbol,tissue,tested) %>% summarise(
  K=sum(status %in% paired), n_WT=sum(status=='WT-associated'),
  n_fer=sum(status=='fer-4-associated'),n_zone=n(),
  n_one_genotype=sum(status=='scored in one genotype only'),
  n_absent=sum(status=='not in network'),
  n_DE_up=sum(current_DE_status=='called_up',na.rm=TRUE),
  n_DE_down=sum(current_DE_status=='called_down',na.rm=TRUE),.groups='drop') %>%
  mutate(f_WT=ifelse(K>0,n_WT/K,NA_real_),f_fer=ifelse(K>0,n_fer/K,NA_real_),
    mixed=n_WT>0 & n_fer>0)
stopifnot(all(d$n_WT+d$n_fer<=d$K),all(d$K+d$n_one_genotype+d$n_absent==d$n_zone))
nominated <- d %>% group_by(gene,symbol) %>% summarise(nom=sum(n_WT+n_fer),.groups='drop') %>%
  filter(nom>0) %>% arrange(desc(nom),symbol)
stopifnot(nrow(nominated)==26,sum(d$mixed)==7,sum(x$de_called)==14)
write.csv(d,file.path(out,'celltype_summary_all37.csv'),row.names=FALSE,na='NA')
write.csv(d[d$mixed,],file.path(out,'opposing_zone_nominations.csv'),row.names=FALSE)
tissues <- c('Atrichoblast','Trichoblast','Cortex','Endodermis','Xylem Pole Pericycle',
 'Phloem Pole Pericycle','Procambium','Protophloem','Metaphloem & Companion Cell',
 'Protoxylem','Metaxylem','Lateral Root Cap','Columella')
tlab <- c('Atri.','Tricho.','Cortex','Endo.','XPP','PPP','Procam.','Proto-\nphloem','Meta-\nphloem/CC',
          'Proto-\nxylem','Meta-\nxylem','LRC','Colum.');names(tlab)<-tissues
tested_order <- c('DOF4.6','OBP1','OBP2','OBP3','DAG1','COG1','CDF4','DOF2.2','DOF6','DOF4.1')
make_plot <- function(ten=FALSE,cvd=NULL) {
  ids <- if(ten) d$gene[d$tested] else nominated$gene
  q <- d[d$gene %in% ids,]
  ord <- if(ten) tested_order else nominated$symbol
  q$xx <- match(q$tissue,tissues);q$yy <- match(q$symbol,rev(ord))
  stopifnot(!anyNA(q$xx),!anyNA(q$yy))
  cols <- c(WT='#0072B2',fer='#D55E00')
  if(!is.null(cvd))cols[]<-get(cvd,asNamespace('colorspace'))(cols)
  a <- q[q$n_WT>0,]; b<-q[q$n_fer>0,]
  title <- if(ten)'Network nominations by cell type: the ten assayed DOFs' else 'Network nominations by cell type: all 26 nominated DOFs'
  lab<-rev(ord);lab[lab=='DOF6']<-'DOF3.2'
  ggplot(q,aes(xx,yy))+
    geom_tile(aes(fill=K>0),width=.94,height=.94,colour='white',size=.3)+
    scale_fill_manual(values=c(`TRUE`='#ECECEC',`FALSE`='#FFFFFF'),guide='none')+
    geom_rect(data=a,aes(xmin=xx-.47,xmax=xx+.47,ymin=yy,ymax=yy+.47),fill=cols['WT'],inherit.aes=FALSE)+
    geom_rect(data=b,aes(xmin=xx-.47,xmax=xx+.47,ymin=yy-.47,ymax=yy),fill=cols['fer'],inherit.aes=FALSE)+
    geom_text(data=a,aes(y=yy+.235,label=paste0(n_WT,'/',K)),size=2.35,colour='white')+
    geom_text(data=b,aes(y=yy-.235,label=paste0(n_fer,'/',K)),size=2.35,colour='white')+
    {if(any(q$K==0))geom_text(data=q[q$K==0,],label='x',size=2.7,colour='#777777')}+
    scale_x_continuous(breaks=seq_along(tissues),labels=tlab,position='top',expand=c(0,0),limits=c(.5,13.5))+
    scale_y_continuous(breaks=seq_along(ord),labels=lab,expand=c(0,0),limits=c(.5,length(ord)+.5))+
    labs(x=NULL,y=NULL,title=title,
      subtitle=paste('Upper blue half: WT-associated nominations. Lower orange half: fer-4-associated nominations.',
       'Labels = nominated zones / zones scored in both genotypes for that gene and cell type.',sep='\n'),
      caption=paste('Gray: paired scores without nomination. x: no zone with paired scores. No DE stars are overlaid.',
       'Both halves colored: opposing nomination directions occur in different zones. These are context counts, not replicate support.',
       'Existing top-100 / pooled 10% tail rule, applied to relative out-degree ranks. No new cell-type network was fitted.',
       'XPP/PPP: xylem-/phloem-pole pericycle; CC: companion cell; LRC: lateral root cap. Full missingness is retained in the CSV.',sep='\n'))+
    theme_minimal(base_size=10,base_family='Arial')+
    theme(panel.grid=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(size=8,colour='black'),
      axis.text.y=element_text(size=9,colour='black',face='italic'),
      plot.title=element_text(face='bold',size=12),plot.subtitle=element_text(size=8.5,lineheight=1.1),
      plot.caption=element_text(size=7.6,hjust=0,lineheight=1.1),plot.title.position='plot',
      plot.caption.position='plot',plot.margin=margin(8,8,8,8))
}
save <- function(p,stem,h) {
  device <- if(capabilities('aqua')) function(filename,width,height,bg,...) grDevices::quartz(type='pdf',file=filename,width=width,height=height,bg=bg) else cairo_pdf
  ggsave(file.path(out,paste0(stem,'.pdf')),p,width=9.5,height=h,device=device,bg='white')
  ggsave(file.path(out,paste0(stem,'.svg')),p,width=9.5,height=h,device=svglite::svglite,bg='white')
  ggsave(file.path(out,paste0(stem,'.png')),p,width=9.5,height=h,dpi=250,bg='white')
}
save(make_plot(),'celltype_all26',8.4)
save(make_plot(TRUE),'celltype_tested10',4.7)
# Equivalent CVD checks with colorspace; colorblindr is unavailable in this runtime.
for(v in c('deutan','protan','tritan'))ggsave(file.path(out,paste0('CVD_',v,'.png')),make_plot(TRUE,v),width=9.5,height=4.7,dpi=150,bg='white')
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('Validated:',nrow(nominated),'nominated DOFs;',length(tissues),'cell identities;',sum(d$mixed),'opposing-zone pairs;',sum(x$de_called),'DE stars.\n')
