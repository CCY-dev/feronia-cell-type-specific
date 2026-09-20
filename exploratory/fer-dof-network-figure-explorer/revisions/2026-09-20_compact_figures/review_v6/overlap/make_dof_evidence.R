#!/usr/bin/env Rscript
# Prototype; not yet carefully vetted.
# Review v6: stack starts derived from counts; input-identity assertions retained.
# Larger v5 labels, source counts and rendered geometry retained.
# Conditional overlap within the DOF family, without unequal-universe area scaling.
suppressPackageStartupMessages({library(ggplot2);library(cowplot);library(colorspace)})
options(device=function(...) grDevices::pdf(file=NULL,...))
args <- commandArgs(FALSE)
out <- dirname(normalizePath(sub('^--file=','',args[grepl('^--file=',args)][1])))
src <- normalizePath(file.path(out,'../../../2026-09-20_fig6_brainstorm/overlap_story/dot_membership.csv'))
setsrc <- normalizePath(file.path(out,'../../../2026-09-20_fig6_brainstorm/fig6_full_intersections.csv'))
d <- read.csv(src,check.names=FALSE)
sets <- read.csv(setsrc)
stopifnot(nrow(d)==37L,!anyDuplicated(d$gene),sum(d$nominated_any)==26L,
  sum(d$called_DE_any)==12L,sum(d$nominated_any&d$called_DE_any)==12L,
  sum(!d$nominated_any&d$called_DE_any)==0L,sum(d$tested_LUC)==10L)
counts <- expand.grid(nominated=c(TRUE,FALSE),DE=c(TRUE,FALSE))
counts$genes <- mapply(function(n,de)sum(d$nominated_any==n & d$called_DE_any==de),counts$nominated,counts$DE)
counts <- counts[order(!counts$nominated,!counts$DE),]
# Derive stack starts within each nomination group; assertions below freeze input identity.
counts$x0 <- ave(counts$genes, counts$nominated, FUN=function(x)c(0,head(cumsum(x),-1)))
counts$x1 <- counts$x0+counts$genes
counts$y <- ifelse(counts$nominated,27,16)
stopifnot(identical(as.integer(counts$genes),c(12L,14L,0L,11L)))
# Cross-check against the independently retained complete intersection table.
for(i in seq_len(nrow(counts))) stopifnot(counts$genes[i] ==
  sets$genes[sets$DOF==1 & sets$GRN_nominated==as.integer(counts$nominated[i]) & sets$called_DEG==as.integer(counts$DE[i])])
write.csv(d[,c('gene','symbol','nominated_any','called_DE_any','tested_LUC')],file.path(out,'plotted_gene_membership.csv'),row.names=FALSE)
write.csv(counts,file.path(out,'conditional_counts.csv'),row.names=FALSE)
pal <- c(DE='#B65F91',notDE='#D6DCE0',text='#28343C',axis='#64717A')
pt <- function(x)1.15*x/(72.27/25.4)
bx <- function(x)30+52*x/30
make_plot <- function(cvd=NULL) {
  p <- if(is.null(cvd))pal else setNames(get(cvd,asNamespace('colorspace'))(pal),names(pal))
  g <- counts[counts$genes>0,]
  g$fill <- ifelse(g$DE,p['DE'],p['notDE'])
  g$text <- ifelse(g$DE,'white',p['text'])
  ggplot()+
    annotate('text',x=1,y=43,label='37 DOF family members',hjust=0,family='Arial',fontface='bold',size=pt(8.2),colour=p['text'])+
    annotate('rect',xmin=c(1,27),xmax=c(3.5,29.5),ymin=35,ymax=37.5,fill=c(p['DE'],p['notDE']),colour=NA)+
    annotate('text',x=c(5,31),y=36.2,label=c('DEG call','No DEG call'),hjust=0,family='Arial',size=pt(7.3),colour=p['text'])+
    annotate('text',x=1,y=c(27,16),label=c('GRN\nnominated','Not GRN\nnominated'),hjust=0,family='Arial',size=pt(7.4),lineheight=.95,colour=p['text'])+
    geom_rect(data=g,aes(xmin=bx(x0),xmax=bx(x1),ymin=y-3.4,ymax=y+3.4,fill=I(fill)),colour=NA)+
    geom_text(data=g,aes(x=bx((x0+x1)/2),y=y,label=genes,colour=I(text)),family='Arial',size=pt(8),fontface='bold')+
    annotate('text',x=bx(11)+2,y=16,label='0 DEGs',hjust=0,family='Arial',size=pt(7.2),colour=p['DE'])+
    annotate('segment',x=bx(0),xend=bx(30),y=9,yend=9,size=.25,colour=p['axis'])+
    annotate('segment',x=bx(c(0,10,20,30)),xend=bx(c(0,10,20,30)),y=9,yend=8.2,size=.25,colour=p['axis'])+
    annotate('text',x=bx(c(0,10,20,30)),y=6,label=c('0','10','20','30'),family='Arial',size=pt(7.2),colour=p['axis'])+
    annotate('text',x=bx(15),y=1.8,label='DOF genes',family='Arial',size=pt(7.4),colour=p['text'])+
    coord_fixed(xlim=c(0,85),ylim=c(0,46),expand=FALSE,clip='off')+
    cowplot::theme_map(font_family='Arial',font_size=7.4)+theme(plot.margin=margin(0,0,0,0))
}
pdfdev <- if(capabilities('aqua'))function(filename,width,height,bg,...)grDevices::quartz(type='pdf',file=filename,width=width,height=height,bg=bg)else grDevices::cairo_pdf
p <- make_plot()
ggsave(file.path(out,'dof_evidence_bars.pdf'),p,width=85/25.4,height=46/25.4,device=pdfdev,bg='white')
ggsave(file.path(out,'dof_evidence_bars.svg'),p,width=85/25.4,height=46/25.4,device=svglite::svglite,bg='white')
ggsave(file.path(out,'dof_evidence_bars.png'),p,width=85/25.4,height=46/25.4,dpi=400,bg='white')
for(cv in c('deutan','protan','tritan'))ggsave(file.path(out,paste0('cvd_',cv,'.png')),make_plot(cv),width=85/25.4,height=46/25.4,dpi=300,bg='white')
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
writeLines(system2('shasum',c('-a','256',shQuote(src),shQuote(setsrc)),stdout=TRUE),file.path(out,'source_sha256.txt'))
writeLines(c('12 nominated + DEG; 14 nominated without DEG; 0 DEG without nomination; 11 neither.',
  '37 unique gene identifiers; 26 nominated; 12 DEG-called; no inference or selection rerun.',
  'Bars share a linear axis starting at 0. Minimum type 8.28 pt; 85 x 46 mm.'),file.path(out,'validation.txt'))
