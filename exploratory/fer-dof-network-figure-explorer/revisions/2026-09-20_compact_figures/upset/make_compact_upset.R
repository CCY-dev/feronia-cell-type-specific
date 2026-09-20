#!/usr/bin/env Rscript
# Prototype; not yet carefully vetted.
# 85 x 55 mm transposed UpSet. Existing audited counts; no new analysis.
suppressPackageStartupMessages({library(ggplot2);library(cowplot);library(colorspace)})
options(device=function(...) grDevices::pdf(file=NULL,...))
a<-commandArgs(FALSE)
out<-dirname(normalizePath(sub('^--file=','',a[grepl('^--file=',a)][1])))
src<-normalizePath(file.path(out,'../../2026-09-20_fig6_brainstorm/fig6_full_intersections.csv'))
all<-read.csv(src,check.names=FALSE)
all$key<-with(all,paste0(DOF,GRN_nominated,called_DEG))
stopifnot(!anyDuplicated(all$key),nrow(all)==7L,sum(all$genes)==5805L,
          sum(all$genes[all$DOF==1])==37L,
          sum(all$genes[all$GRN_nominated==1])==512L,
          sum(all$genes[all$called_DEG==1])==5466L)
d<-all[match(c('111','110','101','011'),all$key),]
stopifnot(identical(as.integer(d$genes),c(12L,14L,0L,172L)))
d$y<-c(37.5,29.5,21.5,13.5)
d$focus<-d$key=='111'
write.csv(d,file.path(out,'plot_input_overlaps.csv'),row.names=FALSE)
write.csv(all,file.path(out,'all_exclusive_intersections.csv'),row.names=FALSE)
sets<-data.frame(label=c('DOFs','GRN-TFs','scDEGs'),total=c(37L,512L,5466L),x=c(9,22,35))
write.csv(sets,file.path(out,'set_totals.csv'),row.names=FALSE)
base_cols<-c(focus='#BA6B98',other='#4E5961',inactive='#DDE2E5',text='#252B30',axis='#777F85',stripe='#F6F7F8')
pt<-function(x)x/(72.27/25.4)
make_plot<-function(cvd=NULL) {
  pal<-if(is.null(cvd))base_cols else setNames(get(cvd,asNamespace('colorspace'))(base_cols),names(base_cols))
  z<-d;z$color<-ifelse(z$focus,pal['focus'],pal['other'])
  dots<-do.call(rbind,lapply(seq_len(nrow(z)),function(i)
    data.frame(x=sets$x,y=z$y[i],member=as.integer(z[i,1:3]),color=z$color[i])))
  lines<-do.call(rbind,lapply(seq_len(nrow(z)),function(i){
    xx<-sets$x[as.integer(z[i,1:3])==1]
    data.frame(x=min(xx),xend=max(xx),y=z$y[i],color=z$color[i])
  }))
  # Horizontal bars: 200 genes occupy 38 mm, with a zero origin at x=43.
  bx<-function(v)43+38*v/200
  ticks<-data.frame(x=bx(c(0,50,100,150,200)),lab=c('0','50','100','150','200'))
  ggplot()+
    annotate('rect',xmin=3,xmax=83,ymin=c(26.5,10.5),ymax=c(32.5,16.5),fill=pal['stripe'],colour=NA)+
    geom_text(data=sets,aes(x=x,y=49,label=label),family='Arial',size=pt(7.7),colour=pal['text'])+
    geom_text(data=sets,aes(x=x,y=44.9,label=format(total,big.mark=',',trim=TRUE)),family='Arial',size=pt(7.2),colour=pal['axis'])+
    annotate('text',x=43,y=49,label='Exclusive overlap',hjust=0,family='Arial',size=pt(7.7),colour=pal['text'])+
    annotate('text',x=43,y=44.9,label='(genes)',hjust=0,family='Arial',size=pt(7.2),colour=pal['axis'])+
    geom_segment(data=lines,aes(x=x,xend=xend,y=y,yend=y,colour=I(color)),size=.45)+
    geom_point(data=dots,aes(x=x,y=y),colour=pal['inactive'],size=2.55)+
    geom_point(data=dots[dots$member==1,],aes(x=x,y=y,colour=I(color)),size=2.9)+
    geom_rect(data=z[z$genes>0,],aes(xmin=43,xmax=bx(genes),ymin=y-1.9,ymax=y+1.9,fill=I(color)),colour=NA)+
    geom_text(data=z,aes(x=bx(genes)+1.25,y=y,label=genes,colour=I(color)),hjust=0,
              family='Arial',size=pt(7.8),fontface=ifelse(z$focus,'bold','plain'))+
    annotate('segment',x=43,xend=81,y=8.2,yend=8.2,colour=pal['axis'],size=.28)+
    geom_segment(data=ticks,aes(x=x,xend=x,y=8.2,yend=7.4),colour=pal['axis'],size=.28)+
    geom_text(data=ticks,aes(x=x,y=4.7,label=lab),family='Arial',size=pt(7.2),colour=pal['axis'])+
    coord_fixed(xlim=c(0,85),ylim=c(0,55),expand=FALSE,clip='off')+
    cowplot::theme_map(font_family='Arial',font_size=7.7)+theme(plot.margin=margin(0,0,0,0))
}
w<-85/25.4;h<-55/25.4
p<-make_plot()
pdf_device<-if(capabilities('aqua')) function(filename,width,height,bg,...) grDevices::quartz(type='pdf',file=filename,width=width,height=height,bg=bg) else grDevices::cairo_pdf
ggsave(file.path(out,'compact_upset.pdf'),p,width=w,height=h,device=pdf_device,bg='white')
ggsave(file.path(out,'compact_upset.svg'),p,width=w,height=h,device=svglite::svglite,bg='white')
ggsave(file.path(out,'compact_upset.png'),p,width=w,height=h,dpi=300,bg='white')
for(cv in c('deutan','protan','tritan'))
  ggsave(file.path(out,paste0('cvd_',cv,'.png')),make_plot(cv),width=w,height=h,dpi=300,bg='white')
status<-system2('python',c(shQuote(file.path(out,'group_svg.py')),shQuote(file.path(out,'compact_upset.svg'))))
stopifnot(status==0L)
writeLines(system2('shasum',c('-a','256',shQuote(src)),stdout=TRUE),file.path(out,'source_sha256.txt'))
writeLines(c('Dimensions:85 x55 mm','Minimum text:7.2pt; main headers7.7pt; counts7.8pt.',
  'Counts verified:12,14,0,172. Set totals37,512,5466. Union5805.',
  'Linear bar axis0-200; no axis breaks or pseudocounts.'),file.path(out,'validation.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
