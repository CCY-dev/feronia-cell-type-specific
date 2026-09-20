#!/usr/bin/env Rscript
# Prototype; not yet carefully vetted.
# Programmatic concept: source dependency + exact DOF-family membership.
suppressPackageStartupMessages({library(ggplot2);library(colorspace)})
pdfFonts(Arial=pdfFonts()$Helvetica); options(device=function(...) grDevices::pdf(file=NULL,...))
arg<-commandArgs(trailingOnly=FALSE)
script<-normalizePath(sub('^--file=','',arg[grepl('^--file=',arg)][1]))
out<-dirname(script);base<-normalizePath(file.path(out,'..'))
m<-read.csv(file.path(base,'fig6_dof_membership.csv'),stringsAsFactors=FALSE)
for(n in c('nominated_any','called_DE_any','tested_LUC')) m[[n]]<-tolower(m[[n]])=='true'
v<-read.csv(file.path(base,'fig6_full_intersections.csv'))
stopifnot(nrow(m)==37,!anyDuplicated(m$gene),sum(m$nominated_any)==26,
 sum(m$called_DE_any)==12,sum(m$tested_LUC)==10,
 !any(m$called_DE_any & !m$nominated_any),!any(m$tested_LUC & !m$nominated_any),
 sum(v$genes[v$GRN_nominated==1])==512,sum(v$genes[v$called_DEG==1])==5466,
 sum(v$genes[v$GRN_nominated==1&v$called_DEG==1])==184)
m$group<-ifelse(m$nominated_any,ifelse(m$called_DE_any,'both','grn_only'),'neither')
m<-m[order(match(m$group,c('both','grn_only','neither')),-m$tested_LUC,m$gene),]
m$within_group<-ave(seq_len(nrow(m)),m$group,FUN=seq_along)
m$center<-c(both=2.5,grn_only=7,neither=11.5)[m$group]
m$x<-m$center+((m$within_group-1)%%7-3)*.43
m$y<-2.8-floor((m$within_group-1)/7)*.55
stopifnot(sum(m$tested_LUC[m$group=='both'])==8,sum(m$tested_LUC[m$group=='grn_only'])==2)
write.csv(m,file.path(out,'dot_membership.csv'),row.names=FALSE)
txt<-function(x,y,label,size=3.2,fontface='plain',colour='#18212B',hjust=.5) annotate('text',x=x,y=y,label=label,size=size,fontface=fontface,colour=colour,hjust=hjust,family='Arial',lineheight=1.1)
box<-function(xmin,xmax,ymin,ymax,fill,colour=NA) annotate('rect',xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill,colour=colour,size=.4)
arr<-function(x,y,xend,yend) annotate('segment',x=x,y=y,xend=xend,yend=yend,colour='#7E8892',size=.5,arrow=grid::arrow(length=grid::unit(1.8,'mm'),type='closed'))
plot_story<-function(cvd=NULL) {
 pal<-c(both='#CC79A7',grn_only='#0072B2',neither='#AFB7BE')
 if(!is.null(cvd))pal[]<-get(cvd,asNamespace('colorspace'))(pal)
 ggplot()+
  box(.15,13.85,0,11.8,'white')+
  txt(.4,11.5,'From single-cell profiles to DOF candidates',5,'bold',hjust=0)+
  box(2.7,8.7,9.65,10.75,'#EEF1F4')+
  txt(5.7,10.31,'Basal WT and fer-4 scRNA-seq',3.9,'bold')+
  txt(5.7,9.91,'Cell-type-resolved expression profiles',2.9)+
  box(10.1,13.5,9.8,10.6,'#F2F4F6')+
  txt(11.8,10.2,'Regulatory prior',3.4,'bold')+
  arr(4.3,9.58,3.0,8.6)+arr(7.3,9.58,10.1,8.6)+arr(11.8,9.7,11.8,8.6)+
  box(.4,5.9,7.1,8.55,'#F7EFF4')+
  txt(3.15,8.08,'5,466 DEGs',4.4,'bold')+
  txt(3.15,7.55,'Basal fer-4 versus WT\nCalled in at least one cell state',2.9)+
  box(8.1,13.6,7.1,8.55,'#EBF4F9')+
  txt(10.85,8.08,'512 GRN-nominated TFs',4.1,'bold')+
  txt(10.85,7.55,'Selected by relative\nout-degree rank changes',2.9)+
  arr(6.0,7.55,7.98,7.55)+txt(7.0,8.04,'DEG-informed\ntarget filtering',2.5)+
  arr(3.15,7.02,5.05,6.56)+arr(10.85,7.02,8.95,6.56)+
  box(4.8,9.2,5.95,6.68,'#F0F2F5')+
  txt(7,6.32,'Focus on the 37 DOFs',3.8,'bold')+
  txt(4.75,5.56,'26 GRN-nominated DOFs',3.4,'bold',colour=pal['grn_only'])+
  annotate('segment',x=.8,xend=8.7,y=5.27,yend=5.27,colour=pal['grn_only'],size=.65)+
  annotate('segment',x=.8,xend=.8,y=5.1,yend=5.27,colour=pal['grn_only'],size=.65)+
  annotate('segment',x=8.7,xend=8.7,y=5.1,yend=5.27,colour=pal['grn_only'],size=.65)+
  box(.4,4.6,1.45,5.0,'#FAF4F8')+box(4.9,9.1,1.45,5.0,'#F0F6FA')+box(9.4,13.6,1.45,5.0,'#F3F4F5')+
  txt(2.5,4.47,'12',7.2,'bold',pal['both'])+
  txt(7,4.47,'14',7.2,'bold',pal['grn_only'])+
  txt(11.5,4.47,'11',7.2,'bold','#606B75')+
  txt(2.5,3.81,'GRN-nominated + DEG',3.2,'bold')+
  txt(7,3.81,'GRN-nominated, no DEG call',3.05,'bold')+
  txt(11.5,3.81,'Not GRN-nominated',3.2,'bold')+
  txt(11.5,3.38,'None has a DEG call',2.7,'plain','#606B75')+
  geom_point(data=m,aes(x,y,colour=group),shape=16,size=4.3)+
  geom_point(data=m[m$tested_LUC,],aes(x,y),shape=21,fill=NA,colour='#17212B',size=4.6,stroke=1.0)+
  scale_colour_manual(values=pal,guide='none')+
  txt(2.5,1.79,'8 assayed',3.25,'bold')+txt(7,1.79,'2 assayed',3.25,'bold')+txt(11.5,1.79,'0 assayed',3.25,'plain','#606B75')+
  txt(.4,.95,'One dot = one DOF. Black outline = assayed with the pYUC5-LUC reporter (10 total).',2.95,hjust=0)+
  txt(.4,.5,'Across all genes: 184 GRN-nominated TFs are DEGs = 12 DOFs + 172 other TFs.',2.95,hjust=0)+
  coord_cartesian(xlim=c(0,14),ylim=c(.1,11.8),expand=FALSE,clip='off')+
  theme_void(base_family='Arial')+theme(plot.margin=margin(6,6,6,6))
}
save_plot<-function(p,stem) {
 dev<-function(filename,width,height,bg,...)grDevices::quartz(type='pdf',file=filename,width=width,height=height,bg=bg)
 ggsave(file.path(out,paste0(stem,'.pdf')),p,width=8.2,height=6.8,device=dev,bg='white')
 ggsave(file.path(out,paste0(stem,'.svg')),p,width=8.2,height=6.8,device=svglite::svglite,bg='white')
 ggsave(file.path(out,paste0(stem,'.png')),p,width=8.2,height=6.8,dpi=250,bg='white')
}
save_plot(plot_story(),'singlecell_dof_evidence_map')
for(v in c('deutan','protan','tritan'))ggsave(file.path(out,paste0('CVD_',v,'.png')),plot_story(v),width=8.2,height=6.8,dpi=150,bg='white')
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('Verified: 37 DOF dots =12+14+11; 10 outlined =8+2; global GRN/DEG overlap184 =12+172.\n')
