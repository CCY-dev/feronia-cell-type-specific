#!/usr/bin/env Rscript
# PROTOTYPE: not yet carefully vetted. Cosmetic derivative of the as-exported union.
# No changes to source pairs, target-family selection, support classes or source analyses.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(colorspace)})
options(device=function(...) grDevices::pdf(file=NULL,...))
pdfFonts(Arial=pdfFonts()$Helvetica)
args<-commandArgs(FALSE)
out<-dirname(normalizePath(sub('^--file=','',args[grepl('^--file=',args)][1])))
src<-normalizePath(file.path(out,'../../review_v2/all_contexts'))
input<-file.path(src,'all_contexts_as_exported_union_edges.csv')
file.copy(input,file.path(out,'source_union_edges.csv'),overwrite=TRUE)
u<-fread(input)
stopifnot(nrow(u)==180,uniqueN(u[,.(source,target)])==180,
          uniqueN(u$source)==10,uniqueN(u$target)==37,
          sum(u$support_class=='supported_in_current_seven')==147,
          sum(u$support_class=='additional_excluded_contexts_only')==33)
tf_order<-c('DOF4.1','CDF4','DOF4.6','OBP2','DAG1','OBP3','COG1','DOF2.2','OBP1','DOF3.2')
t<-unique(u[,.(target,target_name)])
t[,family:=fifelse(grepl('^(YUC|TAA|TAR)',target_name),'YUC/TAA',
                 fifelse(grepl('^GH3',target_name),'GH3',
                 fifelse(grepl('^SAUR',target_name),'SAUR','Aux/IAA')))]
t[,num:=as.numeric(gsub('[^0-9.]','',target_name))]
t[,family_order:=match(family,c('YUC/TAA','GH3','Aux/IAA','SAUR'))]
setorder(t,family_order,num,target_name)
s<-unique(u[,.(source,source_symbol)]);s<-s[match(tf_order,source_symbol)]
# Coordinates are millimetres of display space, not biological measurements.
s[,`:=`(x=seq(14,121,length.out=.N),y=82)]
theta<-seq(pi,2*pi,length.out=nrow(t))
t[,`:=`(x=67.5+49*cos(theta),y=59+34*sin(theta),
        lx=67.5+52*cos(theta),ly=59+37*sin(theta),
        angle=ifelse(theta<=1.5*pi,theta*180/pi-180,theta*180/pi-360),
        hjust=ifelse(theta<=1.5*pi,1,0))]
e<-merge(merge(u,s[,.(source,x0=x,y0=y)],by='source'),
         t[,.(target,x1=x,y1=y)],by='target')
support_levels<-c('supported_in_current_seven','additional_excluded_contexts_only')
e[,support_class:=factor(support_class,levels=support_levels)]
setorder(e,support_class)
stopifnot(nrow(e)==180,all(e$y0>e$y1))
pal<-c(DOF='#253B49',`YUC/TAA`='#E69F00',GH3='#CC79A7',`Aux/IAA`='#70B7D9',SAUR='#009E73')
pt<-72.27/25.4
family_key<-data.table(family=c('YUC/TAA','GH3','Aux/IAA','SAUR'),x=c(20,50,77,108),y=10)
support_key<-data.table(support_class=factor(support_levels,levels=support_levels),
                       x=c(6,73),xend=c(15,82),y=3.3,
                       label=c('147 pairs: seven fits','33 pairs: other fits only'))
plot_net<-function(cvd=NULL){
  cp<-pal
  if(!is.null(cvd))cp<-setNames(get(cvd,asNamespace('colorspace'))(cp),names(cp))
  ggplot()+
    geom_segment(data=e,aes(x0,y0,xend=x1,yend=y1,linetype=support_class,colour=support_class),
                 size=.15,alpha=.50)+
    geom_point(data=t,aes(x,y,fill=family),shape=21,colour='white',size=2.7,stroke=.2)+
    geom_text(data=t,aes(lx,ly,label=target_name,hjust=hjust,angle=angle),
              size=8.5/pt,fontface='italic',family='Arial',colour='#15232B')+
    geom_point(data=s,aes(x,y),shape=21,fill=cp['DOF'],colour='white',size=3.8,stroke=.25)+
    geom_text(data=s,aes(x,y=y+4.4,label=source_symbol),
              size=9/pt,fontface='italic',family='Arial',colour='#15232B')+
    geom_point(data=family_key,aes(x,y,fill=family),shape=21,colour='white',size=2.7,stroke=.2)+
    geom_text(data=family_key,aes(x=x+3,y,label=family),hjust=0,
              size=8.5/pt,family='Arial',colour='#15232B')+
    geom_segment(data=support_key,aes(x,y,xend=xend,yend=y,linetype=support_class,colour=support_class),
                 size=.3,alpha=1)+
    geom_text(data=support_key,aes(x=xend+2,y,label=label),hjust=0,
              size=8.5/pt,family='Arial',colour='#15232B')+
    scale_fill_manual(values=cp,guide='none')+
    scale_linetype_manual(values=c(supported_in_current_seven='solid',additional_excluded_contexts_only='22'),guide='none')+
    scale_colour_manual(values=c(supported_in_current_seven='#B6BFC5',additional_excluded_contexts_only='#87939C'),guide='none')+
    coord_cartesian(xlim=c(0,135),ylim=c(0,90),expand=FALSE,clip='off')+
    theme_void(base_family='Arial',base_size=8.5)+theme(plot.margin=margin(0,0,0,0))
}
p<-plot_net();name<-'all_contexts_panel_d'
ggsave(file.path(out,paste0(name,'.png')),p,width=135/25.4,height=90/25.4,dpi=350,bg='white')
ggsave(file.path(out,paste0(name,'.svg')),p,width=135/25.4,height=90/25.4,device=svglite::svglite,bg='white')
quartz(type='pdf',file=file.path(out,paste0(name,'.pdf')),width=135/25.4,height=90/25.4,family='Arial',bg='white')
print(p);dev.off()
for(cv in c('deutan','protan','tritan'))
  ggsave(file.path(out,paste0('cvd_',cv,'.png')),plot_net(cv),width=135/25.4,height=90/25.4,dpi=180,bg='white')
if(requireNamespace('colorblindr',quietly=TRUE))colorblindr::cvd_grid(p)
fwrite(s,file.path(out,'plotted_tf_nodes.csv'))
fwrite(t,file.path(out,'plotted_target_nodes.csv'))
fwrite(e,file.path(out,'plotted_edges.csv'))
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
cat('Panel d: 135 x 90 mm; 10 DOFs; 37 targets; 180 unchanged pairs.\n')
