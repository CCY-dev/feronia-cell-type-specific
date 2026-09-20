#!/usr/bin/env Rscript
# Prototype; not yet carefully vetted.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(colorspace)})
options(device=function(...) grDevices::pdf(file=NULL,...));pdfFonts(Arial=pdfFonts()$Helvetica)
a<-commandArgs(FALSE);out<-dirname(normalizePath(sub('^--file=','',a[grepl('^--file=',a)][1])))
u<-fread(file.path(out,'all_contexts_as_exported_union_edges.csv'))
stopifnot(nrow(u)==180,uniqueN(u$source)==10,uniqueN(u$target)==37)
tf_order<-c('DOF4.1','CDF4','DOF4.6','OBP2','DAG1','OBP3','COG1','DOF2.2','OBP1','DOF3.2')
t<-unique(u[,.(target,target_name)])
t[,family:=fifelse(grepl('^(YUC|TAA|TAR)',target_name),'YUC/TAA',fifelse(grepl('^GH3',target_name),'GH3',fifelse(grepl('^SAUR',target_name),'SAUR','Aux/IAA')))]
t[,num:=as.numeric(gsub('[^0-9.]','',target_name))]
t[,family_order:=match(family,c('YUC/TAA','GH3','Aux/IAA','SAUR'))]
setorder(t,family_order,num,target_name)
s<-unique(u[,.(source,source_symbol)]);s<-s[match(tf_order,source_symbol)]
s[,`:=`(x=seq(-59,59,length.out=.N),y=50)]
theta<-seq(pi,2*pi,length.out=nrow(t));ry<-48
# Coordinates are display geometry, not scores or effect sizes.
t[,`:=`(x=56*cos(theta),y=19+ry*sin(theta),lx=59*cos(theta),ly=19+(ry+2)*sin(theta),
 angle=ifelse(theta<=1.5*pi,theta*180/pi-180,theta*180/pi-360),hjust=ifelse(theta<=1.5*pi,1,0))]
e<-merge(merge(u,s[,.(source,x0=x,y0=y)],by='source'),t[,.(target,x1=x,y1=y)],by='target')
support_levels<-c('supported_in_current_seven','additional_excluded_contexts_only')
e[,support_class:=factor(support_class,levels=support_levels)]
setorder(e,support_class)
pal<-c(DOF='#253B49',`YUC/TAA`='#E69F00',GH3='#CC79A7',`Aux/IAA`='#70B7D9',SAUR='#009E73')
plot_net<-function(cvd=NULL){
 cp<-pal;if(!is.null(cvd))cp<-setNames(get(cvd,asNamespace('colorspace'))(cp),names(cp))
 ggplot()+
 geom_segment(data=e,aes(x0,y0,xend=x1,yend=y1,linetype=support_class,colour=support_class),size=.19,alpha=.75)+
 geom_point(data=t,aes(x,y,fill=family),shape=21,colour='white',size=2.5,stroke=.25)+
 geom_text(data=t,aes(lx,ly,label=target_name,hjust=hjust,angle=angle),size=2.55,fontface='italic',family='Arial',colour='#15232B')+
 geom_point(data=s,aes(x,y),size=3.0,shape=21,fill=cp['DOF'],colour='white',stroke=.3)+
 geom_text(data=s,aes(x,y=y+5,label=source_symbol),size=2.6,fontface='italic',family='Arial',colour='#15232B')+
 scale_fill_manual(values=cp,breaks=c('YUC/TAA','GH3','Aux/IAA','SAUR'),drop=FALSE,name=NULL)+
 scale_linetype_manual(values=c('supported_in_current_seven'='solid','additional_excluded_contexts_only'='22'),
  labels=c('Supported in seven included fits','Additional pairs from excluded fits'),name=NULL)+
 scale_colour_manual(values=c('supported_in_current_seven'='#ADB5BA','additional_excluded_contexts_only'='#69747C'),guide='none')+
 guides(fill=guide_legend(override.aes=list(size=2.6),nrow=1,order=1),
  linetype=guide_legend(override.aes=list(colour=c('#ADB5BA','#69747C'),size=.3,alpha=1),nrow=1,order=2))+
 coord_cartesian(xlim=c(-71,71),ylim=c(19-ry-17,59),expand=FALSE,clip='off')+
 labs(title='DOF links across all exported cellular contexts',subtitle='Basal WT and fer-4')+
 theme_void(base_family='Arial',base_size=7.4)+
 theme(plot.title=element_text(size=10,face='bold',margin=margin(b=2)),plot.subtitle=element_text(size=7.5,margin=margin(b=5)),
 plot.title.position='plot',plot.margin=margin(5,6,2,6),legend.position='bottom',legend.box='vertical',
 legend.text=element_text(size=7.4),legend.key.width=unit(7,'mm'),legend.key.height=unit(3,'mm'),
 legend.spacing.y=unit(0,'mm'),legend.spacing.x=unit(2,'mm'),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0))
}
p<-plot_net();name<-'all_contexts_overview_compact'
ggsave(file.path(out,paste0(name,'.png')),p,width=170/25.4,height=100/25.4,dpi=350,bg='white')
ggsave(file.path(out,paste0(name,'.svg')),p,width=170/25.4,height=100/25.4,device=svglite::svglite,bg='white')
quartz(type='pdf',file=file.path(out,paste0(name,'.pdf')),width=170/25.4,height=100/25.4);print(p);dev.off()
for(cv in c('deutan','protan','tritan'))ggsave(file.path(out,paste0('cvd_',cv,'.png')),plot_net(cv),width=170/25.4,height=100/25.4,dpi=180,bg='white')
fwrite(s,file.path(out,paste0(name,'_tf_nodes.csv')));fwrite(t,file.path(out,paste0(name,'_target_nodes.csv')))
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
