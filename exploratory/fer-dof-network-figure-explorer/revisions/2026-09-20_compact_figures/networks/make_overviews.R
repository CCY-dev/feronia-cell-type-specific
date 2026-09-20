#!/usr/bin/env Rscript
# Prototype; not yet carefully vetted.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(cowplot);library(colorspace)})
options(device=function(...) grDevices::pdf(file=NULL,...));pdfFonts(Arial=pdfFonts()$Helvetica)
a<-commandArgs(FALSE);out<-dirname(normalizePath(sub('^--file=','',a[grepl('^--file=',a)][1])))
base<-file.path(dirname(dirname(out)),'2026-09-20_dof_subnetworks')
inp<-file.path(base,'inventory/tested10_retained_edges_seven_unaffected_consistent_tissues.csv')
ann<-file.path(base,'annotations/clean_programme_target_edges.tsv')
e<-fread(inp);s<-unique(fread(ann)[,.(target,symbol)])
stopifnot(!anyDuplicated(s$target),!any(e$tissue%in%c('endodermis','phloem','xylem')))
e[,target_export_name:=target_name]
e[target%in%s$target,target_name:=s$symbol[match(target,s$target)]]
e<-e[grepl('^(IAA[0-9]|SAUR[0-9]|GH3[.]|TAA[0-9]|TAR[0-9]|YUC[0-9])',target_name)]
stopifnot(!anyDuplicated(e[,.(source,target,cluster)]))
contexts<-c('Elongation_Phloem Pole Pericycle','Elongation_Xylem Pole Pericycle','Elongation_Procambium')
tf_order<-c('DOF4.1','CDF4','DOF4.6','OBP2','DAG1','OBP3','COG1','DOF2.2','OBP1','DOF3.2')
pal<-c(DOF='#253B49',`YUC/TAA`='#E69F00',GH3='#CC79A7',`Aux/IAA`='#70B7D9',SAUR='#009E73')

make_data<-function(z,prefix,expected) {
 u<-z[,.(n_contexts=uniqueN(cell_type),contexts=paste(sort(unique(cell_type)),collapse='; '),
   contexts_WT=paste(sort(unique(cell_type[genotype=='WT'])),collapse='; '),
   contexts_fer4=paste(sort(unique(cell_type[genotype=='fer4'])),collapse='; '),
   n_context_genotype_links=.N),by=.(source,source_symbol,target,target_name)]
 stopifnot(uniqueN(u$source)==10,uniqueN(u$target)==expected[1],nrow(u)==expected[2],uniqueN(z$cell_type)==expected[3])
 t<-unique(u[,.(target,target_name)])
 t[,family:=fifelse(grepl('^(YUC|TAA|TAR)',target_name),'YUC/TAA',fifelse(grepl('^GH3',target_name),'GH3',fifelse(grepl('^SAUR',target_name),'SAUR','Aux/IAA')))]
 t[,num:=as.numeric(gsub('[^0-9.]','',target_name))]
 t[,family_order:=match(family,c('YUC/TAA','GH3','Aux/IAA','SAUR'))]
 setorder(t,family_order,num,target_name)
 ss<-unique(u[,.(source,source_symbol)]);ss<-ss[match(tf_order,source_symbol)]
 ss[,`:=`(x=seq(-59,59,length.out=.N),y=50)]
 # A U-shaped target arc. Geometry is fixed and encodes no biological quantity.
 theta<-seq(pi,2*pi,length.out=nrow(t)); ry<-if(nrow(t)>25)48 else 35
 t[,`:=`(x=56*cos(theta),y=19+ry*sin(theta))]
 t[,`:=`(lx=60*cos(theta),ly=19+(ry+4)*sin(theta))]
 t[,hjust:=fifelse(x>43,0,fifelse(x< -43,1,.5))]
 t[,angle:=0]
 if(nrow(t)>25) {
  # Radial labels avoid collisions at the denser lower arc at print size.
  t[,`:=`(lx=59*cos(theta),ly=19+(ry+2)*sin(theta),
    angle=ifelse(theta<=1.5*pi,theta*180/pi-180,theta*180/pi-360),
    hjust=ifelse(theta<=1.5*pi,1,0))]
 }
 ed<-merge(merge(u,ss[,.(source,x0=x,y0=y)],by='source'),t[,.(target,x1=x,y1=y)],by='target')
 fwrite(z,file.path(out,paste0(prefix,'_source_context_edges.csv')))
 fwrite(u,file.path(out,paste0(prefix,'_union_edges.csv')))
 fwrite(ss,file.path(out,paste0(prefix,'_tf_nodes.csv')));fwrite(t,file.path(out,paste0(prefix,'_target_nodes.csv')))
 fwrite(z[,.(retained_links=.N,unique_pairs=uniqueN(paste(source,target))),by=.(tissue,cell_type,genotype)],file.path(out,paste0(prefix,'_context_counts.csv')))
 list(s=ss,t=t,e=ed,ry=ry,prefix=prefix,counts=c(TFs=10,targets=nrow(t),pairs=nrow(u),contexts=uniqueN(z$cell_type),source_rows=nrow(z)))
}
v<-make_data(e[cell_type%in%contexts],'vascular_overview_compact',c(23,81,3))
b<-make_data(e,'all_tissues_overview_compact',c(34,147,27))
plot_net<-function(d,title,cvd=NULL) {
 cp<-pal;if(!is.null(cvd))cp<-setNames(get(cvd,asNamespace('colorspace'))(cp),names(cp))
 ggplot()+
  geom_segment(data=d$e,aes(x0,y0,xend=x1,yend=y1),colour='#AAB2B8',size=.19,alpha=.72)+
  geom_point(data=d$t,aes(x,y,fill=family),shape=21,colour='white',size=2.5,stroke=.25)+
  geom_text(data=d$t,aes(lx,ly,label=target_name,hjust=hjust,angle=angle),size=2.55,fontface='italic',family='Arial',colour='#15232B')+
  geom_point(data=d$s,aes(x,y),size=3.0,shape=21,fill=cp['DOF'],colour='white',stroke=.3)+
  geom_text(data=d$s,aes(x,y=y+5,label=source_symbol),size=2.6,fontface='italic',family='Arial',colour='#15232B')+
  scale_fill_manual(values=cp,breaks=c('YUC/TAA','GH3','Aux/IAA','SAUR'),drop=FALSE,name=NULL)+
  guides(fill=guide_legend(override.aes=list(size=2.6),nrow=1))+
  coord_cartesian(xlim=c(-71,71),ylim=c(19-d$ry-ifelse(nrow(d$t)>25,17,10),59),expand=FALSE,clip='off')+
  labs(title=title)+theme_void(base_family='Arial',base_size=7.4)+
  theme(plot.title=element_text(size=9.2,face='bold',margin=margin(b=5)),
   plot.title.position='plot',plot.margin=margin(5,6,2,6),legend.position='bottom',
   legend.text=element_text(size=7.4),legend.key.width=unit(3,'mm'),legend.key.height=unit(3,'mm'),
   legend.spacing.x=unit(3,'mm'),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0))
}
saveplot<-function(p,name,w,h) {
 ggsave(file.path(out,paste0(name,'.png')),p,width=w/25.4,height=h/25.4,dpi=350,bg='white')
 ggsave(file.path(out,paste0(name,'.svg')),p,width=w/25.4,height=h/25.4,device=svglite::svglite,bg='white')
 if(capabilities('aqua')) {quartz(type='pdf',file=file.path(out,paste0(name,'.pdf')),width=w/25.4,height=h/25.4);print(p);dev.off()} else ggsave(file.path(out,paste0(name,'.pdf')),p,width=w/25.4,height=h/25.4,device=cairo_pdf)
}
saveplot(plot_net(v,'DOF links across elongating vascular contexts'),v$prefix,140,75)
saveplot(plot_net(b,'DOF links across the seven included tissue fits'),b$prefix,150,95)
for(cv in c('deutan','protan','tritan')) {
 ggsave(file.path(out,paste0('cvd_',cv,'_vascular.png')),plot_net(v,'DOF links across elongating vascular contexts',cv),width=140/25.4,height=75/25.4,dpi=180,bg='white')
 ggsave(file.path(out,paste0('cvd_',cv,'_all_tissues.png')),plot_net(b,'DOF links across the seven included tissue fits',cv),width=150/25.4,height=95/25.4,dpi=180,bg='white')
}
if(requireNamespace('colorblindr',quietly=TRUE)) colorblindr::cvd_grid(plot_net(v,'DOF links across elongating vascular contexts'))
fwrite(rbindlist(lapply(list(v,b),function(z)data.table(figure=z$prefix,t(z$counts)))),file.path(out,'validation.csv'))
writeLines(c(paste('Input',inp),paste('MD5',tools::md5sum(inp)),paste('Symbols',ann),paste('MD5',tools::md5sum(ann)),
 'No model fitting or new statistical inference. Gene family selection matches previous atlas.',
 'Union retains a source-target pair if present in any included context/genotype. No pooled coefficient or DE value.'),file.path(out,'input_provenance.txt'))
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
