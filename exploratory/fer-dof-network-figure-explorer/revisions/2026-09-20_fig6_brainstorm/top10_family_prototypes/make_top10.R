#!/usr/bin/env Rscript
# Prototype; not yet carefully vetted.
# Run from the repository root. No model fitting or hypothesis testing is rerun.
suppressPackageStartupMessages({library(ggplot2);library(cowplot);library(colorspace)})
pdfFonts(Arial=pdfFonts()$Helvetica)
options(device=function(...) grDevices::pdf(file=NULL,...))
arg <- commandArgs(trailingOnly=FALSE)
script <- normalizePath(sub('^--file=','',arg[grepl('^--file=',arg)][1]))
out <- dirname(script)
repo <- normalizePath(file.path(out,'../../..'))
source_file <- file.path(repo,'.living/tasks/fer-fig6-brainstorm-2026-09-20/outputs/claude-science/data/family_enrichment_paired_universe.csv')
d <- read.csv(source_file,stringsAsFactors=FALSE)
stopifnot(nrow(d)==110, all(d$universe_total==1292),all(d$definition=='PlantTFDB'),
          all(d$universe=='paired_score_KayTF'),!anyDuplicated(paste(d$Family,d$genotype)))
# Symmetric, explicit display selection. BH values retain the source's 55-family correction.
eligible <- d[d$universe_n>=10,]
ranking <- aggregate(odds_ratio~Family,eligible,max)
ranking <- ranking[order(-ranking$odds_ratio,ranking$Family),]
ranking$display_rank <- seq_len(nrow(ranking))
ranking$selected <- ranking$display_rank<=10
names(ranking)[names(ranking)=='odds_ratio'] <- 'max_genotype_OR'
write.csv(ranking,file.path(out,'all_eligible_family_ranking.csv'),row.names=FALSE)
ord <- ranking$Family[ranking$selected]
z <- d[d$Family %in% ord,]
z$display_rank <- match(z$Family,ord)
z <- z[order(z$display_rank,z$genotype),]
z$row <- 11-z$display_rank
z$set <- factor(ifelse(z$genotype=='WT_specific','WT-associated','fer-4-associated'),
                levels=c('WT-associated','fer-4-associated'))
z$evidence <- factor(ifelse(z$p_adj_BH<.05,'BH < 0.05','BH >= 0.05'),
                     levels=c('BH < 0.05','BH >= 0.05'))
z$count <- sprintf('%d/%d',z$foreground_n,z$universe_n)
z$point_x <- pmax(z$odds_ratio,.125)
z$ci_low_x <- pmax(z$odds_lower,.125)
z$log_OR <- ifelse(z$odds_ratio>0,log2(z$odds_ratio),NA_real_)
z$OR_text <- ifelse(z$odds_ratio==0,'0',sprintf('%.2f',z$odds_ratio))
write.csv(z,file.path(out,'plot_input_top10.csv'),row.names=FALSE)
label_family <- function(x) ifelse(x=='Dof','DOF',gsub('_',' ',x))
labels <- label_family(rev(ord))
theme_set(theme_cowplot(font_size=8,font_family='Arial') +
 theme(plot.margin=margin(2,2,2,2),axis.text=element_text(size=8),
       axis.title=element_text(size=8.5),plot.title=element_text(size=9,face='bold')))
cols <- c('WT-associated'='#0072B2','fer-4-associated'='#D55E00')
sim_cols <- function(x,cvd=NULL) if(is.null(cvd)) x else get(cvd,asNamespace('colorspace'))(x)
base_y <- function(label=FALSE) scale_y_continuous(limits=c(.45,10.55),breaks=1:10,
   labels=if(label) labels else NULL,expand=c(0,0))
forest_panel <- function(set,cvd=NULL){
 x <- z[z$set==set,]
 color <- unname(sim_cols(cols,cvd)[match(set,names(cols))])
 p <- ggplot(x,aes(y=row)) +
   geom_vline(xintercept=1,linetype='dashed',colour='grey55',size=.35) +
   geom_segment(aes(x=ci_low_x,xend=odds_upper,yend=row),colour=color,size=.6) +
   geom_point(data=x[x$odds_ratio>0,],aes(x=odds_ratio,shape=evidence),colour=color,size=2.3,stroke=.7) +
   scale_shape_manual(values=c('BH < 0.05'=16,'BH >= 0.05'=1),drop=FALSE) +
   scale_x_log10(limits=c(.125,24),breaks=c(.125,.5,1,2,4,8,16),
      labels=c('0.125','0.5','1','2','4','8','16'),expand=c(0,0)) +
   base_y(FALSE) + labs(x=NULL,y=NULL,title=set) +
   theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),
      axis.text.x=element_text(size=7),plot.title=element_text(colour=color,size=9),
      panel.grid.major.x=element_line(colour='grey94',size=.25),legend.position='none')
 if(any(x$odds_ratio==0)) p <- p + geom_segment(data=x[x$odds_ratio==0,],aes(x=.19,xend=.125,yend=row),colour=color, arrow=grid::arrow(length=grid::unit(1.5,'mm'),type='closed'),size=.6) + geom_text(data=x[x$odds_ratio==0,],aes(x=.145,label='OR = 0'),hjust=0,vjust=-.8,size=2.45)
 p
}
count_panel <- function(set){
 x <- z[z$set==set,]
 ggplot(x,aes(x=1,y=row,label=count))+geom_text(hjust=1,size=2.75,family='Arial')+
   base_y(FALSE)+scale_x_continuous(limits=c(0,1.1),expand=c(0,0))+
   labs(title='n/K',x=NULL,y=NULL)+theme_void(base_family='Arial')+
   theme(plot.title=element_text(size=8,hjust=.9),plot.margin=margin(2,3,2,0))
}
family_panel <- ggplot(data.frame(row=1:10,label=labels),aes(x=1,y=row,label=label))+
 geom_text(hjust=1,size=3,family='Arial')+base_y(FALSE)+
 scale_x_continuous(limits=c(0,1.05),expand=c(0,0))+theme_void()+
 labs(title='TF family')+theme(plot.title=element_text(size=8,hjust=1),plot.margin=margin(2,2,2,0))
forest <- function(cvd=NULL){
 panels <- plot_grid(family_panel,forest_panel('WT-associated',cvd),count_panel('WT-associated'),
                    forest_panel('fer-4-associated',cvd),count_panel('fer-4-associated'),
                    nrow=1,rel_widths=c(1.12,2.46,.62,2.46,.62),align='h',axis='tb')
 # Vertical panel alignment is explicit and all panels share the same continuous y limits.
 ggdraw()+draw_plot(panels,x=0,y=.18,width=1,height=.66)+
   draw_label('Top 10 TF families by enrichment effect size',x=.015,y=.975,hjust=0,vjust=1,
       fontfamily='Arial',fontface='bold',size=11)+
   draw_label('Same families in both sets; ranked by the larger of the two odds ratios',
       x=.015,y=.91,hjust=0,vjust=1,fontfamily='Arial',size=8.3)+
   draw_label('Odds ratio (95% Fisher confidence interval; log scale)',x=.58,y=.14,
       fontfamily='Arial',size=8.5)+
   draw_label('Filled: BH < 0.05    Open: BH >= 0.05    n/K: nominated / eligible family members',
       x=.015,y=.08,hjust=0,fontfamily='Arial',size=7.4)+
   draw_label('>=10 eligible members per family; BH across all 55 families. Prior-conditional nomination enrichment.',
       x=.015,y=.035,hjust=0,fontfamily='Arial',size=7.1)
}
matrix_plot <- function(cvd=NULL){
 # Fixed-size cells; effect size is fill, evidence is border, no significance-based bubble size.
 pal <- sim_cols(colorspace::sequential_hcl(100,palette='Viridis'),cvd)
 x <- z
 x$xx <- as.numeric(x$set)
 x$text_col <- ifelse(!is.na(x$log_OR)&x$log_OR<1.2,'white','black')
 x$text_col[is.na(x$log_OR)] <- 'black'
 g <- ggplot(x,aes(xx,row))+
   geom_tile(aes(fill=log_OR),width=.86,height=.87,colour='white',size=.3)+
   geom_tile(data=x[x$p_adj_BH<.05,],fill=NA,width=.86,height=.87,colour='black',size=.5)+
   geom_text(aes(label=paste(OR_text,count,sep=' | '),colour=text_col),size=3.05,family='Arial')+
   scale_colour_identity()+
   scale_fill_gradientn(colours=pal,limits=c(-1,3),breaks=c(-1,0,1,2,3),
       labels=c('0.5','1','2','4','8'),na.value='grey92',name='Odds ratio')+
   guides(fill=guide_colorbar(title.position='top',title.hjust=.5))+
   scale_x_continuous(breaks=1:2,labels=c('WT-associated','fer-4-associated'),position='top',
       limits=c(.45,2.55),expand=c(0,0))+
   base_y(TRUE)+labs(x=NULL,y=NULL)+
   theme_minimal(base_family='Arial',base_size=9)+
   theme(panel.grid=element_blank(),axis.text.y=element_text(size=9,colour='black'),
      axis.text.x=element_text(size=8.8,colour='black'),legend.position='bottom',
      legend.title=element_text(size=8),legend.text=element_text(size=7.5),
      legend.key.width=grid::unit(.32,'in'),legend.key.height=grid::unit(.08,'in'),
      legend.margin=margin(0,0,0,0),plot.margin=margin(0,3,0,0))
 ggdraw()+draw_plot(g,x=.015,y=.15,width=.97,height=.68)+
   draw_label('Top 10 TF families',x=.02,y=.985,hjust=0,vjust=1,
     fontfamily='Arial',fontface='bold',size=11)+
   draw_label('Ranked by the larger WT or fer-4 odds ratio',x=.02,y=.935,hjust=0,vjust=1,
     fontfamily='Arial',size=8)+
   draw_label('Cell text: odds ratio | nominated / eligible',x=.02,y=.89,hjust=0,vjust=1,
     fontfamily='Arial',size=8)+
   draw_label('Black border: BH < 0.05 across all 55 families',x=.02,y=.105,hjust=0,
     fontfamily='Arial',size=7.5)+
   draw_label('>=10 eligible members per family; zero is grey.',x=.02,y=.066,hjust=0,
     fontfamily='Arial',size=7.3)+
   draw_label('Prior-conditional nomination enrichment; CIs in forest plot.',x=.02,y=.027,hjust=0,
     fontfamily='Arial',size=7.1)
}
save_all <- function(p,name,width,height){
 # Native Quartz PDF is used because local Cairo requires an unavailable XQuartz library.
 grDevices::quartz(type='pdf',file=file.path(out,paste0(name,'.pdf')),width=width,height=height,
                  family='Arial',bg='white');print(p);dev.off()
 svg_path <- file.path(out,paste0(name,'.svg'))
 ggsave(svg_path,p,device=svglite::svglite,width=width,height=height,bg='white')
 # Name generated clipping groups for editable vector handoff.
 svg_lines <- readLines(svg_path,warn=FALSE)
 group_lines <- grep('^<g ',svg_lines)
 for(i in seq_along(group_lines)) svg_lines[group_lines[i]] <- sub('<g ',paste0('<g id=\"',name,'-region-',sprintf('%03d',i),'\" '),svg_lines[group_lines[i]],fixed=TRUE)
 writeLines(svg_lines,svg_path)
 ggsave(file.path(out,paste0(name,'.png')),p,width=width,height=height,dpi=250,bg='white')
}
save_all(forest(),'top10_aligned_forest',6.9,4.3)
save_all(matrix_plot(),'top10_enrichment_matrix',4.2,4.6)
# Mandatory CVD check. colorblindr is not installed locally, so colorspace simulates
# the complete final color encodings under deutan, protan and tritan deficiency.
for(cvd in c('deutan','protan','tritan')){
 ggsave(file.path(out,paste0('cvd_',cvd,'_forest.png')),forest(cvd),width=6.9,height=4.3,dpi=150,bg='white')
 ggsave(file.path(out,paste0('cvd_',cvd,'_matrix.png')),matrix_plot(cvd),width=4.2,height=4.6,dpi=150,bg='white')
}
if(requireNamespace('colorblindr',quietly=TRUE)) colorblindr::cvd_grid(matrix_plot())
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
stopifnot(length(ord)==10,nrow(z)==20,sum(z$p_adj_BH<.05)==6,
  sum(z$odds_ratio==0)==1,all(z$foreground_n<=z$universe_n))
cat('Saved top10 forest and matrix. Selection:',paste(ord,collapse=', '),'\n')
