#!/usr/bin/env Rscript
# Prototype; not yet carefully vetted.
# Descriptive programme membership of fitted edges for the ten assayed DOFs.
# No enrichment tests, no causal edge claims, no inference of TF activity.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE)
repo <- if(length(args)) normalizePath(args[1]) else normalizePath('.')
base <- file.path(repo, 'revisions/2026-09-20_dof_subnetworks')
outdir <- file.path(base, 'annotations')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
write_tsv <- function(x, n) fwrite(x, file.path(outdir, n), sep='\t', na='NA')
gobase <- file.path(base, 'reference/gene_ontology')
obop <- file.path(gobase, 'go-basic.obo')
gafp <- file.path(gobase, 'ARATH-mod.gaf.gz')
q04 <- file.path(repo, 'revisions/2026-09-08_routeA_execution/q04_network/results')
edgep <- file.path(q04, 'DOF_edges_and_all_YUC5_edges.csv')
testedp <- file.path(q04, 'tested_DOF_integrated_evidence_matrix.csv')
# Freeze programme definitions before consulting target overlaps.
programmes <- data.table(
 programme=c('auxin_biosynthesis','auxin_transport','response_to_auxin',
             'auxin_catabolism','auxin_conjugate_metabolism','cell_wall',
             'root_hair_tip_growth','mitotic_cell_cycle'),
 root_go_id=c('GO:0009851','GO:0060918','GO:0009733','GO:0009852',
              'GO:0010249','GO:0071554','GO:0048768','GO:0000278'))
ln <- readLines(obop, warn=FALSE)
starts <- which(ln=='[Term]'); ends <- c(starts[-1]-1L, length(ln))
tr <- list(); ed <- list(); al <- list()
for(i in seq_along(starts)) {
 b <- ln[starts[i]:ends[i]]
 id <- sub('^id: ', '', b[grepl('^id: GO:', b)])[1]
 if(is.na(id)) next
 obsolete <- any(b=='is_obsolete: true')
 tr[[i]] <- data.table(id=id, term_name=sub('^name: ', '', b[grepl('^name: ', b)])[1],
                      namespace=sub('^namespace: ', '', b[grepl('^namespace: ',b)])[1], obsolete=obsolete)
 pa <- b[grepl('^is_a: GO:|^relationship: part_of GO:', b)]
 if(length(pa)&&!obsolete) ed[[i]] <- data.table(child=id, parent=sub('.*(GO:[0-9]{7}).*', '\\1', pa),
                                              relation=ifelse(grepl('^is_a:',pa),'is_a','part_of'))
 aa <- sub('^alt_id: ', '', b[grepl('^alt_id: GO:', b)])
 if(length(aa)) al[[i]] <- data.table(alt_id=aa, id=id)
}
terms <- rbindlist(tr); edges <- unique(rbindlist(ed)); alts <- rbindlist(al)
edges <- edges[child %in% terms[obsolete==FALSE,id] & parent %in% terms[obsolete==FALSE,id]]
programmes <- merge(programmes,terms[,.(root_go_id=id,root_term_name=term_name,namespace,obsolete)],by='root_go_id',sort=FALSE)
stopifnot(nrow(programmes)==8L, !any(programmes$obsolete), all(programmes$namespace=='biological_process'))
write_tsv(programmes,'programme_definitions.tsv')
gaf_con <- gzfile(gafp)
gl <- readLines(gaf_con,warn=FALSE); close(gaf_con)
headers <- gl[startsWith(gl,'!')]
gaf <- fread(text=paste(gl[!startsWith(gl,'!')],collapse='\n'),header=FALSE,sep='\t',quote='',fill=TRUE)
setnames(gaf,paste0('V',1:17),c('database','gene','symbol','qualifier','go_id','reference','evidence',
 'with','aspect','name','synonym','type','taxon','date','assigned_by','extension','form_id'))
# Symbols are data-source labels, not a separate gene-identity inference.
symbols <- unique(gaf[,.(gene,symbol,name,synonym)])
write_tsv(symbols,'gaf_gene_symbols.tsv')
gaf[, original_go_id:=go_id]
altmap <- setNames(alts$id,alts$alt_id)
ix <- which(gaf$go_id %in% names(altmap))
gaf[ix,go_id:=unname(altmap[go_id])]
write_tsv(unique(gaf[!go_id %in% terms[obsolete==FALSE,id],.(go_id,original_go_id)]),'unresolved_or_obsolete_annotations.tsv')
gaf <- gaf[aspect=='P' & !grepl('(^|\\|)NOT($|\\|)',qualifier) &
           go_id %in% terms[obsolete==FALSE,id] & grepl('^AT[1-5CM]G[0-9]{5}$',gene)]
ann <- list(); descs <- list()
for(i in seq_len(nrow(programmes))) {
 desc <- programmes$root_go_id[i]
 repeat { more <- union(desc, edges[parent %in% desc,child]); if(length(more)==length(desc)) break; desc <- more }
 descs[[i]] <- merge(data.table(programme=programmes$programme[i],id=desc),terms,by='id',all.x=TRUE)
 z <- copy(gaf[go_id %in% desc]); z[,`:=`(programme=programmes$programme[i],root_go_id=programmes$root_go_id[i])]
 ann[[i]] <- z
}
annotation <- rbindlist(ann)
write_tsv(rbindlist(descs),'programme_descendant_terms.tsv')
write_tsv(annotation,'programme_annotation_evidence.tsv')
members <- unique(annotation[,.(programme,gene,symbol)])
write_tsv(members,'programme_gene_membership.tsv')
write_tsv(annotation[,.(n_annotations=.N,n_genes=uniqueN(gene)),by=.(programme,evidence)],'annotation_evidence_code_counts.tsv')
tested <- fread(testedp)[,.(source=gene,tested_symbol=symbol)]
dat <- fread(edgep)
dat <- merge(dat,tested,by='source')
stopifnot(uniqueN(dat$source)==10L)
dat[,genotype:=fcase(grepl('_WT_untreated$',cluster),'WT',grepl('_fer4_untreated$',cluster),'fer4',default=NA_character_)]
stopifnot(!anyNA(dat$genotype))
dat[,context:=sub('_(WT|fer4)_untreated$','',cluster)]
dat[,clean_analysis:=!tissue %in% c('phloem','xylem','endodermis') & !target_excluded_by_explicit_map]
write_tsv(dat[,.(edge_rows=.N,n_targets=uniqueN(target)),by=.(tissue,clean_analysis)],'tissue_inclusion_audit.tsv')
clean <- dat[clean_analysis==TRUE]
stopifnot(!anyDuplicated(clean[,.(tissue,source,target,cluster)]),all(is.finite(clean$coef_mean)))
write_tsv(clean[,.(total_unique_targets=uniqueN(target),contexts=uniqueN(context),genotype_contexts=uniqueN(cluster),edge_rows=.N),
 by=.(source,tested_symbol)],'tested_dof_clean_total_targets.tsv')
joined <- merge(clean,members,by.x='target',by.y='gene',allow.cartesian=TRUE)
write_tsv(joined,'clean_programme_target_edges.tsv')
counts <- joined[,.(n_targets=uniqueN(target),n_contexts=uniqueN(context),n_genotype_contexts=uniqueN(cluster),
                    genes=paste(sort(unique(paste0(symbol,' [',target,']'))),collapse='; ')),
                  by=.(source,tested_symbol,programme)]
fullgrid <- CJ(source=tested$source,programme=programmes$programme)
counts <- merge(merge(fullgrid,tested,by='source'),counts,by=c('source','tested_symbol','programme'),all.x=TRUE)
for(v in c('n_targets','n_contexts','n_genotype_contexts')) set(counts,which(is.na(counts[[v]])),v,0L)
write_tsv(counts,'tested_dof_programme_target_counts.tsv')
ctcounts <- joined[,.(n_targets=uniqueN(target),genes=paste(sort(unique(paste0(symbol,' [',target,']'))),collapse='; ')),
 by=.(source,tested_symbol,programme,tissue,context,genotype,cluster)]
setorder(ctcounts,programme,-n_targets,tested_symbol,cluster)
write_tsv(ctcounts,'context_programme_target_counts.tsv')
unioncounts <- joined[,.(n_targets=uniqueN(target),n_tfs=uniqueN(source),n_contexts=uniqueN(context)),by=programme]
unioncounts <- merge(programmes[,.(programme)],unioncounts,by='programme',all.x=TRUE)
for(v in c('n_targets','n_tfs','n_contexts')) set(unioncounts,which(is.na(unioncounts[[v]])),v,0L)
write_tsv(unioncounts,'programme_clean_union_summary.tsv')
targetinfo <- joined[,.(tested_dofs=paste(sort(unique(tested_symbol)),collapse='; '),n_tested_dofs=uniqueN(source),
 tissues=paste(sort(unique(tissue)),collapse='; '),n_contexts=uniqueN(context),n_genotype_contexts=uniqueN(cluster)),
 by=.(programme,target,symbol)]
setorder(targetinfo,programme,-n_tested_dofs,-n_contexts,symbol)
write_tsv(targetinfo,'programme_target_gene_summary.tsv')
# Evidence supporting the observed target annotation only, including all evidence codes and qualifiers.
write_tsv(annotation[gene %in% clean$target],'clean_target_annotation_evidence.tsv')
aux_examples <- joined[tested_symbol %in% c('DOF4.6','DAG1','OBP3','COG1') & programme=='response_to_auxin',
 .(n_targets=uniqueN(target),genes=paste(sort(unique(paste0(symbol,' [',target,']'))),collapse='; ')),
 by=.(tested_symbol,context,genotype)][order(tested_symbol,-n_targets)][,head(.SD,3),by=tested_symbol]
write_tsv(aux_examples,'selected_tf_auxin_response_context_examples.tsv')
# Narrow family lookup is recorded separately from broad GO process annotation.
taa_yuc_symbols <- unique(symbols[grepl('^YUC[0-9]+$|^TAA1$|^TAR[12]$',symbol),.(gene,symbol)])
write_tsv(taa_yuc_symbols,'taa_tar_yuc_family_lookup.tsv')
ty <- merge(clean,taa_yuc_symbols,by.x='target',by.y='gene')
write_tsv(ty,'clean_taa_tar_yuc_edges.tsv')
tymap <- ty[,.(genes=paste(sort(unique(paste0(symbol,' [',target,']'))),collapse='; '),
 tissues=paste(sort(unique(tissue)),collapse='; '),n_genotype_contexts=uniqueN(cluster)),by=.(source,tested_symbol)]
tymap <- merge(tested,tymap,by=c('source','tested_symbol'),all.x=TRUE)
tymap[is.na(genes),`:=`(genes='No retained edge in the seven included tissues',tissues='',n_genotype_contexts=0L)]
write_tsv(tymap,'tested_dof_taa_tar_yuc_summary.tsv')
capture.output(sessionInfo(),file=file.path(outdir,'sessionInfo.txt'))
writeLines(c(ln[grepl('^(format-version|data-version):',ln)],headers,
 'Positive biological-process annotations; all evidence codes retained and explicitly tabulated.',
 'Alternative IDs resolved; obsolete/unresolved annotations excluded; descendants through is_a/part_of only.',
 'No regulates/has_part traversal. Annotation qualifiers are preserved; some describe acts_upstream_of_or_within.',
 'Exclude phloem and xylem for target-map defects, and all endodermis for score-edge mismatch under investigation.',
 'Exclude any row target_excluded_by_explicit_map=TRUE. Fit export edges are descriptive, without newly imposed coefficient threshold.',
 'Programme gene sets overlap. Membership is descriptive, not enrichment and not evidence of direct binding, activity, or causation.'),
 file.path(outdir,'provenance_and_rules.txt'))
paths <- c(obop,gafp,edgep,testedp)
hashes <- vapply(paths,function(p) system2('/usr/bin/shasum',c('-a','256',shQuote(p)),stdout=TRUE),character(1))
writeLines(hashes,file.path(outdir,'inputs.sha256'))
print(unioncounts)
print(dcast(counts,tested_symbol~programme,value.var='n_targets'))
