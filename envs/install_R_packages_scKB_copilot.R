# Re-install the R packages in conda env 'scKB_copilot' that were installed by R itself (not by conda).
# Generated 2026-10-20 from installed.packages() of /resnick/groups/NolanLab/cchen/miniconda3/envs/scKB_copilot.
# Run AFTER recreating the conda env from the matching .linux-64.explicit.txt, using that env's R.
# UNTESTED rebuild recipe: the exact copy of the env is the conda-pack archive deposited with the data.
options(repos = c(CRAN = 'https://cran.r-project.org'))
if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes')
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
pk <- read.csv(text = '
Package,Version,kind,source
BiocManager,1.30.25,cran,CRAN
broom,1.0.7,cran,CRAN
BSgenome.Athaliana.TAIR.TAIR9,1.3.1000,bioc,Bioconductor 3.14 (assumed; no Repository field)
caTools,1.18.3,cran,CRAN
cellranger,1.1.0,cran,CRAN
conflicted,1.2.0,cran,CRAN
cowplot,1.1.3,cran,CRAN
data.table,1.17.0,cran,CRAN
dotCall64,1.2,cran,CRAN
DoubletFinder,2.0.3,github,ekernf01/DoubletFinder@11d18827cf6595da7fd8b9f511e5672beeaf7b84
dtplyr,1.3.1,cran,CRAN
fields,16.3.1,cran,CRAN
fitdistrplus,1.2-2,cran,CRAN
FNN,1.1.4.1,cran,CRAN
forcats,1.0.0,cran,CRAN
future,1.34.0,cran,CRAN
future.apply,1.11.3,cran,CRAN
gargle,1.5.2,cran,CRAN
ggrepel,0.9.6,cran,CRAN
ggridges,0.5.6,cran,CRAN
globals,0.16.3,cran,CRAN
googledrive,2.1.1,cran,CRAN
googlesheets4,1.1.1,cran,CRAN
gplots,3.2.0,cran,CRAN
gridExtra,2.3,cran,CRAN
gtools,3.9.5,cran,CRAN
haven,2.5.4,cran,CRAN
httpuv,1.6.16,cran,CRAN
ica,1.0-3,cran,CRAN
ids,1.0.1,cran,CRAN
igraph,2.1.4,cran,CRAN
IRdisplay,1.1,cran,CRAN
IRkernel,1.3.2,cran,CRAN
irlba,2.3.5.1,cran,CRAN
leiden,0.4.3.1,cran,CRAN
listenv,0.9.1,cran,CRAN
lmtest,0.9-40,cran,CRAN
lubridate,1.9.4,cran,CRAN
maps,3.4.2.1,cran,CRAN
modelr,0.1.11,cran,CRAN
otel,0.2.0,cran,CRAN
parallelly,1.42.0,cran,CRAN
patchwork,1.3.0,cran,CRAN
pbapply,1.7-2,cran,CRAN
pbdZMQ,0.3-13,cran,CRAN
plotly,4.10.4,cran,CRAN
plyr,1.8.9,cran,CRAN
ragg,1.3.3,cran,CRAN
RANN,2.6.2,cran,CRAN
RcppAnnoy,0.0.22,cran,CRAN
RcppEigen,0.3.4.0.2,cran,CRAN
readr,2.1.5,cran,CRAN
readxl,1.4.5,cran,CRAN
rematch,2.0.0,cran,CRAN
repr,1.1.7,cran,CRAN
reprex,2.1.1,cran,CRAN
reshape2,1.4.4,cran,CRAN
ROCR,1.0-11,cran,CRAN
RSpectra,0.16-2,cran,CRAN
rsvd,1.0.5,cran,CRAN
Rtsne,0.17,cran,CRAN
rvest,1.0.4,cran,CRAN
sctransform,0.4.1,cran,CRAN
selectr,0.4-2,cran,CRAN
Seurat,3.1.5,cran,CRAN
shiny,1.13.0,cran,CRAN
sourcetools,0.1.7-1,cran,CRAN
spam,2.11-1,cran,CRAN
systemfonts,1.2.1,cran,CRAN
textshaping,1.0.0,cran,CRAN
tidyverse,2.0.0,cran,CRAN
timechange,0.3.0,cran,CRAN
tsne,0.1-3.1,cran,CRAN
tzdb,0.4.0,cran,CRAN
uuid,1.2-1,cran,CRAN
uwot,0.2.3,cran,CRAN
vroom,1.6.5,cran,CRAN
xtable,1.8-8,cran,CRAN
zoo,1.8-13,cran,CRAN')
bioc_version <- '3.14'
todo <- pk
for (pass in 1:5) {  # retry so packages whose dependencies are installed later succeed on a later pass
  for (i in seq_len(nrow(todo))) {
    p <- todo$Package[i]; v <- todo$Version[i]
    if (requireNamespace(p, quietly = TRUE) && as.character(packageVersion(p)) == v) next
    try(switch(todo$kind[i],
      cran   = remotes::install_version(p, version = v, upgrade = 'never', dependencies = FALSE),
      bioc   = BiocManager::install(p, version = bioc_version, update = FALSE, ask = FALSE, dependencies = FALSE),
      github = remotes::install_github(sub('@.*', '', todo$source[i]), ref = sub('^.*@', '', todo$source[i]), upgrade = 'never', dependencies = FALSE)
    ), silent = TRUE)
  }
}
ip <- installed.packages()[, 'Version']
bad <- pk$Package[is.na(ip[pk$Package]) | ip[pk$Package] != pk$Version]
if (length(bad)) { message('Not at the recorded version: ', paste(bad, collapse = ', ')); quit(status = 1) } else message('All ', nrow(pk), ' R-installed packages at recorded versions.')
