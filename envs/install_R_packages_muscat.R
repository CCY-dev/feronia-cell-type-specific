# Re-install the R packages in conda env 'muscat' that were installed by R itself (not by conda).
# Generated 2026-10-20 from installed.packages() of /resnick/groups/NolanLab/cchen/miniconda3/envs/muscat.
# Run AFTER recreating the conda env from the matching .linux-64.explicit.txt, using that env's R.
# UNTESTED rebuild recipe: the exact copy of the env is the conda-pack archive deposited with the data.
options(repos = c(CRAN = 'https://cran.r-project.org'))
if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes')
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
pk <- read.csv(text = '
Package,Version,kind,source
base64enc,0.1-3,cran,CRAN
BiocManager,1.30.25,cran,CRAN
bluster,1.12.0,bioc,Bioconductor 3.18 (assumed; no Repository field)
brew,1.0-10,cran,CRAN
brio,1.1.5,cran,CRAN
bslib,0.9.0,cran,CRAN
callr,3.7.6,cran,CRAN
car,3.1-3,cran,CRAN
carData,3.0-5,cran,CRAN
cellranger,1.1.0,cran,CRAN
clipr,0.8.0,cran,CRAN
collections,0.3.11,cran,CRAN
commonmark,1.9.2,cran,CRAN
conflicted,1.2.0,cran,CRAN
corrplot,0.95,cran,CRAN
credentials,2.0.2,cran,CRAN
crosstalk,1.2.1,cran,CRAN
dbplyr,2.5.0,cran,CRAN
deldir,2.0-4,cran,CRAN
dendextend,1.19.1,cran,CRAN
desc,1.4.3,cran,CRAN
devtools,2.4.5,cran,CRAN
diffobj,0.3.5,cran,CRAN
dotCall64,1.2,cran,CRAN
downlit,0.4.4,cran,CRAN
DT,0.34.0,cran,CRAN
dtplyr,1.3.1,cran,CRAN
ellipse,0.5.0,cran,CRAN
emmeans,1.11.2-8,cran,CRAN
estimability,1.5.1,cran,CRAN
evaluate,1.0.3,cran,CRAN
factoextra,1.0.7,cran,CRAN
FactoMineR,2.12,cran,CRAN
fastDummies,1.7.5,cran,CRAN
fitdistrplus,1.2-2,cran,CRAN
flashClust,1.01-2,cran,CRAN
fontawesome,0.5.3,cran,CRAN
forcats,1.0.0,cran,CRAN
Formula,1.2-5,cran,CRAN
fs,1.6.5,cran,CRAN
gargle,1.5.2,cran,CRAN
GeneOverlap,1.38.0,bioc,Bioconductor 3.18 (assumed; no Repository field)
gert,2.1.4,cran,CRAN
ggforce,0.5.0,cran,CRAN
ggPlantmap,1.1.0,github,leonardojo/ggPlantmap@e0508172709dc192dfc3e99da631984dbd9a97d2
ggpubr,0.6.0,cran,CRAN
ggridges,0.5.6,cran,CRAN
ggsci,3.2.0,cran,CRAN
ggsignif,0.6.4,cran,CRAN
ggvenn,0.1.10,cran,CRAN
gh,1.4.1,cran,CRAN
gitcreds,0.1.2,cran,CRAN
goftest,1.2-3,cran,CRAN
googledrive,2.1.1,cran,CRAN
googlesheets4,1.1.1,cran,CRAN
gprofiler2,0.2.3,cran,CRAN
harmony,1.2.3,cran,CRAN
haven,2.5.4,cran,CRAN
hdf5r,1.3.12,cran,CRAN
here,1.0.1,cran,CRAN
highr,0.11,cran,CRAN
htmltools,0.5.8.1,cran,CRAN
htmlwidgets,1.6.4,cran,CRAN
httpuv,1.6.15,cran,CRAN
httr2,1.1.1,cran,CRAN
ica,1.0-3,cran,CRAN
ids,1.0.1,cran,CRAN
ini,0.3.1,cran,CRAN
IRdisplay,1.1,cran,CRAN
IRkernel,1.3.2,cran,CRAN
janeaustenr,1.0.0,cran,CRAN
jquerylib,0.1.4,cran,CRAN
knitr,1.49,cran,CRAN
later,1.4.1,cran,CRAN
lazyeval,0.2.2,cran,CRAN
leaps,3.2,cran,CRAN
leidenbase,0.1.32,cran,CRAN
lintr,3.3.0-1,cran,CRAN
lmtest,0.9-40,cran,CRAN
lubridate,1.9.4,cran,CRAN
mathjaxr,1.8-0,cran,CRAN
MatrixModels,0.5-3,cran,CRAN
metap,1.12,cran,CRAN
metapod,1.10.1,bioc,Bioconductor 3.18
miniUI,0.1.1.1,cran,CRAN
mnormt,2.1.1,cran,CRAN
multcomp,1.4-28,cran,CRAN
multcompView,0.1-10,cran,CRAN
multtest,2.58.0,bioc,Bioconductor 3.18 (assumed; no Repository field)
mutoss,0.1-13,cran,CRAN
patchwork,1.3.0,cran,CRAN
pbdZMQ,0.3-13,cran,CRAN
pkgbuild,1.4.6,cran,CRAN
pkgdown,2.1.1,cran,CRAN
pkgload,1.4.0,cran,CRAN
plotly,4.10.4,cran,CRAN
plotrix,3.8-4,cran,CRAN
polyclip,1.10-7,cran,CRAN
polynom,1.4-1,cran,CRAN
praise,1.0.0,cran,CRAN
processx,3.8.6,cran,CRAN
profvis,0.4.0,cran,CRAN
progressr,0.15.1,cran,CRAN
promises,1.3.2,cran,CRAN
proxy,0.4-27,cran,CRAN
ps,1.9.0,cran,CRAN
quantreg,6.1,cran,CRAN
R.cache,0.17.0,cran,CRAN
R.methodsS3,1.8.2,cran,CRAN
R.oo,1.27.1,cran,CRAN
R.utils,2.13.0,cran,CRAN
RANN,2.6.2,cran,CRAN
rappdirs,0.3.3,cran,CRAN
rcmdcheck,1.4.0,cran,CRAN
RcppTOML,0.2.3,cran,CRAN
readr,2.1.5,cran,CRAN
readxl,1.4.5,cran,CRAN
reformulas,0.4.0,cran,CRAN
rematch,2.0.0,cran,CRAN
rematch2,2.1.2,cran,CRAN
remotes,2.5.0,cran,CRAN
repr,1.1.7,cran,CRAN
reprex,2.1.1,cran,CRAN
reticulate,1.41.0.1,cran,CRAN
rex,1.2.1,cran,CRAN
rmarkdown,2.29,cran,CRAN
ROCR,1.0-11,cran,CRAN
roxygen2,7.3.2,cran,CRAN
rprojroot,2.0.4,cran,CRAN
RSpectra,0.16-2,cran,CRAN
rstatix,0.7.2,cran,CRAN
rstudioapi,0.17.1,cran,CRAN
rversions,2.1.2,cran,CRAN
rvest,1.0.4,cran,CRAN
S7,0.2.0,cran,CRAN
sandwich,3.1-1,cran,CRAN
sass,0.4.9,cran,CRAN
scatterplot3d,0.3-44,cran,CRAN
scran,1.30.2,bioc,Bioconductor 3.18
selectr,0.4-2,cran,CRAN
sessioninfo,1.2.3,cran,CRAN
Seurat,5.2.1,cran,CRAN
SeuratData,0.2.2.9002,github,satijalab/seurat-data@3e51f44303069b64f5dc4d68e6a3d4a343f55c39
SeuratDisk,0.0.0.9021,github,mojaveazure/seurat-disk@877d4e18ab38c686f5db54f8cd290274ccdbe295
SeuratObject,5.0.2,cran,CRAN
shiny,1.10.0,cran,CRAN
sn,2.1.1,cran,CRAN
SnowballC,0.7.1,cran,CRAN
sourcetools,0.1.7-1,cran,CRAN
sp,2.2-0,cran,CRAN
spam,2.11-1,cran,CRAN
SparseM,1.84-2,cran,CRAN
spatstat.data,3.1-4,cran,CRAN
spatstat.explore,3.3-4,cran,CRAN
spatstat.geom,3.3-5,cran,CRAN
spatstat.random,3.3-2,cran,CRAN
spatstat.sparse,3.1-0,cran,CRAN
spatstat.univar,3.1-2,cran,CRAN
spatstat.utils,3.1-2,cran,CRAN
tensor,1.5,cran,CRAN
testthat,3.2.3,cran,CRAN
TFisher,0.2.0,cran,CRAN
TH.data,1.1-3,cran,CRAN
tidytext,0.4.2,cran,CRAN
tidyverse,2.0.0,cran,CRAN
timechange,0.3.0,cran,CRAN
tinytex,0.56,cran,CRAN
tokenizers,0.3.0,cran,CRAN
tweenr,2.0.3,cran,CRAN
tzdb,0.4.0,cran,CRAN
UpSetR,1.4.1,cran,CRAN
urlchecker,1.0.1,cran,CRAN
usethis,3.1.0,cran,CRAN
uuid,1.2-1,cran,CRAN
vroom,1.6.5,cran,CRAN
waldo,0.6.1,cran,CRAN
whisker,0.4.1,cran,CRAN
xfun,0.51,cran,CRAN
XML,3.99-0.23,cran,CRAN
xmlparsedata,1.0.5,cran,CRAN
xopen,1.0.1,cran,CRAN
xtable,1.8-4,cran,CRAN
yaml,2.3.10,cran,CRAN
zip,2.3.2,cran,CRAN
zoo,1.8-13,cran,CRAN')
bioc_version <- '3.18'
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
