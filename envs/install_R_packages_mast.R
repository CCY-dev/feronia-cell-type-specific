# Re-install, at the recorded versions, the R packages in conda env 'mast' that conda did not provide
# (packages installed by R itself, and conda packages that R later replaced with another version).
# Generated 2026-10-20 from installed.packages() of /resnick/groups/NolanLab/cchen/miniconda3/envs/mast.
# Run AFTER recreating the conda env from envs/mast.linux-64.explicit.txt, with that env's Rscript.
# Exits with status 1 and lists any package that is not at the recorded version.
options(repos = c(CRAN = 'https://cran.r-project.org'), timeout = 600)
lib <- .libPaths()[1]
if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes', INSTALL_opts = '--no-lock')
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', INSTALL_opts = '--no-lock')
pk <- read.csv(text = 'Package,Version,kind,source
car,3.1-3,cran,CRAN
carData,3.0-5,cran,CRAN
circlize,0.4.16,cran,CRAN
clue,0.3-66,cran,CRAN
ComplexHeatmap,2.22.0,bioc,Bioconductor 3.20
corrplot,0.95,cran,CRAN
Deriv,4.2.0,cran,CRAN
doBy,4.7.0,cran,CRAN
doParallel,1.0.17,cran,CRAN
foreach,1.5.2,cran,CRAN
Formula,1.2-5,cran,CRAN
GetoptLong,1.0.5,cran,CRAN
ggpubr,0.6.2,cran,CRAN
ggsci,4.1.0,cran,CRAN
ggsignif,0.6.4,cran,CRAN
ggvenn,0.1.19,cran,CRAN
GlobalOptions,0.1.2,cran,CRAN
hdf5r,1.3.12,cran,CRAN
inflection,1.3.7,cran,CRAN
IRdisplay,1.1,cran,CRAN
IRkernel,1.3.2,cran,CRAN
iterators,1.0.14,cran,CRAN
janeaustenr,1.0.0,cran,CRAN
leiden,0.4.3.1,cran,CRAN
MatrixModels,0.5-4,cran,CRAN
mclust,6.1.1,cran,CRAN
mgcv,1.9-3,cran,CRAN
microbenchmark,1.5.0,cran,CRAN
nnet,7.3-20,cran,CRAN
pbdZMQ,0.3-14,cran,CRAN
pbkrtest,0.5.5,cran,CRAN
permute,0.9-8,cran,CRAN
polynom,1.4-1,cran,CRAN
quantreg,6.1,cran,CRAN
remotes,2.5.0,cran,CRAN
repr,1.1.7,cran,CRAN
rjson,0.2.23,cran,CRAN
rstatix,0.7.3,cran,CRAN
Seurat,5.0.0,github,satijalab/seurat@a4075e78d52a67e5f901459c6eb024bb8b3c9c44
SeuratDisk,0.0.0.9021,github,mojaveazure/seurat-disk@877d4e18ab38c686f5db54f8cd290274ccdbe295
shape,1.4.6.1,cran,CRAN
SnowballC,0.7.1,cran,CRAN
SparseM,1.84-2,cran,CRAN
tidytext,0.4.3,cran,CRAN
tokenizers,0.3.0,cran,CRAN
vegan,2.7-2,cran,CRAN', stringsAsFactors = FALSE)
bioc_version <- '3.20'
at_version <- function(p, v) { d <- suppressWarnings(packageDescription(p, lib.loc = lib)); is.list(d) && identical(d$Version, v) }
for (pass in 1:8) {
  unlink(Sys.glob(file.path(lib, '00LOCK*')), recursive = TRUE)   # clear locks left by failed installs
  todo <- pk[!mapply(at_version, pk$Package, pk$Version), ]
  if (!nrow(todo)) break
  message('pass ', pass, ': ', nrow(todo), ' packages to install')
  for (i in seq_len(nrow(todo))) {
    p <- todo$Package[i]; v <- todo$Version[i]; s <- todo$source[i]
    try(switch(todo$kind[i],
      cran   = remotes::install_version(p, version = v, upgrade = 'never', dependencies = FALSE, INSTALL_opts = '--no-lock'),
      bioc   = BiocManager::install(p, version = bioc_version, update = FALSE, ask = FALSE, dependencies = FALSE, INSTALL_opts = '--no-lock'),
      github = remotes::install_github(sub('@.*', '', s), ref = sub('^.*@', '', s), upgrade = 'never', dependencies = FALSE, INSTALL_opts = '--no-lock')
    ), silent = TRUE)
  }
}
bad <- pk$Package[!mapply(at_version, pk$Package, pk$Version)]
if (length(bad)) { message('Not at the recorded version: ', paste(bad, collapse = ', ')); quit(status = 1) } else message('All ', nrow(pk), ' packages at recorded versions.')
