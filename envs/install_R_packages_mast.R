# Re-install the R packages in conda env 'mast' that were installed by R itself (not by conda).
# Generated 2026-10-20 from installed.packages() of /resnick/groups/NolanLab/cchen/miniconda3/envs/mast.
# Run AFTER recreating the conda env from the matching .linux-64.explicit.txt, using that env's R.
# UNTESTED rebuild recipe: the exact copy of the env is the conda-pack archive deposited with the data.
options(repos = c(CRAN = 'https://cran.r-project.org'))
if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes')
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
pk <- read.csv(text = '
Package,Version,kind,source
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
microbenchmark,1.5.0,cran,CRAN
pbdZMQ,0.3-14,cran,CRAN
pbkrtest,0.5.5,cran,CRAN
permute,0.9-8,cran,CRAN
polynom,1.4-1,cran,CRAN
quantreg,6.1,cran,CRAN
remotes,2.5.0,cran,CRAN
repr,1.1.7,cran,CRAN
rjson,0.2.23,cran,CRAN
rstatix,0.7.3,cran,CRAN
SeuratDisk,0.0.0.9021,github,mojaveazure/seurat-disk@877d4e18ab38c686f5db54f8cd290274ccdbe295
shape,1.4.6.1,cran,CRAN
SnowballC,0.7.1,cran,CRAN
SparseM,1.84-2,cran,CRAN
tidytext,0.4.3,cran,CRAN
tokenizers,0.3.0,cran,CRAN
vegan,2.7-2,cran,CRAN')
bioc_version <- '3.20'
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
