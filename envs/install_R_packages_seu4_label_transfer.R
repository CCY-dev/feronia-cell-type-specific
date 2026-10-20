# Re-install the R packages in conda env 'seu4_label_transfer' that were installed by R itself (not by conda).
# Generated 2026-10-20 from installed.packages() of /resnick/groups/NolanLab/cchen/miniconda3/envs/seu4_label_transfer.
# Run AFTER recreating the conda env from the matching .linux-64.explicit.txt, using that env's R.
# UNTESTED rebuild recipe: the exact copy of the env is the conda-pack archive deposited with the data.
options(repos = c(CRAN = 'https://cran.r-project.org'))
if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes')
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
pk <- read.csv(text = '
Package,Version,kind,source
bitops,1.0-9,cran,CRAN
caTools,1.18.3,cran,CRAN
deldir,2.0-4,cran,CRAN
dotCall64,1.2,cran,CRAN
fastDummies,1.7.5,cran,CRAN
fitdistrplus,1.2-2,cran,CRAN
ggridges,0.5.6,cran,CRAN
ggvenn,0.1.10,cran,CRAN
goftest,1.2-3,cran,CRAN
gplots,3.2.0,cran,CRAN
gtools,3.9.5,cran,CRAN
ica,1.0-3,cran,CRAN
lmtest,0.9-40,cran,CRAN
patchwork,1.3.0,cran,CRAN
pbapply,1.7-2,cran,CRAN
polyclip,1.10-7,cran,CRAN
RANN,2.6.2,cran,CRAN
RcppHNSW,0.6.0,cran,CRAN
ROCR,1.0-11,cran,CRAN
sctransform,0.4.1,cran,CRAN
Seurat,4.4.0,cran,https://satijalab.r-universe.dev
SeuratDisk,0.0.0.9021,github,mojaveazure/seurat-disk@877d4e18ab38c686f5db54f8cd290274ccdbe295
SeuratObject,4.1.4,cran,https://satijalab.r-universe.dev
sp,2.2-0,cran,CRAN
spam,2.11-1,cran,CRAN
spatstat.data,3.1-6,cran,CRAN
spatstat.explore,3.4-2,cran,CRAN
spatstat.geom,3.3-6,cran,CRAN
spatstat.random,3.3-3,cran,CRAN
spatstat.sparse,3.1-0,cran,CRAN
spatstat.univar,3.1-2,cran,CRAN
spatstat.utils,3.1-3,cran,CRAN
tensor,1.5,cran,CRAN
zoo,1.8-14,cran,CRAN')
bioc_version <- NULL
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
