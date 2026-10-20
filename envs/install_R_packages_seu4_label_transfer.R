# Re-install, at the recorded versions, the R packages in conda env 'seu4_label_transfer' that conda did not provide
# (packages installed by R itself, and conda packages that R later replaced with another version).
# Generated 2026-10-20 from installed.packages() of /resnick/groups/NolanLab/cchen/miniconda3/envs/seu4_label_transfer.
# Run AFTER recreating the conda env from envs/seu4_label_transfer.linux-64.explicit.txt, with that env's Rscript.
# Exits with status 1 and lists any package that is not at the recorded version.
options(repos = c(CRAN = 'https://cran.r-project.org'), timeout = 600)
lib <- .libPaths()[1]
if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes', INSTALL_opts = '--no-lock')
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', INSTALL_opts = '--no-lock')
pk <- read.csv(text = 'Package,Version,kind,source
bitops,1.0-9,cran,CRAN
caTools,1.18.3,cran,CRAN
cluster,2.1.8.1,cran,CRAN
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
zoo,1.8-14,cran,CRAN', stringsAsFactors = FALSE)
bioc_version <- NULL
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
