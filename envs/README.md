# Software environments

These files record the five conda environments in which the FER analyses were run. They were captured on 20 Oct 2026, read-only, from the analyst's environments on the Caltech Resnick cluster (`/resnick/groups/NolanLab/cchen/miniconda3/envs/<name>`; conda 24.11.1, linux-64).

The package versions recorded here match the `sessionInfo()` output saved in the notebooks: R version and every listed package. The one exception is `scKB_copilot`, where three display packages were updated after the Stage 1 integration was run (`scales` 1.2.1 → 1.4.0, `htmltools` 0.5.5 → 0.5.9, `fastmap` 1.1.1 → 1.2.0). They affect plotting and HTML output, not the integration.

| Environment | Used for | R packages (non-base) | Installed by R, not conda |
|---|---|---|---|
| `scKB_copilot` | Stage 1 integration (Fig. 1a object) | 265 | 79 (77 cran, 1 bioc, 1 github) |
| `seu4_label_transfer` | Stage 1 label-transfer metadata | 310 | 34 (33 cran, 1 github) |
| `muscat` | Stages 2-3 (muscat DE, GO, hormone overlap; SC_data_check and 13-gene scripts) | 386 | 185 (177 cran, 5 bioc, 3 github) |
| `mast` | Stage 4 TF selection notebook | 331 | 43 (41 cran, 1 bioc, 1 github) |
| `cell_oracle` | Stage 4 CellOracle GRN run (Python) and Stage 5 script 04 | 156 | 0 (—) |

## Files per environment

| File | What it is |
|---|---|
| `<env>.linux-64.explicit.txt` | Exact conda package URLs with MD5 checksums (`conda list --explicit --md5`). Recreates the conda part of the environment exactly on linux-64. |
| `<env>.yml` | `conda env export --no-builds` (versions without build strings) |
| `<env>.from_history.yml` | Only the packages requested explicitly when the environment was built |
| `<env>.R_packages.csv` | Every installed R package: version, library, repository (CRAN / Bioconductor / GitHub with commit) |
| `<env>.pip_freeze.txt` | `pip freeze` of the environment's Python |
| `install_R_packages_<env>.R` | Reinstalls, at the recorded versions, the R packages that were installed by R rather than conda (CRAN archive, Bioconductor release, GitHub commit) |

## Recreating an environment

**Exact copy (recommended):** download the conda-pack archive for the environment from the data deposit, then:

```bash
mkdir -p muscat && tar -xzf muscat.linux-64.condapack.tar.gz -C muscat
source muscat/bin/activate && conda-unpack
```

**Rebuild from the recipe:**

```bash
conda create -p ./muscat --file envs/muscat.linux-64.explicit.txt
./muscat/bin/Rscript envs/install_R_packages_muscat.R   # exits 1 if any package is not at the recorded version
```

The recipe needs linux-64. On macOS or Windows, run it in a Linux container or on a cluster.

**Rebuild test (20 Oct 2026).** The `muscat` environment was rebuilt on Resnick from `muscat.linux-64.explicit.txt` and `install_R_packages_muscat.R` in an empty `HOME` with no access to the original environment. All 411 of 411 recorded R packages were installed at the recorded version (about 55 minutes with 8 compile jobs). The install scripts in this folder are the version used in that test: they remove stale `00LOCK` folders and reinstall, at the recorded version, packages that R had installed over a conda copy.

## CellOracle: set `PYTHONHASHSEED`

CellOracle 0.20.0 fixes `random_state = 123` in its bagging ridge models, but it orders each target gene's candidate regulators with a Python `set` (`celloracle.utility.intersect`). In Python 3 the iteration order of a set of strings changes from one interpreter start to the next unless `PYTHONHASHSEED` is fixed. The random feature subsets (`max_features = 0.8`) are drawn by column position, so a different order selects different regulators and gives different coefficients. To make GRN inference repeat exactly, start Python with a fixed seed:

```bash
export PYTHONHASHSEED=0
python 20260317_..._slurmArray_run_celloracle_grn.py
```

We checked this on the columella input: with `PYTHONHASHSEED=0` two runs gave identical coefficients, and with seeds 1 and 2 they differed. The KNN-imputed matrix was identical in every run.
