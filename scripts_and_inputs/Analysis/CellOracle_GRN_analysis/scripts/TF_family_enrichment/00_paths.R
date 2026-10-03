## Paths used by every R script in this folder. Run the scripts FROM THIS FOLDER (they `source("00_paths.R")`).
## Defaults follow this repository's layout (CellOracle_GRN_analysis/{inputs,outputs}); each can be overridden
## with an environment variable.
CO_DIR      <- Sys.getenv("CELLORACLE_DIR", file.path("..", "..", "outputs"))   # CellOracle outputs
WIDE_TABLE  <- file.path(CO_DIR,
  "20260317_feronia_salinity_WT_untreated_fer4_untreated_CellOracle_all_metrics_rank_change_wide_TF_GE.txt")  # ranks, Δrank, Label, Is_TF
SCORE_GLOB  <- file.path(CO_DIR, "*_gene_score_iGRN.csv")                       # per-network centrality scores
PLANTTFDB   <- Sys.getenv("PLANTTFDB_LIST", file.path("..", "..", "inputs", "Ath_TF_list.txt"))
              # https://planttfdb.gao-lab.org/download/TF_list/Ath_TF_list.txt.gz
KAY_TF      <- Sys.getenv("KAY_TF_LIST", file.path("..", "..", "inputs", "Kay_TF_thalemine_annotations.csv"))  # TF list used by CellOracle
OUT_DIR     <- Sys.getenv("OUT_DIR", file.path("..", "..", "outputs", "TF_family_enrichment"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
