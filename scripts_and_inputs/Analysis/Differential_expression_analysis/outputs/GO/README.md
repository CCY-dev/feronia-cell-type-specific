# GO enrichment tables

- `Upregulated_genes_GO_table.txt`, `Downregulated_genes_GO_table.txt`: the published tables (−log10 P, GO:BP terms × cell-type/stage clusters) written by the muscat notebook on 17 Mar 2026.
- `gprofiler_e113_eg59_p19_all_results_{up,down}.csv.gz`: complete g:Profiler results for every cluster's up- or down-regulated DEGs (from `../sig_de_fil.csv`), regenerated on 20 Oct 2026 with the archived g:Profiler release `e113_eg59_p19` (Ensembl 113 / Ensembl Genomes 59; organism `athaliana`, FDR correction, domain scope 'annotated', all results).

Reproducibility check (20 Oct 2026): rebuilding both published tables from these results with the notebook's rule (P ≤ 0.01, GO:BP, intersection ≥ 5, top 10 terms per cluster, all terms' values across clusters, missing = 1, −log10) gives the same 72 and 64 terms and 34 clusters, with a largest difference of 5e-14. The same release also reproduces the 68 terms of Extended Data Fig. 3 (top 15 terms per cluster for all DEGs). The live g:Profiler server (release e114 at the time of the check) gives slightly different values.

## The `significant` column

The notebook calls `gost(..., significant = FALSE)`, so g:Profiler returns every term with threshold 1, and its own `significant` flag is then TRUE for nearly every row. In these two files, `significant` has been reset to `p_value ≤ 0.05`, where `p_value` is g:Profiler's FDR-adjusted value (up: 3,218 of 35,738 rows; down: 3,721 of 39,151). No other value was changed. The published tables use their own rule (P ≤ 0.01 and intersection ≥ 5; above).

## Background

No custom background was passed, so the statistical domain is g:Profiler's annotated *A. thaliana* genes (`effective_domain_size` column), not `../../inputs/all_bg.csv`.
