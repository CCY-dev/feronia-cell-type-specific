# GO enrichment tables

- `Upregulated_genes_GO_table.txt`, `Downregulated_genes_GO_table.txt`: the published tables (−log10 P, GO:BP terms × cell-type/stage clusters) written by the muscat notebook.
- `gprofiler_e113_eg59_p19_all_results_{up,down}.csv.gz`: complete g:Profiler results for every cluster's up- or down-regulated DEGs (from `../sig_de_fil.csv`) from the archived g:Profiler release `e113_eg59_p19` (Ensembl 113 / Ensembl Genomes 59; organism `athaliana`, FDR correction, domain scope 'annotated', all results).

The published tables follow the notebook's rule: P ≤ 0.01, GO:BP, intersection ≥ 5, top 10 terms per cluster, every selected term's values across all clusters, missing = 1, −log10 (72 up and 64 down terms, 34 clusters). Later g:Profiler releases give slightly different values.

## The `significant` column

The notebook calls `gost(..., significant = FALSE)`, so g:Profiler returns every term with threshold 1, and its own `significant` flag is then TRUE for nearly every row. In these two files, `significant` has been reset to `p_value ≤ 0.05`, where `p_value` is g:Profiler's FDR-adjusted value (up: 3,218 of 35,738 rows; down: 3,721 of 39,151). No other value was changed. The published tables use their own rule (P ≤ 0.01 and intersection ≥ 5; above).

## Background

No custom background was passed, so the statistical domain is g:Profiler's annotated *A. thaliana* genes (`effective_domain_size` column), not `../../inputs/all_bg.csv`.
