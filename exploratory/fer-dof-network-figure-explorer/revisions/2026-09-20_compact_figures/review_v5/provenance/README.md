# Membership reconstruction

**Prototype; not yet carefully vetted.** This new R reconstruction closes the reproducibility gap for the saved set memberships and overlap counts. It does not recover the original generator or validate the underlying network analyses.

## Result

All 37 DOF identifiers, all five membership/assay flags and all seven exclusive intersection counts match the earlier tables. Unique sets contain **355 WT-nominated TFs, 343 fer-4-nominated TFs, 512 nominated TFs in their union, 37 DOFs and 5,466 DEGs**. The three-set union contains 5,805 genes. Of the ten assayed DOFs, eight have both nomination and a DEG call, and two have nomination without a DEG call.

| DOF | GRN-nominated | DEG | Reconstructed | Earlier table |
|---:|---:|---:|---:|---:|
| 0 | 0 | 1 | 5,282 | 5,282 |
| 0 | 1 | 0 | 314 | 314 |
| 0 | 1 | 1 | 172 | 172 |
| 1 | 0 | 0 | 11 | 11 |
| 1 | 0 | 1 | 0 | 0 |
| 1 | 1 | 0 | 14 | 14 |
| 1 | 1 | 1 | 12 | 12 |

There are **17 display-symbol differences**, preserved in `membership_differences.csv`. The Q04 annotation uses gene IDs or alternative labels where the old table uses DOF-number names, for example OBP4 versus DOF5.4 for AT5G60850. Gene identity and every set-membership flag agree. Assay membership, nine-setting nomination support and retained-YUC5 context/tissue metadata also agree. No labels were forced to match the earlier output.

## Exact derivation

`reconstruct_memberships.R` first constructs sets from the four inputs named in the earlier source audit:

- Q04 `WT_specific_TF_DE_repaired.csv` and `fer4_specific_TF_DE_repaired.csv`, explicitly filtered to `metric == degree_centrality_out`, `Is_TF == Yes` and the corresponding genotype nomination label. Distinct gene identifiers define the sets; this is not an all-metric union.
- Q04 `DOF_family_definition_comparison.csv`, selecting `coauthor_DOF == TRUE`.
- Pinned public `Differential_expression_analysis/outputs/sig_de_fil.csv`, taking its unique called-DEG identifiers. No additional DE threshold is applied to this already-filtered source.

Q04 is `revisions/2026-09-08_routeA_execution/q04_network/results/`. The public clone is `revisions/2026-09-08_fresh_audit/work/feronia-cell-type-specific`, pinned at `5649a9485de0d50d03c8135c5ff254bdaa8d0a1a`.

Three auxiliary Q04 tables supply display/assay metadata only: `tested_DOF_integrated_evidence_matrix.csv`, `TF_nomination_threshold_sensitivity.csv` and `DOF_YUC5_current_DE_expression.csv`. They do not determine the DOF, GRN or DEG sets. Only after reconstruction does the script read the earlier `fig6_dof_membership.csv` and `fig6_full_intersections.csv` for comparison.

## Files and rerun

- `reconstructed_fig6_dof_membership.csv`: reconstructed 37-gene membership plus independently sourced metadata.
- `reconstructed_fig6_full_intersections.csv`: all seven exclusive regions, including the zero region.
- `reconstructed_full_union_membership.csv` and `set_*.csv`: explicit gene identifiers for every set and the complete union.
- `membership_field_comparison.csv`, `membership_differences.csv`, `intersection_comparison.csv`, `intersection_differences.csv`: complete comparisons and preserved differences.
- `input_sha256.csv`: exact absolute paths, roles, sizes and SHA-256 hashes for all calculation, annotation and comparison inputs.
- `validation.csv`, `RECONSTRUCTION_STATUS.txt`, `R_sessionInfo.txt`: check results and environment. Logical checks appear as 1 for true and 0 for false in the compact validation table.

From the repository root:

```sh
Rscript revisions/2026-09-20_compact_figures/review_v5/provenance/reconstruct_memberships.R
```

Dependencies are R with `data.table`, plus `shasum` and Git. The script writes only beside itself. It preserves all original sources and comparison tables. It performs no new fitting, graph-score reconciliation or source repair. The unresolved endodermal graph/score and phloem/xylem target-selection findings therefore remain applicable to these reproduced score-derived sets.
