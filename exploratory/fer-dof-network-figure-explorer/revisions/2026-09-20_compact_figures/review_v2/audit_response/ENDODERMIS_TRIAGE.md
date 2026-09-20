# Endodermal edge/score inconsistency: bounded triage

20 September 2026. Read-only inspection of the pinned public workflow, its source exports, the earlier direct saved-object check, and the currently installed Resnick CellOracle scoring code. No source was repaired, no graph was refitted, and no original/checkpoint hash was changed.

**Decision:** reconcile the three affected fer-4 endodermal graphs and scores before presenting their scores, ranks, or paired networks as final publication results. This is a specific internal-consistency problem, not absence of the underlying data and not evidence that the manuscript's experiments or the whole single-cell dataset are invalid. A new GRN inference run is not yet shown to be necessary; the smallest repair may be recomputation of graph scores from the correct existing filtered graph.

## The mismatch is not a weighted-degree convention

The pinned public workflow explicitly filters links at p < 0.01 and top 20,000 by absolute coefficient, calls `links.get_score()`, saves that object, then exports its `merged_score` and its same `filtered_links`.

- [Public run script, lines 280-293](https://github.com/tmnolan/feronia-cell-type-specific/blob/5649a9485de0d50d03c8135c5ff254bdaa8d0a1a/scripts_and_inputs/Analysis/CellOracle_GRN_analysis/scripts/20260317_feronia_salinity_time_course_WT_untreated_fer4_untreated_CellOracle_newIntegratedSeurat3_NoSuspiciousCluster_slurmArray_run_celloracle_grn.py#L280): filtering, scoring and saving the filtered object.
- Same script, lines 295-307: export the stored merged score table. Lines 314-325: export filtered edges by cluster. The edge export hard-codes 20,000 labels per cluster, but every endodermal cluster actually has 20,000 edges, and the mismatch also occurs inside the HDF5 object, so that labeling assumption does not explain this result.
- Currently installed CellOracle **0.20.0**, verified on Resnick at 2026-09-20 21:58 UTC: [legacy scoring source](installed_get_newtork_scores.R), lines 154-158, constructs a directed graph directly from the supplied source/target columns, without a second edge filter. Lines 89-105 calculate `degree_out = degree(g, mode="out")` and `degree_centrality_out = degree_out/(vcount(g)-1)`. Edge weights are assigned, but outgoing degree is an unweighted count.
- [Installed links-object source](installed_links_object.py), lines 215-233, passes `self.filtered_links` to the legacy scorer and reads its output into `merged_score`. The newer Python implementation, lines 512-535, uses the same degree/count definition, but it is not the method invoked by the public script.

The local source snapshots preserve original line numbers. [Snapshot metadata](installed_celloracle_source_snapshot.json) records current package version, remote paths, SHA-256 hashes and complete text. The historical March execution version was not independently recovered in this bounded check; current installed code is not automatically proof of the exact past environment. The public score columns are consistent with the legacy R scorer, including its cartography fields.

## Small independent reproduction

[check_score_identity.R](check_score_identity.R) rereads the pinned public endodermis edges and score CSVs. All eight edge lists contain 20,000 edges, no duplicate source-target pairs and no self-edges. Every stored score table also sums to 20,000 outgoing edges. The three affected contexts nevertheless describe different node universes and different gene-level outdegrees:

| fer-4 endodermal context | Nodes in filtered edge list | Genes in score table | Stored versus counted DAG1 outdegree | Stored versus counted DOF4.6 outdegree |
|---|---:|---:|---:|---:|
| Maturation | 1,847 | 1,296 | 81 versus 15 | 76 versus 15 |
| Proliferation domain | 1,696 | 1,549 | 20 versus 45 | 21 versus 205 |
| Transition domain | 1,505 | 1,457 | 30 versus 0 | 65 versus 247 |

The other five endodermal genotype-contexts have matching node counts and outgoing-degree comparisons. Exact tables are [graph summary](endodermis_graph_summary.csv), [all gene comparisons](endodermis_gene_degree_comparison.csv) and [ten-DOF comparisons](tested10_endodermis_degree_comparison.csv). Zero counted outgoing links is distinguished from an absent stored score.

The earlier direct Resnick HDF5 check read its `filtered_links` and `merged_score` in the same object and reproduced the same discrepancies: 777, 940 and 971 scored genes differ, with 225, 108 and 65 edge-source genes absent from the corresponding score tables. Evidence: `../../../2026-09-20_dof_subnetworks/inventory/check_endodermis_saved_object_remote.py`, lines 4-21, and `endodermis_saved_object_check.json`, lines 29-62. Those files are under the repository's `revisions/2026-09-20_dof_subnetworks/inventory/`; this is previously recorded direct-object evidence, not a second remote object read today.

Consequently, rounding, weighted degree, duplicate edges, and a CSV-only formatting problem do not explain the discrepancy. At least one of the saved components does not correspond to the other. Which is the correct intended graph remains unresolved.

## A concrete implementation hazard, not an established cause

The public run script is designed for a Slurm array (lines 169-171), uses one common output working directory (lines 164-167), and invokes deprecated `get_score()` (line 287). The installed implementation assigns each cluster a small integer ID and uses the shared relative folder `./network_analysis/` (installed_links_object.py, lines 215-236). The helper writes `network_analysis/<id>/linkList.csv` and launches its R calculation there (installed_use_r_scripts.py, lines 141-164); the results are loaded from those same numbered paths. It prints a failed child-process status without raising immediately (lines 166-172). The caller finally removes the whole shared folder.

**If multiple tissue jobs ran concurrently from the shared directory, they could collide in these temporary inputs/results.** Stale results after a failed subprocess are another possibility worth checking. This is an identifiable code hazard, but concurrency or a collision during the affected run has not been demonstrated. The exact endodermis launch script/log is needed: the currently committed `celltypes` dictionary has endodermis commented out (public script lines 151-161), so that file is not a complete frozen record of the historical endodermis invocation.

A bounded comparison of each affected score table against every other currently exported context found no exact gene-universe/outdegree match to another export. See [cross-context comparison](cross_context_score_identity.csv). This does not support a simple swap with one of the other currently published tables; it does not rule out an earlier/stale/intermediate result. No particular corruption mechanism is claimed as fact.

## Dependency and publication scope

The original nomination notebook reads the per-tissue score CSVs, ranks score values within genotype/context, keeps genes in the top 100 of at least one genotype, computes fer-4 minus WT rank differences, then uses pooled metric-specific 10%/90% tails. The selected TF lists use `degree_centrality_out`. Exact notebook locators are zero-based cell 29 lines 4-8 and 37-43; cell 7 lines 37-55, 75-112; cell 47 lines 2-7; cells 93 and 98. These cells are extracted with line numbers in [notebook evidence](original_nomination_notebook_relevant_cells.txt). The independently reproduced implementation is `revisions/2026-09-08_routeA_execution/q04_network/scripts/02_repair_network_annotations.R`, lines 129-153.

Therefore affected stored scores feed the endodermal nominations and contribute to thresholds pooled across contexts. If scores change, derived pooled TF sets, family enrichment, the 512-member GRN set, and overlaps involving that set require recalculation and a before/after comparison. It has not been shown that those memberships already changed. If the saved scores are correct and only the filtered edge lists are stale, the score-based nomination sets may remain unchanged.

The 84/87 passing exported contexts and the conservative seven-fit target overview can be retained as bounded existing-data results; passing an edge-count check is not comprehensive validation. The all-context overview should remain explicitly exploratory until reconciliation. Independent basal DE tables, the curated DOF-family annotation, promoter assays, EMSA, reporter imaging and growth experiments are not computationally invalidated by this graph/score mismatch. The separate phloem/xylem target-feature selection issue and lack of an endodermal YUC5 target in the model input remain separate concerns.

## Minimal repair gate to propose

1. Freeze the existing raw-link object, filtered-link object, exported scores/edges, input IDs, original launch script, software versions and logs. Record hashes and work in a new directory.
2. On the appropriate compute node, reconstruct the declared p/absolute-coefficient/top-20,000 filter from the existing raw `links_dict`, preserving tie handling, and compare exact source-target keys and coefficients against `filtered_links`. This checks graph stage/provenance without a new regression fit.
3. If the filtered graph is confirmed correct, recompute the required graph scores in a unique per-run/per-tissue directory, preserving the original scoring definitions. At minimum, require exact outgoing degree and degree-centrality agreement with the confirmed graph. Regenerate other dependent metrics needed by the retained ranking workflow, rather than mixing columns from different graphs. Do not silently switch legacy/current scoring implementations because non-degree metrics can differ.
4. If the filtered graph is stale but the raw links and run provenance are sound, reconstruct the correct filtered graph from the existing raw links, then score and export it consistently. Only if raw-model/input provenance cannot be reconciled is a new inference run justified.
5. Require object-to-export identity and edge/score agreement across the regenerated outputs; rerun nominations, enrichment and overlaps using the corrected sources, retaining a membership-change report and the untouched originals. Review an unaffected context as a method-equivalence control.

This repair gate has not been executed. It separates a likely graph-score reconstruction from the larger and currently unproven need to rerun inference from counts.
