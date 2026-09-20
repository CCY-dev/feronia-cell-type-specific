# Underlying FER data and source availability

20 September 2026. Read-only source check; no source repair, model fit, remote write or modification of earlier output/checkpoint hashes.

**Yes, the underlying data and analysis code are available.** The limitation in these figure drafts is reconciliation of particular saved network results, not absence of the source data. Current Resnick access and the GitHub repository were verified again during this check.

## Current access, directly verified

- Read-only SSH to the configured `resnick` host succeeded on `login3.cm.cluster` at **2026-09-20 21:41:30 UTC**. The paths below exist and are readable.
- `git ls-remote https://github.com/tmnolan/feronia-cell-type-specific.git HEAD` succeeded. Current public HEAD is `5649a9485de0d50d03c8135c5ff254bdaa8d0a1a`, identical to the local pinned clone used for this audit.
- Public source: [paper repository](https://github.com/tmnolan/feronia-cell-type-specific), [pinned CellOracle code, inputs and exported tables](https://github.com/tmnolan/feronia-cell-type-specific/tree/5649a9485de0d50d03c8135c5ff254bdaa8d0a1a/scripts_and_inputs/Analysis/CellOracle_GRN_analysis), [pinned differential-expression code and tables](https://github.com/tmnolan/feronia-cell-type-specific/tree/5649a9485de0d50d03c8135c5ff254bdaa8d0a1a/scripts_and_inputs/Analysis/Differential_expression_analysis). Repository/commit access was verified by Git; the linked subdirectories were verified in that commit's local tree. The web page fetcher returned a cache miss, which was not a repository-access failure.

| Available source on Resnick | Current read-only check |
|---|---|
| `/resnick/groups/NolanLab/shared_resources/annotated_objects/rds/integrated/Feronia_salinity_timecourse_4S_20251226_seu3_sc291suspiciousClusterRemoved.rds` | 19,549,746,309 bytes; source integrated Seurat object |
| `/resnick/groups/NolanLab/cchen/Feronia/analysis/CellOracle/CellOracle_input/20250509_4S_feronia_salinity_time_course_untreated_celloracle_input.h5ad` | 5,384,390,704 bytes; CellOracle input object |
| `/resnick/groups/NolanLab/cchen/Feronia/analysis/CellOracle/CellOracle_output/20260317_feronia_salinity_timecourse_WT_untreated_fer4_untreated_newIntegratedSeurat3_NoSuspiciousCluster/endodermis_rawRNA_p01_20k.links.celloracle.links` | 514,497,638 bytes; saved filtered endodermis model object |
| `/resnick/groups/NolanLab/cchen/Feronia/analysis/muscat_result/20260317_salinity_timecourse_useDateT_noSuspiciousClusters/group_timecelltypeannoLi/fer_4-WT/sig_de.csv` | 165,924,398 bytes; full DE table, header read successfully, includes gene, cluster_id, logFC, p_val, p_adj.loc, p_adj.glb and per-library/combined detection fields |

The same CellOracle output directory currently contains **all ten tissues'** `*_rawRNA.links.celloracle.links`, `*_rawRNA_p01_20k.links.celloracle.links` and `*_rawRNA_celloracle.h5ad` files: atrichoblast, columella, cortex, endodermis, lrc, pericycle, phloem, procambium, trichoblast and xylem. Thus raw saved links, filtered saved links and tissue-level inputs are present. These are model/input objects, not a claim that raw sequencing FASTQs were checked today.

This refresh checked metadata and a DE header. It did not load the 19.5 GB source object or recompute its checksum. Prior source-count reconstruction, content checks and hashes remain documented in the existing Route A package and were not redone or replaced.

## Existing local evidence and reusable data

- Public snapshot: `revisions/2026-09-08_fresh_audit/work/feronia-cell-type-specific`, verified HEAD as above. Its CellOracle outputs contain ten edge CSVs, ten score CSVs and the original ranking/nomination tables.
- Full prior provenance: `revisions/2026-09-08_routeA_execution/q04_network/results/source_file_manifest.csv`, `source_sha256_manifest.csv` and `../provenance/public_source_hash_comparison.tsv`. The recorded comparison found all 25 comparable public/shared files byte-identical. This is a prior verification, not a new checksum of every remote file today.
- Source-count reconstruction: `revisions/2026-09-08_routeA_execution/q02_specificity/`; source RDS SHA-256 in `provenance/source_RDS.sha256`. Existing logs report exact reconstruction of 718,502 tested gene-context rows and 349,972 detection-eligible rows from the source counts. Those results are historical recorded analyses, not newly rerun here.
- Current ten-DOF all-fit inventory: `revisions/2026-09-20_dof_subnetworks/inventory/tested10_retained_edges.csv`, 41,726 retained rows and explicit tissue-map and edge/score flags.
- Current matched basal-DE subset: `revisions/2026-09-20_dof_subnetworks/inventory/matching_basal_DE_annotated.csv`, with the read-only extraction source in `basal_DE_extraction_provenance.csv`. This subset contains 96,356 rows for 2,649 queried genes across 37 DE contexts.

## What specifically needs reconciliation

**Endodermis saved edge/score inconsistency.** In three fer-4 contexts, maturation, proliferation and transition, outgoing-edge counts from the saved filtered links disagree with saved `degree_out` scores. This was first observed in the CSV exports and then reproduced inside the original Resnick HDF5 object's own `filtered_links` and `merged_score`, recorded in `revisions/2026-09-20_dof_subnetworks/inventory/endodermis_saved_object_check.json` and its reading script. The corresponding unequal scored-source counts are 777, 940 and 971; 225, 108 and 65 edge sources are absent from the respective score tables. The saved source exists. Which component is stale or requires reconstruction is unresolved. The other five endodermal genotype-contexts pass this count comparison; passing it is not comprehensive model validation.

**Phloem/xylem target-feature selection.** Original substring matching admitted DE-derived target features from phloem-pole/xylem-pole pericycle into the separate phloem/xylem fits. This concerns selected gene features, not inadvertent inclusion of pericycle cells. An explicit tissue map removes 349 of 3,359 phloem-input genes and 650 of 3,635 xylem-input genes. This is documented against actual saved H5AD inputs in `revisions/2026-09-08_routeA_execution/q04_network/README.md` and its source tables. The data are available for a future controlled reconstruction; deleting edges afterward is not equivalent to refitting with the corrected target set.

The compact seven-fit overview conservatively excludes all endodermis plus phloem and xylem. Its scope is a transparent plotting choice while these issues are unresolved, not a statement that data from those tissues are missing.

## Consequences for nominations, family enrichment and overlap counts

The current family plots and overlap diagrams preserve the previously reconstructed all-context score-based nominations. Replacing only the legacy DE annotation previously left network scores and nominations unchanged. That annotation repair is distinct from reconciling the graph/score mismatch or correcting target-feature selection.

A controlled graph/score reconciliation or corrected phloem/xylem fit **could** change scores, ranking thresholds, nominated TF membership, eligible family denominators, enrichment estimates and overlaps involving the GRN set. It has not been shown that they already changed. If the stored scores prove correct and only some saved edge lists are stale, score-based nominations might remain unchanged. Existing eight-fit sensitivity showed DOF enrichment persists after omission of phloem/xylem; it did not resolve the later endodermis finding. The 37-member DOF annotation and 5,466-gene DE set do not inherently change merely because a network is reconciled, whereas the 512-member GRN set and overlaps involving it must be checked after any source correction.

## Requested all-context overview

`all_contexts/` now provides a separate **exploratory as-exported union** using the same ten DOFs, AGI symbol normalization and fixed auxin-family rule as the compact seven-fit overview. All 44 available context labels and 87 genotype-context exports were considered; selected edges occur in 43 context labels and 83 genotype-contexts. The omitted context under this particular target filter is elongating cortex, not a missing network file.

It contains 935 source rows, 180 distinct DOF-target pairs and 37 targets across all ten tissue fits. Of these, 147 pairs are also supported in the seven included fits and 33 occur only in the previously excluded fits. Of the 33 additions, 32 have phloem/xylem support and one has support in another conservatively excluded endodermal context; none is added solely by one of the three endodermal mismatch rows. There are three additional target labels, IAA6, IAA18 and IAA29. All original flags and exact context/genotype provenance are retained. These support classes do not identify which individual phloem/xylem coefficient would survive a corrected refit.

Suggested one-line footer: **Existing Resnick/GitHub data; as-exported network summary, with endodermal edge/score and phloem/xylem target-selection reconciliation pending.**
