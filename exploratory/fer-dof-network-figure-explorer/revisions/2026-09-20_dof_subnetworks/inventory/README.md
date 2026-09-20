# Existing DOF subnetwork inventory

Created 20 September 2026 from existing CellOracle edge exports and current basal differential-expression results. No networks were refitted, no new inference tests were run, and no source files on Resnick were changed.

## Main result

The ten Fig. 6 reporter-assayed DOFs have **41,726 retained source-target-context edges**, involving **2,647 distinct target genes** and **10,743 distinct DOF-target pairs** across the existing networks. These are filtered, prior-constrained fitted associations, not experimentally established targets. Every source and target checked is present in its actual tissue input. No duplicate source-target-cluster keys were found.

The ten canonical labels are OBP2, COG1 (source alias DOF1.5), DOF2.2, CDF4, DOF3.2 (DOF6), OBP1 (source alias DOF3.4), OBP3, DAG1, DOF4.1, and DOF4.6. Canonical labels join by AGI rather than symbol.

## Concrete biological leads

| Lead | Retained links and useful contexts | Interpretation boundary |
|---|---|---|
| YUC5 | 40 links from 8 assayed DOFs. Of these, 27 links from COG1, DAG1, DOF3.2, OBP2, and OBP3 are in pericycle or procambium, unaffected by the target-map issue. | All ten DOF-YUC5 possibilities were supplied in the aggregate prior. These fits provide context-specific retention, not independent discovery of binding. |
| YUC9 | COG1 and DOF4.6 in transition atrichoblast, both fer4. Coefficients +0.01390 and -0.01387, respectively. | Opposite model signs are not proof of causal activation/repression. An absent WT retained edge is not a significant genotype difference. |
| TAA1 | COG1 in distal and proximal lateral root cap, both genotypes, all positive coefficients. | Four fitted edges support a spatial lead, not a demonstrated biosynthetic mechanism. |
| Auxin response | Extensive IAA2, IAA4, IAA11, IAA17, IAA28 and other Aux/IAA targets, plus ARF14/16/2, GH3.3/.5/.6 and SAUR32/50/76/78. | Symbol-screen matches are descriptive leads. No pathway enrichment has been established. |
| Auxin transport | AUX1 links in cortex; LAX3 links occur only in the affected xylem fit. | No PIN-family symbol matches among these fitted targets. This is not evidence that PIN regulation is absent biologically. |
| Cell wall | EXPA11/17/7, EXPB1, XTH15/16/18/19/21/22/23 and pectin methylesterase family targets, among others. | Good candidate program for a manuscript-linked subnetwork, but use a fixed annotation set and eligibility background before claiming enrichment. |
| RALF/FER pathway components | RALF23 receives DOF2.2 and DOF3.2 edges in both genotypes in proliferation cortex, with DAG1 additionally retained in fer4. RALF1, FER, LRX3 and other RALF symbols also occur. | Feedback to the receptor-ligand system is a hypothesis. FER transcript changes in a fer mutant also reflect the perturbed locus itself and do not independently establish feedback. |

The biosynthesis screen finds YUC5, YUC9 and TAA1 among the source-provided target symbols. No other YUC/TAA/TAR-labelled target is retained. Symbol screening can miss alternate or absent annotations, so this statement is about recovered target labels, not an exhaustive curated pathway test.

## Recommended subnetwork settings

1. **Elongating phloem-pole pericycle:** COG1, DAG1 and OBP2 all retain positive YUC5 edges in WT and fer4. Use the same node positions and target set for both genotypes. Elongating xylem-pole pericycle and elongating procambium provide complementary contexts. Pole-pericycle contexts belong to the unaffected pericycle fit, not the affected phloem/xylem fits.
2. **Proliferation cortex:** combine RALF23 with a transparently selected auxin-response and wall-remodelling target program. This can connect ligand biology and growth outputs without forcing every DOF into an auxin-biosynthesis story.
3. **Transition atrichoblast and lateral root cap:** small, explicit YUC9 and TAA1 examples can show spatial specificity, preferably as companion panels or supporting material. Do not combine their edges into an unlabelled pooled network.

Avoid ten unrestricted networks: each DOF links to hundreds of distinct genes when contexts are pooled. A targeted common gene program and 2-3 fixed cell contexts will be easier to interpret. Alternatively, facet one chosen context into ten DOF-centered ego networks, with the target universe identical across genotypes.

Node colour may encode **current basal fer4 versus WT log2 fold change** in the same context. An unfilled or separate status should distinguish tested non-DE, below-detection, and unavailable results. Edge colour/sign should not use the same visual channel as node DE. Edge absence means failure to survive this fit/filter, not proven lack of regulation or differential-network significance.

## Newly detected export inconsistency

The full original CSV audit covers all **87 networks in 10 tissue fits**. In 84 networks, counting exported outgoing edges reproduces saved `degree_out` exactly. Three fer4 endodermal exports do not match their stored scores:

| Context | Scored source genes with unequal outdegree | Exported edge-source genes absent from score table |
|---|---:|---:|
| Maturation endodermis, fer4 | 777 | 225 |
| Proliferation endodermis, fer4 | 940 | 108 |
| Transition endodermis, fer4 | 971 | 65 |

All 30 corresponding comparisons for the ten assayed DOFs fail when zero-edge cases are included. Each export still contains 20,000 edges. This is not caused by duplicate keys or by this inventory's extraction. `full_public_export_degree_audit.csv` records the full check. A read-only check of the original filtered CellOracle HDF5 on Resnick reproduced exactly the same mismatches between its own `filtered_links` and `merged_score`; this is not merely a CSV export problem. The evidence is in `endodermis_saved_object_check.json` and its companion remote-reading script. It remains unresolved whether the saved edge set, the saved scores, or both require correction. Do not attach these three edge lists to the prior ranking claims or use their apparent genotype changes in a biological figure until the saved graph and score provenance are reconciled. WT endodermal and fer4 elongation exports pass the count check, but the paired figure is clearer using other tissues for now.

The pre-existing phloem/xylem issue is separate: substring tissue matching admitted targets from pole-pericycle contexts. Removing those two fits leaves eight target-map-unaffected tissues, but that set still contains the endodermal export inconsistency. `tested10_retained_edges_seven_unaffected_consistent_tissues.csv` excludes phloem, xylem and endodermis completely for conservative paired visualizations (25,617 edges; 1,988 targets). The raw all-ten file preserves both issue flags.

## Files and reproducibility

- `tested10_retained_edges.csv`: all ten assayed DOFs; canonical and source symbols, exact target IDs, tissue/cell state/genotype, signed coefficient, original p, actual input membership, and two independent issue flags.
- `tested10_summary.csv`: per-DOF target counts in all fits and in the eight target-map-unaffected fits.
- `tested10_context_target_counts_and_lists.csv`: exact target lists and count-versus-score checks per DOF and genotype-context.
- `tested10_paired_context_target_counts.csv`: matched WT/fer4 counts, with missing score contexts distinct from zero outgoing edges.
- `auxin_biosynthesis_exact_edges.csv`: all 46 exact YUC5/YUC9/TAA1 edge rows.
- `programme_lead_edges.csv`, `programme_target_summary.csv`, `programme_context_summary.csv`: symbol-based annotation leads. The regex dictionary is explicit; these are not GO results.
- `matching_basal_DE.csv`: 96,356 rows for all 2,649 unique DOF and fitted-target IDs recovered from the full current remote DE output, covering 37 cell-state contexts. The source stayed on Resnick.
- `matching_basal_DE_annotated.csv`: adds detection eligibility and current DE calls. Join `gene + cluster_id` to a network node AGI plus `cell_type`. Missing joins must remain unavailable.
- `extract_matching_de_remote.py`: read-only remote CSV subset extraction, no remote writes. `basal_DE_extraction_provenance.csv` records its source.
- `build_inventory.R`: reproducible inventory and local assertions. Run from project root with `/usr/local/bin/Rscript revisions/2026-09-20_dof_subnetworks/inventory/build_inventory.R .`.

The current DE rule is detection >=10% in either genotype, local adjusted p <=0.05, and absolute log2 fold change >log2(1.5). These results compare basal untreated fer4 and WT, not acute RALF induction.

## Coefficient and filtering provenance

The recovered public Python script calls `oracle.get_links(..., alpha=10)` and then `links.filter_links(p=0.01, weight="coef_abs", threshold_number=20000)`. The saved endodermal HDF5 metadata explicitly records `model_method="bagging_ridge"` and alpha 10. The script exports filtered edges, including `coef_mean`, `coef_abs`, and `p`. In the extracted ten-DOF table, `coef_abs` equals `abs(coef_mean)` exactly and every retained p is below 0.01 (maximum 0.009976515). These p values are model-edge filtering statistics, not replicate-level genotype comparison tests. No BH correction or genotype-edge difference test is asserted here. The original script's comment calling alpha a sparsity control should not be interpreted as a separate hard edge threshold.

Sources are the local frozen q04 outputs under `revisions/2026-09-08_routeA_execution/q04_network/results`, the pinned public-source clone under `revisions/2026-09-08_fresh_audit/work/feronia-cell-type-specific`, and the current full Resnick DE table named in the extraction provenance. `input_md5.csv`, `input_sha256.tsv`, `validation.csv` and `R_sessionInfo.txt` support reproduction.
