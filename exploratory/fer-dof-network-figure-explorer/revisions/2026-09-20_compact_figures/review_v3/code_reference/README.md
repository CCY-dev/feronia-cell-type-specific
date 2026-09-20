# Figure code reference snapshots

**Prototype; not yet carefully vetted.**

This is a code reference bundle for the six figure candidates. It is **not a self-contained portable pipeline**. Inputs, full model objects, package environments and most figure outputs are not bundled. Original relative and absolute path dependencies are preserved in the copied scripts and manifests. No script was run for this packaging step.

The forest and matrix are the original **review v2** code snapshots. They do not represent the separate review v3 plotting edits. Other panel scripts are their original frozen-source versions. The only code change in each copy is the requested prototype comment; shebangs are preserved. The scripts can write derived outputs if run, so these snapshots should be reviewed before execution in a separately configured working directory.

| Panel | Plotting code snapshot |
|---|---|
| TF family enrichment forest, original v2 | [make_review_v2_forest.R](revisions/2026-09-20_compact_figures/review_v2/quantitative/make_review_v2_forest.R) |
| DOF conditional overlap | [make_dof_evidence.R](revisions/2026-09-20_compact_figures/review_v2/overlap/make_dof_evidence.R) |
| Cell-type nomination matrix, original v2 | [make_review_v2_matrix.R](revisions/2026-09-20_compact_figures/review_v2/quantitative/make_review_v2_matrix.R) |
| All exported contexts network | [make_all_context_overview.R](revisions/2026-09-20_compact_figures/review_v2/all_contexts/make_all_context_overview.R) |
| Elongating vascular contexts network | [make_overviews.R](revisions/2026-09-20_compact_figures/networks/make_overviews.R) |
| Exclusive-intersection UpSet | [make_compact_upset.R](revisions/2026-09-20_compact_figures/upset/make_compact_upset.R) |

## What is included

Thirty scripts preserve their original repository-relative hierarchy beneath this directory: the six panel generators, SVG grouping helpers, family/context aggregation and prior figure preparation, the ten-DOF edge inventory, and the Q04 annotation repair, family enrichment, sensitivity, evidence integration, provenance and Slurm launch scripts. The original public inference workflow and full data objects remain at their existing sources; this bundle contains the analysis and plotting references used for these prototypes.

- [Panel code links](panel_code_links.csv): six rows with primary copied-code paths, file URIs for PDF hyperlink construction, hashes, upstream steps and corresponding original figure locations.
- [Source/copy manifest](source_copy_manifest.csv): original paths, copied paths, SHA-256 hashes, roles, reference execution environments and run locations. Commands refer to the original checkout/HPC layout and are documentation, not portable commands for this bundle.
- [Dependency manifest](dependency_manifest.csv): 101 declared source/parameter/generated-artifact dependencies, with hashes where the named local input is a file. Dependencies not supplied by another copied script are not bundled.
- [Literal source excerpts](dependency_source_excerpts.csv): 304 original code locations exposing read calls, path construction, ROOT parameters and Slurm settings. This supplements the curated dependency list and is not an exhaustive dependency resolver.
- [Bundle checks](bundle_validation.json): confirms all original files were unchanged and every code copy differs only by insertion of the prototype header.

## Execution locations and unresolved provenance

The compact plotting scripts were run locally in `/Users/tmn23/projects/fer-sc-manuscript-2026`; several use native macOS Quartz PDF output. The Q04 statistical/provenance workflow used the existing Resnick environments and `/resnick/groups/NolanLab/nolan/FER_resubmission_20260908_RouteA/q04_network`, with full inputs under the original laboratory directories. Exact source expressions and reference commands remain in the manifests. No data or fit is repaired by copying these scripts.

The graph/score and target-selection limitations remain. See the existing [endodermis triage report](../../review_v2/audit_response/ENDODERMIS_TRIAGE.md) and [source availability report](../../review_v2/SOURCE_AVAILABILITY.md). Those reports are linked in place and were not copied or edited. The triage report is outside this code-only archive when distributed separately.

The original paper repository is [feronia-cell-type-specific](https://github.com/tmnolan/feronia-cell-type-specific/tree/5649a9485de0d50d03c8135c5ff254bdaa8d0a1a). This local reference bundle has not been published or pushed.
