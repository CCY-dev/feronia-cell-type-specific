# Validation status

**Prototype; not yet carefully vetted.** The independent audit was delivered and its receipt verified. Its source-analysis findings remain unresolved. Reproducing calculations from saved tables and inspecting the layout do not establish scientific validity or publication readiness.

## Unresolved source issues

- **Endodermis:** filtered graphs and stored outgoing-degree scores disagree in three fer-4 contexts: proliferation, transition and maturation. Verify the filtered graph against the saved raw links before deciding which scores or exports to reconstruct. Dependent nominations, pooled thresholds, family enrichment and GRN-dependent overlaps require reassessment after source correction.
- **Phloem/xylem:** substring matching admitted pole-pericycle DE-derived gene features to the separate tissue fits. This concerns gene-feature selection, not unintended inclusion of pericycle cells. Corrected target sets require controlled refitting; deleting edges afterward is not equivalent.

The current networks compare basal WT and fer-4. They are not RALF-treated networks, and relative out-degree ranks are not direct measurements of TF activity. A union across cell contexts does not establish simultaneous links in one cell or intercellular transport. Missing endodermal YUC5 links are not negative evidence because YUC5 was excluded from that model's input features.

## Evidence to read

1. [Detailed source validation status](revisions/2026-09-20_compact_figures/review_v6/VALIDATION_STATUS.md).
2. [Audit report](.living/tasks/fer-fig6-independent-audit-2026-09-20/outputs/claude-science/AUDIT_REPORT.md) and [findings table](.living/tasks/fer-fig6-independent-audit-2026-09-20/outputs/claude-science/findings.csv).
3. [Receipt verification and wording correction](revisions/2026-09-20_compact_figures/review_v6/AUDIT_RECEIPT_VERIFICATION.md). Read this with the unchanged audit report. Of 33 additional network pairs, 32 have phloem/xylem support. The remaining DAG1 to SAUR78 pair is supported by elongation endodermis in both genotypes and maturation endodermis in WT, outside the three mismatch contexts. The audit's shorter phrase “elongating-endodermis-only” is incomplete.
4. [Endodermal source triage](revisions/2026-09-20_compact_figures/review_v2/audit_response/ENDODERMIS_TRIAGE.md), [source availability](revisions/2026-09-20_compact_figures/review_v2/SOURCE_AVAILABILITY.md) and [spatial-mechanism evidence review](revisions/2026-09-20_compact_figures/review_v2/SPATIAL_MECHANISM_REVIEW.md).
5. [Membership reconstruction and comparisons](revisions/2026-09-20_compact_figures/review_v5/provenance/README.md), [current proof checks](revisions/2026-09-20_compact_figures/review_v7/PDF_VERIFICATION.json) and [panel D checks](revisions/2026-09-20_compact_figures/review_v7/networks/validation.txt).

The saved audit verified calculations and selected reruns. The v7 change reorganized panel D and enlarged nodes and labels while preserving its 180 pairs, ten DOFs and 37 targets. No source repair, new network fit or new biological validation was performed while packaging these files.

The original audit hash list includes historical scratch rasters that are deliberately omitted here. [Audit hash scope](provenance/AUDIT_HASH_SCOPE.csv) distinguishes included outputs from omitted rasters; all listed original output hashes were rechecked in the source checkout during assembly. Some preserved historical documents contain local links or status statements from their original date. Use this package's root navigation and the current sharing-status record for the current entry points.
