# Audit receipt and precision verification

**Prototype; not yet carefully vetted.** Checked 20 September 2026 at 22:41 UTC. Bounded local verification only: saved reports, receipts, hashes and stated source rows were inspected. No auditor output was changed, no statistical analysis was rerun and no remote system was accessed.

The latest observed completed receipt is **e79a0475-81ab-4fa3-b94b-8f01bccd9734**, created 22:35:16 UTC. It explicitly supersedes **c1e0db4a-c32e-4fb4-8e55-8e2ba54eb40b**, created 22:33:21 UTC. Both refer to task `4704c686-0f68-4d94-b928-d8e5540222eb` and frozen manifest `4bd8c3520f2011c8ba503bf3623b9b0713fe96a7b6a17f22d4d356a7e4db36a0`. Completed means the audit was delivered, not publication approval.

All **41/41 files** listed in the auditor's `OUTPUT_SHA256.txt` match their recorded SHA-256 hashes. All **7/7 outputs** explicitly listed in the latest receipt match their hashes and byte sizes, including the hash manifest itself. `FIRST_ASSESSMENT.md` has the same hash in the old receipt, new receipt and current file, preserving the independent first record.

- Latest receipt SHA-256: `b9902a5049cd155fdb7b2b6a0feb43d672b834ffb527a66027c1c983b5548a42`.
- Output hash-list SHA-256: `1931cda6291b428bfcbfce15d39b4d20ff5d0ec345fa531b15405b5b7fe1b4c1`.

## Three requested precision checks

1. **V5 tile width: confirmed.** The saved v5 `quantitative/validation.txt` states an actual tile of **7.094 × 4.762 mm**, hence a half-tile height of 2.381 mm. The corrected audit report, design addendum and F13 agree. This verifies the recorded dimensions; it is not a new physical measurement of the printout.
2. **Dashed additions: count confirmed, developmental-stage wording corrected below.** The saved audit results and original union table agree that **32 of 33 additional pairs have phloem/xylem support**. The remaining pair has only endodermal support outside the three mismatched fer-4 contexts. However, the audit's phrase **“elongating-endodermis-only” is too narrow**: the original source rows include both elongation and WT maturation. This evidence note supersedes that phrase wherever it appears in the audit report, findings, design context or receipt. The auditor's files remain unchanged.
3. **Nomination dependencies: confirmed with the intended interpretation.** The saved audit results report **78 of 512 distinct GRN-nominated TFs** nominated only in the three inconsistent contexts, and **11 of 12 endodermal context nominations among the ten tested DOFs** supplied by those contexts. These are dependencies under the present pooled thresholds. They are **not upper bounds or predictions** of membership changes after correction, because recomputed pooled rank-difference thresholds can also change nominations elsewhere. The latest report, addendum and F03 state this correctly.

## Exact support for the remaining additional pair

**DAG1 (AT3G61850) → SAUR78 (AT1G72430)** has the following three source rows:

| Cell context | Genotype | Recorded edge/score count status |
|---|---|---|
| Elongation_Endodermis | WT | verified_outdegree_match |
| Elongation_Endodermis | fer4 | verified_outdegree_match |
| Maturation_Endodermis | WT | verified_outdegree_match |

The union row reports `n_current_seven = 0`, `n_phloem_xylem_map_issue = 0`, `n_endodermis_score_mismatch = 0`, `n_other_endodermis = 3`, and two distinct cell-context labels. Matching the count check does not establish comprehensive model validity.

**Correct short wording:** “Of 33 additional pairs, 32 have phloem/xylem support; the remaining DAG1–SAUR78 pair is supported in elongating endodermis of both genotypes and WT maturation endodermis, outside the three mismatch contexts.”

## Evidence locations

- [Latest completed receipt](../../../.living/tasks/fer-fig6-independent-audit-2026-09-20/receipts/e79a0475-81ab-4fa3-b94b-8f01bccd9734.json).
- [Audit report](../../../.living/tasks/fer-fig6-independent-audit-2026-09-20/outputs/claude-science/AUDIT_REPORT.md), [design addendum](../../../.living/tasks/fer-fig6-independent-audit-2026-09-20/outputs/claude-science/DESIGN_ADDENDUM_v5.md), [findings](../../../.living/tasks/fer-fig6-independent-audit-2026-09-20/outputs/claude-science/findings.csv) and [output hashes](../../../.living/tasks/fer-fig6-independent-audit-2026-09-20/outputs/claude-science/OUTPUT_SHA256.txt).
- [Saved audit results](../../../.living/tasks/fer-fig6-independent-audit-2026-09-20/outputs/claude-science/checks/results/audit_checks_results.json), keys `sets`, `exposure_to_inconsistent_contexts`, and `inventory_and_union`.
- [V5 recorded dimensions](../review_v5/quantitative/validation.txt).
- [Original all-context union](../review_v2/all_contexts/all_contexts_as_exported_union_edges.csv) and [exact context/genotype edge rows](../review_v2/all_contexts/all_contexts_as_exported_source_context_edges.csv), identified by source `AT3G61850` and target `AT1G72430`.
