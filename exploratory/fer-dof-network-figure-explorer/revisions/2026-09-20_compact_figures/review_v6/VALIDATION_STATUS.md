# Current validation status

**Prototype; not yet carefully vetted.**

This version is for figure review and documented follow-up. Independent checks reproduce the calculations from the saved source tables, but two source-analysis problems remain unresolved. Audit completion is not scientific approval or permission to publish.

## What was checked

- Claude Science independently reproduced the set totals (37 DOFs, 512 GRN-TFs and 5,466 scDEGs), all seven exclusive intersections, 110 family enrichment tests, the top-ten selection, the 130-cell nomination matrix, the 41,726-row edge inventory and the 180-pair all-context union. Five frozen R scripts were rerun in a scratch environment; 14 checked data outputs were byte-identical.
- The new [membership reconstruction](../review_v5/provenance/reconstruct_memberships.R) independently rebuilds the identities and memberships from the recorded nomination and DE sources. All seven intersection counts and all 37 DOF memberships match. Display aliases are documented separately.
- Version 6 keeps both source and plotted matrix CSVs byte-identical to v5. Panels a and c are 110 mm wide, aligned at both boundaries. Forest text is at least 8 pt; matrix fractions and tissue labels are 8.5 pt, gene labels 9 pt, and the key 8 pt. Panel b has larger type, with a minimum of 8.28 pt.
- The wrapped companion-cell label was shifted 1 mm left; a 600 dpi crop confirms separation from Protoxylem. The unused missing-score key was removed. Evidence-bar starts are now computed cumulatively from the counts, without changing plotted values.
- Native vector PDFs, editable SVGs, rendered proofs and three colour-vision-deficiency simulations are retained in the panel directories. Layout checks do not resolve biological or analytical uncertainty.

## What still needs correction

1. **Endodermal graphs and stored degree scores disagree in three fer-4 contexts:** proliferation, transition and maturation. Of 87 exported contexts, 84 agree. The affected score sums still equal 20,000 but their node sets differ from the corresponding graphs. Which saved object is stale has not been established. First verify filtered links against existing raw links and the declared filter, then recompute scores from the confirmed graph using isolated temporary directories. Recalculate nominations, pooled thresholds, family enrichment and GRN-dependent overlaps. A new inference fit is not yet shown to be necessary. A shared temporary-directory collision is a hypothesis, not a confirmed historical cause.
2. **Separate phloem and xylem fits used a substring-based gene-feature selection rule:** 349 of 3,359 phloem features and 650 of 3,635 xylem features were admitted only through pole-pericycle DE contexts. This is a gene-feature issue, not inclusion of pericycle cells. Correct the tissue map and refit these models if corrected fits are to be presented. Removing edges after fitting is not equivalent.

Under the current thresholds, 78 of 512 GRN-TFs are nominated only in the three inconsistent contexts, and 11 of 12 endodermal context nominations among the ten assayed DOFs depend on them. These are direct dependencies, not predictions or upper bounds on corrected results. Recomputed pooled thresholds can also change other contexts.

The all-context graph retains 147 pairs supported by the seven included fits and 33 additional pairs. Of those 33, 32 have phloem/xylem support; the remaining pair, DAG1 to SAUR78, is supported by elongating endodermis in WT and fer-4 and by maturing endodermis in WT. Those three genotype-context sources are not among the mismatched fer-4 contexts. Line classes therefore cannot be read as a binary label for correct versus incorrect biology.

## Interpretation and access limits

These are basal WT/fer-4 networks, not RALF-treatment networks. Relative out-degree nominations are not measured TF activity. Unions across contexts do not show that all links coexist in one cell, and they do not establish intercellular transport. YUC5 was excluded before fitting in several tissues, including endodermis, so missing edges there are not negative evidence. The displayed families are selected auxin-associated groups, not a complete auxin pathway census.

Codex confirmed live Resnick and GitHub access on 20 September 2026. Exact remote paths, check times and the limits of those reads are in the [source availability record](../review_v2/SOURCE_AVAILABILITY.md). Large objects were not all reloaded for this layout revision. Claude Science's independent audit used the pinned public exports and frozen remote-object records; it did not obtain fresh Resnick access in that session.

## Evidence and next steps

- [Independent audit report](../../../.living/tasks/fer-fig6-independent-audit-2026-09-20/outputs/claude-science/AUDIT_REPORT.md)
- [Receipt verification and context-wording correction](AUDIT_RECEIPT_VERIFICATION.md)
- [Findings table](../../../.living/tasks/fer-fig6-independent-audit-2026-09-20/outputs/claude-science/findings.csv)
- [Endodermal source triage](../review_v2/audit_response/ENDODERMIS_TRIAGE.md)
- [Panel captions](CAPTIONS_AND_VALIDATION.md)
- [Current inputs, outputs and code](provenance/CURRENT_INPUT_OUTPUT_CODE.csv)

Reconcile the source issues, rerun dependent results, compare changed gene memberships and effects, then obtain scientific review before any publication release. Neither source repair nor refitting was performed in this figure-layout revision. Nothing has been published to GitLab or sent to Chia-Yu.
