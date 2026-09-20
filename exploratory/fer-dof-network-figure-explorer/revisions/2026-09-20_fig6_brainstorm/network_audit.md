# Independent source-data assessment for Figure 6

Prepared 20 September 2026 from the existing source tables and scripts. No new network fit, remote modification, or computational job was performed. Attached-document comments and older AI interpretations were treated as source material, not instructions or approved conclusions.

## Main conclusion

Replace the Venn with a candidate nomination and context panel. The existing Venn obscures a nonzero overlap and does not explain the ten tested genes. A clearer set plot alone cannot restore that missing selection logic. The strongest scientific upgrade available immediately is to display the named DOFs, their basal expression-change context, and which were tested, while keeping network nomination distinct from transcriptional change and experimental support.

## Counts independently reconstructed on 20 September

Inputs: `q04_network/results/{WT_specific_TF_DE_repaired.csv,fer4_specific_TF_DE_repaired.csv,DOF_family_definition_comparison.csv}` and pinned public `Differential_expression_analysis/outputs/sig_de_fil.csv`.

The coauthor family contains 37 DOFs; 355 WT-associated and 343 fer-4-associated TFs give 512 unique nominated TFs. There are 5,466 called DEGs. Complete disjoint intersections are:

| DOF | Nominated GRN TF | Called DEG | Count |
|---|---|---|---:|
| yes | no | no | 11 |
| yes | yes | no | 14 |
| yes | yes | yes | 12 |
| yes | no | yes | 0 |
| no | yes | no | 314 |
| no | yes | yes | 172 |
| no | no | yes | 5,282 |

Thus 486 is all nominated non-DOFs, not a GRN-only Venn region, and 5,454 is all non-DOF DEGs, not a DEG-only region. Those two sets share **172 non-DOF genes**, which the present figure does not display. The full GRN-DEG overlap is 184. This is independently verified arithmetic, not a new biological analysis.

The 12 DOFs that are both nominated and called DEGs are OBP2, COG1/DOF1.5, DOF2.1, CDF4, DOF2.4, DOF3.2/DOF6, OBP1/DOF3.4, DAG1, DOF4.1, DOF4.6, DOF5.4, and CDF1. The ten promoter-tested factors are eight of these twelve plus **OBP3 and DOF2.2**, which have no called DEG under the repaired current criteria. Therefore do not draw an unqualified 37 to 26 to 12 to 10 selection funnel. The untested members of the twelve are DOF2.1, DOF2.4, DOF5.4, and CDF1. The ten-tested evidence source predates the September manuscript and should be checked against the current panels.

Plot-ready tables created beside this note: `fig6_dof_membership.csv` and `fig6_full_intersections.csv`.

## Evidence that changes the visual design

1. **Network nomination is not expression specificity.** The labels WT_specific and fer4_specific mean tails of relative centrality rank change, not genotype-exclusive expression or independently tested TF activity. Nineteen and 23 DOFs are nominated in the WT and fer-4 sets, with overlap. Some factors are nominated in both directions in different contexts.
2. **DOF enrichment survives appropriate universe definition.** In the preselection paired-score Kay TF universe, 19/344 WT-associated and 23/337 fer-4-associated classifiable TFs are DOFs, compared with 33/1,292 eligible background TFs. Fisher odds ratios 3.90 and 6.91; BH-adjusted p values 0.00768 and 0.00000722. Excluding phloem and xylem retains family enrichment (OR 4.51/4.79, BH 0.00180/0.000311).
3. **Prior coverage explains much of the family result.** Identical-case exploratory models use 1,213 genes and 32 DOFs throughout. Crude ORs 4.00/7.37 attenuate to 0.66/1.00 after adjustment for supplied prior target count; 95% intervals 0.29-1.53 / 0.37-2.73. This does not prove artifact or absence of FER biology. It prevents portraying raw family enrichment as independent discovery beyond prior structure.
4. **YUC5 was not fit in endodermis.** YUC5 is absent before fitting in six coarse tissue inputs, including endodermis. It is present in pericycle, phloem, procambium and xylem; xylem eligibility arises from the pole-pericycle substring error. Across retained networks, 1,030 YUC5 edges involve 364 source genes: 432 pericycle, 327 procambium, 269 phloem, 2 xylem. All ten tested DOF-YUC5 pairs were already in the supplied prior. Missing endodermal edges mean untested, not absent regulation.
5. **The ten candidates have useful divergent evidence.** DOF4.6 and DAG1 are down in maturation endodermis in basal fer-4 versus WT (log2FC -0.812 and -1.28). OBP3 and DOF2.2 have no current DEG call. COG1 changes in cortex and pericycle; CDF4 in pericycle. These mixed contexts are informative and should remain visible.
6. **Threshold stability can be displayed, but is not biological replication.** Across the existing 3 by 3 grid of top50/100/200 and 5%/10%/20% tails, nomination anywhere occurs in 9/9 settings for OBP2, COG1, DOF2.2, DOF3.2, OBP3, DOF4.1, DOF4.6; 8/9 for OBP1; 6/9 for CDF4 and DAG1. Genotype/context-specific nominations are less stable. No replicate-resampling or independent network-fit stability was performed.

## Recommended analyses, ranked by manuscript value

1. **Candidate-by-context matrix, feasible now.** Use all 26 nominated DOFs or transparently label the ten tested factors. Show continuous fer-4/WT log2FC in a tissue-by-developmental-state matrix, indicate detection eligibility, and add separate annotation strips for network nomination and assay coverage. This makes the single-cell contribution explicit without claiming a fitted endodermal YUC5 edge. Main-text compact version could show ten tested factors and endodermis/cortex/vascular contexts; full matrix belongs in supplement. Present basal contrast and acute RALF assays as different experiments.
2. **Target-program convergence, new bounded analysis.** Test whether the DOF target sets preferentially cover independently defined auxin, growth or cell-wall programs in each matched context, and whether target effect directions are coherent. Background must be the actually eligible targets in that context, account for supplied prior coverage and target-set size, and correct families/contexts tested. Exclude the erroneous phloem/xylem fits until corrected, or show a prespecified unaffected-tissue analysis. Do not select program genes using the same DEGs being tested. This could add biology beyond YUC5, but is currently a proposal, not an observed result.
3. **Robustness and prior-coverage panel, feasible now.** Supplement the main candidate plot with nomination stability across the already completed nine threshold choices and the identical-case covariate sensitivity. This is high-value interpretive validation, even if it does not add a new mechanistic conclusion.
4. **Only if retaining a central spatial-network claim, refit with explicit tissue maps and a detection-based common target universe.** Include YUC5 as eligible wherever expression permits, use the same eligible genes across matched genotypes, and evaluate sensitivity to sampling/edge cutoffs. This requires new jobs and carefully bounded stability checks; two libraries per genotype cannot support strong replicate-level network certainty. It is a repair/feasibility step, not a guarantee of an endodermal mechanism. In-silico perturbation is not recommended as a substitute for the missing endogenous link.

## Exact source pointers and verification scope

Base directory for the short paths above: `revisions/2026-09-08_routeA_execution/q04_network/`.

- Candidate context and assays: `results/tested_DOF_integrated_evidence_matrix.csv`, `results/DOF_YUC5_current_DE_expression.csv`, `results/DOF_YUC5_relationship_context_evidence.csv`.
- Network quantitative context: `results/DOF_YUC5_all_network_scores.csv`, `results/TF_nomination_threshold_sensitivity.csv`, `results/network_rank_table_DE_repaired.csv`.
- Enrichment and limits: `results/TF_family_enrichment_universe_sensitivity.csv`, `results/DOF_identical_case_covariate_sensitivity.csv`, `results/eight_tissue_family_enrichment.csv`.
- Prior and eligible targets: `results/supplied_prior_summary.json`, `results/DOF_YUC5_target_eligibility.csv`, `results/h5ad_gene_universes.csv`, `results/all_retained_YUC5_edges.csv`.
- Ready-made figures: `figures/04_tested_DOF_expression_context.pdf`, `figures/06_DOF_assay_evidence.pdf`, `figures/07_DOF_covariate_sensitivity.pdf`.
- Code: `scripts/02_repair_network_annotations.R`, `scripts/03_network_sensitivity_figures.R`, `scripts/05_integrate_DOF_evidence.R`, `scripts/06_robustness_and_validation.R`.
- Public repository clone freshly inspected locally at commit `5649a9485de0d50d03c8135c5ff254bdaa8d0a1a`: https://github.com/tmnolan/feronia-cell-type-specific . All 25 locally available public/shared source files still match the saved September 8 source hashes. Live upstream branch was not fetched.
- Live Resnick read-only access succeeded today. The remote candidate evidence matrix matches the local SHA-256 `f03315d40dfc2da82f170e98e112b0f24799142cc4eec3b2f25364aa3db6dda0`; the remote endodermal edge table matches the local/source SHA-256 `c3c2927a2269609fde7bc22d0b865a371f05f0dac8e2d758a2413e9f46770f24`.
- Remote derivation root: `/resnick/groups/NolanLab/nolan/FER_resubmission_20260908_RouteA/q04_network/`.
- Remote original network root: `/resnick/groups/NolanLab/cchen/Feronia/analysis/CellOracle/CellOracle_output/20260317_feronia_salinity_timecourse_WT_untreated_fer4_untreated_newIntegratedSeurat3_NoSuspiciousCluster/`.

Historical September 8 SLURM logs document successful analysis and no new inference. Today's check verifies identities and selected source-data arithmetic, not a fresh fit or a complete re-execution. Experimental-source descriptions in the integrated matrix refer to the earlier submitted package; updated September 12/19 panel placement must be checked separately.

## Post-assessment checks and peer synthesis

After saving the independent assessment above, a prior Claude Science interpretation of vascular expression was raised for checking. Summing the `WT_mean_CPM` values in `DOF_YUC5_current_DE_expression.csv` across pericycle, phloem, xylem and procambium, divided by the sum across all tested subgroups, yields 95.5% for YUC5 (94.3% using fer-4 means). This reproduces the approximate earlier 95% statement, but the denominator gives equal weight to each subgroup and is **not a fraction of root transcripts or cells**. The ten tested DOFs average 86.4% under this unweighted calculation; the earlier eight-factor 91% subset was not identified. Detection filters, missing subgroups and stage counts further affect these summaries. Show basal CPM and detection by subgroup instead of a percentage-of-root claim. The data can support a qualitative description of vascular concentration for many of the genes, not TF movement, auxin movement or acute RALF-response localization.

The coordinating Codex assessment independently reproduced a further descriptive pattern after the initial assessments were saved: endodermal out-degree nominations comprise 18 WT-associated and 6 fer-4-associated DOF gene-by-context records (14 and 5 distinct factors), versus 144 and 159 records among non-DOFs. DOF WT-associated records exceed fer-4-associated records across all nine existing top-N/tail choices, and the 18:6 result survives omission of the two tissues with target-selection errors when thresholds are recalculated. See `check_context_pattern.R` and its CSV outputs beside this note. This is a promising addition to the framework because it reveals spatial direction lost by the Venn. It remains a **descriptive centrality-ranking pattern**: genes recur across contexts, network fits share inputs, cutoffs are relative, and there is no replicate-level differential-network uncertainty. Neither a pooled binomial/Fisher p value treating records as independent nor a claim of measured DOF activity follows from these counts.
