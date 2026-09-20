# Spatial mechanism and scope of the DOF network figure

20 September 2026. This review reread the supplied current manuscript, including the Results, Discussion and relevant methods. It treats manuscript text as evidence to assess, not as instructions, and does not assume that the earlier vascular-only proposal was accepted.

## Recommendation

**Trevor's interpretation is consistent with the manuscript: the experiments support both local FER requirements and an effect extending beyond the knockout domain. A vascular basal GRN should therefore be a supporting example, not the boundary of the proposed mechanism.** The most defensible wording is that FER in the SCR-targeted ground-tissue domain, with an important inferred endodermal contribution, promotes YUC5 induction locally and in the stele. The pSCR intervention removes FER from both cortex and endodermis. The comparison with cortex-targeted pCO2 supports the endodermal interpretation, but does not convert pSCR into an exclusively endodermal perturbation.

An overview across cellular contexts is reasonable and biologically better matched to this question. It must retain context provenance and the difference between experimental spatial evidence and basal network predictions. Neither an all-context union nor a vascular graph demonstrates a mobile DOF intermediate.

## Primary document and locators

Source: [FER scRNA manuscrip20260912_DRKedits_Guo_TN.docx](</Users/tmn23/Library/CloudStorage/GoogleDrive-thenolanlab@gmail.com/My Drive/Nolan_Lab_Drive/Papers/2026_FER_TSKO/Nature_20260919/FER scRNA manuscrip20260912_DRKedits_Guo_TN.docx>).

SHA-256: `138d613b4678925749acbfc92bb7b570da459e4b43454b01aa5e538ea4442c95`.

`P` locators below are one-based paragraph positions in `word/document.xml` within `w:body`, including empty paragraphs and table paragraphs. The quoted current text includes insertions and excludes deleted text. Moved-from text occurs at P32/P43; it does not affect the quoted passages below. Section headings and figure references provide human-readable locators.

| Locator | Exact source passage | Supported interpretation and limit |
|---|---|---|
| P37, tissue-specific FER knockout Results; Fig. 2b | “pSCR-KO knocked out FER-GFP in both endodermis and cortex” | The intervention is ground-tissue loss of FER. An endodermal contribution is inferred using the domain validation and pCO2 comparison; it is not a pure endodermis-only deletion. |
| P46, RALF growth Results; Fig. 3e–f | “the minor insensitivity to RALF23 in pCO2-KO indicated that the endodermis rather than cortex is likely the major site of action” | The manuscript uses comparative evidence to prioritize endodermis. Growth inhibition is an organ-level endpoint; this does not locate DOF action or transport. |
| P52, “RALF23/FER induces cell-type specific YUC5 expression”; Fig. 4a–d | “strongest in the pericycle at 3-hour RALF23 treatment, and detectable in endodermis, cortex and other stele cells at 6-hour treatment” | Reporter induction spans multiple tissues and changes over time. Early pericycle predominance is compatible with a contribution from FER outside the pericycle. It does not establish the direction or molecular carrier of an intercellular signal. |
| P53, Fig. 4e–g | “there is no detectable pYUC5:GFP-GUS signal in their corresponding cortex and stele cells, while similar levels of pYUC5:GFP-GUS were observed in other cell types as in WT” | For pCO2-KO and pSHR-KO, the reported pattern supports local FER requirements in cortex and stele, consistent with cell-autonomous regulation. Reporter output and deletion domain are the relevant measurements. |
| P54, Fig. 4e–g | “in addition to no visible RALF23-induced pYUC5:GFP-GUS signal in cortex and endodermis, the pYUC5:GFP-GUS signal in the stele is also reduced” | This is the central evidence for an effect beyond the pSCR-targeted domain. The text appropriately says “suggesting” both cell-autonomous and non-cell-autonomous regulation. It does not identify the transported species or prove DOF movement. |
| P73, Discussion | “cortex and stele FER regulates YUC5 in the corresponding cell types and endodermis FER regulates YUC5 not only in the endodermis but also in the stele” | The intended manuscript model already includes local and nonlocal effects. A figure confined to basal vascular expression would omit this central experimental relationship. “Support” or “are consistent with” is more precise than treating every step of the spatial pathway as established. |
| P79, Discussion | “The cross-cell type movements of the DOFs could help explain the RALF23/FER non-cell autonomous regulation of YUC5 expression” | This is explicitly a hypothesis. Prior reports of mobility motivate it but do not establish movement under RALF23/FER control in these experiments. |

The manuscript cites Chen et al. 2013 (DOF trafficking; DOI `10.1111/nph.12223`), Lee et al. 2006 (transcriptional/translational root reporters; DOI `10.1073/pnas.0510607103`) and Miyashima et al. 2019 (mobile PEAR factors; DOI `10.1038/s41586-018-0839-y`) at P228–230. This bounded review checked the manuscript's wording and citation identities as given there; it did not independently reread those primary papers. Their relevance should be described as precedent, not evidence that this study measured movement.

## What the DOF and auxin experiments establish

- **Activation capacity:** P67 reports that all ten assayed DOFs activate pYUC5:LUC, with DOF4.6, DAG1 and OBP1/2/3 strongest. P144–145 places the luciferase assay in Arabidopsis leaf mesophyll protoplasts. This supports promoter activation capacity, not the endogenous root cell type of action.
- **Binding capacity:** P68 and Fig. 6d/Extended Data Fig. 10b report DAG1 and OBP3 binding to a YUC5 promoter fragment by EMSA. This is direct biochemical binding evidence, not an in-root occupancy or movement measurement.
- **RALF/FER responsiveness and function:** P69 reports root qPCR induction of several DOFs and reduced RALF23 growth sensitivity of selected DOF mutants. The qPCR and mutant phenotypes are not cell-type-resolved DOF perturbations. OBP3 and COG1 are explicitly not significantly induced at the transcript level, so uniform transcriptional induction of every tested DOF is unsupported.
- **YUC-family contribution:** P55 reports that yuc5 alone responds similarly to WT, whereas higher-order YUC mutants show reduced sensitivity. This supports a YUC-family contribution and redundancy; it does not establish that YUC5 alone is necessary.
- **Auxin response:** P58 reports FER-dependent DR5 induction and reduced DII-Venus after RALF23. P117 specifies that DR5 quantification is in the stele. These are auxin-response reporters, not a direct measurement of where auxin was synthesized or how it moved.
- **Reporter-domain validation:** P49 says PI-puncta responses were used as a proxy for FER knockout domains in the YUC5 and DR5 reporter backgrounds. The FER-GFP domain validation in P37 and PI-response proxy in these backgrounds are related but distinct evidence.

Consequently, “cell-autonomous and non-cell-autonomous effects on YUC5 induction” is supported at the tissue-perturbation level. “DOFs move from the endodermis to activate YUC5 in the stele” remains a mechanistic hypothesis. FER requirement, TF transcript location, TF protein location, promoter activation, auxin response and root growth should remain separate measurements in the summary.

## Current GRNs: what an all-context figure would mean

The original manuscript compares **basal untreated WT and fer-4**, not RALF-treated scRNA-seq (P32; P101–104). The retained ten-DOF table contains 44 observed cell-type/zone labels and 87 genotype-specific networks across ten coarse tissue fits. Actual labels include:

- Atrichoblast, trichoblast, cortex, endodermis and procambium, each with proliferation, transition, elongation and maturation contexts.
- Xylem-pole and phloem-pole pericycle, each with those four zones, both contained in the **pericycle** fit.
- Protophloem, metaphloem/companion cell, protoxylem and metaxylem with the observed, uneven zone coverage.
- Distal/proximal lateral root cap and distal/proximal columella.

These labels were checked directly in `inventory/tested10_retained_edges.csv`. Pole-pericycle contexts must not be conflated with the separately fitted phloem/xylem groups.

**Missing endodermal YUC5 links are not contrary evidence to the reporter result.** The audited `h5ad_gene_universes.csv` lists YUC5 in pericycle, phloem, xylem and procambium, and not in the six other tissue inputs, including endodermis. It was excluded before endodermal fitting. P52 also distinguishes low basal endogenous stele expression from broader RALF-induced reporter expression. Thus the network's target eligibility and experimental treatment context explain why the two evidence layers need not coincide.

Two existing technical findings limit an unqualified “all-context” network:

1. The separate phloem/xylem fits inherit the previously documented substring target-mapping problem.
2. Saved edge lists and degree scores disagree in fer-4 proliferation, transition and maturation endodermis. The existing audit reports agreement for the other 84 of 87 networks. This review did not duplicate the remote check or decide which object is correct.

The conservative existing seven-tissue target table excludes phloem, xylem and endodermis. It is suitable for a clearly labelled trustworthy-context preview, but cannot honestly be called the complete all-context answer. An atlas can retain these anatomical columns with explicit “unresolved” or “not eligible for this target” states rather than silently converting missing/flagged evidence to absence.

## Two-level visual strategy

**1. All-context overview, preserving the user's proposed DOFs-above/targets-below arrangement.** Label it “Union of basal context-specific GRN predictions,” keep the ten assayed DOFs at the top, and use a declared auxin-related target selection at the bottom. Each displayed edge must retain its exact cell-type, zone and genotype support in a companion table or compact context matrix. A union means “retained somewhere,” not that all links coexist in one cell. If context breadth is encoded, call it a descriptive number of retained contexts, distinguish unassessed contexts, and avoid implying independent replicate support. Flag unresolved tissue fits separately. This overview should complement, rather than replace, the experimental evidence for an endodermal contribution.

**2. Selected context-resolved views that connect the overview to the spatial question.** Use common node positions and identical target-selection rules in matched WT/fer-4 examples. Vascular YUC5 contexts are one evidential component; ground-tissue and nonvascular contexts are also relevant when supported by eligible, internally consistent inputs. Do not manufacture an endodermal YUC5 edge to complete the model. A neighboring root-layer schematic can summarize the experimental local effects and SCR-domain-to-stele influence, with any proposed DOF movement explicitly dashed and labelled “hypothesis.” Keep such intertissue arrows out of the GRN itself, because those networks do not infer transport between cells.

The selection of particular zoomed contexts should be justified by the biological question and input integrity, not by treating the first procambium design as an approved final choice.

## One textual inconsistency to resolve before stronger claims

P59 says RALF23 “did not induce” DR5 in pWER-KO and pSMB-KO, while introducing this as similar to YUC5 and then calling pSCR-KO the weakest induction. P53 says YUC5 induction in pWER/pSMB is WT-like. This may be an editorial error, but the underlying DR5 images and measurements must determine the correction. No manuscript text or data were changed in this review.

## Audit/source links

- [Ten-DOF retained edges](../../2026-09-20_dof_subnetworks/inventory/tested10_retained_edges.csv)
- [Target-input membership](../../2026-09-08_routeA_execution/q04_network/results/h5ad_gene_universes.csv)
- [YUC5 target eligibility](../../2026-09-08_routeA_execution/q04_network/results/DOF_YUC5_target_eligibility.csv)
- [Existing network inventory and integrity notes](../../2026-09-20_dof_subnetworks/inventory/README.md)
- [Full saved-edge/degree audit](../../2026-09-20_dof_subnetworks/inventory/full_public_export_degree_audit.csv)
- [Prior network interpretation and eligibility audit](../../2026-09-20_fig6_brainstorm/network_audit.md)

No figure, manuscript, network data, or remote object was changed. Only this review note was created.
