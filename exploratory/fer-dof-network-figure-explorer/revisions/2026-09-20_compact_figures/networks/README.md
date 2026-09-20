# Compact DOF network overviews

Primary preview: `vascular_overview_compact.png`, PDF and SVG, 140 × 75 mm. Supporting broader overview: `all_tissues_overview_compact.png`, PDF and SVG, 150 × 95 mm. Both use the same ten experimentally examined DOFs and the previously declared auxin-family selection rule. No inference was rerun.

| Overview | DOFs | Targets | Distinct source-target pairs | Contexts with selected links | Original context/genotype edge rows |
|---|---:|---:|---:|---:|---:|
| Elongating vascular contexts | 10 | 23 | 81 | 3 | 207 |
| Seven included tissue fits | 10 | 34 | 147 | 27 | 610 |

The primary overview is the union of Elongation_Phloem Pole Pericycle, Elongation_Xylem Pole Pericycle and Elongation_Procambium, across basal WT and fer-4 networks. This choice relates to the manuscript's vascular/pericycle YUC5 response and late-transition/elongation focus. It does not identify procambium as the principal site of DOF activity, and the manuscript reporter data do not distinguish PPP from XPP. The underlying basal networks are not acute RALF23-response networks.

The supporting overview uses all contexts represented by selected links from the seven included tissue fits: atrichoblast, trichoblast, cortex, pericycle, procambium, lateral root cap and columella. The 27 contexts count only contexts containing a selected edge. Endodermis, phloem and xylem are excluded owing to the unresolved integrity issues documented in `../../2026-09-20_dof_subnetworks/SUBNETWORK_DESIGN.md` and its inventory. Passing the existing checks for the seven remaining fits is not a comprehensive model audit.

## Selection and encoding

Source is `../../2026-09-20_dof_subnetworks/inventory/tested10_retained_edges_seven_unaffected_consistent_tissues.csv`. Gene symbols are normalized by AGI using `../../2026-09-20_dof_subnetworks/annotations/clean_programme_target_edges.tsv`; the original export symbol is retained as `target_export_name`.

Target symbols match `^(IAA[0-9]|SAUR[0-9]|GH3[.]|TAA[0-9]|TAR[0-9]|YUC[0-9])`. This transparent gene-family filter is reused from the single-context atlas. It is not a pathway-completeness claim or enrichment result. There is no ranking by coefficient, no selection for differential expression, and no top-target truncation.

Every gray line is one retained DOF-source to target pair in at least one included context and genotype. All lines have the same appearance because the figure is a context/genotype union. The union does not imply that all links coexist in one cell state, that the TF activates the target, or that links differ significantly between genotypes. Layout and node size encode no quantitative statistic. Target colors identify gene families, not expression or TF activity. Arrowheads are omitted for compactness; the top row contains the modeled sources and the lower arc contains targets. DOF-YUC5 connections were present in the supplied prior, so retention is not independent discovery of binding.

`*_source_context_edges.csv` preserves every selected original edge, including its genotype, cell type, tissue and coefficients. `*_union_edges.csv` provides one plotted source-target pair with exact WT/fer-4 context lists and counts. `*_context_counts.csv` reports context/genotype counts, and the node tables record exact labels and coordinates.

## Suggested captions

**Primary:** Predicted connections from ten tested DOF transcription factors to selected auxin-related gene families, combined across elongating phloem-pole pericycle, xylem-pole pericycle and procambium basal networks. Gray lines indicate a retained source-target link in at least one included WT or fer-4 context; colors identify target gene families. The union contains 81 distinct pairs connecting ten DOFs and 23 targets. It summarizes modeled associations across contexts and does not imply simultaneous interactions or activation.

**Supporting overview:** Predicted DOF connections to the same selected auxin-related gene families, combined across the seven included tissue fits and both basal genotypes. Gray lines indicate a retained source-target link in at least one context. The union contains 147 distinct pairs connecting ten DOFs and 34 targets, with selected edges in 27 contexts. Endodermis, phloem and xylem fits are excluded pending the source-integrity review.

## Reproduction and checks

Run `/usr/local/bin/Rscript make_overviews.R`, then `python3 tag_svg_groups.py`. The R script resolves input/output locations relative to itself. `input_provenance.txt` records exact inputs and MD5 hashes, `R_sessionInfo.txt` records package versions, and `validation.csv` contains assertions verified during generation. SVG text remains editable and logical layers have named groups. PDFs use native Quartz on this Mac because Cairo's XQuartz library is unavailable.

Both plot renders were visually inspected. The dense supporting arc uses radial labels to resolve lower-arc overlap. Target labels are 7.25 pt; TF labels are 7.4 pt. Deutan, protan and tritan simulations are supplied using `colorspace`; symbols and contiguous gene-family blocks provide redundancy where simulated hues approach one another. This compact network remains useful as a relationship overview, while exact context-specific statements should be made from the source tables or matched WT/fer-4 atlas.

Further reading: [Nolan et al., Science (2023), Fig. 6](https://doi.org/10.1126/science.adf4721), the design reference for the TF-to-functional-target subnetwork presentation. The present union is explicitly broader than the single biological context used in that reference.
