# PROTOTYPE: not yet carefully vetted

This additive presentation revision increases type size and matches the widths of panels A and C. Source values, row/column order, denominators, selection rules and statistical estimates are unchanged. Previous versions and frozen inputs remain untouched.

## Dimensions and typography

- `top10_aligned_forest_review_v5.pdf/.svg/.png`: 110 × 80 mm, minimum text 8 pt; family labels 9 pt, counts 8.5 pt. The DOF family label remains bold without added color meaning.
- `celltype_tested10_review_v5.pdf/.svg/.png`: 110 × 90 mm, fractions and gene labels 9 pt, tissue labels 8.5 pt, two-line key 8 pt. The actual matrix tiles measure 7.094 × 4.762 mm. A fraction such as 1/4 is 4.413 mm wide, leaving 1.34 mm on each side.
- `four_families_review_v5.pdf/.svg/.png`: optional 110 × 43 mm companion, minimum text 8 pt.

Full tissue names appear below the matrix at 60 degrees, as explicitly requested. Metaphloem and companion cell wraps across two lines. Place A and C at their native 110 mm widths with identical left and right bounds in the assembled proof. Do not scale down and thereby reduce typography below these stated sizes.

## Checks

`source_manifest.json` verifies that all four plotting/source CSVs are byte-identical to review v3. There are 130 gene/tissue entries for ten DOFs and thirteen tissues, and twenty forest estimates for ten families. There are six BH-significant estimates and one zero odds ratio in the forest. `svg_validation.json` confirms editable text, named SVG groups, dimensions and minimum font sizes. `font_check.txt` records embedded Arial TrueType fonts without Type 3. `validation.txt` summarizes the checks.

PDF rasters were inspected at 240 dpi and at 96 dpi for an approximate native display-size check. Angled labels, keys, confidence intervals, count headings and zero-OR arrow were inspected without visible clipping or text collisions. Deuteranopia, protanopia and tritanopia simulations were inspected for both primary figures. Panel position and upper/lower tile position retain genotype information independently of hue. Screen-size inspection cannot replace a physical print proof.

## Data and interpretation

The matrix reads the frozen `revisions/2026-09-20_fig6_brainstorm/celltype_prototypes/celltype_summary_all37.csv`. Its plotted subset is written to `plot_source_celltype_tested10.csv`; `source_celltype_summary_all37.csv` preserves all 37 DOFs. Fractions are nominated developmental zones divided by zones with paired genotype scores for each DOF and tissue. They are context counts, not biological replicates or TF activity measurements.

The forest reads `revisions/2026-09-20_compact_figures/quantitative/source_top10_family.csv`, copied here unchanged. `plot_source_four_significant_families.csv` contains the optional companion subset. The eligible universe contains 1,292 paired-score TFs and 33 DOFs; it differs from the supplied 37-member DOF list. See `CAPTION.md` for full definitions.

These remain exploratory figures. Saved edges and degree scores are inconsistent in three fer-4 endodermal contexts, and phloem/xylem fits have separate target-mapping issues. This cosmetic revision does not repair those upstream limitations or add new analyses.

## Reproduction

Run from the repository root:

```sh
Rscript revisions/2026-09-20_compact_figures/review_v5/quantitative/make_review_v5_matrix.R
Rscript revisions/2026-09-20_compact_figures/review_v5/quantitative/make_review_v5_forest.R
python revisions/2026-09-20_compact_figures/review_v5/quantitative/tag_svg.py
python revisions/2026-09-20_compact_figures/review_v5/quantitative/validate_output.py
```

The lab dataviz skill and guide were applied. The user's request for angled full tissue names overrides the guide's default preference for unrotated labels.
