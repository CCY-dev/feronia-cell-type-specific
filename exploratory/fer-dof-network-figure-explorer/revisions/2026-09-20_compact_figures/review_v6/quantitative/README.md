# PROTOTYPE: not yet carefully vetted

This additive matrix-only revision addresses independent audit findings F11, F12 and F13. The forest remains the unchanged v5 file. Both figures retain a native width of 110 mm for aligned assembly.

`celltype_tested10_review_v6.pdf/.svg/.png` is **110 × 90 mm**. Only the wrapped “Metaphloem and companion cell” label shifts left by 1 mm; full tissue names remain angled at 60 degrees. The adjustment is made to the generated R grid text grob, not by editing the exported SVG. Counts decrease from 9 to 8.5 pt to improve half-tile clearance. Tissue labels remain 8.5 pt, gene labels 9 pt and key text 8 pt. The unused x key is omitted because all 130 plotted entries have at least one paired-score zone; conditional missingness code remains for future inputs.

The actual cell width remains **7.094 mm**, not the 5.4 mm quoted in the initial audit F13. A full tile is 4.762 mm high, and each colored half is 2.381 mm high. The 8.5 pt fraction 1/4 is 4.170 mm wide, leaving approximately 1.46 mm horizontal padding on each side.

## Validation

- `source_manifest.json` verifies byte-identical source and plotted CSVs against v5: 130 entries, ten DOFs and thirteen tissues.
- `svg_validation.json` verifies a minimum text size of 8 pt, editable text and named SVG groups.
- `font_check.txt` records embedded Arial TrueType fonts without Type 3.
- `check_pdf_600dpi.png` and `check_long_label_crop_600dpi.png` show the rendered PDF. The long-label crop was inspected against `check_v5_long_label_crop_600dpi.png`; the adjusted label has clear separation from Protoxylem and Protophloem.
- The native-size display preview and deutan, protan and tritan previews were inspected. Upper/lower position provides redundant genotype coding. A display preview is not a physical print proof.
- The PDF has one page and a 110 × 90 mm page box. Previous versions, forest outputs and frozen source tables remain unchanged.

The v5 layout check missed the crowded long-label pair. The v6 check therefore uses a 600 dpi close crop in addition to the full figure and native-size preview. `validation.txt` summarizes the present checks.

## Inputs and scope

The R generator reads `revisions/2026-09-20_fig6_brainstorm/celltype_prototypes/celltype_summary_all37.csv`, preserving a copy as `source_celltype_summary_all37.csv`. The ten-assayed-DOF subset is `plot_source_celltype_tested10.csv`. Data, gene/tissue order, denominator definitions and nomination rules remain unchanged. This is a presentation revision, not an analysis repair. Saved edges and degree scores remain inconsistent in three fer-4 endodermal contexts, and phloem/xylem fits have separate target-mapping issues. See `CAPTION.md` for interpretation and source limitations.

Reproduce from the repository root:

```sh
Rscript revisions/2026-09-20_compact_figures/review_v6/quantitative/make_review_v6_matrix.R
python revisions/2026-09-20_compact_figures/review_v6/quantitative/tag_svg.py
python revisions/2026-09-20_compact_figures/review_v6/quantitative/validate_output.py
pdftoppm -png -r 600 -x 1410 -y 1200 -W 700 -H 650 -singlefile revisions/2026-09-20_compact_figures/review_v6/quantitative/celltype_tested10_review_v6.pdf revisions/2026-09-20_compact_figures/review_v6/quantitative/check_long_label_crop_600dpi
```

The dataviz skill and guide were applied; angled full labels follow the user's explicit request. **Further reading:** [Small axis labels](https://clauswilke.com/dataviz/small-axis-labels.html).
