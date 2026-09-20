# PROTOTYPE: not yet carefully vetted

Panel d is an additive **135 × 90 mm** presentation derivative of `review_v2/all_contexts/all_contexts_overview_compact`. It preserves all ten DOFs, 37 targets and 180 distinct source-target pairs. The assembler supplies the external title and panel letter, so no large internal title is included.

## Design

Ten DOFs occupy the top row; targets follow the same family/numeric order around the lower arc. The figure uses 9 pt DOF labels and 8.5 pt target and key labels. TF and target point-size settings are 3.8 and 2.7; their actual SVG diameters are approximately 2.985 and 2.131 mm. Edge strokes export at 0.32 pt with 0.50 opacity. Legend strokes are stronger at 0.64 pt and full opacity.

The compact key explicitly counts pairs: 147 also supported in the seven-fit reference subset, and 33 additional pairs found only in the other, previously excluded fits. Solid and dashed lines preserve the original support classes. Node colors identify the four target families; gene-name prefixes and family ordering provide redundant cues.

The dense all-context union conveys the extent of modeled connections. It does not make every crossing edge individually traceable at this size. The exact edge table is retained for pair-level inspection. No edges are bundled, removed or weighted by expression, effect size, confidence or context count.

## Files and validation

- `all_contexts_panel_d.pdf/.svg/.png`: final display files.
- `make_panel_d.R`: R generator with the prototype header.
- `source_union_edges.csv`: byte-identical copy of the v2 union table.
- `plotted_edges.csv`, `plotted_tf_nodes.csv`, `plotted_target_nodes.csv`: source values with display coordinates.
- `finalize_and_validate.py`: SVG group tagging and independent identity/render checks.
- `source_manifest.json`: source hash, equality results, exact counts and dimensions.
- `font_check.txt`: embedded Arial TrueType fonts, no Type 3.
- `check_pdf_600dpi.png`, `check_native_size.png`, and two 600 dpi label crops: visual checks.
- `cvd_deutan.png`, `cvd_protan.png`, `cvd_tritan.png`: color-vision simulations.

All original edge columns match the source table after key-based comparison. The SVG contains exactly 147 solid and 33 dashed network segments, ten DOF labels and 37 target labels. All TF nodes are above all targets. The PDF is one page at exactly 135 × 90 mm. The minimum text size is 8.5 pt. Full preview, native-size display preview, 600 dpi label crops and all three CVD previews were inspected without visible label collisions or clipping. This is a screen-based production check, not a physical print proof.

No prior figure, source table or analysis code was changed. The unresolved endodermal edge/score inconsistencies and phloem/xylem target-mapping issues remain. These outputs do not repair the source analyses or establish TF activity, simultaneous interactions, activation/repression or a genotype difference.

## Reproduction

Run from the repository root:

```sh
Rscript revisions/2026-09-20_compact_figures/review_v7/networks/make_panel_d.R
python revisions/2026-09-20_compact_figures/review_v7/networks/finalize_and_validate.py
```

The lab dataviz skill and guide were applied. **Further reading:** [Multi-panel figures](https://clauswilke.com/dataviz/multi-panel-figures.html).
