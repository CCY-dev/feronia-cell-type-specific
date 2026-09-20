# Design addendum for review_v5 (claude-science)

Prototype; not yet carefully vetted. Design notes only; no data change is proposed here and none is approved.

Inspected: `review_v5/review_at_print_size_v5.pdf` (3 pages, 200 dpi), `quantitative/celltype_tested10_review_v5.pdf` at 600 dpi actual size, `quantitative/make_review_v5_{forest,matrix}.R`, `overlap/make_dof_evidence.R`, `provenance/*`. Plot-source CSVs in v5 are byte-identical to v2; v5 changes layout only, as stated.

Confirmed in v5: prototype label on every proof page and in code; panels a and c both 110 mm and aligned; forest minimum 8.0 pt (svg_validation.json); DOF row bold; matrix counts 9 pt; evidence-bar text enlarged; full cell-type labels below the matrix; page-3 note now points to the source record.

Coordinator status received after this pass (v6): item 5 applied in the caption; hard-coded stack starts replaced by cumulative counts; items 1 and 3 being fixed without abbreviating full names.

## Material remaining issues (in priority order)

1. **Label collision, panel c.** At actual size the two-line 60° label "Metaphloem and companion cell" runs into "Protoxylem". Fix options: single-line "Metaphloem & CC" (define CC in the key), or rotate the column labels 90°, or widen only that column's gap. (findings F12)
2. **Vertical footprint, panel c.** The angled full labels take ~35 of the 90 mm; the 10-row matrix itself is ~45 mm. If page space matters, 90° labels or abbreviations with a key recover ~15 mm without changing the data. Trade-off against Trevor's request for full labels; report only.
3. **Dead key entry, panel c.** "x: no paired scores" describes a state that never occurs in the ten-DOF matrix (all 130 tiles have K>0). Drop it here; keep it only in the 26-DOF supplementary version where K=0 tiles exist. (F11)
4. **Tile-count fit, panel c.** 9 pt counts span the full ~2.38 mm half-tile height with no vertical margin (v5 tile 7.094 × 4.762 mm, count text 4.413 mm wide, 1.34 mm side padding, per `quantitative/validation.txt`; width is not the constraint). Legible; 8.5 pt or a slightly taller row would restore a vertical margin. (F13)
5. **Statistical notation, panel a.** Add "one-sided Fisher enrichment test; two-sided 95 % conditional-MLE CI; BH across 55 families within each genotype" to the caption or key line. (F06)
6. **Exposure statement, panels a and c.** Wherever endodermal nominations appear (panel c endodermis column; DOF4.6/DAG1 tiles), the caption should say that 11 of the 12 endodermal entries for the ten assayed DOFs, and 78 of the 512 GRN-TFs, rest on the three fer-4 endodermal contexts whose saved graphs and scores disagree. These are direct dependencies under the current pooled thresholds — not a prediction that they will change, and not an upper bound, because recomputed pooled rank tails can also move nominations in other contexts. (F03)
7. **Panel d hairball.** Unchanged from v2; legible only as edge density. A 10 × 37 TF-by-target matrix (fill = support class, glyph colour = family) carries the same 180 pairs traceably at the same footprint. Preference, not a defect. (F15)

Not issues: colour choices (Okabe–Ito, CVD-simulated), CI axis range, the 37 vs 33 denominator note, the reconstruction's 17 symbol aliases (all valid: DOF1=DOF1.7, DAG2=DOF2.5, DOF2=DOF3.1, TMO6=DOF5.3, OBP4=DOF5.4, HCA2=DOF5.6, and AGI-only names for genes lacking a Kay symbol).
