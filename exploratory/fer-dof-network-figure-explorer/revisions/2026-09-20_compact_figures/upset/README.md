# Compact UpSet candidate

**85 × 55 mm.** This revision keeps only a tightly aligned membership matrix and horizontal count bars. Set totals appear under the column names. No oversized title, marginal bar chart or in-panel prose is retained. Minimum text is 7.2 points at the stated output dimensions; headers are 7.7 points and count labels 7.8 points.

Files: `compact_upset.pdf`, `compact_upset.svg`, and `compact_upset.png`. SVG retains editable text and named groups. R source: `make_compact_upset.R`; SVG grouping helper: `group_svg.py`. Run the R script from any working directory. Previous candidates remain unchanged.

## Proposed caption

**Overlap of DOF family membership, GRN nomination and single-cell differential expression.** Filled dots identify set membership and pale dots indicate absence. Each row is an exclusive overlap involving at least two sets; horizontal bars give gene counts. Twelve DOFs are both GRN nominated and called differentially expressed, 14 are GRN nominated without a DE call, and none has a DE call without GRN nomination. An additional 172 genes belong to the GRN and DEG sets but are not DOFs. Set totals are 37 DOFs, 512 GRN-nominated TFs and 5,466 scDEGs. Singleton intersections, omitted from the matrix, contain 11 DOFs, 314 GRN-TFs and 5,282 scDEGs.

Membership is pooled across cell contexts. A gene's GRN nomination and DE call need not occur in the same context. GRN and DE analyses share scRNA-seq data, and differential expression informed target selection; this overlap is not independent experimental validation. Set definitions are inherited from the audited source tables and have not been recalculated after subsequent network-integrity findings.

## Inputs and checks

- Audited count source: `../../2026-09-20_fig6_brainstorm/fig6_full_intersections.csv`.
- `plot_input_overlaps.csv` holds the four plotted rows: 12, 14, 0 and 172. The zero remains explicit.
- `all_exclusive_intersections.csv` preserves all seven rows. `set_totals.csv` preserves set totals.
- Assertions verify all counts, every set total, seven unique membership combinations, and union size 5,805.
- Every bar has a zero origin and shares the same linear 0–200 axis; there are no broken axes or pseudocounts.
- CVD simulations using colorspace are supplied for deuteranopia, protanopia and tritanopia. Membership position and count labels provide redundant interpretation independently of color.
- Native Quartz PDF avoids the unavailable local XQuartz/Cairo dependency. SVG is the editable vector master.
- PDF page dimensions were verified as 240.945 × 155.905 points (85 × 55 mm), with embedded TrueType Arial fonts. The PDF was rendered at 300 dpi and at 96 dpi (about 321 × 208 pixels), then visually checked for readability, label collisions, clipping and alignment. The main PNG and deuteranopia simulation were also inspected.
- The plotting code follows the lab dataviz skill. At the original 85 mm width, all text remains at least 7.2 points; do not rescale smaller without redesign.

Further reading: [Wilke, balance between data and context](https://clauswilke.com/dataviz/balance-data-context.html).
