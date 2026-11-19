# Figure kit: the main figures and the Extended Data charts as published

The authors' figure code (figure kit by Trevor Nolan) that draws the paper's figures in the journal style:

- **Computational panels:** Fig. 1 (a–d) and Fig. 5a, b, from the outputs of the analyses in this repository.
- **Wet-lab charts:** every quantification chart of Figs. 2–6 and of the Extended Data figures, from the paper's
  **Source Data** workbook. The workbook is published with the paper and is not part of this repository.

The analysis scripts elsewhere in this repository and the Source Data produce the numbers; this kit turns them into
tables and lays them out.

## Quick start (Python 3.13)
```
python -m pip install -r requirements.txt

# computational panels
python tools/fig1_tables.py --check     # recompute data/fig1/ from this repository and compare (11 checks)
python tools/fig5ab_tables.py --check   # recompute data/fig5/ (a, b) from this repository and compare
python build.py fig1                    # -> output/figures/Fig1.{pdf,png,svg}
python build.py fig5ab                  # -> output/figures/Fig5ab.{pdf,png,svg}

# wet-lab charts (put the Source Data workbook in source_data/ first)
python tools/source_data_to_long.py "source_data/<Source Data>.xlsx"   # -> source_data/long/ (tidy tables)
python tools/fig2_tables.py             # -> data/fig2/ measurement tables + output/stats_checks/STATS_CHECK_fig2.csv
python tools/fig3_tables.py
python tools/fig4_tables.py
python tools/fig5_tables.py "source_data/<Source Data>.xlsx"   # Fig. 5d LUC tables (the xlsx adds the workbook P labels to the check)
python tools/fig6_tables.py
python tools/ed_charts.py               # -> output/ed_charts/<chart>.{pdf,png,svg} + STATS_CHECK_ED.csv
python tools/make_placeholder_tiles.py  # blank stand-ins for the images (see below)
python build.py fig2                    # likewise fig3, fig4, fig5, fig6 -> output/figures/
```
Outputs go to `output/`; `source_data/`, the measurement tables and the images are not tracked.

**Font.** The figures use Arial (journal requirement). Arial is licensed and is **not** included: install it on your
system (or point `FER_FONT_DIR` at a folder with ARIAL.TTF, ARIALBD.TTF, ARIALI.TTF, ARIALBI.TTF) before building.
Without it the kit falls back to Liberation Sans / DejaVu Sans and warns; text widths then differ slightly.

## What is and is not in this repository
- **Source Data:** not included (published with the paper). `tools/source_data_to_long.py` reads it; every number
  on a chart is computed from it at build time.
- **Images** (micrographs, seedling photographs, the EMSA gel): not included. `tools/make_placeholder_tiles.py`
  writes a uniformly grey tile of the same pixel size for every image listed in `data/figN/figN_image_sizes.csv`, so
  the figures build with every chart, label and scale bar at its published position. To rebuild the complete
  figures, put the original images at the paths in the manifests `data/figN/figN_images.csv` (crop, orientation and
  µm per pixel are given there). The three colour-scale keys (`data/fig3|4|6/images/fig*_lut_key.png`) are included.
- **Statistics labels** (letters and asterisks as printed by the authors, and the printed treatment/mock ratios):
  `data/figN/*_stats.csv`, `data/fig2/fig2h_ratios.csv`, `data/fig4/fig4[ef]_*ratios.csv`, the `ratio` column of `data/fig6/fig6b_seedlings.csv`,
  `data/fig5/fig5d_luc_displayed_values.csv` and `data/ed/ED_author_labels.csv`. Each table script re-runs the test
  named in the figure legend on the Source Data and compares it with these labels (`output/stats_checks/`). The
  public GitHub copy does not contain these label files; it carries `data/LABELS_NOT_INCLUDED`, and the charts are
  then drawn without labels (take them from the published figures) and the checks list the computed labels only.

## Inputs → tables → figure
| Figure | Tables (`data/`) | Made by | From |
|---|---|---|---|
| 1a | `fig1/fig1a_umap_cells.csv.gz` (per-cell UMAP + annotations), `fig1a_palette.csv` | `tools/hpc/extract_fig1a_umap.R` (needs the integrated Seurat object, data deposit); palette from the integration notebook | Stage 1 object; `Preprocessing/scripts/202501226_Ref_integration_…_seu3.ipynb` |
| 1b, c | `fig1/fig1b_genes_by_n_subpopulations.csv`, `fig1c_deg_per_subpopulation.csv` | `tools/fig1_tables.py` | `Differential_expression_analysis/outputs/sig_de_fil.csv` |
| 1d | `fig1/fig1d_hormone_overlap_q_authors.csv` (drawn), `fig1d_hormone_overlap.csv`, `fig1d_colour_scale.csv` | `tools/fig1_tables.py` (the q table is copied unchanged) | `Differential_expression_analysis/outputs/hormone_overlap/Fig1d_overlap_q.csv` |
| 2d, f, g (+ h ratios) | `fig2/fig2d_meristem_cell_number.csv`, `fig2f_root_hair_length.csv`, `fig2g_root_length_ralf23.csv` | `tools/fig2_tables.py` | Source Data sheets `Fig_2d&Ext_Fig_4a`, `Fig_2f&Ext_Fig_4d`, `Fig_2g` |
| 3b, d, g | `fig3/fig3b_yuc5_reporter_timecourse.csv`, `fig3d_yuc5_reporter_col0_fer4.csv`, `fig3g_yuc5_reporter_tsko.csv` | `tools/fig3_tables.py` | sheets `Fig_3b`, `Fig_3d`, `Fig_3g` |
| 4b, d, e, g | `fig4/fig4b_dr5_col0_fer4.csv`, `fig4d_dr5_tsko.csv`, `fig4e_root_length_peo_iaa_L3.csv`, `fig4g_root_length_afb1_tir1.csv` | `tools/fig4_tables.py` | sheets `Fig_4b`, `Fig_4d`, `Fig_4e`, `Fig_4g` |
| 5a, b | `fig5/fig5a_TF_family_enrichment_all_tests.csv`, `fig5b_network_{dof_nodes,target_nodes,edges,layout}.csv` | `tools/fig5ab_tables.py` | `CellOracle_GRN_analysis/outputs/TF_family_enrichment/` (Stage 5) |
| 5c | `fig5/fig5c_promoter_elements.csv`, `fig5c_probe_sequences.csv`, `fig5/promoter/` | `tools/fig5_promoter/derive_coordinates.py`, then `make_promoter_table.py` | TAIR10 / Araport11 / Ensembl / NCBI records in `tools/fig5_promoter/sources/` (re-fetch: `fetch_sources.sh`) |
| 5d | `fig5/fig5d_luc.csv`, `fig5d_luc_wells.csv`, `fig5d_luc_rep_map.csv` | `tools/fig5_tables.py` | sheet `Fig_5d LUC assay` |
| 5e | `fig5/fig5e_emsa_lanes.csv` (lane conditions) | – | gel image (not included) |
| 6a, d, e, g | `fig6/fig6a_root_length.csv`, `fig6d_pyuc5_dof.csv`, `fig6e_qpcr_{dof46,cdf4,dof22,dag1}.csv`, `fig6g_{pdag1,pdof22}.csv` | `tools/fig6_tables.py` | sheets `Fig_6a`, `Fig_6d`, `Fig_6e`, `Fig_6g` |
| ED 2b, c, e; 4a, b, d, f; 5g; 6c; 7b, d, f; 9d, e; 10c | drawn directly (chart sizes in `ed/ED_chart_sizes.csv`) | `tools/ed_charts.py` | the `Ext_Fig_*` sheets |

The Extended Data figures themselves combine these charts with images and schematics on the authors' page layouts;
the kit draws the charts (`output/ed_charts/`), not the full pages. `fig1d_hormone_overlap_q_proposed.csv` is an
alternative statistic that is not drawn (the config uses `statistic: q_authors`).

## Layout of the kit
- `build.py` – entry point; `ferfigs/` – the kit (page model, house style, journal checks, panel drawing code);
  `figures/fig1.py` … `fig6.py`, `figures/fig5ab.py` – one module per figure; `config/` – every word, colour and
  position (`style.yaml` = house style, one YAML per figure); `ferfigs/profiles/nature_portfolio.yaml` – journal limits.
- `tools/` – the table scripts above, `make_placeholder_tiles.py`, and `hpc/extract_fig1a_umap.R`.
- `python build.py check` runs the journal checks on the built figures.

## Verification (19 Nov 2026)
- `tools/fig1_tables.py --check` and `tools/fig5ab_tables.py --check`: every table equals a fresh derivation from this
  repository. With Arial installed, `Fig1.pdf` is pixel-identical to the published Fig. 1; `Fig5ab.pdf` has the same
  text and positions as the top 62 mm of the published Fig. 5.
- From the Source Data workbook: with the original images in place, Figs. 2–6 and all 22 Extended Data charts are
  pixel-identical to the published versions. With the blank tiles, every vector element and text span of Figs. 2–6
  (charts, axes, labels, scale bars) and every image position are identical; only the image pixels differ.
- Statistics checks (legend test re-run on the Source Data vs the printed labels): Fig. 2 30 of 31 (the remaining row
  is an alternative experimental unit, per-root means, shown for information), Fig. 3 25/25, Fig. 4 48/48,
  Fig. 5 20/20, Fig. 6 36/36, Extended Data 135/135.
