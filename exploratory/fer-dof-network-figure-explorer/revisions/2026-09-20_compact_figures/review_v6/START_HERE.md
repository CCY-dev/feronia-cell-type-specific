# FER Figure 6 code and provenance reference

**Prototype; not yet carefully vetted.** The endodermal graph/score and phloem/xylem target-selection discrepancies remain unresolved. The independent audit has been delivered and its receipt verified, but its findings remain unresolved. Calculation and layout checks do not establish scientific validity or publication readiness.

The current code inventory is `revisions/2026-09-20_compact_figures/review_v6/provenance/CURRENT_INPUT_OUTPUT_CODE.csv` inside the archive. Panel A retains the v5 forest. Panel B uses the v6 DOF evidence bars, with stack starts calculated from counts. Panel C uses the v6 matrix, which reads only `celltype_summary_all37.csv`. The other panel source versions remain recorded in that inventory.

The v5 membership reconstruction and its results remain under `revisions/2026-09-20_compact_figures/review_v5/provenance/`. Its nine table inputs are included at their original repository-relative paths. It reproduces all membership flags and seven intersection counts; 17 display-symbol differences remain explicitly recorded. This new reconstruction does not claim that the historical generator was recovered.

This ZIP is **not a self-contained portable analysis pipeline**. It contains marked code references, current plotting scripts, selected derived tables, the v5 forest's immediate input and membership-reconstruction table dependencies. Raw network objects, the full public source checkout and its Git metadata, all upstream source data, installed libraries and remote environments are not bundled. Scripts can require original project paths or Resnick paths and can write alongside themselves. Inspect those dependencies and use an isolated project copy before running them. The membership reconstruction also inspects the pinned public clone's Git revision.

Original repository-relative paths are preserved. Historical reference scripts use the marked copies recorded in `reference_metadata/source_copy_manifest.csv`. Where needed, a prototype comment is added only to archived script bytes, preserving shebangs and executable content; original files remain unchanged. `BUNDLE_FILES.csv` records source and archived hashes and whether a header was added.

Plot recreation and summary reconstruction are distinct from rebuilding the underlying networks. No network fit, graph repair, GitLab publication or external collaborator transfer is performed by creating this archive.
