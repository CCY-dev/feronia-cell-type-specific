"""Figure 5: network analysis identified the DOF family TFs that regulate YUC5 expression.

Reads config/fig5.yaml and data/fig5/:
  a  TF-family enrichment forest plot     fig5a_TF_family_enrichment_all_tests.csv (CellOracle analysis, Stage 5)
  b  DOF -> auxin-gene network            fig5b_network_* tables (tools/fig5ab_tables.py)
  c  pYUC5 promoter with DOF sites        coordinate tables from TAIR10/Araport11 (tools/fig5_promoter/)
  d  pYUC5-LUC heatmap                    fig5d_luc.csv (fold change per DOF x replicate from the Source Data,
                                          tools/fig5_tables.py); the displayed-values table is used only without it
  e  EMSA                                 gel image + lane table; words from the config
Every number drawn comes from a table in data/fig5/; every word from the config or config/labels.yaml.
"""
from ferfigs.panels import fig6_forest, fig6_gel, fig6_heatmap, fig6_network, fig6_promoter
from ferfigs.panels.fig6_common import MM

import pandas as pd


def build(ctx):
    cfg = ctx.cfg
    P = cfg['panels']
    page = ctx.page()
    mm = MM(page)
    fig_dir = ctx.path('data/fig5')
    man = ctx.table(cfg['manifest'], 'image_manifest')
    coloured = bool(cfg.get('coloured_text', False))

    # ---- a: TF-family enrichment
    pa = dict(P['a'], coloured_text=coloured)
    df_a = fig6_forest.read(ctx.path(pa['table']))
    rec = fig6_forest.draw(page, mm, pa['box'], df_a, pa)
    rec['source_tables'] = {'TF-family enrichment, all 110 tests (CellOracle rerun 2026-10-21, unchanged)': pa['table']}
    ctx.record('a', rec)

    # ---- b: network
    pb = P['b']
    ed = pd.read_csv(ctx.path(pb['edges']))
    dn = pd.read_csv(ctx.path(pb['dof_nodes']))
    tn = pd.read_csv(ctx.path(pb['target_nodes']))
    lay = pd.read_csv(ctx.path(pb['layout'])) if ctx.path(pb['layout']).exists() else None
    rec = fig6_network.draw(page, mm, pb['box'], ed, dn, tn, lay, pb)
    if rec['layout'] != 'layout_table':
        ctx.warnings.append('fig5 b: node set differs from fig5b_network_layout.csv; automatic layout used')
    rec['source_tables'] = {'edges drawn': pb['edges'], 'DOF nodes': pb['dof_nodes'], 'target nodes': pb['target_nodes'],
                            'node layout (y, mm before scaling)': pb['layout']}
    ctx.record('b', rec)

    # ---- c: promoter
    pc = dict(P['c'])
    if coloured:
        pc['highlight'] = 'text'
    el = pd.read_csv(ctx.path(pc['elements']))
    sq = pd.read_csv(ctx.path(pc['sequences']), dtype=str, keep_default_na=False)
    rec = fig6_promoter.draw(page, mm, pc['box'], el, sq, pc, cfg['fonts']['sequence'])
    rec['source_tables'] = {'promoter elements (TSS-relative)': pc['elements'], 'probe sequences': pc['sequences']}
    ctx.record('c', rec)

    # ---- d: LUC heatmap (measurements when supplied, else displayed values)
    pd_ = P['d']
    meas = ctx.path(pd_['values_measured'])
    if meas.exists():
        vals, how = fig6_heatmap.from_measurements(meas, pd_['effectors'], pd_['replicates'])
        from ferfigs import provenance as PV          # status from the registry, not from the file's existence
        st_d, w_d = PV.data_status(pd_['values_measured'], pd_.get('data_status_measured', 'author_supplied'),
                                   ctx.kit, ctx.registry())
        P['d']['data_status'] = st_d
        if w_d:
            ctx.data_warn(f'd: {w_d}')
        src = pd_['values_measured']
    else:
        vals, how = fig6_heatmap.from_displayed(ctx.path(pd_['values_displayed']), pd_['effectors'], pd_['replicates'])
        src = pd_['values_displayed']
    rec = fig6_heatmap.draw(page, mm, pd_['box'], vals, pd_)
    rec.update(values_source=how, source_tables={'LUC values drawn': src})
    ctx.record('d', rec)

    # ---- e: EMSA
    pe = P['e']
    lanes = pd.read_csv(ctx.path(pe['lanes']))
    row_e = man[man['panel'] == 'e'].iloc[0]
    rec, _ = fig6_gel.draw(page, mm, fig_dir, row_e, lanes, pe, pe['gel_x'], pe['gel_y'], pe['gel_w'])
    rec['source_tables'] = {'lane centres (gel pixels)': pe['lanes']}
    ctx.record('e', rec)

    for letter, (x, y_top) in cfg['letters'].items():
        page.panel_letter(letter, x, y_top)
    return page


def validate(ctx):
    """python build.py validate: the Fig. 5 tables that are not one of the standard schemas (plain-language
    errors: missing file, missing columns, non-numeric values). Optional inputs are checked only when present."""
    P = ctx.cfg['panels']
    hint = 'Regenerate with ROOT/scripts/25_fig5_tables.py (a, b, d) or see docs/notes/fig6_NOTES.md (c, e).'
    checks = [
        (P['a']['table'], ['genotype', 'family', 'nominated', 'eligible', 'OR', 'lo', 'hi', 'q'],
         ['nominated', 'eligible', 'OR', 'lo', 'hi', 'q']),
        (P['b']['edges'], ['source', 'target'], []),
        (P['b']['dof_nodes'], ['source', 'outdeg', 'label', 'group'], ['outdeg']),
        (P['b']['target_nodes'], ['target', 'target_sym', 'target_class', 'indeg'], ['indeg']),
        (P['c']['elements'], ['element', 'strand', 'tss_rel_start', 'tss_rel_end'], ['tss_rel_start', 'tss_rel_end']),
        (P['c']['sequences'], ['name', 'sequence_5to3', 'highlight_1based'], []),
        (P['d']['values_displayed'], ['effector', 'replicate', 'fold_change', 'label'], ['fold_change']),
        (P['e']['lanes'], ['lane', 'x_px', 'group'] + [r['column'] for r in P['e']['rows']], ['lane', 'x_px']),
    ]
    if ctx.path(P['b']['layout']).exists():
        checks.append((P['b']['layout'], ['node', 'side', 'y_mm'], ['y_mm']))
    for rel, req, num in checks:
        ctx.custom_table(rel, req, num, hint=hint)
    if ctx.path(P['d']['values_measured']).exists():     # the same reader the build uses: numbers, ids, complete design
        fig6_heatmap.from_measurements(ctx.path(P['d']['values_measured']), P['d']['effectors'], P['d']['replicates'])
        return len(checks) + 1
    return len(checks)
