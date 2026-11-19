"""Figure 5a, b: TF-family enrichment forest plot (a) and the DOF -> auxin-gene network (b).

The computational panels of the paper's Fig. 5, drawn with the authors' figure kit (Trevor Nolan); panels c-e are
wet-lab results and are not part of this repository. Inputs: data/fig5/fig5a_*, fig5b_* (written by
tools/fig5ab_tables.py from the CellOracle outputs of this repository). Rebuild: python build.py fig5ab
"""
from ferfigs.panels import fig6_forest, fig6_network
from ferfigs.panels.fig6_common import MM

import pandas as pd


def build(ctx):
    cfg = ctx.cfg
    P = cfg['panels']
    page = ctx.page()
    mm = MM(page)
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

    for letter, (x, y_top) in cfg['letters'].items():
        page.panel_letter(letter, x, y_top)
    return page


def validate(ctx):
    """python build.py validate: the Fig. 5a,b tables (plain-language errors: missing file, columns, non-numeric values)."""
    P = ctx.cfg['panels']
    hint = 'Regenerate with tools/fig5ab_tables.py.'
    checks = [
        (P['a']['table'], ['genotype', 'family', 'nominated', 'eligible', 'OR', 'lo', 'hi', 'q'],
         ['nominated', 'eligible', 'OR', 'lo', 'hi', 'q']),
        (P['b']['edges'], ['source', 'target'], []),
        (P['b']['dof_nodes'], ['source', 'outdeg', 'label', 'group'], ['outdeg']),
        (P['b']['target_nodes'], ['target', 'target_sym', 'target_class', 'indeg'], ['indeg']),
    ]
    if ctx.path(P['b']['layout']).exists():
        checks.append((P['b']['layout'], ['node', 'side', 'y_mm'], ['y_mm']))
    for rel, req, num in checks:
        ctx.custom_table(rel, req, num, hint=hint)
    return len(checks)
