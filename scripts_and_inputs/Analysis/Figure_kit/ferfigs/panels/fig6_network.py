"""Fig. 6b: predicted DOF -> auxin-gene target network (bipartite: DOFs left, targets right), ported from the R
panel (compact_panels.R) to the kit.

Inputs (unchanged upstream tables, derived from Chia-Yu Chen's CellOracle edge export by the ten-DOF run):
  ten_tie_inclusive_dof_nodes.csv     source (AGI), outdeg, label, group
  ten_tie_inclusive_target_nodes.csv  target (AGI), target_sym, target_class, indeg
  ten_tie_inclusive_edges.csv         source, target, target_sym, ... (one row per drawn edge)
Layout:
  fig6b_network_layout.csv            node, side, label, y_mm: vertical positions of the current figure (from the
                                      R code), used when every node of the node tables has a row there;
  otherwise a deterministic automatic layout (groups / classes in config order, then degree, then name) is used
  and the build record says so.
DOF node AREA is proportional to the number of displayed targets (outdeg): d = node_d_per_sqrt * sqrt(outdeg).
"""
from __future__ import annotations

import math

import pandas as pd

from ..schemas import ValidationError
from .fig6_common import MM, cap_mm


def auto_layout(dn, tn, cfg):
    """Deterministic layout for new data: y in the same units as the layout table (mm before scaling)."""
    al = cfg['auto_layout']
    rows, y = [], 0.0
    prev = None
    tn = tn.assign(_c=tn['target_class'].map({c: i for i, c in enumerate(cfg['class_order'])}).fillna(99))
    for _, r in tn.sort_values(['_c', 'indeg', 'target_sym'], ascending=[True, False, True], kind='stable').iterrows():
        if prev is not None:
            y += al['target_pitch'] + (al['class_gap'] if r['target_class'] != prev else 0.0)
        rows.append(dict(node=r['target'], side='target', y_mm=y))
        prev = r['target_class']
    span = y
    dn = dn.assign(_g=dn['group'].map({g: i for i, g in enumerate(cfg['group_order'])}).fillna(99))
    dn = dn.sort_values(['_g', 'outdeg', 'label'], ascending=[True, False, True], kind='stable')
    n, ng = len(dn), dn['group'].nunique()
    pitch = (span - (ng - 1) * al['group_gap']) / max(n - 1, 1) if n > 1 else 0.0
    pitch = max(pitch, 0.5)
    y, prev = 0.0, None
    for _, r in dn.iterrows():
        if prev is not None:
            y += pitch + (al['group_gap'] if r['group'] != prev else 0.0)
        rows.append(dict(node=r['source'], side='dof', y_mm=y))
        prev = r['group']
    return pd.DataFrame(rows)


def draw(page, mm: MM, box, ed, dn, tn, lay, cfg, gid='fig6-b'):
    st = page.st
    x, y, w, h = box
    fs_head = cfg.get('header_pt', st['font']['size_pt'])
    fs_head_t = cfg.get('header_targets_pt', fs_head)
    fs_node = cfg.get('node_label_pt', st['font']['small_pt'])
    fs_dof = cfg.get('dof_label_pt', fs_node)                # DOF names
    fs_tgt = cfg.get('target_label_pt', fs_node)             # target gene names
    fs_grp = cfg.get('group_label_pt', fs_node)              # nomination groups (WT only, ...)
    fs_cls = cfg.get('class_label_pt', fs_node)              # target classes (GH3, ...)
    fs_hi = cfg.get('highlight_label_pt', st['font']['size_pt'])
    # ---- checks
    probs = []
    for col in ('source', 'target'):
        if col not in ed.columns:
            probs.append(f'edges table: missing column {col!r}')
    if probs:
        raise ValidationError(probs)
    unknown = sorted(set(ed['source']) - set(dn['source'])) + sorted(set(ed['target']) - set(tn['target']))
    if unknown:
        raise ValidationError([f'edges refer to nodes missing from the node tables: {unknown[:6]}'])
    # ---- layout: the stored table when it covers every node, else automatic
    need = set(dn['source']) | set(tn['target'])
    have = set(lay['node']) if lay is not None else set()
    mode = 'layout_table' if need <= have else 'automatic'
    L = lay if mode == 'layout_table' else auto_layout(dn, tn, cfg)
    ymap = dict(zip(L['node'], L['y_mm'].astype(float)))
    span = max(ymap[n] for n in need)
    y_first, y_last = y + cfg['y_first'], y + cfg['y_last']
    scale = (y_last - y_first) / span if span > 0 else 1.0

    def Y(node):
        return y_first + ymap[node] * scale
    X = {k: x + v for k, v in cfg['x'].items()}         # group_bracket, group_label, dof_label, dof_node, ...
    hi = cfg['highlight_target']
    pal = cfg['class_colours']
    # ---- edges (highlighted edges last, on top)
    tsym = dict(zip(tn['target'], tn['target_sym']))
    e = ed.assign(_hi=[tsym.get(t) == hi for t in ed['target']]).sort_values('_hi', kind='stable')
    ew = cfg['edge']
    for _, r in e.iterrows():
        sty = ew['highlight'] if r['_hi'] else ew['other']
        mm.line(X['dof_node'], Y(r['source']), X['target_node'], Y(r['target']), lw=sty['lw_pt'],
                colour=sty['colour'], alpha=sty.get('alpha'), gid=f'{gid}-edge{"-hi" if r["_hi"] else ""}',
                z=2 if r['_hi'] else 1, cap='round')
    # ---- DOF nodes and labels
    k = cfg['node_d_per_sqrt']
    for _, r in dn.iterrows():
        yy = Y(r['source'])
        mm.circle(X['dof_node'], yy, k * math.sqrt(float(r['outdeg'])), face=cfg['dof_colour'], gid=f'{gid}-dof', z=4)
        page.text(X['dof_label'], yy + cap_mm(fs_dof) / 2, r['label'], size=fs_dof, ha='right', fontstyle='italic',
                  gid=f'{gid}-doflabel')
    # ---- target nodes and labels
    for _, r in tn.iterrows():
        yy = Y(r['target'])
        col = pal.get(r['target_class'], '#808080')
        mm.circle(X['target_node'], yy, cfg['target_d'], face=col, edge='#ffffff', lw=cfg['target_edge_pt'],
                  gid=f'{gid}-target', z=4)
        is_hi = r['target_sym'] == hi
        sz = fs_hi if is_hi else fs_tgt
        page.text(X['target_label'], yy + cap_mm(sz) / 2, r['target_sym'], size=sz, fontstyle='italic',
                  weight='bold' if is_hi else None, gid=f'{gid}-targetlabel')
    # ---- nomination groups (left) and target classes (right)
    gl = cfg['group_labels']
    pad = cfg['bracket_pad']
    for g in cfg['group_order']:
        ys = [Y(n) for n in dn[dn['group'] == g]['source']]
        if not ys:
            continue
        mm.line(X['group_bracket'], min(ys) - pad, X['group_bracket'], max(ys) + pad, lw=cfg['group_bracket_pt'],
                colour=cfg['group_bracket_colour'], gid=f'{gid}-group')
        page.text(X['group_label'], (min(ys) + max(ys)) / 2 + cap_mm(fs_grp) / 2, gl.get(g, g), size=fs_grp,
                  ha='right', gid=f'{gid}-group')
    cl = cfg['class_labels']
    for c in cfg['class_order']:
        ys = [Y(n) for n in tn[tn['target_class'] == c]['target']]
        if not ys:
            continue
        mm.line(X['class_bracket'], min(ys) - pad, X['class_bracket'], max(ys) + pad, lw=cfg['class_bracket_pt'],
                colour=pal.get(c, '#808080'), gid=f'{gid}-class')
        page.text(X['class_label'], (min(ys) + max(ys)) / 2 + cap_mm(fs_cls) / 2, cl.get(c, c), size=fs_cls,
                  gid=f'{gid}-class')
    # ---- headers
    yh = y + cfg['y_head']
    page.text(X['dof_label'], yh, cfg['header_dofs'], size=fs_head, ha='right', weight='bold', gid=f'{gid}-header')
    if 'header_targets_x' in cfg:                          # centre of the target header (mm from the box left)
        xc = x + cfg['header_targets_x']
    else:
        tw = max(page.text_width(s, fs_tgt) for s in tn['target_sym'])
        xc = (X['target_node'] - cfg['target_d'] / 2 + X['target_label'] + tw) / 2
    page.text(xc, yh, cfg['header_targets'], size=fs_head_t, ha='center', weight='bold', gid=f'{gid}-header')
    return dict(type='network', layout=mode, n_dofs=int(len(dn)), n_targets=int(len(tn)), n_edges=int(len(ed)),
                n_highlight_edges=int(e['_hi'].sum()), highlight=hi, y_scale=round(scale, 4),
                nodes={n: round(Y(n), 3) for n in sorted(need)})
