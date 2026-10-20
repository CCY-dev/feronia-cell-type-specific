"""Figure layouts. Every position is in mm from the top-left corner; sizes come from config.yaml.

Layouts
  si_dot      Supplementary option A: a = UMAP grid (rows = genes, columns = WT | fer-4), b = dot plot by cluster.
  si_violin   Supplementary option B: a = UMAP grid (rows = genotype, columns = genes, aligned with b's columns),
              b = split violins by cluster.
  ed7_dot     Extended Data Fig. 7 d,e option A: d = YUC5 WT | fer-4 maps, e = dot plot by cluster (two row blocks).
  ed7_violin  Extended Data Fig. 7 d,e option B: same d, e = split violins by cluster (two row blocks).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import panels as P
from .colour import Gradient2, col2rgb
from .page import MM_PER_PT, Page


class Ctx:
    """Shared state for one build: config, cells, cluster rows, summaries, colour scales, names."""

    def __init__(self, cfg, cells, rows, summ):
        self.cfg, self.cells, self.rows, self.summ = cfg, cells, rows, summ
        st = cfg['style']
        self.gens = cfg['genotype_order']
        self.gname = cfg['names']['genotypes']
        self.names = cfg['names']
        self.frame = P.umap_frame(cells['umap_1'].to_numpy(float), cells['umap_2'].to_numpy(float),
                                  st['umap']['pad_frac'])
        self.gene_scales = {s: P.GeneScale(cells[g['column']].to_numpy(float), st['expression_colours'])
                            for s, g in cfg['genes'].items()}
        # dot colour = mean SCT per cluster; one scale shared by every dot plot (ED7 and both Supplementary figures)
        mode = st['dot']['mean']
        self.mean_col = {'arithmetic': 'mean_sct', 'seurat_dotplot': 'mean_seurat_dotplot',
                         'seurat_scaled': 'mean_seurat_scaled'}[mode]
        if mode == 'seurat_scaled':
            # Seurat DotPlot(scale = TRUE): per gene, z-score of mean(expm1(x)) over the plotted groups
            # (cluster x genotype), clipped to [col.min, col.max] = [-2.5, 2.5]
            avg = np.expm1(summ['mean_seurat_dotplot'].astype(float))
            z = avg.groupby(summ['gene']).transform(lambda v: (v - v.mean()) / v.std(ddof=1))
            summ['mean_seurat_scaled'] = z.clip(-2.5, 2.5)
            vmin, vmax = -2.5, 2.5
        else:
            vmin = 0.0
            vmax = st['dot']['mean_max'] or float(np.nanmax(summ[self.mean_col]))
        self.dot_vmin, self.dot_vmax = vmin, vmax
        self.dot_scale = Gradient2(st['expression_colours'][0], st['expression_colours'][1], vmin, vmax)
        self.dot_key_title = st['dot'].get('key_title', 'Mean expression')
        self.counts = {(r.cluster, r.genotype): int(r.n_cells)
                       for r in summ[summ['gene'] == summ['gene'].iloc[0]].itertuples()}
        self.vstyles = {g: dict(fill=P.rgb255(st['violin']['fill'][g]), edge=P.rgb255(st['violin']['edge'][g]))
                        for g in self.gens}
        self.palette = None
        if cfg.get('celltype_palette'):
            pal = pd.read_csv(cfg['celltype_palette'])
            pal = pal[pal['map'] == 'celltype']
            self.palette = {r['name']: tuple(c / 255 for c in col2rgb(r['hex'])) for _, r in pal.iterrows()}

    def gene_label(self, sym):
        return f"*{self.cfg['genes'][sym].get('label', sym)}*"


def _letter(page, x, y, s, st):
    page.text(x, y, s, size=st['letter_pt'], weight='bold', gid=f'letter-{s}')


def _map_h(ctx, w):
    xr, yr = ctx.frame
    return w * (yr[1] - yr[0]) / (xr[1] - xr[0])


def _gene_key_v(ctx, page, sym, x, y_bottom, geo, size):
    sc = ctx.gene_scales[sym]
    h = geo['key_h']
    ticks = P.nice_ticks(sc.vmin, sc.vmax, geo.get('key_step'))
    return P.colour_bar(page, x, y_bottom - h, geo['key_w'], h, sc.rgb01, sc.vmin, sc.vmax, ticks,
                        geo['ppi'], 'v', size, gid=f'key-{sym}')


def _gene_key_h(ctx, page, sym, x, y, w, geo, size):
    sc = ctx.gene_scales[sym]
    ticks = P.nice_ticks(sc.vmin, sc.vmax, geo.get('key_step'))
    return P.colour_bar(page, x, y, w, geo['key_w'], sc.rgb01, sc.vmin, sc.vmax, ticks, geo['ppi'], 'h', size,
                        gid=f'key-{sym}')


def _rows(ctx, subset=None):
    rows = ctx.rows if subset is None else [r for r in ctx.rows if r['celltype'] in subset]
    return rows


def _dot_keys(ctx, page, x, y, st):
    d = st['dot']
    size = st['small_pt']
    page.text(x, y, ctx.dot_key_title, size=size, gid='dot-colour-key')
    ticks = P.nice_ticks(ctx.dot_vmin, ctx.dot_vmax, d.get('mean_step'))
    w = d['colour_key_w']
    bottom = P.colour_bar(page, x, y + 1.3, w, st['umap']['key_w'], ctx.dot_scale.rgb01, ctx.dot_vmin, ctx.dot_vmax, ticks,
                          st['umap']['ppi'], 'h', size, gid='dot-colour-key')
    yb = P.size_key(page, x, bottom + 3.0, d['size_key'], d['dot_max_d'], size, 'Cells expressing')
    return yb


# ---- Supplementary, option A ---------------------------------------------------------------------------------------
def si_dot(ctx, fig):
    cfg, st = ctx.cfg, ctx.cfg['style']
    geo = dict(st['umap'], **fig.get('umap', {}))
    genes = fig['genes']
    W = fig['width_mm']
    page = Page(W, fig['height_mm'], st['font'], st['label_pt'])
    mw = geo['map_w']
    mh = _map_h(ctx, mw)
    x_maps, y_maps = fig['a_x'], fig['top'] + 7.0
    _letter(page, 0.0, fig['top'] + 2.9, 'a', st)
    for j, g in enumerate(ctx.gens):
        page.text(x_maps + j * (mw + geo['gap_x']) + mw / 2, y_maps - 2.0, ctx.gname[g], size=st['header_pt'],
                  ha='center', weight='bold', gid='a-genotype')
    for i, sym in enumerate(genes):
        y = y_maps + i * (mh + geo['gap_y'])
        col = cfg['genes'][sym]['column']
        for j, g in enumerate(ctx.gens):
            sub = ctx.cells[ctx.cells['genotype'] == g]
            P.draw_map(page, sub, sub[col].to_numpy(float), ctx.gene_scales[sym],
                       x_maps + j * (mw + geo['gap_x']), y, ctx.frame, geo, st['umap']['order'], gid=f'a-{sym}-{g}')
        page.text(x_maps - geo['row_label_dx'], y + mh / 2, ctx.gene_label(sym), size=st['header_pt'],
                  rotation=90, ha='center', gid='a-gene')
        kx = x_maps + len(ctx.gens) * (mw + geo['gap_x']) - geo['gap_x'] + geo['key_gap']
        _gene_key_v(ctx, page, sym, kx, y + mh - 1.0, geo, st['small_pt'])
    y_bot_maps = y_maps + len(genes) * (mh + geo['gap_y']) - geo['gap_y']
    P.umap_glyph(page, x_maps + 1.0, y_bot_maps + 1.2 + 4.5, size=st['small_pt'])
    # b
    rg = dict(st['rows'], **fig.get('rows', {}))
    xb = fig['b_x']
    _letter(page, xb - 2.0, fig['top'] + 2.9, 'b', st)
    header_y = fig['top'] + rg['header_h']
    rows = _rows(ctx)
    ys, spans, y_end = P.row_positions(rows, header_y + 1.6, rg['row_pitch'], rg['group_gap'], min_h=P.label_min_heights(ctx.names, st['label_pt']))
    dg = dict(st['dot'], **fig.get('dot', {}))
    dg['small_pt'] = st['small_pt']
    n_g = len(ctx.gens)
    chart_w = len(genes) * (n_g * dg['dot_pitch'] + dg['gene_gap']) - dg['gene_gap']
    x_chart = xb + rg['celltype_w'] + rg['stage_w'] + n_g * rg['count_w'] + rg['chart_gap']
    rg2 = dict(rg, label_pt=st['label_pt'], small_pt=st['small_pt'])
    P.draw_row_block(page, rows, ys, spans, xb, x_chart + chart_w + 0.8, rg2, ctx.names, ctx.counts, ctx.gens,
                     ctx.gname, header_y)
    x_end, _ = P.draw_dots(page, rows, ys, ctx.summ, genes, x_chart, dict(dg, gene_pt=st['header_pt'], mean_col=ctx.mean_col),
                           ctx.dot_scale, ctx.gens, ctx.gname, ctx.gene_label, header_y)
    _dot_keys(ctx, page, x_end + fig['key_dx'], header_y + 1.6, st)
    return page


# ---- Supplementary, option B ---------------------------------------------------------------------------------------
def si_violin(ctx, fig):
    cfg, st = ctx.cfg, ctx.cfg['style']
    genes = fig['genes']
    page = Page(fig['width_mm'], fig['height_mm'], st['font'], st['label_pt'])
    rg = dict(st['rows'], **fig.get('rows', {}))
    vg = dict(st['violin'], **fig.get('violin', {}))
    vg.update(small_pt=st['small_pt'], gene_pt=st['header_pt'], row_pitch=rg['row_pitch'],
              half_h=rg['row_pitch'] * st['violin']['half_frac'])
    n_g = len(ctx.gens)
    xb = fig['b_x']
    x_chart = xb + rg['celltype_w'] + rg['stage_w'] + n_g * rg['count_w'] + rg['chart_gap']
    cw, cg = vg['col_w'], vg['col_gap']
    # a: maps aligned with the violin columns
    geo = dict(st['umap'], map_w=cw, **fig.get('umap', {}))
    mh = _map_h(ctx, cw)
    y_maps = fig['top'] + 7.0
    _letter(page, 0.0, fig['top'] + 2.9, 'a', st)
    for i, sym in enumerate(genes):
        x = x_chart + i * (cw + cg)
        col = cfg['genes'][sym]['column']
        page.text(x + cw / 2, y_maps - 2.0, ctx.gene_label(sym), size=st['header_pt'], ha='center', gid='a-gene')
        for j, g in enumerate(ctx.gens):
            sub = ctx.cells[ctx.cells['genotype'] == g]
            P.draw_map(page, sub, sub[col].to_numpy(float), ctx.gene_scales[sym], x,
                       y_maps + j * (mh + geo['gap_y']), ctx.frame, geo, st['umap']['order'], gid=f'a-{sym}-{g}')
        ky = y_maps + n_g * (mh + geo['gap_y']) - geo['gap_y'] + 1.2
        kb = _gene_key_h(ctx, page, sym, x + (cw - geo['key_h_w']) / 2, ky, geo['key_h_w'], geo, st['small_pt'])
    for j, g in enumerate(ctx.gens):
        page.text(x_chart - geo['row_label_dx'], y_maps + j * (mh + geo['gap_y']) + mh / 2, ctx.gname[g],
                  size=st['header_pt'], weight='bold', rotation=90, ha='center', gid='a-genotype')
    y_a_bottom = kb
    page.text(x_chart - geo['row_label_dx'] - 1.0, ky + 1.5, 'Expression', size=st['small_pt'], ha='right',
              gid='a-key-title')
    # optional reference map (cell types, Fig. 1a colours) in the space left of the maps
    if fig.get('reference_map') and ctx.palette:
        rw = fig['reference_map_w']
        rx = x_chart - geo['row_label_dx'] - 4.5 - rw
        sub = ctx.cells
        xr, yr = ctx.frame
        w_px, h_px, w, h = P.map_size(rw, xr, yr, geo['ppi'])
        rgb = np.array([ctx.palette.get(c, (0.85, 0.85, 0.85)) for c in sub['celltype.anno.Li']])
        arr = P.render_map(sub['umap_1'].to_numpy(float), sub['umap_2'].to_numpy(float), rgb, xr, yr, w_px, h_px,
                           geo['ppi'], geo['point_size_pt2'] * (rw / 28.0) ** 2)
        page.image(arr, rx, y_maps, w, h, gid='a-reference')
        page.text(rx + w / 2, y_maps - 2.0, 'Cell type (WT and *fer-4*)', size=st['small_pt'], ha='center',
                  gid='a-reference')
        P.umap_glyph(page, rx + 1.0, y_maps + h + 1.2 + 4.5, size=st['small_pt'])
    else:
        P.umap_glyph(page, xb + 1.0, y_maps + mh * 2 + geo['gap_y'] - 1.0, size=st['small_pt'])
    # b
    top_b = y_a_bottom + fig['b_gap']
    _letter(page, 0.0, top_b + 2.9, 'b', st)
    header_y = top_b + rg['header_h']
    rows = _rows(ctx)
    ys, spans, y_end = P.row_positions(rows, header_y + 1.6, rg['row_pitch'], rg['group_gap'], min_h=P.label_min_heights(ctx.names, st['label_pt']))
    chart_w = len(genes) * (cw + cg) - cg
    rg2 = dict(rg, label_pt=st['label_pt'], small_pt=st['small_pt'])
    if ctx.palette and fig.get('reference_map'):
        rg2['swatch'] = ctx.palette
    P.draw_row_block(page, rows, ys, spans, xb, x_chart + chart_w + 0.8, rg2, ctx.names, ctx.counts, ctx.gens,
                     ctx.gname, header_y)
    x_end, y_lab, _ = P.draw_violins(page, rows, ys, ctx.cells, genes, cfg['genes'], x_chart, vg, ctx.gens,
                                     ctx.gene_label, header_y, ctx.vstyles)
    page.text((x_chart + x_end) / 2, y_lab, st['violin']['axis_title'], size=st['label_pt'], ha='center',
              gid='vln-axis-title')
    P.split_key(page, xb + 1.0, top_b + 2.2, ctx.vstyles, ctx.gens, ctx.gname, st['small_pt'])
    return page


# ---- Extended Data Fig. 7 d,e ------------------------------------------------------------------------------------
def _ed7_d(ctx, page, fig, st):
    sym = fig['genes'][0]
    geo = dict(st['umap'], **fig.get('umap', {}))
    mw = geo['map_w']
    mh = _map_h(ctx, mw)
    col = ctx.cfg['genes'][sym]['column']
    x0, y0 = fig['d_x'], fig['top'] + 7.0
    _letter(page, 0.0, fig['top'] + 2.9, 'd', st)
    for j, g in enumerate(ctx.gens):
        x = x0 + j * (mw + geo['gap_x'])
        sub = ctx.cells[ctx.cells['genotype'] == g]
        P.draw_map(page, sub, sub[col].to_numpy(float), ctx.gene_scales[sym], x, y0, ctx.frame, geo,
                   st['umap']['order'], gid=f'd-{sym}-{g}')
        page.text(x + mw / 2, y0 - 2.0, ctx.gname[g], size=st['header_pt'], ha='center', weight='bold',
                  gid='d-genotype')
    kx = x0 + 2 * mw + geo['gap_x'] + geo['key_gap']
    page.text(kx, y0 + mh - 1.0 - geo['key_h'] - 1.6, ctx.gene_label(sym), size=st['label_pt'], gid='d-key-title')
    right = _gene_key_v(ctx, page, sym, kx, y0 + mh - 1.0, geo, st['small_pt'])
    page.text(x0 - geo['row_label_dx'], y0 + mh / 2, f"{ctx.gene_label(sym)} expression", size=st['header_pt'],
              rotation=90, ha='center', gid='d-rowlabel')
    yb = P.umap_glyph(page, x0 + 1.0, y0 + mh + 1.2 + 4.5, size=st['small_pt'])
    return yb, right


def _ed7_blocks(ctx, fig):
    split = fig['split_before']
    rows = ctx.rows
    k = next(i for i, r in enumerate(rows) if r['celltype'] == split)
    return [rows[:k], rows[k:]]


def ed7(ctx, fig, kind):
    st = ctx.cfg['style']
    sym = fig['genes'][0]
    page = Page(fig['width_mm'], fig['height_mm'], st['font'], st['label_pt'])
    y_d_bottom, d_right = _ed7_d(ctx, page, fig, st)
    rg = dict(st['rows'], **fig.get('rows', {}))
    rg2 = dict(rg, label_pt=st['label_pt'], small_pt=st['small_pt'])
    n_g = len(ctx.gens)
    top_e = y_d_bottom + fig['e_gap']
    _letter(page, 0.0, top_e + 2.9, 'e', st)
    header_y = top_e + rg['header_h']
    lab_w = rg['celltype_w'] + rg['stage_w'] + n_g * rg['count_w'] + rg['chart_gap']
    x = fig['e_x']
    y_lab_max = 0
    for bi, rows in enumerate(_ed7_blocks(ctx, fig)):
        ys, spans, y_end = P.row_positions(rows, header_y + 1.6, rg['row_pitch'], rg['group_gap'], min_h=P.label_min_heights(ctx.names, st['label_pt']))
        if kind == 'dot':
            dg = dict(st['dot'], **fig.get('dot', {}))
            dg.update(small_pt=st['small_pt'], gene_pt=st['header_pt'], mean_col=ctx.mean_col)
            chart_w = n_g * dg['dot_pitch']
            P.draw_row_block(page, rows, ys, spans, x, x + lab_w + chart_w + 0.8, rg2, ctx.names, ctx.counts,
                             ctx.gens, ctx.gname, header_y)
            x_end, _ = P.draw_dots(page, rows, ys, ctx.summ, [sym], x + lab_w, dg, ctx.dot_scale, ctx.gens,
                                   ctx.gname, ctx.gene_label, header_y)
        else:
            vg = dict(st['violin'], **fig.get('violin', {}))
            vg.update(small_pt=st['small_pt'], gene_pt=st['header_pt'], row_pitch=rg['row_pitch'],
                      half_h=rg['row_pitch'] * st['violin']['half_frac'])
            chart_w = vg['col_w']
            P.draw_row_block(page, rows, ys, spans, x, x + lab_w + chart_w + 0.8, rg2, ctx.names, ctx.counts,
                             ctx.gens, ctx.gname, header_y)
            x_end, y_lab, _ = P.draw_violins(page, rows, ys, ctx.cells, [sym], ctx.cfg['genes'], x + lab_w, vg,
                                             ctx.gens, ctx.gene_label, header_y, ctx.vstyles)
            page.text(x + lab_w + chart_w / 2, y_lab, st['violin']['axis_title'], size=st['small_pt'],
                      ha='center', gid='vln-axis-title')
        x = x_end + fig['block_gap']
    if kind == 'dot':
        _dot_keys(ctx, page, fig['key_x'], fig['top'] + 7.0, st)
    else:
        P.split_key(page, fig['key_x'], fig['top'] + 9.0, ctx.vstyles, ctx.gens, ctx.gname, st['small_pt'])
    return page


# ---- option A grouped by genotype ----------------------------------------------------------------------------------
def _geno_geo(ctx, fig):
    st = ctx.cfg['style']
    dg = {**st['dot'], **st.get('dot_by_genotype', {}), **fig.get('dot', {})}
    dg.update(small_pt=st['small_pt'], gene_pt=st['header_pt'], mean_col=ctx.mean_col,
              band_grey=st['rows'].get('band_grey'))
    return dg


def si_dot_geno(ctx, fig):
    """As si_dot, but panel b is grouped by genotype: WT block | fer-4 block, each with its own cell counts."""
    cfg, st = ctx.cfg, ctx.cfg['style']
    geo = dict(st['umap'], **fig.get('umap', {}))
    genes = fig['genes']
    page = Page(fig['width_mm'], fig['height_mm'], st['font'], st['label_pt'])
    mw = geo['map_w']; mh = _map_h(ctx, mw)
    x_maps, y_maps = fig['a_x'], fig['top'] + 7.0
    _letter(page, 0.0, fig['top'] + 2.9, 'a', st)
    for j, g in enumerate(ctx.gens):
        page.text(x_maps + j * (mw + geo['gap_x']) + mw / 2, y_maps - 2.0, ctx.gname[g], size=st['header_pt'],
                  ha='center', weight='bold', gid='a-genotype')
    for i, sym in enumerate(genes):
        y = y_maps + i * (mh + geo['gap_y'])
        col = cfg['genes'][sym]['column']
        for j, g in enumerate(ctx.gens):
            sub = ctx.cells[ctx.cells['genotype'] == g]
            P.draw_map(page, sub, sub[col].to_numpy(float), ctx.gene_scales[sym],
                       x_maps + j * (mw + geo['gap_x']), y, ctx.frame, geo, st['umap']['order'], gid=f'a-{sym}-{g}')
        page.text(x_maps - geo['row_label_dx'], y + mh / 2, ctx.gene_label(sym), size=st['header_pt'],
                  rotation=90, ha='center', gid='a-gene')
        kx = x_maps + len(ctx.gens) * (mw + geo['gap_x']) - geo['gap_x'] + geo['key_gap']
        _gene_key_v(ctx, page, sym, kx, y + mh - 1.0, geo, st['small_pt'])
    y_bot_maps = y_maps + len(genes) * (mh + geo['gap_y']) - geo['gap_y']
    P.umap_glyph(page, x_maps + 1.0, y_bot_maps + 1.2 + 4.5, size=st['small_pt'])
    # b
    rg = dict(st['rows'], **fig.get('rows', {}))
    dg = _geno_geo(ctx, fig)
    xb = fig['b_x']
    _letter(page, xb - 2.0, fig['top'] + 2.9, 'b', st)
    # genotype header baseline level with panel a's genotype headers
    header_y = (y_maps - 2.0) + dg['geno_dy'] + dg['col_label_h']
    rows = _rows(ctx)
    ys, spans, y_end = P.row_positions(rows, header_y + 1.6, rg['row_pitch'], rg['group_gap'], min_h=P.label_min_heights(ctx.names, st['label_pt']))
    x_chart = xb + rg['celltype_w'] + rg['stage_w'] + rg['chart_gap']
    rg2 = dict(rg, label_pt=st['label_pt'], small_pt=st['small_pt'], show_counts=False, band_grey=None)
    P.draw_row_block(page, rows, ys, spans, xb, x_chart, rg2, ctx.names, ctx.counts, ctx.gens, ctx.gname, header_y)
    x_end, _ = P.draw_genotype_blocks(page, rows, ys, spans, x_chart, ctx.summ, genes, dg, ctx.dot_scale, ctx.gens,
                                      ctx.gname, ctx.gene_label, ctx.counts, header_y, band_x0=xb)
    kx = fig.get('key_x') or (x_end + fig['key_dx'])
    _dot_keys(ctx, page, kx, fig.get('key_y', header_y + 1.6), st)
    return page


def ed7_geno(ctx, fig):
    """ED7 d,e option A grouped by genotype: e = two row blocks, each [WT: cells, YUC5 | fer-4: cells, YUC5]."""
    st = ctx.cfg['style']
    sym = fig['genes'][0]
    page = Page(fig['width_mm'], fig['height_mm'], st['font'], st['label_pt'])
    y_d_bottom, d_right = _ed7_d(ctx, page, fig, st)
    rg = dict(st['rows'], **fig.get('rows', {}))
    dg = _geno_geo(ctx, fig)
    top_e = y_d_bottom + fig['e_gap']
    _letter(page, 0.0, top_e + 2.9, 'e', st)
    header_y = top_e + 2.9 + dg['geno_dy'] + dg['col_label_h']
    rg2 = dict(rg, label_pt=st['label_pt'], small_pt=st['small_pt'], show_counts=False, band_grey=None)
    lab_w = rg['celltype_w'] + rg['stage_w'] + rg['chart_gap']
    x = fig['e_x']
    for rows in _ed7_blocks(ctx, fig):
        ys, spans, y_end = P.row_positions(rows, header_y + 1.6, rg['row_pitch'], rg['group_gap'], min_h=P.label_min_heights(ctx.names, st['label_pt']))
        P.draw_row_block(page, rows, ys, spans, x, x + lab_w, rg2, ctx.names, ctx.counts, ctx.gens, ctx.gname,
                         header_y)
        x_end, _ = P.draw_genotype_blocks(page, rows, ys, spans, x + lab_w, ctx.summ, [sym], dg, ctx.dot_scale,
                                          ctx.gens, ctx.gname, ctx.gene_label, ctx.counts, header_y, band_x0=x,
                                          show_gene_labels=dg.get('gene_labels', True))
        x = x_end + fig['block_gap']
    _dot_keys(ctx, page, fig['key_x'], fig['top'] + 7.0, st)
    return page


# ---- Extended Data Fig. 7d and 7e sized to the deck slots ------------------------------------------------------------
def ed7_d_slot(ctx, fig):
    """d only: WT | fer-4 YUC5 maps with the gene key, sized to the current ED7 d slot."""
    st = ctx.cfg['style']
    page = Page(fig['width_mm'], fig['height_mm'], st['font'], st['label_pt'])
    _ed7_d(ctx, page, fig, st)
    return page


def _stage_short(ctx, stage):
    return ctx.cfg['names'].get('stages_short', {}).get(stage, ctx.names['stages'].get(stage, stage))


def cluster_strip(ctx, page, genes, top, sg, x0, x1, gene_rows):
    """Clusters as columns grouped by cell type (wrapped cell-type names over a bracket, rotated stage names), one
    dot row per (genotype, gene). gene_rows=False: one row per genotype labelled with the genotype (single gene).
    gene_rows=True: a block of gene rows per genotype, genotype name rotated beside its block.
    Returns (y_bottom_of_dots, x_end, pitch)."""
    st = ctx.cfg['style']
    lp, sp = sg['label_pt'], sg['stage_pt']
    lead = lp * 0.3528 * 1.15
    rows = _rows(ctx)
    groups = []
    for r in rows:
        if not groups or groups[-1][0] != r['celltype']:
            groups.append((r['celltype'], []))
        groups[-1][1].append(r)
    lines = {ct: sg.get('celltype_lines', {}).get(ct, ctx.names['celltype_lines'].get(ct, [ctx.names['celltypes'].get(ct, ct)]))
             for ct, _ in groups}
    need = {ct: max(page.text_width(l, lp) for l in lines[ct]) + sg['label_pad'] for ct, _ in groups}
    avail = x1 - x0 - (len(groups) - 1) * sg['group_gap']
    lo, hi = 1.0, sg['pitch_max']
    for _ in range(40):
        mid = (lo + hi) / 2
        tot = sum(max(len(rr) * mid, need[ct]) for ct, rr in groups)
        lo, hi = (mid, hi) if tot <= avail else (lo, mid)
    pitch = lo
    widths = [max(len(rr) * pitch, need[ct]) for ct, rr in groups]
    n_lines = max(len(v) for v in lines.values())
    y_ct_last = top + sg['ct_top'] + (n_lines - 1) * lead
    y_br = y_ct_last + sg['bracket_dy']
    stage_w = max(page.text_width(_stage_short(ctx, r['stage']), sp) for r in rows)
    y_st_top = y_br + sg['stage_dy']
    y = y_st_top + stage_w + sg['dots_dy']
    row_y = {}                                            # (genotype, gene) -> y
    for j, g in enumerate(ctx.gens):
        if gene_rows and j:
            y += sg['block_gap']
        y_block0 = y
        for sym in (genes if gene_rows else genes[:1]):
            row_y[(g, sym)] = y
            y += sg['row_pitch']
        if gene_rows:
            yc = (y_block0 + y - sg['row_pitch']) / 2
            page.text(x0 - sg['geno_label_dx'], yc, ctx.gname[g], size=st['header_pt'], weight='bold', rotation=90,
                      ha='center', gid='strip-genotype')
            for sym in genes:
                page.text(x0 - sg['row_label_gap'], row_y[(g, sym)], ctx.gene_label(sym), size=sg['gene_pt'],
                          ha='right', va='center', gid='strip-gene')
        else:
            page.text(x0 - sg['row_label_gap'], row_y[(g, genes[0])], ctx.gname[g], size=st['header_pt'],
                      weight='bold', ha='right', va='center', gid='strip-genotype')
    y_last = y - sg['row_pitch']
    band_t, band_b = y_br + 0.5, y_last + sg['row_pitch'] / 2
    s = ctx.summ.set_index(['gene', 'cluster', 'genotype'])
    x = x0
    for k, ((ct, rr), w) in enumerate(zip(groups, widths)):
        if k % 2 == 0 and st['rows'].get('band_grey'):
            page.rect(x - sg['group_gap'] / 2, band_t, w + sg['group_gap'], band_b - band_t,
                      fill=P.grey(st['rows']['band_grey']), gid='strip-band', zorder=0)
        for i, ln in enumerate(lines[ct]):
            page.text(x + w / 2, y_ct_last - (len(lines[ct]) - 1 - i) * lead, ln, size=lp, ha='center',
                      gid='strip-celltype')
        page.line(x + 0.3, y_br, x + w - 0.3, y_br, lw=0.4, colour=P.grey(0.3), gid='strip-bracket')
        cx0 = x + (w - len(rr) * pitch) / 2
        for i, r in enumerate(rr):
            xc = cx0 + (i + 0.5) * pitch
            page.text(xc + 0.716 * sp * 0.3528 / 2, y_st_top, _stage_short(ctx, r['stage']), size=sp, ha='right',
                      rotation=90, gid='strip-stage')
            for (g, sym), yy in row_y.items():
                rec = s.loc[(sym, r['cluster'], g)]
                if rec['n_cells'] > 0 and rec['pct_expressing'] > 0:
                    page.dots(np.array([xc]), np.array([yy]), P.dot_diameter([rec['pct_expressing']], sg['dot_max_d']),
                              ctx.dot_scale.rgb01(np.array([rec[ctx.mean_col]])),
                              edge=P.grey(st['dot']['dot_edge_grey']), lw=st['dot']['dot_edge_pt'], gid=f'strip-dot-{g}')
        x += w + sg['group_gap']
    return band_b, x - sg['group_gap'], pitch


def _strip_keys(ctx, page, x0, yk, sg):
    st = ctx.cfg['style']; d = st['dot']
    page.text(x0, yk, ctx.dot_key_title, size=st['small_pt'], gid='dot-colour-key')
    ticks = P.nice_ticks(ctx.dot_vmin, ctx.dot_vmax, d.get('mean_step'))
    P.colour_bar(page, x0 + page.text_width(ctx.dot_key_title, st['small_pt']) + 1.5, yk - 1.6, d['colour_key_w'],
                 st['umap']['key_w'], ctx.dot_scale.rgb01, ctx.dot_vmin, ctx.dot_vmax, ticks, st['umap']['ppi'], 'h',
                 st['small_pt'], gid='dot-colour-key')
    xs = x0 + sg['size_key_dx']
    page.text(xs, yk, 'Cells expressing', size=st['small_pt'], gid='dot-size-key')
    xs += page.text_width('Cells expressing', st['small_pt']) + 2.0
    for pct in d['size_key']:
        dd = float(P.dot_diameter([pct], sg['dot_max_d'])[0])
        page.dots(np.array([xs + 1.2]), np.array([yk - 0.7]), np.array([dd]), np.array([[0.6, 0.6, 0.6]]),
                  edge=P.grey(d['dot_edge_grey']), lw=d['dot_edge_pt'], gid='dot-size-key')
        lab = f'{pct}%'
        page.text(xs + 2.8, yk, lab, size=st['small_pt'], gid='dot-size-key')
        xs += 2.8 + page.text_width(lab, st['small_pt']) + 2.2


def ed7_e_strip(ctx, fig):
    """e as a full-width strip: clusters as columns grouped by cell type; one row of dots per genotype."""
    st = ctx.cfg['style']
    sg = fig['strip']
    page = Page(fig['width_mm'], fig['height_mm'], st['font'], st['label_pt'])
    _letter(page, 0.0, fig['top'] + 2.9, 'e', st)
    yb, _, _ = cluster_strip(ctx, page, fig['genes'], fig['top'], sg, sg['x_dots'],
                             fig['width_mm'] - sg['right_margin'], gene_rows=False)
    _strip_keys(ctx, page, sg['x_dots'], yb + sg['key_dy'], sg)
    return page


def si_dof_strip(ctx, fig):
    """All ten DOFs. a = one row of maps per genotype (genes as columns) with a labelled key per gene;
    b = cluster strip with a block of gene rows per genotype. Shared left label column."""
    cfg, st = ctx.cfg, ctx.cfg['style']
    genes = fig['genes']
    geo = dict(st['umap'], **fig.get('umap', {}))
    x0, x1 = fig['x_left'], fig['width_mm'] - fig['right_margin']
    n = len(genes)
    mw = (x1 - x0 - (n - 1) * geo['gap_x']) / n
    geo['map_w'] = mw
    mh = _map_h(ctx, mw)
    page = Page(fig['width_mm'], fig['height_mm'], st['font'], st['label_pt'])
    _letter(page, 0.0, fig['top'] + 2.9, 'a', st)
    y_maps = fig['top'] + 6.0
    for i, sym in enumerate(genes):
        x = x0 + i * (mw + geo['gap_x'])
        col = cfg['genes'][sym]['column']
        page.text(x + mw / 2, y_maps - 1.5, ctx.gene_label(sym), size=st['header_pt'], ha='center', gid='a-gene')
        for j, g in enumerate(ctx.gens):
            sub = ctx.cells[ctx.cells['genotype'] == g]
            P.draw_map(page, sub, sub[col].to_numpy(float), ctx.gene_scales[sym], x,
                       y_maps + j * (mh + geo['gap_y']), ctx.frame, geo, st['umap']['order'], gid=f'a-{sym}-{g}')
        ky = y_maps + len(ctx.gens) * (mh + geo['gap_y']) - geo['gap_y'] + fig['key_dy']
        kw = mw * fig['key_frac']
        kb = _gene_key_h(ctx, page, sym, x + (mw - kw) / 2, ky, kw, geo, st['small_pt'])
    for j, g in enumerate(ctx.gens):
        page.text(x0 - fig['geno_label_gap'], y_maps + j * (mh + geo['gap_y']) + mh * fig.get('geno_label_frac', 0.5), ctx.gname[g],
                  size=st['header_pt'], weight='bold', ha='right', va='center', gid='a-genotype')
    # colour-key title in the label column, level with the keys
    for k, t in enumerate(fig['key_title']):
        page.text(x0 - fig['geno_label_gap'], ky + 1.0 + k * 2.2, t, size=st['small_pt'], ha='right', gid='a-key-title')
    # UMAP axes glyph tucked under the label column, against the first map's lower-left corner
    P.umap_glyph(page, fig['glyph_x'], y_maps + 2 * mh + geo['gap_y'] - fig['glyph_up'], L=fig['glyph_L'],
                 size=st['small_pt'])
    # b
    top_b = kb + fig['b_gap']
    _letter(page, 0.0, top_b + 2.9, 'b', st)
    sg = dict(fig['strip'])
    yb, _, _ = cluster_strip(ctx, page, genes, top_b, sg, x0, x1, gene_rows=True)
    _strip_keys(ctx, page, x0, yb + sg['key_dy'], sg)
    return page


# ---- one Supplementary figure for all ten DOFs ----------------------------------------------------------------------
def si_dof_combined(ctx, fig):
    """a = maps in blocks of n genes (columns) x genotype (rows) with one key per gene; b = dot plot by genotype."""
    cfg, st = ctx.cfg, ctx.cfg['style']
    genes = fig['genes']
    per = fig['genes_per_block']
    geo = dict(st['umap'], **fig.get('umap', {}))
    mw = geo['map_w']; mh = _map_h(ctx, mw)
    page = Page(fig['width_mm'], fig['height_mm'], st['font'], st['label_pt'])
    _letter(page, 0.0, fig['top'] + 2.9, 'a', st)
    xa = fig['a_x']
    y = fig['top'] + 5.6
    for b0 in range(0, len(genes), per):
        blk = genes[b0:b0 + per]
        y_maps = y + 2.0
        for i, sym in enumerate(blk):
            x = xa + i * (mw + geo['gap_x'])
            col = cfg['genes'][sym]['column']
            page.text(x + mw / 2, y_maps - 1.6, ctx.gene_label(sym), size=st['header_pt'], ha='center', gid='a-gene')
            for j, g in enumerate(ctx.gens):
                sub = ctx.cells[ctx.cells['genotype'] == g]
                P.draw_map(page, sub, sub[col].to_numpy(float), ctx.gene_scales[sym], x,
                           y_maps + j * (mh + geo['gap_y']), ctx.frame, geo, st['umap']['order'], gid=f'a-{sym}-{g}')
            ky = y_maps + len(ctx.gens) * (mh + geo['gap_y']) - geo['gap_y'] + 1.0
            kb = _gene_key_h(ctx, page, sym, x + (mw - geo['key_h_w']) / 2, ky, geo['key_h_w'], geo, st['small_pt'])
        for j, g in enumerate(ctx.gens):
            page.text(xa - geo['row_label_dx'], y_maps + j * (mh + geo['gap_y']) + mh / 2, ctx.gname[g],
                      size=st['header_pt'], weight='bold', rotation=90, ha='center', gid='a-genotype')
        y = kb + fig['block_gap']
    xg = xa + per * (mw + geo['gap_x']) - geo['gap_x'] + 3.0
    P.umap_glyph(page, xg, fig['top'] + 5.6 + 2.0 + 2 * mh + geo['gap_y'] - 1.0, size=st['small_pt'])
    # b
    top_b = y - fig['block_gap'] + fig['b_gap']
    rg = dict(st['rows'], **fig.get('rows', {}))
    dg = _geno_geo(ctx, fig)
    xb = fig['b_x']
    _letter(page, 0.0, top_b + 2.9, 'b', st)
    header_y = top_b + 2.9 + dg['geno_dy'] + dg['col_label_h']
    rows = _rows(ctx)
    ys, spans, y_end = P.row_positions(rows, header_y + 1.6, rg['row_pitch'], rg['group_gap'],
                                       min_h=P.label_min_heights(ctx.names, st['label_pt']))
    x_chart = xb + rg['celltype_w'] + rg['stage_w'] + rg['chart_gap']
    rg2 = dict(rg, label_pt=st['label_pt'], small_pt=st['small_pt'], show_counts=False, band_grey=None)
    P.draw_row_block(page, rows, ys, spans, xb, x_chart, rg2, ctx.names, ctx.counts, ctx.gens, ctx.gname, header_y)
    x_end, _ = P.draw_genotype_blocks(page, rows, ys, spans, x_chart, ctx.summ, genes, dg, ctx.dot_scale, ctx.gens,
                                      ctx.gname, ctx.gene_label, ctx.counts, header_y, band_x0=xb)
    _dot_keys(ctx, page, x_end + fig['key_dx'], header_y + 1.6, st)
    return page


# ---- DE heatmap (gene x cluster), same cluster columns as the strip ------------------------------------------------
def _column_header(ctx, page, top, sg, x0, x1, reserve_right=0.0, draw=True):
    """Cell-type names over brackets and rotated stage names; returns (cols [(cluster, xc)], pitch, groups, widths,
    y_below_header)."""
    lp, sp = sg['label_pt'], sg['stage_pt']
    lead = lp * 0.3528 * 1.15
    rows = _rows(ctx)
    groups = []
    for r in rows:
        if not groups or groups[-1][0] != r['celltype']:
            groups.append((r['celltype'], []))
        groups[-1][1].append(r)
    lines = {ct: sg.get('celltype_lines', {}).get(ct, ctx.names['celltype_lines'].get(ct, [ctx.names['celltypes'].get(ct, ct)]))
             for ct, _ in groups}
    need = {ct: max(page.text_width(l, lp) for l in lines[ct]) + sg['label_pad'] for ct, _ in groups}
    avail = x1 - reserve_right - x0 - (len(groups) - 1) * sg['group_gap']
    lo, hi = 1.0, sg['pitch_max']
    for _ in range(40):
        mid = (lo + hi) / 2
        tot = sum(max(len(rr) * mid, need[ct]) for ct, rr in groups)
        lo, hi = (mid, hi) if tot <= avail else (lo, mid)
    pitch = lo
    widths = [max(len(rr) * pitch, need[ct]) for ct, rr in groups]
    n_lines = max(len(v) for v in lines.values())
    y_ct_last = top + sg['ct_top'] + (n_lines - 1) * lead
    y_br = y_ct_last + sg['bracket_dy']
    stage_w = max(page.text_width(_stage_short(ctx, r['stage']), sp) for r in rows)
    y_st_top = y_br + sg['stage_dy']
    cols, x = [], x0
    for (ct, rr), w in zip(groups, widths):
        if draw:
            for i, ln in enumerate(lines[ct]):
                page.text(x + w / 2, y_ct_last - (len(lines[ct]) - 1 - i) * lead, ln, size=lp, ha='center', gid='hm-celltype')
            page.line(x + 0.3, y_br, x + w - 0.3, y_br, lw=0.4, colour=P.grey(0.3), gid='hm-bracket')
        cx0 = x + (w - len(rr) * pitch) / 2
        for i, r in enumerate(rr):
            xc = cx0 + (i + 0.5) * pitch
            if draw:
                page.text(xc + 0.716 * sp * 0.3528 / 2, y_st_top, _stage_short(ctx, r['stage']), size=sp, ha='right',
                          rotation=90, gid='hm-stage')
            cols.append((r['cluster'], xc, ct))
        x += w + sg['group_gap']
    return cols, pitch, groups, widths, y_st_top + stage_w + sg['dots_dy'] - sg['row_pitch'] / 2


def _diverging(c_neg, c_pos, lim):
    neg = Gradient2('white', c_neg, 0.0, lim)
    pos = Gradient2('white', c_pos, 0.0, lim)
    def rgb01(v):
        v = np.clip(np.asarray(v, float), -lim, lim)
        out = np.where(v[..., None] < 0, neg.rgb01(-v), pos.rgb01(v))
        return out
    return rgb01


def draw_de_heatmap(ctx, page, fig, top, x0, x1, letter=None):
    st = ctx.cfg['style']; hg = fig['heatmap']; sg = dict(fig['strip'])
    de = pd.read_csv(ctx.cfg['de_table'])
    genes = fig['genes'] + (fig.get('extra_genes') or [])
    if letter:
        _letter(page, 0.0, top + 2.9, letter, st)
    cols, pitch, groups, widths, y0 = _column_header(ctx, page, top, sg, x0, x1, reserve_right=hg['count_w'])
    rp = hg['row_pitch']
    lim = hg['lim']
    rgb = _diverging(hg['down'], hg['up'], lim)
    key = de.set_index(['symbol', 'cluster'])
    y = y0
    row_y = {}
    for g in genes:
        if g in (fig.get('extra_genes') or []) and g == genes[len(fig['genes'])]:
            y += hg['extra_gap']
        row_y[g] = y
        y += rp
    y_end = y
    # bands per cell-type group (as in panel b) behind the cells
    x = x0
    for k, ((ct, rr), w) in enumerate(zip(groups, widths)):
        if k % 2 == 0 and st['rows'].get('band_grey'):
            page.rect(x - sg['group_gap'] / 2, y0 - 0.6, w + sg['group_gap'], y_end - y0 + 1.2,
                      fill=P.grey(st['rows']['band_grey']), gid='hm-band', zorder=0)
        x += w + sg['group_gap']
    cw, ch = pitch - hg['cell_gap'], rp - hg['cell_gap']
    n_de = {}
    for g in genes:
        yy = row_y[g]
        page.text(x0 - 1.2, yy + rp / 2, ctx.gene_label(g), size=hg['gene_pt'], ha='right', va='center', gid='hm-gene')
        n_de[g] = 0
        for cl, xc, ct in cols:
            xl, yt = xc - cw / 2, yy + hg['cell_gap'] / 2
            if (g, cl) not in key.index:                       # not tested by muscat
                page.rect(xl, yt, cw, ch, fill='white', edge=P.grey(0.75), lw=0.25, gid='hm-nottested', zorder=2)
                continue
            rec = key.loc[(g, cl)]
            if not bool(rec['DE']):
                page.rect(xl, yt, cw, ch, fill=P.grey(hg['ns_grey']), gid='hm-ns', zorder=2)
                continue
            n_de[g] += 1
            lf = float(rec['logFC'])
            c = tuple(rgb(np.array([lf]))[0])
            page.rect(xl, yt, cw, ch, fill=c, edge=P.grey(0.15), lw=0.3, gid='hm-de', zorder=3)
            fc = 2 ** abs(lf)
            lab = ('\u2191' if lf > 0 else '\u2193') + f"{fc:.1f}"      # direction printed too (greyscale/CVD safe)
            if page.text_width(lab, hg['value_pt']) > cw - 0.15:
                raise ValueError(f'heatmap value {lab} does not fit its cell')
            lum = 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]
            page.text(xc, yt + ch / 2, lab, size=hg['value_pt'], ha='center', va='center',
                      colour='white' if lum < 0.5 else 'k', gid='hm-value')
        if hg.get('count_w', 0) > 0:
            page.text(x1, yy + rp / 2, str(n_de[g]), size=hg['gene_pt'], ha='right', va='center', gid='hm-count')
    # key row
    yk = y_end + hg['key_dy']
    sz = st['small_pt']
    t1 = 'Fold change, fer-4 vs WT (DE only):'
    t1 = 'Fold change, *fer-4* vs WT (DE only):'
    page.text(x0, yk, t1, size=sz, gid='hm-key')
    xk = x0 + page.text_width(t1, sz) + 1.5
    page.text(xk, yk, 'down', size=sz, gid='hm-key')
    xk += page.text_width('down', sz) + 1.0
    ticks = [-1, 0, 1]
    P.colour_bar(page, xk, yk - 1.6, hg['key_w'], st['umap']['key_w'], rgb, -lim, lim, ticks, st['umap']['ppi'], 'h', sz,
                 gid='hm-key')
    page.text(xk + hg['key_w'] + 1.0, yk, 'up', size=sz, gid='hm-key')
    xs = xk + hg['key_w'] + 1.0 + page.text_width('up', sz) + 4.0
    page.text(xs, yk, 'log2 scale; \u2191/\u2193 number = fold change up/down', size=sz, colour=P.grey(0.25), gid='hm-key')
    xs += page.text_width('log2 scale; \u2191/\u2193 number = fold change up/down', sz) + 4.0
    page.rect(xs, yk - 1.9, 2.6, 2.2, fill=P.grey(hg['ns_grey']), gid='hm-key')
    page.text(xs + 3.3, yk, 'tested, not DE', size=sz, gid='hm-key')
    xs += 3.3 + page.text_width('tested, not DE', sz) + 3.0
    page.rect(xs, yk - 1.9, 2.6, 2.2, fill='white', edge=P.grey(0.75), lw=0.25, gid='hm-key')
    page.text(xs + 3.3, yk, 'not tested (too few expressing cells)', size=sz, gid='hm-key')
    return yk + 6.0, n_de


def si_de_heatmap(ctx, fig):
    st = ctx.cfg['style']
    page = Page(fig['width_mm'], fig['height_mm'], st['font'], st['label_pt'])
    draw_de_heatmap(ctx, page, fig, fig['top'], fig['x_left'], fig['width_mm'] - fig['right_margin'],
                    letter=fig.get('letter'))
    return page


def si_dof_strip_de(ctx, fig):
    """Supplementary figure a (maps), b (dot strip), c (DE heatmap on the same cluster columns)."""
    page = si_dof_strip(ctx, fig)
    st = ctx.cfg['style']
    yb = max(it[2][3] for it in page.items) + fig['c_gap']
    draw_de_heatmap(ctx, page, fig, yb, fig['x_left'], fig['width_mm'] - fig['right_margin'], letter='c')
    return page


def _nom_mark(page, m, x, y, d, colours, gid):
    """WT = filled circle, fer-4 = filled diamond (shape + colour, so greyscale and CVD still separate them)."""
    fill = tuple(c / 255 for c in col2rgb(colours[m]))
    if m == 'fer-4':
        r = d * 0.62
        page.polygon([x, x + r, x, x - r], [y - r, y, y + r, y], fill=fill, edge=fill, lw=0.3, gid=gid, zorder=4)
    else:
        page.dots(np.array([x]), np.array([y]), np.array([d]), np.array([fill]), edge=fill, lw=0.3, gid=gid, zorder=4)


def draw_nomination_grid(ctx, page, fig, top, x0, x1, letter='d'):
    """Same columns as the DE heatmap (no header repeated): where each DOF was nominated by the network analysis."""
    st = ctx.cfg['style']; hg = fig['heatmap']; sg = dict(fig['strip']); ng = fig['nominations']
    nom = pd.read_csv(ctx.cfg['nomination_table'])
    _letter(page, 0.0, top + 2.9, letter, st)
    cols, pitch, groups, widths, _ = _column_header(ctx, page, top, sg, x0, x1, reserve_right=hg['count_w'], draw=False)
    rp = hg['row_pitch']; y0 = top + ng['top_pad']
    genes = fig['genes']
    y_end = y0 + len(genes) * rp
    x = x0
    for k, ((ct, rr), w) in enumerate(zip(groups, widths)):
        if k % 2 == 0 and st['rows'].get('band_grey'):
            page.rect(x - sg['group_gap'] / 2, y0 - 0.6, w + sg['group_gap'], y_end - y0 + 1.2,
                      fill=P.grey(st['rows']['band_grey']), gid='nom-band', zorder=0)
        x += w + sg['group_gap']
    look = {(r.symbol, r.cluster): r.more_central_in for r in nom.itertuples()}
    de = pd.read_csv(ctx.cfg['de_table']).set_index(['symbol', 'cluster'])
    cw, ch = pitch - hg['cell_gap'], rp - hg['cell_gap']
    shown = both_de = 0
    for i, g in enumerate(genes):
        yy = y0 + i * rp
        page.text(x0 - 1.2, yy + rp / 2, ctx.gene_label(g), size=hg['gene_pt'], ha='right', va='center', gid='nom-gene')
        for cl, xc, ct in cols:
            page.rect(xc - cw / 2, yy + hg['cell_gap'] / 2, cw, ch, fill='white', edge=P.grey(0.82), lw=0.25,
                      gid='nom-cell', zorder=2)
            m = look.get((g, cl))
            if not isinstance(m, str):
                continue
            shown += 1
            _nom_mark(page, m, xc, yy + rp / 2, ng['mark_d'], ng['colours'], f'nom-{m}')   # Fig. 6a colours
            if (g, cl) in de.index and bool(de.loc[(g, cl), 'DE']):
                both_de += 1
                page.rect(xc - cw / 2, yy + hg['cell_gap'] / 2, cw, ch, fill='none', edge=(0, 0, 0), lw=0.9,
                          gid='nom-also-de', zorder=5)
    yk = y_end + hg['key_dy']; sz = st['small_pt']
    page.text(x0, yk, 'Network nomination (CellOracle out-degree rank change):', size=sz, gid='nom-key')
    xs = x0 + page.text_width('Network nomination (CellOracle out-degree rank change):', sz) + 2.0
    for m, lab in [('WT', 'more central in WT'), ('fer-4', 'more central in *fer-4*')]:
        _nom_mark(page, m, xs + 1.0, yk - 0.75, ng['mark_d'], ng['colours'], 'nom-key')
        page.text(xs + 2.4, yk, lab, size=sz, gid='nom-key')
        xs += 2.4 + page.text_width(lab, sz) + 3.5
    page.rect(xs, yk - 1.9, 2.6, 2.2, fill='none', edge=(0, 0, 0), lw=0.9, gid='nom-key')
    page.text(xs + 3.3, yk, 'also DE in that cluster (a)', size=sz, gid='nom-key')
    n_hidden = int((nom.status == 'cluster not DE-eligible').sum())
    hidden_ctx = nom.loc[nom.status == 'cluster not DE-eligible', 'cluster'].nunique()
    page.text(x0, yk + 3.0, f"{n_hidden} further nominations in {hidden_ctx} clusters with fewer than 5 cells in a sample "
              f"are not shown (Source Data).", size=sz, colour=P.grey(0.25), gid='nom-note')
    return yk + 5.0, shown, both_de


def si_de_and_nominations(ctx, fig):
    st = ctx.cfg['style']
    page = Page(fig['width_mm'], fig['height_mm'], st['font'], st['label_pt'])
    x0, x1 = fig['x_left'], fig['width_mm'] - fig['right_margin']
    yb, _ = draw_de_heatmap(ctx, page, fig, fig['top'], x0, x1, letter='a')
    _, shown, both = draw_nomination_grid(ctx, page, fig, yb + fig['d_gap'], x0, x1, letter='b')
    print(f"  nominations drawn {shown}; also DE {both}")
    return page


LAYOUTS = {'si_dot': si_dot, 'si_violin': si_violin,
           'ed7_dot': lambda c, f: ed7(c, f, 'dot'), 'ed7_violin': lambda c, f: ed7(c, f, 'violin'),
           'si_dot_geno': si_dot_geno, 'ed7_dot_geno': ed7_geno,
           'ed7_d_slot': ed7_d_slot, 'ed7_e_strip': ed7_e_strip, 'si_dof_combined': si_dof_combined, 'si_dof_strip': si_dof_strip, 'si_de_heatmap': si_de_heatmap,
           'si_dof_strip_de': si_dof_strip_de, 'si_de_and_nominations': si_de_and_nominations}
