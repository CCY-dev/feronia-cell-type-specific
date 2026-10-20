"""Panel builders: UMAP maps (as Fig. 1a), cluster row labels, dot plot, split violins, keys."""
from __future__ import annotations

import math

import numpy as np

from .colour import Gradient2, col2rgb
from .data import violin_density
from .page import MINUS, MM_PER_PT


def grey(v):
    return (v, v, v)


def rgb255(c):
    return tuple(x / 255 for x in c)


def fmt_tick(t):
    s = f'{t:g}'
    return s.replace('-', MINUS)


# ---- UMAP --------------------------------------------------------------------------------------------------------
def umap_frame(x, y, pad_frac):
    """Data limits shared by all maps (Fig. 1a rule): range of all cells plus pad_frac on each side."""
    xr = (x.min() - pad_frac * np.ptp(x), x.max() + pad_frac * np.ptp(x))
    yr = (y.min() - pad_frac * np.ptp(y), y.max() + pad_frac * np.ptp(y))
    return xr, yr


def map_size(w_mm, xr, yr, ppi):
    w_px = int(round(w_mm / 25.4 * ppi))
    h_px = int(round(w_px * (yr[1] - yr[0]) / (xr[1] - xr[0])))
    return w_px, h_px, w_px * 25.4 / ppi, h_px * 25.4 / ppi


def render_map(xs, ys, rgb01, xr, yr, w_px, h_px, ppi, point_pt2):
    """Scatter -> RGB uint8 (h_px, w_px, 3) on white, equal x/y scaling (same code path as Fig. 1a)."""
    import matplotlib
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    # +0.25 px guards against float truncation in Agg (int(inches * dpi)); the canvas is exactly w_px x h_px
    fig = Figure(figsize=((w_px + 0.25) / ppi, (h_px + 0.25) / ppi), dpi=ppi)
    FigureCanvasAgg(fig)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_axis_off()
    ax.set_xlim(*xr)
    yc = (yr[0] + yr[1]) / 2
    ys_span = (xr[1] - xr[0]) * h_px / w_px
    ax.set_ylim(yc - ys_span / 2, yc + ys_span / 2)
    with matplotlib.rc_context({'path.simplify': False}):
        ax.scatter(xs, ys, s=point_pt2, c=rgb01, linewidths=0, marker='o')
    fig.canvas.draw()
    arr = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
    assert arr.shape[:2] == (h_px, w_px), (arr.shape, h_px, w_px)
    return arr


class GeneScale:
    """Seurat FeaturePlot colours for one gene: lightgrey -> blue from the minimum to the maximum over all cells."""

    def __init__(self, values, colours):
        self.vmin, self.vmax = float(np.min(values)), float(np.max(values))
        self.g = Gradient2(colours[0], colours[1], self.vmin, self.vmax)

    def rgb01(self, v):
        return self.g.rgb01(v)


def draw_map(page, sub, values, scale, x, y, frame, geo, order, gid):
    """One map of the cells in `sub` (rows of the per-cell table) coloured by `values`."""
    xr, yr = frame
    w_px, h_px, w, h = map_size(geo['map_w'], xr, yr, geo['ppi'])
    xs, ys = sub['umap_1'].to_numpy(float), sub['umap_2'].to_numpy(float)
    rgb = scale.rgb01(values)
    if order == 'value':                       # Seurat FeaturePlot(order = TRUE): expressing cells drawn last
        idx = np.argsort(values, kind='stable')
    else:                                      # 'object': table order (FeaturePlot order = FALSE, as Fig. 1a)
        idx = np.arange(len(values))
    pt2 = geo['point_size_pt2'] * (geo['map_w'] / 28.0) ** 2
    arr = render_map(xs[idx], ys[idx], rgb[idx], xr, yr, w_px, h_px, geo['ppi'], pt2)
    page.image(arr, x, y, w, h, gid=gid)
    return w, h


def colour_bar(page, x, y, w, h, scale_rgb, vmin, vmax, ticks, ppi, orient, size, lw=0.3, tick_len=0.7,
               label_gap=0.9, gid=None):
    """Continuous colour bar (raster at ppi) with vector frame, ticks and labels. orient 'v' (low at bottom, labels
    right) or 'h' (low at left, labels below). Returns the right/bottom extent in mm."""
    if orient == 'v':
        npx_w, npx_h = max(1, int(round(w / 25.4 * ppi))), int(round(h / 25.4 * ppi))
        v = vmax - (np.arange(npx_h) + 0.5) / npx_h * (vmax - vmin)
        grad = (scale_rgb(v) * 255).round().astype(np.uint8)
        arr = np.repeat(grad[:, None, :], npx_w, axis=1)
    else:
        npx_w, npx_h = int(round(w / 25.4 * ppi)), max(1, int(round(h / 25.4 * ppi)))
        v = vmin + (np.arange(npx_w) + 0.5) / npx_w * (vmax - vmin)
        grad = (scale_rgb(v) * 255).round().astype(np.uint8)
        arr = np.repeat(grad[None, :, :], npx_h, axis=0)
    page.image(arr, x, y, w, h, gid=gid)
    page.rect(x, y, w, h, edge=grey(0.4), lw=lw, gid=gid, zorder=3)
    ext = 0.0
    for t in ticks:
        f = (t - vmin) / (vmax - vmin)
        if orient == 'v':
            yy = y + h - f * h
            page.line(x + w, yy, x + w + tick_len, yy, lw=lw, colour=grey(0.3), gid=gid)
            b = page.text(x + w + tick_len + 0.4, yy, fmt_tick(t), size=size, va='center', gid=gid)
            ext = max(ext, b[2])
        else:
            xx = x + f * w
            page.line(xx, y + h, xx, y + h + tick_len, lw=lw, colour=grey(0.3), gid=gid)
            b = page.text(xx, y + h + tick_len + label_gap + 0.716 * size * MM_PER_PT, fmt_tick(t), size=size,
                          ha='center', gid=gid)
            ext = max(ext, b[3])
    return ext


def nice_ticks(vmin, vmax, step=None):
    span = vmax - vmin
    if step is None:
        step = 1.0 if span >= 2 else 0.5
    t0 = math.ceil(vmin / step - 1e-9) * step
    return [round(t, 6) for t in np.arange(t0, vmax + 1e-9, step)]


def umap_glyph(page, ox, oy, L=4.5, head=0.8, head_half=0.45, size=5.5, lw=0.5, labels=('UMAP 1', 'UMAP 2')):
    """Two arrows from one corner (ox, oy) as in Fig. 1a: 'UMAP 1' under the horizontal arrow, 'UMAP 2' rotated."""
    hw = head_half * head
    page.line(ox, oy, ox + L - head * 0.5, oy, lw=lw, gid='umap-glyph')
    page.line(ox, oy, ox, oy - L + head * 0.5, lw=lw, gid='umap-glyph')
    page.polygon([ox + L, ox + L - head, ox + L - head], [oy, oy - hw, oy + hw], fill='k', gid='umap-glyph', zorder=3)
    page.polygon([ox, ox - hw, ox + hw], [oy - L, oy - L + head, oy - L + head], fill='k', gid='umap-glyph', zorder=3)
    page.text(ox + L / 2 + 1.2, oy + 2.3, labels[0], size=size, ha='center', gid='umap-glyph')
    page.text(ox - 0.9, oy - L / 2 - 1.2, labels[1], size=size, ha='center', rotation=90, gid='umap-glyph')
    return oy + 2.2 + 0.5


# ---- cluster rows ------------------------------------------------------------------------------------------------
def row_positions(rows, y0, pitch, group_gap, min_h=None):
    """Centre y of each row; a gap of group_gap between cell types. A cell type whose label needs more height than
    its rows (min_h[celltype], e.g. a two-line label on a single row) gets equal padding above and below its rows.
    Returns ys, spans {celltype: (y_top, y_bot)}, y_end."""
    min_h = min_h or {}
    groups = []
    for r in rows:
        if not groups or groups[-1][0] != r['celltype']:
            groups.append((r['celltype'], []))
        groups[-1][1].append(r)
    ys, spans, y = [], {}, y0
    for gi, (ct, rr) in enumerate(groups):
        if gi:
            y += group_gap
        h = max(len(rr) * pitch, min_h.get(ct, 0.0))
        pad = (h - len(rr) * pitch) / 2
        for i in range(len(rr)):
            ys.append(y + pad + (i + 0.5) * pitch)
        spans[ct] = (y, y + h)
        y += h
    return np.array(ys), spans, y


def label_min_heights(names, label_pt, pad=0.5):
    """Height each cell-type label needs (lines x leading + pad), from the configured line breaks."""
    lead = label_pt * MM_PER_PT * 1.15
    return {ct: len(lines) * lead + pad for ct, lines in names['celltype_lines'].items()}


def draw_row_block(page, rows, ys, spans, x0, x_right, geo, names, counts, gens, gname, header_y):
    """Bands, cell-type labels (wrapped as configured), stage labels and cells-per-genotype columns.
    Returns x where the chart area starts."""
    sz, szs = geo['label_pt'], geo['small_pt']
    xs_stage = x0 + geo['celltype_w']
    xn = xs_stage + geo['stage_w']
    band = geo.get('band_grey')
    for k, (ct, (t, b)) in enumerate(spans.items()):
        if band and k % 2 == 0:
            page.rect(x0, t, x_right - x0, b - t, fill=grey(band), gid='rows-band', zorder=0)
        lines = names['celltype_lines'].get(ct, [names['celltypes'].get(ct, ct)])
        lead = sz * MM_PER_PT * 1.15
        yc = (t + b) / 2 + 0.716 * sz * MM_PER_PT / 2 - (len(lines) - 1) * lead / 2
        xl = x0 + geo['pad']
        sw = geo.get('swatch')
        if sw:                                   # cell-type colour of the reference map (Fig. 1a palette)
            s = 1.3
            page.rect(xl, (t + b) / 2 - s / 2, s, s, fill=sw.get(ct, grey(0.85)), edge=grey(0.5), lw=0.2,
                      gid='rows-swatch', zorder=2)
            xl += s + 0.8
        for i, ln in enumerate(lines):
            page.text(xl, yc + i * lead, ln, size=sz, gid='rows-celltype')
    for r, yc in zip(rows, ys):
        page.text(xs_stage, yc, names['stages'].get(r['stage'], r['stage']), size=szs, va='center',
                  gid='rows-stage')
    # cells per genotype (right-aligned numbers)
    if geo.get('show_counts', True):
        cw = geo['count_w']
        for j, g in enumerate(gens):
            xr = xn + (j + 1) * cw - 0.6
            page.text(xr, header_y, gname[g], size=szs, ha='right', gid='rows-count-header')
            for r, yc in zip(rows, ys):
                n = counts[(r['cluster'], g)]
                page.text(xr, yc, f'{n:,}', size=szs, ha='right', va='center', colour=grey(0.25),
                          gid='rows-count')
        page.text(xn + len(gens) * cw / 2, header_y - geo['count_title_dy'], 'Cells', size=szs, ha='center',
                  gid='rows-count-header')
        xn += len(gens) * cw
    return xn + geo['chart_gap']


# ---- dot plot ----------------------------------------------------------------------------------------------------
def dot_diameter(pct, dmax):
    """Area proportional to the percentage of expressing cells; 100 % = dmax mm."""
    return dmax * np.sqrt(np.clip(np.asarray(pct, float), 0, 100) / 100.0)


def draw_dots(page, rows, ys, summ, genes, x0, geo, scale, gens, gname, gene_label, header_y):
    dp, gap, dmax = geo['dot_pitch'], geo['gene_gap'], geo['dot_max_d']
    s = summ.set_index(['gene', 'cluster', 'genotype'])
    xs_all = {}
    for i, sym in enumerate(genes):
        xg = x0 + i * (len(gens) * dp + gap)
        xc = xg + len(gens) * dp / 2
        page.text(xc, header_y - geo['gene_dy'], gene_label(sym), size=geo['gene_pt'], ha='center', gid='dot-gene')
        for j, g in enumerate(gens):
            x = xg + (j + 0.5) * dp
            xs_all[(sym, g)] = x
            page.text(x, header_y, gname[g], size=geo['small_pt'], ha='left', va='center', rotation=90,
                      gid='dot-genotype')
            pct, mu, n = [], [], []
            for r in rows:
                rec = s.loc[(sym, r['cluster'], g)]
                n.append(rec['n_cells'])
                pct.append(rec['pct_expressing'])
                mu.append(rec[geo.get('mean_col', 'mean_sct')])
            pct, mu, n = np.array(pct, float), np.array(mu, float), np.array(n)
            ok = (n > 0) & (pct > 0)
            if ok.any():
                page.dots(np.full(ok.sum(), x), ys[ok], dot_diameter(pct[ok], dmax), scale.rgb01(mu[ok]),
                          edge=grey(geo['dot_edge_grey']), lw=geo['dot_edge_pt'], gid=f'dot-{sym}-{g}')
            for yy in ys[n == 0]:
                page.text(x, yy, '\u2013', size=geo['small_pt'], ha='center', va='center', colour=grey(0.5),
                          gid='dot-nd')
    x_end = x0 + len(genes) * (len(gens) * dp + gap) - gap
    return x_end, xs_all


def size_key(page, x, y, pcts, dmax, size, title, gid='dot-size-key', pitch=None):
    """Vertical size key: title then one dot per value with its label. Returns bottom y."""
    page.text(x, y, title, size=size, gid=gid)
    yy = y + 1.2
    pitch = pitch or dmax + 0.9
    for p in pcts:
        yy += pitch
        d = float(dot_diameter(p, dmax))
        page.dots([x + dmax / 2], [yy - pitch / 2], [d], [grey(0.6)], edge=grey(0.45), lw=0.2, gid=gid)
        page.text(x + dmax + 1.0, yy - pitch / 2, f'{p:g}%', size=size, va='center', gid=gid)
    return yy


# ---- violins -----------------------------------------------------------------------------------------------------
def draw_violins(page, rows, ys, df, genes, cfg_genes, x0, geo, gens, gene_label, header_y, styles):
    """Split violins, one column per gene: WT upper half, fer-4 lower half, value axis = SCT expression shared by all
    clusters of the gene (0 to its maximum over all cells). Width-scaled per half (ggplot2 scale = 'width')."""
    cw, gap, halfh = geo['col_w'], geo['col_gap'], geo['half_h']
    y_axis = ys[-1] + geo['row_pitch'] / 2 + geo['axis_gap']
    by_cluster = {c: d for c, d in df.groupby('time.celltype.anno.Li')}
    cols = {}
    for i, sym in enumerate(genes):
        col = cfg_genes[sym]['column']
        xl = x0 + i * (cw + gap)
        vmax_data = float(df[col].max())
        step = geo.get('tick_step') or (1.0 if vmax_data >= 2 else 0.5)
        xmax = math.ceil(vmax_data / step - 1e-9) * step
        sx = lambda v: xl + v / xmax * cw  # noqa: E731
        cols[sym] = (xl, xmax)
        page.text(xl + cw / 2, header_y, gene_label(sym), size=geo['gene_pt'], ha='center', gid='vln-gene')
        for r, yc in zip(rows, ys):
            sub = by_cluster[r['cluster']]
            for j, g in enumerate(gens):
                v = sub.loc[sub['genotype'] == g, col].to_numpy(float)
                kd = violin_density(v)
                sign = -1 if j == 0 else 1             # WT above the centre line, fer-4 below
                st = styles[g]
                if kd is None:
                    continue
                grid, dens = kd
                if len(grid) == 1:                     # all cells equal (almost always all zero): flat line
                    page.line(sx(grid[0]), yc, sx(grid[0]), yc + sign * halfh, lw=geo['flat_pt'],
                              colour=st['edge'], gid=f'vln-{sym}-{g}', zorder=4)
                    continue
                wv = dens / dens.max() * halfh
                xs = np.concatenate([sx(grid), sx(grid[::-1])])
                yv = np.concatenate([yc + sign * wv, np.full(len(grid), yc)])
                page.polygon(xs, yv, fill=st['fill'], edge=st['edge'], lw=geo['edge_pt'], gid=f'vln-{sym}-{g}',
                             zorder=3 + j)
        # value axis
        lw = geo['axis_pt']
        page.line(xl, y_axis, xl + cw, y_axis, lw=lw, gid='vln-axis')
        for t in nice_ticks(0, xmax, step):
            page.line(sx(t), y_axis, sx(t), y_axis + geo['tick_len'], lw=lw, gid='vln-axis')
            page.text(sx(t), y_axis + geo['tick_len'] + 0.6 + 0.716 * geo['small_pt'] * MM_PER_PT, fmt_tick(t),
                      size=geo['small_pt'], ha='center', gid='vln-axis')
        if geo.get('zero_line_grey') is not None:
            page.line(xl, ys[0] - geo['row_pitch'] / 2, xl, y_axis, lw=0.25, colour=grey(geo['zero_line_grey']),
                      gid='vln-zero', zorder=1)
    x_end = x0 + len(genes) * (cw + gap) - gap
    y_lab = y_axis + geo['tick_len'] + 0.6 + 0.716 * geo['small_pt'] * MM_PER_PT + geo['axis_title_dy']
    return x_end, y_lab, cols


def split_key(page, x, y, styles, gens, gname, size, w=6.0, halfh=1.5, edge_pt=0.3):
    """Genotype key for split violins: a sample split shape with WT (upper half) and fer-4 (lower half)."""
    t = np.linspace(0, 1, 60)
    prof = np.sin(np.pi * t) ** 0.8
    for j, g in enumerate(gens):
        sign = -1 if j == 0 else 1
        xs = np.concatenate([x + t * w, x + t[::-1] * w])
        yv = np.concatenate([y + sign * prof * halfh, np.full(len(t), y)])
        page.polygon(xs, yv, fill=styles[g]['fill'], edge=styles[g]['edge'], lw=edge_pt, gid='vln-key', zorder=3 + j)
        page.text(x + w + 1.0, y + sign * 0.95, gname[g], size=size, va='center', gid='vln-key')
    return x + w + 1.0 + max(page.text_width(gname[g], size) for g in gens)


# ---- dot plot grouped by genotype ----------------------------------------------------------------------------------
def draw_genotype_blocks(page, rows, ys, spans, x0, summ, genes, geo, scale, gens, gname, gene_label, counts,
                         header_y, band_x0, show_gene_labels=True):
    """One block per genotype: [cells | gene 1 ... gene n], blocks separated by geo['block_gap'] with the row bands
    broken in the gap. Header per block: genotype name (bold) over a rule, then rotated column labels.
    Returns x_end and the x-centre of every (gene, genotype) column."""
    dp, dmax, cw, bg = geo['dot_pitch'], geo['dot_max_d'], geo['count_w'], geo['block_gap']
    szs = geo['small_pt']
    s = summ.set_index(['gene', 'cluster', 'genotype'])
    show_n = geo.get('show_counts', True)
    if not show_n:
        cw = 0.0
    block_w = cw + (geo['count_gap'] if show_n else 0.0) + len(genes) * dp
    band = geo.get('band_grey')
    xs_all, x = {}, x0
    y_rule = header_y - geo['col_label_h']
    for j, g in enumerate(gens):
        xb0, xb1 = x, x + block_w
        # row bands: the first block continues the label bands; later blocks start fresh after a white gap
        if band:
            for k, (ct, (t, b)) in enumerate(spans.items()):
                if k % 2 == 0:
                    bx0 = band_x0 if j == 0 else xb0 - 0.6
                    page.rect(bx0, t, xb1 + 0.6 - bx0, b - t, fill=grey(band), gid='rows-band', zorder=0)
        # header: genotype name over a rule spanning the block
        page.text((xb0 + xb1) / 2, y_rule - geo['geno_dy'], gname[g], size=geo['gene_pt'], ha='center',
                  weight='bold', gid='dot-genotype-header')
        page.line(xb0, y_rule, xb1, y_rule, lw=geo['rule_pt'], colour=grey(0.25), gid='dot-genotype-rule')
        # cell counts
        xr = xb0 + cw - 0.4
        if show_n:
          page.text(xr - 1.2, header_y, 'Cells', size=szs, ha='left', va='center', rotation=90, colour=grey(0.25),
                  gid='rows-count-header')
          for r, yc in zip(rows, ys):
            page.text(xr, yc, f"{counts[(r['cluster'], g)]:,}", size=szs, ha='right', va='center',
                      colour=grey(0.25), gid='rows-count')
        # dots
        for i, sym in enumerate(genes):
            xc = xb0 + cw + (geo['count_gap'] if show_n else 0.0) + (i + 0.5) * dp
            xs_all[(sym, g)] = xc
            if show_gene_labels:
                page.text(xc, header_y, gene_label(sym), size=geo['gene_label_pt'], ha='left', va='center',
                          rotation=90, gid='dot-gene')
            pct, mu, n = [], [], []
            for r in rows:
                rec = s.loc[(sym, r['cluster'], g)]
                n.append(rec['n_cells']); pct.append(rec['pct_expressing'])
                mu.append(rec[geo.get('mean_col', 'mean_sct')])
            pct, mu, n = np.array(pct, float), np.array(mu, float), np.array(n)
            ok = (n > 0) & (pct > 0)
            if ok.any():
                page.dots(np.full(ok.sum(), xc), ys[ok], dot_diameter(pct[ok], dmax), scale.rgb01(mu[ok]),
                          edge=grey(geo['dot_edge_grey']), lw=geo['dot_edge_pt'], gid=f'dot-{sym}-{g}')
            for yy in ys[n == 0]:
                page.text(xc, yy, '\u2013', size=szs, ha='center', va='center', colour=grey(0.5), gid='dot-nd')
        x = xb1 + bg
    return x - bg, xs_all
