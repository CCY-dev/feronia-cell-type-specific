"""UMAP maps for Fig. 1a: one shared frame for every map, points rasterized at a fixed resolution inside the vector
page, every label vector text.

Each map is rendered once into an RGB array at `ppi` pixels per inch at its printed size (Agg, deterministic) and
placed with Page.image, so the PDF holds one image per map (no resampling: the array has exactly the printed size in
pixels) while all text, keys and the axis glyph stay vector. All maps share the same data limits (the coordinate
range of all cells plus `pad_frac` on each side) and equal x/y scaling.
"""
from __future__ import annotations

import math

import numpy as np

from .fig1_colour import Gradient2, col2rgb

MINUS = '\u2212'


def frame(x, y, pad_frac):
    """Data limits shared by all maps: full range of all cells plus pad_frac of the range on each side."""
    xr = (x.min() - pad_frac * np.ptp(x), x.max() + pad_frac * np.ptp(x))
    yr = (y.min() - pad_frac * np.ptp(y), y.max() + pad_frac * np.ptp(y))
    return xr, yr


def map_size_px(w_mm, xr, yr, ppi):
    """Pixel size of a map that is w_mm wide at ppi, with equal x/y scaling. Returns (w_px, h_px, w_mm, h_mm) with the
    mm sizes recomputed from whole pixels so the placed image is exactly ppi at print size."""
    w_px = int(round(w_mm / 25.4 * ppi))
    h_px = int(round(w_px * (yr[1] - yr[0]) / (xr[1] - xr[0])))
    return w_px, h_px, w_px * 25.4 / ppi, h_px * 25.4 / ppi


def render(xs, ys, rgb01, xr, yr, w_px, h_px, ppi, point_pt2):
    """Scatter of the given points (drawn in the given order) -> RGB uint8 array (h_px, w_px, 3) on white."""
    import matplotlib
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    fig = Figure(figsize=(w_px / ppi, h_px / ppi), dpi=ppi)
    FigureCanvasAgg(fig)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_axis_off()
    ax.set_xlim(*xr)
    # equal scaling: the y span shown is the x span times the pixel aspect, centred on the data's y range
    yc = (yr[0] + yr[1]) / 2
    ys_span = (xr[1] - xr[0]) * h_px / w_px
    ax.set_ylim(yc - ys_span / 2, yc + ys_span / 2)
    with matplotlib.rc_context({'path.simplify': False}):
        ax.scatter(xs, ys, s=point_pt2, c=rgb01, linewidths=0, marker='o')
    fig.canvas.draw()
    arr = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
    assert arr.shape[:2] == (h_px, w_px), (arr.shape, h_px, w_px)
    return arr


def colours_for(df, kind, spec, palette, fer_scale):
    """Per-cell colours (0-1 RGB) and drawing order for one map kind: 'celltype', 'stage' or 'fer'."""
    if kind == 'fer':
        rgb = fer_scale.rgb01(df[spec['fer_column']].astype(float).values)
        if spec.get('fer_order', 'object') == 'value':                       # FeaturePlot order = TRUE: expressing cells on top
            return rgb, np.argsort(df[spec['fer_column']].astype(float).values, kind='stable')
        return rgb, np.arange(len(df))                                      # FeaturePlot order = FALSE: object order
    col = spec['columns'][kind]
    pal = palette[palette['map'] == kind].set_index('name')
    missing = sorted(set(df[col]) - set(pal.index))
    if missing:
        from ..schemas import ValidationError
        raise ValidationError([f"fig1a_palette.csv: no {kind} colour for {missing}. Re-run tools/fig1_tables.py or add "
                               f"the name to the palette table."])
    rgb = np.array([col2rgb(h) for h in pal.loc[df[col], 'hex']], float) / 255
    if (pal['draw_order'] == 'object').all():
        return rgb, np.arange(len(df))
    rank = pal['draw_order'].astype(int).reindex(df[col]).values               # Seurat DimPlot(order = ...)
    return rgb, np.argsort(rank, kind='stable')


def draw_grid(page, df, palette, spec, geo, ctx_names, lines=None):
    """Six maps: rows = spec['rows'] (cell type, stage, FER), columns = genotypes. geo = layout.a of the config (mm);
    lines = line widths (style lines with the figure's overrides). Returns geometry + record."""
    lines = lines or page.st['lines']
    x, y = df['umap_1'].astype(float).values, df['umap_2'].astype(float).values
    xr, yr = frame(x, y, geo['pad_frac'])
    w_px, h_px, w, h = map_size_px(geo['map_w'], xr, yr, geo['ppi'])
    fer_vals = df[spec['fer_column']].astype(float).values
    fer_scale = Gradient2(spec['fer_colours'][0], spec['fer_colours'][1], fer_vals.min(), fer_vals.max())
    gens = list(spec['genotypes'])
    x0, y0 = geo['x_maps'], geo['y_maps']
    # positions follow the nominal map size (map_w x map_w * aspect); each image is placed top-left at its slot with
    # its whole-pixel size (at most one pixel, 0.04 mm, smaller), so rounding never accumulates down the page
    w_nom = geo['map_w']
    h_nom = w_nom * (yr[1] - yr[0]) / (xr[1] - xr[0])
    col_x = [x0 + i * (w_nom + geo['gap_x']) for i in range(len(gens))]
    row_y = [y0 + j * (h_nom + geo['gap_y']) for j in range(len(spec['rows']))]
    placed, n_cells = {}, {}
    for j, row in enumerate(spec['rows']):
        for i, g in enumerate(gens):
            sub = df[df[spec['columns']['genotype']] == g]
            n_cells[g] = len(sub)
            rgb, order = colours_for(sub, row['map'], spec, palette, fer_scale)
            arr = render(x[sub.index.values][order], y[sub.index.values][order], rgb[order], xr, yr, w_px, h_px,
                         geo['ppi'], geo['point_size_pt2'])
            placed[(row['map'], g)] = page.image(arr, col_x[i], row_y[j], w=w, gid=f'fig1-a-{row["map"]}-{g}',
                                                 name=f'fig1a_{row["map"]}_{g}')
    # genotype headers (bold), sample header with a bracket over both columns
    gh, hd = geo['genotype'], geo['header']
    for i, g in enumerate(gens):
        page.text(col_x[i] + w_nom / 2, gh['baseline'], ctx_names['genotype'](g), size=gh['size_pt'], ha='center',
                  weight='bold', gid='fig1-a-header')
    head = spec['header'].format(n_cells=f'{len(df):,}')
    bx0, bx1 = col_x[0] + hd['bracket_inset'], col_x[-1] + w_nom - hd['bracket_inset']
    page.text((bx0 + bx1) / 2, hd['baseline'], head, size=hd['size_pt'], ha='center', gid='fig1-a-header')
    yb, tk = hd['baseline'] + hd['bracket_gap'], hd['bracket_tick']
    page.polyline([(bx0, yb + tk), (bx0, yb), (bx1, yb), (bx1, yb + tk)], lw=lines['bracket_pt'],
                  colour=(hd['bracket_grey'],) * 3, gid='fig1-a-header')
    # row labels, rotated, centred on each row
    rl = geo['row_label']
    for j, row in enumerate(spec['rows']):
        page.text(rl['x'], row_y[j] + h_nom / 2, row['label'], size=rl['size_pt'], rotation=90, ha='center',
                  gid='fig1-a-rowlabel')
    # FER colour key right of the last FER map: continuous bar (raster at ppi), ticks and label vector
    fk = geo['fer_key']
    j_fer = [r['map'] for r in spec['rows']].index('fer')
    kx = col_x[-1] + w_nom + fk['gap']
    kh = fk['h'] if fk.get('h') else h * fk['h_frac']
    ky1 = row_y[j_fer] + h_nom - fk['bottom_inset']
    kpx_w, kpx_h = int(round(fk['w'] / 25.4 * geo['ppi'])), int(round(kh / 25.4 * geo['ppi']))
    kw, kh = kpx_w * 25.4 / geo['ppi'], kpx_h * 25.4 / geo['ppi']
    ky0 = ky1 - kh
    vmin, vmax = float(fer_vals.min()), float(fer_vals.max())
    v = vmax - (np.arange(kpx_h) + 0.5) / kpx_h * (vmax - vmin)
    grad = (fer_scale.rgb01(v) * 255).round().astype(np.uint8)
    key_arr = np.repeat(grad[:, None, :], kpx_w, axis=1)
    page.image(key_arr, kx, ky0, w=kw, gid='fig1-a-ferkey', name='fig1a_fer_key')
    page.rect(kx, ky0, kw, kh, lw=fk['frame_pt'], edge=(0.4, 0.4, 0.4), gid='fig1-a-ferkey')
    ticks = np.arange(math.ceil(vmin), math.floor(vmax) + 1e-9, fk['tick_step'])
    for t in ticks:
        yy = ky1 - (t - vmin) / (vmax - vmin) * kh
        page.line(kx + kw, yy, kx + kw + fk['tick_len'], yy, lw=fk['frame_pt'], colour=(0.3, 0.3, 0.3),
                  gid='fig1-a-ferkey')
        page.text(kx + kw + fk['label_dx'], yy + fk['label_dy'], f'{t:g}'.replace('-', MINUS), size=fk['size_pt'],
                  gid='fig1-a-ferkey')
    page.text(kx, ky0 - fk['title_gap'], spec['fer_key_label'], size=fk['size_pt'], gid='fig1-a-ferkey')
    key_right = kx + kw + fk['label_dx'] + max(page.text_width(f'{t:g}', fk['size_pt']) for t in ticks)
    # UMAP axis glyph under the lower-left map: two arrows from one corner, 'UMAP 1' centred under the horizontal
    # arrow, 'UMAP 2' rotated beside the vertical arrow (as on the approved page)
    gl = geo['glyph']
    L, hd_ = gl['arrow_mm'], gl['head_mm']
    hw = gl['head_half'] * hd_
    bottom = row_y[-1] + h_nom
    ox, oy = col_x[0] + gl['dx'], bottom + gl['dy'] + L
    lw = lines['bracket_pt']
    page.line(ox, oy, ox + L - hd_ * 0.5, oy, lw=lw, gid='fig1-a-glyph')
    page.line(ox, oy, ox, oy - L + hd_ * 0.5, lw=lw, gid='fig1-a-glyph')
    page.polyline([(ox + L, oy), (ox + L - hd_, oy - hw), (ox + L - hd_, oy + hw)], lw=0,
                  colour=None, fill='k', closed=True, gid='fig1-a-glyph')
    page.polyline([(ox, oy - L), (ox - hw, oy - L + hd_), (ox + hw, oy - L + hd_)], lw=0,
                  colour=None, fill='k', closed=True, gid='fig1-a-glyph')
    page.text(ox + L / 2 + gl['x_label_dx'], oy + gl['x_label_dy'], spec['glyph_labels'][0], size=gl['size_pt'],
              ha='center', gid='fig1-a-glyph')
    page.text(ox + gl['y_label_dx'], oy - L / 2 + gl['y_label_dy'], spec['glyph_labels'][1], size=gl['size_pt'],
              rotation=90, ha='center', gid='fig1-a-glyph')
    a_bottom = oy + gl['bottom_pad']
    rec = dict(type='umap_grid', n_cells=int(len(df)), n_cells_by_genotype=n_cells, map_mm=[round(w, 3), round(h, 3)],
               map_px=[w_px, h_px], ppi=geo['ppi'], xlim=[round(v_, 4) for v_ in xr], ylim=[round(v_, 4) for v_ in yr],
               fer_range=[vmin, vmax], point_size_pt2=geo['point_size_pt2'],
               maps=[f'{k[0]}:{k[1]}' for k in placed])
    return dict(col_x=col_x, row_y=row_y, w=w_nom, h=h_nom, bottom=a_bottom, right=key_right, rec=rec,
                fer_key=(kx, ky0, kw, kh), glyph=(ox, oy, L))
