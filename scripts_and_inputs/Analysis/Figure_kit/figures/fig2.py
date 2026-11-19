"""Figure 2: cell type-specific FER knockout (TSKO) lines.

a schematic | b FER-GFP images; c meristems | d meristem cell number; e root hairs | f root hair length;
g root length under 1 uM RALF23 | h seedlings with the printed treatment/mock ratios.
Every block is placed at the position given in config/fig2.yaml. Every chart number is computed from the
measurement tables (tools/fig2_tables.py, from the Source Data) at build time; letters/asterisks and printed ratios
are drawn exactly as given in their tables; the build compares the h ratios with mean(1 uM)/mean(Mock) of the g
table and records the result. Image pixels are never changed (orientation and crop only). The build stops with a
message when labels collide or content leaves the page.
"""
from ferfigs import style as S
from ferfigs.panels import dotbar
from ferfigs.panels import fig23_tiles as T
from ferfigs.panels import fig2_schematic as FS
from ferfigs.schemas import ValidationError

WHITE, BLACK = (1.0, 1.0, 1.0), (0.0, 0.0, 0.0)


def build(ctx):
    cfg = ctx.cfg
    L = cfg['layout']
    page = ctx.page()
    fig_dir = ctx.path('data/fig2')
    man = ctx.table(cfg['manifest'], 'image_manifest')
    probs = []
    extents = {}

    # ---------------------------------------------------------------- a: schematic
    sch = cfg['schematic']
    geo = FS.draw(page, sch['x'], sch['y'], sch['w'], sch, gid='fig2-a')
    ctx.record('a', dict(type='schematic', **geo))
    extents['a'] = geo['panel_box_mm']

    # ---------------------------------------------------------------- b, c, e, g: image panels
    def image_panel(p):
        im = cfg['images'][p]
        rows = man[man['panel'] == p]
        if rows.empty:
            raise ValidationError([f"{cfg['manifest']}: no tiles for panel {p}"])
        fit = tuple(im.get('fit_row_height') or ())
        placed, geo = T.place_grid(page, fig_dir, rows, im['x'], im['y'], im['um_per_mm'],
                                   im.get('col_gap', L['tile_gap']), im.get('row_gap', L['tile_gap']),
                                   gid=f'fig2-{p}', valign=im.get('valign', 'top'), fit_row_height=fit)
        bar = WHITE if im.get('bar_colour') == 'white' else BLACK
        for _, r in rows.iterrows():
            T.draw_marks(page, placed[(r['row'], r['col'])], r, cfg['marks'], bar_colour=bar, gid=f'fig2-{p}')
        if im.get('frame_lw'):
            T.frames(page, placed, lw=im['frame_lw'], gid=f'fig2-{p}-frame')
        # labels: the first tile row above, the last tile row below (2 x 4 grids)
        lab = cfg['tile_labels'][im['labels']]
        size, pitch = im.get('label_pt'), im.get('label_pitch')
        first, last = rows['row'].min(), rows['row'].max()
        low = max(pl.y0 + pl.h_mm for pl in placed.values())
        for ri, where in ((first, 'above'), (last, 'below')):
            if where == 'below' and ri == first:
                continue
            sel = rows[rows['row'] == ri]
            centres = {r['genotype']: placed[(ri, r['col'])].x0 + placed[(ri, r['col'])].w_mm / 2
                       for _, r in sel.iterrows()}
            missing = [g for g in sel['genotype'] if g not in lab]
            if missing:
                raise ValidationError([f"config/fig2.yaml tile_labels.{im['labels']}: no label for {missing}"])
            lines = {g: lab[g] for g in sel['genotype']}
            probs.extend(T.check_label_spacing(page, centres, lines, L['min_label_gap'], size, where=f'fig2 {p}'))
            if where == 'above':
                top_row = min(placed[(ri, c)].y0 for c in sel['col'])
                T.labels_above(page, centres, lines, top_row - im['label_gap_above'], size, pitch,
                               gid=f'fig2-{p}-labels')
            else:
                bottom_row = max(placed[(ri, c)].y0 + placed[(ri, c)].h_mm for c in sel['col'])
                low = T.labels_below(page, centres, lines, bottom_row + im['label_gap_below'], size, pitch,
                                     gid=f'fig2-{p}-labels')
        tb, box = T.boxes(rows, placed)
        mag = {r['tile']: round(float(r['um_per_px']) / placed[(r['row'], r['col'])].mm_per_px, 3)
               for _, r in rows.iterrows()}
        ctx.record(p, dict(type='images', tiles=list(rows.sort_values(['row', 'col'])['file']),
                           um_per_mm=im['um_per_mm'], fit_row_height=list(fit), um_per_mm_by_tile=mag,
                           scale_bar=[dict(tile=r['tile'], um=float(r['scale_bar_um']),
                                           mm=round(float(r['scale_bar_um']) / im['um_per_mm'], 3))
                                      for _, r in rows.iterrows() if r.get('scale_bar_um')],
                           tile_mm=[[round(pl.w_mm, 2), round(pl.h_mm, 2)] for pl in placed.values()],
                           tile_boxes_mm=tb, panel_box_mm=box))
        extents[p] = [box[0], box[1], box[2], round(low - box[1], 3)]
        return box

    tops = {p: image_panel(p)[1] for p in ('b', 'c', 'e')}

    # ---------------------------------------------------------------- d, f: charts beside c, e
    ch = cfg['charts']
    for cht, img in ch['rows'].items():
        spec = cfg['panels'][cht]
        df = ctx.table(spec['data'], 'measurements')
        st = ctx.table(spec['stats'], 'annotations') if spec.get('stats') else None
        with T.chart_text(page, ch['font_pt']):
            ax, rec = dotbar.draw(page, (ch['x'], tops[img], ch['w'], ch['h']), df, spec, st, gid=f'fig2-{cht}')
            ax.tick_params(labelsize=ch['font_pt'])
        rec['panel_box_mm'] = [round(v, 3) for v in rec['box_mm']]
        ctx.record(cht, rec)
        extents[cht] = [ch['x'], tops[img], ch['w'], round(rec['x_label_bottom_mm'] - tops[img], 3)]

    # ---------------------------------------------------------------- g: root length under RALF23
    spec = cfg['panels']['g']
    x, y, w, hh = spec['frame']
    g_df = ctx.table(spec['data'], 'measurements')
    st = ctx.table(spec['stats'], 'annotations') if spec.get('stats') else None
    with T.chart_text(page, ch['font_pt']):
        ax, rec = dotbar.draw(page, (x, y, w, hh), g_df, spec, st, gid='fig2-g')
        ax.tick_params(labelsize=ch['font_pt'])
    rec['panel_box_mm'] = [round(v, 3) for v in rec['box_mm']]
    k_ = spec['key']
    T.key(page, x + w / 2, y + k_['dy'], [(t['text'], S.series_colour(t['series'])) for t in k_['items']],
          title=k_.get('title'), marker_gap=0.9, item_gap=2.4, title_gap=1.6, size=ch['font_pt'], gid='fig2-g-key')
    ctx.record('g', rec)
    extents['g'] = [x, y, w, round(rec['x_label_bottom_mm'] - y, 3)]

    # ---------------------------------------------------------------- h: seedlings with ratios
    hs = cfg['h_seedlings']
    im = cfg['images']['h']
    size = hs['text_pt']
    rows = man[man['panel'] == 'h']
    if rows.empty:
        raise ValidationError([f"{cfg['manifest']}: no tiles for panel h"])
    placed, _ = T.place_grid(page, fig_dir, rows, im['x'], im['y'], im['um_per_mm'], im['col_gap'], im['row_gap'],
                             gid='fig2-h')
    for _, r in rows.iterrows():
        T.draw_marks(page, placed[(r['row'], r['col'])], r, bar_colour=WHITE, gid='fig2-h')
    ncol = rows['col'].max()
    right_edge = 0.0
    for _, r in rows.iterrows():
        pl = placed[(r['row'], r['col'])]
        if r['genotype'] not in hs['row_labels']:
            raise ValidationError([f"config/fig2.yaml h_seedlings.row_labels: no label for {r['genotype']}"])
        word = hs['row_labels'][r['genotype']]
        if r['col'] == 0:
            page.text(im['label_left'], pl.y0 + pl.h_mm / 2, word, size=size, ha='right', va='center',
                      gid='fig2-h-rowlabels')
        elif r['col'] == ncol:
            x0 = pl.x0 + pl.w_mm + im['label_right_gap']
            wd = page.text(x0, pl.y0 + pl.h_mm / 2, word, size=size, ha='left', va='center', gid='fig2-h-rowlabels')
            right_edge = max(right_edge, x0 + (wd or page.text_width(word, size)))
    yh = min(pl.y0 for pl in placed.values()) - im['header_gap']
    page.text(im['label_left'], yh, hs['title'], size=size, ha='right', gid='fig2-h-header')
    for c in sorted(rows['col'].unique()):
        pl = placed[(0, c)]
        for col in hs['columns']:
            page.text(pl.mm(col['x_px'], 0)[0], yh, col['text'], size=size, ha='center', gid='fig2-h-header')
    printed = ctx.table(hs['ratios'], 'annotations')
    m = g_df.groupby(['group', 'series'])['value'].mean()
    check = []
    for _, r in printed.iterrows():
        ratio = float(m[(r['group'], r['series'])] / m[(r['group'], 'Mock')])
        check.append(dict(group=r['group'], series=r['series'], printed=r['label'], g_table_ratio=round(ratio, 4),
                          g_table_ratio_printed_precision=hs['ratio_format'].format(ratio),
                          equal=hs['ratio_format'].format(ratio) == r['label']))
    src = hs.get('ratio_source', 'table')
    if src not in ('table', 'computed'):
        probs.append(f"fig2 h: ratio_source '{src}' must be 'table' or 'computed'")
    shown = {(c['group'], c['series']): (c['printed'] if src == 'table' else c['g_table_ratio_printed_precision'])
             for c in check}
    ser = hs['columns'][1]['series']
    for _, r in rows.iterrows():
        pl = placed[(r['row'], r['col'])]
        if (r['genotype'], ser) not in shown:
            if printed.empty and src == 'table':    # public copy without the printed ratios: nothing to draw
                continue
            probs.append(f"{hs['ratios']}: no ratio for {r['genotype']} / {ser}")
            continue
        xx, yb = pl.mm(hs['ratio_x_px'], hs['ratio_baseline_px'])
        page.text(xx, yb, shown[(r['genotype'], ser)], size=size, ha='center', colour=WHITE, gid='fig2-h-ratios')
    tb, box = T.boxes(rows, placed)
    how = ('as printed in the published figure' if src == 'table'
           else f"mean(1 uM RALF23) / mean(Mock) of the panel g table, {hs['ratio_format']}")
    drawn = dict(columns=['group', 'series', 'drawn', 'printed', 'g_table_ratio_of_means', 'drawn_from'],
                 rows=[[c['group'], c['series'], shown[(c['group'], c['series'])], c['printed'], c['g_table_ratio'], how]
                       for c in check])
    ctx.record('h', dict(type='images',
                         tiles=[f"{r['tile']}: {r['file']} [{r['crop_px']}]" for _, r in
                                rows.sort_values(['row', 'col']).iterrows()],
                         um_per_mm=im['um_per_mm'],
                         scale_bar=[dict(tile=r['tile'], um=float(r['scale_bar_um']),
                                         mm=round(float(r['scale_bar_um']) / im['um_per_mm'], 3))
                                    for _, r in rows.iterrows() if r.get('scale_bar_um')],
                         tile_boxes_mm=tb, panel_box_mm=box, ratio_source=src, ratio_check=check,
                         ratios_equal_g_table=all(c['equal'] for c in check),
                         source_rows={'Ratios drawn on h': drawn},
                         source_tables={'Ratios printed on h': hs['ratios']}))
    x_left = im['label_left'] - max(page.text_width(v, size) for v in hs['row_labels'].values())
    extents['h'] = [x_left, yh - 2.5, max(right_edge, box[0] + box[2]) - x_left, box[1] + box[3] - (yh - 2.5)]

    # ---------------------------------------------------------------- letters and checks
    for k, (x, y_top) in cfg['letters'].items():
        page.panel_letter(k, x, y_top)
    low = max(e[1] + e[3] for e in extents.values())
    if low > page.H + 0.01:
        probs.append(f'fig2: content reaches {low:.1f} mm but size_mm height is {page.H:g} mm')
    right = max(e[0] + e[2] for e in extents.values())
    if right > page.W + 0.01:
        probs.append(f'fig2: content reaches x = {right:.1f} mm, beyond the {page.W:g}-mm width')
    if probs:
        raise ValidationError(probs)
    ctx.note('extents_mm', {k: [round(v, 3) for v in e] for k, e in sorted(extents.items())})
    ctx.note('content_bottom_mm', round(low, 2))
    return page
