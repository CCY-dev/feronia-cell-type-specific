"""Figure 4: DR5rev::GFP images (a, c) and their quantification (b, d), root length with RALF23 and RALF23 &
PEO-IAA in FER/fer-4 and the five _L3 TSKO lines (e), seedlings of Col-0, afb1 and tir1-1/afb2-3 (f) and their
root length (g).

Rows: a + b, c + d, then e, f, g. Every chart value is read from data/fig4/ (measurement tables written from the
Source Data by tools/fig4_tables.py); labels and the printed ratios are drawn exactly as given in their tables.
Image pixels are never changed.
"""
from ferfigs import style as S
from ferfigs.layout import ASCENDER
from ferfigs.panels import dotbar, imagerow
from ferfigs.panels import images as I
from ferfigs.panels.imagerow import tile_px
from ferfigs.schemas import ValidationError

MINUS = '−'


def _check_tiles(rows, groups, pm, panel):
    rows = rows.sort_values('col')
    want = [(g, t) for g in groups for t in pm]
    have = list(zip(rows['genotype'], rows['treatment']))
    if have != want:
        raise ValidationError([f"fig4_images.csv panel {panel}: tiles (by col) are {have}; the config expects {want}."])
    return rows


def _box(placed):
    x0 = min(p.x0 for p in placed); y0 = min(p.y0 for p in placed)
    x1 = max(p.x0 + p.w_mm for p in placed); y1 = max(p.y0 + p.h_mm for p in placed)
    return [round(x0, 3), round(y0, 3), round(x1 - x0, 3), round(y1 - y0, 3)]


def _ratio_check(ratios, table):
    out = []
    m = table.groupby(['group', 'series'])['value'].mean()
    for _, r in ratios.iterrows():
        g, s, lab = r['group'], r['series'], r['label']
        dec = len(lab.split('.')[1]) if '.' in lab else 0
        v = float(m[(g, s)] / m[(g, 'Mock')])
        out.append(dict(group=g, series=s, printed=lab, ratio_of_means=round(v, 4), matches=f'{v:.{dec}f}' == lab))
    return out


def validate(ctx):
    for p in ('e', 'f'):
        ctx.table(ctx.cfg['panels'][p]['ratios'], 'annotations')
    return 2


def build(ctx):
    cfg = ctx.cfg
    im = cfg['images']
    page = ctx.page()
    st = page.st
    size = st['font']['size_pt']
    fig_dir = ctx.path('data/fig4')
    man = ctx.table(cfg['manifest'], 'image_manifest')
    boxes = {}

    def letter(p, base=None):
        x, b0 = cfg['letters'][p]
        base = b0 if base is None else base
        top = base - ASCENDER * st['font']['panel_label_pt'] / S.PT_PER_MM
        page.panel_letter(p, x, top)
        boxes[f'letter {p}'] = [round(x, 3), round(base, 3), None, None]

    # ---- layout: row 1 a + b, row 2 c + d, row 3 e, f, g. Each image row has its own
    # magnification (its own vector scale bar from its first tile) and its chart beside it; the chart frame top is
    # level with the tile tops and its x axis with the tile bottoms. One Fire LUT key (same LUT for a and c) at the
    # right of row 1.
    pm = im['pm']
    ra = _check_tiles(man[man['panel'] == 'a'], cfg['panels']['a']['groups'], pm, 'a')
    rc = _check_tiles(man[man['panel'] == 'c'], cfg['panels']['c']['groups'], pm, 'c')
    ch = cfg['charts']
    tables, recs = {}, {}

    def chart(p, x, y, w, h):
        s = cfg['panels'][p]
        df = ctx.table(s['data'], 'measurements')
        tables[p] = df
        stt = ctx.table(s['stats'], 'annotations') if s.get('stats') else None
        ax, rec = dotbar.draw(page, (x, y, w, h), df, s, stt, gid=f'fig4-{p}')
        recs[p] = rec
        boxes[p] = [round(x, 3), round(y, 3), round(w, 3), round(h, 3)]
        return rec

    def pairs(rows, groups):
        return [(S.display(g), [2 * i, 2 * i + 1]) for i, g in enumerate(groups)], \
               {int(c): pm[t] for c, t in zip(rows['col'], rows['treatment'])}

    def image_row(p, rows, y_rep, pm_label):
        r = im['rows'][p]
        upp = {float(v) for v in rows['um_per_px']}
        widths = {tile_px(fig_dir, q)[0] for _, q in rows.iterrows()}
        if len(upp) != 1 or len(widths) != 1:
            raise ValidationError([f'fig4_images.csv panel {p}: one um_per_px and one crop width needed '
                                   f'(found um/px {sorted(upp)}, widths {sorted(widths)} px).'])
        mpp = r['tile_w_mm'] / widths.pop()
        g, pmd = pairs(rows, cfg['panels'][p]['groups'])
        off = im['header_offsets']
        out = imagerow.draw(page, fig_dir, rows, r['x0'], y_rep + off['tiles'], mm_per_px=mpp, groups=g,
                            within=im['within_gap'], between=im['between_gap'], pm=pmd, pm_label=pm_label,
                            gid=f'fig4-{p}', y_groups=y_rep + off['groups'], y_pm=y_rep + off['pm'])
        pls = list(out['placed'].values())
        x1 = max(q.x0 + q.w_mm for q in pls)
        page.header(r['x0'], x1, y_rep, im['reporter'], gap=im.get('reporter_gap', 0.9), gid=f'fig4-{p}-reporter')
        boxes[p] = _box(pls)
        ctx.record(p, dict(type='images', tiles=list(rows['file']), crop_px=list(rows['crop_px']),
                           mm_per_px=round(mpp, 6), um_per_mm=round(upp.pop() / mpp, 3),
                           tile_mm=[round(pls[0].w_mm, 3), round(max(q.h_mm for q in pls), 3)], ppi=pls[0].ppi))
        y_top = min(q.y0 for q in pls); y_bot = max(q.y0 + q.h_mm for q in pls)
        return x1, y_top, y_bot

    def letter_at(p, y_rep):
        letter(p, y_rep + im.get('letter_dy', 0.3))

    # row 1: a + b
    y1 = im['rows']['a']['y_reporter']
    a_x1, a_top, a_bot = image_row('a', ra, y1, im['pm_label'])
    bx, bw = ch['boxes_x']['b']
    if bx < a_x1 + 12.0:
        raise ValidationError(['fig4: chart b starts too close to the a images (room for its y label).'])
    chart('b', bx, a_top, bw, a_bot - a_top)
    ctx.record('b', recs['b'])
    letter_at('a', y1); letter_at('b', y1)
    k = im['key']
    kx = bx + bw + k['gap']
    kp = page.image(fig_dir / k['file'], kx, a_bot - k['bottom_offset'] - k['h'], h=k['h'], gid='fig4-key',
                    name=k['file'])
    page.text(kx - 0.2, kp.y0 - 0.7, k['labels'][0], gid='fig4-key')
    page.text(kx - 0.2, kp.y0 + kp.h_mm + 2.46, k['labels'][1], gid='fig4-key')
    if kx + max(page.text_width(t) for t in k['labels']) > page.W - 0.3:
        raise ValidationError(['fig4: the colour key ends beyond the page; narrow chart b.'])
    boxes['key'] = _box([kp])

    # row 2: c alone (full width)
    y2 = max(recs['b']['x_label_bottom_mm'], kp.y0 + kp.h_mm + 2.46) + ch['row2_gap']
    c_x1, c_top, c_bot = image_row('c', rc, y2, im['pm_label'])
    if c_x1 > page.W - 0.5:
        raise ValidationError(['fig4: the c image row ends beyond the page; lower images.rows.c.tile_w_mm.'])
    letter_at('c', y2)

    # row 3: d + e (same frame top and height, so their x axes are level)
    y3 = c_bot + ch['row3_gap']
    for p in ('d', 'e'):
        chart(p, ch['boxes_x'][p][0], y3, ch['boxes_x'][p][1], ch['h'])
        letter(p, y3 + 0.4)
    ctx.record('d', recs['d'])

    # e: treatment/mock ratios printed inside the treatment columns (as given in the table), checked vs the data
    se = cfg['panels']['e']
    rat_e = ctx.table(se['ratios'], 'annotations')
    bx = {(b['group'], b['series']): b for b in recs['e']['bars']}
    ybase = y3 + ch['h'] - se.get('ratio_y_mm', 0.7)
    drawn_e = []
    for _, r in rat_e.iterrows():
        if (r['group'], r['series']) not in bx:
            raise ValidationError([f"{se['ratios']}: no bar for {r['group']} / {r['series']}."])
        b = bx[(r['group'], r['series'])]
        page.text(b['x_mm'], ybase, r['label'], size=se.get('ratio_size_pt', 5.5), rotation=90, ha='left',
                  va='center', gid='fig4-e-ratio')
        drawn_e.append([r['group'], r['series'], r['label']])
    chk_e = _ratio_check(rat_e, tables['e'])
    recs['e']['ratios_drawn'] = drawn_e
    recs['e']['ratio_check'] = chk_e
    recs['e']['ratio_matches'] = f"{sum(c['matches'] for c in chk_e)} of {len(chk_e)}"
    recs['e']['source_rows'] = {'Root length relative to mock, as printed in the treatment columns': dict(
        columns=['group', 'series', 'printed', 'ratio_of_means_from_e'],
        rows=[[c['group'], c['series'], c['printed'], c['ratio_of_means']] for c in chk_e])}
    recs['e']['source_tables'] = {'Root length relative to mock (as printed, transcription)': se['ratios']}
    ctx.record('e', recs['e'])

    # ---- row 4: f seedling composite (genotype labels at the left, RALF23 -/+ over it) + g, whose frame spans the
    # image height (frame top = image top, x axis = image bottom)
    sd = cfg['seedlings']
    y4 = max(recs['d']['x_label_bottom_mm'], recs['e']['x_label_bottom_mm']) + ch['row4_gap']
    rf = man[man['panel'] == 'f']
    if len(rf) != 1:
        raise ValidationError(['fig4_images.csv: panel f needs exactly one row (the seedling composite).'])
    r = rf.iloc[0]
    pl = I.place(page, fig_dir, r, sd['x'], y4, h=sd['h'], gid='fig4-f-image')
    chart('g', ch['boxes_x']['g'][0], pl.y0, ch['boxes_x']['g'][1], pl.h_mm)
    k_ = pl.mm_per_px
    ysg = y4 - sd['pm_gap']
    for t, xp in sd['pm_x_px'].items():
        page.text(pl.x0 + xp * k_, ysg, MINUS if pm[t] in ('-', MINUS) else '+', size=7.0, ha='center',
                  gid='fig4-f-pm')
    page.text(pl.x0 + min(sd['pm_x_px'].values()) * k_ - 3.0, ysg, sd['pm_label'], ha='right', gid='fig4-f-pm')
    for g, (r0, r1) in zip(sd['genotypes'], sd['rows_px']):
        page.text(pl.x0 - sd['label_gap'], pl.y0 + (r0 + r1) / 2 * k_, S.display(g), ha='right', va='center',
                  gid='fig4-f-genotype')
    rat_f = ctx.table(cfg['panels']['f']['ratios'], 'annotations')
    white = S.rgb(sd.get('ratio_colour', [255, 255, 255]))
    drawn_f = []
    for _, rr in rat_f.iterrows():
        if rr['group'] not in sd['genotypes'] or rr['series'] != 'RALF23':
            raise ValidationError([f"{cfg['panels']['f']['ratios']}: {rr['group']} / {rr['series']} has no block."])
        i = sd['genotypes'].index(rr['group'])
        yb = pl.y0 + sd['ratio_baseline_px'][i] * k_
        if not sd['rows_px'][i][0] < sd['ratio_baseline_px'][i] < sd['rows_px'][i][1]:
            raise ValidationError([f"seedlings: ratio baseline of {rr['group']} is outside its block."])
        page.text(pl.x0 + sd['ratio_x_px'] * k_, yb, rr['label'], ha='center', colour=white, gid='fig4-f-ratio')
        drawn_f.append([rr['group'], rr['series'], rr['label']])
    chk_f = _ratio_check(rat_f, tables['g'])
    lab_left = min(pl.x0 - sd['label_gap'] - page.text_width(S.display(g)) for g in sd['genotypes'])
    xf = cfg['letters']['f'][0]
    if lab_left < xf + 1.0:
        raise ValidationError([f'fig4: the f genotype labels start at {lab_left:.1f} mm, too close to letter f.'])
    if pl.x0 + pl.w_mm + 1.0 > cfg['letters']['g'][0]:
        raise ValidationError(['fig4: the f image runs into letter g; lower seedlings.h or move g.'])
    letter('f', y4 + 0.4)
    letter('g', y4 + 0.4)
    ctx.record('g', recs['g'])
    boxes['f'] = _box([pl])
    ctx.record('f', dict(type='seedlings', file=r['file'], ppi=pl.ppi, box_mm=boxes['f'], ratios_drawn=drawn_f,
                         ratio_check_vs_g=chk_f, ratio_matches_vs_g=f"{sum(c['matches'] for c in chk_f)} of {len(chk_f)}",
                         scale_bar='1-cm bar burned into the pixels, shown as drawn (no calibration known)',
                         source_rows={'Root length relative to mock, as printed': dict(
                             columns=['group', 'series', 'printed', 'ratio_of_means_from_g'],
                             rows=[[c['group'], c['series'], c['printed'], c['ratio_of_means']] for c in chk_f])},
                         source_tables={'Root length relative to mock (as printed, transcription)':
                                        cfg['panels']['f']['ratios']}))
    bottom = max(recs['e']['x_label_bottom_mm'], recs['g']['x_label_bottom_mm'], pl.y0 + pl.h_mm)
    if bottom > page.H - 0.5:
        raise ValidationError([f'fig4: content ends at {bottom:.1f} mm but size_mm gives {page.H:g} mm; raise the '
                               f'height in config/fig4.yaml.'])
    ctx.note('panel_boxes_mm', boxes)
    ctx.note('content_bottom_mm', round(bottom, 2))
    return page
