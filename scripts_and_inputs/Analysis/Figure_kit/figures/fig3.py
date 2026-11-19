"""Figure 3: pYUC5:GFP-GUS reporter images (a, c, e), cross sections (f) and their quantification (b, d, g).

Panel slots come from config/fig3.yaml (image_rows a, c, e; section_row f; chart_panels b, d, g).
Reads config/fig3.yaml, data/fig3/fig3_images.csv (tiles) and the b, d, g measurement and stats tables. Every
number on the charts is computed from the measurement tables at build time; the letters are drawn exactly as
given in the stats tables. Pixels of the micrographs are never changed (orientation and crop only).

Layout: full page (a | b, c, d, e | f | g; 180 x 231.6 mm). Every position and size comes from the config
(`letters`, `images`, `charts`); the module only does layout arithmetic. Header words come from the manifest (time
for a, genotype for b-d, treatment for the -/+ signs) and from the config (reporter, treatment label), so a new
tile set is relabelled by its manifest.
"""
import copy

from ferfigs import style as S
from ferfigs.layout import ASCENDER
from ferfigs.panels import dotbar, imagerow


def _runs(rows, column):
    """[(markup, [cols])]: consecutive tiles (by col) that share `column` (e.g. time or genotype) form one group."""
    out = []
    for _, r in rows.sort_values('col').iterrows():
        lab = S.display(str(r[column]))
        if out and out[-1][0] == lab:
            out[-1][1].append(r['col'])
        else:
            out.append((lab, [r['col']]))
    return out


def _key(page, fig_dir, k, x, y_tiles_bottom, gid):
    """Colour key as approved: the LUT strip in a thin frame, 'Max' centred above it and 'Min' centred below it with
    its baseline level with the tile bottoms. Returns the right edge (mm)."""
    y1 = y_tiles_bottom - k['bottom_offset']
    y0 = y1 - k['h']
    kp = page.image(fig_dir / k['file'], x, y0, h=k['h'], gid=gid, name=k['file'])
    page.rect(x, y0, kp.w_mm, k['h'], lw=0.25, edge='k', gid=gid)
    hi, lo = k.get('labels', ('Max', 'Min'))
    cx = x + kp.w_mm / 2
    page.text(cx, y0 - k['max_gap'], hi, ha='center', gid=gid)
    page.text(cx, y_tiles_bottom, lo, ha='center', gid=gid)
    return x + kp.w_mm


def build(ctx):
    cfg = ctx.cfg
    im = cfg['images']
    page = ctx.page()
    # Fig. 3 text size (config font_pt; the approved page used 7-7.25 pt). A private copy of the style, so the
    # shared house style seen by the other figures is unchanged.
    page.st = copy.deepcopy(page.st)
    page.st['font']['size_pt'] = float(cfg.get('font_pt', page.st['font']['size_pt']))
    fs = page.st['font']['size_pt']
    fig_dir = ctx.path('data/fig3')
    man = ctx.table(cfg['manifest'], 'image_manifest')
    signs = im['signs']
    mags = {}

    # ---- image rows (a, c, e): tile rows with reporter / group / -/+ headers and a colour key
    for p in cfg.get('image_rows', ['a', 'c', 'e']):
        rc, rows = im[p], man[man['panel'] == p].sort_values('col')
        first = rows.iloc[0]
        mm_per_px = rc['tile_h'] / imagerow.tile_px(fig_dir, first)[1]
        out = imagerow.draw(page, fig_dir, rows, rc['x0'], rc['y'], mm_per_px=mm_per_px,
                            groups=_runs(rows, rc['group_by']), within=im['within_gap'], between=rc['between'],
                            reporter=im['reporter'], pm={r['col']: signs[r['treatment']] for _, r in rows.iterrows()},
                            pm_label=im['treatment_label'], y_reporter=rc['y_reporter'], y_groups=rc['y_groups'],
                            y_pm=rc['y_pm'], gid=f'fig3-{p}')
        _key(page, fig_dir, im['key'], out['x1'] + im['key']['gap'], out['y_bottom'], gid=f'fig3-{p}-key')
        upp = float(first['um_per_px'])
        mags[p] = dict(mm_per_px=round(mm_per_px, 6), um_per_mm=round(upp / mm_per_px, 3),
                       tile_mm=[round(imagerow.tile_px(fig_dir, first)[0] * mm_per_px, 3), round(rc['tile_h'], 3)])
        ctx.record(p, dict(type='images', tiles=list(rows['file']), mm_per_px=round(mm_per_px, 6)))

    # ---- f: optical cross sections, one display width (um/px confirmed for Col-0 only; see the manifest notes)
    sp = cfg.get('section_row', 'f')
    rd, rows = im[sp], man[man['panel'] == sp].sort_values('col')
    imagerow.draw(page, fig_dir, rows, rd['x0'], rd['y'], tile_w=rd['tile_w'], within=rd['gap'], between=rd['gap'],
                  labels_above={r['col']: S.display(r['genotype']) for _, r in rows.iterrows()},
                  y_groups=rd['y_labels'], gid=f'fig3-{sp}')
    ctx.record(sp, dict(type='images', tiles=list(rows['file'])))

    # ---- charts (b, d, g) on one row; frames (top tick, y axis, width, height) from the config
    ch = cfg['charts']
    for p in cfg.get('chart_panels', ['b', 'd', 'g']):
        s = cfg['panels'][p]
        df = ctx.table(s['data'], 'measurements')
        st = ctx.table(s['stats'], 'annotations') if s.get('stats') else None
        ax, rec = dotbar.draw(page, (ch['x'][p], ch['y'], ch['w'][p], ch['h']), df, s, st, gid=f'fig3-{p}')
        ax.tick_params(axis='y', labelsize=fs)          # tick labels at the figure's text size
        ctx.record(p, rec)

    # ---- panel letters: 8 pt bold at the approved x and baseline
    lpt = page.st['font']['panel_label_pt']
    for letter, (x, base) in cfg['letters'].items():
        page.panel_letter(letter, x, base - ASCENDER * lpt / S.PT_PER_MM)
    ctx.note('magnification_ace', mags)
    return page
