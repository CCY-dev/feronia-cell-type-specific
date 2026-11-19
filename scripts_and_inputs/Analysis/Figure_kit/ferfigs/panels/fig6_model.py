"""Fig. 6h: root image with a callout onto the working model, port of draw_model_h.py.

Root image: manifest panel h (native pixels, mirrored left-right as published, crop from the manifest);
the callout box comes from the manifest overlays ('box x0 y0 x1 y1' in pixels of the oriented, cropped image).
The provisional inset label sits at the top of the image exactly as in the candidate (white italic, two lines, on a
black box at 65 % opacity); the build stops if it does not fit inside the image (widen the crop or make h taller).
Model: pill-shaped nodes (20 % tint, 1-pt border in the node colour), arrows for 'activates', a T-bar for
'inhibits'; words, order and colours from the config (labels black, journal: no coloured text).
"""
from __future__ import annotations

from ..schemas import ValidationError
from . import images as I
from .fig6_common import MM, arrow_head, cap_mm, tint


def draw(page, mm: MM, fig_dir, row, cfg, box, gid='fig6-h'):
    st = page.st
    fs = cfg.get('text_pt', st['font']['size_pt'])
    x, y, w, h = box
    pl = I.place(page, fig_dir, row, x, y, h=h, gid=f'{gid}-root')
    # provisional inset label (candidate: 6 pt italic white, two lines, black 65 % box inset 1.2 pt)
    lab = cfg['root_label']
    inset = lab['inset_pt'] * 25.4 / 72
    lw_ = max(page.text_width(f"*{t}*", lab['size_pt']) for t in lab['lines'])
    bw = pl.w_mm - 2 * inset
    if lw_ > bw - 2 * lab['min_pad']:
        raise ValidationError([f"h: the root-image label ({lw_:.1f} mm) does not fit the {pl.w_mm:.1f}-mm image; widen the "
                               f"h crop in data/fig6/fig6_images.csv or make h taller (config h.box)."])
    lh = lab['line_pt'] * 25.4 / 72
    bh = lab['box_pad_pt'] * 25.4 / 72 + lh * len(lab['lines'])
    mm.rect(pl.x0 + inset, pl.y0 + inset, bw, bh, face='#000000', alpha=lab['alpha'], gid=f'{gid}-label', z=5)
    for i, t in enumerate(lab['lines']):
        page.text(pl.x0 + pl.w_mm / 2, pl.y0 + inset + lh * (i + 1), t, size=lab['size_pt'], ha='center',
                  fontstyle='italic', colour=(1, 1, 1), gid=f'{gid}-label')
    # callout box on the image and lines to the model panel
    boxes = [o for o in I.parse_overlays(row.get('overlays')) if o[0] == 'box']
    if not boxes:
        raise ValidationError(['fig6_images.csv: panel h needs a callout overlay "box x0 y0 x1 y1"'])
    _, n, _ = boxes[0]
    bx0, by0 = pl.mm(n[0], n[1])
    bx1, by1 = pl.mm(n[2], n[3])
    cc = cfg['callout_colour']
    mm.rect(bx0, by0, bx1 - bx0, by1 - by0, face='none', edge=cc, lw=cfg['callout_pt'], gid=f'{gid}-callout', z=4)
    ins = cfg.get('panel_inset', 0.0)                  # model panel inset from the box top, bottom and right
    mx0, mx1 = pl.x0 + pl.w_mm + cfg['model_gap'], x + w - ins
    my0, my1 = y + ins, y + h - ins
    mm.rbox(mx0, my0, mx1 - mx0, my1 - my0, cfg['panel_radius'], face=cfg['panel_colour'], gid=f'{gid}-panel', z=0)
    ce = cfg.get('callout_end_inset', 0.1)
    mm.line(bx1, by0, mx0 + ce, my0 + ce, lw=cfg['callout_pt'], colour=cc, gid=f'{gid}-callout', z=1)
    mm.line(bx1, by1, mx0 + ce, my1 - ce, lw=cfg['callout_pt'], colour=cc, gid=f'{gid}-callout', z=1)
    # nodes
    nodes, edges = cfg['nodes'], cfg['edges']
    if len(edges) != len(nodes) - 1:
        raise ValidationError(['h: need one edge kind between each pair of consecutive nodes'])
    nh = cfg['node_h']
    nw = mx1 - mx0 - 2 * cfg['node_side_pad']
    xc = (mx0 + mx1) / 2
    pad = cfg['node_end_pad']
    pitch = (my1 - my0 - 2 * pad - nh) / (len(nodes) - 1)
    ys = [my0 + pad + i * pitch for i in range(len(nodes))]
    widest = 0.0
    for nd, yy in zip(nodes, ys):
        col = nd['colour']
        mm.rbox(xc - nw / 2, yy, nw, nh, nh / 2, face=tint(col, cfg['node_tint']), edge=col, lw=cfg['node_edge_pt'],
                gid=f'{gid}-node', z=3)
        tc = tint(col, 0.55, base=(0, 0, 0)) if cfg.get('coloured_text') else None   # candidate: dark shade
        widest = max(widest, page.text(xc, yy + nh / 2 + cap_mm(fs) / 2, nd['label'], size=fs, ha='center',
                                       weight='bold', colour=tc, gid=f'{gid}-nodelabel'))
    if widest > nw - 1.5:
        raise ValidationError([f'h: node label ({widest:.1f} mm) wider than the node ({nw:.1f} mm)'])
    ac = cfg['arrow_colour']
    hl, hw = cfg.get('arrow_head', [1.5, 1.2])          # head length, full width (mm)
    tb = cfg.get('tbar_half', 1.6)
    for k, kind in enumerate(edges):
        ya, yb = ys[k] + nh + cfg['arrow_gap'], ys[k + 1] - cfg['arrow_gap']
        if kind == 'activates':
            mm.line(xc, ya, xc, yb - hl * 0.8, lw=cfg['arrow_pt'], colour=ac, gid=f'{gid}-arrow', z=3)
            arrow_head(mm, xc, yb, 1.5707963267948966, hl, hw, ac, gid=f'{gid}-arrow', z=3)
        elif kind == 'inhibits':
            mm.line(xc, ya, xc, yb, lw=cfg['arrow_pt'], colour=ac, gid=f'{gid}-inhibit', z=3)
            mm.line(xc - tb, yb, xc + tb, yb, lw=cfg['tbar_pt'], colour=ac, gid=f'{gid}-inhibit', z=3)
        else:
            raise ValidationError([f'h: unknown edge kind {kind!r} (activates | inhibits)'])
    return dict(type='model', root_image_mm=[round(pl.w_mm, 2), round(pl.h_mm, 2)], root_ppi=pl.ppi,
                nodes=[nd['label'] for nd in nodes], edges=edges, label=lab['lines'], label_status=lab['status'])
