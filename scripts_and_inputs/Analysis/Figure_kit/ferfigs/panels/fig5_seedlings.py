"""Seedling photographs shown as genotype blocks (Fig. 5f): rows of blocks (4 + 3 in the proposed layout), each block
cut from a native photo composite at its white separator rows (manifest crop_px), every block at one width and one
magnification, with a genotype header, treatment labels over the first row, the numbers printed on the seedlings
(drawn exactly as given in a table) and a vector scale bar.

Pixels are never changed. A block shorter than block_px can be extended downwards (pad: true) with a flat VECTOR
rectangle in the block's own background colour (median of its bottom rows) placed behind the photograph, so every
block has one height; the earlier candidate did the same extension in pixels. Set pad: false to show the blocks at
their native heights.

spec (config 'seedlings' block):
  block_w: 35.2              block width in mm (all blocks share one width, so one magnification)
  block_px: 258              shown height of every block in native px (crop_px may be taller only if pad is false)
  pad: true                  extend shorter blocks to block_px with their background colour (vector, behind)
  pad_rows: 12               rows used for the background colour
  x_left, x_right            span of the widest row (mm); blocks spread across it, shorter rows centred
  y: 103.0                   baseline of the first row's genotype headers (mm)
  first_header_to_image: 8.4 header baseline to image top for the first row (room for the treatment labels)
  header_to_image: 2.4       the same for the other rows
  row_gap: 4.2               image bottom to the next row's header baseline
  treatments: [{lines: [Mock], x_px: 54.8}, ...]   labels over the first row, centred at x_px of the block
  treatment_line_mm: 2.7     line pitch of two-line treatment labels
  treatment_to_image: 1.3    last treatment baseline to image top
  ratio_x_px: {RALF23: 203.8, ...}                  ratio label centre (native px from the block's left edge)
  ratio_baseline_px: {default: 226, fer-4: 251}     ratio baseline (native px from the block's top)
  ratio_colour: [255, 255, 255]
  ratio_label: {block: FER/fer-4, lines: [..], size_pt: 5.5, x_px: 6, line_mm: 2.3, min_clear_mm: 1.0}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image

from .. import style as S
from ..schemas import ValidationError
from . import images as I


def block_pixels(fig_dir, row):
    """Pixels of one block exactly as placed (orientation and crop from the manifest)."""
    from ..layout import oriented
    with Image.open(Path(fig_dir) / row['file']) as im:
        arr = np.asarray(im.convert('RGB'))
    arr = oriented(arr, row.get('orient') or 'identity')
    c = I.parse_crop(row.get('crop_px'))
    if c:
        arr = arr[int(round(c[1])):int(round(c[3])), int(round(c[0])):int(round(c[2]))]
    return arr


def draw(page, fig_dir, rows, ratios, spec, gid='seedlings'):
    """rows: manifest rows of the panel (row = grid row, col = position in the row). ratios: annotations table
    (group = genotype, series = treatment column, label = text printed as given) or None.
    Returns dict(placed={genotype: Placed}, y_top, y_bottom, blocks=[...], ratios_drawn=[...], x0, x1)."""
    st = page.st
    size = st['font']['size_pt']
    white = S.rgb(spec.get('ratio_colour', [255, 255, 255]))
    bw, block_px = float(spec['block_w']), int(spec['block_px'])
    rows = rows.sort_values(['row', 'col'])
    grid = {}
    for _, r in rows.iterrows():
        grid.setdefault(int(r['row']), []).append(r)
    n_max = max(len(v) for v in grid.values())
    span = float(spec['x_right']) - float(spec['x_left'])
    gap = (span - n_max * bw) / (n_max - 1) if n_max > 1 else 0.0
    if gap < st['gaps_mm']['tile'] - 1e-9:
        raise ValidationError([f"seedlings: {n_max} blocks of {bw} mm do not fit in {span:.1f} mm with the "
                               f"{st['gaps_mm']['tile']} mm minimum gap; lower seedlings.block_w in the config."])
    ratio_rows = {}
    if ratios is not None:
        for _, r in ratios.iterrows():
            ratio_rows.setdefault(r['group'], []).append(r)
        known = {r['genotype'] for _, r in rows.iterrows()}
        bad = sorted(set(ratio_rows) - known)
        if bad:
            raise ValidationError([f"{ratios.attrs.get('path', 'ratios')}: genotype(s) {bad} have no seedling block "
                                   f"in the manifest (known: {sorted(known)})."])
    placed, blocks, drawn = {}, [], []
    y_gen = float(spec['y'])
    y_bottom = y_gen
    for ri in sorted(grid):
        blks = grid[ri]
        x_row = float(spec['x_left']) + (n_max - len(blks)) * (bw + gap) / 2
        y_img = y_gen + float(spec['first_header_to_image'] if ri == min(grid) else spec['header_to_image'])
        bottoms = []
        for ci, r in enumerate(blks):
            x = x_row + ci * (bw + gap)
            arr = block_pixels(fig_dir, r)
            hp, wp = arr.shape[:2]
            k = bw / wp                                     # mm per native px of this block
            if hp > block_px:
                raise ValidationError([f"seedlings: block {r['genotype']} is {hp} px tall after crop_px; the shown "
                                       f"height is block_px = {block_px}. Crop it in the manifest."])
            h_mm = (block_px if spec.get('pad', True) else hp) * k
            pad_rgb = None
            if spec.get('pad', True) and hp < block_px:
                pad_rgb = [int(v) for v in np.median(arr[-int(spec.get('pad_rows', 12)):].reshape(-1, 3), axis=0)]
                page.rect(x, y_img, bw, h_mm, lw=0, edge='none', face=pad_rgb, gid=f'{gid}-pad', zorder=-1)
            pl = I.place(page, fig_dir, r, x, y_img, w=bw, gid=f'{gid}-block')
            placed[r['genotype']] = pl
            I.draw_scalebar(page, pl, r, gid=f'{gid}-scalebar')
            page.header(x, x + bw, y_gen, S.display(r['genotype']), gid=f'{gid}-header')
            if ri == min(grid):
                pitch = float(spec.get('treatment_line_mm', 2.7))
                y_tr = y_img - float(spec.get('treatment_to_image', 1.3))
                for t in spec.get('treatments', []):
                    lines = t['lines']
                    for i, line in enumerate(lines):
                        page.text(x + float(t['x_px']) * k, y_tr - (len(lines) - 1 - i) * pitch, line, ha='center',
                                  gid=f'{gid}-treatment')
            # ratios: one baseline per block (default, or a genotype-specific one), centred on the treatment column
            base_px = spec['ratio_baseline_px'].get(r['genotype'], spec['ratio_baseline_px']['default'])
            for rr in ratio_rows.get(r['genotype'], []):
                if rr['series'] not in spec['ratio_x_px']:
                    raise ValidationError([f"seedlings: no ratio_x_px for series '{rr['series']}' in the config."])
                if base_px > (block_px if spec.get('pad', True) else hp):
                    raise ValidationError([f"seedlings: ratio baseline {base_px} px is below block {r['genotype']}."])
                xc = x + float(spec['ratio_x_px'][rr['series']]) * k
                page.text(xc, y_img + base_px * k, rr['label'], ha='center', colour=white, gid=f'{gid}-ratio')
                drawn.append(dict(group=r['genotype'], series=rr['series'], label=rr['label'], x_mm=round(xc, 3),
                                  y_mm=round(y_img + base_px * k, 3)))
            lab = spec.get('ratio_label')
            if lab and r['genotype'] == lab['block']:
                lsize = float(lab.get('size_pt', size))
                yb = y_img + base_px * k
                lines = lab['lines']
                right = 0.0
                for i, line in enumerate(lines):
                    w = page.text(x + float(lab['x_px']) * k, yb - (len(lines) - 1 - i) * float(lab['line_mm']),
                                  line, size=lsize, colour=white, gid=f'{gid}-ratio-label')
                    right = max(right, x + float(lab['x_px']) * k + w)
                firsts = [d for d in drawn if d['group'] == r['genotype']]
                if firsts:
                    left = min(d['x_mm'] - page.text_width(d['label']) / 2 for d in firsts)
                    if right + float(lab.get('min_clear_mm', 1.0)) > left + 1e-6:
                        raise ValidationError([f"seedlings: the ratio label ends {left - right:.2f} mm before the first "
                                               f"ratio (minimum {lab.get('min_clear_mm', 1.0)} mm); shorten its lines or "
                                               f"lower its size in the config."])
            blocks.append(dict(genotype=r['genotype'], x_mm=round(x, 3), y_mm=round(y_img, 3), w_mm=round(bw, 3),
                               h_mm=round(h_mm, 3), px=[wp, hp], ppi=round(wp / (bw / 25.4), 1), padded_rgb=pad_rgb,
                               file=r['file'], crop_px=r.get('crop_px', '')))
            bottoms.append(y_img + h_mm)
        y_bottom = max(bottoms)
        y_gen = y_bottom + float(spec['row_gap'])
    y_top = float(spec['y']) - 0.716 * size / S.PT_PER_MM
    xs = [b['x_mm'] for b in blocks]
    return dict(placed=placed, y_top=y_top, y_bottom=y_bottom, blocks=blocks, ratios_drawn=drawn, gap_mm=round(gap, 3),
                x0=min(xs), x1=max(xs) + bw)
