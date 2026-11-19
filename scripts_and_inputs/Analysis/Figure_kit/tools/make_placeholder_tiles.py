#!/usr/bin/env python
"""Write blank stand-in tiles for the photographs and micrographs of Figs. 2-6.

The figures combine charts (drawn from the Source Data) with original images. The images are not part of this
repository, but the figure layout depends on their pixel size: tiles are scaled to a fixed magnification and the
charts are placed relative to them. This script writes, for every image listed in data/figN/figN_image_sizes.csv, a
uniformly grey PNG with exactly that width, height and colour mode. With these tiles the figures build with every
chart, label and scale bar in its published position; only the image areas are blank.

To rebuild the complete figures, put the original image files at the same paths (data/figN/images/...); the
manifests data/figN/figN_images.csv give the crop, orientation and pixel size used for each tile.

usage (from the kit folder):
  python tools/make_placeholder_tiles.py            # write missing tiles only (never overwrites a real image)
  python tools/make_placeholder_tiles.py --force    # overwrite every listed tile with a blank one
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import pandas as pd
from PIL import Image

KIT = Path(__file__).resolve().parents[1]
GREY = 230


def main(argv=None):
    ap = argparse.ArgumentParser(description='Write blank stand-in image tiles for Figs. 2-6')
    ap.add_argument('--force', action='store_true', help='overwrite existing image files')
    a = ap.parse_args(argv)
    wrote = kept = 0
    for sizes in sorted(glob.glob(str(KIT / 'data' / 'fig*' / 'fig*_image_sizes.csv'))):
        fig_dir = Path(sizes).parent
        for r in pd.read_csv(sizes).itertuples(index=False):
            out = fig_dir / r.file
            if out.exists() and not a.force:
                kept += 1
                continue
            out.parent.mkdir(parents=True, exist_ok=True)
            fill = (GREY,) * len(r.mode) if len(r.mode) > 1 else GREY
            if r.mode == 'RGBA':
                fill = (GREY, GREY, GREY, 255)
            Image.new(r.mode, (int(r.width_px), int(r.height_px)), fill).save(out)
            wrote += 1
    print(f'{wrote} blank tiles written, {kept} existing images kept')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
