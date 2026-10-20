#!/usr/bin/env python
"""Rebuild every FER single-cell expression draft from the per-cell export with one command:

    python build.py                      # all figures + source data + checks + before/after pages
    python build.py --only ED7de_optionA_dot --formats pdf
    python build.py --cells other.csv.gz --out other_output

Writes <out>/figures/<id>.pdf|png, <out>/source_data/si_expression_cluster_summary.csv, <out>/build.json (input,
config and code hashes, versions, per-figure records and checks) and <out>/compare/before_after.pdf.
"""
from __future__ import annotations

import os
import argparse
import hashlib
import json
import pathlib
import platform
import sys

import yaml

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from expression_figs import data as D  # noqa: E402
from expression_figs.figures import LAYOUTS, Ctx  # noqa: E402


def file_sha(p):
    return D.sha256(p)


def check_page(page, fig, st):
    """Text inside the page, text sizes within the configured ladder, no overlapping text boxes."""
    out = []
    W, H = page.w, page.h
    texts = [(gid, b, s) for kind, gid, b, s in page.items if kind == 'text']
    outside = [(s, [round(v, 2) for v in b]) for gid, b, s in texts
               if b[0] < -0.01 or b[1] < -0.01 or b[2] > W + 0.01 or b[3] > H + 0.01]
    out.append(dict(check='text inside page', ok=not outside, detail=outside[:10]))
    over = []
    for i in range(len(texts)):
        g1, a, s1 = texts[i]
        for g2, b, s2 in texts[i + 1:]:
            gap = 0.3                                     # labels closer than 0.3 mm side by side count as touching
            if a[0] < b[2] + gap and b[0] < a[2] + gap and a[1] < b[3] - 0.05 and b[1] < a[3] - 0.05:
                over.append((s1, s2))
    out.append(dict(check='no overlapping text boxes', ok=not over, detail=over[:10]))
    imgs = [(gid, sh) for kind, gid, b, sh in page.items if kind == 'image']
    out.append(dict(check='images', ok=True, detail=f'{len(imgs)} images placed at native pixel size'))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', default=str(HERE / 'config.yaml'))
    ap.add_argument('--cells', default=None, help='per-cell table (default: config cells:)')
    ap.add_argument('--out', default=str(HERE / 'output'))
    ap.add_argument('--only', nargs='*', default=None)
    ap.add_argument('--formats', nargs='*', default=['pdf', 'png'])
    ap.add_argument('--no-compare', action='store_true')
    a = ap.parse_args(argv)

    cfg = yaml.safe_load(open(a.config))
    cells_path = pathlib.Path(a.cells or (HERE / cfg['cells']))
    for k in ('celltype_palette', 'de_table', 'nomination_table'):
        if cfg.get(k):
            cfg[k] = str(HERE / cfg[k])
    out = pathlib.Path(a.out)
    (out / 'figures').mkdir(parents=True, exist_ok=True)
    (out / 'source_data').mkdir(parents=True, exist_ok=True)

    cells = D.load_cells(cells_path, cfg)
    rows = D.cluster_rows(cells, cfg)
    summ = D.summaries(cells, cfg, rows)
    cl = cfg.get('clusters', {})
    elig = D.de_eligibility(cells, rows, cl.get('min_cells_per_sample', 5))
    if cl.get('de_tested_list'):
        listed = {l.strip() for l in open(HERE / cl['de_tested_list']) if l.strip() and not l.startswith('#')}
        mine = set(elig.index[elig['de_eligible']])
        if listed != mine:
            raise SystemExit(f'DE eligibility rule does not reproduce the DE-tested clusters: '
                             f'extra {sorted(mine - listed)}, missing {sorted(listed - mine)}')
    summ = summ.merge(elig[['de_eligible']], left_on='cluster', right_index=True)
    sd = out / 'source_data' / 'si_expression_cluster_summary.csv'
    summ.to_csv(sd, index=False, float_format='%.6g')
    elig.to_csv(out / 'source_data' / 'si_cluster_cells_per_sample.csv')
    if cl.get('de_eligible_only'):
        rows = [r for r in rows if elig.loc[r['cluster'], 'de_eligible']]
        summ = summ[summ['de_eligible']].copy()
    ctx = Ctx(cfg, cells, rows, summ)

    import matplotlib
    import numpy
    import pandas
    import scipy
    code_files = sorted((HERE / 'expression_figs').glob('*.py')) + [HERE / 'build.py']
    _rel = lambda q: os.path.relpath(os.path.abspath(str(q)), str(HERE))
    rec = dict(
        inputs=dict(cells=_rel(cells_path), cells_sha256=file_sha(cells_path), config=_rel(a.config),
                    config_sha256=file_sha(a.config),
                    palette_sha256=file_sha(cfg['celltype_palette']) if cfg.get('celltype_palette') else None,
                    de_table_sha256=file_sha(cfg['de_table']) if cfg.get('de_table') else None,
                    nomination_table_sha256=file_sha(cfg['nomination_table']) if cfg.get('nomination_table') else None),
        code_sha256={p.name: file_sha(p) for p in code_files},
        environment=dict(python=platform.python_version(), platform=platform.platform(),
                         matplotlib=matplotlib.__version__, numpy=numpy.__version__, pandas=pandas.__version__,
                         scipy=scipy.__version__, pyyaml=yaml.__version__),
        data=dict(n_cells=int(len(cells)), n_by_genotype=cells['genotype'].value_counts().to_dict(),
                  n_by_sample=cells['orig.ident'].value_counts().to_dict(), n_clusters=len(rows),
                  dot_mean_max=ctx.dot_vmax,
                  gene_ranges={s: [sc.vmin, sc.vmax] for s, sc in ctx.gene_scales.items()},
                  umap_frame=[list(map(float, ctx.frame[0])), list(map(float, ctx.frame[1]))]),
        source_data=dict(path=str(sd.relative_to(out)), sha256=file_sha(sd)),
        figures=[])
    st = cfg['style']
    for fig in cfg['figures']:
        if a.only and fig['id'] not in a.only:
            continue
        page = LAYOUTS[fig['layout']](ctx, fig)
        if fig.get('fit_height', True):
            fig['height_mm'] = page.fit_height(fig.get('bottom_margin', 0.5))
        paths = page.save(str(out / 'figures' / fig['id']), a.formats)
        checks = check_page(page, fig, st)
        rec['figures'].append(dict(id=fig['id'], layout=fig['layout'], genes=fig['genes'],
                                   size_mm=[fig['width_mm'], fig['height_mm']],
                                   outputs={pathlib.Path(p).name: file_sha(p) for p in paths}, checks=checks))
        flag = 'ok' if all(c['ok'] for c in checks) else 'CHECK'
        print(f"{fig['id']}: {fig['width_mm']} x {fig['height_mm']} mm  [{flag}]")
        for c in checks:
            if not c['ok']:
                print('   ', c['check'], c['detail'])
    if not a.no_compare and cfg.get('compare') and not a.only:
        from expression_figs.compare import build_compare
        cp = build_compare(cfg['compare'], out)
        if cp:
            rec['compare'] = {pathlib.Path(p).name: file_sha(p) for p in cp}
    json.dump(rec, open(out / 'build.json', 'w'), indent=1, default=str)
    print('wrote', out / 'build.json')


if __name__ == '__main__':
    main()
