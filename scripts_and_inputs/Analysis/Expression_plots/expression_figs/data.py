"""Per-cell table -> cluster rows, per-cluster summaries and violin densities.

Nothing here is a new analysis: every number is a count, a fraction or a mean of the SCT data-slot values that the
original Seurat VlnPlot / FeaturePlot calls drew, grouped by the same cluster column (time.celltype.anno.Li).
"""
from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd

REQUIRED = ['cell', 'orig.ident', 'genotype', 'celltype.anno.Li', 'time.anno.Li', 'time.celltype.anno.Li',
            'umap_1', 'umap_2']


class DataError(ValueError):
    pass


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def load_cells(path, cfg):
    df = pd.read_csv(path)
    miss = [c for c in REQUIRED if c not in df.columns]
    genes = cfg['genes']
    miss += [g['column'] for g in genes.values() if g['column'] not in df.columns]
    if miss:
        raise DataError(f'{path}: missing columns {miss}')
    samples = cfg['samples']
    bad = sorted(set(df['orig.ident']) - set(samples))
    if bad:
        raise DataError(f'unexpected samples {bad}; expected {list(samples)}')
    exp = df['orig.ident'].map(samples)
    if not (exp == df['genotype']).all():
        raise DataError('genotype column disagrees with the sample -> genotype map in config')
    if df['cell'].duplicated().any():
        raise DataError('duplicated cell barcodes')
    for g in genes.values():
        v = df[g['column']]
        if v.isna().any() or (v < 0).any():
            raise DataError(f"{g['column']}: NA or negative values")
    return df


def cluster_rows(df, cfg):
    """Ordered list of row dicts (one per cluster present in the data): tissue then stage, as in config."""
    order = cfg['row_order']
    tissues, stages = order['celltypes'], order['stages']
    pairs = df.groupby(['celltype.anno.Li', 'time.anno.Li', 'time.celltype.anno.Li']).size().reset_index()
    # one cluster name per (cell type, stage) and vice versa
    if pairs['time.celltype.anno.Li'].duplicated().any() or pairs[['celltype.anno.Li', 'time.anno.Li']].duplicated().any():
        raise DataError('time.celltype.anno.Li is not one-to-one with (celltype.anno.Li, time.anno.Li)')
    unknown_t = sorted(set(pairs['celltype.anno.Li']) - set(tissues))
    unknown_s = sorted(set(pairs['time.anno.Li']) - set(stages))
    if unknown_t or unknown_s:
        raise DataError(f'add to config row_order: cell types {unknown_t}, stages {unknown_s}')
    pairs['t'] = pairs['celltype.anno.Li'].map({t: i for i, t in enumerate(tissues)})
    pairs['s'] = pairs['time.anno.Li'].map({s: i for i, s in enumerate(stages)})
    pairs = pairs.sort_values(['t', 's'])
    return [dict(cluster=r['time.celltype.anno.Li'], celltype=r['celltype.anno.Li'], stage=r['time.anno.Li'])
            for _, r in pairs.iterrows()]


def summaries(df, cfg, rows):
    """Tidy per cluster x genotype x gene: n_cells, n_expressing (SCT > 0), pct_expressing, mean_sct (arithmetic mean
    of the SCT data values, zeros included), mean_seurat_dotplot (Seurat DotPlot's 'average expression' without
    scaling: log1p(mean(expm1(x)))), median_sct, max_sct."""
    out = []
    gens = cfg['genotype_order']
    for sym, g in cfg['genes'].items():
        col = g['column']
        for r in rows:
            sub = df[df['time.celltype.anno.Li'] == r['cluster']]
            for gt in gens:
                v = sub.loc[sub['genotype'] == gt, col].to_numpy(float)
                n = len(v)
                out.append(dict(gene=sym, agi=g['agi'], cluster=r['cluster'], celltype=r['celltype'],
                                stage=r['stage'], genotype=gt, n_cells=n, n_expressing=int((v > 0).sum()),
                                pct_expressing=100.0 * (v > 0).mean() if n else np.nan,
                                mean_sct=v.mean() if n else np.nan,
                                mean_seurat_dotplot=np.log1p(np.expm1(v).mean()) if n else np.nan,
                                median_sct=np.median(v) if n else np.nan,
                                max_sct=v.max() if n else np.nan))
    return pd.DataFrame(out)


def bw_nrd0(x):
    """R stats::bw.nrd0 (the bandwidth ggplot2's geom_violin uses by default)."""
    x = np.asarray(x, float)
    n = len(x)
    hi = np.std(x, ddof=1)
    q75, q25 = np.quantile(x, [0.75, 0.25])   # R quantile type 7 = numpy 'linear'
    lo = min(hi, (q75 - q25) / 1.34)
    if not lo > 0:
        lo = hi if hi > 0 else (abs(x[0]) if abs(x[0]) > 0 else 1.0)
    return 0.9 * lo * n ** (-0.2)


def violin_density(x, n_grid=512, adjust=1.0):
    """ggplot2 stat_ydensity(trim = TRUE, bw = 'nrd0', kernel = 'gaussian'): Gaussian KDE on a 512-point grid over
    [min(x), max(x)], evaluated exactly (R's density() uses an FFT approximation; differences are far below print
    resolution). Returns (grid, density) or None when there are fewer than two cells (ggplot2 drops such groups).
    A group whose values are all identical returns (array([v]), array([1.0])) - drawn as a flat line."""
    x = np.asarray(x, float)
    if len(x) < 2:
        return None
    if np.all(x == x[0]):
        return np.array([x[0]]), np.array([1.0])
    bw = bw_nrd0(x) * adjust
    grid = np.linspace(x.min(), x.max(), n_grid)
    # exact Gaussian sum, chunked; ties collapsed with weights for speed
    vals, counts = np.unique(x, return_counts=True)
    dens = np.zeros_like(grid)
    for i in range(0, len(vals), 2048):
        v, c = vals[i:i + 2048], counts[i:i + 2048]
        dens += (c[None, :] * np.exp(-0.5 * ((grid[:, None] - v[None, :]) / bw) ** 2)).sum(1)
    dens /= len(x) * bw * np.sqrt(2 * np.pi)
    return grid, dens


def de_eligibility(df, rows, min_cells):
    """Per cluster: cells per sample and whether every sample has >= min_cells cells. This is the muscat rule used
    for differential expression (pbDS(min_cells = 5) with two samples per genotype): a cluster was tested only if all
    four samples passed. Returns a DataFrame indexed by cluster."""
    ct = pd.crosstab(df['time.celltype.anno.Li'], df['orig.ident'])
    ct = ct.reindex([r['cluster'] for r in rows]).fillna(0).astype(int)
    out = ct.add_prefix('n_')
    out['de_eligible'] = (ct >= min_cells).all(axis=1)
    return out
