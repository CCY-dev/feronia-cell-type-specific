#!/usr/bin/env python
"""Derive the Fig. 1 tables in data/fig1/ from this repository's DE outputs, exactly as the notebooks do,
and verify them against every saved output that exists.

Inputs: the files of this repository listed below (paths relative to the repository root).
Her code is not changed and not re-run here; each step below names the notebook cell it mirrors.
  M = scripts_and_inputs/Analysis/Differential_expression_analysis/scripts/
      20260317_feronia_salinity_timecourse_muscat_newIntegratedSeurat3_NoSuspiciousCluster.ipynb   (Fig. 1b, c)
  H = .../scripts/20260317_feronia_salinity_DE_hormone_geneoverlap_newIntegratedSeurat3_NoSuspiciousCluster.ipynb
                                                                                                   (Fig. 1d)
  R = scripts_and_inputs/Preprocessing/scripts/202501226_Ref_integration_feronia_salinity_time_course_seu3.ipynb
                                                                                                   (Fig. 1a colours)

Usage (from the kit folder):
  python tools/fig1_tables.py             # write data/fig1/*.csv and provenance/fig1/verification.json
  python tools/fig1_tables.py --check     # recompute and compare with data/fig1/ without writing (exit 1 if not equal)

Tables written (all small; the build reads only these):
  fig1a_umap_cells.csv.gz          per-cell UMAP export (copied unchanged; see tools/hpc/extract_fig1a_umap.R)
  fig1a_palette.csv                cell-type and stage colours and drawing order from her notebook R (cells 44, 46)
  fig1b_genes_by_n_subpopulations.csv   M cells 62-64
  fig1c_deg_per_subpopulation.csv  M cells 56-58 (rank = her order, descending total)
  fig1d_hormone_overlap.csv        H cells 7-28: overlap counts and her one-sided Fisher P (as coded, the default)
  fig1d_colour_scale.csv           H cell 29: colorRamp2(c(0, -log10(0.05), 10), c("white", "beige", "red"))
  fig1d_hormone_overlap_q_proposed.csv  OPT-IN alternative (not used in the paper); never drawn unless the config
                                   sets panels.d.statistic: q_proposed
  fig1_summary.csv                 counts quoted in the figure or the legend, with their source
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

KIT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(KIT))
sys.path.insert(0, str(KIT / 'tools'))
from ferfigs.panels.fig1_colour import ColorRamp2, srgb_to_lab, hex_to_rgb01  # noqa: E402

COMMIT = 'd9d516eab98c6d043915c26328d3d4b83f35a3ff'
UP = KIT.parents[2]          # this repository's root (the kit sits in scripts_and_inputs/Analysis/); was a pinned copy
DE = UP / 'scripts_and_inputs' / 'Analysis' / 'Differential_expression_analysis'
NB_M = DE / 'scripts' / '20260317_feronia_salinity_timecourse_muscat_newIntegratedSeurat3_NoSuspiciousCluster.ipynb'
NB_H = DE / 'scripts' / '20260317_feronia_salinity_DE_hormone_geneoverlap_newIntegratedSeurat3_NoSuspiciousCluster.ipynb'
NB_R = UP / 'scripts_and_inputs' / 'Preprocessing' / 'scripts' / \
    '202501226_Ref_integration_feronia_salinity_time_course_seu3.ipynb'
SIG = DE / 'outputs' / 'sig_de_fil.csv'
ALL_BG = DE / 'inputs' / 'all_bg.csv'
# H cell 12: object name -> file (JA is Hickman et al. 2017)
HORMONE_FILES = {'ABA': 'ABA.csv', 'ACC': 'ACC.csv', 'BR': 'BR.csv', 'CK': 'CK.csv', 'GA': 'GA.csv',
                 'IAA': 'IAA.csv', 'JA': 'Hickman_2017_JA.csv', 'SA': 'SA.csv'}
# H cells 10-11 and M cells 53-54 (the same filter in both notebooks)
P_VAL_THRES, FC_THRES = 0.05, 1.5
# H cell 28-29
PSEUDO = 1e-200
CUTOFF = -math.log10(0.05)
RAMP_BREAKS, RAMP_COLOURS = [0.0, CUTOFF, 10.0], ['white', 'beige', 'red']
Q_AUTHORS = DE / 'outputs' / 'hormone_overlap' / 'Fig1d_overlap_q.csv'
OUT = KIT / 'data' / 'fig1'
PROV = KIT / 'provenance' / 'fig1'
UMAP_EXPORT_NAME = 'fig1a_umap_cells.csv.gz'   # per-cell export, shipped unchanged; its SHA-256 goes into fig1_summary.csv


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ------------------------------------------------------------------ reading (as R reads the files)
def read_sig():
    """H cell 7 / M cell 52: read.delim(sep = ',', header = T); then the DEG filter of H cell 11 / M cell 54."""
    sig = pd.read_csv(SIG, dtype={'gene': str, 'cluster_id': str, 'up_dn_label': str}, keep_default_na=True)
    keep = (sig['p_adj.loc'] <= P_VAL_THRES) & (sig['logFC'].abs() > math.log2(FC_THRES))
    return sig, sig[keep].copy()


def read_hormone_lists():
    """H cells 12, 14, 23: read.csv(header = T) (R drops the UTF-8 byte-order mark; values are kept as written, so
    'AT1G56130 ' keeps its trailing space), then as.character(unique(GeneID)) per hormone."""
    out = {}
    for h, fn in HORMONE_FILES.items():
        df = pd.read_csv(DE / 'inputs' / fn, dtype=str, keep_default_na=False, encoding='utf-8-sig')
        assert list(df.columns) == ['GeneID'], (fn, list(df.columns))
        out[h] = list(dict.fromkeys(df['GeneID']))
    return out


def genome_size_from_all_bg():
    """H cell 15: length(unique(all_bg$gene))."""
    bg = pd.read_csv(ALL_BG, usecols=['gene', 'cluster_id'], dtype=str)
    return bg, int(bg['gene'].nunique())


# ------------------------------------------------------------------ Fig. 1b and 1c (M cells 56-64)
def derive_b(deg):
    n_clusters_de = deg.groupby('gene').size()                      # M cell 62: tally per gene
    gene_de_n = n_clusters_de.value_counts()                         # M cell 63: tally per n_clusters_de
    binned = {}
    for n, k in gene_de_n.items():                                   # M cell 64: >= 15 -> '15+'
        lab = '15+' if n >= 15 else str(int(n))
        binned[lab] = binned.get(lab, 0) + int(k)
    levels = [str(i) for i in range(1, 16)] + ['15+']
    rows = [(lab, binned[lab]) for lab in levels if lab in binned]
    b = pd.DataFrame(rows, columns=['n_subpopulations', 'n_genes'])
    total = int(deg['gene'].nunique())                               # M cell 63: total_de_gene_number
    assert int(b['n_genes'].sum()) == total
    return b, total


def derive_c(deg):
    nb = deg.groupby(['cluster_id', 'up_dn_label']).size().unstack(fill_value=0)   # M cell 56
    c = pd.DataFrame({'subpopulation': nb.index, 'up': nb.get('Up', 0).astype(int).values,
                      'down': nb.get('Down', 0).astype(int).values})
    c['total'] = c['up'] + c['down']                                 # M cell 58: total_n
    assert not c['total'].duplicated().any(), 'ties in total_n: fct_reorder tie order would need R collation'
    c = c.sort_values('total', ascending=False, kind='stable').reset_index(drop=True)   # .desc = TRUE
    c['rank'] = np.arange(1, len(c) + 1)
    return c[['rank', 'subpopulation', 'up', 'down', 'total']]


# ------------------------------------------------------------------ Fig. 1d as coded (H cells 16-29)
def gene_overlap_p(n_a, n_b, n_i, n_u, genome):
    """GeneOverlap::testGeneOverlap: fisher.test(matrix(c(genome - |union|, |B| - |I|, |A| - |I|, |I|), ncol = 2),
    alternative = 'greater')$p.value."""
    tbl = [[genome - n_u, n_a - n_i], [n_b - n_i, n_i]]
    return float(fisher_exact(tbl, alternative='greater').pvalue)


def derive_d(deg, hormones, genome):
    markers = {k: list(dict.fromkeys(g['gene'])) for k, g in deg.groupby('cluster_id', sort=True)}   # H 16-20
    rows = []
    for s, A in markers.items():
        sa = set(A)
        for h, B in hormones.items():
            sb = set(B)
            n_i, n_u = len(sa & sb), len(sa | sb)
            p = gene_overlap_p(len(sa), len(sb), n_i, n_u, genome)
            rows.append(dict(subpopulation=s, hormone=h, overlap=n_i, n_deg=len(sa), n_set=len(sb), genome_size=genome,
                             p_value=p, mlog10p=-math.log10(p + PSEUDO)))
    d = pd.DataFrame(rows)
    return d, markers


def derive_q_proposed(d, markers, hormones, bg):
    """OPT-IN alternative (not used in the paper): universe of subpopulation s =
    genes tested in s (all_bg rows of s); hormone sets whitespace-stripped and restricted to that universe; one-sided
    Fisher exact test (greater); Benjamini-Hochberg across the 272 tests. Never used unless the config selects it."""
    from statsmodels.stats.multitest import multipletests
    elig = {k: set(g['gene']) for k, g in bg.groupby('cluster_id')}
    hs = {h: {x.strip() for x in B} for h, B in hormones.items()}
    out = []
    for r in d.itertuples():
        E, A = elig[r.subpopulation], set(markers[r.subpopulation])
        assert A <= E, f'DEGs outside the tested genes of {r.subpopulation}'
        S = hs[r.hormone] & E
        k = len(A & S)
        assert k == r.overlap, (r.subpopulation, r.hormone, k, r.overlap)
        p = float(fisher_exact([[len(E) - len(A | S), len(S) - k], [len(A) - k, k]], alternative='greater').pvalue)
        out.append(dict(subpopulation=r.subpopulation, hormone=r.hormone, overlap=k, n_deg=len(A),
                        n_universe=len(E), n_set_in_universe=len(S), p_proposed=p))
    q = pd.DataFrame(out)
    q['q_proposed'] = multipletests(q['p_proposed'], method='fdr_bh')[1]
    q['mlog10q'] = -np.log10(q['q_proposed'] + PSEUDO)
    return q


# ------------------------------------------------------------------ Fig. 1a colours (R cells 44, 46)
def _r_vector(src, name):
    m = re.search(rf'{re.escape(name)}\s*(?:<-|=)\s*c\((.*?)\)', src, re.S)
    if not m:
        raise ValueError(f'no vector {name} in notebook source')
    return re.findall(r'["\']([^"\']*)["\']', m.group(1))


def derive_palette():
    nb = json.loads(NB_R.read_text())
    src44 = ''.join(nb['cells'][44]['source'])
    src46 = ''.join(nb['cells'][46]['source'])
    order, palette = _r_vector(src44, 'order'), _r_vector(src44, 'palette')
    assert len(order) == len(palette) == 17
    rows = []
    # cell 44: factor levels = order[sort(match(unique(x), order))]; colour of a level = palette at the same position
    for i, (name, col) in enumerate(zip(order, palette)):
        rows.append(dict(map='celltype', name=name, hex=col.upper(), level=i + 1, draw_order='object',
                         source='R notebook cell 44 (order, palette)'))
    # cell 46: DimPlot(group.by = 'time.anno.Li', order = o, cols = k). Seurat SingleDimPlot sets the factor levels to
    # rev(o) and sorts the cells by level, so level 1 is drawn first (bottom); scale_color_manual(values = k) gives
    # level i the colour k[i].
    dm = re.search(r'DimPlot\((.*)\)', src46, re.S).group(1)
    o = re.findall(r'["\']([^"\']*)["\']', re.search(r'order\s*=\s*c\((.*?)\)', dm, re.S).group(1))
    k = re.findall(r'["\']([^"\']*)["\']', re.search(r'cols\s*=\s*c\((.*?)\)', dm, re.S).group(1))
    assert len(o) == len(k) == 8
    for i, (name, col) in enumerate(zip(o[::-1], k)):
        rows.append(dict(map='stage', name=name, hex=col.upper(), level=i + 1, draw_order=str(i + 1),
                         source='R notebook cell 46 (DimPlot order, cols)'))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ saved outputs of her notebooks
def _plain(nb, cell, k=None):
    outs = nb['cells'][cell].get('outputs', [])
    texts = [''.join(o.get('data', {}).get('text/plain', '')) for o in outs]
    return texts if k is None else texts[k]


def parse_r_list(txt):
    """R print() of a named list of character vectors -> {name: [values]} (quoted strings kept exactly)."""
    out, cur = {}, None
    for line in txt.splitlines():
        m = re.match(r'^\$`?(.*?)`?\s*$', line)
        if m and line.startswith('$'):
            cur = m.group(1); out[cur] = []
            continue
        if cur is not None:
            out[cur] += re.findall(r'"((?:[^"\\]|\\.)*)"', line)
    return out


def parse_r_table(txt):
    """R print() of a named integer vector (table) -> {name: int}."""
    lines = [l for l in txt.splitlines() if l.strip()]
    out = {}
    for names, vals in zip(lines[0::2], lines[1::2]):
        # names are right-aligned over their values: split the value line, then cut the name line at the same ends
        ends = [m.end() for m in re.finditer(r'\d+', vals)]
        starts = [0] + ends[:-1]
        nums = [int(v) for v in re.findall(r'\d+', vals)]
        for s, e, v in zip(starts, ends, nums):
            out[names[s:e].strip()] = v
    return out


def saved_output_checks(b, total, c, d, markers, hormones, genome, umap, deg):
    """Compare the derived tables with every saved output of her notebooks that holds the same numbers."""
    checks = []

    def add(name, ok, detail):
        checks.append(dict(check=name, passed=bool(ok), detail=detail))
    nbh = json.loads(NB_H.read_text())
    add('H cell 15: genome_size', int(_plain(nbh, 15, 0).split()[-1]) == genome, f'printed {_plain(nbh, 15, 0).strip()}')
    add('H cell 17: length(markers_list) = 34', _plain(nbh, 17, 0).split()[-1] == str(len(markers)), str(len(markers)))
    hs = parse_r_table(_plain(nbh, 25, 0))
    add('H cell 25: hormone set sizes', hs == {h: len(v) for h, v in hormones.items()}, hs)
    ml = parse_r_list(_plain(nbh, 21, 0))
    same_m = set(ml) == set(markers) and all(ml[k] == markers[k] for k in ml)
    add('H cell 21: printed markers_list equals the DEG lists from sig_de_fil.csv (names, genes and order)', same_m,
        f'{len(ml)} lists, {sum(map(len, ml.values()))} genes')
    hl = parse_r_list(_plain(nbh, 24, 0))
    same_h = set(hl) == set(hormones) and all(hl[k] == hormones[k] for k in hl)
    add('H cell 24: printed hormone_list equals the vendored gene sets (exact strings, order)', same_h,
        f'{len(hl)} sets; SA keeps "AT1G56130 " with a trailing space: {"AT1G56130 " in hl.get("SA", [])}')
    # overlaps and P recomputed from the printed lists alone (her session's objects)
    mism = []
    for r in d.itertuples():
        A, B = set(ml[r.subpopulation]), set(hl[r.hormone])
        n_i = len(A & B)
        p = gene_overlap_p(len(A), len(B), n_i, len(A | B), int(_plain(nbh, 15, 0).split()[-1]))
        if n_i != r.overlap or abs(p - r.p_value) > 1e-12 * max(p, 1e-300):
            mism.append((r.subpopulation, r.hormone, n_i, r.overlap, p, r.p_value))
    add('H cells 21+24+15: all 272 overlaps and P values recomputed from her printed lists', not mism,
        f'{len(d) - len(mism)} of {len(d)} identical' + (f'; first mismatches {mism[:3]}' if mism else ''))
    nbm = json.loads(NB_M.read_text())
    add('M cell 63: total_de_gene_number', _plain(nbm, 63, 0).split()[-1] == str(total), _plain(nbm, 63, 0).strip())
    head_b = _plain(nbm, 65, 0)
    exp = [(m.group(1), int(m.group(2))) for m in re.finditer(r'^\d+\s+(\S+)\s+(\d+)\s*$', head_b, re.M)]
    bmap = dict(zip(b['n_subpopulations'], b['n_genes']))
    add('M cell 65: head(gene_de_binned, 3)', all(bmap.get(k) == v for k, v in exp) and len(exp) == 3, exp)
    head_c = _plain(nbm, 57, 0)
    expc = [(m.group(1).strip(), m.group(2), int(m.group(3)))
            for m in re.finditer(r'^\d+\s+(.+?)\s+(Up|Down)\s+(\d+)\s*$', head_c, re.M)]
    cm = c.set_index('subpopulation')
    add('M cell 57: head(nb_genes_per_clust, 3)',
        len(expc) == 3 and all(int(cm.loc[s, 'up' if ud == 'Up' else 'down']) == v for s, ud, v in expc), expc)
    tu = re.findall(r'"(\d+) \\n Total DE genes"', ''.join(_plain(nbm, 73)))
    deg_up = int(deg['gene'][deg['up_dn_label'] == 'Up'].nunique())
    deg_dn = int(deg['gene'][deg['up_dn_label'] == 'Down'].nunique())
    add('M cell 73: unique genes ever up / ever down', tu == [str(deg_up), str(deg_dn)], f'printed {tu}; derived {[deg_up, deg_dn]}')
    # cell counts of the per-cell UMAP export against her saved tables (same object)
    t40 = parse_r_table(_plain(nbm, 40, 0)); t41 = parse_r_table(_plain(nbm, 41, 0)); t44 = parse_r_table(_plain(nbm, 44, 0))
    t36 = parse_r_table(_plain(nbm, 36, 0))
    u_stage = umap['time.anno.Li'].value_counts().to_dict()
    u_ct = umap['celltype.anno.Li'].value_counts().to_dict()
    u_gen = umap['genotype'].map({'WT': 'WT', 'fer4': 'fer_4'}).value_counts().to_dict()
    rootcap = {'Distal Columella', 'Proximal Columella', 'Distal Lateral Root Cap', 'Proximal Lateral Root Cap'}
    tc = np.where(umap['time.anno.Li'].isin(rootcap), umap['time.anno.Li'],
                  umap['time.anno.Li'] + '_' + umap['celltype.anno.Li'])
    u_sub = pd.Series(tc).value_counts().to_dict()
    add('M cells 36, 40, 41, 44: cells per subpopulation, stage, genotype and cell type in the UMAP export',
        u_stage == t40 and u_gen == t41 and u_ct == t44 and u_sub == t36,
        f'{len(umap)} cells; stage {u_stage == t40}, genotype {u_gen == t41}, cell type {u_ct == t44}, '
        f'subpopulation ({len(t36)}) {u_sub == t36}')
    return checks


def transcription_checks(b, c, d):
    """Numbers read off the plots saved in her notebooks (provenance/fig1/saved_plot_transcriptions.json: one manual
    reading and one independent vision-model reading of each plot)."""
    p = PROV / 'saved_plot_transcriptions.json'
    if not p.exists():
        return []
    T = json.loads(p.read_text())
    checks = []
    bl = list(zip(b['n_subpopulations'], b['n_genes']))
    cm = c.set_index('subpopulation')
    dm = d.set_index(['subpopulation', 'hormone'])['overlap']
    for reader in ('manual', 'vision_model'):
        R = T.get(reader, {})
        if 'M_cell66_bins' in R:
            got = [tuple(x) for x in R['M_cell66_bins']]
            checks.append(dict(check=f'M cell 66 plot (Fig. 1b counts), {reader} reading', passed=got == bl,
                               detail=f'{sum(1 for a, b_ in zip(got, bl) if a == b_)} of {len(bl)} bars equal'))
        for key, col in (('M_cell70_up', 'up'), ('M_cell70_down', 'down')):
            if key in R:
                rr = R[key]
                ok = [int(cm.loc[n, col]) == v for n, v in rr if n in cm.index]
                checks.append(dict(check=f'M cell 70 plot ({col}), {reader} reading',
                                   passed=len(rr) == len(cm) and all(ok) and len(ok) == len(rr),
                                   detail=f'{sum(ok)} of {len(cm)} bars equal'))
        if 'H_cell30_heatmap' in R:
            H = R['H_cell30_heatmap']
            ok = n = 0
            bad = []
            for row in H['rows']:
                for h, v in zip(H['columns'], row['values']):
                    n += 1
                    if dm.get((row['subpopulation'], h)) == v:
                        ok += 1
                    else:
                        bad.append((row['subpopulation'], h, v, dm.get((row['subpopulation'], h))))
            checks.append(dict(check=f'H cell 30 heatmap (Fig. 1d counts), {reader} reading', passed=ok == n == len(d),
                               detail=f'{ok} of {n} cells equal' + (f'; differences {bad[:5]}' if bad else '')))
    return checks


def heatmap_colour_check(d):
    """Tile colours of her saved heatmap (H cell 30 PNG) against the colour her col_fun gives for the derived P.
    Rows are located by the counts (each row's 8 counts identify it); colours are the median of tile pixels away
    from the printed number."""
    nb = json.loads(NB_H.read_text())
    import base64
    from PIL import Image
    png = [o['data']['image/png'] for o in nb['cells'][30]['outputs'] if 'image/png' in o.get('data', {})][0]
    im = np.asarray(Image.open(io.BytesIO(base64.b64decode(png))).convert('RGB')).astype(int)
    T = json.loads((PROV / 'saved_plot_transcriptions.json').read_text())['vision_model']['H_cell30_heatmap']
    # grid: the 1-px tile borders (rect_gp col = 'lightgrey', blended with the fill) differ from both neighbours
    # along most of the heatmap body; columns first, then rows within the body
    def lines(a):
        dd = np.abs(np.diff(a, axis=1)).sum(-1)
        s = ((dd[:, :-1] > 20) & (dd[:, 1:] > 20)).mean(0)
        return [int(i + 1) for i in np.where(s > 0.6)[0]]
    xe = lines(im)
    ye = lines(im[:, xe[0]:xe[-1] + 1].transpose(1, 0, 2)) if len(xe) >= 2 else []
    if len(xe) != 9 or len(ye) != 35:
        return dict(check='H cell 30 heatmap tile colours', passed=False,
                    detail=f'grid not found ({len(xe)} x edges, {len(ye)} y edges)')
    ramp = ColorRamp2(RAMP_BREAKS, RAMP_COLOURS)
    dd = d.set_index(['subpopulation', 'hormone'])
    dE, exact, n = [], 0, 0
    for i, row in enumerate(T['rows']):
        for j, h in enumerate(T['columns']):
            x0, x1, y0, y1 = xe[j] + 3, xe[j + 1] - 2, ye[i] + 3, ye[i + 1] - 2
            w = x1 - x0
            patch = np.concatenate([im[y0:y1, x0:x0 + w // 5].reshape(-1, 3), im[y0:y1, x1 - w // 5:x1].reshape(-1, 3)])
            obs = np.median(patch, 0)
            exp = np.array(hex_to_rgb01(ramp([dd.loc[(row['subpopulation'], h), 'mlog10p']])[0])) * 255
            e = float(np.linalg.norm(srgb_to_lab(obs / 255) - srgb_to_lab(exp / 255)))
            dE.append(e); n += 1
            exact += int(np.all(np.abs(obs - exp) <= 1))
    dE = np.array(dE)
    return dict(check='H cell 30 heatmap tile colours vs her col_fun applied to the derived P (as coded)',
                passed=bool(dE.max() < 2.0), detail=f'{n} tiles; {exact} equal within 1/255 per channel; CIE76 dE '
                f'median {np.median(dE):.2f}, max {dE.max():.2f}')


def r_rerun_check(d):
    p = PROV / 'r_hormone_overlap_as_coded.csv'
    if not p.exists():
        return []
    r = pd.read_csv(p, dtype={'pval': str, 'p_val_log': str})
    m = d.merge(r, on=['subpopulation', 'hormone'], how='outer', indicator=True)
    both = m['_merge'] == 'both'
    pr = m['pval'].astype(float)
    rel = (np.abs(m['p_value'] - pr) / np.maximum(pr, 1e-300))[both]
    ramp = ColorRamp2(RAMP_BREAKS, RAMP_COLOURS)
    py_col = ramp(m.loc[both, 'mlog10p'].values)
    return [dict(check='R re-run of H cells 7-29 (tools/fig1_verify_r.R): overlaps, n_deg, n_set',
                 passed=bool(both.all() and (m['overlap'] == m['intersection']).all() and (m['n_deg_x'] == m['n_deg_y']).all()
                             and (m['n_set_x'] == m['n_set_y']).all()),
                 detail=f'{int(both.sum())} of {len(d)} cells matched'),
            dict(check='R re-run: one-sided Fisher P (scipy vs R fisher.test)', passed=bool(rel.max() < 1e-9),
                 detail=f'max relative difference {rel.max():.2e}'),
            dict(check='R re-run: heatmap colours (Python ColorRamp2 vs circlize colorRamp2 hex)',
                 passed=list(py_col) == list(m.loc[both, 'colour']),
                 detail=f'{sum(a == b for a, b in zip(py_col, m.loc[both, "colour"]))} of {int(both.sum())} identical')]


def candidate_checks(b, total, c, d):
    """Numbers printed in the current candidate figure (provenance/fig1/candidate_reference_values.json, extracted from
    the vector text of Fig1_candidate.pdf by tools/seed_fig1_reference_values.py)."""
    p = PROV / 'candidate_reference_values.json'
    if not p.exists():
        return []
    V = json.loads(p.read_text())
    out = []
    out.append(dict(check='current figure b: 15 bar labels', passed=V['b_counts'] == b['n_genes'].tolist(),
                    detail=f"printed {V['b_counts']}"))
    one_two = int(b['n_genes'].iloc[:2].sum())
    out.append(dict(check='current figure b: total and one-or-two-subpopulation share',
                    passed=V['b_total'] == total and V['b_one_two'] == one_two and V['b_pct'] == round(100 * one_two / total),
                    detail=f"printed {V['b_total']:,} / {V['b_one_two']:,} ({V['b_pct']}%)"))
    out.append(dict(check='current figure c: 34 totals in the same order', passed=V['c_totals'] == c['total'].tolist(),
                    detail=f"{sum(a == b_ for a, b_ in zip(V['c_totals'], c['total']))} of {len(c)} equal"))
    out.append(dict(check='current figure c/d: row order = her ranking by total DEGs (M cell 58)',
                    passed=V['c_order'] == c['subpopulation'].tolist(),
                    detail=f"{sum(a == b_ for a, b_ in zip(V['c_order'], c['subpopulation']))} of {len(c)} rows in place"))
    dm = d.set_index(['subpopulation', 'hormone'])['overlap']
    cells = V['d_counts']
    ok = sum(1 for s, row in cells.items() for h, v in row.items() if dm.get((s, h)) == v)
    n = sum(len(r) for r in cells.values())
    out.append(dict(check='current figure d: 272 printed counts', passed=ok == n == len(d), detail=f'{ok} of {n} equal'))
    if 'd_q_reference' in V:
        q = V['d_q_reference']
        out.append(dict(check='opt-in q: identical to the candidate run (source_data/fig1/fig1d_hormone_overlap.csv)',
                        passed=bool(q.get('identical')), detail=q.get('detail', '')))
    return out


# ------------------------------------------------------------------ main
def derive():
    """All tables (as CSV text) and all checks. Needs every pinned file, including the fetch-only ones."""
    sig, deg = read_sig()
    hormones = read_hormone_lists()
    bg, genome = genome_size_from_all_bg()
    b, total = derive_b(deg)
    c = derive_c(deg)
    d, markers = derive_d(deg, hormones, genome)
    q = derive_q_proposed(d, markers, hormones, bg)
    pal = derive_palette()
    umap_path = OUT / UMAP_EXPORT_NAME
    umap = pd.read_csv(umap_path, dtype=str, keep_default_na=False)
    scale = pd.DataFrame(dict(value=RAMP_BREAKS, colour=RAMP_COLOURS,
                              hex=[f'#{r:02X}{g:02X}{b_:02X}' for r, g, b_ in
                                   [tuple(int(v) for v in np.round(np.array(hex_to_rgb01(h)) * 255))
                                    for h in ColorRamp2(RAMP_BREAKS, RAMP_COLOURS)(RAMP_BREAKS)]],
                              note=['0', '-log10(0.05) (cutoff)', 'cap: values above 10 take this colour']))
    n_wt = int((umap['genotype'] == 'WT').sum()); n_fer = int((umap['genotype'] == 'fer4').sum())
    one_two = int(b['n_genes'].iloc[:2].sum())
    summary = pd.DataFrame([
        ('n_cells', len(umap), 'rows of fig1a_umap_cells.csv.gz (= M cell 41: 21,655 WT + 23,009 fer-4)'),
        ('n_cells_WT', n_wt, 'genotype == WT'), ('n_cells_fer4', n_fer, 'genotype == fer4'),
        ('n_deg_rows', len(deg), 'sig_de_fil.csv rows passing p_adj.loc <= 0.05 and |logFC| > log2(1.5)'),
        ('n_deg_rows_removed_by_filter', len(sig) - len(deg), 'the filter of H cell 11 / M cell 54 removes nothing'),
        ('n_deg_unique', total, 'M cell 63 total_de_gene_number'),
        ('n_deg_one_or_two_subpopulations', one_two, 'bins 1 and 2 of fig1b'),
        ('n_subpopulations_with_deg', len(c), 'clusters in sig_de_fil.csv (H cell 17)'),
        ('n_subpopulations_tested', int(bg['cluster_id'].nunique()), 'clusters in all_bg.csv'),
        ('genome_size', genome, 'H cell 15: length(unique(all_bg$gene))'),
        ('n_tests_d', len(d), 'subpopulations x hormone sets'),
        ('n_p_lt_0.05_as_coded', int((d['p_value'] < 0.05).sum()), 'unadjusted one-sided Fisher P (as coded)'),
        ('n_mlog10p_capped', int((d['mlog10p'] >= 10).sum()), 'cells at the top colour (-log10 P >= 10)'),
        ('n_q_lt_0.05_proposed', int((q['q_proposed'] < 0.05).sum()), 'OPT-IN proposal only'),
        ('n_mlog10q_capped_proposed', int((q['mlog10q'] >= 10).sum()), 'OPT-IN proposal only'),
        ('upstream_commit', COMMIT, 'CCY-dev/feronia-cell-type-specific'),
        ('umap_export_sha256', sha256(umap_path), 'data/fig1/fig1a_umap_cells.csv.gz (tools/hpc/extract_fig1a_umap.R)'),
    ], columns=['key', 'value', 'source'])
    tables = {
        'fig1a_palette.csv': pal,
        'fig1b_genes_by_n_subpopulations.csv': b,
        'fig1c_deg_per_subpopulation.csv': c,
        'fig1d_hormone_overlap.csv': d,
        'fig1d_colour_scale.csv': scale,
        'fig1d_hormone_overlap_q_proposed.csv': q,
        'fig1_summary.csv': summary,
    }
    texts = {}
    for name, df in tables.items():
        buf = io.StringIO()
        # colour breaks at full precision (her cutoff is -log10(0.05) exactly); statistics at 10 significant digits
        df.to_csv(buf, index=False, float_format='%.17g' if name == 'fig1d_colour_scale.csv' else '%.10g',
                  lineterminator='\n')
        texts[name] = buf.getvalue()
    checks = []
    checks += saved_output_checks(b, total, c, d, markers, hormones, genome, umap, deg)
    checks += transcription_checks(b, c, d)
    if (PROV / 'saved_plot_transcriptions.json').exists():
        checks.append(heatmap_colour_check(d))
    checks += r_rerun_check(d)
    checks += candidate_checks(b, total, c, d)
    return texts, checks


def main(argv=None):
    ap = argparse.ArgumentParser(description='Derive and verify the Fig. 1 tables')
    ap.add_argument('--check', action='store_true', help='compare with data/fig1/ without writing')
    a = ap.parse_args(argv)
    need = [SIG, ALL_BG, NB_M, NB_H, NB_R, Q_AUTHORS] + [DE / 'inputs' / f for f in HORMONE_FILES.values()]
    probs = [str(p) for p in need if not p.exists()]
    if probs:
        print('Repository files missing:\n  ' + '\n  '.join(probs), file=sys.stderr)
        return 2
    texts, checks = derive()
    # Fig. 1d colour statistic (q_authors) = this repository's BH table, copied unchanged
    # (scripts/20261020_fig1d_hormone_overlap_tested_universe_BH.R; universe 17,554 genes, BH over 272 tests)
    texts['fig1d_hormone_overlap_q_authors.csv'] = Q_AUTHORS.read_text()
    diff = []
    for name, txt in texts.items():
        p = OUT / name
        if a.check:
            if not p.exists() or p.read_text() != txt:
                diff.append(name)
        else:
            p.write_text(txt)
    PROV.mkdir(parents=True, exist_ok=True)
    rep = dict(upstream_commit=COMMIT, tables={n: hashlib.sha256(t.encode()).hexdigest() for n, t in texts.items()},
               checks=checks, n_passed=sum(c_['passed'] for c_ in checks), n_checks=len(checks))
    if not a.check:
        (PROV / 'verification.json').write_text(json.dumps(rep, indent=1, default=str) + '\n')
    for c_ in checks:
        print(('PASS  ' if c_['passed'] else 'FAIL  ') + c_['check'] + '  |  ' + str(c_['detail'])[:160])
    print(f"{rep['n_passed']} of {rep['n_checks']} checks passed")
    if a.check and diff:
        print('data/fig1 differs from a fresh derivation: ' + ', '.join(diff), file=sys.stderr)
        return 1
    return 0 if rep['n_passed'] == rep['n_checks'] else 1


if __name__ == '__main__':
    sys.exit(main())
