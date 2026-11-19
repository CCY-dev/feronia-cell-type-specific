"""Reshape the paper's Source Data workbook (one sheet per wet-lab panel) into tidy long tables.

One row per measurement: sheet, panel, group, series, replicate, value, sample_id, cell.
`cell` is the source cell (e.g. B3), so every plotted number traces back to the workbook.
Nothing is typed in and nothing is imputed (empty cells are skipped). Derived columns in the workbook ("Ave",
"Firefly/Renilla", per-experiment means / normalised rows) are handled explicitly per sheet below.

usage (from the kit folder):
  python tools/source_data_to_long.py "<Source data.xlsx>" [out_dir]      (default out_dir: source_data/long)
writes <out_dir>/<slug>.csv per sheet + ALL_long.csv + summary.csv (n, mean, sd per group x series).
The Source Data workbook is published with the paper; it is not part of this repository.
"""
import csv, os, re, sys, statistics as stt
from collections import defaultdict
import openpyxl
from openpyxl.utils import get_column_letter as L

# sheet -> (panel id(s), unit, parser)
SHEETS = {
    'Ext_Fig_2b Up_qPCR': ('ED2b', 'relative expression', 'qpcr_genes'),
    'Ext_Fig_2c Down_qPCR': ('ED2c', 'relative expression', 'qpcr_genes'),
    'Ext_Fig_2e pWOL_AUX1-Venus': ('ED2e', 'a.u.', 'aux1_wol'),
    'Fig_2d&Ext_Fig_4a Meristem#': ('Fig2d;ED4a', 'cells', 'left_right_ave'),
    'Ext_Fig_4b Meristem length(uM)': ('ED4b', 'µm', 'left_right_ave'),
    'Fig_2f&Ext_Fig_4d Root hair(mm)': ('Fig2f;ED4d', 'mm', 'root_hair'),
    'Fig_2g RALF_CRISPR-L1(mm_': ('Fig2g', 'mm', 'two_row'),
    'Ext_Fig_4f RALF_CRISPR-L3': ('ED4f', 'mm', 'two_row'),
    'Fig_3b pYUC5-GFP_TimeCourse': ('Fig3b', 'a.u.', 'two_row'),
    'Fig_3d pYUC5-GFP_fer4': ('Fig3d', 'a.u.', 'two_row'),
    'Fig_3g pYUC5-GFP_CRISPR': ('Fig3g', 'a.u.', 'two_row'),
    'Ext_Fig_5g RALF_yuc (mm)': ('ED5g', 'mm', 'two_row'),
    'Fig_4b DR5-GFP_fer4': ('Fig4b', 'a.u.', 'two_row'),
    'Fig_4d DR5-GFP_CRISPR': ('Fig4d', 'a.u.', 'two_row'),
    'Fig_4e PEO-IAA_CRISPR-L3': ('Fig4e', 'mm', 'two_row'),
    'Fig_4g RALF_afb1(mm)': ('Fig4g', 'mm', 'two_row'),
    'Ext_Fig_6c PEO-IAA fer-4': ('ED6c', 'mm', 'two_row'),
    'Ext_Fig_7b MultiRALFs_fer4': ('ED7b', 'mm', 'two_row'),
    'Ext_Fig_7d pYUC5-GFP_MultiRALFS': ('ED7d', 'a.u.', 'one_row_series'),
    'Ext_Fig_7f DR5_MultiRALFs': ('ED7f', 'a.u.', 'one_row_series'),
    'Fig_5d LUC assay': ('Fig5d', 'Firefly/Renilla', 'luc'),
    'Fig_6a RALF_dof': ('Fig6a', 'mm', 'two_row'),
    'Fig_6d pYUC5_dof 6h': ('Fig6d', 'a.u.', 'two_row'),
    'Fig_6e DOFs qPCR': ('Fig6e', 'relative expression', 'qpcr_biorep_blocks'),   # 3 biological replicates per gene
    'Fig_6g DOFs-GFP': ('Fig6g', 'a.u.', 'dofs_gfp'),
    'Ext_Fig_9d 35S-DOF RALF': ('ED9d', 'mm', 'two_row'),
    'Ext_Fig_9e DOF qPCR': ('ED9e', 'relative expression', 'qpcr_geno_blocks'),
    'Ext Fig 10c PI puncta': ('ED10c', 'puncta', 'pi_puncta'),
    'Supplementary Fig_3': ('SF3', 'relative expression', 'qpcr_geno_blocks'),
}
# sheets that are not wet-lab measurements (the Fig. 5a TF-family enrichment tables): skipped
SKIP = {'Fig_5a metadata', 'Fig_5a background_1301_TFs', 'Fig_5a enrichment_results'}

# spelling normalisation to the kit's labels (config/labels.yaml uses these)
NORM = {'WT': 'Col-0', 'fer4': 'fer-4', 'pCO2-KO-L1': 'pCO2-KO_L1', 'GUSGFP': 'GUS-GFP', 'GUSGFP ': 'GUS-GFP'}


def norm(s):
    s = str(s).strip()
    s = NORM.get(s, s)
    s = s.replace('μ', 'µ').replace('uM', 'µM')
    s = re.sub(r'^1µM RALF23$', '1 µM RALF23', s)
    return s


def series_norm(s):
    s = norm(s)
    # treatment names: "1 µM RALF23" -> RALF23 (dose is constant within every sheet; kept in notes)
    m = {'1 µM RALF23': 'RALF23', '1 µM RALF': 'RALF23', 'RALF': 'RALF23', '1 µM RALF23 & 4µM PEO-IAA': 'RALF23 & PEO-IAA',
         'RALF23 & 4µM PEO-IAA': 'RALF23 & PEO-IAA', 'Mock ': 'Mock',
         '1 µM RALF1': 'RALF1', '1 µM RALF31': 'RALF31', '1 µM RALF33': 'RALF33'}
    s = s.replace('1µM', '1 µM').replace('4 µM', '4µM')
    return m.get(s, s).strip()


def num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def is_summary(ws, r, c=1):
    """True for the workbook's own summary rows ('Ave', 'Ratio to Col-0', ...): a text label in the index column."""
    a = ws.cell(r, c).value
    return isinstance(a, str) and not a.startswith('Root')


def p_two_row(ws, start_col=2):
    """row1 = group (first column of a block), row2 = series; data from row 3; column A = replicate index."""
    out = []
    grp = None
    for c in range(start_col, ws.max_column + 1):
        g, s = ws.cell(1, c).value, ws.cell(2, c).value
        if g is not None:
            grp = norm(g)
        if s is None:
            if g is None:
                grp = None  # spacer column ends a block
            continue
        for r in range(3, ws.max_row + 1):
            v = ws.cell(r, c).value
            if num(v) and not is_summary(ws, r):
                out.append(dict(group=grp, series=series_norm(s), replicate=str(ws.cell(r, 1).value or r - 2),
                                value=float(v), cell=f'{L(c)}{r}'))
    return out


def p_one_row_series(ws):
    """row1 = series only (single genotype, Col-0)."""
    out = []
    for c in range(2, ws.max_column + 1):
        s = ws.cell(1, c).value
        if s is None:
            continue
        for r in range(2, ws.max_row + 1):
            v = ws.cell(r, c).value
            if num(v) and not is_summary(ws, r):
                out.append(dict(group='Col-0', series=series_norm(s), replicate=str(ws.cell(r, 1).value or r - 1),
                                value=float(v), cell=f'{L(c)}{r}'))
    return out


def p_left_right_ave(ws):
    """Meristem sheets: row1 group per 3 columns (Left, Right, Ave). value = Ave (mean of the two sides of one root),
    as plotted; Left/Right kept in notes. Ave is recomputed and checked against the workbook."""
    out = []
    for c in range(2, ws.max_column + 1):
        if ws.cell(2, c).value != 'Left':
            continue
        g = norm(ws.cell(1, c).value)
        for r in range(3, ws.max_row + 1):
            lv, rv, av = ws.cell(r, c).value, ws.cell(r, c + 1).value, ws.cell(r, c + 2).value
            if not num(av) or is_summary(ws, r):
                continue
            if num(lv) and num(rv):
                assert abs((lv + rv) / 2 - av) < 1e-6, (ws.title, r, c)
            out.append(dict(group=g, series='all', replicate=str(ws.cell(r, 1).value or r - 2), value=float(av),
                            cell=f'{L(c + 2)}{r}', notes=f'left={lv}; right={rv}'))
    return out


def p_aux1_wol(ws):
    out = []
    # AUX1-Venus: Col-0 B..D (Left Right Ave), fer4 E..G; value = Ave
    for c, g in ((2, 'Col-0'), (5, 'fer-4')):
        for r in range(3, ws.max_row + 1):
            av = ws.cell(r, c + 2).value
            if num(av) and not is_summary(ws, r):
                out.append(dict(group=g, series='pAUX1:AUX1-Venus', replicate=str(ws.cell(r, 1).value), value=float(av),
                                cell=f'{L(c + 2)}{r}', notes=f'left={ws.cell(r, c).value}; right={ws.cell(r, c + 1).value}'))
    # pWOL-erGFP: I = index, J = Col-0, K = fer-4, data from row 2
    for c, g in ((10, 'Col-0'), (11, 'fer-4')):
        for r in range(2, ws.max_row + 1):
            v = ws.cell(r, c).value
            if num(v) and not is_summary(ws, r, 9):
                out.append(dict(group=g, series='pWOL:erGFP', replicate=str(ws.cell(r, 9).value), value=float(v),
                                cell=f'{L(c)}{r}'))
    return out


def p_root_hair(ws):
    """Two blocks (L1 lines: B..I; L3 lines: K..R). Column A marks 'Root k' at the first hair of each root;
    rows below belong to that root until the next label. value = one root hair; sample_id = root."""
    out = []
    root = None
    roots = {}
    for r in range(2, ws.max_row + 1):
        a = ws.cell(r, 1).value
        if a is not None and str(a).startswith('Root'):
            root = str(a).replace(' ', '')
        roots[r] = root
    for c in range(2, ws.max_column + 1):
        g = ws.cell(1, c).value
        if g is None:
            continue
        block = 'L1 set' if c < 10 else 'L3 set'
        k = 0
        for r in range(2, ws.max_row + 1):
            v = ws.cell(r, c).value
            if num(v) and not is_summary(ws, r):
                k += 1
                out.append(dict(group=norm(g), series=block, replicate=f'{roots[r]}_h{k}', value=float(v),
                                sample_id=roots[r], cell=f'{L(c)}{r}'))
    return out


def p_qpcr_genes(ws):
    """row2 = header (Gene ID, Col-0-rep1..3, fer-4-rep1..3); one row per gene."""
    out = []
    hdr = [ws.cell(2, c).value for c in range(1, ws.max_column + 1)]
    for r in range(3, ws.max_row + 1):
        gene = ws.cell(r, 1).value
        if gene is None:
            continue
        for c in range(2, ws.max_column + 1):
            v = ws.cell(r, c).value
            if num(v):
                g, rep = hdr[c - 1].rsplit('-rep', 1)
                out.append(dict(group=str(gene), series=norm(g), replicate=f'rep{rep}', value=float(v), cell=f'{L(c)}{r}'))
    return out


def p_qpcr_geno_blocks(ws):
    """Blocks: a header row 'GENE, Mock-rep1.., RALF23-rep1..' then rows Col-0 / fer-4."""
    out = []
    gene, hdr = None, None
    for r in range(1, ws.max_row + 1):
        a = ws.cell(r, 1).value
        b = ws.cell(r, 2).value
        if isinstance(b, str) and b.startswith('Mock-rep'):
            gene = str(a)
            hdr = [ws.cell(r, c).value for c in range(1, ws.max_column + 1)]
            continue
        # older 'Supplementary Fig_3' layout: gene headers 'CDF4_rep2' (= repeat experiment), rows 'Col-0_rep3' / 'fer-4-rep3'
        if gene and isinstance(a, str) and re.sub(r'[_-]rep\d+$', '', a) in ('Col-0', 'fer-4'):
            m = re.match(r'^(.*?)_rep(\d+)$', gene)
            g0, expt = (m.group(1), f'expt{m.group(2)}') if m else (gene, '')
            a = re.sub(r'[_-]rep\d+$', '', a)
            for c in range(2, ws.max_column + 1):
                v = ws.cell(r, c).value
                if num(v):
                    tr, rep = hdr[c - 1].rsplit('-rep', 1)
                    out.append(dict(group=a, series=series_norm(tr), replicate=f'rep{rep}', value=float(v),
                                    sample_id=(f'{g0}|{expt}' if expt else g0), cell=f'{L(c)}{r}',
                                    notes=f'gene={g0}' + (f'; experiment={expt}' if expt else '')))
    return out


def p_qpcr_biorep_blocks(ws):
    """'Fig_6e DOFs qPCR' (3 biological replicates): per gene a row 'Mock ... RALF23' (treatment
    headers), a row 'Bio-rep1 .. Bio-rep3' (each spans 3 columns), then rows [gene, 'Col-0', values] / [None, 'fer-4',
    values]. Each value is one sample; replicate = 'bio<k>_<i>' (i = position within the bio rep)."""
    out = []
    trt, bio, gene = {}, {}, None
    for r in range(1, ws.max_row + 1):
        row = [ws.cell(r, c).value for c in range(1, ws.max_column + 1)]
        if 'Mock' in row and 'RALF23' in row:
            cur = None
            for c, v in enumerate(row, 1):
                if v in ('Mock', 'RALF23'):
                    cur = v
                trt[c] = cur if c >= 3 else None
            continue
        if any(isinstance(v, str) and v.startswith('Bio-rep') for v in row):
            cur = None
            for c, v in enumerate(row, 1):
                if isinstance(v, str) and v.startswith('Bio-rep'):
                    cur = v.replace('Bio-rep', '')
                bio[c] = cur if c >= 3 else None
            continue
        if row[1] in ('Col-0', 'fer-4'):
            gene = row[0] or gene
            pos = {}
            for c in range(3, len(row) + 1):
                if trt.get(c) is None or bio.get(c) is None:
                    continue
                key = (trt[c], bio[c]); pos[key] = pos.get(key, 0) + 1   # position counts empty cells too
                v = row[c - 1]
                if num(v):
                    out.append(dict(group=row[1], series=series_norm(trt[c]), replicate=f'bio{bio[c]}_{pos[key]}',
                                    value=float(v), sample_id=gene, cell=f'{L(c)}{r}',
                                    notes=f'gene={gene}; bio_rep={bio[c]}'))
    return out


def p_dofs_gfp(ws):
    out = []
    for c0, rep in ((1, 'pDAG1:erGFP'), (5, 'pDOF2.2:erGFP')):
        for k, s in ((1, 'Mock'), (2, 'RALF23')):
            for r in range(2, ws.max_row + 1):
                v = ws.cell(r, c0 + k).value
                if num(v):
                    out.append(dict(group=rep, series=s, replicate=str(ws.cell(r, c0).value), value=float(v),
                                    cell=f'{L(c0 + k)}{r}'))
    return out


def p_pi_puncta(ws):
    out = []
    head = None
    for c in range(2, ws.max_column + 1):
        h = ws.cell(1, c).value
        if h is not None:
            head = str(h)
        t = ws.cell(2, c).value
        if t is None or head is None:
            continue
        g, s = head.rsplit(' ', 1)
        for r in range(3, ws.max_row + 1):
            v = ws.cell(r, c).value
            if num(v) and not is_summary(ws, r):
                out.append(dict(group=norm(g), series=series_norm(s), replicate=f'{ws.cell(r, 1).value}_{t}',
                                value=float(v), sample_id=str(ws.cell(r, 1).value), cell=f'{L(c)}{r}', notes=f'tissue={t}'))
    return out


def p_luc(ws):
    """Experiment blocks: row with a date in column A and construct names; next row = 'Firefly LUC, Renilla LUC,
    Firefly/Renilla' triplets; then well rows. value = Firefly/Renilla per well (recomputed from the raw counts and
    checked). The workbook's mean / 'Normalized' (fold over GUS-GFP) / P-value rows are read into notes per construct.
    series = experiment date."""
    import datetime
    out = []
    r = 1
    while r <= ws.max_row:
        a = ws.cell(r, 1).value
        if isinstance(a, datetime.datetime):
            exp = a.strftime('%Y-%m-%d')
            names = {c: ws.cell(r, c).value for c in range(2, ws.max_column + 1, 3) if ws.cell(r, c).value}
            rr = r + 2
            wells = []
            while rr <= ws.max_row and any(num(ws.cell(rr, c).value) for c in range(2, ws.max_column + 1)) \
                    and ws.cell(rr, 3).value != 'Normalized' and not isinstance(ws.cell(rr, 1).value, datetime.datetime):
                wells.append(rr)
                rr += 1
            # wells: rows whose Firefly column is numeric; summary rows have Firefly empty
            for c, nm in names.items():
                k = 0
                summ = {}
                for w in wells:
                    f, rn, ratio = ws.cell(w, c).value, ws.cell(w, c + 1).value, ws.cell(w, c + 2).value
                    if num(f) and num(rn):
                        k += 1
                        calc = f / rn
                        note = '' if not num(ratio) else ('' if abs(calc - ratio) / ratio < 1e-6 else f'workbook ratio {ratio} != F/R {calc:.6g}')
                        if not num(ratio):
                            note = 'ratio cell empty in workbook; F/R computed'
                        out.append(dict(group=norm(nm), series=exp, replicate=f'{exp}_w{k}', value=calc,
                                        sample_id=exp, cell=f'{L(c)}{w}:{L(c + 2)}{w}', notes=note))
                    elif num(ratio):
                        summ['mean_row'] = ratio
                # normalised / p rows
                for w in range(rr, min(rr + 3, ws.max_row + 1)):
                    if ws.cell(w, c + 1).value == 'Normalized' or ws.cell(w, 3).value == 'Normalized':
                        summ['normalized'] = ws.cell(w, c + 2).value
                    v = ws.cell(w, c + 2).value
                    if isinstance(v, str) and v.startswith('P'):
                        summ['p'] = v
                if summ:
                    out.append(dict(group=norm(nm), series=exp, replicate=f'{exp}_summary', value=float('nan'),
                                    sample_id=exp, cell=f'{L(c + 2)}', notes='; '.join(f'{k}={v}' for k, v in summ.items())))
            r = rr
        else:
            r += 1
    return out


PARSERS = {k[2:] if False else k: v for k, v in dict(two_row=p_two_row, one_row_series=p_one_row_series,
           left_right_ave=p_left_right_ave, aux1_wol=p_aux1_wol, root_hair=p_root_hair, qpcr_genes=p_qpcr_genes,
           qpcr_geno_blocks=p_qpcr_geno_blocks, qpcr_biorep_blocks=p_qpcr_biorep_blocks, dofs_gfp=p_dofs_gfp, pi_puncta=p_pi_puncta, luc=p_luc).items()}


def slug(s):
    return re.sub(r'[^A-Za-z0-9]+', '_', s).strip('_')


def main(xlsx, out):
    os.makedirs(out, exist_ok=True)
    wb = openpyxl.load_workbook(xlsx, data_only=True)
    OPTIONAL = {'Supplementary Fig_3'}   # older layout only; its experiments are part of Fig. 6e
    sheets = {k: v for k, v in SHEETS.items() if k in wb.sheetnames or k not in OPTIONAL}
    assert set(wb.sheetnames) - SKIP == set(sheets), (set(wb.sheetnames) - SKIP) ^ set(sheets)
    allrows = []
    cols = ['sheet', 'panel', 'unit', 'group', 'series', 'replicate', 'value', 'sample_id', 'cell', 'notes']
    for name, (panel, unit, parser) in sheets.items():
        rows = PARSERS[parser](wb[name])
        for r in rows:
            r.update(sheet=name, panel=panel, unit=unit)
            r.setdefault('sample_id', '')
            r.setdefault('notes', '')
        with open(os.path.join(out, slug(name) + '.csv'), 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader(); w.writerows(rows)
        allrows += rows
    with open(os.path.join(out, 'ALL_long.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader(); w.writerows(allrows)
    agg = defaultdict(list)
    for r in allrows:
        if r['value'] == r['value']:  # not NaN
            agg[(r['panel'], r['group'], r['series'])].append(r['value'])
    with open(os.path.join(out, 'summary.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['panel', 'group', 'series', 'n', 'mean', 'sd'])
        for (p, g, s), v in agg.items():
            w.writerow([p, g, s, len(v), round(stt.mean(v), 5), round(stt.stdev(v), 5) if len(v) > 1 else ''])
    print(f'{len(allrows)} rows from {len(SHEETS)} sheets; {len(agg)} groups')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'source_data', 'long'))
