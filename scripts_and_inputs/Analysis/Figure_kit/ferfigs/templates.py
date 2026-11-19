"""Collaborator templates: one Excel workbook per figure (plus plain CSV versions) that asks for source data in
exactly the tidy format the kit reads, and the importer that turns a filled workbook back into kit tables.

    python build.py templates                      -> templates/FER_FigN_source_data_template.xlsx, templates/csv/
    python build.py import-template FILE.xlsx --out data/figN/
                                                   -> CSV tables named as in config/figN.yaml, validated

Workbook layout (per figure):
  README                 how to fill the sheets (one row per measurement, units, no merged cells or colours)
  REQUESTS               what is needed for each panel (docs/DATA_REQUESTS.csv)
  <panel>_values         blank table with drop-down lists for group and series (names as in the figure)
  <panel>_stats          blank statistics table (labels exactly as they should be printed + the test used)
  <panel>_current        the values the current draft is drawn from, with their status (for checking only)
  <figure>_images        the image manifest, pre-filled with the current tiles: confirm um/px, give raw file
                         names and acquisition settings, or replace files
"""
from __future__ import annotations

import csv
import io
from pathlib import Path

import pandas as pd
import yaml

from .schemas import SCHEMAS, read_table
from .style import KIT
from .sourcedata import STATUS_TEXT, _normalize_zip

README_LINES = [
    'FER manuscript: source data template',
    '',
    'Please fill one row per individual measurement (one root, seedling, plate or image), keep the column names, and do',
    'not merge cells, add colours or summary rows. The figure kit reads these sheets directly and redraws the figure,',
    'so the plotted means, error bars and n are always computed from your rows.',
    '',
    'Columns (values sheets):',
    '  group      x-axis category exactly as in the drop-down list (genotype or time point)',
    '  series     treatment or condition exactly as in the drop-down list',
    '  replicate  an identifier unique within group and series (1, 2, 3 ... or your sample id)',
    '  value      the measured number (numbers only), in the figure\'s unit shown in the example row',
    '  unit       the figure\'s unit (or leave empty); the kit never converts units, so a row in another unit',
    '             (e.g. cm in a mm chart) is refused: convert it first',
    '  sample_id  optional: your own sample or image id, so values can be traced to raw files',
    '  experiment optional: biological replicate / experiment date when an experiment was repeated',
    '  notes      optional',
    '',
    'Statistics sheets: write the labels exactly as they should be printed (letters such as a, ab; asterisks; ns),',
    'one row per bar (group + series), or one row per bracket (group/series and group2/series2). Name the test in',
    'the test column (e.g. "two-way ANOVA with Tukey\'s multiple comparisons, P < 0.05") and give exact P values when',
    'you have them. The kit draws these labels as given and computes no statistics. Statistics belong to the values',
    'they were computed from: when you send new values, send their statistics sheet too (labels made for earlier',
    'values are flagged as stale until they are re-supplied or confirmed).',
    'In experimental_unit say what n counts in that test (roots, seedlings, plates or independent experiments), so the',
    'legend can state whether n counts independent biological replicates.',
    '',
    'Images sheet: one row per image tile. Please confirm micrometres per pixel (um_per_px) from the microscope',
    'metadata, give the raw file name and acquisition settings, and send tiles without burned-in scale bars when',
    'possible (the kit draws vector scale bars). Apply any brightness/contrast change to the whole image and equally to',
    'controls (Nature image-integrity policy); full-length unprocessed gels/blots are needed as PDF for any gel panel.',
    '',
    'The *_current sheets show the values the current draft is drawn from and their status. Values marked as recovered',
    'from the figure deck were read from the vector charts and are not your original measurements: please replace them.',
    '',
    'Import (python build.py import-template): every filled sheet is checked before anything is written. If any sheet',
    'has a problem (a name not in the drop-down list, a unit other than the figure\'s, a missing column, text where a',
    'number should be), nothing is imported and every problem is listed with its sheet and row. Imported tables are',
    'recorded in data/PROVENANCE.json with this workbook\'s name and checksum.',
]


def _dv(ws, col_letter, n_rows, options, sheet_name, wb):
    """Drop-down list validation from a hidden list sheet (works for long lists and names with commas)."""
    from openpyxl.worksheet.datavalidation import DataValidation
    lists = wb['_lists'] if '_lists' in wb.sheetnames else wb.create_sheet('_lists')
    lists.sheet_state = 'hidden'
    col = lists.max_column + 1 if lists.max_row > 1 or lists['A1'].value else 1
    for i, o in enumerate(options, start=1):
        lists.cell(i, col, o)
    from openpyxl.utils import get_column_letter
    L = get_column_letter(col)
    dv = DataValidation(type='list', formula1=f"=_lists!${L}$1:${L}${len(options)}", allow_blank=True,
                        showErrorMessage=True, errorTitle='Not in the list',
                        error='Use a name from the list (the names printed in the figure).')
    ws.add_data_validation(dv)
    dv.add(f'{col_letter}2:{col_letter}{n_rows + 1}')


def _header(ws, cols):
    from openpyxl.styles import Font, PatternFill
    ws.append(cols)
    for c in ws[1]:
        c.font = Font(bold=True)
        c.fill = PatternFill('solid', fgColor='EEEEEE')
    ws.freeze_panes = 'A2'
    for i, c in enumerate(cols, start=1):
        ws.column_dimensions[ws.cell(1, i).column_letter].width = max(12, len(c) + 4)


def _requests():
    p = KIT / 'docs' / 'DATA_REQUESTS.csv'
    if not p.exists():
        return pd.DataFrame(columns=['figure', 'panel', 'request', 'format'])
    return pd.read_csv(p, dtype=str, keep_default_na=False)


def write_figure(fig, out_dir=None):
    from openpyxl import Workbook
    from openpyxl.styles import Font
    from datetime import datetime
    cfg = yaml.safe_load((KIT / 'config' / f'{fig}.yaml').read_text())
    label = cfg.get('label', fig)
    n = label.replace('Fig. ', '').replace(' ', '')
    out_dir = Path(out_dir or KIT / 'templates')
    (out_dir / 'csv').mkdir(parents=True, exist_ok=True)
    wb = Workbook()
    ws = wb.active
    ws.title = 'README'
    for line in README_LINES:
        ws.append([line.replace('FER manuscript', f'FER manuscript, {label}') if line.startswith('FER') else line])
    ws['A1'].font = Font(bold=True)
    ws.column_dimensions['A'].width = 120
    rq = _requests()
    rq = rq[rq['figure'].str.lower() == fig.lower()]
    wr = wb.create_sheet('REQUESTS')
    _header(wr, ['panel', 'what we need', 'format'])
    for _, r in rq.iterrows():
        wr.append([r['panel'], r['request'], r['format']])
    wr.column_dimensions['B'].width = 100
    wr.column_dimensions['C'].width = 40
    manifests = set()
    custom = {}                                              # panel -> [{name, columns, file?, note?}]
    for panel, ts in (cfg.get('templates') or {}).items():   # figure-level block (custom tables)
        custom.setdefault(panel, []).extend(ts or [])
    for panel, pc in (cfg.get('panels') or {}).items():
        if pc.get('manifest'):
            manifests.add(pc['manifest'])
        if pc.get('type') != 'dotbar':
            custom.setdefault(panel, []).extend(pc.get('templates', []) or [])
            continue
        unit = pc.get('unit', '')
        cols = list(SCHEMAS['measurements']['required']) + ['unit', 'sample_id', 'experiment', 'notes']
        wv = wb.create_sheet(f'Fig{n}{panel}_values'[:31])
        _header(wv, cols)
        wv.append(['(example, delete)', '', 'R1', '', unit, '', '', ''])
        _dv(wv, 'A', 2000, list(pc['groups']), wv.title, wb)
        _dv(wv, 'B', 2000, list(pc['series']), wv.title, wb)
        with open(out_dir / 'csv' / f'Fig{n}{panel}_values.csv', 'w', newline='') as f:
            csv.writer(f, lineterminator='\n').writerow(cols)
        scols = list(SCHEMAS['annotations']['required']) + list(SCHEMAS['annotations']['optional'])
        wsS = wb.create_sheet(f'Fig{n}{panel}_stats'[:31])
        _header(wsS, scols)
        _dv(wsS, 'A', 200, list(pc['groups']), wsS.title, wb)
        _dv(wsS, 'B', 200, list(pc['series']), wsS.title, wb)
        _dv(wsS, 'D', 200, ['letters', 'stars', 'ns', 'text'], wsS.title, wb)
        with open(out_dir / 'csv' / f'Fig{n}{panel}_stats.csv', 'w', newline='') as f:
            csv.writer(f, lineterminator='\n').writerow(scols)
        # current values (for checking only)
        data_rel = pc.get('data') or pc.get('data_file')
        wc = wb.create_sheet(f'Fig{n}{panel}_current'[:31])
        if data_rel and (KIT / data_rel).exists():
            cur = pd.read_csv(KIT / data_rel, dtype=str, keep_default_na=False)
            wc.append([f'{label}{panel}: values the current draft is drawn from. ' +
                       STATUS_TEXT.get(pc.get('data_status', ''), pc.get('data_status', ''))])
            wc.append(list(cur.columns))
            for row in cur.itertuples(index=False):
                wc.append(list(row))
        else:
            wc.append([f'{label}{panel}: no values yet; the panel is drawn as a placeholder until this data arrives.'])
        wc['A1'].font = Font(bold=True, color='9E0029')
    for panel in sorted(custom):
        for t in custom[panel]:
            wsx = wb.create_sheet(f'Fig{n}{panel}_{t["name"]}'[:31])
            _header(wsx, t['columns'])
            if t.get('note'):
                wsx.cell(row=1, column=len(t['columns']) + 2, value=t['note']).font = Font(italic=True)
            if t.get('prefill') and t.get('file') and (KIT / t['file']).exists():     # current values, to edit
                cur = pd.read_csv(KIT / t['file'], dtype=str, keep_default_na=False)
                for row in cur.reindex(columns=t['columns'], fill_value='').itertuples(index=False):
                    wsx.append(list(row))
            with open(out_dir / 'csv' / f'Fig{n}{panel}_{t["name"]}.csv', 'w', newline='') as f:
                csv.writer(f, lineterminator='\n').writerow(t['columns'])
    for m in sorted(manifests):
        man = pd.read_csv(KIT / m, dtype=str, keep_default_na=False)
        wi = wb.create_sheet(f'Fig{n}_images'[:31])
        cols = list(SCHEMAS['image_manifest']['required']) + list(SCHEMAS['image_manifest']['optional'])
        _header(wi, cols)
        for _, r in man.iterrows():
            wi.append([r.get(c, '') for c in cols])
        man.reindex(columns=cols).to_csv(out_dir / 'csv' / f'Fig{n}_images.csv', index=False, lineterminator='\n')
    wb.properties.creator = 'ferfigs'
    wb.properties.created = datetime(2026, 1, 1)
    wb.properties.modified = datetime(2026, 1, 1)
    p = out_dir / f'FER_Fig{n}_source_data_template.xlsx'
    wb.save(p)
    _normalize_zip(p)
    return p


def write_all(out_dir=None):
    from .context import figures
    return [write_figure(f, out_dir) for f in figures()]


def _sheet_rows(msg, name, rows):
    """A validation message about the staged CSV, with its line numbers translated to the workbook's row numbers."""
    import re

    def sub(m):
        word = {'line ': 'row ', 'line(s) ': 'row(s) ', 'lines ': 'rows '}[m.group(1)]
        return word + ', '.join(str(rows.get(int(n), n)) for n in re.findall(r'\d+', m.group(2)))
    msg = re.sub(r'(line\(s\) |lines |line )(\d+(?:, \d+)*)', sub, msg)
    return 'sheet' + msg[len(name):] if msg.startswith(name) else msg


def _same_table(p, df) -> bool:
    """True when the CSV at p holds the same cells as df (columns missing on one side count as empty)."""
    if not Path(p).exists():
        return False
    cur = pd.read_csv(p, dtype=str, keep_default_na=False)
    cols = list(dict.fromkeys(list(cur.columns) + list(df.columns)))
    a = cur.reindex(columns=cols, fill_value='').astype(str).apply(lambda c: c.str.strip()).reset_index(drop=True)
    b = df.reindex(columns=cols, fill_value='').astype(str).apply(lambda c: c.str.strip()).reset_index(drop=True)
    return a.shape == b.shape and a.equals(b)


def _table_bytes(df, name) -> bytes:
    """CSV bytes as the kit writes them; a .gz target is gzip-compressed with a fixed header (byte-reproducible)."""
    import gzip
    b = df.to_csv(index=False, lineterminator='\n').encode('utf-8')
    return gzip.compress(b, mtime=0) if str(name).endswith('.gz') else b


def _inside(p, kit) -> bool:
    try:
        Path(p).resolve().relative_to(Path(kit).resolve())
        return True
    except ValueError:
        return False


def _dry_run(kit, out_dir, staged):
    """Build the affected figures once in a scratch copy of the kit with the staged tables in place.
    Returns (problems, notes). Nothing in the kit is touched."""
    import json
    import os
    import shutil
    import subprocess
    import sys
    import tempfile
    rel_out = Path(out_dir).resolve().relative_to(Path(kit).resolve())
    figs = sorted({s['fig'] for s in staged})
    problems, notes = [], []
    with tempfile.TemporaryDirectory() as td:
        tk = Path(td) / 'kit'
        shutil.copytree(kit, tk, ignore=shutil.ignore_patterns('output', '__pycache__', '.pytest_cache', '*.pyc',
                                                               'tests', 'templates', '.tmp', 'scratch'))
        for s in staged:
            dst = tk / rel_out / s['name']
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_bytes(s['bytes'])
        r = subprocess.run([sys.executable, 'build.py', *figs], cwd=tk, capture_output=True, text=True)
        rep_p = tk / 'output' / 'report.json'
        rep = json.loads(rep_p.read_text()) if rep_p.exists() else {}
        # messages name kit files relative to the kit (data/fig6/luc.csv), not the scratch copy
        prefixes = sorted({str(tk), os.path.abspath(tk), str(tk.resolve()), os.path.relpath(tk),
                           './' + os.path.relpath(tk)}, key=len, reverse=True)

        def _kit_rel(msg):
            for pre in prefixes:
                msg = msg.replace(pre + os.sep, '').replace(pre, '.')
            return msg
        for fig, kind, probs in rep.get('failures', []):
            problems += [f'dry-run build of {fig} with the new tables: {_kit_rel(str(q))}' for q in probs]
        if r.returncode != 0 and not rep.get('failures'):
            problems.append(f'dry-run build failed: {_kit_rel((r.stderr or r.stdout).strip()[-600:])}')
        reads = {}
        for res in rep.get('results', []):
            reads[res['figure'].lower()] = {i['path'] for i in res.get('inputs', [])}
        if not problems:
            for s in staged:
                rel = (rel_out / s['name']).as_posix()
                if rel in reads.get(s['fig'], set()):
                    continue
                if s.get('used_when'):
                    notes.append(f"{s['sheet']}: imported to {rel}; it is drawn only when {s['used_when']}.")
                else:
                    problems.append(f"{s['sheet']}: {s['fig'].replace('fig', 'Fig. ')} does not read {rel}, so the "
                                    f"values would be ignored. Check the target file (config/{s['fig']}.yaml).")
    return problems, notes


def import_workbook(path, out_dir, kit=None, register=True, dry_run=True):
    """Filled template -> CSV tables named exactly as the config of that figure expects (e.g. the Fig4e_values sheet
    becomes the file named in config/fig4.yaml panels: e: data:).

    Validate first, then write. Every filled sheet is checked before anything is written: group and series names
    against the figure's lists, the unit column against the panel's unit, the columns of custom tables, and
    statistics marked as proposed by the helper are refused. If anything is wrong a ValidationError lists every
    problem and no file is changed. Then each table is written through a temporary file and an atomic rename, and
    tables written inside the kit are registered in data/PROVENANCE.json (status author_supplied, workbook name and
    checksum; statistics are bound to the values they were supplied with). Sheets left blank keep the current
    tables; values imported without a statistics sheet leave the old labels bound to the old values, so the build
    warns that they are stale.

    Before anything is written, the affected figures are also built once in a scratch copy of the kit with the new
    tables in place (dry run): the import is refused if that build fails (e.g. a custom table the drawing code cannot
    use) or if the figure does not read a file the workbook would write (a table that would be silently ignored),
    unless the template entry says when the file is used (`used_when:`). A statistics sheet that is filled in is
    bound to the values on disk after the import even when its cells did not change (re-supplying unchanged labels
    confirms them); a blank statistics sheet never does.

    Returns dict(written=[paths], unchanged=[paths], ignored=[sheet names], registered=[kit-relative paths],
    notes=[text])."""
    import os
    import re
    import tempfile
    from . import provenance as P
    from .schemas import ValidationError
    kit = Path(kit or KIT)
    path = Path(path)
    xl = pd.ExcelFile(path)
    out_dir = Path(out_dir)
    staged, problems, ignored, notes = [], [], [], []
    for sh in xl.sheet_names:
        m = re.fullmatch(r'Fig(\d+)([a-z]?)_(\w+)', sh)
        if not m:
            if sh not in ('README', 'REQUESTS') and not sh.startswith('_'):    # _lists holds the drop-down lists
                ignored.append(sh)
            continue
        if sh.endswith('_current'):
            continue
        fig, panel, what = f'fig{m.group(1)}', m.group(2), m.group(3)
        cfgp = kit / 'config' / f'{fig}.yaml'
        if not cfgp.exists():
            problems.append(f'{sh}: this kit has no config/{fig}.yaml')
            continue
        cfg = yaml.safe_load(cfgp.read_text())
        df = pd.read_excel(xl, sh, dtype=str).fillna('')
        df = df.loc[:, [c for c in df.columns if not str(c).startswith('Unnamed')]]
        if len(df.columns):
            df = df[~df.iloc[:, 0].astype(str).str.startswith('(example')]
        df = df[(df != '').any(axis=1)]
        if df.empty:                 # sheets left blank keep the kit's current tables
            continue
        if what not in ('values', 'stats', 'images'):            # custom table declared under templates:
            ts = list((cfg.get('templates') or {}).get(panel, []) or []) + \
                 list(((cfg.get('panels') or {}).get(panel) or {}).get('templates', []) or [])
            t = next((t for t in ts if t['name'] == what), None)
            if t is None:
                ignored.append(sh)
                continue
            missing = [c for c in t['columns'] if c not in df.columns]
            if missing:
                problems.append(f'{sh}: missing column(s) {missing}; keep the template header.')
                continue
            staged.append(dict(sheet=sh, df=df[t['columns']], name=Path(t.get('file') or f'{fig}{panel}_{what}.csv').name,
                               kind='custom', levels=None, unit=None, data_name=None, fig=fig,
                               used_when=t.get('used_when'), role=t.get('role', 'custom'),
                               rows={k + 2: int(ix) + 2 for k, ix in enumerate(df.index)}))
            continue
        kind = {'values': 'measurements', 'stats': 'annotations', 'images': 'image_manifest'}[what]
        keep = list(SCHEMAS[kind]['required']) + [c for c in SCHEMAS[kind]['optional'] if c in df.columns]
        df = df[[c for c in keep if c in df.columns]]
        levels = unit = data_name = None
        if what == 'images':
            target = cfg.get('manifest') or next((pc['manifest'] for pc in (cfg.get('panels') or {}).values()
                                                  if pc.get('manifest')), None)
            if target is None:
                problems.append(f'{sh}: config/{fig}.yaml names no image manifest')
                continue
        else:
            pc = (cfg.get('panels') or {}).get(panel) or {}
            data_target = pc.get('data') or pc.get('data_file')
            if pc.get('type') != 'dotbar' or not data_target:
                problems.append(f'{sh}: config/{fig}.yaml has no dot-bar panel {panel!r} with a data table')
                continue
            target = data_target if what == 'values' else \
                (pc.get('stats') or pc.get('stats_file') or Path(data_target).stem + '_stats.csv')
            levels = {'group': pc['groups'], 'series': pc['series']}
            if what == 'stats':
                levels.update({'group2': pc['groups'], 'series2': pc['series']})
            unit = pc.get('unit') if what == 'values' else None
            data_name = Path(data_target).name
        staged.append(dict(sheet=sh, df=df, name=Path(target).name, kind=kind, levels=levels, unit=unit,
                           data_name=data_name, fig=fig, used_when=None,
                           role={'measurements': 'measurements', 'annotations': 'statistics',
                                 'image_manifest': 'image_manifest'}[kind],
                           rows={k + 2: int(ix) + 2 for k, ix in enumerate(df.index)}))
    names = [s['name'] for s in staged]
    for n in sorted({n for n in names if names.count(n) > 1}):
        problems.append(f"sheets {[s['sheet'] for s in staged if s['name'] == n]} would all write {n}")
    with tempfile.TemporaryDirectory() as td:          # validate every staged table before writing any
        for s in staged:
            tp = Path(td) / s['name']
            s['bytes'] = _table_bytes(s['df'], s['name'])
            tp.write_bytes(s['bytes'])
            if s['kind'] == 'custom':
                continue
            try:
                read_table(tp, s['kind'], s['levels'], unit=s['unit'], allow_proposed=False)
            except ValidationError as e:
                problems += [f"{s['sheet']}: {_sheet_rows(q, s['name'], s['rows'])}" for q in e.problems]
    if problems:
        raise ValidationError(problems + ['Nothing was written. Correct the workbook and import it again.'])
    changed = [s for s in staged if not _same_table(out_dir / s['name'], s['df'])]
    if changed and dry_run and _inside(out_dir, kit):
        probs, dnotes = _dry_run(kit, out_dir, changed)
        if probs:
            raise ValidationError(probs + ['Nothing was written. Correct the workbook and import it again.'])
        notes += dnotes
    out_dir.mkdir(parents=True, exist_ok=True)
    written, unchanged = [], []
    for s in staged:
        p = out_dir / s['name']
        if _same_table(p, s['df']):            # e.g. the pre-filled images sheet returned untouched
            unchanged.append(p)
            continue
        tmp = p.with_name(p.name + '.importing')
        tmp.write_bytes(s['bytes'])
        os.replace(tmp, p)                             # atomic: a table is either the old or the new version
        written.append(p)

    def _rel(p):
        try:
            return Path(p).resolve().relative_to(kit.resolve()).as_posix()
        except ValueError:
            return None
    registered = []
    if register:
        reg = P.load(kit)
        wsha = P.sha256(path)[:12]
        for s in sorted(staged, key=lambda s: s['kind'] == 'annotations'):     # values before their statistics
            p = out_dir / s['name']
            r = _rel(p)
            if r is None:
                continue
            src = f"import-template {path.name} (sha256 {wsha}), sheet {s['sheet']}"
            d = out_dir / s['data_name'] if s['data_name'] else None
            if s['kind'] == 'annotations' and d is not None and d.exists() and _rel(d):
                e = reg.get(r) or {}
                if p in unchanged and e.get('for_data') == _rel(d) and e.get('sha256') == P.sha256(p) \
                        and e.get('for_data_sha256') == P.sha256(d):
                    continue                           # already bound to exactly these values
                P.register(r, kit, source=src, for_data=_rel(d), reg=reg)   # filled-in sheet = confirmed labels
            else:
                if p in unchanged:
                    continue                           # same content as the current table: keep its record
                P.register(r, kit, status='author_supplied', source=src, reg=reg, role=s['role'])
            registered.append(r)
        if registered:
            P.save(reg, kit)
    return dict(written=written, unchanged=unchanged, ignored=ignored, registered=registered, notes=notes)
