"""Source Data workbooks in the Nature Portfolio format: one Excel file per figure, one sheet per panel, the
linked figure noted in the file. Built from exactly the tables the figures were drawn from.

Each chart sheet holds the individual values (tidy), the per-bar summary actually drawn (n, mean, SD, SEM) and
the statistics labels with the test named in the stats table. Panels whose values were recovered from figure
drawings say so in the first rows of the sheet; replace them with original measurements before submission.
"""
from __future__ import annotations

import io
import re
import zipfile
from pathlib import Path

import pandas as pd

STATUS_TEXT = {
    'original': 'Original measurements supplied by the authors.',
    'author_supplied': 'Values supplied by the authors through the collaborator template (workbook and checksum in '
                       'data/PROVENANCE.json).',
    'unregistered': 'VALUES OF UNKNOWN PROVENANCE: the table was added or changed without being registered '
                    '(python build.py import-template or register). Check before use.',
    'layout': 'Drawing geometry or label tables (positions, lane annotations, colours); not measured values.',
    'recovered_from_deck_vector': 'VALUES RECOVERED FROM THE VECTOR CHART IN THE ORIGINAL FIGURE DECK (marker positions, '
                                  'drawing precision). Not the original measurement files: replace before submission.',
    'digitized_from_raster': 'APPROXIMATE VALUES DIGITIZED FROM A RASTER IMAGE. Not for submission.',
    'computed_from_upstream': 'Computed by the kit from the pinned upstream analysis outputs (see docs/PROVENANCE.md).',
    'upstream_table': 'Table produced by the upstream analysis code (see docs/PROVENANCE.md).',
    'placeholder': 'PLACEHOLDER: no source data yet.',
    'mockup': 'MOCK-UP FOR LAYOUT ONLY: bars, error bars and labels traced from the original chart; individual '
              'points simulated. Not data; never for submission.',
    'original_images_from_deck': 'Native image pixels from the original figure deck (cropped and oriented as placed; '
                                 'no pixel values changed).',
    'schematic': 'Line art drawn by code from words in the figure config (no data).',
    'displayed_values_from_deck': 'VALUES AS PRINTED IN THE ORIGINAL FIGURE DECK (rounded display values, not the '
                                  'measurements): replace with the measurement table before submission.',
}


def _normalize_zip(path):
    """Fixed timestamps inside the .xlsx zip so identical content gives identical bytes."""
    src = Path(path).read_bytes()
    zin = zipfile.ZipFile(io.BytesIO(src))
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as zout:
        for item in sorted(zin.infolist(), key=lambda i: i.filename):
            zi = zipfile.ZipInfo(item.filename, date_time=(2026, 1, 1, 0, 0, 0))
            zi.compress_type = zipfile.ZIP_DEFLATED
            zi.external_attr = 0o644 << 16
            data = zin.read(item.filename)
            if item.filename == 'docProps/core.xml':     # openpyxl writes the save time here
                data = re.sub(rb'(<dcterms:(?:created|modified)[^>]*>)[^<]*(</dcterms:(?:created|modified)>)',
                              rb'\g<1>2026-01-01T00:00:00Z\g<2>', data)
            zout.writestr(zi, data)
    Path(path).write_bytes(buf.getvalue())


def write(fig_label, records, out_path, kit):
    """records: {panel: record from ctx.record}. A record may carry 'data' (measurements CSV), 'stats',
    'source_tables' ({sheet_suffix: csv path}), 'source_rows' ({sheet_suffix: {'columns': [...], 'rows': [[...]]}},
    tables computed during the build, e.g. the ratios actually drawn) and 'data_status'."""
    from openpyxl import Workbook
    from openpyxl.styles import Font
    from datetime import datetime
    wb = Workbook()
    ws0 = wb.active
    ws0.title = 'README'
    ws0.append([f'Source Data for {fig_label}'])
    ws0['A1'].font = Font(bold=True)
    ws0.append(['One sheet per panel. Values are exactly those drawn by the figure kit (ferfigs).'])
    ws0.append([])
    ws0.append(['Panel', 'Sheet', 'Data status', 'Table(s)'])
    for c in ws0[4]:
        c.font = Font(bold=True)
    for panel, rec in sorted(records.items()):
        sheet = f'{fig_label.replace(" ", "").replace(".", "")}{panel}'[:31]
        status = rec.get('data_status', 'unspecified')
        tables = [rec[k] for k in ('data', 'stats') if rec.get(k)] + list((rec.get('source_tables') or {}).values())
        tables += [f'{n} (computed by the build)' for n in (rec.get('source_rows') or {})]
        if not tables and rec.get('type') != 'dotbar':
            continue
        ws0.append([panel, sheet, STATUS_TEXT.get(status, status), '; '.join(map(str, tables))])
        ws = wb.create_sheet(sheet)
        ws.append([f'{fig_label}{panel}'])
        ws['A1'].font = Font(bold=True)
        ws.append([STATUS_TEXT.get(status, status)])
        if rec.get('inputs'):                    # every input file of the panel with its own role and status
            ws.append(['Input file', 'role', 'status (data/PROVENANCE.json)'])
            for c in ws[ws.max_row]:
                c.font = Font(italic=True)
            for path, v in rec['inputs'].items():
                ws.append([path, v.get('role', ''), v.get('status', '')])
        ws.append([])
        if rec.get('data'):
            d = pd.read_csv(Path(kit) / rec['data'], dtype=str, keep_default_na=False)
            ws.append(['Individual values'])
            ws.cell(ws.max_row, 1).font = Font(bold=True)
            ws.append(list(d.columns))
            for row in d.itertuples(index=False):
                ws.append([_num(v) for v in row])
            ws.append([])
        if rec.get('bars'):
            ws.append(['Summary drawn (bar = mean; error bar = ' + rec.get('error', 'sd').upper() + ')'])
            ws.cell(ws.max_row, 1).font = Font(bold=True)
            ws.append(['group', 'series', 'n', 'mean', 'SD', 'error bar'])
            for b in rec['bars']:
                ws.append([b['group'], b['series'], b['n'], b['mean'], b['sd'], b['err']])
            ws.append([])
        if rec.get('stats') and (Path(kit) / rec['stats']).exists():   # absent in the public copy (no author labels)
            s = pd.read_csv(Path(kit) / rec['stats'], dtype=str, keep_default_na=False)
            ws.append(['Statistics labels as printed'])
            ws.cell(ws.max_row, 1).font = Font(bold=True)
            ws.append(list(s.columns))
            for row in s.itertuples(index=False):
                ws.append(list(row))
            ws.append([])
        for name, t in (rec.get('source_rows') or {}).items():
            ws.append([name])
            ws.cell(ws.max_row, 1).font = Font(bold=True)
            ws.append(list(t['columns']))
            for row in t['rows']:
                ws.append([_num(v) for v in row])
            ws.append([])
        for name, p in (rec.get('source_tables') or {}).items():
            if not (Path(kit) / p).exists():      # label tables are absent in the public copy
                continue
            t = pd.read_csv(Path(kit) / p, dtype=str, keep_default_na=False)
            ws.append([name])
            ws.cell(ws.max_row, 1).font = Font(bold=True)
            ws.append(list(t.columns))
            for row in t.itertuples(index=False):
                ws.append([_num(v) for v in row])
            ws.append([])
    wb.properties.creator = 'ferfigs'
    wb.properties.created = datetime(2026, 1, 1)
    wb.properties.modified = datetime(2026, 1, 1)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)
    _normalize_zip(out_path)
    return out_path


def _num(v):
    try:
        f = float(v)
        return int(f) if f.is_integer() and '.' not in str(v) else f
    except (TypeError, ValueError):
        return v
