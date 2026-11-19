"""Command line:  python build.py <command> [options]

  python build.py                       build every figure, check it against the journal profile, write
                                        output/REPORT.md and the Source Data workbooks
  python build.py fig4 fig6             build only these figures
  python build.py validate              check every table named in the configs (no drawing)
  python build.py check                 re-check the PDFs already in output/
  python build.py stats-propose FILE --test twoway_tukey|tukey_all|tukey_within|vs_control
                                        write FILE_stats_PROPOSED.csv for author review (never drawn)
  python build.py templates             regenerate the collaborator templates in templates/
  python build.py import-template FILE.xlsx --out data/figN/
                                        convert a filled template workbook back into the kit's CSV tables
                                        (validates every sheet first; writes nothing if any sheet has a problem;
                                        registers what it writes in data/PROVENANCE.json)
  python build.py register TABLE --status STATUS --source "who/what"
  python build.py register STATS --for DATA --source "who/how"
                                        record a table edited by hand (or confirm statistics for these values)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import traceback
import warnings
from pathlib import Path

from .style import KIT


def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def environment():
    """Versions and font files that determine the bytes of the outputs. Rebuilding with identical entries gives
    byte-identical files; other versions or font files can change bytes (compare with the review PDF and `check`)."""
    import platform
    import matplotlib
    import matplotlib.font_manager as fm
    import numpy
    import openpyxl
    import pandas
    import PIL
    import yaml
    from .style import resolve_font
    try:
        import fitz
        pymupdf = str(fitz.VersionBind)
    except Exception:                 # noqa: BLE001
        pymupdf = 'not installed'
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fam = resolve_font()
    files = sorted({f.fname for f in fm.fontManager.ttflist if f.name == fam})
    return dict(python=platform.python_version(), platform=f'{platform.system()} {platform.machine()}',
                matplotlib=matplotlib.__version__, numpy=numpy.__version__, pandas=pandas.__version__,
                pillow=PIL.__version__, pyyaml=yaml.__version__, openpyxl=openpyxl.__version__, pymupdf=pymupdf,
                font_family=fam, font_files={Path(f).name: _sha(f)[:16] for f in files})


def _scale_note(r):
    """Smallest text and the smallest print scale at which it stays at the journal's 5-pt floor (scaled placements,
    e.g. a full-page figure shrunk to share its page with the legend, also shrink the 8-pt panel letters)."""
    sizes = [float(k) for k in (r.get('info', {}).get('text_sizes_pt') or {})]
    if not sizes:
        return 'smallest text: none'
    from .style import load_profile
    floor = float(load_profile()['text_pt']['min'])
    lo = min(sizes)
    smin = min(1.0, floor / lo)
    w = r.get('info', {}).get('size_mm', [180, 0])[0]
    return (f"smallest text {lo:g} pt: stays at or above {floor:g} pt when printed at {smin:.0%} or more "
            f"({w * smin:.0f} mm wide or wider); panel letters are 8 pt only at 100%")


def _panels_with_bad_images(ctx, bad):
    """A panel whose image (listed for it in its manifest) is unregistered or changed is itself 'unregistered'."""
    if not bad:
        return
    import pandas as pd
    for panel, rec in ctx.panel_records.items():
        man = rec.get('manifest')
        if not man or not ctx.path(man).exists():
            continue
        m = pd.read_csv(ctx.path(man), dtype=str, keep_default_na=False)
        if not {'panel', 'file'} <= set(m.columns):
            continue
        base = Path(man).parent
        files = {(base / f).as_posix() for f in m.loc[m['panel'] == panel, 'file']}
        if files & bad:
            rec['data_status'] = 'unregistered'
            for f in sorted(files & bad):
                rec.setdefault('inputs', {})[f] = dict(role='image', status='unregistered')


def _inputs_note(inputs):
    from collections import Counter
    if not inputs:
        return 'none recorded'
    roles = Counter(i['role'] for i in inputs)
    bad = [i['path'] for i in inputs if i['status'] == 'unregistered']
    return (f"{len(inputs)} files under data/ ({', '.join(f'{k} {v}' for k, v in sorted(roles.items()))}); "
            f"unregistered or changed: {bad or 'none'}")


def build(figs=None, out=None, check=True, source_data=True):
    from .context import Context, load_module, figures
    from .schemas import ValidationError
    from . import check as C, sourcedata as SD, provenance as PV
    out = Path(out or KIT / 'output')
    figs = figs or figures()
    results, failures = [], []
    for fig in figs:
        t0 = time.time()
        try:
            with warnings.catch_warnings(record=True) as wlist, PV.recording() as reads:
                warnings.simplefilter('always')
                ctx = Context(fig)
                mod = load_module(fig)
                page = mod.build(ctx)
                paths, rec = page.save(out / 'figures', stem=ctx.cfg.get('output_stem', fig))
            for w in wlist:
                if 'font cache' not in str(w.message):
                    ctx.warnings.append(str(w.message))
        except ValidationError as e:
            failures.append((fig, 'data problem', e.problems))
            print(f'[{fig}] DATA PROBLEM:\n  - ' + '\n  - '.join(e.problems), file=sys.stderr)
            continue
        except Exception as e:      # noqa: BLE001
            failures.append((fig, 'error', [repr(e)]))
            traceback.print_exc()
            continue
        r = C.check_pdf(paths['pdf'], rec) if check and 'pdf' in paths else dict(figure=fig, errors=[], warnings=[],
                                                                                 info={})
        r = C.apply_accepted(r, ctx.cfg.get('accepted'))
        r['inputs'] = []                          # every file under data/ the build opened, with its registry status
        for rel in PV.data_inputs(reads):
            st, role = ctx.input_status(rel) if ctx.registry() else ('unchecked', 'unknown')
            r['inputs'].append(dict(path=rel, role=role, status=st))
        _panels_with_bad_images(ctx, {i['path'] for i in r['inputs'] if i['status'] == 'unregistered'})
        r['build_warnings'] = ctx.warnings
        r['data_warnings'] = ctx.data_warnings
        secs = round(time.time() - t0, 1)         # printed only, so report.json is identical between builds
        r['outputs'] = {k: dict(path=v.relative_to(out).as_posix(), sha256=_sha(v)) for k, v in paths.items()}
        r['data_status'] = {p: v.get('data_status') for p, v in ctx.panel_records.items()}
        results.append(r)
        if source_data:
            label = ctx.cfg.get('label', fig)
            SD.write(label, ctx.panel_records, out / 'source_data' / f"Source_Data_{label.replace(' ', '_').replace('.', '')}.xlsx",
                     KIT)
        n_acc, n_dw = len(r.get('accepted', [])), len(r['data_warnings'])
        print(f"[{fig}] built in {secs} s; {len(r['errors'])} journal error(s), {len(r['warnings'])} warning(s)"
              + (f", {n_acc} accepted on purpose (reasons in output/REPORT.md)" if n_acc else '')
              + f"; {n_dw} data warning(s)")
        for w in r['data_warnings']:
            print(f'  data warning: {w}')
    rep = out / 'REPORT.md'
    lines = ['# Build report', '', C.report_markdown(results), '']
    for r in results:
        lines += [f"## {r['figure']}", '', f"- outputs: {json.dumps(r['outputs'])}",
                  f"- data status by panel: {json.dumps(r['data_status'])}",
                  f"- text sizes (pt -> characters): {json.dumps(r['info'].get('text_sizes_pt', {}))}",
                  f"- stroke widths (pt -> count): {json.dumps(r['info'].get('stroke_widths_pt', {}))}",
                  f"- image ppi at print size: {r['info'].get('image_ppi', [])}",
                  f"- {_scale_note(r)}",
                  f"- inputs read: {_inputs_note(r['inputs'])}",
                  f"- build warnings: {r['build_warnings'] or 'none'}",
                  f"- data warnings (provenance): {r['data_warnings'] or 'none'}",
                  f"- accepted on purpose (author decisions, not journal waivers): "
                  f"{json.dumps(r.get('accepted', [])) if r.get('accepted') else 'none'}",
                  f"- accepted rules that matched nothing: {r.get('accepted_unused') or 'none'}", '']
    for fig, kind, probs in failures:
        lines += [f'## {fig}: {kind}', ''] + [f'- {p}' for p in probs] + ['']
    env = environment()
    lines += ['## Build environment', '',
              'Rebuilding with these versions and font files gives byte-identical outputs; other versions or font '
              'files can change the bytes (and text widths slightly). Compare against the review PDF and `check`.', '',
              f'- {json.dumps(env)}', '']
    out.mkdir(parents=True, exist_ok=True)
    rep.write_text('\n'.join(lines))
    (out / 'report.json').write_text(json.dumps(dict(results=results, failures=failures, environment=env), indent=1,
                                                default=str))
    return results, failures


def _panel_levels(pc, key):
    if pc.get('groups') and pc.get('series') and key in ('data', 'stats', 'data_file', 'stats_file'):
        lv = {'group': pc['groups'], 'series': pc['series']}
        if key.startswith('stats'):
            lv.update({'group2': pc['groups'], 'series2': pc['series']})
        return lv
    return None


def validate():
    """Check every table named in the configs (schemas, names, units, proposed statistics) without drawing, and list
    provenance warnings (unregistered tables, statistics made for other values)."""
    import yaml
    from .schemas import read_table, ValidationError
    from . import provenance as PV
    probs, n, dwarn = [], 0, []
    reg = PV.load(KIT)
    for cfg in sorted((KIT / 'config').glob('fig*.yaml')):
        c = yaml.safe_load(cfg.read_text())
        for panel, pc in (c.get('panels') or {}).items():
            for key, default_kind in (('data', 'measurements'), ('stats', 'annotations'), ('manifest', 'image_manifest'),
                                      ('data_file', 'measurements'), ('stats_file', 'annotations')):
                rel = pc.get(key)
                if not rel:
                    continue
                if (pc.get('optional') or key.endswith('_file')) and not (KIT / rel).exists():
                    continue            # optional table not supplied yet (e.g. awaiting collaborator data)
                kind = pc.get(f'{key}_kind', default_kind)
                if kind == 'custom':
                    continue
                unit = pc.get('unit') if kind == 'measurements' else None
                try:
                    read_table(KIT / rel, kind, _panel_levels(pc, key), unit=unit,
                               allow_proposed=bool(pc.get('allow_proposed_stats')))
                    n += 1
                except ValidationError as e:
                    probs += [f'{cfg.stem} {panel}: {p}' for p in e.problems]
            d, st = pc.get('data') or pc.get('data_file'), pc.get('stats') or pc.get('stats_file')
            if pc.get('type') == 'dotbar' and d and (KIT / d).exists():
                _, w = PV.data_status(d, pc.get('data_status', 'unspecified'), KIT, reg)
                w2 = PV.stats_check(st, d, KIT, reg) if st and (KIT / st).exists() else None
                dwarn += [f'{cfg.stem} {panel}: {x}' for x in (w, w2) if x]
    from .context import Context, load_module, figures
    for fig in figures():                     # figure-specific tables (e.g. Fig. 1 tables: block)
        mod = load_module(fig)
        if hasattr(mod, 'validate'):
            try:
                n += int(mod.validate(Context(fig)) or 0)
            except ValidationError as e:
                probs += [f'{fig}: {p}' for p in e.problems]
    print(f'{n} tables valid' if not probs else '\n'.join(probs))
    for w in dwarn:
        print(f'data warning: {w}')
    return probs


def register_cmd(rest):
    """python build.py register TABLE --status S --source TEXT | register STATS --for DATA --source TEXT"""
    import yaml
    from .schemas import read_table, ValidationError
    from .sourcedata import STATUS_TEXT
    from . import provenance as PV
    ap = argparse.ArgumentParser(prog='build.py register')
    ap.add_argument('table')
    ap.add_argument('--status', choices=sorted(k for k in STATUS_TEXT if k != 'unregistered'))
    ap.add_argument('--for', dest='for_data')
    ap.add_argument('--source', required=True, help='who supplied or computed the table, and how')
    a = ap.parse_args(rest)
    rel = Path(a.table).as_posix()
    if not (KIT / rel).exists():
        print(f'{rel}: file not found (give the path relative to the kit folder)')
        return 1
    if bool(a.status) == bool(a.for_data):
        print('give --status for a data table, or --for DATA for a statistics table (one of the two)')
        return 1
    if a.for_data and not (KIT / a.for_data).exists():
        print(f'{a.for_data}: file not found')
        return 1
    for cfg in sorted((KIT / 'config').glob('fig*.yaml')):     # validate like the build would before recording
        c = yaml.safe_load(cfg.read_text())
        for panel, pc in (c.get('panels') or {}).items():
            for key, kind in (('data', 'measurements'), ('data_file', 'measurements'), ('stats', 'annotations'),
                              ('stats_file', 'annotations')):
                if pc.get(key) == rel and pc.get(f'{key}_kind', kind) != 'custom':
                    try:
                        read_table(KIT / rel, kind, _panel_levels(pc, key),
                                   unit=pc.get('unit') if kind == 'measurements' else None,
                                   allow_proposed=bool(pc.get('allow_proposed_stats')))
                    except ValidationError as e:
                        print(f'{rel} was not registered:\n  - ' + '\n  - '.join(e.problems))
                        return 1
    PV.register(rel, KIT, status=a.status, source=a.source, for_data=a.for_data)
    print(f'registered {rel}' + (f' for {a.for_data}' if a.for_data else f' as {a.status}'))
    return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    cmds = {'validate', 'check', 'stats-propose', 'templates', 'import-template', 'build', 'register'}
    if not argv or argv[0] not in cmds:
        res, fails = build(argv or None)
        return 1 if fails else 0
    cmd, rest = argv[0], argv[1:]
    if cmd == 'build':
        res, fails = build(rest or None)
        return 1 if fails else 0
    if cmd == 'validate':
        return 1 if validate() else 0
    if cmd == 'check':
        from . import check as C
        pdfs = sorted(p for p in (KIT / 'output' / 'figures').glob('*.pdf') if p.name.lower().startswith('fig'))
        if not pdfs:
            print('no figures in output/figures; run  python build.py  first')
            return 1
        import yaml
        res = []
        for p in pdfs:                                  # uses FigN.build.json next to each PDF when present
            cfgp = KIT / 'config' / f'{p.stem.lower()}.yaml'
            rules = (yaml.safe_load(cfgp.read_text()) or {}).get('accepted') if cfgp.exists() else None
            res.append(C.apply_accepted(C.check_pdf(p), rules))
        print(C.report_markdown(res))
        return 1 if any(r['errors'] for r in res) else 0
    if cmd == 'stats-propose':
        ap = argparse.ArgumentParser(prog='build.py stats-propose')
        ap.add_argument('file'); ap.add_argument('--test', default='twoway_tukey')
        ap.add_argument('--control', default=None); ap.add_argument('--alpha', type=float, default=0.05)
        a = ap.parse_args(rest)
        from .stats_helper import propose
        p, df = propose(a.file, a.test, a.control, a.alpha)
        print(f'wrote {p} ({len(df)} labels). Review it; the build does not read this file.')
        return 0
    if cmd == 'templates':
        from .templates import write_all
        for p in write_all():
            print('wrote', p)
        return 0
    if cmd == 'import-template':
        ap = argparse.ArgumentParser(prog='build.py import-template')
        ap.add_argument('file'); ap.add_argument('--out', required=True)
        a = ap.parse_args(rest)
        from .templates import import_workbook
        from .schemas import ValidationError
        try:
            res = import_workbook(a.file, a.out)
        except ValidationError as e:
            print('IMPORT REFUSED:\n  - ' + '\n  - '.join(e.problems), file=sys.stderr)
            return 1
        for p in res['written']:
            print('wrote', p)
        for p in res['unchanged']:
            print('unchanged', p)
        for sh in res['ignored']:
            print(f'ignored sheet {sh!r} (not a template sheet of this kit)')
        for r in res['registered']:
            print('registered', r, 'in data/PROVENANCE.json')
        for n in res.get('notes', []):
            print('note:', n)
        return 0
    if cmd == 'register':
        return register_cmd(rest)
    return 0
