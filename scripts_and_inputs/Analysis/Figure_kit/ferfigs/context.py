"""Build context handed to every figures/figN.py module.

A figure module defines   build(ctx) -> Page   and reads everything through ctx:
  ctx.cfg                 the figure config (config/figN.yaml) as a dict
  ctx.table(path, kind)   a validated table (path relative to the kit folder; kind from schemas.SCHEMAS)
  ctx.path(rel)           absolute path of a kit file
  ctx.page()              a Page of the configured size
  ctx.record(panel, rec)  store what a panel drew (bars, n, labels) for the checks and the Source Data file
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import yaml

from .layout import Page
from .schemas import read_table, ValidationError
from .style import KIT


class Context:
    def __init__(self, fig: str, kit: Path = KIT):
        self.fig = fig
        self.kit = Path(kit)
        cfg_path = self.kit / 'config' / f'{fig}.yaml'
        if not cfg_path.exists():
            raise FileNotFoundError(f'no config for {fig}: {cfg_path}')
        with open(cfg_path) as f:
            self.cfg = yaml.safe_load(f)
        self.tables = {}
        self.records = {}
        self.warnings = []
        self.data_warnings = []          # provenance problems (unregistered tables, stale statistics)
        self._page = None
        self._registry = None
        self._warned = set()             # input files already warned about (one warning per file)

    def path(self, rel):
        p = Path(rel)
        return p if p.is_absolute() else self.kit / p

    def registry(self):
        from . import provenance as P
        if self._registry is None:
            self._registry = P.load(self.kit)
        return self._registry

    def data_warn(self, text):
        if text not in self.data_warnings:
            self.data_warnings.append(text)

    def panel_for(self, rel):
        """The dotbar panel config whose data or statistics table is rel (data/stats or data_file/stats_file)."""
        for pc in (self.cfg.get('panels') or {}).values():
            if pc.get('type') == 'dotbar' and str(rel) in (pc.get('data'), pc.get('stats'), pc.get('data_file'),
                                                            pc.get('stats_file')):
                return pc
        return None

    def labels_omitted(self, rel):
        """True when rel is absent and this copy of the kit ships without the authors' label tables (the public copy
        carries data/LABELS_NOT_INCLUDED). The labels are then not drawn; the published figure shows them."""
        if self.path(rel).exists() or not self.path('data/LABELS_NOT_INCLUDED').exists():
            return False
        self.data_warn(f'{rel}: author labels are not in this copy of the repository; drawn without them '
                       f'(the labels are in the published figure).')
        return True

    def table(self, rel, kind, levels=None):
        if kind == 'annotations' and self.labels_omitted(rel):
            import pandas as pd
            df = pd.DataFrame(columns=['group', 'series', 'label', 'kind', 'group2', 'series2'], dtype=str)
            self.tables[str(rel)] = dict(kind=kind, rows=0, path=str(rel))
            return df
        pc = self.panel_for(rel)
        if levels is None and pc is not None:   # a dotbar panel's tables: names checked against its config lists
            levels = {'group': pc['groups'], 'series': pc['series']}
            if kind == 'annotations':
                levels.update({'group2': pc['groups'], 'series2': pc['series']})
        unit = pc.get('unit') if (pc is not None and kind == 'measurements') else None
        allow = bool(pc.get('allow_proposed_stats')) if pc is not None else False
        df = read_table(self.path(rel), kind, levels, unit=unit, allow_proposed=allow)
        if allow and kind == 'annotations' and 'source' in df.columns and (df['source'] == 'proposed_by_helper').any():
            self.data_warn(f'{rel}: statistics proposed by the helper are drawn on purpose (allow_proposed_stats: '
                           f'true); they are not the authors\' analysis.')
        self.tables[str(rel)] = dict(kind=kind, rows=len(df), path=str(rel))
        self.warnings += df.attrs.get('warnings', [])
        return df

    def custom_table(self, rel, required, numeric=(), hint=''):
        """Read a figure-specific CSV that is not one of the three standard schemas: the file must exist, carry the
        required columns, and the numeric columns must parse. Recorded in the build report like any table."""
        import numpy as np
        import pandas as pd
        from .schemas import ValidationError
        p = self.path(rel)
        if self.labels_omitted(rel):
            return pd.DataFrame(columns=list(required), dtype=str)
        if not p.exists():
            raise ValidationError([f'{rel}: file not found.' + (f' {hint}' if hint else '')])
        df = pd.read_csv(p, dtype=str, keep_default_na=False)
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValidationError([f'{rel}: missing column(s) {missing}; required: {list(required)}.'])
        probs = []
        for col in numeric:
            bad = pd.to_numeric(df[col], errors='coerce').isna()
            if bad.any():
                probs.append(f'{rel} line {int(np.where(bad)[0][0]) + 2}: {col} = {df[col][bad].iloc[0]!r} is not a number.')
            else:
                df[col] = pd.to_numeric(df[col])
        if probs:
            raise ValidationError(probs)
        self.tables[str(rel)] = dict(kind='custom', rows=len(df), path=str(rel))
        return df

    def page(self):
        w, h = self.cfg['size_mm']
        self._page = Page(w, h, self.fig)
        return self._page

    def note(self, key, value):
        """Non-panel build information (layout checks, row bottoms) stored in the build record under '_notes';
        not a panel, so it never appears in data_status or the Source Data workbook."""
        if self._page is not None:
            self._page.records.setdefault('_notes', {})[key] = value

    def record(self, panel, rec):
        rec = dict(rec)
        pc = (self.cfg.get('panels') or {}).get(panel, {})
        rec['data_status'] = rec.get('data_status') or pc.get('data_status', 'unspecified')
        for k in ('data', 'stats', 'manifest'):
            if k in pc:
                rec[k] = pc[k]
        self._provenance(panel, rec)
        self.records[panel] = rec
        if self._page is not None:
            self._page.records[panel] = rec

    def input_status(self, rel, panel='inputs'):
        """(status, role) of one input file, warning once per file when it is unregistered or changed."""
        from . import provenance as P
        st, role, w = P.file_status(rel, self.kit, self.registry())
        if w and rel not in self._warned:
            self._warned.add(rel)
            self.data_warn(f'{panel}: {w}')
        return st, role

    def _provenance(self, panel, rec):
        """Status of every input a panel names (data, statistics, image manifest, source tables) from the registry.
        Any input that is unregistered or changed makes the panel's status 'unregistered'; statistics made for other
        values give a stale-statistics warning. The per-file roles and statuses go into rec['inputs']."""
        from . import provenance as P
        if not self.registry():            # a kit without a registry: configured statuses, as before
            return
        inputs, unknown = {}, False
        data = rec.get('data') if rec.get('data') and self.path(rec['data']).exists() else None
        if data:
            st, w = P.data_status(data, rec['data_status'], self.kit, self.registry())
            rec['data_status'] = st
            if w and data not in self._warned:
                self._warned.add(data)
                self.data_warn(f'{panel}: {w}')
            inputs[data] = dict(role=(self.registry().get(data) or {}).get('role', 'measurements'), status=st)
            unknown |= st == 'unregistered'
            if rec.get('stats') and self.path(rec['stats']).exists():
                srel = rec['stats']
                w2 = P.stats_check(srel, data, self.kit, self.registry())
                if w2 and srel not in self._warned:
                    self._warned.add(srel)
                    self.data_warn(f'{panel}: {w2}')
                inputs[srel] = dict(role='statistics', status='not confirmed for these values' if w2 else 'statistics')
        others = [rec.get('manifest')] + list((rec.get('source_tables') or {}).values())
        for rel in others:
            if not rel or rel in inputs or not str(rel).startswith('data/') or not self.path(rel).exists():
                continue
            st, role = self.input_status(str(rel), panel)
            inputs[str(rel)] = dict(role=role, status=st)
            unknown |= st == 'unregistered'
        if unknown:
            rec['data_status'] = 'unregistered'
        elif any(v['status'] == 'author_supplied' for v in inputs.values()):
            rec['data_status'] = 'author_supplied'    # values replaced through a template; per-file list says which
        if inputs:
            rec['inputs'] = dict(sorted(inputs.items()))

    @property
    def panel_records(self):
        return {k: v for k, v in self.records.items() if not k.startswith('_')}


def load_module(fig: str, kit: Path = KIT):
    p = Path(kit) / 'figures' / f'{fig}.py'
    spec = importlib.util.spec_from_file_location(f'ferfigs_fig_{fig}', p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def figures(kit: Path = KIT):
    return sorted(p.stem for p in (Path(kit) / 'config').glob('fig*.yaml'))
