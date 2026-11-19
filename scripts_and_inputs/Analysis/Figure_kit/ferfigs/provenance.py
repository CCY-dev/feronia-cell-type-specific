"""Data provenance registry: which version of each table the figures may draw, where it came from, and which values
each statistics table was made for.

data/PROVENANCE.json maps a table path (relative to the kit) to
    {"sha256": ..., "status": <data status>, "source": <plain words>}            for measurement tables
    {"sha256": ..., "for_data": <table>, "for_data_sha256": ..., "source": ...}  for statistics tables
`python build.py import-template` registers what it writes; `python build.py register` registers a table edited by
hand. At build time a table that differs from its registered version is drawn with the status 'unregistered' and a
data warning, and a statistics table made for other values than the ones drawn gives a 'stale statistics' warning.
The registry records; it never changes a table.

Every build also records which files it actually opened (an audit hook on file opens, active only while a figure is
built), so every input under data/ - measurement tables, statistics, printed values, upstream tables, layout tables,
image manifests and image files - is checked against the registry, not only the tables a panel names. Each entry
carries a role (measurements, statistics, printed_values, upstream_table, layout, image_manifest, image, custom) so
mixed evidence within one panel stays visible.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from contextlib import contextmanager
from pathlib import Path

from .style import KIT

REGISTRY = Path('data') / 'PROVENANCE.json'
ROLES = ('measurements', 'statistics', 'printed_values', 'upstream_table', 'layout', 'image_manifest', 'image', 'custom')


def _audit(event, args):
    """Audit hook: while a build is being recorded, remember every file opened for reading."""
    if event != 'open':
        return
    sink = getattr(sys, '_ferfigs_reads', None)
    if sink is None:
        return
    path, mode, flags = (list(args) + [None, None, None])[:3]
    if path is None or isinstance(path, int):
        return
    if isinstance(mode, str):
        if any(c in mode for c in 'wax+'):
            return
    elif isinstance(flags, int) and flags & (os.O_WRONLY | os.O_RDWR):
        return
    try:
        sink.add(os.path.realpath(os.fsdecode(path)))
    except (TypeError, ValueError):
        pass


if not getattr(sys, '_ferfigs_audit_installed', False):      # one hook per process, shared by every kit copy
    sys.addaudithook(_audit)
    sys._ferfigs_audit_installed = True


@contextmanager
def recording():
    """with recording() as reads: ...  -> the set of absolute paths opened for reading inside the block."""
    prev = getattr(sys, '_ferfigs_reads', None)
    sink = set()
    sys._ferfigs_reads = sink
    try:
        yield sink
    finally:
        sys._ferfigs_reads = prev


def data_inputs(reads, kit: Path = KIT):
    """Kit-relative paths under data/ among the recorded reads (the registry itself excluded), sorted."""
    root = Path(os.path.realpath(Path(kit)))
    out = set()
    for r in reads:
        try:
            rel = Path(r).relative_to(root).as_posix()
        except ValueError:
            continue
        if rel.startswith('data/') and rel != REGISTRY.as_posix():
            out.add(rel)
    return sorted(out)


def file_status(rel, kit: Path = KIT, reg=None):
    """(status, role, warning or None) for any input file. Statistics entries report status 'statistics'."""
    reg = load(kit) if reg is None else reg
    rel = Path(rel).as_posix()
    e = reg.get(rel)
    if e is None:
        return 'unregistered', 'unknown', (f'{rel} is not in {REGISTRY.as_posix()}: its provenance is unknown. Import '
                                           f'it with the template, or register it: python build.py register {rel} '
                                           f'--status ... --source "who supplied it".')
    role = e.get('role') or ('statistics' if e.get('for_data') else 'measurements')
    if sha256(Path(kit) / rel) != e['sha256']:
        return 'unregistered', role, (f'{rel} has changed since it was registered ({e.get("source", "")}). Register '
                                      f'the new version (python build.py register {rel} ...) or import it with the '
                                      f'template.')
    return e.get('status', 'statistics' if e.get('for_data') else 'unspecified'), role, None


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(kit: Path = KIT) -> dict:
    p = Path(kit) / REGISTRY
    return json.loads(p.read_text()) if p.exists() else {}


def save(reg: dict, kit: Path = KIT):
    p = Path(kit) / REGISTRY
    tmp = p.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(dict(sorted(reg.items())), indent=1, ensure_ascii=False) + '\n')
    tmp.replace(p)


def register(rel, kit: Path = KIT, status=None, source='', for_data=None, reg=None, role=None):
    """Record the current bytes of a table. Data tables carry a status; statistics tables name the data table (and
    its current checksum) they were made for. role: one of ROLES (kept from an earlier entry when not given).
    Returns the updated registry (saved unless reg is passed in)."""
    own = reg is None
    reg = load(kit) if own else reg
    rel = Path(rel).as_posix()
    old = reg.get(rel) or {}
    entry = dict(sha256=sha256(Path(kit) / rel), source=source)
    if for_data:
        fd = Path(for_data).as_posix()
        entry.update(for_data=fd, for_data_sha256=sha256(Path(kit) / fd), role='statistics')
    else:
        entry['status'] = status or 'unspecified'
        entry['role'] = role or old.get('role') or 'measurements'
    if role and role not in ROLES:
        raise ValueError(f'role {role!r} is not one of {ROLES}')
    reg[rel] = entry
    if own:
        save(reg, kit)
    return reg


def data_status(rel, configured, kit: Path = KIT, reg=None):
    """(status, warning or None) for a data table about to be drawn."""
    reg = load(kit) if reg is None else reg
    rel = Path(rel).as_posix()
    p = Path(kit) / rel
    e = reg.get(rel)
    if e is None:
        if not reg:                      # kit without a registry: config status, as before
            return configured, None
        return 'unregistered', (f'{rel} is not in {REGISTRY.as_posix()}: its provenance is unknown. Import it with '
                                f'the template, or register it: python build.py register {rel} --status original '
                                f'--source "who supplied it".')
    if not p.exists():
        return configured, None
    if sha256(p) != e['sha256']:
        return 'unregistered', (f'{rel} has changed since it was registered ({e.get("source", "")}). Register the '
                                f'new version (python build.py register {rel} --status ... --source ...) or import '
                                f'it with the template.')
    return e.get('status', configured), None


def stats_check(stats_rel, data_rel, kit: Path = KIT, reg=None):
    """Warning text when a statistics table was made for other values than the data table drawn, else None."""
    reg = load(kit) if reg is None else reg
    if not reg:
        return None
    s, d = Path(stats_rel).as_posix(), Path(data_rel).as_posix()
    e = reg.get(s)
    if e is None or e.get('for_data') != d:
        return (f'{s} is not registered for {d}: the statistics labels cannot be traced to these values. Supply '
                f'them with the values (template stats sheet), or confirm: python build.py register {s} --for {d} '
                f'--source "who computed them and how".')
    if sha256(Path(kit) / s) != e['sha256']:
        return (f'{s} has changed since it was registered; confirm it: python build.py register {s} --for {d} '
                f'--source "...".')
    if sha256(Path(kit) / d) != e['for_data_sha256']:
        return (f'stale statistics: {s} was made for a different version of {d}. Re-supply the labels and the test '
                f'with the new values, or confirm they still apply: python build.py register {s} --for {d} '
                f'--source "...".')
    return None
