#!/usr/bin/env python3
# Prototype; not yet carefully vetted.
"""Verify copied-source bytes, current figures, inputs and entry-point links.

Local file-identity and packaging checks only. Does not validate scientific claims.
The optional --source-root checks hashes against the originating checkout.
"""
from pathlib import Path
import argparse, csv, hashlib, json, re
from urllib.parse import unquote, urlparse

PKG=Path(__file__).resolve().parents[1]
C='revisions/2026-09-20_compact_figures'
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def pieces(value): return [p.strip() for p in value.split(';') if p.strip()]

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source-root',type=Path)
    args=ap.parse_args();source=args.source_root.resolve() if args.source_root else None
    rows=list(csv.DictReader((PKG/'provenance/SOURCE_PROVENANCE.csv').open()))
    for row in rows:
        path=PKG/row['package_path'];assert path.is_file(),path
        assert sha(path)==row['package_sha256'],path
        assert path.stat().st_size==int(row['bytes']),path
        if source and row['original_repository_path']:
            original=source/row['original_repository_path']
            assert original.is_file(),original
            assert sha(original)==row['original_source_sha256'],original
    panels=list(csv.DictReader((PKG/C/'review_v7/provenance/CURRENT_INPUT_OUTPUT_CODE.csv').open()))
    assert len(panels)==6
    inputs=set();outputs=[]
    for row in panels:
        for raw in pieces(row['input_tables_read_by_plot_script']):
            assert (PKG/raw).is_file(),raw;inputs.add(raw)
        for raw in pieces(row['outputs']):
            assert (PKG/raw).is_file(),raw;outputs.append(raw)
        for raw in pieces(row['plotting_code']): assert (PKG/raw).is_file(),raw
    assert len(inputs)==7 and len(outputs)==18
    assert sum(p.endswith('.svg') for p in outputs)==6
    linkcount=0
    for name in ['README.md','SOURCES.md','VALIDATION_STATUS.md']:
        for target in re.findall(r'\]\(([^)]+)\)',(PKG/name).read_text()):
            parsed=urlparse(target)
            if parsed.scheme or target.startswith('#'): continue
            path=PKG/unquote(parsed.path);assert path.exists(),(name,target)
            linkcount+=1
    forbidden=[]
    for path in PKG.rglob('*'):
        if not path.is_file(): continue
        if '.git' in path.relative_to(PKG).parts: continue
        assert not path.is_symlink(),path
        if path.suffix.lower() in {'.rds','.h5ad','.h5','.hdf5','.loom','.docx','.pptx'} or path.name in {'.env','credentials.json'}:
            forbidden.append(str(path.relative_to(PKG)))
    assert not forbidden,forbidden
    result={'status':'Packaging checks passed; not scientific validation','copied_source_files':len(rows),
        'original_source_hashes_checked':bool(source),'panels':len(panels),'unique_immediate_input_tables':len(inputs),
        'panel_output_files':len(outputs),'editable_SVGs':6,'root_navigation_local_links':linkcount,
        'forbidden_large_objects_or_document_types':forbidden}
    (PKG/'provenance/PACKAGE_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
