#!/usr/bin/env python3
# Prototype; not yet carefully vetted.
"""Build a v6 code/provenance reference after explicit --create-zip invocation.

Default invocation checks the proposed contents without creating a ZIP. Original
files and earlier bundles remain untouched. Run only after proof sources settle.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import io
import zipfile

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
V5 = HERE.parent/'review_v5'
REFERENCE = HERE.parent/'review_v3/code_reference'
INVENTORY = HERE/'provenance/CURRENT_INPUT_OUTPUT_CODE.csv'
SUFFIXES = {'.R','.py','.csv','.md','.txt','.json','.sbatch'}
CODE_SUFFIXES = {'.R','.py','.sbatch'}
MARKER = b'# Prototype; not yet carefully vetted.\n'

def digest(data):
    return hashlib.sha256(data).hexdigest()

def marked_bytes(data, suffix):
    if suffix not in CODE_SUFFIXES or b'not yet carefully vetted' in data[:700].lower():
        return data,False
    if data.startswith(b'#!'):
        first,sep,rest=data.partition(b'\n')
        return first+sep+MARKER+rest,True
    return MARKER+data,True

def collect():
    entries={};manifest=[];sources={}
    def add(source, archive_path=None, role='reference'):
        source=Path(source)
        assert source.is_file(),source
        key=str(archive_path if archive_path is not None else source.relative_to(REPO))
        assert not key.startswith('/') and '..' not in Path(key).parts,key
        data=source.read_bytes()
        archived,header_added=marked_bytes(data,source.suffix)
        if key in entries:
            assert entries[key]==archived,f'Different contents for duplicate archive path: {key}'
            return
        entries[key]=archived;sources[source]=digest(data)
        manifest.append({'archive_path':key,'source_file':str(source),'role':role,
                         'source_sha256':digest(data),'archived_sha256':digest(archived),
                         'source_bytes':len(data),'archived_bytes':len(archived),
                         'prototype_header_added_in_archive':header_added})

    # Existing marked upstream reference copies, preserving their original paths.
    with (REFERENCE/'source_copy_manifest.csv').open() as f:
        for row in csv.DictReader(f):
            copy=Path(row['copy_absolute_path'])
            assert digest(copy.read_bytes())==row['copy_sha256'],copy
            add(copy,row['original_repo_relative_path'],'historical marked reference copy')

    # Current selected plotting/calculation code, including dependencies omitted
    # from the old reference snapshots. Reference copies take precedence where
    # the original path was already included with a prototype comment.
    with INVENTORY.open() as f:
        panel_rows=list(csv.DictReader(f))
    assert len(panel_rows)==6
    for row in panel_rows:
        assert row['status']=='Prototype; not yet carefully vetted'
        for field in ('calculation_code','plotting_code'):
            for rel in row[field].split(' ; '):
                if rel not in entries:
                    add(REPO/rel,role='current inventory code')
    matrix=next(row for row in panel_rows if row['panel'].startswith('c:'))
    assert matrix['input_tables_read_by_plot_script']=='revisions/2026-09-20_fig6_brainstorm/celltype_prototypes/celltype_summary_all37.csv'

    # Preserve v5 forest code, selected data, metadata and complete reconstruction.
    for source in sorted((V5/'quantitative').iterdir()):
        if source.suffix in SUFFIXES and ('matrix' not in source.name and 'celltype' not in source.name):
            add(source,role='v5 forest dependency or metadata')
    for source in sorted((V5/'provenance').iterdir()):
        if source.is_file() and source.suffix in SUFFIXES:
            add(source,role='v5 membership reconstruction or provenance')
    # All nine reconstruction table inputs, including comparison-only tables.
    with (V5/'provenance/input_sha256.csv').open() as f:
        for row in csv.DictReader(f):
            source=Path(row['path'])
            assert digest(source.read_bytes())==row['sha256'],source
            add(source,role='membership reconstruction input: '+row['role'])

    # Immediate current plot inputs include the forest's older compact-table path.
    for row in panel_rows[:3]:
        for rel in row['input_tables_read_by_plot_script'].split(' ; '):
            add(REPO/rel,role='current panel immediate input')
    for dirname in ('quantitative','overlap','provenance'):
        for source in sorted((HERE/dirname).iterdir()):
            if source.is_file() and source.suffix in SUFFIXES:
                add(source,role='v6 script, table or documentation')

    # These generators/docs are owned by the proof coordinator. Read only after
    # completion signal; absence is a failed check, never a partial final bundle.
    for name in ('make_review_pdf.py','make_review_with_sources.py','make_code_bundle.py',
                 'CAPTIONS_AND_VALIDATION.md','START_HERE.md','README.md',
                 'VALIDATION_STATUS.md','AUDIT_RECEIPT_VERIFICATION.md'):
        add(HERE/name,role='v6 proof or bundle documentation')
    for name in ('source_copy_manifest.csv','dependency_manifest.csv',
                 'dependency_source_excerpts.csv','README.md'):
        add(REFERENCE/name,Path('reference_metadata')/name,'historical reference metadata')
    add(HERE/'START_HERE.md','START_HERE.md','archive entry point')

    # Ensure the snapshot was not assembled across concurrent source changes.
    for source,expected in sources.items():
        assert digest(source.read_bytes())==expected,f'Source changed during collection: {source}'
    buf=io.StringIO();w=csv.DictWriter(buf,fieldnames=list(manifest[0]));w.writeheader();w.writerows(manifest)
    entries['BUNDLE_FILES.csv']=buf.getvalue().encode()
    return entries

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--create-zip',action='store_true',help='Create ZIP after the coordinator signals that proof sources are stable.')
    args=parser.parse_args()
    entries=collect()
    print(f'Checked {len(entries)} proposed archive entries; original files unchanged.')
    if not args.create_zip:
        print('Check only. ZIP not created. Await the final bundle signal.')
        return
    output=HERE/'FER_fig6_prototype_code_v6.zip'
    if output.exists():
        raise FileExistsError(f'Preserve the existing archive; choose a new version before recreating: {output}')
    with zipfile.ZipFile(output,'w',compression=zipfile.ZIP_DEFLATED) as archive:
        for name,data in sorted(entries.items()):
            archive.writestr(name,data)
    with zipfile.ZipFile(output) as archive:
        assert archive.testzip() is None
        for name,data in entries.items():
            assert digest(archive.read(name))==digest(data),name
    print(output)
    print(f'{output.stat().st_size:,} bytes; SHA-256 {digest(output.read_bytes())}')

if __name__=='__main__':
    main()
