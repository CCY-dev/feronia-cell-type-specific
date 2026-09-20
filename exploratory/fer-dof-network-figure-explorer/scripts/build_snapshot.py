#!/usr/bin/env python3
# Prototype; not yet carefully vetted.
"""Assemble the bounded figure-sharing snapshot from an existing manuscript checkout.

Requires --source-root. This copies saved files, never reruns analysis or publishes.
Preserves original bytes and the already marked code snapshots from the v6 ZIP.
Do not use this on a package after its coordinator has edited or published it.
"""
from pathlib import Path, PurePosixPath
import argparse, csv, hashlib, io, json, stat, zipfile

PKG=Path(__file__).resolve().parents[1]
C='revisions/2026-09-20_compact_figures'
AUDIT='.living/tasks/fer-fig6-independent-audit-2026-09-20'
D='revisions/2026-09-20_dof_subnetworks'
Q='revisions/2026-09-08_routeA_execution/q04_network'
B='revisions/2026-09-20_fig6_brainstorm'

def sha(data): return hashlib.sha256(data).hexdigest()
def write_csv(path,rows,fields=None):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-root',required=True,type=Path)
    args=ap.parse_args(); source=args.source_root.resolve()
    assert (source/'CLAUDE.md').is_file(),source
    archive=source/C/'review_v6/FER_fig6_prototype_code_v6.zip'
    assert sha(archive.read_bytes())=='77e33ed99c67a7d9f89768b166762689e9feadca22764f69b70eaa7c4b125ae3'
    provenance={}
    def save(rel,data,origin,original_sha,role,preservation='Exact source bytes',container=''):
        rel=PurePosixPath(rel)
        assert not rel.is_absolute() and '..' not in rel.parts and '.git' not in rel.parts,rel
        dest=PKG/str(rel);dest.parent.mkdir(parents=True,exist_ok=True)
        if dest.exists(): assert dest.read_bytes()==data,f'Refusing to overwrite changed package file: {rel}'
        else: dest.write_bytes(data)
        provenance[str(rel)]={'package_path':str(rel),'original_repository_path':origin,
            'original_source_sha256':original_sha,'package_sha256':sha(data),'bytes':len(data),
            'role':role,'preservation':preservation,'source_container':container}
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        manifest={r['archive_path']:r for r in csv.DictReader(io.StringIO(z.read('BUNDLE_FILES.csv').decode()))}
        references={r['original_repo_relative_path']:r for r in csv.DictReader(io.StringIO(z.read('reference_metadata/source_copy_manifest.csv').decode()))}
        for item in z.infolist():
            rel=PurePosixPath(item.filename)
            assert not item.is_dir() and not rel.is_absolute() and '..' not in rel.parts,rel
            assert not stat.S_ISLNK(item.external_attr>>16),rel
            # The original v6 entry point is retained under its repository path.
            if str(rel)=='START_HERE.md': continue
            data=z.read(item);row=manifest.get(str(rel));ref=references.get(str(rel))
            if row:
                assert sha(data)==row['archived_sha256'],rel
                original_sha=ref['original_sha256'] if ref else row['source_sha256']
                preservation='Already marked code snapshot' if original_sha!=sha(data) else 'Exact source bytes'
                origin=str(rel) if not str(rel).startswith('reference_metadata/') else str(Path(row['source_file']).relative_to(source))
                save(rel,data,origin,original_sha,'V6 base: '+row['role'],preservation,str(archive.relative_to(source)))
            else:
                assert str(rel)=='BUNDLE_FILES.csv'
                save(rel,data,'','', 'Original v6 bundle manifest','Exact archive entry',str(archive.relative_to(source)))
    def add(rel,role):
        rel=str(rel);p=source/rel;assert p.is_file(),p
        data=p.read_bytes()
        if rel in provenance:
            old=provenance[rel]
            assert old['original_source_sha256']==sha(data) or old['package_sha256']==sha(data),rel
            return
        save(rel,data,rel,sha(data),role)
    def add_tree(rel,role,suffixes=None):
        for p in sorted((source/rel).rglob('*')):
            if p.is_file() and (suffixes is None or p.suffix in suffixes): add(p.relative_to(source),role)

    inventory=f'{C}/review_v7/provenance/CURRENT_INPUT_OUTPUT_CODE.csv'
    rows=list(csv.DictReader((source/inventory).open()));assert len(rows)==6
    for row in rows:
        for field in ['input_tables_read_by_plot_script','archived_plot_inputs','calculation_code','plotting_code','outputs']:
            for rel in row[field].split(' ; '): add(rel,f"{row['panel']}: {field}")
    add(inventory,'Current six-panel input/output/code map')
    for name in ['README.md','PDF_VERIFICATION.json','make_review_pdf.py','make_review_with_sources.py',
                 'review_at_print_size_v7.pdf','review_with_sources_v7.pdf','sources_and_validation.pdf','review-1.png']:
        add(f'{C}/review_v7/{name}','V7 proof, preview and source generator')
    add_tree(f'{C}/review_v7/networks','V7 panel D code, tables and validation',
             {'.R','.py','.md','.csv','.json','.txt'})
    for name in ['PDF_VERIFICATION.json']:
        add(f'{C}/review_v6/{name}','Recorded v6 proof verification')
    for name in ['SOURCE_AVAILABILITY.md','SPATIAL_MECHANISM_REVIEW.md']:
        add(f'{C}/review_v2/{name}','Source availability or manuscript mechanism review')
    add_tree(f'{C}/review_v2/audit_response','Endodermis source-triage evidence')
    # Exact explanatory tables behind the network's support classes.
    for rel in [f'{C}/review_v2/all_contexts/all_contexts_as_exported_source_context_edges.csv',
                f'{C}/review_v2/all_contexts/README.md',f'{C}/review_v2/all_contexts/validation.csv',
                f'{C}/review_v2/all_contexts/all_contexts_overview_compact.pdf',
                f'{C}/networks/all_tissues_overview_compact_union_edges.csv',
                f'{C}/networks/CAPTIONS.md',f'{C}/networks/README.md',
                f'{C}/upset/README.md',
                f'{D}/inventory/tested10_retained_edges.csv',f'{D}/inventory/README.md',
                f'{D}/inventory/endodermis_saved_object_check.json',
                f'{D}/inventory/check_endodermis_saved_object_remote.py',
                f'{D}/inventory/full_public_export_degree_audit.csv',
                f'{D}/inventory/tested10_edge_score_count_mismatches.csv',
                f'{D}/inventory/tested10_identities.csv',
                f'{D}/inventory/observed_auxin_biosynthesis_target_input_eligibility.csv',
                f'{D}/inventory/basal_DE_extraction_provenance.csv',
                f'{Q}/README.md',f'{Q}/results/h5ad_gene_universes.csv',
                f'{Q}/results/DOF_YUC5_target_eligibility.csv',
                f'{Q}/results/source_file_manifest.csv',f'{Q}/results/source_sha256_manifest.csv',
                f'{B}/network_audit.md']:
        if not (source/rel).is_file():
            if rel.endswith('/CAPTIONS.md'): continue
            raise FileNotFoundError(rel)
        add(rel,'Selected source evidence and context provenance')
    # Selected independent-audit deliverables, not task histories or prompts.
    for name in ['AUDIT_REPORT.md','FIRST_ASSESSMENT.md','DESIGN_ADDENDUM_v5.md','findings.csv','OUTPUT_SHA256.txt']:
        add(f'{AUDIT}/outputs/claude-science/{name}','Independent audit deliverable')
    add_tree(f'{AUDIT}/outputs/claude-science/checks','Independent audit code and saved results')
    add_tree(f'{AUDIT}/outputs/claude-science/scratch_rerun','Independent audit scratch code and session record',{'.R','.py','.txt','.csv','.json'})
    for rel in [f'{AUDIT}/run-manifest.json',
                f'{AUDIT}/receipts/e79a0475-81ab-4fa3-b94b-8f01bccd9734.json',
                f'{AUDIT}/receipts/c1e0db4a-c32e-4fb4-8e55-8e2ba54eb40b.json']:
        add(rel,'Independent audit frozen-input manifest or receipt')
    # Include metadata at the exact original paths linked from the source appendix.
    for name in ['README.md','source_copy_manifest.csv','dependency_manifest.csv','dependency_source_excerpts.csv']:
        add(f'{C}/review_v3/code_reference/{name}','Historical code-reference provenance')
    # The v3 provenance guide is required by the source-guide links.
    for name in ['README.md','INPUT_OUTPUT_CODE.csv','LOCAL_FILES_SHA256.csv','validation.json']:
        add(f'{C}/review_v3/provenance/{name}','Historical v2 source trace')
    # All auditor hashes are rechecked in the source checkout; omitted audit
    # rasters remain explicitly outside this bounded sharing package.
    audit_scope=[]; auditdir=source/AUDIT/'outputs/claude-science'
    for line in (auditdir/'OUTPUT_SHA256.txt').read_text().splitlines():
        if not line.strip(): continue
        expected,rel=line.split(maxsplit=1);rel=rel.lstrip('*')
        p=auditdir/rel
        if not p.exists(): p=source/rel
        assert p.is_file(),p
        observed=sha(p.read_bytes());assert observed==expected,p
        repo_rel=str(p.relative_to(source))
        audit_scope.append({'audit_file':rel,'original_repository_path':repo_rel,'recorded_sha256':expected,
                            'source_hash_matches':True,'included_in_package':repo_rel in provenance,
                            'omission_reason':'' if repo_rel in provenance else 'Historical audit scratch raster omitted; source hash checked'})
    write_csv(PKG/'provenance/AUDIT_HASH_SCOPE.csv',audit_scope)
    write_csv(PKG/'provenance/SOURCE_PROVENANCE.csv',[provenance[k] for k in sorted(provenance)])
    (PKG/'provenance/SNAPSHOT_STATUS.json').write_text(json.dumps({
        'title':'FER-DOF Network Figure Explorer','status':'Prototype; not yet carefully vetted',
        'source_zip_sha256':sha(archive.read_bytes()),'copied_files':len(provenance),
        'copied_bytes':sum(r['bytes'] for r in provenance.values()),'panel_types':len(rows),
        'audit_output_hashes_checked':len(audit_scope),
        'audit_outputs_included':sum(r['included_in_package'] for r in audit_scope),
        'scope':'Saved outputs and selected sources; not a self-contained portable raw-analysis pipeline',
        'omitted':['raw HPC objects','complete public source checkout','Git metadata','full manuscript',
                   'task histories and prompts','older proof rasters','historical audit scratch rasters'],
        'source_changes':False,'analysis_rerun':False,'published':False},indent=2)+'\n')
    print(f'Assembled {len(provenance)} source files, {sum(r["bytes"] for r in provenance.values()):,} bytes at {PKG}')

if __name__=='__main__': main()
