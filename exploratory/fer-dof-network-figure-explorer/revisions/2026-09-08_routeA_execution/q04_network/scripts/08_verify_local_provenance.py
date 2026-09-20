#!/usr/bin/env python3
# Prototype; not yet carefully vetted.
"""Local artifact/source hash verification only; all statistical analysis is remote."""
import csv,hashlib,sys
from pathlib import Path
root=Path(sys.argv[1]).resolve()
public=Path(sys.argv[2]).resolve()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 return h.hexdigest()
rows=[]
with (root/'results/source_sha256_manifest.csv').open() as f:
 for r in csv.DictReader(f):
  name=Path(r['path']).name
  local=next((p for p in [public/'outputs'/name,public/'inputs'/name] if p.exists()),None)
  if local is None:continue
  digest=sha(local)
  rows.append(dict(source_path=r['path'],public_path=str(local),source_sha256=r['sha256'],public_sha256=digest,identical=digest==r['sha256']))
with (root/'provenance/public_source_hash_comparison.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print('Public/shared source byte-identities:',sum(r['identical'] for r in rows),'/',len(rows))
for r in rows:
 if not r['identical']:print('DIFFERENT',r['source_path'])
rows=[]
for p in sorted(root.rglob('*')):
 if not p.is_file() or p.name in ['local_artifact_sha256.tsv']:continue
 rows.append(dict(path=str(p.relative_to(root)),size_bytes=p.stat().st_size,sha256=sha(p)))
with (root/'provenance/local_artifact_sha256.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['path','size_bytes','sha256'],delimiter='\t');w.writeheader();w.writerows(rows)
print('Local artifacts hashed:',len(rows))
