#!/usr/bin/env python
# Prototype; not yet carefully vetted.
"""Read the existing supplied prior once; summarize without running CellOracle."""
import os,sys,csv,json,pickle,collections,hashlib,resource,time
ROOT=sys.argv[1]
SRC='/resnick/groups/NolanLab/shared_resources/source_files/celloracle_inputs'
PRIOR=SRC+'/base_GRN_TG_to_TF_dictionary.pkl'
def write(name,fields,rows):
 with open(ROOT+'/results/'+name,'w') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
start=time.time()
with open(PRIOR,'rb') as f:prior=pickle.load(f)
assert isinstance(prior,dict),type(prior)
counts=collections.Counter();n_entries=0;n_pairs=0
for target,regulators in prior.items():
 assert isinstance(target,str)
 unique=set(regulators)
 assert all(isinstance(tf,str) for tf in unique)
 counts.update(unique);n_entries+=len(regulators);n_pairs+=len(unique)
yuc5=sorted(set(prior.get('AT5G43890',[])))
write('supplied_prior_TF_target_counts.csv',['gene','prior_target_count'],[{'gene':tf,'prior_target_count':n} for tf,n in sorted(counts.items())])
write('supplied_prior_YUC5_regulators.csv',['source','target'],[{'source':tf,'target':'AT5G43890'} for tf in yuc5])
info={'source':PRIOR,'size_bytes':os.stat(PRIOR).st_size,'object_type':str(type(prior)),'target_genes':len(prior),'TFs':len(counts),'entries_including_duplicates':n_entries,'unique_target_TF_pairs':n_pairs,'YUC5_target_present':'AT5G43890' in prior,'YUC5_prior_regulators':len(yuc5),'elapsed_seconds':time.time()-start,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'python':sys.version,'executable':sys.executable,'interpretation':'The aggregate prior was supplied to inference. Membership does not identify ATAC/DAP/iGRN component provenance or prove binding in the assayed root cells.'}
with open(ROOT+'/results/supplied_prior_summary.json','w') as f:json.dump(info,f,indent=2)
print(json.dumps(info,indent=2),flush=True)
del prior
# Hash the exact small/medium inputs used for numerical postprocessing. Large
# H5AD count files were read only for metadata; identity is recorded separately.
paths=[]
with open(ROOT+'/results/source_file_manifest.csv') as f:
 paths.extend(row['path'] for row in csv.DictReader(f))
paths.extend([PRIOR,SRC+'/features.tsv.gz'])
rows=[]
for p in sorted(set(paths)):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 rows.append({'path':p,'size_bytes':os.stat(p).st_size,'sha256':h.hexdigest()})
write('source_sha256_manifest.csv',['path','size_bytes','sha256'],rows)
