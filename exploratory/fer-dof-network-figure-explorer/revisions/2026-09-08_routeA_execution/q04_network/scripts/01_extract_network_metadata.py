#!/usr/bin/env python
# Prototype; not yet carefully vetted.
"""Scheduled, read-only source extraction. Never load counts or infer networks."""
import os, sys, csv, json, hashlib, collections
import h5py
ROOT=sys.argv[1]
FER='/resnick/groups/NolanLab/cchen/Feronia'
CO=FER+'/analysis/CellOracle/CellOracle_output/20260317_feronia_salinity_timecourse_WT_untreated_fer4_untreated_newIntegratedSeurat3_NoSuspiciousCluster'
SOURCE=FER+'/analysis/CellOracle/CellOracle_input/20250509_4S_feronia_salinity_time_course_untreated_celloracle_input.h5ad'
def dec(x):return x.decode() if isinstance(x,bytes) else str(x)
def write(name, fields, rows):
 with open(ROOT+'/results/'+name,'w') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
rows=[]
for tissue,path in [('all',SOURCE)]+[(t,CO+'/'+t+'_rawRNA_celloracle.h5ad') for t in ['cortex','endodermis','trichoblast','atrichoblast','pericycle','phloem','xylem','procambium','lrc','columella']]:
 with h5py.File(path,'r') as f:
  index=dec(f['var'].attrs['_index'])
  ids=[dec(x) for x in f['var'][index][:]]
  rows.extend({'tissue':tissue,'gene':x} for x in ids)
  print(tissue,len(ids),path,flush=True)
write('h5ad_gene_universes.csv',['tissue','gene'],rows)
with h5py.File(SOURCE,'r') as f:
 obj=f['obs']['time.celltype.anno.Li'];cats=[dec(x) for x in obj['categories'][:]];counts=collections.Counter(cats[int(i)] for i in obj['codes'][:])
write('source_subgroup_counts.csv',['cell_type','n_cells'],[{'cell_type':k,'n_cells':v} for k,v in sorted(counts.items())])
with open(ROOT+'/logs/python_metadata_environment.json','w') as f:json.dump({'python':sys.version,'h5py':h5py.__version__,'executable':sys.executable,'conda_env':os.getenv('CONDA_DEFAULT_ENV')},f,indent=2)
