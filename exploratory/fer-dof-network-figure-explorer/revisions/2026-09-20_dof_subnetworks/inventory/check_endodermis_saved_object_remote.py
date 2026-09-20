# Read the lab-generated CellOracle object in place. No remote writes or jobs.
import json, pickle, zlib
import h5py
path = '/resnick/groups/NolanLab/cchen/Feronia/analysis/CellOracle/CellOracle_output/20260317_feronia_salinity_timecourse_WT_untreated_fer4_untreated_newIntegratedSeurat3_NoSuspiciousCluster/endodermis_rawRNA_p01_20k.links.celloracle.links'
with h5py.File(path, 'r') as f:
    def read(key):
        return pickle.loads(zlib.decompress(f['&' + key][()].tobytes()))
    links = read('filtered_links')
    score = read('merged_score')
    meta = {'path':path, 'model_method':read('model_method'), 'alpha':read('ALPHA_used')}
    out = []
    for cluster, edges in links.items():
        s = score[score.cluster == cluster]
        counts = edges.groupby('source').size()
        observed = s['degree_out']
        union = counts.index.union(observed.index)
        a = counts.reindex(union).fillna(0)
        b = observed.reindex(union)
        out.append({'cluster':cluster, 'edges':len(edges), 'score_genes':len(s),
                    'edge_sources_absent_from_scores':int(b.isna().sum()),
                    'scored_sources_count_mismatch':int(((a != b) & b.notna()).sum())})
    print(json.dumps({'metadata':meta, 'checks':out}, indent=2, default=str))
