#!/usr/bin/env python
"""Fig. 5c: write the two tidy tables the figure draws, from the outputs of derive_coordinates.py.

    python tools/fig5_promoter/derive_coordinates.py      # TAIR10/Araport11 sources -> data/fig5/promoter/
    python tools/fig5_promoter/make_promoter_table.py     # -> data/fig5/fig5c_promoter_elements.csv
                                                          #    data/fig5/fig5c_probe_sequences.csv

Coordinates are TSS-relative (Araport11 AT5G43890.1 5' end = +1, no position 0), read on the YUC5 sense strand,
5' (distal promoter) on the left. Nothing here is typed by hand: every number is copied from
data/fig5/promoter/YUC5_promoter_coordinates.json and YUC5_DOF_motifs_predicted.tsv.
"""
import csv
import json
from pathlib import Path

KIT = Path(__file__).resolve().parents[2]
SRC = KIT / 'data' / 'fig5' / 'promoter'
OUT = KIT / 'data' / 'fig5'


def blocks(pos):
    """[12, 13, 14, 15, 16, 29, ...] -> '12-16;29-33' (1-based, inclusive)."""
    out, run = [], []
    for p in sorted(pos):
        if run and p == run[-1] + 1:
            run.append(p)
        else:
            if run:
                out.append(run)
            run = [p]
    if run:
        out.append(run)
    return ';'.join(f'{r[0]}-{r[-1]}' for r in out)


def main():
    c = json.loads((SRC / 'YUC5_promoter_coordinates.json').read_text())
    with open(SRC / 'YUC5_DOF_motifs_predicted.tsv') as fh:
        sites = list(csv.DictReader(fh, delimiter='\t'))
    ref = c['selected_reference_interval']
    pr = c['emsa_probe']
    rows = [
        dict(element='reference_interval', id='pYUC5_2011bp', strand='', tss_rel_start=ref['tss_rel'][0],
             tss_rel_end=ref['tss_rel'][1], note=ref['definition']),
        dict(element='five_prime_utr', id='5UTR', strand='', tss_rel_start=c['five_prime_utr_tss_rel'][0],
             tss_rel_end=c['five_prime_utr_tss_rel'][1], note='Araport11 AT5G43890.1'),
        dict(element='tss', id='TSS', strand='', tss_rel_start=1, tss_rel_end=1, note=c['tss_definition']),
        dict(element='atg', id='ATG', strand='', tss_rel_start=c['atg_tss_rel'][0], tss_rel_end=c['atg_tss_rel'][1],
             note='start codon'),
        dict(element='emsa_probe', id='EMSA_probe', strand='sense', tss_rel_start=pr['tss_rel'][0],
             tss_rel_end=pr['tss_rel'][1], note=f"chr5:{pr['chr5'][0]}-{pr['chr5'][1]}; 45 bp"),
    ]
    for s in sites:
        rows.append(dict(element='dof_site', id=s['site_id'], strand=s['strand'], tss_rel_start=int(s['tss_rel_start']),
                         tss_rel_end=int(s['tss_rel_end']),
                         note=f"{s['motif_read_5to3']} ({s['sense_strand_text']} on the sense strand); {s['region']}"))
    assert sum(r['element'] == 'dof_site' for r in rows) == c['predicted_dof_motif']['n_sites_in_reference_with_atg']
    cols = ['element', 'id', 'strand', 'tss_rel_start', 'tss_rel_end', 'note']
    with open(OUT / 'fig5c_promoter_elements.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    hl = blocks(pr['substituted_positions_1based'])
    seq = [dict(name='probe', sequence_5to3=pr['wt_sequence_5to3'], highlight_1based=hl,
                note='EMSA probe (YUC5 sense strand); highlighted = the two TAAAG sites substituted in the mutant'),
           dict(name='mutant', sequence_5to3=pr['mutant_sequence_5to3'], highlight_1based=hl,
                note='mutant probe; highlighted = substituted bases')]
    with open(OUT / 'fig5c_probe_sequences.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['name', 'sequence_5to3', 'highlight_1based', 'note'], lineterminator='\n')
        w.writeheader()
        w.writerows(seq)
    print(f'wrote {len(rows)} elements and 2 sequences to data/fig5/')


if __name__ == '__main__':
    main()
