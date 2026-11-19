#!/usr/bin/env python3
"""Derive verified YUC5 (AT5G43890) promoter coordinates from fetched genome-browser records.

Inputs (all in ../sources, fetched by fetch_sources.sh; see sources/fetch_log.tsv):
  ensembl_lookup_AT5G43890.json              Ensembl Plants REST gene/transcript model (Araport11 on TAIR10)
  ensembl_chr5_17648376-17653875_plus.fa     Ensembl REST genomic sequence, chr5 plus strand
  ncbi_NC_003076.8_17648376-17653875_plus.fa NCBI RefSeq genomic sequence, same window
  ncbi_NC_003076.8_17648376-17653875_plus.gb NCBI RefSeq GenBank features (NM_123756.2 / NP_199202.1)
  TAIR10_GFF3_genes_AT5G43890.gff            TAIR10 gene model lines (TAIR download API)
  ../../../../inputs/probe_sequences.json    frozen probe and mutant sequences (task input)

Outputs (../data): coordinates JSON/TSV, predicted motif table, probe alignment,
FASTA of the selected reference interval, and a GFF3 of the derived features.

Coordinate conventions
  Genomic: TAIR10 chr5 = RefSeq NC_003076.8, 1-based inclusive. YUC5 is on the minus strand.
  Displayed sequence: YUC5 sense (coding) strand, 5'->3' = reverse complement of chr5 plus strand.
  TSS-relative: Araport11 transcript 5' end = +1, the base upstream = -1 (no position 0).
  ATG-relative: A of the start codon = +1, the base upstream = -1 (no position 0).
Predicted DOF motif (legend definition "T/AAAAG or CTTTT/A"): WAAAG read on the sense strand
  (sense-strand site) or CTTTW read on the sense strand (= WAAAG on the antisense strand).
"""
import hashlib
import json
import re
from pathlib import Path

# Figure kit copy (Fig. 5c): sources and the frozen probe sequences live beside this script; outputs go to
# data/fig5/promoter/ in the kit. Run from anywhere:
#     python tools/fig5_promoter/derive_coordinates.py && python tools/fig5_promoter/make_promoter_table.py
HERE = Path(__file__).resolve().parent
SRC = HERE / "sources"
OUT = HERE.parents[1] / "data" / "fig5" / "promoter"
BUNDLE = HERE                          # BUNDLE / "inputs/probe_sequences.json" -> tools/fig5_promoter/inputs/
OUT.mkdir(parents=True, exist_ok=True)

WIN_START, WIN_END = 17648376, 17653875  # fetched window, chr5 plus-strand coordinates
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def fasta(p):
    lines = Path(p).read_text().splitlines()
    return lines[0], "".join(x.strip() for x in lines[1:] if not x.startswith(">")).upper()


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ---------------------------------------------------------------- sequence sources
h_ens, plus_ens = fasta(SRC / "ensembl_chr5_17648376-17653875_plus.fa")
h_ncbi, plus_ncbi = fasta(SRC / "ncbi_NC_003076.8_17648376-17653875_plus.fa")
assert len(plus_ens) == WIN_END - WIN_START + 1, len(plus_ens)
assert plus_ens == plus_ncbi, "Ensembl TAIR10 and NCBI NC_003076.8 sequences differ"
assert h_ens.startswith(">chromosome:TAIR10:5:17648376:17653875:1"), h_ens
PLUS = plus_ens


def plus_base(g):
    return PLUS[g - WIN_START]


def sense_seq(g_hi, g_lo):
    """Sense-strand (minus-strand) sequence 5'->3' from genomic g_hi down to g_lo inclusive."""
    assert g_hi >= g_lo
    return rc(PLUS[g_lo - WIN_START:g_hi - WIN_START + 1])


# ---------------------------------------------------------------- annotation sources
ens = json.loads((SRC / "ensembl_lookup_AT5G43890.json").read_text())
assert ens["id"] == "AT5G43890" and ens["display_name"] == "YUC5" and ens["strand"] == -1
assert ens["assembly_name"] == "TAIR10" and ens["seq_region_name"] == "5"
assert len(ens["Transcript"]) == 1, "expected a single Araport11 transcript"
tx = ens["Transcript"][0]
utr5 = [u for u in tx["UTR"] if u["type"] == "five_prime_utr"]
assert len(utr5) == 1
ENS = dict(tx_id=tx["id"], tx_start=tx["start"], tx_end=tx["end"], cds_lo=tx["Translation"]["start"],
           cds_hi=tx["Translation"]["end"], utr5_lo=utr5[0]["start"], utr5_hi=utr5[0]["end"],
           n_exon=len(tx["Exon"]), protein_len=tx["Translation"]["length"], canonical=tx["is_canonical"])

gb = (SRC / "ncbi_NC_003076.8_17648376-17653875_plus.gb").read_text()
m_mrna = re.search(r"\n {5}mRNA {12}complement\((\d+)\.\.(\d+)\)", gb)
m_cds = re.search(r"\n {5}CDS {13}complement\((\d+)\.\.(\d+)\)", gb)
assert m_mrna and m_cds, "unexpected GenBank layout (multi-exon or plus strand?)"
NCBI = dict(mrna_lo=WIN_START + int(m_mrna.group(1)) - 1, mrna_hi=WIN_START + int(m_mrna.group(2)) - 1,
            cds_lo=WIN_START + int(m_cds.group(1)) - 1, cds_hi=WIN_START + int(m_cds.group(2)) - 1,
            transcript_id=re.search(r'/transcript_id="([^"]+)"', gb).group(1),
            protein_id=re.search(r'/protein_id="([^"]+)"', gb).group(1),
            assembly=re.search(r"Assembly: (\S+)", gb).group(1),
            version=re.search(r"\nVERSION\s+(\S+)", gb).group(1))

tair = [l.split("\t") for l in (SRC / "TAIR10_GFF3_genes_AT5G43890.gff").read_text().splitlines() if l.strip()]
tair_mrna = [f for f in tair if f[2] == "mRNA"]
tair_cds = [f for f in tair if f[2] == "CDS"]
assert len(tair_mrna) == 1 and len(tair_cds) == 1
TAIR10 = dict(mrna_lo=int(tair_mrna[0][3]), mrna_hi=int(tair_mrna[0][4]), strand=tair_mrna[0][6],
              cds_lo=int(tair_cds[0][3]), cds_hi=int(tair_cds[0][4]),
              five_prime_utr=[f for f in tair if f[2] == "five_prime_UTR"])

# Concordance: Ensembl (Araport11) == NCBI RefSeq (TAIR/Araport curated)
assert (ENS["tx_start"], ENS["tx_end"]) == (NCBI["mrna_lo"], NCBI["mrna_hi"])
assert (ENS["cds_lo"], ENS["cds_hi"]) == (NCBI["cds_lo"], NCBI["cds_hi"])
assert (TAIR10["cds_lo"], TAIR10["cds_hi"]) == (ENS["cds_lo"], ENS["cds_hi"]), "TAIR10 CDS differs"

TSS = ENS["tx_end"]           # minus strand: transcript 5' end is the higher coordinate
ATG_A = ENS["cds_hi"]         # first base of the start codon
assert ENS["utr5_hi"] == TSS and ENS["utr5_lo"] == ATG_A + 1
UTR5_LEN = TSS - ATG_A
start_codon = sense_seq(ATG_A, ATG_A - 2)
stop_codon = sense_seq(ENS["cds_lo"] + 2, ENS["cds_lo"])
cds_len = ENS["cds_hi"] - ENS["cds_lo"] + 1
assert start_codon == "ATG", start_codon
assert stop_codon in ("TAA", "TAG", "TGA"), stop_codon
assert cds_len % 3 == 0 and cds_len // 3 - 1 == ENS["protein_len"]


def tss_rel(g):
    d = TSS - g
    return d + 1 if d >= 0 else d


def atg_rel(g):
    d = ATG_A - g
    return d + 1 if d >= 0 else d


def g_from_tss_rel(p):
    return TSS - (p - 1) if p > 0 else TSS - p


# ---------------------------------------------------------------- reference interval
# Legend: "YUC5 gene promoter including 5' UTR (2011bp)". Selected reference interval:
# the 2011 bp immediately upstream of the ATG (ends at the last 5'UTR base). The distal
# endpoint is therefore set by the legend length, not by a construct sequence.
LEGEND_LEN, GRAPHIC_LEN = 2011, 2014
REF_HI = ATG_A + LEGEND_LEN           # distal (5') end, genomic
REF_LO = ATG_A + 1                    # last 5'UTR base
assert REF_HI - REF_LO + 1 == LEGEND_LEN
ref_seq = sense_seq(REF_HI, REF_LO)
ref_atg_seq = sense_seq(REF_HI, ATG_A - 2)
assert len(ref_atg_seq) == GRAPHIC_LEN and ref_atg_seq.endswith("ATG")

# ---------------------------------------------------------------- probe
probe = json.loads((BUNDLE / "inputs/probe_sequences.json").read_text())
WT, MUT = probe["WT"].upper(), probe["mutant"].upper()
assert set(WT) <= set("ACGT") and set(MUT) <= set("ACGT")
win_sense = rc(PLUS)   # sense strand over the whole fetched window, 5'->3'
hits_sense = [m.start() for m in re.finditer("(?=%s)" % WT, win_sense)]
hits_anti = [m.start() for m in re.finditer("(?=%s)" % rc(WT), win_sense)]
assert len(hits_sense) == 1 and not hits_anti, (hits_sense, hits_anti)
j = hits_sense[0]
PROBE_HI = WIN_END - j                 # genomic coordinate of probe 5' base (sense strand)
PROBE_LO = PROBE_HI - len(WT) + 1
assert sense_seq(PROBE_HI, PROBE_LO) == WT
assert REF_LO <= PROBE_LO and PROBE_HI <= REF_HI, "probe outside selected interval"
mut_hits = len(re.findall(MUT, win_sense)) + len(re.findall(rc(MUT), win_sense))
subs = [i for i, (a, b) in enumerate(zip(WT, MUT)) if a != b]
assert len(WT) == len(MUT)
sub_blocks = []
for i in subs:
    if sub_blocks and i == sub_blocks[-1][1] + 1:
        sub_blocks[-1][1] = i
    else:
        sub_blocks.append([i, i])

# ---------------------------------------------------------------- predicted DOF motifs
MOTIF_SENSE = re.compile(r"(?=([AT]AAAG))")
MOTIF_ANTI = re.compile(r"(?=(CTTT[AT]))")    # WAAAG on the antisense strand


def scan(seq_sense, g_first):
    """seq_sense: sense-strand string whose first base is genomic g_first (descending)."""
    sites = []
    for strand, pat in (("sense", MOTIF_SENSE), ("antisense", MOTIF_ANTI)):
        for m in pat.finditer(seq_sense):
            g5 = g_first - m.start()       # sense-strand 5' end of the 5-mer
            g3 = g5 - 4
            sites.append(dict(strand=strand, sense_strand_text=m.group(1),
                              motif_read_5to3=m.group(1) if strand == "sense" else rc(m.group(1)),
                              chr5_start=g3, chr5_end=g5,
                              tss_rel_start=tss_rel(g5), tss_rel_end=tss_rel(g3),
                              atg_rel_start=atg_rel(g5), atg_rel_end=atg_rel(g3)))
    sites.sort(key=lambda s: -s["chr5_end"])
    return sites


sites = scan(ref_atg_seq, REF_HI)     # includes the ATG; no site can contain ATG though
for s in sites:
    s["region"] = "5'UTR" if s["tss_rel_start"] > 0 else ("spans TSS" if s["tss_rel_end"] > 0 else "upstream of TSS")
    s["in_emsa_probe"] = PROBE_LO <= s["chr5_start"] and s["chr5_end"] <= PROBE_HI
    s["destroyed_in_mutant"] = False
    if s["in_emsa_probe"]:
        k0 = PROBE_HI - s["chr5_end"]
        s["destroyed_in_mutant"] = any(k0 <= i <= k0 + 4 for i in subs)
for i, s in enumerate(sites, 1):
    s["site_id"] = "DOF%02d" % i
probe_sites = [s for s in sites if s["in_emsa_probe"]]
mut_left = len(MOTIF_SENSE.findall(MUT)) + len(MOTIF_ANTI.findall(MUT))
aaag_core_probe = [m.start() for m in re.finditer("(?=AAAG)", WT)] + [m.start() for m in re.finditer("(?=CTTT)", WT)]

# ---------------------------------------------------------------- write outputs
coords = {
    "gene": "YUC5", "agi": "AT5G43890", "aliases_ncbi": "YUCCA5, SUPER1 (NCBI Gene 834411)",
    "assembly": "TAIR10 (GCF_000001735.4; chromosome 5 = NC_003076.8)",
    "annotation": "Araport11 (Ensembl Plants release 63 / Ensembl 116 REST; NCBI RefSeq NM_123756.2, curated by TAIR and Araport)",
    "chromosome": "5", "gene_strand": "-",
    "displayed_orientation": "YUC5 sense strand 5'->3' (reverse complement of chr5 plus strand); upstream at left",
    "transcript": ENS["tx_id"], "transcript_chr5": [ENS["tx_start"], ENS["tx_end"]], "exons": ENS["n_exon"],
    "refseq_mrna": NCBI["transcript_id"], "refseq_protein": NCBI["protein_id"],
    "tss_chr5": TSS, "tss_definition": "5' end of Araport11 AT5G43890.1 (= RefSeq NM_123756.2 5' end)",
    "five_prime_utr_chr5": [ATG_A + 1, TSS], "five_prime_utr_length_bp": UTR5_LEN,
    "five_prime_utr_tss_rel": [1, UTR5_LEN],
    "atg_chr5": [ATG_A - 2, ATG_A], "atg_tss_rel": [tss_rel(ATG_A), tss_rel(ATG_A - 2)],
    "cds_chr5": [ENS["cds_lo"], ENS["cds_hi"]], "cds_length_bp": cds_len, "protein_length_aa": ENS["protein_len"],
    "stop_codon": stop_codon,
    "tair10_model": {"mrna_chr5": [TAIR10["mrna_lo"], TAIR10["mrna_hi"]], "cds_chr5": [TAIR10["cds_lo"], TAIR10["cds_hi"]],
                     "five_prime_utr_annotated": bool(TAIR10["five_prime_utr"]),
                     "note": "TAIR10 AT5G43890.1 had no UTRs (mRNA = CDS); the TSS and 5'UTR come from Araport11"},
    "selected_reference_interval": {
        "definition": "2011 bp immediately upstream of the ATG (distal promoter through the last 5'UTR base); length from the legend",
        "chr5": [REF_LO, REF_HI], "strand": "-", "length_bp": LEGEND_LEN,
        "tss_rel": [tss_rel(REF_HI), tss_rel(REF_LO)], "atg_rel": [atg_rel(REF_HI), atg_rel(REF_LO)],
        "with_atg_chr5": [ATG_A - 2, REF_HI], "with_atg_length_bp": GRAPHIC_LEN,
        "status": "reference candidate; assayed construct endpoints not confirmed by any supplied construct sequence"},
    "emsa_probe": {
        "wt_sequence_5to3": WT, "mutant_sequence_5to3": MUT, "length_bp": len(WT), "mutant_length_bp": len(MUT),
        "strand_written": "YUC5 sense strand (chr5 minus strand)",
        "exact_matches_in_5500bp_window": {"sense": len(hits_sense), "antisense": len(hits_anti)},
        "chr5": [PROBE_LO, PROBE_HI], "tss_rel": [tss_rel(PROBE_HI), tss_rel(PROBE_LO)],
        "atg_rel": [atg_rel(PROBE_HI), atg_rel(PROBE_LO)],
        "substituted_positions_1based": [i + 1 for i in subs], "n_substitutions": len(subs),
        "substitution_blocks_1based": [[a + 1, b + 1] for a, b in sub_blocks],
        "substitution_blocks_text": ["%s>%s" % (WT[a:b + 1], MUT[a:b + 1]) for a, b in sub_blocks],
        "predicted_sites_in_probe": [s["site_id"] for s in probe_sites],
        "predicted_sites_remaining_in_mutant": mut_left,
        "mutant_exact_matches_in_window": mut_hits,
        "aaag_or_cttt_cores_in_wt_probe_1based": sorted(i + 1 for i in aaag_core_probe),
        "experimental_status": "region tested by EMSA (DAG1); the predicted sites themselves are motif matches"},
    "predicted_dof_motif": {"definition": "(A/T)AAAG on either strand (legend: T/AAAAG or CTTTT/A)",
                            "n_sites_in_reference_with_atg": len(sites),
                            "n_sense": sum(s["strand"] == "sense" for s in sites),
                            "n_antisense": sum(s["strand"] == "antisense" for s in sites),
                            "n_in_5utr": sum(s["region"] == "5'UTR" for s in sites)},
    "source_files_sha256": {p.name: sha(p) for p in sorted(SRC.iterdir()) if p.is_file()},
}
(OUT / "YUC5_promoter_coordinates.json").write_text(json.dumps(coords, indent=2) + "\n")

cols = ["site_id", "strand", "motif_read_5to3", "sense_strand_text", "chr5_start", "chr5_end", "tss_rel_start",
        "tss_rel_end", "atg_rel_start", "atg_rel_end", "region", "in_emsa_probe", "destroyed_in_mutant"]
with open(OUT / "YUC5_DOF_motifs_predicted.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for s in sites:
        fh.write("\t".join(str(s[c]) for c in cols) + "\n")

feat_rows = [
    ("reference_interval_2011bp", REF_LO, REF_HI), ("reference_interval_with_ATG_2014bp", ATG_A - 2, REF_HI),
    ("five_prime_UTR", ATG_A + 1, TSS), ("TSS", TSS, TSS), ("ATG", ATG_A - 2, ATG_A),
    ("CDS", ENS["cds_lo"], ENS["cds_hi"]), ("transcript_AT5G43890.1", ENS["tx_start"], ENS["tx_end"]),
    ("EMSA_probe", PROBE_LO, PROBE_HI)]
with open(OUT / "YUC5_promoter_features.tsv", "w") as fh:
    fh.write("feature\tchr\tchr5_start\tchr5_end\tstrand\tlength_bp\ttss_rel_5prime\ttss_rel_3prime\tatg_rel_5prime\tatg_rel_3prime\n")
    for name, lo, hi in feat_rows:
        fh.write(f"{name}\t5\t{lo}\t{hi}\t-\t{hi - lo + 1}\t{tss_rel(hi)}\t{tss_rel(lo)}\t{atg_rel(hi)}\t{atg_rel(lo)}\n")

with open(OUT / "YUC5_promoter_features.gff3", "w") as fh:
    fh.write("##gff-version 3\n##sequence-region 5 1 26975502\n# TAIR10 chr5 (NC_003076.8); derived by derive_coordinates.py\n")
    for name, lo, hi in feat_rows:
        fh.write(f"5\tderived\tregion\t{lo}\t{hi}\t.\t-\t.\tID={name};Name={name}\n")
    for s in sites:
        st = "-" if s["strand"] == "sense" else "+"
        fh.write(f"5\tmotif_scan\tTF_binding_site\t{s['chr5_start']}\t{s['chr5_end']}\t.\t{st}\t.\t"
                 f"ID={s['site_id']};Name=WAAAG;Note=predicted DOF motif on YUC5 {s['strand']} strand\n")


def wrap(s, n=60):
    return "\n".join(s[i:i + n] for i in range(0, len(s), n))


(OUT / "YUC5_reference_2011bp_sense.fa").write_text(
    f">YUC5_AT5G43890_reference_2011bp TAIR10 chr5:{REF_LO}-{REF_HI} minus strand, sense 5'->3', "
    f"TSS-relative {tss_rel(REF_HI)}..+{tss_rel(REF_LO)}; ends at last 5'UTR base before ATG\n{wrap(ref_seq)}\n")
(OUT / "YUC5_reference_with_ATG_2014bp_sense.fa").write_text(
    f">YUC5_AT5G43890_reference_with_ATG_2014bp TAIR10 chr5:{ATG_A - 2}-{REF_HI} minus strand, sense 5'->3', "
    f"2011 bp + ATG\n{wrap(ref_atg_seq)}\n")

ruler = "".join("|" if i in subs else " " for i in range(len(WT)))
(OUT / "YUC5_probe_alignment.txt").write_text(
    f"EMSA probe vs mutant, YUC5 sense strand 5'->3' (chr5 minus strand)\n"
    f"WT probe  chr5:{PROBE_LO}-{PROBE_HI}, TSS-relative {tss_rel(PROBE_HI)} to {tss_rel(PROBE_LO)}, "
    f"ATG-relative {atg_rel(PROBE_HI)} to {atg_rel(PROBE_LO)}; {len(WT)} bp; exact single match (sense) in the 5.5-kb window\n\n"
    f"WT      5'-{WT}-3'\n           {ruler}\nMutant  5'-{MUT}-3'\n\n"
    f"Substituted positions (1-based within probe): {[i + 1 for i in subs]}\n"
    f"Predicted (A/T)AAAG sites in WT probe: {[(s['site_id'], s['motif_read_5to3'], s['tss_rel_start'], s['tss_rel_end']) for s in probe_sites]}\n"
    f"Predicted sites remaining in mutant: {mut_left}\n")

print(json.dumps({k: coords[k] for k in ("tss_chr5", "five_prime_utr_chr5", "five_prime_utr_length_bp", "atg_chr5",
                                          "atg_tss_rel", "cds_length_bp", "protein_length_aa", "stop_codon")}))
print(json.dumps(coords["selected_reference_interval"]))
print(json.dumps(coords["emsa_probe"]))
print(json.dumps(coords["predicted_dof_motif"]))
