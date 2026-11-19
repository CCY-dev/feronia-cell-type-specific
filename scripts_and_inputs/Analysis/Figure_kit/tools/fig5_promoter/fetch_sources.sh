#!/usr/bin/env bash
# Record of the genome-browser fetches used for this panel (run once on 2026-10-04 UTC;
# see sources/fetch_log.tsv for the exact times and HTTP codes). Re-running refetches live
# database records into sources_refetch/ beside this script (never over the frozen sources/), so a later
# Ensembl/NCBI/TAIR release can be compared with diff before anything is replaced.
#   Ensembl Plants REST (TAIR10 assembly, Araport11 genebuild): https://rest.ensembl.org
#     browser view: https://plants.ensembl.org/Arabidopsis_thaliana/Gene/Summary?g=AT5G43890
#   NCBI RefSeq (NC_003076.8, NM_123756.2, Gene 834411): https://eutils.ncbi.nlm.nih.gov
#     browser view: https://www.ncbi.nlm.nih.gov/gene/834411
#   TAIR10 gene models (TAIR download API): https://www.arabidopsis.org
#     locus page: https://www.arabidopsis.org/locus?name=AT5G43890
set -euo pipefail
OUT="$(cd "$(dirname "$0")" && pwd)/sources_refetch"   # kit copy: beside this script
mkdir -p "$OUT"; cd "$OUT"
LOG=fetch_log.tsv; printf "utc_time\turl\thttp_code\tfile\n" > $LOG
get () { local code; code=$(curl -s -L -m 300 -H "Accept: $3" -o "$2" -w "%{http_code}" "$1")
  printf "%s\t%s\t%s\t%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$1" "$code" "$2" >> $LOG; echo "$code $2"; }
E=https://rest.ensembl.org; N=https://eutils.ncbi.nlm.nih.gov/entrez/eutils
get "$E/lookup/id/AT5G43890?species=arabidopsis_thaliana;expand=1;utr=1;content-type=application/json" ensembl_lookup_AT5G43890.json application/json
get "$E/info/assembly/arabidopsis_thaliana?content-type=application/json" ensembl_assembly_arabidopsis_thaliana.json application/json
get "$E/info/software?content-type=application/json" ensembl_info_software.json application/json
get "$E/info/eg_version?content-type=application/json" ensembl_info_eg_version.json application/json
get "$E/sequence/region/arabidopsis_thaliana/5:17648376..17653875:1?content-type=text/x-fasta" ensembl_chr5_17648376-17653875_plus.fa text/x-fasta
get "$E/overlap/region/arabidopsis_thaliana/5:17646376-17655875?feature=gene;feature=transcript;feature=exon;feature=cds;content-type=application/json" ensembl_overlap_chr5_17646376-17655875.json application/json
get "$E/sequence/region/arabidopsis_thaliana/5:17650129..17652142:-1?content-type=text/x-fasta" ensembl_chr5_17650129-17652142_minus.fa text/x-fasta
get "$E/sequence/region/arabidopsis_thaliana/5:17651949..17651993:-1?content-type=text/x-fasta" ensembl_chr5_17651949-17651993_minus_probe.fa text/x-fasta
get "$N/esearch.fcgi?db=gene&term=AT5G43890%5BGene%20Name%5D%20AND%20Arabidopsis%20thaliana%5BOrganism%5D&retmode=json" ncbi_gene_esearch_AT5G43890.json application/json
get "$N/esummary.fcgi?db=gene&id=834411&retmode=json" ncbi_gene_esummary_834411.json application/json
get "$N/efetch.fcgi?db=nuccore&id=NC_003076.8&seq_start=17648376&seq_stop=17653875&strand=1&rettype=gbwithparts&retmode=text" ncbi_NC_003076.8_17648376-17653875_plus.gb text/plain
sleep 1
get "$N/efetch.fcgi?db=nuccore&id=NC_003076.8&seq_start=17648376&seq_stop=17653875&strand=1&rettype=fasta&retmode=text" ncbi_NC_003076.8_17648376-17653875_plus.fa text/plain
get "$N/efetch.fcgi?db=nuccore&id=NC_003076.8&seq_start=17650129&seq_stop=17652142&strand=2&rettype=fasta&retmode=text" ncbi_NC_003076.8_17650129-17652142_minus.fa text/plain
TMP=$(mktemp)
get "https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff" "$TMP" application/octet-stream
shasum -a 256 "$TMP" | awk -v n="$(wc -c < "$TMP")" '{print $1"  TAIR10_GFF3_genes.gff ("n" bytes)"}' > tair10_gff_fullfile_sha256.txt
grep "AT5G43890" "$TMP" > TAIR10_GFF3_genes_AT5G43890.gff; rm -f "$TMP"
echo "Compare with: diff -r tools/fig5_promoter/sources tools/fig5_promoter/sources_refetch"
