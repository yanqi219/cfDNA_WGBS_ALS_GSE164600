#!/usr/bin/env bash
set -euo pipefail

# Download STAR Protocol resources into work/files (hg19-based).
# This follows the protocol’s “Download resources and supporting files” section:
# https://pmc.ncbi.nlm.nih.gov/articles/PMC11489062/
#
# Output directory:
#   cfDNA_WGBS_ALS_GSE164600/star_protocol_phtscfdna/work/files

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK_DIR="${ROOT_DIR}/work"
FILES_DIR="${WORK_DIR}/files"
mkdir -p "${FILES_DIR}"

cd "${FILES_DIR}"

echo "[resources] writing to: ${FILES_DIR}"

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: missing required command: $1" >&2
    exit 1
  }
}

# Prefer curl; allow wget.
if command -v curl >/dev/null 2>&1; then
  DL="curl -L -o"
elif command -v wget >/dev/null 2>&1; then
  DL="wget -O"
else
  echo "ERROR: need curl or wget to download resources." >&2
  exit 1
fi

need_cmd gunzip
need_cmd awk

# (a) Human reference genome hg19
if [[ ! -f hg19.fa ]]; then
  echo "[resources] downloading hg19.fa.gz"
  ${DL} hg19.fa.gz "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz"
  gunzip -f hg19.fa.gz
else
  echo "[resources] hg19.fa exists; skipping"
fi

# (b) ENCODE blacklist (hg19)
if [[ ! -f hg19-blacklist.v2.bed ]]; then
  echo "[resources] downloading hg19-blacklist.v2.bed.gz"
  ${DL} hg19-blacklist.v2.bed.gz "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz"
  gunzip -f hg19-blacklist.v2.bed.gz
else
  echo "[resources] hg19-blacklist.v2.bed exists; skipping"
fi

# (c) Genomic gaps for hg19 → gaps.hg19.bed (protocol awk reformat)
if [[ ! -f gaps.hg19.bed ]]; then
  echo "[resources] downloading gap.txt.gz"
  ${DL} gap.txt.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz"
  gunzip -f gap.txt.gz
  # Protocol transformation:
  # awk 'BEGIN {OFS="\t"} { $1=""; print $0 }' gap.txt | awk 'BEGIN {OFS="\t"} {$1=$1; print}' > gaps.hg19.bed
  awk 'BEGIN {OFS="\t"} { $1=""; print $0 }' gap.txt | awk 'BEGIN {OFS="\t"} {$1=$1; print}' > gaps.hg19.bed
else
  echo "[resources] gaps.hg19.bed exists; skipping"
fi

# (d) 5-Mb bins with mappability + GC content → bins5mb_filtered.bed
# Protocol:
# wget bins_5mb.csv
# awk -F, 'NR==1 || ($5 >= 0.3 && $6 >= 0.9) {print $1,$2,$3,$4}' OFS='\t' bins_5mb.csv > bins5mb_filtered.bed
if [[ ! -f bins_5mb.csv ]]; then
  echo "[resources] downloading bins_5mb.csv"
  ${DL} bins_5mb.csv "https://raw.githubusercontent.com/cancer-genomics/reproduce_lucas_wflow/master/code/preprocessing/bins_5mb.csv"
fi

if [[ ! -f bins5mb_filtered.bed ]]; then
  echo "[resources] creating bins5mb_filtered.bed"
  awk -F, 'NR==1 || ($5 >= 0.3 && $6 >= 0.9) {print $1,$2,$3,$4}' OFS='\t' bins_5mb.csv > bins5mb_filtered.bed
else
  echo "[resources] bins5mb_filtered.bed exists; skipping"
fi

echo "[resources] done"

