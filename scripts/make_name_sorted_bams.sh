#!/usr/bin/env bash
set -euo pipefail

# Create read-name-sorted BAMs for bismark_methylation_extractor (paired-end).
#
# Why: bismark_methylation_extractor expects paired-end BAMs to be read-name sorted
# (coordinate-sorted BAMs trigger an error).
#
# Default input:
#   /Users/qyan/Desktop/bam/*.bam
# Default output:
#   cfDNA_WGBS_ALS_GSE164600/data/processed/bam_name_sorted/*.name_sorted.bam
#
# You can override paths via environment variables:
#   IN_DIR=/path/to/bams OUT_DIR=/path/to/output THREADS=8 bash make_name_sorted_bams.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PROJ_DIR="${ROOT_DIR}/cfDNA_WGBS_ALS_GSE164600"

IN_DIR_DEFAULT="/Users/qyan/Desktop/bam"
OUT_DIR_DEFAULT="${PROJ_DIR}/data/processed/bam_name_sorted"

IN_DIR="${IN_DIR:-$IN_DIR_DEFAULT}"
OUT_DIR="${OUT_DIR:-$OUT_DIR_DEFAULT}"

THREADS="${THREADS:-4}"

if ! command -v samtools >/dev/null 2>&1; then
  echo "ERROR: samtools not found in PATH." >&2
  exit 127
fi

if [[ ! -d "${IN_DIR}" ]]; then
  echo "ERROR: IN_DIR not found: ${IN_DIR}" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"

shopt -s nullglob
bams=( "${IN_DIR}"/*.bam )
shopt -u nullglob

if [[ ${#bams[@]} -eq 0 ]]; then
  echo "ERROR: No BAM files found in: ${IN_DIR}" >&2
  exit 1
fi

echo "Found ${#bams[@]} BAM(s)"
echo "IN : ${IN_DIR}"
echo "OUT: ${OUT_DIR}"
echo "THREADS: ${THREADS}"

for bam in "${bams[@]}"; do
  base="$(basename "${bam}")"
  sample="${base%.bam}"
  out_bam="${OUT_DIR}/${sample}.name_sorted.bam"

  if [[ -f "${out_bam}" ]]; then
    echo "Skipping (exists): ${out_bam}"
    continue
  fi

  echo ""
  echo "==> ${sample}"
  echo "IN : ${bam}"
  echo "OUT: ${out_bam}"

  # Name sort (required for paired-end methylation extraction)
  samtools sort -n -@ "${THREADS}" -o "${out_bam}" "${bam}"
done

echo ""
echo "Done. Name-sorted BAMs are under: ${OUT_DIR}"

