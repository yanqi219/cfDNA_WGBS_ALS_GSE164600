#!/usr/bin/env bash
set -euo pipefail

# Run Bismark methylation extraction for all BAMs in a folder.
#
# Defaults are set to avoid whitespace-in-path issues by writing outputs to:
#   /Users/qyan/primamente_outputs/methylation_extractor/
#
# Environment overrides:
#   BAM_DIR=... OUT_DIR=... GENOME_FOLDER=... THREADS=8 BISMARK_EXTRACTOR=...

THREADS="${THREADS:-4}"
OUT_DIR="${OUT_DIR:-/Users/qyan/Desktop/methy/methylation_extractor}"

# Prefer name-sorted BAMs produced by scripts/make_name_sorted_bams.sh
if [[ -d "/Users/qyan/primamente_task/cfDNA_WGBS_ALS_GSE164600/data/processed/bam_name_sorted" ]]; then
  DEFAULT_BAM_DIR="/Users/qyan/primamente_task/cfDNA_WGBS_ALS_GSE164600/data/processed/bam_name_sorted"
else
  ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -L)"
  DEFAULT_BAM_DIR="${ROOT_DIR}/cfDNA_WGBS_ALS_GSE164600/data/processed/bam_name_sorted"
fi
BAM_DIR="${BAM_DIR:-$DEFAULT_BAM_DIR}"

BISMARK_EXTRACTOR="${BISMARK_EXTRACTOR:-$HOME/Bismark/bismark_methylation_extractor}"
if [[ ! -x "${BISMARK_EXTRACTOR}" ]]; then
  if command -v bismark_methylation_extractor >/dev/null 2>&1; then
    BISMARK_EXTRACTOR="$(command -v bismark_methylation_extractor)"
  else
    echo "ERROR: bismark_methylation_extractor not found (PATH or \$HOME/Bismark)." >&2
    exit 127
  fi
fi

[[ -d "${BAM_DIR}" ]] || { echo "ERROR: BAM_DIR not found: ${BAM_DIR}" >&2; exit 1; }

mkdir -p "${OUT_DIR}"

shopt -s nullglob
bams=( "${BAM_DIR}"/*.bam )
shopt -u nullglob
[[ ${#bams[@]} -gt 0 ]] || { echo "ERROR: No *.bam found in: ${BAM_DIR}" >&2; exit 1; }

echo "BAM_DIR:        ${BAM_DIR}"
echo "OUT_DIR:        ${OUT_DIR}"
echo "THREADS:        ${THREADS}"
echo "EXTRACTOR:      ${BISMARK_EXTRACTOR}"
echo "BAM count:      ${#bams[@]}"

for bam in "${bams[@]}"; do
  sample="$(basename "${bam%.bam}")"
  out="${OUT_DIR}/${sample}"
  mkdir -p "${out}"

  echo ""
  echo "==> ${sample}"
  (
    cd "${out}"
    "${BISMARK_EXTRACTOR}" \
      --paired-end \
      --no_overlap \
      --bedGraph \
      --counts \
      --gzip \
      --parallel "${THREADS}" \
      "${bam}"
  )
done

echo ""
echo "Done. Outputs under: ${OUT_DIR}"

