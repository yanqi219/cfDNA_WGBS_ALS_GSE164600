#!/usr/bin/env bash
set -euo pipefail

# Run STAR Protocol end-motif preprocessing (verbatim upstream scripts) up to:
#   - 21-combine_motif.R → endmotif_4bp_gc_summary.rds
#
# This wrapper does NOT modify vendored upstream scripts.
# It creates a working copy under:
#   cfDNA_WGBS_ALS_GSE164600/star_protocol_phtscfdna/work/
# and substitutes PROJECT_DIR=/user/project → PROJECT_DIR=<absolute_work_dir>
# inside the working copy only.
#
# Protocol: https://pmc.ncbi.nlm.nih.gov/articles/PMC11489062/
# Upstream: https://github.com/dxl668/PHTScfDNA

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UPSTREAM_DIR="${ROOT_DIR}/upstream"
WORK_DIR="${ROOT_DIR}/work"

mkdir -p "${WORK_DIR}"

if [[ ! -d "${UPSTREAM_DIR}/scripts" ]]; then
  echo "ERROR: missing vendored scripts at: ${UPSTREAM_DIR}/scripts" >&2
  exit 1
fi

echo "[STAR wrapper] work dir: ${WORK_DIR}"

# Prepare working copy of scripts (so upstream scripts remain unchanged)
mkdir -p "${WORK_DIR}/scripts"
cp -f "${UPSTREAM_DIR}/scripts/"* "${WORK_DIR}/scripts/"

# Ensure work contains expected directories (scripts will mkdir too; this is just clarity)
mkdir -p "${WORK_DIR}/bam_files"
mkdir -p "${WORK_DIR}/files"

# Substitute PROJECT_DIR placeholder in working copy only
for f in "${WORK_DIR}/scripts/"*.sh; do
  # Only replace exact placeholder line used by upstream scripts.
  # Keep the rest of the file identical.
  perl -0777 -pe "s/^PROJECT_DIR=\\/user\\/project\\s*\$/PROJECT_DIR=${WORK_DIR//\//\\/}/m" -i "$f"
done

# Also patch the R script paths referenced by the .sh scripts by ensuring they live in work/scripts/
# (Upstream scripts already set R_SCRIPT=\"${PROJECT_DIR}/scripts/<name>.R\")
# so copying into work/scripts/ is sufficient.

echo "[STAR wrapper] (optional) download resources into work/files/"
echo "             run: bash \"${ROOT_DIR}/files/get_hg19_resources.sh\""

echo "[STAR wrapper] starting pipeline (01 → 07 → 20 → 21)"

cd "${WORK_DIR}"

# Run in protocol order
bash "scripts/01-filter_bam.sh"
bash "scripts/02-bamtobed.sh"
bash "scripts/03-filter_frags_gc.sh"
bash "scripts/04-bins5mb.sh"
bash "scripts/05-filter_bedpe.sh"
bash "scripts/06-endmotif_bed.sh"
bash "scripts/07-motif_fasta.sh"
bash "scripts/20-motif_gc.sh"
bash "scripts/21-combine_motif.sh"

echo "[STAR wrapper] done"
echo "[STAR wrapper] combined output should be at:"
echo "  ${WORK_DIR}/21-combine_motif/endmotif_4bp_gc_summary.rds"

