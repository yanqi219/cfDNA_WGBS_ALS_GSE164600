## STAR Protocol (verbatim) end-motif pipeline

This folder vendors the **original** bash + R scripts from the upstream repository
`dxl668/PHTScfDNA` (as referenced by the STAR Protocols paper) and provides a thin
wrapper to run them without editing the vendored files.

- **Protocol paper**: [Liu et al., STAR Protocols 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11489062/)
- **Upstream code**: `https://github.com/dxl668/PHTScfDNA`

### What is included

- `upstream/scripts/`: scripts copied **verbatim** from `dxl668/PHTScfDNA/scripts/`:
  - `01-filter_bam.sh`
  - `02-bamtobed.sh`
  - `03-filter_frags_gc.sh`
  - `04-bins5mb.sh`
  - `05-filter_bedpe.sh`
  - `06-endmotif_bed.sh`
  - `07-motif_fasta.sh`
  - `20-motif_gc.sh` + `20-motif_gc.R`
  - `21-combine_motif.sh` + `21-combine_motif.R`

- `files/get_hg19_resources.sh`: downloads the protocol resources into `files/`
  (hg19 FASTA, ENCODE blacklist, gaps BED, and bins5mb filtered BED).

- `run_endmotif_pipeline.sh`: **wrapper** that:
  - creates a working directory under `work/`
  - copies `upstream/scripts/` into `work/scripts/`
  - substitutes the placeholder `PROJECT_DIR=/user/project` **in the working copy only**
    (the vendored upstream scripts remain unchanged)
  - runs the STAR Protocol steps up to and including `21-combine_motif.R`

### Requirements

This is the STAR Protocol’s original tooling, so you need:

- `samtools` (v1.16.1 in protocol)
- `bedtools` (v2.29.0 in protocol)
- `R` (v4.2.3 in protocol) with packages: `tidyverse`, `ggplot2`

Note: the upstream scripts use `module load ...` and `#SBATCH` headers (HPC/Slurm).
They can still be executed on a workstation **if** `samtools`, `bedtools`, and `Rscript`
are on `PATH` (the `module` commands may need to be available or ignored by your shell
environment).

### Resource download (hg19 + blacklist + gaps + bins)

From the repository root:

```bash
bash "cfDNA_WGBS_ALS_GSE164600/star_protocol_phtscfdna/files/get_hg19_resources.sh"
```

This writes into:

`cfDNA_WGBS_ALS_GSE164600/star_protocol_phtscfdna/work/files/`

### Running the STAR Protocol end-motif preprocessing (steps 01–07, 20–21)

1. Put your BAMs here (or adjust after creating the work dir):

`cfDNA_WGBS_ALS_GSE164600/star_protocol_phtscfdna/work/bam_files/*.bam`

2. Run:

```bash
bash "cfDNA_WGBS_ALS_GSE164600/star_protocol_phtscfdna/run_endmotif_pipeline.sh"
```

Outputs will be under:

`cfDNA_WGBS_ALS_GSE164600/star_protocol_phtscfdna/work/`

matching the upstream script directory names:

- `01-filter_bam/`
- `02-bamtobed/`
- `03-filter_frags_gc/`
- `04-bins5mb/`
- `05-filter_bedpe/`
- `06-endmotif_bed/`
- `07-motif_fasta/`
- `08-motif_merge/`
- `20-motif_gc/`
- `21-combine_motif/`

### Scope (per your request)

This wrapper runs everything **before** the protocol’s “End motif analysis and visualization”
step (i.e., it stops after generating `endmotif_4bp_gc_summary.rds`).

