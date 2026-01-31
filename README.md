# cfDNA WGBS Analysis: ALS vs Control (GSE164600)

Analysis of cell-free DNA whole-genome bisulfite sequencing data to identify epigenetic markers distinguishing ALS patients from healthy controls. This pipeline performs cfDNA characterization including quality control, insert size analysis, end motif profiling, differential methylation analysis, and classification.

## Data

- **Source:** GEO GSE164600 (downsampled to chromosome 21)
- **Samples:** 22 total (12 ALS, 10 Control)
- **Format:** Bismark-aligned BAM files

## Project Structure

```
cfDNA_WGBS_ALS_GSE164600/
├── data/
│   ├── sample_metadata.csv
│   ├── raw/                          # Original BAM files
│   └── processed/
│       ├── fragments/                # Fragment BED files
│       ├── methylation_extractor/    # Bismark methylation calls
│       └── methylation/rds/          # Processed methylation objects
├── scripts/
│   ├── 01_qc_processing.R            # QC metrics and visualization
│   ├── 02_fragmentation_analysis.R   # Fragment size analysis
│   ├── 03_end_motif_analysis.R       # 4-mer end motif profiling
│   ├── 03A_end_motif_protocol_preprocess.R
│   ├── 04_DNAm_analysis_DMRseq.R     # Differential methylation
│   └── 06_classification_analysis.R  # ML classification
├── results/
│   ├── figures/                      # Publication-ready figures
│   └── tables/                       # Summary statistics and results
├── cfDNA_analysis_report.qmd         # Quarto analysis report
├── cfDNA_WGBS_analysis_report/       # Rendered HTML report
├── single_sample_cfdna_analysis/     # Single-sample analysis module
├── references.bib
└── nature.csl
```
