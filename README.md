# cfDNA WGBS Analysis - ALS vs Control (GSE164600)

Analysis of cell-free DNA whole-genome bisulfite sequencing data to identify markers distinguishing ALS from healthy controls.

## Data

- Source: GSE164600 (downsampled to chr21)
- Samples: 12 ALS, 10 Control
- Format: BAM files (Bismark aligned)

## Project Structure

```
cfDNA_WGBS_ALS_GSE164600/
├── data/
│   ├── raw/           # Original BAM files (not tracked by git)
│   └── processed/     # Processed data outputs
├── scripts/           # Analysis scripts (numbered for order)
│   └── 00_setup.R     # Package installation and setup
├── R/                 # Reusable functions
├── results/
│   ├── figures/       # Plots and visualizations
│   └── tables/        # Summary tables and metrics
├── renv/              # Package management (auto-generated)
└── renv.lock          # Package versions lock file
```

## Setup

1. Open the project in RStudio (or set working directory)
2. Run the setup script:
   ```r
   source("scripts/00_setup.R")
   ```
   This will initialize renv and install required packages.

## Analysis Pipeline

Scripts are numbered in execution order:

- `00_setup.R` - Install packages
- `01_qc_processing.R` - Comprehensive QC and initial statistics
- `02_analysis.R` - Statistical analysis and visualization
- `03_classification.R` - ALS vs Control classification

## Key Analyses

### 1. QC & Processing (`01_qc_processing.R`)

Generates publication-quality figures covering:

| Analysis | Description | cfDNA Relevance |
|----------|-------------|-----------------|
| **Fragment Length** | Insert size distribution, 167bp peak, 10.4bp periodicity | Nucleosome footprint |
| **End Motifs** | 4-mer frequencies at 5' fragment ends | Tissue-specific cleavage patterns |
| **Methylation QC** | CpG methylation levels | Disease marker |
| **Bisulfite Conversion** | CHH methylation (should be <1%) | Data quality |
| **Sample Metrics** | Read counts, MAPQ, fragment ratios | Technical quality |

**Output figures:**
- `fig1_fragment_length.pdf` - Fragment size distribution
- `fig2_end_motifs.pdf` - End motif analysis
- `fig3_methylation_qc.pdf` - Methylation quality control
- `fig4_sample_qc_overview.pdf` - Sample quality summary
- `fig5_statistical_comparison.pdf` - ALS vs Control comparison

### 2. Downstream Analysis

- **Statistical Analysis** - Differential features between groups
- **Classification** - Binary classification with precision, sensitivity, F1

## Notes

- Data is bisulfite-treated (affects motif analysis)
- Methylation status in XM tag (Bismark convention)
- Consider scalability for full-size samples (~50GB per BAM)
