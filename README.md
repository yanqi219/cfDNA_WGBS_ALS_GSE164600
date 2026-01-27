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
- `01_extract_features.R` - Extract fragment features from BAM files
- `02_analysis.R` - Statistical analysis and visualization
- `03_classification.R` - ALS vs Control classification

## Key Analyses

1. **Fragment Length Distribution** - Insert size analysis
2. **End Motif Analysis** - 4-mer motif frequencies at fragment ends
3. **Methylation Analysis** - CpG methylation from Bismark XM tag
4. **Classification** - Binary classification with performance metrics

## Notes

- Data is bisulfite-treated (affects motif analysis)
- Methylation status in XM tag (Bismark convention)
- Consider scalability for full-size samples (~50GB per BAM)
