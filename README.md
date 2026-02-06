# cd4_tcells_hpv_abseq

This repository contains the scripts related to the analysis of HPV Abseq  CD4+ T cells dataset as part of manuscript titled "_Unique naive CD4+ T cell subsets are associated with containment of high-risk human papillomavirus infection_"

# Overview

This repository provides a complete analysis pipeline for single-cell RNA sequencing and antibody sequencing (AbSeq) data from CD4+ T cells.

Analysis workflow:
- 1_main - initial Seurat objects creation, quality control and integration of RNA and AbSeq assays (scripts are not provided)
- 2_clustering - cell typing analysis with master object (scripts are not provided)
- 3_plotting & de - graphical exploration, differential expression analysis and T-cell fitness (TCF) signature analysis
- 4_pseudotime - trajectory analysis with slinghsot and tradeSeq for pseudotime interference and differential expression along trajectories

All functions are documented with roxygen2-style comments for clarity.

# Repository structure
```{bash}
.
├── README.md
├── scripts
│   ├── 3_plotting_and_de
│   │   ├── 3_abseq_DE.R
│   │   ├── 3_Figure_3B_3C.R
│   │   └── 3_Figure6.R
│   ├── 4_pseudotime
│   │   ├── 4_Figure5A_5D.R
│   │   ├── 4_Figure5B_5C_5E_5F.R
│   │   ├── 4_Figure5G.R
│   │   └── 4_trajectory_analysis.R
│   └── plotting_setup.R
└── utils
    ├── 3_utils.R
    └── 4_utils.R

```

# Dependencies

## CRAN Packages
- **Seurat** (v4.3.0.1) - Single-cell RNA-seq analysis
- **tidyverse** (includes: dplyr, tidyr, readr, purrr, tibble, stringr, forcats, ggplot2)
- **ggpubr** - Data visualisation
- **scCustomize** - Custom Seurat visualisations
- **EnhancedVolcano** - Enhanced volcano plots
- **RColorBrewer** - Colour palettes
- **showtext** - Text rendering
- **readxl** - Read Excel files (only for section 1)
- **ggarchery** - Arrow annotations for plots
- **Matrix** - Sparse and dense matrix classes

## Bioconductor Packages
- **tradeSeq** (v1.20.0)- Trajectory inference and differential expression
- **slingshot** (v2.14.0)- Trajectory inference
- **SingleCellExperiment** - Single-cell data container
- **scran** - Single-cell RNA-seq analysis methods
- **BiocParallel** - Parallel processing for Bioconductor
- **clusterProfiler** (v4.10.1)- Functional enrichment analysis
- **enrichplot** - Visualisation for enrichment results
- **ReactomePA** (v1.46.0)- Reactome pathway analysis
- **org.Hs.eg.db** - Human organism annotation database
- **celldex** - Reference datasets for cell type annotation (only for section 1 & 2)
- **SingleR** - Single-cell recognition (only for section 1)
- **MAST** (v1.28.0)- Model-based Analysis of Single-cell Transcriptomics

# Data Availability

GEO/SRA: TBD
Zenodo: TBD

## Input data

**Raw data files:**
- Sample Tag Calls CSV files: `*Sample_Tag_Calls.csv` (per cartridge/experiment)
- DBEC Molecules Per Cell CSV files: `*DBEC_MolsPerCell.csv` (per cartridge/experiment)

**Section 3 & 4 input file:**
- `cd4_tcells_hpv_abseq.rds` - Master Seurat object (required for scripts 3 and 4)