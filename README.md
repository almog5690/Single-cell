# DOSCEVIA - Determinants Of Single-Cell Expression Variability In Aging 
Analysis of gene-expression variablity and its determinants in single-cell RNA-seq (scRNA-seq) data. 


## Installation
Clone the repository into a directory of your choice (denoted here '/SC/'), and start an R session within this directory. 
You may need to install several R packages that are required for the analysis of scRNA-seq data, such as Seurat, rhdf5, BASiCS. 

# Datasets
The code perfoms analysis of the Tabula Maris senis scRNA-seq data. Please download and put the dataset in a 'data' sub-directory: '/SC/data/'. 

In addition, the code uses gene features data. 
Download the gene selection data from gnomad called 'gnomad.v2.1.1.lof_metrics.by_gene.txt', 
and the gene length data from "mart_export.txt". 


## Code Structure
The main code is in the `Code/` directory:
- `run_scRNA_seq_variability_determinants_analysis.R` - Main entry point for the analysis pipeline
- `scRNA_seq_config.R` - Configuration (paths, parameters, data types)
- `scRNA_seq_utilities.R` - Shared utility functions
- `scRNA_seq_preprocess.R` - Data preprocessing
- `ExpressionRegressionAnalysis.R` - Regression analysis for both mean and overdispersion
- `BASiCSAnalysis.R` - BASiCS statistical modeling
- `MeanFigures.R`, `OverdispersionFigures.R` - Figure generation

**Note:** The `Code/old_src/` directory contains deprecated files that are no longer used. These are kept for reference only - do not use or modify them.

## Environment Configuration
The code supports multiple running environments. Set the `run.env` variable before sourcing `scRNA_seq_config.R`:

| Environment | Description |
|-------------|-------------|
| `or.windows` | Or Zuk, Windows local (default) |
| `or.wsl` | Or Zuk, WSL (Windows Subsystem for Linux) |
| `or.cluster` | Or Zuk, Unix cluster |
| `almog.windows` | Almog Yair, Windows local |

Example usage:
```r
run.env <- "or.wsl"  # Set before sourcing config
source("scRNA_seq_config.R")
```

## Workflow
First data pre-processing is performed on read counts data files (filtering of cells and genes, normalization, etc.), by running the preprocessing pipeline.

Next, analysis of the correlation between expression statistics (mean, overdispersion) and gene features (selection, gene length) is performed using `ExpressionRegressionAnalysis.R`.

Estimation of over-dispersion for each gene in each cell uses BASiCS (NOTE: the BASiCS analysis can take a lot of time: about 1 week for the entire dataset on a single laptop).

Running the entire pipeline produces figures showing the correlation between the different covariates (determinants) and gene-expression over-dispersion and mean expression. Run the file `run_scRNA_seq_variability_determinants_analysis.R` to run the pipeline. 


# Contact: 
This repository was developed by Almog Yair. For questions about the package or analysis, please email Almog Yair (almog.yair@mail.huji.ac.il) or Or Zuk (or.zuk@mail.huji.ac.il). 

# Ref
The analysis performed using this code is described in the paper:
“Determinant of Cell-to-Cell Transcriptional Variability During Mammalian Aging”, A. Yair, T. Landsberger and O. Zuk. 
