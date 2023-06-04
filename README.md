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


## Workflow
First data pre-processing is performed on read counts data files (filtering of cells and genes, normalization, atc.), by running the commands in the files 
'data pre-processing facs.R' and 'data pre-processing droplet.R'. Both files contaiins a for loop which perform the pre-processing to every tissue individually.

Next, analysis of the correlation between mean expression and selection for both old and young groups is performed by running the files 'Mean expression facs.R' and 'Mean expression droplet.R'.

Estimation of over-dispersion for each gene in each cell using BASiCS (NOTE: the BASiCS analysis can take a lot of time: about 1 week for the entire dataset on a single laptop).

After running the BASiCS analysis we can perform our over-dispersion analysis, by running the files commands in the files 'Over dispersion facs.R' and 'Over dispersion droplet.R'. 

Running the entire pipeline produces figures showing the correlation between the different covaraites (determinents) and gene-expression over-dispersion and mean expression. Run the file 'scRNA_seq_variability_determinants_analysis.R' to run the pipeline. 


# Contact: 
This repository was developed by Almog Yair. For questions about the package or analysis, please email Almog Yair (almog.yair@mail.huji.ac.il) or Or Zuk (or.zuk@mail.huji.ac.il). 

# Ref
The analysis performed using this code is described in the paper:
“Determinant of Cell-to-Cell Transcriptional Variability During Mammalian Aging”, A. Yair, T. Landsberger and O. Zuk. 
