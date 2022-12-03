# Single-cell
Single cell RNA analysis about paper “Determinant of Cell-to-Cell Transcriptional Variability During
Mammalian Aging”

## Workflow
After obtaining the Tabula Maris senis raw (or processed) data 
we first perform data pre-processing (filtering of cells and genes, normalization, atc.).
Next, we analyze the correlation between mean expression and selection for bot old and young.
Before we can examine the correlation between over-dispersion and selection and transcript length we first have 
estimate the over-dispersion for each gene in each cell using BASiCS (NOTE: BASiCS analysis can take alot of time).
After the BASics analysis we can perform our over-dispersion analysis.
