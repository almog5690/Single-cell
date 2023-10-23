#!/bin/sh

#SBATCH --time=168:00:00
#SBATCH --ntasks=4
#SBATCH --mem=48G
module load R4
#!/usr/bin/env Rscript

module load R4
Rscript --vanilla ../run_BASiCS_commandline.R Blood_SC Blood CD8-Tem
