module load R4

#SBATCH --ntasks=2
#SBATCH --mem=4G
#!/usr/bin/env Rscript

Rscript --vanilla ../run_BASiCS_commandline.R Blood_SC Blood pDC
