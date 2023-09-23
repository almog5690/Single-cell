module load R4

#SBATCH --time=

#SBATCH -n 

#SBATCH --mem=

#SBATCH --time=10:00:00
#SBATCH --ntasks=2
#SBATCH --mem=1G

#!/usr/bin/env Rscript

# R Script for running BASiCS on large files 
--------------------------------------------

Rscript --vanilla run_BASiCS_commandline.R Blood_SC Blood NK-FCER1G
