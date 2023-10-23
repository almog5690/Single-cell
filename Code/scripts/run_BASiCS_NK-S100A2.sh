<<<<<<< Updated upstream
#!/bin/sh

module load R4

#SBATCH --ntasks=2
#SBATCH --mem=16G
#!/usr/bin/env Rscript

Rscript --vanilla ../run_BASiCS_commandline.R Blood_SC Blood NK-S100A2
=======
#!/bin/sh

module load R4

#SBATCH --ntasks=2
#SBATCH --mem=16G
#!/usr/bin/env Rscript

Rscript --vanilla ../run_BASiCS_commandline.R Blood_SC Blood NK-S100A2
>>>>>>> Stashed changes
