<<<<<<< Updated upstream
#!/bin/sh

module load R4

#SBATCH --ntasks=2
#SBATCH --mem=16G
#!/usr/bin/env Rscript

Rscript --vanilla ../run_BASiCS_commandline.R Blood_SC Blood Memory-B-IL32
=======
#!/bin/sh

module load R4

#SBATCH --ntasks=2
#SBATCH --mem=16G
#!/usr/bin/env Rscript

Rscript --vanilla ../run_BASiCS_commandline.R Blood_SC Blood Memory-B-IL32
>>>>>>> Stashed changes
