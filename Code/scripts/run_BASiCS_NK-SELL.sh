<<<<<<< Updated upstream
#!/bin/bash

#SBATCH --time=168:00:00
#SBATCH --ntasks=4
#SBATCH --mem=48G
module load R4

#!/usr/bin/env Rscript
=======
#!/bin/sh
>>>>>>> Stashed changes

#SBATCH --time=168:00:00
#SBATCH --ntasks=4
#SBATCH --mem=48G
module load R4
<<<<<<< Updated upstream

=======
#!/usr/bin/env Rscript

module load R4
>>>>>>> Stashed changes
Rscript --vanilla ../run_BASiCS_commandline.R Blood_SC Blood NK-SELL
