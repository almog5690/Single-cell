#!/usr/bin/env bash
#SBATCH --job-name=download_aging
#SBATCH --time=12:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/sci/labs/orzuk/orzuk/projects/SingleCell/Data/_download_logs/%x_%j.out
#SBATCH --error=/sci/labs/orzuk/orzuk/projects/SingleCell/Data/_download_logs/%x_%j.err

set -euo pipefail
mkdir -p /sci/labs/orzuk/orzuk/projects/SingleCell/Data/_download_logs

DATASET="${1:-}"
if [[ -z "$DATASET" ]]; then
  echo "Usage: sbatch download_aging_data.sbatch <dataset>" >&2
  exit 1
fi

/sci/labs/orzuk/orzuk/github/Single-cell/bin/download_aging_data.sh "$DATASET"
SBATCH
