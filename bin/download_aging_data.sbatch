#!/usr/bin/env bash
#SBATCH --job-name=download_aging
#SBATCH --time=12:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/sci/labs/orzuk/orzuk/projects/SingleCell/Data/_download_logs/%x_%j.out
#SBATCH --error=/sci/labs/orzuk/orzuk/projects/SingleCell/Data/_download_logs/%x_%j.err

set -euo pipefail
mkdir -p /sci/labs/orzuk/orzuk/projects/SingleCell/Data/_download_logs

DATASET="${1:-}"
shift || true
if [[ -z "$DATASET" ]]; then
  echo "Usage: sbatch download_aging_data.sbatch <dataset> [--base DIR] [--no-extract] [--ts-max-files N]" >&2
  exit 1
fi

# Prefer local copy if submitting from the repo; otherwise fallback to the repo path.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -x "${SCRIPT_DIR}/download_aging_data.sh" ]]; then
  "${SCRIPT_DIR}/download_aging_data.sh" "$DATASET" "$@"
else
  /sci/labs/orzuk/orzuk/github/Single-cell/bin/download_aging_data.sh "$DATASET" "$@"
fi
