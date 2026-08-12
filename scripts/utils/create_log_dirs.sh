#!/bin/bash
# scripts/utils/create_log_dirs.sh
#
# Create the log directory tree that the SLURM wrappers write to.
# SLURM does not create directories, so missing directories cause logs
# to vanish silently. Run this once before submitting any jobs.
#
# Usage:
#   bash scripts/utils/create_log_dirs.sh

set -euo pipefail

CODEDIR="${EMSEQ_CODE_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"

dirs=(
    "$CODEDIR/logs"
    "$CODEDIR/logs/02b_features"
    "$CODEDIR/logs/sections"
    "$CODEDIR/logs/sections/10_chromatin"
    "$CODEDIR/logs/sections/20_chip_integration"
    "$CODEDIR/logs/sections/30_hic"
    "$CODEDIR/logs/sections/40_permutation"
    "$CODEDIR/logs/sections/50_features"
    "$CODEDIR/logs/sections/60_mecp2"
    "$CODEDIR/logs/sections/70_neuronal"
)

for d in "${dirs[@]}"; do
    mkdir -p "$d"
    echo "  $d"
done

echo ""
echo "Log directories ready."
