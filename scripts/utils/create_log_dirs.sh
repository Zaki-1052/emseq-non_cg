#!/bin/bash
# scripts/utils/create_log_dirs.sh
#
# Create the log directory tree that the SLURM wrappers write to, and place a
# .gitkeep in each leaf so git tracks the empty directories. Run once from
# anywhere; the script finds the repo root from its own location.
#
# Usage:
#   bash scripts/utils/create_log_dirs.sh

set -euo pipefail

CODEDIR="$(cd "$(dirname "$0")/../.." && pwd)"

dirs=(
    logs
    logs/02b_features
    logs/sections
    logs/sections/10_chromatin
    logs/sections/20_chip_integration
    logs/sections/30_hic
    logs/sections/40_permutation
    logs/sections/50_features
    logs/sections/60_mecp2
    logs/sections/70_neuronal
    logs/sections/80_cross_modality
)

for d in "${dirs[@]}"; do
    full="$CODEDIR/$d"
    mkdir -p "$full"
    touch "$full/.gitkeep"
    echo "  $d/"
done

echo ""
echo "Log directories ready (with .gitkeep files)."
