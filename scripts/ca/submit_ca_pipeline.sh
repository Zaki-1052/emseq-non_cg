#!/bin/bash
# scripts/ca/submit_ca_pipeline.sh
#
# Submits the CA-only mCH analysis pipeline as a dependency chain.
# Each step waits for the previous one to complete successfully.
#
# Usage: bash scripts/submit_ca_pipeline.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "=== CA-Only mCH Pipeline Submission ==="
echo "Script dir: $SCRIPT_DIR"
echo ""

mkdir -p logs

JOB_A=$(sbatch --parsable "${SCRIPT_DIR}/01_ca_filter.sb")
echo "Step 01: CA-filter (array 0-15)        -> Job $JOB_A"

JOB_B=$(sbatch --parsable --dependency=afterok:${JOB_A} "${SCRIPT_DIR}/02_ca_combine.sb")
echo "Step 02: Combine CA (array 0-7)        -> Job $JOB_B  (after $JOB_A)"

JOB_C=$(sbatch --parsable --dependency=afterok:${JOB_B} "${SCRIPT_DIR}/03_ca_aggregate.sb")
echo "Step 03: Aggregate per gene (array 0-7) -> Job $JOB_C  (after $JOB_B)"

JOB_D=$(sbatch --parsable --dependency=afterok:${JOB_C} "${SCRIPT_DIR}/04_ca_differential.sb")
echo "Step 04: Differential testing           -> Job $JOB_D  (after $JOB_C)"

echo ""
echo "=== All jobs submitted ==="
echo "Monitor: squeue -u \$USER"
echo "Cancel all: scancel $JOB_A $JOB_B $JOB_C $JOB_D"
