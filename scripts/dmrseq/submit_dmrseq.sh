#!/bin/bash
# scripts/dmrseq/submit_dmrseq.sh
#
# Submit the dmrseq pipeline with SLURM dependency chaining.
#
#   Phase 1:  Split samples by chromosome          (single job)
#   Phase 2a: Local DMR discovery per chromosome    (array 0-20, after Phase 1)
#   Phase 2b: Block detection per chromosome        (array 0-20, after Phase 1)
#   Phase 3a: Gather local DMR results              (after Phase 2a)
#   Phase 3b: Gather block results                  (after Phase 2b)
#
# Usage: bash submit_dmrseq.sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASEDIR="/expanse/lustre/projects/csd940/zalibhai/emseq"
CH_DIR="${BASEDIR}/combined_ch"

echo "========================================="
echo "dmrseq Pipeline"
echo "========================================="
echo "Script dir: $SCRIPT_DIR"
echo "Base dir:   $BASEDIR"
echo ""

# --- Validate inputs ---

SAMPLES=(ctrl_M1 ctrl_M2 ctrl_F1 ctrl_F2 mut_M1 mut_M2 mut_F1 mut_F2)
MISSING=0
for s in "${SAMPLES[@]}"; do
    f="${CH_DIR}/${s}_CH.methylKit.gz"
    if [[ ! -s "$f" ]]; then
        echo "ERROR: Combined CH file not found: $f"
        MISSING=1
    fi
done
if [[ $MISSING -eq 1 ]]; then
    echo ""
    echo "Fix missing files before running."
    exit 1
fi
echo "All 8 combined CH files verified."
echo ""

# --- Create output directories ---

mkdir -p "${BASEDIR}/results/dmrseq" \
         "${BASEDIR}/.dmrseq_splits" \
         "${BASEDIR}/logs/dmrseq"

# --- Submit jobs ---

SPLIT_JID=$(sbatch --parsable "${SCRIPT_DIR}/01_split_samples.sb")
echo "  Phase 1  - Split samples:      job $SPLIT_JID"

LOCAL_JID=$(sbatch --parsable --dependency=afterok:${SPLIT_JID} \
    "${SCRIPT_DIR}/02_dmrseq_local.sb")
echo "  Phase 2a - Local DMRs (array): job $LOCAL_JID (after $SPLIT_JID)"

BLOCK_JID=$(sbatch --parsable --dependency=afterok:${SPLIT_JID} \
    "${SCRIPT_DIR}/02_dmrseq_blocks.sb")
echo "  Phase 2b - Blocks (array):     job $BLOCK_JID (after $SPLIT_JID)"

GATHER_L_JID=$(sbatch --parsable --dependency=afterok:${LOCAL_JID} \
    "${SCRIPT_DIR}/03_gather_local.sb")
echo "  Phase 3a - Gather local:       job $GATHER_L_JID (after $LOCAL_JID)"

GATHER_B_JID=$(sbatch --parsable --dependency=afterok:${BLOCK_JID} \
    "${SCRIPT_DIR}/03_gather_blocks.sb")
echo "  Phase 3b - Gather blocks:      job $GATHER_B_JID (after $BLOCK_JID)"

echo ""
echo "========================================="
echo "Dependency graph:"
echo ""
echo "  Phase 1 ($SPLIT_JID)"
echo "    ├── Phase 2a ($LOCAL_JID, array 0-20)"
echo "    │     └── Phase 3a ($GATHER_L_JID)"
echo "    └── Phase 2b ($BLOCK_JID, array 0-20)"
echo "          └── Phase 3b ($GATHER_B_JID)"
echo ""
echo "Monitor: squeue -u \$USER"
echo "Logs:    ${BASEDIR}/logs/dmrseq/"
echo "========================================="
