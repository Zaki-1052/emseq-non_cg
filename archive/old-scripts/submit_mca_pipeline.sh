#!/bin/bash
# scripts/submit_mca_pipeline.sh
#
# Submit the gene-body mCA analysis pipeline with SLURM dependency chaining.
# Steps run strictly in sequence: the spike-in QC gate in step 2 must pass
# before gene-body aggregation starts.
#
#   1. CA-context filtering        (array 0-15, one task per methylKit file)
#   2. Merge CA stats + spike-in QC
#   3. Gene-body aggregation       (array 0-7, one task per sample)
#   4. DSS differential testing
#
# Step 0 (extract_gene_bodies.py) is run interactively beforehand.
#
# Usage: bash submit_mca_pipeline.sh
#        (run from the emseq/scripts directory on Expanse)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASEDIR="/expanse/lustre/projects/csd940/zalibhai/emseq"
LOGDIR="${BASEDIR}/logs"
GENE_BED="${BASEDIR}/references/gene_bodies.protein_coding.bed"

mkdir -p "$LOGDIR"

echo "========================================="
echo "Gene-body mCA Pipeline"
echo "========================================="
echo "Script dir: $SCRIPT_DIR"
echo "Base dir:   $BASEDIR"
echo "Log dir:    $LOGDIR"
echo ""

if [[ ! -s "$GENE_BED" ]]; then
    echo "ERROR: gene body BED not found or empty: $GENE_BED"
    echo ""
    echo "Run step 0 first:"
    echo "  conda activate emseq_env"
    echo "  python ${SCRIPT_DIR}/extract_gene_bodies.py \\"
    echo "      --gtf /expanse/lustre/projects/csd940/zalibhai/abc/reference/gencode.vM25.annotation.gtf.gz \\"
    echo "      --output ${GENE_BED}"
    exit 1
fi

echo "Gene body BED: $GENE_BED ($(wc -l < "$GENE_BED") genes)"
echo ""

# Step 1: CA-context filtering
STEP1_JID=$(sbatch --parsable "${SCRIPT_DIR}/run_ca_filter.sb")
echo "  Step 1 - CA filter:       job $STEP1_JID (array 0-15)"

# Step 2: Merge stats and gate on spike-in QC
STEP2_JID=$(sbatch --parsable --dependency=afterok:${STEP1_JID} \
    "${SCRIPT_DIR}/run_merge_ca_stats.sb")
echo "  Step 2 - Merge + QC gate: job $STEP2_JID (after $STEP1_JID)"

# Step 3: Gene-body aggregation
STEP3_JID=$(sbatch --parsable --dependency=afterok:${STEP2_JID} \
    "${SCRIPT_DIR}/run_genebody_aggregate.sb")
echo "  Step 3 - Gene-body agg:   job $STEP3_JID (array 0-7, after $STEP2_JID)"

# Step 4: Differential testing
STEP4_JID=$(sbatch --parsable --dependency=afterok:${STEP3_JID} \
    "${SCRIPT_DIR}/run_mca_differential.sb")
echo "  Step 4 - DSS testing:     job $STEP4_JID (after $STEP3_JID)"

echo ""
echo "========================================="
echo "Dependency chain:"
echo "  $STEP1_JID -> $STEP2_JID -> $STEP3_JID -> $STEP4_JID"
echo ""
echo "Monitor: squeue -u \$USER"
echo "Logs:    $LOGDIR"
echo "========================================="
