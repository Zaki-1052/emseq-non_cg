#!/bin/bash
# scripts/submit_mch_pipeline.sh
#
# Submit the gene-body mCH analysis pipeline with SLURM dependency chaining.
#
#   1. Combine CHH + CHG            (array 0-7)
#   2a. Gene-body aggregation       (array 0-7, after 1)
#   2b. edgeR differential testing  (after 2a)
#   3. dmrseq DMR discovery         (after 1, parallel with 2a/2b)
#   4. Integration analyses         (after 2b)
#
# Usage: bash submit_mch_pipeline.sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASEDIR="/expanse/lustre/projects/csd940/zalibhai/emseq"
LOGDIR="${BASEDIR}/logs"
MARINER_BASE="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c"

mkdir -p "$LOGDIR" \
         "${BASEDIR}/combined_ch" \
         "${BASEDIR}/results/02_aggregate/aggregated" \
         "${BASEDIR}/results/03_differential" \
         "${BASEDIR}/results/04_integration" \
         "${BASEDIR}/results/dmrseq" \
         "${BASEDIR}/results/plots"

echo "========================================="
echo "Gene-body mCH Pipeline"
echo "========================================="
echo "Script dir: $SCRIPT_DIR"
echo "Base dir:   $BASEDIR"
echo "Log dir:    $LOGDIR"
echo ""

GENE_BED="${BASEDIR}/references/gene_bodies.protein_coding.bed"
K119UB_CONSENSUS="${BASEDIR}/references/K119ub_consensus_v3.bed"
MECP2_CONSENSUS="${BASEDIR}/references/mecp2-peaks/MeCP2_adult_concensus_peakset_Conc4.txt"
MECP2_DIFFBIND="${BASEDIR}/references/mecp2-peaks/260804_MeCP2_adult_diffbind_results_filtered_Conc4.txt"
K119UB_DIFFBIND="${MARINER_BASE}/peaks/diffbind/K119ub_diffbind_results_summit_appended_ap.txt"
NEURONAL_GENES="${MARINER_BASE}/biomodal/downstream/plots/visualizations/tables/72_neuronal_gene_set_go_derived.tsv"
MC_DMR="${MARINER_BASE}/biomodal/downstream/modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260402_191818/DMR_mc_control__mutant_20260402_191818.bed"
HMC_DMR="${MARINER_BASE}/biomodal/downstream/modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260402_191818/DMR_hmc_control__mutant_20260402_191818.bed"

MISSING=0
for f in "$GENE_BED" "$K119UB_CONSENSUS" "$MECP2_CONSENSUS" "$MECP2_DIFFBIND" \
         "$K119UB_DIFFBIND" "$NEURONAL_GENES" "$MC_DMR" "$HMC_DMR"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: Required file not found or empty: $f"
        MISSING=1
    fi
done

SAMPLES=(ctrl_M1 mut_M1 ctrl_F1 mut_F1 ctrl_M2 mut_M2 ctrl_F2 mut_F2)
declare -A EMSEQ_PREFIX
EMSEQ_PREFIX=(
    [ctrl_M1]="260629_Emseq_LV_index_E1_Bap1_Math1_adult_ctrl_M1_S1"
    [mut_M1]="260629_Emseq_LV_index_F1_Bap1_Math1_adult_mut_M1_S2"
    [ctrl_F1]="260629_Emseq_LV_index_G1_Bap1_Math1_adult_ctrl_F1_S3"
    [mut_F1]="260629_Emseq_LV_index_H1_Bap1_Math1_adult_mut_F1_S4"
    [ctrl_M2]="260629_Emseq_LV_index_A2_Bap1_Math1_adult_ctrl_M2_S5"
    [mut_M2]="260629_Emseq_LV_index_B2_Bap1_Math1_adult_mut_M2_S6"
    [ctrl_F2]="260629_Emseq_LV_index_C2_Bap1_Math1_adult_ctrl_F2_S7"
    [mut_F2]="260629_Emseq_LV_index_D2_Bap1_Math1_adult_mut_F2_S8"
)
EXTRACT_DIR="${BASEDIR}/em-seq_output/methylDackelExtracts"
for s in "${SAMPLES[@]}"; do
    for ctx in CHH CHG; do
        f="${EXTRACT_DIR}/${EMSEQ_PREFIX[$s]}_${ctx}.methylKit.gz"
        if [[ ! -s "$f" ]]; then
            echo "ERROR: Input file not found: $f"
            MISSING=1
        fi
    done
done

if [[ $MISSING -eq 1 ]]; then
    echo ""
    echo "Fix missing files before running."
    exit 1
fi

echo "Gene body BED: $GENE_BED ($(wc -l < "$GENE_BED") genes)"
echo "All 16 input methylKit files verified."
echo "All reference files verified."
echo ""

# Step 1 (combine CH) already completed — verify outputs exist
CH_DIR="${BASEDIR}/combined_ch"
for s in "${SAMPLES[@]}"; do
    if [[ ! -s "${CH_DIR}/${s}_CH.methylKit.gz" ]]; then
        echo "ERROR: Combined CH file missing: ${CH_DIR}/${s}_CH.methylKit.gz"
        MISSING=1
    fi
done
if [[ $MISSING -eq 1 ]]; then
    echo "Re-run with Step 1 enabled."
    exit 1
fi
echo "Step 1 skipped — all 8 combined CH files verified."
echo ""

STEP2A_JID=$(sbatch --parsable "${SCRIPT_DIR}/02_mch_aggregate_sample.sb")
echo "  Step 02 - Aggregate:     job $STEP2A_JID (array 0-7)"

STEP2B_JID=$(sbatch --parsable --dependency=afterok:${STEP2A_JID} \
    "${SCRIPT_DIR}/03_mch_differential.sb")
echo "  Step 03 - edgeR:         job $STEP2B_JID (after $STEP2A_JID)"

echo "  Step 06 - dmrseq:        submitting via scripts/dmrseq/submit_dmrseq.sh"
bash "${SCRIPT_DIR}/dmrseq/submit_dmrseq.sh"

STEP4_JID=$(sbatch --parsable --dependency=afterok:${STEP2B_JID} \
    "${SCRIPT_DIR}/04_mch_integration.sb")
echo "  Step 04 - Integration:   job $STEP4_JID (after $STEP2B_JID)"

echo ""
echo "========================================="
echo "Dependency chain:"
echo "  Step2a($STEP2A_JID) -> Step2b($STEP2B_JID) -> Step4($STEP4_JID)"
echo "  Step06 (dmrseq) — submitted via scripts/dmrseq/submit_dmrseq.sh (parallel)"
echo ""
echo "Monitor: squeue -u \$USER"
echo "Logs:    $LOGDIR"
echo "========================================="
