#!/bin/bash
# scripts/utils/copy_reference_data.sh
#
# Copy the reference data the section pipeline needs from the companion
# mariner_hi-c repository into this repository's data/ directory.
#
# Run once during setup, before any section job. Every path below was verified
# against the mariner_hi-c layout. The script stops at the first missing source
# file so that a partial copy never looks like a complete one.
#
# BigWig tracks are NOT copied. They are large and _shared_config.R reads them
# in place through the bigwigs_dir root in paths.yaml.
#
# Usage:
#   bash scripts/utils/copy_reference_data.sh <hic_dir> <data_dir>
#
# Example (Expanse):
#   bash scripts/utils/copy_reference_data.sh \
#       /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c \
#       /expanse/lustre/projects/csd940/zalibhai/emseq
#
# Example (local):
#   bash scripts/utils/copy_reference_data.sh \
#       /Users/zakiralibhai/Documents/GitHub/mariner_hi-c \
#       /Users/zakiralibhai/emseq

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <hic_dir> <data_dir>" >&2
    exit 1
fi

HIC_DIR="$1"
DATA_DIR="$2"
BIO_DIR="${HIC_DIR}/biomodal/downstream"

if [[ ! -d "$HIC_DIR" ]]; then
    echo "ERROR: mariner_hi-c directory not found: $HIC_DIR" >&2
    exit 1
fi
if [[ ! -d "$DATA_DIR" ]]; then
    echo "ERROR: data directory not found: $DATA_DIR" >&2
    exit 1
fi

echo "========================================="
echo "Copying section pipeline reference data"
echo "========================================="
echo "Source (mariner_hi-c): $HIC_DIR"
echo "Target (emseq data):   $DATA_DIR"
echo ""

for sub in chip_peaks chromatin diffbind hic mecp2 neuronal features; do
    mkdir -p "${DATA_DIR}/data/${sub}"
done

N_COPIED=0

# copy_one <source> <destination>
copy_one() {
    local src="$1"
    local dst="$2"

    if [[ ! -s "$src" ]]; then
        echo "ERROR: source file missing or empty: $src" >&2
        exit 1
    fi

    cp "$src" "$dst"
    local size
    size=$(wc -c < "$dst" | tr -d ' ')
    printf '  %-46s %12s bytes\n' "$(basename "$dst")" "$size"
    N_COPIED=$((N_COPIED + 1))
}

echo "--- Histone mark consensus peaks -> data/chip_peaks/ ---"
copy_one "${HIC_DIR}/peaks/CTCF.bed"                                  "${DATA_DIR}/data/chip_peaks/CTCF.bed"
copy_one "${HIC_DIR}/peaks/beds/H3K27acCerebellumLate2.bed"           "${DATA_DIR}/data/chip_peaks/H3K27ac.bed"
copy_one "${HIC_DIR}/peaks/beds/H3K27me3CerebellumLate1.bed"          "${DATA_DIR}/data/chip_peaks/H3K27me3.bed"
copy_one "${HIC_DIR}/peaks/beds/H3K4me1CerebellumLate1.bed"           "${DATA_DIR}/data/chip_peaks/H3K4me1.bed"
copy_one "${HIC_DIR}/peaks/beds/H3K4me3CerebellumLate2.bed"           "${DATA_DIR}/data/chip_peaks/H3K4me3.bed"
copy_one "${HIC_DIR}/peaks/beds/Bivalent_Cerebellum_Late.bed"         "${DATA_DIR}/data/chip_peaks/Bivalent.bed"
echo ""

echo "--- chromHMM emission states -> data/chromatin/ ---"
copy_one "${HIC_DIR}/peaks/251230-Challana-EmissionState12-activeenhancer.bed" "${DATA_DIR}/data/chromatin/activeenhancer.bed"
copy_one "${HIC_DIR}/peaks/251230-Challana-EmissionState12-activepromoter.bed" "${DATA_DIR}/data/chromatin/activepromoter.bed"
copy_one "${HIC_DIR}/peaks/251230-Challana-EmissionState3-bivalent.bed"        "${DATA_DIR}/data/chromatin/bivalent.bed"
echo ""

echo "--- Quantitative DiffBind results -> data/diffbind/ ---"
copy_one "${HIC_DIR}/peaks/diffbind/ATAC_allATAC_diffbind_results_summit_appended_ap.txt" "${DATA_DIR}/data/diffbind/ATAC_diffbind.txt"
copy_one "${HIC_DIR}/peaks/diffbind/K27ac_diffbind_results_summit_appended_ap.txt"        "${DATA_DIR}/data/diffbind/K27ac_diffbind.txt"
copy_one "${HIC_DIR}/peaks/diffbind/K27me3_diffbind_results_summit_appended_ap.txt"       "${DATA_DIR}/data/diffbind/K27me3_diffbind.txt"
echo ""

echo "--- ATAC peak sets -> data/diffbind/ ---"
copy_one "${HIC_DIR}/peaks/atac_seq/ATAC_up.bed"           "${DATA_DIR}/data/diffbind/ATAC_up.bed"
copy_one "${HIC_DIR}/peaks/atac_seq/ATAC_down.bed"         "${DATA_DIR}/data/diffbind/ATAC_down.bed"
copy_one "${HIC_DIR}/peaks/atac_seq/consensus_control.bed" "${DATA_DIR}/data/diffbind/ATAC_consensus_ctrl.bed"
copy_one "${HIC_DIR}/peaks/atac_seq/consensus_mutant.bed"  "${DATA_DIR}/data/diffbind/ATAC_consensus_mut.bed"
echo ""

echo "--- Condition-specific peak sets -> data/diffbind/ ---"
copy_one "${HIC_DIR}/peaks/intersect/P51_K119ub_ctrl_intersect.bed" "${DATA_DIR}/data/diffbind/K119ub_ctrl.bed"
copy_one "${HIC_DIR}/peaks/intersect/P51_K119ub_mut_intersect.bed"  "${DATA_DIR}/data/diffbind/K119ub_mut.bed"
copy_one "${HIC_DIR}/peaks/intersect/P60_K27ac_ctrl_intersect.bed"  "${DATA_DIR}/data/diffbind/K27ac_ctrl.bed"
copy_one "${HIC_DIR}/peaks/intersect/P60_K27ac_mut_intersect.bed"   "${DATA_DIR}/data/diffbind/K27ac_mut.bed"
echo ""

echo "--- MeCP2 peak-level and aging data -> data/mecp2/ ---"
copy_one "${HIC_DIR}/peaks/mecp2/MeCP2_annotated.txt" "${DATA_DIR}/data/mecp2/MeCP2_annotated.txt"
copy_one "${HIC_DIR}/peaks/mecp2/MeCP2_up.bed"        "${DATA_DIR}/data/mecp2/MeCP2_up.bed"
copy_one "${HIC_DIR}/peaks/mecp2/MeCP2_down.bed"      "${DATA_DIR}/data/mecp2/MeCP2_down.bed"
copy_one "${BIO_DIR}/peaks/MeCP2_ctrl_adultvsyoung_diffbind_results.txt" "${DATA_DIR}/data/mecp2/MeCP2_ctrl_aging_diffbind.txt"
copy_one "${BIO_DIR}/peaks/MeCP2_mut_adultvsyoung_diffbind_results.txt"  "${DATA_DIR}/data/mecp2/MeCP2_mut_aging_diffbind.txt"
echo ""

echo "--- Hi-C loops, compartments, subcompartments -> data/hic/ ---"
copy_one "${HIC_DIR}/peaks/loop_annotation_extended/late/extended_characterized_loops.tsv" "${DATA_DIR}/data/hic/characterized_loops.tsv"
copy_one "${HIC_DIR}/tads/tad-pc-analysis/inputs/late/diffPC/diffcompartments.txt"         "${DATA_DIR}/data/hic/diffcompartments.txt"
copy_one "${HIC_DIR}/ML/cmpts/outputs/calder2/late/250402_subcompartment_labels_100kb.tsv" "${DATA_DIR}/data/hic/calder2_subcompartments_100kb.tsv"
echo ""

echo "--- Per-gene K119ub signal -> data/ ---"
copy_one "${BIO_DIR}/data/k119ub_gene_signal.tsv" "${DATA_DIR}/data/k119ub_gene_signal.tsv"
echo ""

echo "--- Neuronal gene sets -> data/neuronal/ ---"
copy_one "${BIO_DIR}/plots/visualizations/tables/72_neuronal_gene_set_go_derived.tsv" "${DATA_DIR}/data/neuronal/neuronal_gene_set_go_derived.tsv"
copy_one "${BIO_DIR}/plots/visualizations/tables/76_synapse_axon_gene_set.tsv"        "${DATA_DIR}/data/neuronal/synapse_axon_gene_set.tsv"
echo ""

echo "========================================="
echo "Copied ${N_COPIED} files."
echo ""
echo "Still required before section jobs can run:"
echo "  data/features/*.bed  — run scripts/utils/generate_feature_beds.R"
echo "========================================="
