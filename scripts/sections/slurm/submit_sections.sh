#!/bin/bash
# scripts/sections/slurm/submit_sections.sh
#
# Submits the 23 mCH section jobs to SLURM in dependency order.
#
# Every section has one wrapper in this directory named <section>.sb. Sections
# that read another section's output are chained with --dependency=afterok, so
# a consumer starts only after its producer exits 0.
#
# Usage:
#   bash scripts/sections/slurm/submit_sections.sh
#   bash scripts/sections/slurm/submit_sections.sh --dry-run
#   bash scripts/sections/slurm/submit_sections.sh 70_02_chromatin_remodeling
#   bash scripts/sections/slurm/submit_sections.sh --ntimes 20000 40_04_gene_level
#   bash scripts/sections/slurm/submit_sections.sh --force-rerun 40_01_dmr_marks
#
# The wrappers write to logs/sections/<name>_%j.out, a path relative to the
# code directory, so this script changes to the code directory before it calls
# sbatch.

set -euo pipefail

# ---------------------------------------------------------------------------
# Fixed locations
# ---------------------------------------------------------------------------

CODEDIR="${EMSEQ_CODE_DIR:-/expanse/lustre/projects/csd940/zalibhai/emseq-repo}"
SLURM_SUBDIR="scripts/sections/slurm"
LOG_SUBDIR="logs/sections"

# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------

ALL_SECTIONS=(
    10_01_chromatin_state
    10_02_ab_compartment
    10_03_polycomb_enrichment
    10_04_subcompartment
    20_01_mecp2_correlation
    20_02_multi_mark_diffbind
    20_03_quadrant_scatters
    20_04_mch_mecp2_by_mark
    30_01_loop_anchor_methylation
    30_02_mecp2_loop_anchors
    40_01_dmr_marks
    40_02_atac_loops
    40_03_domains
    40_04_gene_level
    50_01_feature_methylation
    60_01_methylation_scale
    60_02_k119ub_unmethylated
    60_03_reconciliation
    60_04_aging
    70_01_k119ub_neuronal
    70_02_chromatin_remodeling
    70_03_geneset_overlap
    70_04_synapse_chromatin
)
N_ALL_SECTIONS=${#ALL_SECTIONS[@]}

# Sections whose Fisher tests section 40_04 validates by label shuffle.
# 40_04 reads results/sections/fisher_test_registry.tsv, so it waits for every
# section named here.
FISHER_SECTIONS=(
    10_01_chromatin_state
    10_02_ab_compartment
    10_03_polycomb_enrichment
    20_01_mecp2_correlation
    20_02_multi_mark_diffbind
    30_01_loop_anchor_methylation
    30_02_mecp2_loop_anchors
    60_01_methylation_scale
    60_02_k119ub_unmethylated
    60_04_aging
    70_01_k119ub_neuronal
    70_03_geneset_overlap
)

# Sections that read NTIMES and FORCE_RERUN from the environment.
PERMUTATION_PREFIX="40_"

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------

DRY_RUN=0
FORCE_RERUN=0
NTIMES=""
TARGET_SECTION=""
AFTER_02B=""

usage() {
    cat <<'USAGE'
submit_sections.sh -- submit the mCH section pipeline to SLURM

  bash scripts/sections/slurm/submit_sections.sh [options] [section]

Options
  --dry-run          Print every sbatch command; submit nothing. Job ids in the
                     printed commands and in the summary are placeholders.
  --ntimes N         Permutation count for the 40_* sections. Exported as
                     NTIMES. Without it each section uses its own default
                     (5000 for 40_01, 40_02, 40_03; 10000 for 40_04).
  --force-rerun      Export FORCE_RERUN=1 to the 40_* sections, which makes
                     them ignore their RDS caches and recompute.
  --after-02b JOBID  Chain 50_01_feature_methylation to wait on this step 02b
                     job. Without it 50_01 runs immediately and fails if the
                     feature tables are not ready.
  -h, --help         Print this text.

Arguments
  section            Submit only this section, with no dependency. Give the
                     full name (70_02_chromatin_remodeling) or a unique
                     prefix (70_02). Without it, all 23 sections are submitted
                     in dependency order.

Environment
  EMSEQ_CODE_DIR     Code directory holding scripts/ and logs/.
                     Default: /expanse/lustre/projects/csd940/zalibhai/emseq-repo
USAGE
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --force-rerun)
            FORCE_RERUN=1
            shift
            ;;
        --after-02b)
            [[ $# -ge 2 ]] || die "--after-02b needs a SLURM job id."
            AFTER_02B="$2"
            shift 2
            ;;
        --after-02b=*)
            AFTER_02B="${1#--after-02b=}"
            shift
            ;;
        --ntimes)
            [[ $# -ge 2 ]] || die "--ntimes needs a value."
            NTIMES="$2"
            shift 2
            ;;
        --ntimes=*)
            NTIMES="${1#--ntimes=}"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            usage >&2
            die "Unknown option: $1"
            ;;
        *)
            [[ -z "$TARGET_SECTION" ]] || die "Only one section can be named. Got '$TARGET_SECTION' and '$1'."
            TARGET_SECTION="$1"
            shift
            ;;
    esac
done

if [[ -n "$NTIMES" ]] && [[ ! "$NTIMES" =~ ^[0-9]+$ ]]; then
    die "--ntimes must be a whole number, got '$NTIMES'."
fi

# ---------------------------------------------------------------------------
# Section-name helpers
# ---------------------------------------------------------------------------

# "10_01_chromatin_state" -> "10_01"
short_name() {
    local full="$1"
    local head="${full%%_*}"
    local rest="${full#*_}"
    printf '%s_%s' "$head" "${rest%%_*}"
}

# Accept a full section name or a unique prefix. Print the full name.
resolve_section() {
    local want="$1"
    local i match=""
    local n_match=0

    want="${want%.sb}"

    for (( i = 0; i < N_ALL_SECTIONS; i++ )); do
        if [[ "${ALL_SECTIONS[$i]}" == "$want" ]]; then
            printf '%s' "$want"
            return 0
        fi
    done

    for (( i = 0; i < N_ALL_SECTIONS; i++ )); do
        case "${ALL_SECTIONS[$i]}" in
            "$want"*)
                match="${ALL_SECTIONS[$i]}"
                n_match=$(( n_match + 1 ))
                ;;
        esac
    done

    if [[ $n_match -eq 1 ]]; then
        printf '%s' "$match"
        return 0
    fi

    {
        if [[ $n_match -eq 0 ]]; then
            echo "ERROR: no section matches '$want'."
        else
            echo "ERROR: '$want' matches $n_match sections."
        fi
        echo "Sections:"
        for (( i = 0; i < N_ALL_SECTIONS; i++ )); do
            echo "  ${ALL_SECTIONS[$i]}"
        done
    } >&2
    return 1
}

# ---------------------------------------------------------------------------
# Submission bookkeeping
# ---------------------------------------------------------------------------

N_ROWS=0
SUMMARY_SECTION=()
SUMMARY_JOBID=()
SUMMARY_WAIT=()
DRY_NEXT_JOBID=9000001

jobid_of() {
    local want="$1"
    local i
    for (( i = 0; i < N_ROWS; i++ )); do
        if [[ "${SUMMARY_SECTION[$i]}" == "$want" ]]; then
            printf '%s' "${SUMMARY_JOBID[$i]}"
            return 0
        fi
    done
    echo "ERROR: no job id recorded for section '$want'. Submission order is wrong." >&2
    return 1
}

# submit_section <section> [dependency-section ...]
submit_section() {
    local section="$1"
    shift

    local script_rel="${SLURM_SUBDIR}/${section}.sb"
    local dep_ids=""
    local dep_names=""
    local d dep_id

    for d in "$@"; do
        dep_id=$(jobid_of "$d")
        if [[ -z "$dep_ids" ]]; then
            dep_ids="$dep_id"
            dep_names="$(short_name "$d")"
        else
            dep_ids="${dep_ids}:${dep_id}"
            dep_names="${dep_names}, $(short_name "$d")"
        fi
    done

    local export_spec="ALL"
    case "$section" in
        "${PERMUTATION_PREFIX}"*)
            if [[ -n "$NTIMES" ]]; then
                export_spec="${export_spec},NTIMES=${NTIMES}"
            fi
            if [[ $FORCE_RERUN -eq 1 ]]; then
                export_spec="${export_spec},FORCE_RERUN=1"
            fi
            ;;
    esac

    local -a cmd
    cmd=(sbatch --parsable --export="$export_spec")
    if [[ -n "$dep_ids" ]]; then
        cmd+=(--dependency="afterok:${dep_ids}" --kill-on-invalid-dep=yes)
    fi
    cmd+=("$script_rel")

    local jid=""
    if [[ $DRY_RUN -eq 1 ]]; then
        echo "${cmd[*]}"
        jid="$DRY_NEXT_JOBID"
        DRY_NEXT_JOBID=$(( DRY_NEXT_JOBID + 1 ))
    else
        if ! jid=$("${cmd[@]}"); then
            die "sbatch failed for section ${section}. Command: ${cmd[*]}"
        fi
        jid="${jid%%;*}"
        if [[ ! "$jid" =~ ^[0-9]+$ ]]; then
            die "Could not read a job id from sbatch for ${section}. Output: '${jid}'"
        fi
        printf '  %-30s job %-10s waits on: %s\n' \
            "$section" "$jid" "${dep_names:--}"
    fi

    SUMMARY_SECTION[$N_ROWS]="$section"
    SUMMARY_JOBID[$N_ROWS]="$jid"
    SUMMARY_WAIT[$N_ROWS]="${dep_names:--}"
    N_ROWS=$(( N_ROWS + 1 ))
}

print_summary() {
    local i all_ids=""
    echo ""
    echo "========================================================================"
    if [[ $DRY_RUN -eq 1 ]]; then
        echo "Dry run: nothing was submitted. Job ids below are placeholders."
    else
        echo "Submitted ${N_ROWS} job(s)."
    fi
    echo "========================================================================"
    printf '%-10s  %-30s  %s\n' "JOB ID" "SECTION" "WAITS ON"
    printf '%-10s  %-30s  %s\n' "----------" "------------------------------" "--------"
    for (( i = 0; i < N_ROWS; i++ )); do
        printf '%-10s  %-30s  %s\n' \
            "${SUMMARY_JOBID[$i]}" "${SUMMARY_SECTION[$i]}" "${SUMMARY_WAIT[$i]}"
        if [[ -z "$all_ids" ]]; then
            all_ids="${SUMMARY_JOBID[$i]}"
        else
            all_ids="${all_ids} ${SUMMARY_JOBID[$i]}"
        fi
    done
    echo "========================================================================"
    echo "Logs:    ${CODEDIR}/${LOG_SUBDIR}/<section>_<jobid>.out"
    if [[ $DRY_RUN -eq 0 ]]; then
        echo "Monitor: squeue -u \$USER"
        echo "Cancel:  scancel ${all_ids}"
    fi
}

# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------

[[ -d "$CODEDIR" ]] || die "Code directory not found: ${CODEDIR}. Set EMSEQ_CODE_DIR."
cd "$CODEDIR"
[[ -d "$SLURM_SUBDIR" ]] || die "Wrapper directory not found: ${CODEDIR}/${SLURM_SUBDIR}"

PLANNED_SECTIONS=()
N_PLANNED=0

if [[ -n "$TARGET_SECTION" ]]; then
    TARGET_SECTION=$(resolve_section "$TARGET_SECTION")
    PLANNED_SECTIONS[0]="$TARGET_SECTION"
    N_PLANNED=1
else
    for (( i = 0; i < N_ALL_SECTIONS; i++ )); do
        PLANNED_SECTIONS[$i]="${ALL_SECTIONS[$i]}"
    done
    N_PLANNED=$N_ALL_SECTIONS
fi

for (( i = 0; i < N_PLANNED; i++ )); do
    wrapper="${SLURM_SUBDIR}/${PLANNED_SECTIONS[$i]}.sb"
    [[ -f "$wrapper" ]] || die "Wrapper not found: ${CODEDIR}/${wrapper}"
done

if [[ -n "$NTIMES" || $FORCE_RERUN -eq 1 ]]; then
    n_permutation=0
    for (( i = 0; i < N_PLANNED; i++ )); do
        case "${PLANNED_SECTIONS[$i]}" in
            "${PERMUTATION_PREFIX}"*) n_permutation=$(( n_permutation + 1 )) ;;
        esac
    done
    if [[ $n_permutation -eq 0 ]]; then
        die "--ntimes and --force-rerun only reach the ${PERMUTATION_PREFIX}* sections, and none is being submitted."
    fi
fi

if [[ $DRY_RUN -eq 0 ]]; then
    command -v sbatch >/dev/null 2>&1 || die "sbatch is not on PATH. Submit from an Expanse login node."
    mkdir -p "$LOG_SUBDIR"
fi

echo "========================================================================"
echo "EM-seq mCH sections: SLURM submission"
echo "========================================================================"
echo "Code directory:  ${CODEDIR}"
echo "Wrappers:        ${CODEDIR}/${SLURM_SUBDIR}"
echo "Job logs:        ${CODEDIR}/${LOG_SUBDIR}"
echo "Sections:        ${N_PLANNED}"
echo "Dry run:         $( [[ $DRY_RUN -eq 1 ]] && echo yes || echo no )"
echo "NTIMES:          ${NTIMES:-section default}"
echo "FORCE_RERUN:     $( [[ $FORCE_RERUN -eq 1 ]] && echo 1 || echo "not set" )"
echo "Start time:      $(date)"
echo "========================================================================"
echo ""

# ---------------------------------------------------------------------------
# Single section: no dependency
# ---------------------------------------------------------------------------

if [[ -n "$TARGET_SECTION" ]]; then
    echo "Single section, dependencies skipped."
    echo ""
    submit_section "$TARGET_SECTION"
    print_summary
    exit 0
fi

# ---------------------------------------------------------------------------
# Wave 1: no dependencies
# ---------------------------------------------------------------------------

echo "Wave 1: independent sections"
submit_section 10_01_chromatin_state
submit_section 10_02_ab_compartment
submit_section 10_04_subcompartment
submit_section 20_01_mecp2_correlation
submit_section 30_01_loop_anchor_methylation
submit_section 30_02_mecp2_loop_anchors
submit_section 40_01_dmr_marks
if [[ -n "$AFTER_02B" ]]; then
    # 50_01 reads the feature tables step 02b writes. Chain it.
    SUMMARY_SECTION[$N_ROWS]="step_02b"
    SUMMARY_JOBID[$N_ROWS]="$AFTER_02B"
    SUMMARY_WAIT[$N_ROWS]="(external)"
    N_ROWS=$(( N_ROWS + 1 ))
    submit_section 50_01_feature_methylation step_02b
else
    submit_section 50_01_feature_methylation
fi
submit_section 60_01_methylation_scale
submit_section 60_03_reconciliation
submit_section 60_04_aging
submit_section 70_01_k119ub_neuronal
echo ""

# ---------------------------------------------------------------------------
# Wave 2: one producer each
# ---------------------------------------------------------------------------

echo "Wave 2: consumers of 10_01, 60_01, 70_01"
submit_section 10_03_polycomb_enrichment    10_01_chromatin_state
submit_section 20_02_multi_mark_diffbind    10_01_chromatin_state
submit_section 40_02_atac_loops             10_01_chromatin_state
submit_section 40_03_domains                10_01_chromatin_state
submit_section 60_02_k119ub_unmethylated    60_01_methylation_scale
submit_section 70_03_geneset_overlap        70_01_k119ub_neuronal
echo ""

# ---------------------------------------------------------------------------
# Wave 3: consumers of 20_02
# ---------------------------------------------------------------------------

echo "Wave 3: consumers of 20_02"
submit_section 20_03_quadrant_scatters      20_02_multi_mark_diffbind
submit_section 20_04_mch_mecp2_by_mark      20_02_multi_mark_diffbind
submit_section 70_02_chromatin_remodeling   70_01_k119ub_neuronal 20_02_multi_mark_diffbind
submit_section 70_04_synapse_chromatin      70_01_k119ub_neuronal 20_02_multi_mark_diffbind
echo ""

# ---------------------------------------------------------------------------
# Wave 4: Fisher registry validation, last
# ---------------------------------------------------------------------------

echo "Wave 4: Fisher registry permutation"
submit_section 40_04_gene_level "${FISHER_SECTIONS[@]}"
echo ""

EXPECTED_ROWS=$N_ALL_SECTIONS
[[ -n "$AFTER_02B" ]] && EXPECTED_ROWS=$(( EXPECTED_ROWS + 1 ))
if [[ $N_ROWS -ne $EXPECTED_ROWS ]]; then
    die "Submitted ${N_ROWS} sections, expected ${EXPECTED_ROWS}."
fi

print_summary
echo "End time: $(date)"
