# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Non-CG methylation analysis in BAP1-knockout mouse cerebellum. Tests whether BAP1 loss (via H2AK119ub accumulation) reshapes gene-body mCH that MeCP2 reads — the first test of the Polycomb/H2AK119ub → DNMT3A → mCH → MeCP2 axis. Eight EM-seq samples (4 ctrl, 4 mut; balanced M/F), aligned to mm10. The upstream EM-seq Nextflow pipeline and MethylDackel extraction are complete and not rerun from here.

The BAP1-KO genome-wide mCH shift is expected biology, not a confounding artifact. High lambda_GC values in QQ plots are consistent with this.

## Compute Environment

All pipeline scripts run on SDSC Expanse (SLURM) via conda environment `emseq_env`. Account `csd940`. Base directory on Expanse: `/expanse/lustre/projects/csd940/zalibhai/emseq`.

This local working directory mirrors the remote scripts and results. Scripts are developed locally and transferred to Expanse.

## Pipeline Architecture

Two parallel pipelines for gene-body differential mCH testing, both using SLURM dependency chaining:

### All-CH pipeline (`scripts/submit_mch_pipeline.sh`)

1. `01_mch_combine_ch.sb` — merge CHH + CHG methylKit files per sample (array 0-7)
2. `02_mch_aggregate_sample.sb` → `02_mch_aggregate_sample.R` — gene-body aggregation using GenomicRanges overlap (array 0-7, after step 1)
3. `03_mch_differential.sb` → `03_mch_differential.R` — edgeR QL F-test with `~ genotype + sex + lambda_ch_rate`, `log(total_coverage)` offset (after step 2)
4. `04_mch_integration.sb` → `04_mch_integration.R` — K119ub overlap, MeCP2 binding, neuronal gene set/gene-length interaction, CG×mCH cross-modality (after step 3)
5. `05_mch_volcano.R` — volcano plots with direction, FDR stringency, and gene-length coloring (run locally or interactively)
6. `06_mch_dmrseq.sb` → `06_mch_dmrseq.R` — genome-wide CH DMR discovery, parallel with steps 2–4

### CA-only pipeline (`scripts/ca/submit_ca_pipeline.sh`)

1. `ca/01_ca_filter.sb` → `ca/01_ca_filter.py` — reference-genome lookup for CA trinucleotide context (pysam, array 0-15)
2. `ca/02_ca_combine.sb` — merge CA-filtered CHH + CHG per sample
3. `ca/03_ca_aggregate.sb` — reuses `02_mch_aggregate_sample.R` on CA-only input
4. `ca/04_ca_differential.sb` — reuses `03_mch_differential.R` on CA-only aggregates

### Shared components

- `02_mch_aggregate_sample.R` and `03_mch_differential.R` are shared between the all-CH and CA-only pipelines. The CA pipeline passes different input/output paths but uses the same R scripts.
- `utils/multi_format_output.R` — saves every plot in four formats (PDF, SVG, PNG, JPEG) into per-figure subdirectories. Sourced by all R scripts that produce plots.

## Key Methodological Decisions

- **All-CH, not CA-only**: CHH + CHG combined gives ~3x more sites per gene body than CA alone. The CA pipeline exists for comparison with published mCA atlases but is secondary.
- **Lambda as covariate, not subtracted**: per-sample lambda spike-in conversion noise rates enter the edgeR model as a covariate rather than being pre-subtracted from counts. The noise floor (~0.5–0.9%) is close to signal (~1%), so subtraction zeroes too many genes.
- **Gene-body aggregation, not site-level**: site-level differential mCH is underpowered. Aggregating ~3,050 CA sites per 50kb gene gives ~76,000 observations per gene per sample at 25x.

## Samples

Eight samples with IDs `ctrl_M1`, `ctrl_M2`, `ctrl_F1`, `ctrl_F2`, `mut_M1`, `mut_M2`, `mut_F1`, `mut_F2`. Coverage ranges from ~38x to ~118x with two deep outliers (mut_M2 ~118x, mut_F1 ~101x), handled by the edgeR offset term.

File prefix mapping (e.g., `ctrl_M1` → `260629_Emseq_LV_index_E1_Bap1_Math1_adult_ctrl_M1_S1`) is defined in the `EMSEQ_PREFIX` associative array in `submit_mch_pipeline.sh` and individual `.sb` scripts.

## Reference Files and Cross-Project Dependencies

The integration step (`04_mch_integration.R`) pulls reference files from both this project and the companion Hi-C project at `/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c`:
- Gene body BED: `references/gene_bodies.protein_coding.bed`
- K119ub consensus peaks: `references/K119ub_consensus_v3.bed`
- K119ub DiffBind: from `mariner_hi-c/peaks/diffbind/`
- MeCP2 consensus and DiffBind: `references/mecp2-peaks/`
- Neuronal gene set: from `mariner_hi-c/biomodal/downstream/`
- CG mC/hmC DMRs: from `mariner_hi-c/biomodal/downstream/modality/outputs/`

Local copies of K119ub and MeCP2 data files live in `data/`.

## Directory Layout

- `scripts/` — pipeline scripts (numbered execution order)
- `scripts/ca/` — CA-only pipeline variant
- `scripts/utils/` — shared R utilities
- `results/` — pipeline output (aggregated counts, differential tables, integration TSVs, plots)
- `results/ca/` — CA-only pipeline output (mirrors `results/` structure)
- `data/` — local copies of reference peak files
- `docs/` — project aims, abstract, handoff context
- `plans/` — analysis plans
- `logs/` — SLURM log files (organized by step)
- `adult_mecp2_vs_h2aub/` — separate MeCP2 vs H2AK119ub quadrant analysis (standalone R scripts)
- `archive/` — superseded scripts and prior results
- `txt/` — meeting notes, preliminary results text

## R Dependencies

All packages listed in `scripts/00_install_r_packages.R`. Key packages: edgeR, GenomicRanges, data.table, bsseq, dmrseq, fgsea, clusterProfiler, ggplot2, patchwork, svglite. Python: pysam (CA-filter step only).
