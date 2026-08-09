# emseq

Non-CG methylation in BAP1-knockout mouse cerebellum. EM-seq, 8 samples, mm10.

## What this project tests

BAP1 is a deubiquitinase. It removes H2AK119Ub, a Polycomb repressive mark. When BAP1 is lost, H2AK119Ub accumulates across the genome. The question: does that accumulation change gene-body mCH, and does MeCP2 binding follow?

The predicted axis: Polycomb/H2AK119Ub -> DNMT3A -> mCH -> MeCP2.

This repo tests the first link. BAP1-KO in Math1-Cre conditional knockout mice (adult cerebellum). 4 controls, 4 mutants, balanced male/female.

## Samples

| ID | Genotype | Sex | Coverage |
|----|----------|-----|----------|
| ctrl_M1 | control | M | ~38x |
| ctrl_M2 | control | M | ~57x |
| ctrl_F1 | control | F | ~54x |
| ctrl_F2 | control | F | ~68x |
| mut_M1 | mutant | M | ~55x |
| mut_M2 | mutant | M | ~118x |
| mut_F1 | mutant | F | ~101x |
| mut_F2 | mutant | F | ~53x |

Two deep outliers (mut_M2 and mut_F1). The edgeR model handles this with a `log(total_coverage)` offset term.

File prefix mapping (for example, `ctrl_M1` maps to `260629_Emseq_LV_index_E1_Bap1_Math1_adult_ctrl_M1_S1`) is in the `EMSEQ_PREFIX` associative array in `submit_mch_pipeline.sh` and individual `.sb` scripts.

## Pipeline architecture

Two parallel pipelines for gene-body differential mCH. Both use SLURM dependency chaining on SDSC Expanse.

### All-CH pipeline (`scripts/submit_mch_pipeline.sh`)

This is the primary analysis. CHH + CHG combined gives about 3x more sites per gene body than CA alone.

1. `01_mch_combine_ch.sb` -- merge CHH + CHG methylKit files per sample (array 0-7)
2. `02_mch_aggregate_sample.sb` / `02_mch_aggregate_sample.R` -- gene-body aggregation with GenomicRanges overlap (array 0-7, after step 1)
3. `03_mch_differential.sb` / `03_mch_differential.R` -- edgeR QL F-test with `~ genotype + sex + lambda_ch_rate`, `log(total_coverage)` offset (after step 2)
4. `04_mch_integration.sb` / `04_mch_integration.R` -- K119ub overlap, MeCP2 binding, neuronal gene set/gene-length interaction, CG x mCH cross-modality (after step 3)
5. `05_mch_volcano.R` -- volcano plots with direction, FDR stringency, gene-length coloring (run locally)
6. `06_mch_dmrseq.sb` / `06_mch_dmrseq.R` -- genome-wide CH DMR discovery (parallel with steps 2-4)

### CA-only pipeline (`scripts/ca/submit_ca_pipeline.sh`)

Exists for comparison with published mCA atlases. Secondary to the all-CH pipeline.

1. `ca/01_ca_filter.sb` / `ca/01_ca_filter.py` -- reference-genome lookup for CA trinucleotide context (pysam, array 0-15)
2. `ca/02_ca_combine.sb` -- merge CA-filtered CHH + CHG per sample
3. `ca/03_ca_aggregate.sb` -- reuses `02_mch_aggregate_sample.R` on CA-only input
4. `ca/04_ca_differential.sb` -- reuses `03_mch_differential.R` on CA-only aggregates

### Shared components

`02_mch_aggregate_sample.R` and `03_mch_differential.R` are shared between both pipelines. The CA pipeline passes different input/output paths but uses the same R scripts.

`scripts/utils/multi_format_output.R` saves every plot in 4 formats (PDF, SVG, PNG, JPEG) into per-figure subdirectories. All R scripts that produce plots source this utility.

## Key methodological decisions

**All-CH over CA-only.** CHH + CHG combined gives about 3x more sites per gene body than CA alone. The CA pipeline exists for comparison with published mCA atlases.

**Lambda as covariate, not subtracted.** Per-sample lambda spike-in conversion noise rates enter the edgeR model as a covariate. The noise floor (about 0.5-0.9%) is close to signal (about 1%), so pre-subtraction zeroes too many genes.

**Gene-body aggregation, not site-level.** Site-level differential mCH is underpowered. Aggregation gives about 3,050 CA sites per 50 kb gene and about 76,000 observations per gene per sample at 25x.

**FDR as binary gatekeeper.** In the downstream ChIP analyses (steps 07-08), FDR is a threshold (pass/fail at 0.05). It is not a continuous weight or confidence score. This follows the ASA 2016 statement on P-values (Wasserstein and Lazar, doi:10.1080/00031305.2016.1154108).

## Key results

8,767 genes are differentially methylated at FDR < 0.05. The shift is genome-wide: lambda_GC = 8.457. This is expected biology in a BAP1-KO, not a confounding artifact.

Direction: 7,466 genes lose mCH in the mutant. 1,301 gain it. The loss is 5.7x more common.

## Downstream analyses

### Step 07: MeCP2 vs H2AK119Ub quadrant analysis

Script: `scripts/07_mecp2_h2aub_quadrant.R`. Output: `results/07_quadrant/`.

Takes DiffBind peak output for MeCP2 and H2AK119Ub from `data/`. Annotates peaks to genes with ChIPseeker/TxDb. Filters to gene-body peaks (Promoter, Exon, Intron, 5'/3' UTR). Collapses multiple peaks per gene into one fold change: median of significant peaks, with FDR as a binary gatekeeper.

Produces `gene_level_results_v2.csv`. Step 08 reads this file.

7,371 genes have peaks in both marks. 1,403 are significant (at least one significant peak in both MeCP2 and H2AK119Ub). Physical overlap analysis (50% or more reciprocal overlap) is integrated.

GO enrichment runs on all 4 quadrants with corrected background (genes-with-peaks-in-both-marks) and goseq gene-length correction. The enrichGO terms do not survive goseq. The synapse/cognition enrichment from the original analysis was a compound artifact of wrong background + gene-length bias.

### Step 08: Three-way mCH integration

Script: `scripts/08_mch_mecp2_h2aub_integration.R`. Output: `results/08_three_way/`. Run step 07 first.

Reads the step 07 gene-level results and merges with mCH edgeR output. Deduplicates 23 mCH genes that have multiple ENSMUSG IDs (keeps the highest |logFC| per gene name).

6,351 genes have data across all 3 modalities. Produces 3 scatter plots:

1. mCH vs MeCP2 (Spearman rho = 0.38)
2. mCH vs H2AK119Ub (Spearman rho = 0.37)
3. MeCP2 vs H2AK119Ub, colored by mCH direction (blue = mCH down, red = mCH up)

The three-way scatter shows a clean directional split. Q3 (both ChIP marks down) has 990 of 992 mCH-significant genes losing mCH. Q1 (both marks up) is the only quadrant with a mixed mCH response.

Three-way GO enrichment (genes significant in both ChIP marks AND mCH) also does not survive goseq. Diagnostics confirm that both gene length and coverage are associated with mCH significance (both p < 2.2e-16).

## Reference files and cross-project dependencies

Step 04 (`04_mch_integration.R`) pulls reference files from both this project and the companion Hi-C project at `/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c`:

- Gene body BED: `references/gene_bodies.protein_coding.bed`
- K119ub consensus peaks: `references/K119ub_consensus_v3.bed`
- K119ub DiffBind: from `mariner_hi-c/peaks/diffbind/`
- MeCP2 consensus and DiffBind: `references/mecp2-peaks/`
- Neuronal gene set: from `mariner_hi-c/biomodal/downstream/`
- CG mC/hmC DMRs: from `mariner_hi-c/biomodal/downstream/modality/outputs/`

Local copies of K119ub and MeCP2 data files live in `data/`.

## Directory layout

```
scripts/              Pipeline scripts (numbered order)
  01-06               mCH pipeline (combine, aggregate, edgeR, integration, volcano, DMR)
  07_mecp2_h2aub_quadrant.R      MeCP2 vs H2AK119Ub gene-level analysis
  08_mch_mecp2_h2aub_integration.R  Three-way mCH + ChIP integration
  ca/                 CA-only pipeline variant
  utils/              Shared R utilities (multi_format_output.R)
results/              Pipeline output (numbered to match scripts)
  02_aggregate/       Gene-body aggregated counts
  03_differential/    edgeR differential results
  04_integration/     K119ub/MeCP2 integration
  07_quadrant/        MeCP2 vs H2AK119Ub quadrant analysis output
  08_three_way/       Three-way mCH + ChIP integration output
  ca/                 CA-only pipeline output
data/                 DiffBind peak files, reference BEDs, gencode gene bodies
docs/                 Project aims, abstract, context
logs/                 SLURM log files
archive/              Superseded scripts and prior results
```

## Compute environment

SDSC Expanse. Conda environment `emseq_env`. Account `csd940`.

Scripts are developed locally and transferred to Expanse. This local directory mirrors the remote layout.

## Dependencies

R: edgeR, GenomicRanges, data.table, bsseq, dmrseq, clusterProfiler, goseq, ggplot2, patchwork, svglite, ChIPseeker.

Python: pysam (CA-filter step only).

Full list in `scripts/00_install_r_packages.R`.
