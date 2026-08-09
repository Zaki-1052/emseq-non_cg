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

## Pipeline

Two parallel pipelines test gene-body differential mCH. Both run on SDSC Expanse (SLURM).

**All-CH pipeline** (`scripts/submit_mch_pipeline.sh`): merges CHH + CHG per sample, aggregates to gene bodies, runs edgeR QL F-test with `~ genotype + sex + lambda_ch_rate`. This is the primary analysis. About 3x more sites per gene than CA alone.

**CA-only pipeline** (`scripts/ca/submit_ca_pipeline.sh`): filters to CA trinucleotide context first, then runs the same aggregation and edgeR steps. Exists for comparison with published mCA atlases.

Both pipelines share `02_mch_aggregate_sample.R` and `03_mch_differential.R`. The CA pipeline passes different input/output paths but uses the same R scripts.

## Key results

8,767 genes are differentially methylated at FDR < 0.05. The shift is genome-wide: lambda_GC = 8.457. This is expected biology in a BAP1-KO, not a confounding artifact.

Direction: 7,466 genes lose mCH in the mutant. 1,301 gain it. The loss is 5.7x more common.

## Downstream analyses

### MeCP2 vs H2AK119Ub quadrant analysis

Directory: `adult_mecp2_vs_h2aub/`

Takes DiffBind peak output for MeCP2 and H2AK119Ub. Annotates peaks to genes with ChIPseeker. Collapses multiple peaks per gene into one fold change (median of significant peaks, FDR as binary gatekeeper). Plots gene-level MeCP2 FC vs H2AK119Ub FC as a quadrant scatter.

Active script: `adult_MeCP2vsH2AUb_v2.R`. Produces `gene_level_results_v2.csv`.

GO enrichment runs on all 4 quadrants with corrected background and goseq gene-length correction. The enrichGO terms do not survive goseq. The synapse/cognition enrichment from the original analysis was a compound artifact of wrong background + gene-length bias.

### Three-way mCH integration

Directory: `mecp2_h2aub_mch/`

Reads the v2 quadrant output, merges with mCH edgeR results, and produces 3 scatter plots:

1. mCH vs MeCP2 (Spearman rho = 0.38)
2. mCH vs H2AK119Ub (Spearman rho = 0.37)
3. MeCP2 vs H2AK119Ub, colored by mCH direction (blue = mCH down, red = mCH up)

Active script: `mecp2_h2aub_mch_v2.R`. Depends on the v2 quadrant analysis. Run that first.

The three-way scatter shows a clean directional split. Q3 (both ChIP marks down) has 990 of 992 mCH-significant genes losing mCH. Q1 (both marks up) is the only quadrant with a mixed mCH response.

Three-way GO enrichment (genes significant in both ChIP marks AND mCH) also does not survive goseq.

## Directory layout

```
scripts/           Pipeline scripts (numbered order)
scripts/ca/        CA-only pipeline variant
scripts/utils/     Shared R utilities (multi_format_output.R)
results/           Pipeline output
results/ca/        CA-only output
adult_mecp2_vs_h2aub/  MeCP2 vs H2AK119Ub quadrant analysis
mecp2_h2aub_mch/       Three-way mCH + ChIP integration
data/              Local copies of reference peak files
docs/              Project aims, abstract, context
logs/              SLURM log files
archive/           Superseded scripts
```

## Compute environment

SDSC Expanse. Conda environment `emseq_env`. Account `csd940`.

Scripts are developed locally and transferred to Expanse. This local directory mirrors the remote layout.

## Dependencies

R: edgeR, GenomicRanges, data.table, bsseq, dmrseq, clusterProfiler, goseq, ggplot2, patchwork, svglite, ChIPseeker.

Python: pysam (CA-filter step only).

Full list in `scripts/00_install_r_packages.R`.
