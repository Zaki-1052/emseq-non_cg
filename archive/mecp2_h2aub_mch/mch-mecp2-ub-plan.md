# Plan: mCH + MeCP2 + H2AK119Ub Three-Way Integration (v2)

## Context

The `mecp2_h2aub_mch/` directory contains Janvi's second iteration of the MeCP2 vs H2AK119Ub analysis, extended with non-CG methylation (mCH) integration. Three scripts exist:

1. `mecp2_h2aub_all_elements_txdb.R` (886 lines) — ChIPseeker/TxDb annotation approach
2. `runfirst_mecp2_h2aub_genebodies_gencode.R` (789 lines) — gencode BED annotation approach
3. `runsecond_mecp2_h2aub_noncg_gencode.R` (490 lines) — merges ChIP results with mCH edgeR output

A code review (`review-mch-mecp2-ub.md`) identified 11 issues. The upstream MeCP2/H2AK119Ub gene-level analysis has already been fixed in `adult_mecp2_vs_h2aub/adult_MeCP2vsH2AUb_v2.R`, which produces `gene_level_results_v2.csv` with corrected median-of-significant-peaks aggregation. The new script reads that output instead of redoing the upstream analysis, and focuses on the mCH integration.

## Issues from review and how each is resolved

| # | Issue | Resolution |
|---|-------|------------|
| 1 | `runfirst` references undefined `gene_bodies_encode` variable | Eliminated — script reads v2 output, no gencode BED needed |
| 2 | Duplicate gene entries in mCH edgeR (23 genes) | Deduplicate: keep highest \|logFC\| per gene_name |
| 3 | FDR used as continuous weight (both scripts) | Eliminated — v2 upstream already uses median aggregation |
| 4 | FDR floor of 1e-300 | Eliminated — v2 upstream does not use FDR as weight |
| 5 | No goseq / gene-length correction | Add goseq alongside enrichGO for three-way GO tests |
| 6 | Only Q1 gets GO analysis | Run GO on all 4 MeCP2/H2AK119Ub quadrants (filtered to mCH-significant genes) |
| 7 | mCH correlations computed against FDR-weighted FCs | Fixed — reads v2 output which uses median FCs |
| 8 | Q1 label position wrong (0.10\*ymax) | Fix to 0.90\*ymax |
| 9 | Hardcoded paths | here::here() relative paths |
| 10 | Promoter excluded from TxDb filter | N/A — v2 upstream already includes Promoter |
| 11 | ggsave instead of save_multiformat_ggplot | Use multi_format_output.R |

## Decisions (from user)

| Question | Decision |
|----------|----------|
| Upstream source | Read v2 output CSV (no code duplication) |
| Duplicate mCH genes (23 entries) | Keep highest \|logFC\| per gene_name |
| Three-way GO enrichment | Yes — per MeCP2/H2AK119Ub quadrant, restricted to mCH-significant genes |
| Gencode BED | Use for gene body coordinates (cross-reference after join) |
| GO gene sets | Per MeCP2/H2A quadrant + mCH significance (FDR < 0.05) |
| Third plot | Yes — MeCP2 vs H2AK119Ub scatter colored by mCH status/direction |
| mCH significance definition | FDR only (edger_fdr < 0.05), no FC threshold |
| mCH diagnostics | Both gene-length and coverage diagnostics |

## Deliverables

One new R script: `mecp2_h2aub_mch/mecp2_h2aub_mch_v2.R`

New file — additive. The three original scripts remain for reference.

One new file: `mecp2_h2aub_mch/README.txt` — documents methodology.

## Script Structure

### 1. Setup (~lines 1-30)

```r
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(GenomicRanges)
library(IRanges)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(goseq)

source(here("scripts", "utils", "multi_format_output.R"))

out_dir <- here("mecp2_h2aub_mch")
```

Note: Does not load ChIPseeker or TxDb — the upstream annotation is handled by v2.

### 2. Load upstream v2 output (~lines 31-50)

Read the corrected gene-level MeCP2/H2AK119Ub results:

```r
gene_df <- read_csv(here("adult_mecp2_vs_h2aub", "gene_level_results_v2.csv"))
```

This table has columns: SYMBOL, n_mecp2_peaks, n_mecp2_sig_peaks, mecp2_has_sig, mecp2_Fold, mecp2_min_FDR, n_h2a_peaks, n_h2a_sig_peaks, h2a_has_sig, h2a_Fold, h2a_min_FDR, sig, quadrant, go_eligible, has_physical_overlap, n_total_peaks, gene_length.

Report the number of genes and significant genes.

### 3. Load and deduplicate mCH edgeR output (~lines 51-90)

```r
mch <- read_tsv(here("results", "03_differential", "mch_differential_results.tsv"))
```

The mCH output has 21,097 rows (20,369 testable + header) with 23 duplicate gene_names (same gene_name, different ENSMUSG IDs from overlapping gene models). Deduplicate by keeping the row with the highest absolute edger_logFC per gene_name:

```r
mch <- mch %>%
  group_by(gene_name) %>%
  slice_max(abs(edger_logFC), n = 1, with_ties = FALSE) %>%
  ungroup()
```

Report how many genes before and after dedup, and list the dropped duplicates to console.

Key mCH columns used downstream: gene_name, gene_id, edger_logFC, edger_fdr, sig_fdr005, total_coverage_mean, gene_length (as mch_gene_length to avoid collision with v2's gene_length).

### 4. Load gencode gene body coordinates (~lines 91-110)

```r
gencode_bed <- read.table(
  here("mecp2_h2aub_mch", "gencode.vM25.mouse.genes.annotation.bed"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
```

Columns: Chromosome (numeric, no "chr"), Start, End, Annotation, Name. Convert to GRanges with `chr` prefix for cross-reference.

### 5. Merge (~lines 111-140)

Inner join v2 output with deduplicated mCH output:

```r
merged <- inner_join(gene_df, mch, by = c("SYMBOL" = "gene_name"))
```

Left join gencode coordinates for gene body positions:

```r
merged <- left_join(merged, gencode_coords, by = c("SYMBOL" = "Name"))
```

Report: number of genes with all three modalities, overlap breakdown.

### 6. Significance flags (~lines 141-170)

Define combined significance for the two scatter plots:

- `mch_sig`: edger_fdr < 0.05 (boolean)
- `mch_sig_vs_mecp2`: mch_sig AND mecp2_has_sig — "Significant" / "Not significant" factor
- `mch_sig_vs_h2a`: mch_sig AND h2a_has_sig — same
- `mch_direction`: "mCH down" / "mCH up" / "mCH not sig" — for three-way scatter coloring

### 7. Quadrant assignment (~lines 171-200)

For the mCH scatter plots (mCH on y-axis):
- `quadrant_mecp2`: sign(mecp2_Fold) x sign(edger_logFC)
- `quadrant_h2a`: sign(h2a_Fold) x sign(edger_logFC)

Use a shared helper function `assign_quadrant(x_fold, y_fold)` returning standard labels.

### 8. Correlation tests (~lines 201-230)

Spearman correlations:
1. mCH logFC vs MeCP2 Fold (all merged genes)
2. mCH logFC vs H2AK119Ub Fold (all merged genes)

Print results and format rho labels for plot subtitles.

### 9. Quadrant statistics and label placement (~lines 231-280)

Shared function `build_quadrant_counts(df, quadrant_col, sig_col, x_fold_col, y_fold_col)` — same structure as runsecond but with corrected Q1 label position:

```r
quadrant_counts$y <- c(
  0.90 * ymax,   # Q1 — FIXED from 0.10
  0.95 * ymax,   # Q2
  -0.90 * ymax,  # Q3
  -0.90 * ymax   # Q4
)
```

### 10. Scatter plot: mCH vs MeCP2 (~lines 281-340)

X-axis: MeCP2 log2(FC, BAP1-KO / ctrl)
Y-axis: mCH log2FC (edgeR, BAP1-KO / ctrl)
Color: mch_sig_vs_mecp2 (grey/deeppink)
Labels: top 10 significant genes from Q1 + top 10 from Q3 (by |mecp2_Fold| + |edger_logFC|)
Quadrant count annotations with n and % per group, plus median annotated.
Subtitle: Spearman rho, total genes, significant gene count.

### 11. Scatter plot: mCH vs H2AK119Ub (~lines 341-400)

Same structure as section 10, with h2a_Fold on x-axis and mch_sig_vs_h2a for coloring.

### 12. Three-way scatter: MeCP2 vs H2AK119Ub colored by mCH (~lines 401-460)

X-axis: MeCP2 log2(FC, BAP1-KO / ctrl)
Y-axis: H2AK119Ub log2(FC, BAP1-KO / ctrl)
Color: three-level factor — grey = mCH not sig, blue = mCH sig and down (logFC < 0), red = mCH sig and up (logFC > 0)
Quadrant labels show total genes AND mCH-sig gene count per quadrant.
Subtitle: total merged genes, how many mCH-sig, how many mCH-sig-down vs up.

### 13. Diagnostic: gene length vs mCH significance (~lines 461-500)

Violin + boxplot of gene_length (from mCH output) by mCH significance status (sig_fdr005).
Wilcoxon test, annotate with n per group and median values on the plot.

### 14. Diagnostic: coverage vs mCH significance (~lines 501-540)

Violin + boxplot of total_coverage_mean by mCH significance status.
Wilcoxon test, annotate with n per group and median values on the plot.

### 15. Three-way GO enrichment (~lines 541-650)

Universe: all genes in merged table (genes with MeCP2 + H2AK119Ub peaks + testable mCH).

For each MeCP2/H2AK119Ub quadrant with >= 10 mCH-significant genes:

**Gene set**: genes that are (a) significant in both MeCP2 and H2AK119Ub (from v2's `sig` column) AND (b) mCH-significant (edger_fdr < 0.05), in that quadrant.

**enrichGO** (corrected background):
```r
universe_entrez <- bitr(merged$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
enrichGO(gene = q_entrez, universe = universe_entrez$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)
```
Run BP, MF, CC for each qualifying quadrant.

**goseq** (gene-length correction):
Use gene_length from the merged table (mCH pipeline's gene body length).
```r
gene_vector <- as.integer(merged$SYMBOL %in% q_symbols)
names(gene_vector) <- merged_entrez$ENTREZID
pwf <- nullp(gene_vector, bias.data = gene_lengths_matched)
goseq_result <- goseq(pwf, "mm10", "knownGene")
```

Save enrichGO dotplots and goseq results per quadrant.

### 16. Output (~lines 651-730)

All output to `mecp2_h2aub_mch/` via `here::here()`.

**CSVs:**
- `mch_mecp2_h2aub_merged_v2.csv` — full merged table
- `quadrant_counts_mecp2_v2.csv` — mCH vs MeCP2 quadrant stats
- `quadrant_counts_h2a_v2.csv` — mCH vs H2AK119Ub quadrant stats
- `three_way_quadrant_summary.csv` — MeCP2/H2AK119Ub quadrants with mCH breakdown
- `GO_three_way_{Q1,...}_{BP,MF,CC}.csv` — enrichGO results per quadrant
- `GO_goseq_three_way_{Q1,...}.csv` — goseq results per quadrant

**Plots (all via save_multiformat_ggplot → PDF/SVG/PNG/JPEG per figure):**
- `mch_vs_mecp2_v2/` — mCH vs MeCP2 scatter
- `mch_vs_h2aub_v2/` — mCH vs H2AK119Ub scatter
- `three_way_scatter_v2/` — MeCP2 vs H2AK119Ub colored by mCH status
- `diagnostic_mch_genelength/` — gene length vs mCH significance
- `diagnostic_mch_coverage/` — coverage vs mCH significance
- `GO_three_way_{Q1,...}_{ont}_dotplot/` — GO dotplots

### 17. README

New `mecp2_h2aub_mch/README.txt` documenting:
- That upstream MeCP2/H2AK119Ub data comes from `adult_mecp2_vs_h2aub/adult_MeCP2vsH2AUb_v2.R` (corrected aggregation)
- mCH edgeR deduplication method (highest |logFC| per gene_name)
- Three-way significance definition
- GO methodology (enrichGO + goseq, corrected universe)
- Gencode BED cross-reference for coordinates

## Critical files

| File | Action |
|------|--------|
| `mecp2_h2aub_mch/mecp2_h2aub_mch_v2.R` | CREATE — new unified script |
| `mecp2_h2aub_mch/README.txt` | CREATE — methodology documentation |
| `mecp2_h2aub_mch/mecp2_h2aub_all_elements_txdb.R` | UNCHANGED — kept for reference |
| `mecp2_h2aub_mch/runfirst_mecp2_h2aub_genebodies_gencode.R` | UNCHANGED — kept for reference |
| `mecp2_h2aub_mch/runsecond_mecp2_h2aub_noncg_gencode.R` | UNCHANGED — kept for reference |

## Dependencies

| Input file | Source |
|------------|--------|
| `adult_mecp2_vs_h2aub/gene_level_results_v2.csv` | v2 script output (must be run first) |
| `results/03_differential/mch_differential_results.tsv` | edgeR pipeline step 03 |
| `mecp2_h2aub_mch/gencode.vM25.mouse.genes.annotation.bed` | Gencode vM25 gene bodies |

## Reusable patterns

- `scripts/utils/multi_format_output.R` — `save_multiformat_ggplot(plot, base_path, width, height)` → PDF/SVG/PNG/JPEG per figure
- `here::here()` for all paths
- All comparison plots annotated with n per group and summary values (median) per memory feedback

## Verification

After writing the script, hand the user a complete `Rscript` command with `2>&1 | tee`. The user runs it and checks:

1. No duplicate genes in merged output — row count matches unique gene count
2. Spearman correlations computed against corrected median FCs (not FDR-weighted)
3. Q1 label sits in upper-right corner of all scatter plots
4. Three-way scatter shows mCH direction coloring
5. GO universe size matches merged gene count, not full genome
6. goseq results exist alongside enrichGO
7. Diagnostic plots show n and median per group
8. All plots saved via multi_format_output in per-figure subdirectories
9. Axis labels say "BAP1-KO / ctrl"
