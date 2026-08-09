# Plan: Rewrite MeCP2 vs H2AK119Ub Quadrant Analysis

## Context

Janvi's MeCP2 vs H2AK119Ub quadrant analysis (`adult_mecp2_vs_h2aub/`) has 15 issues identified in Jai's code review (`jai-review.md`). The analysis takes DiffBind peak output for MeCP2 and H2AK119Ub ChIP-seq, annotates peaks to genes, collapses multiple peaks per gene into a single fold-change value, builds a quadrant scatter plot of gene-level MeCP2 FC vs H2AK119Ub FC, and runs GO enrichment on significant concordant genes.

The most serious problems: (1) a dplyr sequential-evaluation bug making `n_*_sig_peaks` always 0 or 1, (2) gene-body filtering described in the README but never implemented, (3) FDR used as a continuous weight (statistically unsound — FDR is a binary gatekeeper per the ASA 2016 statement, not a quantitative score), (4) wrong GO background inflating enrichment, and (5) no gene-length bias correction.

This is a BAP1-KO model (Math1-Cre conditional knockout in cerebellum), not a MeCP2 mutant. The axis labels should say "MeCP2 FC in BAP1 KO / ctrl."

## Decisions (from user)

| Issue | Decision |
|-------|----------|
| Gene-body filter (#2) | Implement: keep only Promoter, Intron, Exon, 5'/3' UTR annotations |
| Weighting scheme (#4) | Remove FDR from weights entirely. FDR is a binary filter, not a score |
| Aggregation method | Median FC of significant peaks per gene (robust to outlier peaks) |
| FDR floor (#5) | N/A — FDR no longer used as weight |
| FC threshold (#8) | Require \|log2FC\| >= 0.5 in both marks for GO entry |
| GO background (#6) | Fix to genes-with-peaks-in-both-marks |
| Gene-length bias (#9) | Run both enrichGO (corrected background) AND goseq (length correction) |
| GO scope (#10) | All four quadrants |
| Script integration (#11) | Single unified script merging gene-level + physical-overlap analyses |
| Mutant identity (#12) | BAP1-KO — update axis labels |
| File paths (#13) | here::here() relative paths |
| Peak-count bias (#7) | Diagnostic plot (n_peaks vs sig status), no correction |
| Q1 label position (#15) | Fix to upper-right quadrant |
| Stale comment (#3) | Fix (comment says 80%, code uses 50%) |
| merged$ in pipe (#14) | Fix to bare column names |

## Deliverables

One new R script: `adult_mecp2_vs_h2aub/adult_MeCP2vsH2AUb_v2.R`

New file — additive, does not modify or delete the original script. The original remains as-is for reference.

One updated file: `adult_mecp2_vs_h2aub/README.txt` — rewritten to match the new methodology.

## Script Structure

The unified script has these sections in order:

### 1. Setup and input (~lines 1–40)

- `library()` calls: GenomicRanges, IRanges, ggplot2, dplyr, tidyr, ChIPseeker, TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db, ggrepel, clusterProfiler, enrichplot, DOSE, goseq
- Use `here::here("adult_mecp2_vs_h2aub", "mecp2.txt")` and same for `h2aub.txt`
- `read.table()` as current, unchanged column expectations (seqnames, start, end, width, strand, Conc, Fold, p.value, FDR)
- Output directory: `here::here("adult_mecp2_vs_h2aub")`

### 2. GRanges conversion and peak annotation (~lines 41–90)

- Convert both peak tables to GRanges, attach Fold/p.value/FDR as mcols (unchanged from current)
- `annotatePeak()` against mm10 knownGene TxDb (unchanged)
- Build annotation data frames with SYMBOL, geneID, annotation, distanceToTSS, Fold, FDR, peak width

### 3. Gene-body filtering (NEW — fixes #2)

- Define gene-body annotations: `c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron")`
- ChIPseeker annotations include parenthetical detail (e.g., "Promoter (<=1kb)", "Intron (ENST...)"). Use `grepl()` to match the leading term
- Filter both `mecp2_df` and `h2a_df` to rows where annotation matches any gene-body category
- Also filter `!is.na(SYMBOL)` as before
- Report how many peaks survive filtering (cat to console)

### 4. Gene-level collapse (fixes #1, #4, #5)

**The aggregation change is the most important fix in the rewrite.**

Current (broken): `weight = |FC| * -log10(FDR)`, weighted mean of FC, then `n_sig_peaks = sum(FDR < 0.05)` after FDR has been overwritten to a scalar.

New:
```
gene_level <- peaks_df %>%
  group_by(SYMBOL) %>%
  summarise(
    n_total_peaks = n(),
    n_sig_peaks = sum(FDR < 0.05),           # computed FIRST, before any aggregation
    has_sig_peak = any(FDR < 0.05),
    gene_Fold = median(Fold[FDR < 0.05]),     # median FC of significant peaks only
    min_FDR = min(FDR),                        # renamed to avoid overwrite bug
    .groups = "drop"
  )
```

Key changes:
- `n_sig_peaks` is computed before any column is overwritten (fixes #1)
- `gene_Fold` uses `median()` of significant peaks only — no FDR weighting at all
- For genes with zero significant peaks, `gene_Fold` will be `NA` — these genes won't pass the significance filter anyway, so this is correct
- `min_FDR` uses a distinct name from the per-peak `FDR` column to prevent the dplyr sequential-evaluation trap
- No FDR floor needed — FDR is not used as a weight

### 5. Merge and significance (~lines ~130–170)

- Inner join MeCP2 and H2AUb gene-level tables by SYMBOL (unchanged logic)
- Significance: gene has `has_sig_peak == TRUE` in BOTH marks (same definition, just using the pre-computed boolean)
- Report total genes and significant gene counts

### 6. Quadrant assignment and correlation (~lines ~170–210)

- Same quadrant logic (sign of FC in each mark)
- Spearman correlation on gene-level FCs (unchanged)
- Quadrant statistics: total and significant genes per quadrant with percentages

### 7. FC threshold for GO analysis (NEW — addresses #8)

- After quadrant assignment, flag genes for GO eligibility: `go_eligible = sig == "Significant" & abs(mecp2_Fold) >= 0.5 & abs(h2a_Fold) >= 0.5`
- Report how many significant genes pass vs. fail the FC filter

### 8. Physical overlap analysis (merged from overlap script — addresses #11)

- Convert MeCP2 and H2AUb peaks to GRanges (already done in step 2 — reuse)
- `findOverlaps()` between the two peak sets
- Compute overlap fractions for both directions (unchanged from overlap script)
- Filter to >= 50% reciprocal overlap (fixes #3: remove stale "80%" comment, or just write correct comment)
- Fix #14: use bare column names in dplyr pipes, not `merged$column`
- Assign significance and quadrants to overlapping loci (same as overlap script)
- Build the overlap-specific scatter plot
- Add a cross-reference column to gene_df: for each gene, does it have at least one physically overlapping peak pair?

### 9. Diagnostic: peak count vs. significance (NEW — addresses #7)

- Boxplot or violin plot: `n_total_peaks` (combined MeCP2 + H2AUb) on y-axis, significance status on x-axis
- Wilcoxon test for the difference
- This goes into the output as a diagnostic figure, not used for filtering

### 10. Diagnostic: gene length vs. significance (NEW — addresses #9)

- Retrieve gene lengths from TxDb (`genes()` and `width()`)
- Boxplot: gene length by significance status
- Wilcoxon test
- This supports interpretation of GO results

### 11. Quadrant scatter plot (fixes #15)

- Same ggplot structure as current but:
  - Fix Q1 label position: place at `(0.90*xmax, 0.90*ymax)` instead of `(0.90*xmax, 0.10*ymax)`
  - Axis labels: "MeCP2 log2(FC, BAP1-KO / ctrl)" and "H2AK119Ub log2(FC, BAP1-KO / ctrl)"
  - Subtitle includes Spearman rho, total genes, significant genes, and the FC threshold note
- Top gene labels: same logic (top 10 overall + top 5 per non-Q1 quadrant)

### 12. GO enrichment — all four quadrants (fixes #6, #9, #10)

For each quadrant with >= 10 GO-eligible significant genes:

**enrichGO (corrected background):**
```
universe_entrez <- bitr(unique(gene_df$SYMBOL), ...)$ENTREZID
enrichGO(gene = q_entrez, universe = universe_entrez, ...)
```
- Runs BP, MF, CC for each qualifying quadrant
- Background is all genes in gene_df (genes with peaks in both marks), not the full org.Mm.eg.db

**goseq (gene-length correction):**
```
gene_lengths <- # from TxDb
pwf <- nullp(gene_vector, bias.data = gene_lengths, ...)
goseq_result <- goseq(pwf, "mm10", "knownGene")
```
- Runs on the same gene sets as enrichGO
- Reports both results side by side for comparison

### 13. Output

All output goes to `adult_mecp2_vs_h2aub/` via `here::here()`.

**CSVs:**
- `gene_level_results_v2.csv` — full gene-level table with corrected n_sig_peaks, median FC, overlap flag
- `quadrant_counts_v2.csv` — quadrant statistics
- `top_labeled_genes_v2.csv` — genes labeled on plot
- `GO_enrichGO_{Q1,Q2,Q3,Q4}_{BP,MF,CC}.csv` — enrichGO results per quadrant per ontology (only for quadrants with enough genes)
- `GO_goseq_{Q1,Q2,Q3,Q4}.csv` — goseq results per quadrant
- `overlap_loci_results.csv` — physical overlap analysis results

**Plots:**
- `quadrant_plot_v2.{png,svg}` — corrected gene-level scatter
- `overlap_plot_v2.{png,svg}` — physical overlap scatter (from merged overlap script)
- `diagnostic_peak_count_vs_sig.{png,svg}` — peak count bias diagnostic
- `diagnostic_gene_length_vs_sig.{png,svg}` — gene length bias diagnostic
- `GO_enrichGO_{Q1,...}_dotplot.png` — dotplots per quadrant
- `GO_goseq_{Q1,...}_dotplot.png` — goseq dotplots per quadrant

### 14. Updated README

Rewrite `README.txt` to match the new methodology:
- Document the median-of-significant-peaks aggregation (no FDR weighting) and cite the ASA 2016 statement as the rationale
- Document gene-body filtering (ChIPseeker annotation categories retained)
- Document the FC threshold for GO eligibility
- Document the corrected GO background and goseq length correction
- Document the physical overlap integration

## Critical files

| File | Action |
|------|--------|
| `adult_mecp2_vs_h2aub/adult_MeCP2vsH2AUb_v2.R` | CREATE — new unified script |
| `adult_mecp2_vs_h2aub/README.txt` | EDIT — rewrite to match new methods |
| `adult_mecp2_vs_h2aub/adult_MeCP2vsH2AUb.R` | UNCHANGED — original kept for reference |
| `adult_mecp2_vs_h2aub/log2log2plot_physicaloverlap.R` | UNCHANGED — merged into v2, original kept |

## Reusable patterns from this repo

- `scripts/utils/multi_format_output.R` — `save_multiformat_ggplot(plot, base_path, width, height)` saves plots in PDF/SVG/PNG/JPEG into per-figure subdirectories. Source via `here::here("scripts/utils/multi_format_output.R")`. All plots use this utility instead of inline ggsave calls.

## Verification

After writing the script, hand the user a complete run command. The user runs it in RStudio or terminal and checks:

1. `n_mecp2_sig_peaks` in the output CSV has values > 1 for multi-peak genes (bug #1 fixed)
2. Gene count is lower than 8,276 (gene-body filter is removing distal intergenic peaks)
3. GO results use a BgRatio denominator matching the gene_df size, not ~28,900
4. goseq results exist alongside enrichGO for comparison
5. Q1 label sits in the upper-right corner of the plot (not near the x-axis)
6. Axis labels say "BAP1-KO / ctrl"
7. Diagnostic plots show peak-count and gene-length distributions by significance status
8. Overlap analysis results appear in the output with a cross-reference to gene-level data
