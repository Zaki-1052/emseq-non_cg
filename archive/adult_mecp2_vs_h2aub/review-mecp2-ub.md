Good, I've now verified the data outputs against the code. Here's the full critique.

---

## Chat

### Confirmed Bugs

**1. `n_mecp2_sig_peaks` and `n_h2a_sig_peaks` are broken — always 0 or 1**

This is the most concrete code bug. In `adult_MeCP2vsH2AUb.R:177–198`, the `summarise()` call redefines `mecp2_FDR` as a scalar (`min(mecp2_FDR)`) on line 186 before `n_mecp2_sig_peaks = sum(mecp2_FDR < 0.05)` on line 195. In dplyr, `summarise` evaluates sequentially and the later expression sees the already-computed scalar, not the original per-peak vector. So `sum(scalar < 0.05)` is always 0 or 1.

I confirmed this in the actual output: `n_mecp2_sig_peaks` has only values {0, 1} across all 8,276 genes, and `n_h2a_sig_peaks` is identical ({0, 1} only). Look at Etv1 in the top-genes CSV — it has 36 MeCP2 peaks with a gene-level min FDR of 3.55e-39, yet `n_mecp2_sig_peaks = 1`. It's nearly certain that many of those 36 peaks are individually significant.

`mecp2_sig = any(mecp2_FDR < 0.05)` happens to give the correct answer by coincidence (`any(min < 0.05)` ≡ `any(original_vector < 0.05)`), but the `n_*_sig_peaks` columns are wrong and should not be reported.

**Fix:** Either reorder the columns so `n_*_sig_peaks` is computed before `*_FDR` is overwritten, or give the aggregated FDR a different name (e.g., `mecp2_minFDR`).

**2. README says gene-body filtering was applied; code doesn't implement it**

README line 3: *"If not in gene body, it was not used."* The code stores the `annotation` and `distanceToTSS` columns from ChIPseeker but only filters `!is.na(SYMBOL)` (lines 149–154). No annotation-based filter is applied. Peaks assigned to "Distal Intergenic", "Downstream (>300kb)", etc. are all retained.

Either the README is aspirational documentation that wasn't implemented, or the filtering was meant to happen in the code and was missed. If the input files were pre-curated to gene body only, the code should assert that rather than silently assuming it.

**3. Stale comment in overlap script**

`log2log2plot_physicaloverlap.R:40` says *"require at least 80% of the MeCP2 peak"* but the actual filter on line 41 uses `>= 0.5` (50%). The plot title correctly says 50%. The comment misleads anyone reading just the code.

---

### Statistical Methodology Issues

**4. Self-reinforcing weight creates a multiplicative bias toward extreme peaks**

The core aggregation formula:

```
weight = |FC| × −log10(FDR)
Gene FC = Σ(FC × weight) / Σ(weight)
```

FC appears in both the value and the weight. Expanding the numerator for positive FC peaks: `Σ(FC² × −log10(FDR))`. This is not a weighted mean of fold changes — it's a ratio where extreme peaks contribute quadratically. A peak with FC=5 gets 5× the weight of FC=1 *and* contributes a 5× larger value. The result is dominated by whichever single peak has the largest |FC|, making the multi-peak aggregation largely cosmetic.

Concrete example with two peaks for one gene:
|             | FC  | FDR   | weight | FC × weight          |
| ----------- | --- | ----- | ------ | -------------------- |
| Peak A      | 5.0 | 0.001 | 15.0   | 75.0                 |
| Peak B      | 0.5 | 0.01  | 1.0    | 0.5                  |
| **Gene FC** |     |       |        | **75.5 / 16 = 4.72** |

Simple average would be 2.75. Weighting by `−log10(FDR)` alone would give 3.2. The self-reinforcing weight pushes it to 4.72 — essentially recovering Peak A's value. The aggregation is performing "pick the most extreme peak" with extra steps.

**Recommendation:** Weight by `−log10(FDR)` alone (statistical confidence without circularity), or by peak width, or by concentration. Any of these decouple the weight from the value being averaged.

**5. FDR floor of 1e-300 creates enormous dynamic range in weights**

Line 166: `pmax(mecp2_FDR, 1e-300)` — this floors FDR to avoid `log10(0)`. But `−log10(1e-300) = 300`, while `−log10(0.05) ≈ 1.3`. A single peak with FDR ≈ 0 gets a weight ~230× larger than a barely-significant peak, completely dominating the gene-level value. Combined with issue #4, this means one very significant peak with a large FC overwhelms all others.

**Recommendation:** Floor at the smallest observed non-zero FDR in the dataset, or at something biologically reasonable like 1e-10 or 1e-20.

**6. GO enrichment uses the wrong background**

`enrichGO()` is called without the `universe` argument (lines 647–656). This defaults to all annotated genes in org.Mm.eg.db (~28,900 genes per the BgRatio column in the output). The correct background is the 8,276 genes that could have appeared in the analysis (genes with peaks in both datasets). Using all genome genes inflates enrichment because:

- Genes with peaks tend to be expressed, active, and longer
- Neuronal/synaptic genes are among the longest in the genome and disproportionately represented among genes-with-peaks
- The enrichment is testing "are Q1 genes enriched for synapse terms vs. all genes?" when it should be testing "vs. genes that had peaks in both marks"

The GO results (synapse organization, axon guidance, cognition, learning) are biologically expected for MeCP2, which makes this harder to catch — but the p-values and fold enrichments are likely overstated.

**Recommendation:** Pass `universe = unique(gene_df$SYMBOL)` converted to Entrez IDs. Better yet, use `goseq` or `GREAT`, which can correct for gene length bias — a critical confound here since longer genes accumulate more peaks.

**7. No gene-level multiple testing correction for the "significant" label**

A gene is labeled "Significant" if it has ≥1 peak with FDR < 0.05 in both MeCP2 *and* H2AUb. But genes with more peaks get more chances. A gene with 36 peaks (like Etv1) has 36 independent shots at having one pass the threshold; a gene with 1 peak gets one shot. This inflates significance for long, peak-dense genes — the same class that drives the GO enrichment.

This isn't a standard multiple-testing problem with a clean correction, but it should be acknowledged and tested for. One approach: check whether `n_peaks` correlates with `sig` status after controlling for fold change.

**8. No fold-change threshold — genes with trivial FC are included**

A gene with FC = 0.001 in both marks lands in Q1 as "MeCP2 up / H2AK119Ub up." 2,724 genes are in Q1 (32.9%), but many of those likely have near-zero fold changes clustered at the origin. The quadrant counts and GO analysis don't distinguish between biologically meaningful changes and noise near zero. Adding a minimum |FC| filter (e.g., 0.5 or 1.0) or using a distance-from-origin metric would focus the GO analysis on genes with real effect sizes.

---

### Biological Analysis Design

**9. Gene length bias is unaddressed and likely confounds the GO results**

This is related to #6 but worth stating as a separate biological concern. ChIP-seq peak density scales with gene length. Longer genes are more likely to have peaks in both datasets (surviving the inner join), more likely to have ≥1 significant peak in each (meeting the sig threshold), and more likely to have extreme weighted FCs (more peaks feeding the self-reinforcing weight). Neuronal genes (cadherins, contactins, SLITRKs, NRXNs, etc.) are among the longest protein-coding genes in mammals. The entire analysis pipeline — from gene assignment through significance through GO — preferentially surfaces long genes. The GO enrichment may be partly or largely a gene-length artifact.

This is testable: plot gene length vs. significance status, and run `goseq` (which models length bias) as a comparison to `enrichGO`.

**10. Only Q1 gets GO analysis**

76% of significant genes fall in Q1 (both marks up in mutant), which makes it the obvious focus. But Q3 (both down, 285 sig genes / 14%) and Q2 (MeCP2 down / H2A up, 121 sig genes / 5.9%) may reveal distinct biology. Q3 genes — where both marks decrease in the mutant — could identify loci losing both MeCP2 and Polycomb silencing. Q2 genes — MeCP2 down but H2AK119Ub up — could point to a compensatory Polycomb response. Even a supplementary GO analysis on these quadrants would strengthen the story.

**11. The two scripts answer different questions but aren't connected**

`adult_MeCP2vsH2AUb.R` asks: "for genes with peaks in both marks, do the fold changes correlate?" `log2log2plot_physicaloverlap.R` asks: "at loci where the peaks physically overlap, do the fold changes correlate?" These are complementary views — gene-level convergence vs. locus-level co-occupancy — but neither script references the other's output. Understanding which Q1 *genes* also have physically overlapping peaks would directly test whether the concordance is from co-occupancy or from independent regulation at different positions within the same gene.

**12. Biological question about the experimental design**

The axis labels say "mutant/control." If this is a MeCP2 loss-of-function mutant (Rett model), then the x-axis is measuring binding of the *mutant* MeCP2 protein, which complicates interpretation — gain of binding at some loci could reflect altered DNA-binding specificity of the mutant protein rather than increased abundance. If it's an overexpression model (MECP2 duplication), the Q1 enrichment (more MeCP2 → more H2AK119Ub) has a cleaner interpretation. The manuscript should address what the mutant is and what "MeCP2 fold change in a MeCP2 mutant" means mechanistically.

---

### Code Quality

**13. Hardcoded paths to another user's machine**

Lines 36, 43, and 705 reference `/Users/janviramchandra/Desktop/`. This will fail for anyone else without modification. Use `here::here()` or `file.path()` with a configurable base directory.

**14. `merged$column` inside dplyr pipe**

`log2log2plot_physicaloverlap.R:41`:
```r
merged <- merged %>% filter(merged$overlap_fraction_of_mecp2 >= 0.5 ...)
```
Using `merged$` inside a pipe bypasses the dplyr data mask. It works here but is fragile — if `merged` were updated earlier in a longer pipe chain, the `$` would still reference the old version. Should be just `filter(overlap_fraction_of_mecp2 >= 0.5 & overlap_fraction_of_h2a >= 0.5)`.

**15. Q1 label position is misplaced on the plot**

In line 431, `quadrant_counts$y[1] = 0.10*ymax` places the Q1 label near the x-axis rather than in the upper-right quadrant where Q1 data lives. The plot confirms this — the Q1 stats text sits at the bottom of Q1, overlapping the dense data region. Compare to the overlap script, which places Q1 at `0.97*ymax`.

---

### Summary — what I'd prioritize for the grad student

| Priority               | Issue                                         | Impact                                                  |
| ---------------------- | --------------------------------------------- | ------------------------------------------------------- |
| Fix now                | #1: `n_*_sig_peaks` bug                       | Reported values are wrong                               |
| Fix now                | #2: gene-body filter not implemented          | Results may include peaks the README says were excluded |
| Fix before publication | #4: self-reinforcing weight                   | Gene-level FCs are biased toward extremes               |
| Fix before publication | #6: GO background                             | Enrichment p-values likely inflated                     |
| Investigate            | #9: gene-length bias                          | Could confound the entire GO story                      |
| Fix before publication | #5: FDR floor                                 | Weights have absurd dynamic range                       |
| Discuss                | #7, #8: sig criteria, FC threshold            | Affects which genes enter GO                            |
| Discuss                | #10, #11: other quadrants, script integration | Missed biology                                          |
| Quick fix              | #3, #13, #14, #15                             | Code hygiene                                            |

The GO results (synapse, cognition, learning) are biologically plausible for MeCP2 — which is exactly what makes the statistical concerns harder to catch. The results may well be real, but the current pipeline can't distinguish "real enrichment" from "long-gene artifact feeding a self-reinforcing weight into an improperly-backgrounded GO test." Fixing #4, #6, and testing for #9 would answer that.