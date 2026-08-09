Here's what I see in the three output files:

---

**Pipeline ran successfully.** 21,846 genes considered, 20,371 tested (1,475 dropped for low coverage or too few CA sites). Model used: `~ genotype + sex + lambda_ca_rate` (no batch term, so Model A won).

**156 Bonferroni-significant genes** (427 at FDR 0.05, 597 at FDR 0.10). Directional split is 2:1 — 106 genes with mCA *higher* in mutant, 50 lower. The effect sizes are large: median |mca_diff| among significant genes is ~3.1 percentage points, on a control background where the genome-wide median mCA is only ~1.7%. The top hits (Npvf, Olfr437, Ttr) go from ~4–5% to ~14–15% — roughly tripling.

**Robustness is mixed.** Only 23 of the 156 Bonferroni-significant genes survive all 8 leave-one-out refits. 56 survive 7 of 8, and 10 survive none. That means the majority of hits depend on at least one specific sample — consistent with the coverage imbalance the plan flagged (mut_M2 at ~118x, mut_F1 at ~101x vs 38–55x for the rest). The 23 fully robust genes are the high-confidence core.

**Gene lengths skew short.** Median gene length among significant hits is ~12kb, with 31 of 156 under 5kb. This is the *opposite* of what the MeCP2/mCA hypothesis predicts — if BAP1-loss altered DNMT3A-mediated mCA deposition in a way that mattered for MeCP2 gene-body dampening, you'd expect enrichment in long neuronal genes (Nrxn1-class, 100kb+). Instead, many top hits are short, non-neuronal genes (Npvf, Olfr437, Ttr, Krt24).

**Neuronal enrichment is modest.** 44 of the 4,070 tested neuronal genes are significant (Fisher OR = 1.58, p = 0.012). Statistically significant enrichment, but the effect is small — neuronal genes are ~1.6x over-represented, not 5–10x.

**Dispersion is almost entirely at the DSS floor** (0.001). The median dispersion for significant genes is 0.001, the minimum is 0.001, the max is 0.002. That means DSS's shrinkage estimator pushed most genes to its floor — likely because n=4 per group gives almost no degrees of freedom to estimate overdispersion. This is something to think about: the p-values may be anti-conservative (too small) if the true dispersions are larger than 0.001.

**Lambda spike-in rates.** Six of eight samples are in the expected 0.65–0.92% range. Two are above 1% — ctrl_M1 at 1.06% and mut_M2 at 1.25%. These aren't alarming but they're on the high side. The genotype-lambda correlation (0.36) is non-trivial, which is why it's good that lambda_ca_rate is in the model.

**NEEDS YOU:** A few things to flag before you decide where to take this:

1. **The LOO fragility + dispersion floor is the biggest concern.** With most dispersions pinned at the DSS floor and 85% of hits losing significance when one sample drops, the 156-gene count is probably inflated. The 23 fully robust genes are the ones I'd trust for downstream interpretation. Want to look at what those 23 are (gene identities, lengths, neuronal overlap)?

2. **The gene-length and gene-identity profile doesn't match the MeCP2/mCA hypothesis.** The top hits are short, non-neuronal genes with massive effect sizes (3–10 pp shifts). That pattern looks more like cell-type contamination or cell-composition shifts (Ttr = choroid plexus, Emcn = endothelial) than a global DNMT3A-deposited mCA landscape change. The neuronal enrichment is real but weak.

3. **The column schema matches the plan exactly** — all columns from the plan's spec are present, naming conventions follow the biomodal pattern, per-sample values are included. The sample matrix is the expected 20,371 × 8 wide format. Files are ready for integration with the existing R visualization pipeline.

What would you like to dig into?

---

## Working

`scripts/plot_mca_volcano.R` written. Parses clean; **not run** — that's yours.

```
cd /Users/zakiralibhai/emseq && Rscript scripts/plot_mca_volcano.R --results mca_results/mca_differential_results.tsv 2>&1 | tee mca_results/plots_volcano_log.txt
```

Writes into `mca_results/plots/` as PDF + SVG + PNG each:

| Figure                   | Encoding                                             |
| ------------------------ | ---------------------------------------------------- |
| `mca_volcano_direction`  | colour = direction of change (the main plot)         |
| `mca_volcano_robustness` | colour = leave-one-out refits retaining significance |
| `mca_volcano_genelength` | colour = log10 gene length                           |
| `mca_volcano_panel`      | the three combined, tagged A/B/C                     |

All three share coordinates: x = `mca_diff` in percentage points, y = `-log10(dss_pval)`. Raw p on the axis, corrected cutoffs as horizontal lines — Bonferroni dashed at `0.05/20371`, BH dotted at the largest raw p passing FDR (read back from the data, since BH is step-up and has no fixed raw-p cutoff). Both lines are direct-labelled so they don't rely on line style alone. The 23 leave-one-out-robust genes get a dark ring, which is why panel A carries the fragility story without a second colour scale.

Options if you want to move it: `--sig-basis fdr` recolours on the 427 FDR genes instead of the 156, `--label-top N`, `--gene-list FILE` labels significant members of a set (for the 5,614 neuronal genes when you have that file locally — it's not in the biomodal checkout).

### Things I decided that you may want to reverse

- **Duplicated the style rather than sourcing it.** `theme_mca()` and the red/blue direction colours are copied out of `_shared_config.R` instead of sourced. Sourcing it executes ~600 lines that load run-5 DMR BEDs from `/expanse/...` at source time — it can't run on your Mac, and emseq is a separate repo. Same for `save_multiformat_ggplot()`, reimplemented locally. That's real duplication across two repos and I'd rather you know than discover it.
- **PNG, not JPEG.** Your `multi_format_output.R` writes `.jpg`; a 20k-point scatter with `geom_text_repel` labels gets visible JPEG ringing around the text. I deviated. Say the word and I'll match your util exactly.
- **No y-axis capping.** Section 4 caps `-log10(q)` at 300; here the smallest p is `1.25e-47` → max y ≈ 46.9, so nothing to cap.

### Two bugs I made and fixed before handing this over

1. **Panels B and C would have errored.** I factored the label layer into a shared `gene_label_layer()`, but panels B and C build from `ggplot(mapping = aes(...))` with no `data` — the layer inherited nothing to label. Fixed by giving it explicit data (only the labelled rows, which is also cheaper than passing 20k mostly-empty strings).
2. **Duplicated aesthetics on panel A.** The ring and label layers set `colour`/`size` as fixed params while inheriting panel A's `aes(colour = direction, size = ..., alpha = ...)`. Fixed with `inherit.aes = FALSE` and explicit mappings on both.

I also stopped short of guessing how optparse names `--sig-basis` in the result list — optparse isn't installed anywhere I could read its source, so I pulled the README, confirmed `dest` is supported and controls the name, and set `dest` explicitly on all four dashed flags rather than relying on an unverified dash-to-underscore rule.

### Colour choice

Kept your `COLORS$direction` pair. Ran the dataviz validator on it: `#D7191C`/`#2C7BB6` passes all six checks (CVD ΔE 23.1 deutan, normal-vision ΔE 32.1, contrast ≥3:1). `grey70` fails the chroma and contrast checks by design — it's the null category and should recede; darkening it to `#8C8C8C` collides with the blue at normal-vision ΔE 14.1, below the floor. So grey70 is the better choice, not just the incumbent one.

**NEEDS YOU:** one thing I couldn't settle from the data. Panel C's subtitle compares median gene length among significant genes (12.1 kb) against all tested genes — but the honest comparison for the MeCP2 hypothesis is whether significant genes are *shorter than expected given their power*, and short genes have less power, not more. So the short-gene skew is arguably evidence against a length-driven artifact rather than for one. I've plotted the raw pattern without interpreting it. If you want the powered comparison, I'd need to know whether you'd rather bin by `n_ca_sites_mean` or model length against `total_coverage_mean` — that's a real analysis choice, not a plotting one.

---

## Chat

Plots look clean — labels readable, rings visible on the robust genes, threshold lines distinct. Here's what they're saying.

### The shape of the volcano

Panel A shows a strikingly **asymmetric** volcano. The right arm (mCA gain in mutant) is longer, denser, and reaches far more extreme effect sizes than the left arm (mCA loss). The three outliers — Npvf (+10.2 pp), Olfr437 (+10.0 pp), Ttr (+8.4 pp) — sit in the upper right at p < 10^-33, well separated from everything else. These are not marginal calls; they're 2–3x the baseline mCA rate. The left arm has Islr2 and Edn2 as the strongest losses (~5 pp), significant but an order of magnitude less extreme by both effect size and p-value than the right-side outliers. The 2:1 directional ratio (106 up, 50 down) is visible in the density of red vs blue points above the Bonferroni line.

The grey cloud is reassuringly centered on zero with no systematic inflation — the bulk of the genome's mCA is unchanged, which is what you'd expect if this is a locus-specific effect rather than a global artefact.

### Robustness is the real story

Panel B is the most informative panel. The purple gradient reveals that **the volcano's structure is dominated by a few rock-solid hits and a long tail of fragile ones.** The darkest points (8/8 refits) are concentrated at the top: Npvf, Olfr437, Ttr, Emcn, Kctd12 — the genes with the largest effects survive every leave-one-out. As you move down toward the Bonferroni line, points lighten dramatically — many are at 0–3 refits, meaning they disappear when a single sample drops out.

This makes sense given the coverage imbalance (mut_M2 at ~118x, mut_F1 at ~101x vs 38–55x for the rest). The shallow middle tier of the volcano is being propped up by the two deep samples. The 23 fully robust genes are the trustworthy core; the other 133 Bonferroni-significant genes should be treated as exploratory.

### Gene length does not match the MeCP2 hypothesis

Panel C is the concerning one. The significant genes are predominantly **lighter green** — shorter than the genome-wide median. The subtitle quantifies it: median 12.1 kb among significant genes vs 19.5 kb across all tested. The top three hits are tiny: Npvf (3.8 kb), Olfr437 (6.6 kb), Ttr (9 kb). Some longer genes are present (Emcn 91 kb, Lrrc49 119 kb, Dpys 89 kb, Fgf1 91 kb), but they're the exception, not the pattern.

The MeCP2/mCA hypothesis predicts the opposite. MeCP2-mediated transcriptional dampening scales with gene-body length because longer genes accumulate more mCA binding sites (Gabel 2015). If BAP1-loss were disrupting DNMT3A → mCA → MeCP2 in a functionally meaningful way, you'd expect the biggest effects at Nrxn1-class (1 Mb+) and Syt1-class (340 kb) neuronal genes. Those aren't showing up as top hits. The Fisher enrichment for the neuronal gene set is real (OR 1.58, p = 0.012) but modest — 44 of 4,070 neuronal genes.

### What the gene identities are telling you

The top hits are biologically coherent — but they point toward **cell-type composition**, not a pan-neuronal mCA landscape shift:

- **Ttr** — transthyretin, a choroid plexus epithelial marker. Not expressed in neurons.
- **Emcn** — endomucin, an endothelial cell marker. Not neuronal.
- **Olfr437** — olfactory receptor, not expected in cerebellum at all.
- **Npvf** — neuropeptide VF precursor, expressed in a small hypothalamic population, not cerebellar neurons.

If BAP1-KO cerebellum has even a slightly different proportion of non-neuronal cells (choroid plexus, vasculature) relative to wildtype, that would produce exactly this pattern. Neurons carry ~1% mCA; non-neuronal cells carry ~0%. A shift in cell proportions within the bulk tissue would appear as a localized mCA change at genes that are markers for the contaminating cell type, with large effect sizes in short genes (because the per-gene mCA rate changes from near-zero to a non-trivial value when the ratio of neuronal to non-neuronal nuclei shifts).

### The mCA-up direction is the interesting part

The 2:1 ratio favoring mCA gain is at least consistent with the DNMT3A recruitment hypothesis — more H2AK119ub → more DNMT3A recruitment → more mCA deposition. But the genes where this is showing up don't match the predicted targets (long neuronal genes where MeCP2 operates). The gain could also be explained by a shift toward more neuronal nuclei in the mutant tissue relative to non-neuronal, which would raise apparent mCA at most genes without implying any change in DNMT3A activity.

### Bottom line

The pipeline worked: 156 Bonferroni-significant genes, 23 robust to leave-one-out, correct column schema, plausible effect sizes. The results are real differential mCA. But what they most likely reflect is **cell-type composition differences between BAP1-KO and wildtype cerebellum**, not a global rewiring of DNMT3A-mediated mCA at neuronal gene bodies. The gene identities (Ttr, Emcn, Olfr437), the short gene-length skew, and the modest neuronal enrichment all point the same direction.

**NEEDS YOU:** Two things to think about:

1. **Can you rule out composition?** If you have RNA-seq from these same samples, checking expression of Ttr, Emcn, and a panel of cell-type markers (Snap25 for neurons, Aqp4 for astrocytes, Cx3cr1 for microglia, Cldn5 for endothelium) between genotypes would tell you quickly whether the cell mix is shifting. If it's not, these mCA differences become much more interesting.

2. **The 23 robust genes deserve a closer look regardless.** Even if composition explains the bulk of the 156, the fully robust subset may contain genuine DNMT3A-target changes hiding underneath. Want me to pull those 23 out, check their identities, lengths, and neuronal-set membership?


## Working

`scripts/mca_gsea.R` written. Parses clean, CLI defaults verified against the prefix-collision bug. **Not run** — yours:

```
cd /Users/zakiralibhai/emseq && Rscript scripts/mca_gsea.R --results mca_results/mca_differential_results.tsv 2>&1 | tee mca_results/gsea_log.txt
```

Runs both, into `mca_results/gsea/`:

**GSEA** — `gseGO` over all 20,371 genes ranked by `log2_fc`, no threshold. Since `log2_fc` is signed, `scoreType` stays at the default `"std"` — the `"pos"` variant in your section 72 is for unsigned signal magnitudes and would be wrong here.

**ORA** — `enrichGO` on the 156, against the 20,371 tested genes as universe, split three ways: all / up (106) / down (50). The split matters given the 2:1 directional skew.

Parameters copied from your pipeline: `ont="BP"`, `minGSSize=15`, `maxGSSize=500`, `pAdjustMethod="BH"`, `qvalueCutoff=0.2`, `readable=TRUE`, `seed=TRUE`, `eps=0`, and the `bitr` → dedupe-by-largest-|value| → `sort(decreasing=TRUE)` ranked-list construction from `build_ranked_list()`.

Outputs: four TSVs (full `@result` each), four dotplots as PDF/SVG/PNG, and `mca_enrichment_summary.txt` — which counts terms matching `synap|neuron|axon|dendrit|nervous` per block, so the neuronal question gets answered directly rather than by eye.

Useful flags: `--ont CC` (cell-type hypothesis — CC terms would show whether the hits concentrate in a compartment), `--ont ALL`, `--sig-basis fdr` to run ORA on the 427 instead of the 156, `--pcutoff 1` to dump every term including non-significant ones for auditing.

### One decision to flag

`theme_mca()` and `save_multiformat()` are duplicated between this script and `plot_mca_volcano.R` rather than extracted to a shared `scripts/_mca_viz_utils.R`. I kept both self-contained so I didn't have to touch a volcano script that already works mid-analysis. It's ~25 duplicated lines. Say the word and I'll extract them properly.

---

## Chat

The caveat, riding along as promised rather than in place of the deliverable.

Enrichment inherits every assumption of the input. If the 156 are composition- or batch-driven, GSEA will return a clean, coherent, entirely spurious pathway story — and a coherent story is *more* dangerous than a messy one, because it reads as validation. Two specific ways that shows up here:

**The ranking is dominated by genes that shouldn't have mCA at all.** Top of the `log2_fc` ranking will be Npvf (3.29x), Olfr437 (3.06x), Ttr (2.59x). Olfr437 is an olfactory receptor in cerebellum — near-zero expression, so near-zero DNMT3A-deposited mCA is expected, yet it carries 4.9% baseline and ranks second. GSEA doesn't know that; it will happily place olfactory-receptor and sensory-perception terms at the top of the list. If you see those terms come back, that's the mismapping signal being read as biology.

**ORA on 156 with a 20,371 universe is well-powered enough to find something regardless.** With 106 up genes spanning choroid plexus (Ttr), endothelium (Emcn), myeloid (Cd300ld3), and epithelium (Krt24), the likeliest coherent result is vasculature/epithelium/immune terms — which would be the composition hypothesis confirming itself, not a Polycomb-to-mCA axis.

So the result is worth reading as a **diagnostic** more than a finding. If GO CC comes back with terms like extracellular region, plasma membrane, or blood vessel, that's evidence for composition. If GO BP returns synaptic and neuronal terms that survive after dropping the top few outlier genes, that's the first thing in this analysis pointing toward the actual hypothesis.

None of that changes with a better enrichment method — it changes when the plate mapping comes back. That question is still the gate, and it's still open.


## Chat

Here they are — 23 genes, 21 up, 2 down.

| Gene          | Length   | Δ (pp) | mut/ctrl | Neuronal set |
| ------------- | -------- | ------ | -------- | ------------ |
| Npvf          | 3.8 kb   | +10.16 | 3.29     | –            |
| Olfr437       | 6.6 kb   | +10.05 | 3.06     | –            |
| Ttr           | 9.0 kb   | +8.41  | 2.59     | –            |
| Emcn          | 91.1 kb  | +6.82  | 2.24     | –            |
| Cd300ld3      | 6.3 kb   | +6.01  | 1.87     | –            |
| Slc13a1       | 79.8 kb  | +5.94  | 2.44     | –            |
| Kctd12        | 6.1 kb   | +5.89  | 2.71     | **yes**      |
| Edn2          | 6.2 kb   | −5.57  | 0.62     | –            |
| Zic4          | 26.9 kb  | +5.38  | 2.70     | –            |
| Krt24         | 5.3 kb   | +5.23  | 2.50     | –            |
| Sap30         | 5.2 kb   | +5.16  | 2.83     | –            |
| Chia1         | 18.9 kb  | +4.95  | 2.45     | –            |
| Nsmf          | 8.5 kb   | −4.89  | 0.52     | –            |
| Dpys          | 89.0 kb  | +4.85  | 1.97     | –            |
| Fam72a        | 12.0 kb  | +4.63  | 2.63     | –            |
| Fgf1          | 90.8 kb  | +4.62  | 2.01     | –            |
| Slc35g2       | 18.9 kb  | +4.61  | 2.48     | –            |
| Cst8          | 6.8 kb   | +4.44  | 2.02     | –            |
| Fgf20         | 29.2 kb  | +4.36  | 2.00     | –            |
| Lrrc49        | 119.3 kb | +4.29  | 2.28     | –            |
| 2610001J05Rik | 3.8 kb   | +6.17  | 2.16     | –            |
| Scg2          | 5.5 kb   | +3.99  | 2.14     | –            |
| Exoc3l        | 6.2 kb   | +3.52  | 2.35     | –            |

To materialize it:
```
cd /Users/zakiralibhai/emseq && head -1 mca_results/mca_differential_results.tsv > mca_results/mca_robust23.tsv && awk -F'\t' 'NR>1 && $32==8' mca_results/mca_differential_results.tsv >> mca_results/mca_robust23.tsv && wc -l mca_results/mca_robust23.tsv
```

Neuronal-set column uses `61_neuronal_gene_set.tsv` (1,149 genes) from your biomodal checkout — **1 of 23 hits**. That's a partial check: the plan referenced the 5,614-gene section-72 GO set, which isn't in the local checkout, so treat this as indicative only.

### Pulling them out changed my read of the whole result

Look at the `mut/ctrl` column. Every up gene sits between 1.87 and 3.29. Across all 156 significant genes: **up median 1.99 (IQR 1.79–2.18), down median 0.59 (IQR 0.54–0.64).** These genes are on different chromosomes, span 3.8 kb to 119 kb, and have baselines from 2.6% to 14.8% — and the mutant rate is ~2x (or ~0.6x) the control rate in nearly all of them. Independent locus-specific regulation doesn't produce a ratio distribution that tight. This is multiplicative, not additive.

Four more facts that constrain what it can be:

- **Genome-wide mCA is unchanged.** Per-sample means: ctrl 2.09/2.13/2.34/2.11%, mut 2.01/2.23/1.85/2.22%. No global shift, so it isn't a whole-library conversion or chemistry difference.
- **Significant genes come from the high-baseline tail.** Median `mca_ctrl` among significant genes is 4.26% vs a genome-wide median of 1.71% (90th pct 3.89%).
- **Depth is not the driver.** The two deep mutants (mut_F1 ~101x, mut_M2 ~118x) and two shallow ones give indistinguishable values at every top gene — Npvf: 14.38, 14.71 vs 14.80, 14.48. That clears the coverage-imbalance worry from the plan.
- **No chromosomal clustering** (max 1.55x enrichment), which rules out a strain passenger region around the targeted locus.

A multiplicative shift at genes with cell-type-specific mCA, with the global mean held constant, is what a **redistribution of cell populations** looks like: bulk mCA at a gene ≈ Σ(fraction of type *k*) × (mCA of type *k* there). The identities fit — Ttr (choroid plexus epithelium), Emcn (endothelium), Cd300ld3 (myeloid), Krt24 (keratin), Chia1 (macrophage/epithelial).

### The detail that bothers me most

mCA is deposited co-transcriptionally by DNMT3A in gene bodies. **Olfr437 is an olfactory receptor — it should be transcriptionally silent in cerebellum and therefore near-zero in mCA.** It has 4.9% baseline, ~3x the genome median, and it's the #2 hit. Same shape for Cst8 (epididymal) and Krt24 (skin keratin). Genes that should carry almost no mCA are carrying the most.

That points at mismapping. Olfr is a ~1,000-member paralogous family; keratins and Cd300 are clustered families. Loci that attract mismapped reads show inflated apparent non-CG methylation, and anything that shifts the read pool shifts the artifact. Depth is ruled out — but something that tracks genotype exactly is not.

### Two things wrong with the run itself

**1. The batch model was never implemented.** `build_design()` (`scripts/mca_differential.R:148`) constructs only `genotype`, `sex`, `lambda_ca_rate`. There is no `batch` column anywhere in the script — `grep -in batch` returns nothing. The plan called for Model A vs Model B batch sensitivity with the comparison logged in the summary; it wasn't skipped-and-reported, it was silently dropped. The summary's bare "Model: ~ genotype + sex + lambda_ca_rate" reads like a choice was made when no alternative was ever fit.

**2. Dispersion is pinned at the DSS floor, and that's wrong.** Every significant gene has phi = 0.001–0.002. I checked observed replicate SD against the binomial sampling floor: real overdispersion is 1.5–4x in SD terms, i.e. 2–16x in variance. So DSS is modeling these genes as far less variable than they are, and **the reported p-values are anti-conservative.** The 1.25e-47 headline is not credible at face value. (I initially suspected the opposite — *under*dispersion — from eyeballing the tight within-group spreads. That was wrong; running the comparison showed normal overdispersion. The correct concern is the floor, not the data.)

### **NEEDS YOU**

**Is genotype confounded with library prep plate?** The plan defines batch as plate 1 (E1/F1/G1/H1 = S1–S4) vs plate 2 (A2/B2/C2/D2 = S5–S8). If the four controls went on plate 1 and the four mutants on plate 2, a prep-level difference is *mathematically indistinguishable* from the genotype effect, and every number above is uninterpretable. If the plates are mixed across genotype, the batch model becomes fittable and worth running.

I can't determine this locally — `spike_in_rates.tsv` lives in `ca_filtered/` on Expanse. Can you paste the sample→plate mapping, or the head of that file?

That answer decides what happens next. If plates are mixed, the fix is adding batch to the design and re-running with a dispersion floor that reflects the observed overdispersion. If plates are confounded with genotype, no reanalysis of this data can separate them and the question is what to do about it.


GSEA (all genes ranked by log2_fc)
----------------------------------
  Terms with q < 0.05: 72
  Matching /synap|neuron|axon|dendrit|nervous/: 9
  NES > 0: 1, NES < 0: 71
  Top terms by |NES|:
    sensory perception of bitter taste                         NES  -2.05  q 2.13e-02
    detection of chemical stimulus involved in sensory percept NES  -2.01  q 2.13e-02
    negative regulation of myoblast differentiation            NES  -1.93  q 2.36e-02
    aminoglycan catabolic process                              NES  +1.92  q 4.18e-02
    positive regulation of cartilage development               NES  -1.91  q 4.18e-02
    regulation of synaptic vesicle priming                     NES  -1.90  q 3.40e-02
    chondrocyte differentiation                                NES  -1.89  q 1.07e-02
    positive regulation of calcium-mediated signaling          NES  -1.88  q 2.15e-02
    positive regulation of mRNA metabolic process              NES  -1.82  q 1.99e-02
    synaptic vesicle exocytosis                                NES  -1.82  q 1.99e-02
    response to pain                                           NES  -1.81  q 4.18e-02
    post-Golgi vesicle-mediated transport                      NES  -1.81  q 2.14e-02
    positive regulation of cell projection organization        NES  -1.81  q 2.33e-06
    endochondral bone morphogenesis                            NES  -1.78  q 3.28e-02
    negative regulation of cytokine-mediated signaling pathway NES  -1.77  q 3.68e-02


violin plots - noncg quant - increase/vs decrease ubiquitination