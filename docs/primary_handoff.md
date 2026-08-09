## Non-CG methylation analysis in BAP1-KO cerebellum — context and goals

### The system

BAP1 is a deubiquitinase that removes H2AK119ub (a Polycomb-deposited histone mark). In BAP1-knockout mouse cerebellum, H2AK119ub accumulates genome-wide. Prior work in this lab (Hi-C) identified ~2,910 dysregulated chromatin loops, and dual-modality DUET evoC sequencing revealed a coordinated pattern of 5mC increase / 5hmC decrease at ~6,750 gene-body loci (85% concordance among co-significant genes). Gene bodies are the primary affected compartment (42% of genes significant for 5mC, 47% for 5hmC), while promoters show almost no change (<0.2%).

### Why non-CG methylation matters here

Neurons are the exception to the rule that mammalian DNA methylation lives at CpG. During early postnatal life, DNMT3A deposits non-CG methylation (mostly mCA) across transcribed gene bodies, and in mature neurons it becomes ~25% of all methylation. MeCP2 is the reader: it binds mCA and dampens transcription, with the effect scaling with gene length (longer genes accumulate more mCA and are more MeCP2-sensitive). This is settled biology (Guo 2014, Chen 2015, Gabel 2015, Lavery 2020, Moore 2025).

The connection to BAP1 is the open question. No published paper links H2AK119ub to non-CG methylation or to MeCP2 chromatin-level function. The mechanistic pieces exist separately — DNMT3A recognizes H2AK119ub via its ubiquitin-recognition domain (Gretarsson et al. 2024), DNMT3A writes mCA in neurons (Stroud 2017), and MeCP2 reads mCA (Guo 2014) — but no one has connected "H2AK119ub → DNMT3A → mCA → MeCP2" as a single axis. This lab's data (H2AK119ub CUT&RUN + dual 5mC/5hmC + expression + loops in the same genetic model) can test whether BAP1 loss reshapes the mCA landscape that MeCP2 reads.

### The analytical approach: gene-body aggregation

Site-level differential mCA analysis is not feasible with whole-genome sequencing. At ~1% methylation with a ~0.5% conversion noise floor, detecting a difference at a single CA site would require thousands of reads covering that site — genome-wide, that's trillions of reads. Nobody does this.

What works instead is aggregating all CA sites within a gene body. A 50kb gene contains roughly 3,050 CA dinucleotides (CA occurs at ~6.1% of positions in the mouse genome). At 25x genome coverage, those 3,050 sites × 25 reads each give ~76,000 total observations per gene per sample, which is more than enough to estimate a gene-level mCA rate with high precision. This trades spatial resolution (individual sites) for statistical power (gene-level rates from thousands of pooled sites).

This approach is biologically appropriate because mCA's functional unit is the gene body, not the individual site. MeCP2 binding and transcriptional dampening scale with the density of mCA across the gene body as a whole (Gabel 2015). It's also the approach the field uses — Luo et al. 2017 built their neuronal mCA atlas this way.

Neuronal genes are famously long (Syt1 ~340kb, Nrxn1 1Mb+), which means the genes most biologically relevant to mCA/MeCP2 function are also the best-powered for detection. This works in the project's favor.

### Detection thresholds

The minimum detectable difference (MDD) for a ctrl-vs-mut comparison depends on gene length, sequencing depth, and the unknown biological variance between replicates. Two variance components combine: sampling noise (from finite read counts) and biological noise (replicate-to-replicate variability in mCA rates). Calculations below assume n = 4 per group, 80% power, Bonferroni correction over ~5,000 testable genes (z ≈ 5.0), and baseline mCA ≈ 1%.

At **25x genome coverage (~400M reads/sample)** with low biological variance (σ_bio = 0.03%):

|Gene length|CA sites|MDD|
|---|---|---|
|10kb|610|0.30%|
|20kb|1,220|0.22%|
|50kb|3,050|0.17%|
|100kb|6,100|0.14%|
|340kb|20,700|0.12%|

At **37.5x (~600M reads/sample)**, long-gene MDD improves modestly (50kb: 0.15%), because biological variance is already the dominant noise source. At **50x (~800M reads/sample)** and beyond, gains are marginal for 50kb+ genes. The biggest beneficiaries of deeper sequencing are shorter genes (10–20kb), where sampling noise is still large relative to biological noise.

The biological variance estimate (σ_bio) is the most consequential unknown. If σ_bio is 0.05% instead of 0.03%, MDDs increase by 30–40% across the board. The first batch of replicates will reveal which regime applies.

Published neuronal mCA differences between cell types run 1–5%, and within-genotype perturbations show 0.2–0.5% shifts (Luo et al. 2017). Detection of 0.15–0.20% differences at 50kb+ genes is therefore biologically meaningful — it sits at or below the scale of known perturbation effects.

### The conversion noise floor

Lambda spike-in controls from prior EM-seq runs showed 0.5–0.9% apparent non-CG "methylation" on unmethylated DNA — the conversion noise floor. Biological non-CG signal is ~1%, so the noise floor is roughly half the signal. For gene-body aggregation this is a systematic per-sample bias, not random noise: it can be measured from the spike-in and either subtracted as a correction or included as a covariate/offset in the statistical model. Every library must include the spike-in controls for this to work.

### The dataset

Eight EM-seq samples from BAP1-KO and wildtype adult mouse cerebellum (n = 4 per genotype: 2 male, 2 female per group). These are a new sequencing batch, separate from and deeper than the initial DUET cohort. Sex and batch are available as covariates.

### What the non-CG analysis would test

Whether BAP1 loss — through H2AK119ub accumulation and its downstream effects on DNMT3A recruitment and/or TET access — alters gene-body mCA rates at neuronal and synaptic genes. If it does, this would be the first evidence connecting the Polycomb/H2AK119ub axis to the mCA/MeCP2 gene-regulation system, and would link the lab's existing 5mC/5hmC findings and Hi-C loop remodeling to a third, functionally distinct epigenetic layer.
