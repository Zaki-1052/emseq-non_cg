# TODO

## Port section 20: RNA-seq expression integration

Biomodal `section_20_coordinated_rnaseq.R` was missed during the port triage.
It was not classified as Not Portable, Portable, or Adapts — it fell through
the triage entirely, likely because its title references "coordinated mC/hmC"
which made it look dual-modality. The analytical core is modality-independent.

### What the biomodal version does

1. Loads RNA-seq DESeq2 results (`adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx`)
2. Joins coordinated mC-up/hmC-down genes with RNA-seq log2FC and padj
3. Classifies expression outcomes (Up / Down / Unchanged) at two thresholds
4. Figures:
   - 20a: Stacked bar — expression direction breakdown (coordinated vs other DMR vs all)
   - 20b: Scatter — combined methylation effect vs log2FC (Spearman correlation)
   - 20c: Violin — log2FC distributions (coordinated vs other DMR genes)
   - 20d: Heatmap — 2x2 enrichment (methylation direction x expression direction)
5. Fisher's exact tests comparing expression direction ratios between gene sets

### What the mCH port needs

- Gene selection: mCH-significant genes (FDR < 0.05) stratified by mCH direction
  (hypo in mutant = 4,915 genes, hyper = 642) instead of coordinated mC/hmC
- Same RNA-seq data file — same BAP1-KO model, same animals
- Same figures adapted: stacked bar (mCH-hypo vs mCH-hyper vs all), scatter
  (mCH logFC vs RNA-seq log2FC), violin, 2x2 enrichment
- Drop the "combined_effect = |mC| + |hmC|" metric — use edger_logFC directly

### Data files needed

- Source: `mariner_hi-c/tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx` (6.1 MB)
- Copy to: `data/rnaseq/adult_timepoint_rnaseq_deseq2.xlsx` (or similar)
- BigWigs: `RNActrl.bw` and `RNAmut.bw` already synced to local bigwigs dir

### Where it goes

- Script: `scripts/sections/80_cross_modality/80_03_rnaseq_integration.R`
  (or new group if 80 is full)
- SLURM wrapper: `scripts/sections/slurm/80_03_rnaseq_integration.sb`
- Results: `results/sections/80_cross_modality/`

### Why this matters

The abstract says "transcriptional downregulation" and "dysregulating neuronal
and synaptic gene expression." The aims discuss expression outcomes of BAP1 loss.
This section is the only place in the pipeline that tests whether methylation
changes correspond to expression changes. Without it, the paper's claims about
transcription are not supported by the EM-seq analysis.

## Audit other missed biomodal sections

The port plan triaged 35 of 93+ biomodal sections. Section 20 was missed.
A focused re-triage of remaining un-categorized sections should check for
other analyses that are modality-independent but were skipped because their
code references mC/hmC variables. Candidates to check:

- Sections 51, 53 (MeCP2 non-CG methylation) — ask mCH questions the EM-seq
  data can answer better than biomodal evoC, but need substantial adaptation
- Sections 57, 58 (Ecker non-CG validation, dose-response) — use external WGBS
  data to validate non-CG findings, may be relevant
- Section 18 (K119ub BigWig signal) — gene-level K119ub signal quantification
- Section 42 (max significance gene list) — could produce a useful gene table

Priority: section 20 (RNA-seq) first, audit second.
