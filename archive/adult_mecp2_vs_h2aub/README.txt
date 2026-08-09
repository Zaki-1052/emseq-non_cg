MeCP2 vs H2AK119Ub GENE-LEVEL QUADRANT ANALYSIS (v2)

Model: BAP1-KO (Math1-Cre conditional knockout, adult cerebellum)
Axes: MeCP2 and H2AK119Ub log2(FC, BAP1-KO / ctrl)

Script: adult_MeCP2vsH2AUb_v2.R
Original scripts kept for reference: adult_MeCP2vsH2AUb.R, log2log2plot_physicaloverlap.R

GENE-LEVEL QUADRANT PLOT

1. MeCP2 and H2AK119Ub DiffBind output filtered for Conc >= 4.
2. Peaks annotated to mm10 genes via ChIPseeker.
   Only gene-body peaks retained (Promoter, Exon, Intron, 5' UTR, 3' UTR).
   Distal intergenic, downstream (>300kb), and other non-gene-body peaks removed.
3. Multiple peaks per gene collapsed to one FC value per mark per gene:
    * Gene FC = median(FC) of significant peaks only (FDR < 0.05)
    * No FDR-based weighting. FDR is used as a binary filter (pass/fail at 0.05),
      not as a continuous weight. P-values and FDR are threshold statistics,
      not measures of effect importance.
      Ref: Wasserstein & Lazar 2016, ASA Statement on P-Values
      (doi:10.1080/00031305.2016.1154108)
4. Significance (pink label): gene has at least one sig peak (FDR < 0.05)
   in BOTH MeCP2 and H2AK119Ub.

PHYSICAL OVERLAP ANALYSIS (integrated from separate script)

1. MeCP2 and H2AK119Ub peaks tested for physical overlap via GenomicRanges.
2. Reciprocal overlap filter: >= 50% of each peak must be covered by the other.
3. Significance: both overlapping peaks must have FDR < 0.05.
4. Cross-referenced to gene-level table (has_physical_overlap column).

GO ANALYSIS

1. GO-eligible genes: significant in both marks AND |log2FC| >= 0.5 in both marks.
   The FC threshold removes genes near the origin with trivial effect sizes.
2. GO enrichment run on all four quadrants (Q1-Q4), each with >= 10 eligible genes.
3. Two methods run in parallel for comparison:
    a. enrichGO (clusterProfiler) with CORRECTED background:
       universe = all genes with peaks in both marks, not the full genome.
    b. goseq with gene-length bias correction:
       models the relationship between gene length and selection probability.
       If enrichGO and goseq agree, the enrichment is robust to length bias.
       If they diverge, long-gene confounding may inflate the enrichGO result.
4. BP, MF, and CC ontologies run for enrichGO; goseq runs BP.
5. BH multiple testing correction on all results.

DIAGNOSTICS

1. Peak count vs. significance: violin/boxplot of total peaks per gene by sig status.
   Tests whether genes with more peaks are more likely to be called significant
   (expected if significance is defined as "any peak passes FDR < 0.05").
2. Gene length vs. significance: violin/boxplot of gene length by sig status.
   TODO: Check this diagnostic after running. If significant genes are substantially
   longer, the goseq results should be preferred over enrichGO, and the gene-length
   confound should be discussed in any manuscript text interpreting the GO results.
   Neuronal genes (cadherins, contactins, SLITRKs, NRXNs) are among the longest
   protein-coding genes in the mammalian genome, so a length bias could explain
   part of the synapse/cognition GO enrichment independent of MeCP2/H2AK119Ub biology.

OUTPUT FILES

CSVs:
  gene_level_results_v2.csv         Full gene-level table
  quadrant_counts_v2.csv            Quadrant statistics
  top_labeled_genes_v2.csv          Genes labeled on quadrant plot
  overlap_loci_results.csv          Physical overlap analysis
  GO_enrichGO_{Q}__{ont}.csv        enrichGO results per quadrant per ontology
  GO_goseq_{Q}.csv                  goseq results per quadrant

Plots (each saved as PDF/SVG/PNG/JPEG via multi_format_output.R):
  quadrant_plot_v2/                 Gene-level scatter
  overlap_plot_v2/                  Physical overlap scatter
  diagnostic_peak_count_vs_sig/     Peak count bias diagnostic
  diagnostic_gene_length_vs_sig/    Gene length bias diagnostic
  GO_enrichGO_{Q}_{ont}_dotplot/    enrichGO dotplots
