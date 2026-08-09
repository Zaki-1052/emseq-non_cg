mCH + MeCP2 + H2AK119Ub THREE-WAY INTEGRATION (v2)

Model: BAP1-KO (Math1-Cre conditional knockout, adult cerebellum)
Axes: MeCP2, H2AK119Ub log2(FC, BAP1-KO / ctrl), mCH log2FC (edgeR)

Script: mecp2_h2aub_mch_v2.R
Original scripts kept for reference:
  mecp2_h2aub_all_elements_txdb.R
  runfirst_mecp2_h2aub_genebodies_gencode.R
  runsecond_mecp2_h2aub_noncg_gencode.R

UPSTREAM DEPENDENCY

This script reads gene_level_results_v2.csv from the corrected
MeCP2 vs H2AK119Ub analysis (adult_mecp2_vs_h2aub/adult_MeCP2vsH2AUb_v2.R).
That script must be run first. The upstream v2 analysis uses:
  - ChIPseeker/TxDb annotation with gene-body filter
  - Median FC of significant peaks (no FDR weighting)
  - Corrected GO background
  - Physical overlap analysis

mCH DATA

Non-CG methylation edgeR results from results/03_differential/mch_differential_results.tsv.
21,097 rows, 23 duplicate gene names (same gene_name, different ENSMUSG IDs from
overlapping gene models). Deduplicated by keeping the entry with the highest
absolute edger_logFC per gene_name.

GENCODE COORDINATES

Gene body coordinates from gencode.vM25.mouse.genes.annotation.bed (mm10, vM25).
Used for positional cross-reference. Duplicate gene names (21 entries) resolved
by keeping the longest gene body.

THREE-WAY MERGE

Inner join of v2 gene-level results (SYMBOL) with deduplicated mCH output (gene_name).
Left join of gencode coordinates for gene body positions.
Only genes present in all three sources are retained.

SIGNIFICANCE DEFINITIONS

For mCH scatter plots:
  mCH vs MeCP2: significant if edger_fdr < 0.05 AND mecp2_has_sig (any MeCP2 peak FDR < 0.05)
  mCH vs H2AK119Ub: significant if edger_fdr < 0.05 AND h2a_has_sig

For three-way scatter:
  mCH direction coloring: grey = mCH not sig, blue = mCH sig and down, red = mCH sig and up

mCH significance uses FDR only (no FC threshold). The genome-wide mCH shift in BAP1-KO
is expected biology, so FDR < 0.05 captures the relevant genes without additional filtering.

PLOTS

1. mCH vs MeCP2 scatter: mCH logFC on y-axis, MeCP2 FC on x-axis.
   Top 10 significant genes labeled from Q1 and Q3.

2. mCH vs H2AK119Ub scatter: mCH logFC on y-axis, H2AK119Ub FC on x-axis.
   Top 10 significant genes labeled from Q1 and Q3.

3. Three-way scatter: MeCP2 FC on x-axis, H2AK119Ub FC on y-axis (same axes as v2
   quadrant plot), colored by mCH significance and direction. Shows where the mCH
   signal concentrates across the MeCP2/H2AK119Ub quadrants.

THREE-WAY GO ENRICHMENT

For each MeCP2/H2AK119Ub quadrant with >= 10 three-way significant genes:
  Gene set: genes significant in BOTH ChIP marks (from v2 sig column) AND
            mCH-significant (edger_fdr < 0.05).

Two methods run for comparison:
  a. enrichGO (clusterProfiler) with corrected background:
     universe = all genes in merged table, not the full genome.
  b. goseq with gene-length bias correction:
     uses mCH pipeline gene body length as bias data.
     If enrichGO and goseq agree, the enrichment is robust to length bias.
     If they diverge, long-gene confounding may inflate the enrichGO result.

BP, MF, CC ontologies for enrichGO; goseq runs BP.
BH multiple testing correction on all results.

DIAGNOSTICS

1. Gene length vs mCH significance: violin/boxplot of gene body length by mCH
   significance status. Checks whether longer genes are more likely to show
   significant mCH changes (expected, given that longer genes have more CH sites
   and therefore more statistical power).

2. Coverage vs mCH significance: violin/boxplot of total coverage by mCH
   significance status. Checks whether coverage drives significance.

Both plots annotated with n per group and median values.

OUTPUT FILES

CSVs:
  mch_mecp2_h2aub_merged_v2.csv         Full merged table
  quadrant_counts_mecp2_v2.csv           mCH vs MeCP2 quadrant stats
  quadrant_counts_h2a_v2.csv             mCH vs H2AK119Ub quadrant stats
  three_way_quadrant_summary.csv         MeCP2/H2AK119Ub quadrants with mCH breakdown
  GO_three_way_{Q}__{ont}.csv            enrichGO results per quadrant per ontology
  GO_goseq_three_way_{Q}.csv             goseq results per quadrant

Plots (each saved as PDF/SVG/PNG/JPEG via multi_format_output.R):
  mch_vs_mecp2_v2/                       mCH vs MeCP2 scatter
  mch_vs_h2aub_v2/                       mCH vs H2AK119Ub scatter
  three_way_scatter_v2/                   MeCP2 vs H2AK119Ub colored by mCH status
  diagnostic_mch_genelength/             Gene length vs mCH significance
  diagnostic_mch_coverage/               Coverage vs mCH significance
  GO_three_way_{Q}_{ont}_dotplot/        GO dotplots
