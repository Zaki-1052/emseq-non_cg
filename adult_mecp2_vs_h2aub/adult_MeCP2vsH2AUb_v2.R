# adult_mecp2_vs_h2aub/adult_MeCP2vsH2AUb_v2.R
#
# MeCP2 vs H2AK119Ub gene-level quadrant analysis (v2)
#
# Rewrites the original analysis with fixes from code review:
#   - Median FC of significant peaks (no FDR weighting)
#   - Gene-body annotation filter
#   - Corrected GO background + goseq gene-length correction
#   - GO enrichment on all four quadrants
#   - Physical overlap analysis integrated
#   - Peak-count and gene-length bias diagnostics
#
# Input: DiffBind peak tables (mecp2.txt, h2aub.txt) with Conc >= 4
# Output: gene-level CSVs, quadrant plots, overlap plots, GO results, diagnostics

############################################################
# 1. Setup
############################################################

library(here)
library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(goseq)

source(here("scripts", "utils", "multi_format_output.R"))

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

out_dir <- here("adult_mecp2_vs_h2aub")

mecp2_raw <- read.table(
  here("adult_mecp2_vs_h2aub", "mecp2.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

h2aub_raw <- read.table(
  here("adult_mecp2_vs_h2aub", "h2aub.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

cat("MeCP2 peaks loaded:", nrow(mecp2_raw), "\n")
cat("H2AK119Ub peaks loaded:", nrow(h2aub_raw), "\n")

############################################################
# 2. GRanges conversion and peak annotation
############################################################

gr_mecp2 <- GRanges(
  seqnames = mecp2_raw$seqnames,
  ranges = IRanges(start = mecp2_raw$start, end = mecp2_raw$end)
)
mcols(gr_mecp2) <- mecp2_raw[, c("Fold", "p.value", "FDR")]

gr_h2a <- GRanges(
  seqnames = h2aub_raw$seqnames,
  ranges = IRanges(start = h2aub_raw$start, end = h2aub_raw$end)
)
mcols(gr_h2a) <- h2aub_raw[, c("Fold", "p.value", "FDR")]

peakAnno_mecp2 <- annotatePeak(gr_mecp2, TxDb = txdb, annoDb = "org.Mm.eg.db")
anno_mecp2 <- as.data.frame(peakAnno_mecp2)

peakAnno_h2a <- annotatePeak(gr_h2a, TxDb = txdb, annoDb = "org.Mm.eg.db")
anno_h2a <- as.data.frame(peakAnno_h2a)

mecp2_df <- data.frame(
  SYMBOL       = anno_mecp2$SYMBOL,
  geneID       = anno_mecp2$geneId,
  annotation   = anno_mecp2$annotation,
  distanceToTSS = anno_mecp2$distanceToTSS,
  peak_width   = mecp2_raw$width,
  mecp2_Fold   = mecp2_raw$Fold,
  mecp2_FDR    = mecp2_raw$FDR,
  stringsAsFactors = FALSE
)

h2a_df <- data.frame(
  SYMBOL       = anno_h2a$SYMBOL,
  geneID       = anno_h2a$geneId,
  annotation   = anno_h2a$annotation,
  distanceToTSS = anno_h2a$distanceToTSS,
  peak_width   = h2aub_raw$width,
  h2a_Fold     = h2aub_raw$Fold,
  h2a_FDR      = h2aub_raw$FDR,
  stringsAsFactors = FALSE
)

############################################################
# 3. Gene-body filtering
############################################################

gene_body_pattern <- "^(Promoter|5' UTR|3' UTR|Exon|Intron)"

mecp2_df <- mecp2_df %>%
  filter(
    !is.na(SYMBOL),
    grepl(gene_body_pattern, annotation)
  )

h2a_df <- h2a_df %>%
  filter(
    !is.na(SYMBOL),
    grepl(gene_body_pattern, annotation)
  )

cat("MeCP2 peaks after gene-body filter:", nrow(mecp2_df), "\n")
cat("H2AK119Ub peaks after gene-body filter:", nrow(h2a_df), "\n")

############################################################
# 4. Gene-level collapse
#
# Aggregation: median FC of significant peaks per gene.
# FDR is a binary gatekeeper (in/out at 0.05), not a weight.
# See: Wasserstein & Lazar 2016, ASA Statement on P-Values
# (doi:10.1080/00031305.2016.1154108)
############################################################

fdr_cutoff <- 0.05

mecp2_gene <- mecp2_df %>%
  group_by(SYMBOL) %>%
  summarise(
    n_mecp2_peaks     = n(),
    n_mecp2_sig_peaks = sum(mecp2_FDR < fdr_cutoff),
    mecp2_has_sig     = any(mecp2_FDR < fdr_cutoff),
    mecp2_Fold        = if (sum(mecp2_FDR < fdr_cutoff) > 0)
                            median(mecp2_Fold[mecp2_FDR < fdr_cutoff])
                          else median(mecp2_Fold),
    mecp2_min_FDR     = min(mecp2_FDR),
    .groups = "drop"
  )

h2a_gene <- h2a_df %>%
  group_by(SYMBOL) %>%
  summarise(
    n_h2a_peaks     = n(),
    n_h2a_sig_peaks = sum(h2a_FDR < fdr_cutoff),
    h2a_has_sig     = any(h2a_FDR < fdr_cutoff),
    h2a_Fold        = if (sum(h2a_FDR < fdr_cutoff) > 0)
                          median(h2a_Fold[h2a_FDR < fdr_cutoff])
                        else median(h2a_Fold),
    h2a_min_FDR     = min(h2a_FDR),
    .groups = "drop"
  )

############################################################
# 5. Merge and significance
############################################################

gene_df <- inner_join(mecp2_gene, h2a_gene, by = "SYMBOL")

cat("Genes with peaks in both marks:", nrow(gene_df), "\n")

gene_df$sig <- ifelse(
  gene_df$mecp2_has_sig & gene_df$h2a_has_sig,
  "Significant",
  "Not significant"
)

gene_df$sig <- factor(
  gene_df$sig,
  levels = c("Not significant", "Significant")
)

total_genes <- nrow(gene_df)
total_sig_genes <- sum(gene_df$sig == "Significant")

cat("Total genes:", total_genes, "\n")
cat("Significant genes (sig peak in both marks):", total_sig_genes, "\n")

############################################################
# 6. Quadrant assignment and correlation
############################################################

gene_df$quadrant <- case_when(
  gene_df$mecp2_Fold >= 0 & gene_df$h2a_Fold >= 0 ~ "Q1 (MeCP2 up / H2A up)",
  gene_df$mecp2_Fold <  0 & gene_df$h2a_Fold >= 0 ~ "Q2 (MeCP2 down / H2A up)",
  gene_df$mecp2_Fold <  0 & gene_df$h2a_Fold <  0 ~ "Q3 (MeCP2 down / H2A down)",
  gene_df$mecp2_Fold >= 0 & gene_df$h2a_Fold <  0 ~ "Q4 (MeCP2 up / H2A down)"
)

cor_test <- cor.test(
  gene_df$mecp2_Fold,
  gene_df$h2a_Fold,
  method = "spearman"
)

rho_label <- paste0(
  "Spearman ρ = ",
  round(cor_test$estimate, 2),
  " | p = ",
  format.pval(cor_test$p.value, digits = 2, eps = 1e-10)
)

print(cor_test)

quadrant_all <- gene_df %>%
  count(quadrant) %>%
  complete(
    quadrant = c(
      "Q1 (MeCP2 up / H2A up)",
      "Q2 (MeCP2 down / H2A up)",
      "Q3 (MeCP2 down / H2A down)",
      "Q4 (MeCP2 up / H2A down)"
    ),
    fill = list(n = 0)
  ) %>%
  rename(n_total_genes = n) %>%
  mutate(percent_total_genes = round(100 * n_total_genes / total_genes, 1))

quadrant_sig <- gene_df %>%
  filter(sig == "Significant") %>%
  count(quadrant) %>%
  complete(
    quadrant = c(
      "Q1 (MeCP2 up / H2A up)",
      "Q2 (MeCP2 down / H2A up)",
      "Q3 (MeCP2 down / H2A down)",
      "Q4 (MeCP2 up / H2A down)"
    ),
    fill = list(n = 0)
  ) %>%
  rename(n_sig_genes = n) %>%
  mutate(percent_sig_genes = round(100 * n_sig_genes / total_sig_genes, 1))

quadrant_counts <- left_join(quadrant_all, quadrant_sig, by = "quadrant")

quadrant_counts <- quadrant_counts %>%
  mutate(
    label = paste0(
      "Total: n=", n_total_genes, " (", percent_total_genes, "%)\n",
      "Sig: n=", n_sig_genes, " (", percent_sig_genes, "%)"
    )
  )

print(quadrant_counts)

############################################################
# 7. FC threshold for GO eligibility
############################################################

fc_threshold <- 0

gene_df <- gene_df %>%
  mutate(
    go_eligible = sig == "Significant" &
      abs(mecp2_Fold) >= fc_threshold &
      abs(h2a_Fold) >= fc_threshold
  )

n_go_eligible <- sum(gene_df$go_eligible, na.rm = TRUE)
n_sig_below_fc <- total_sig_genes - n_go_eligible

cat("GO-eligible significant genes (|FC| >= ", fc_threshold, " in both marks): ",
    n_go_eligible, "\n", sep = "")
cat("Significant genes below FC threshold: ", n_sig_below_fc, "\n", sep = "")

############################################################
# 8. Physical overlap analysis
############################################################

hits <- findOverlaps(gr_mecp2, gr_h2a)

cat("\nOverlapping peak pairs (any overlap):", length(hits), "\n")

overlap_df <- data.frame(
  mecp2_Fold = mcols(gr_mecp2)$Fold[queryHits(hits)],
  mecp2_FDR  = mcols(gr_mecp2)$FDR[queryHits(hits)],
  h2a_Fold   = mcols(gr_h2a)$Fold[subjectHits(hits)],
  h2a_FDR    = mcols(gr_h2a)$FDR[subjectHits(hits)]
)

overlap_widths <- width(pintersect(gr_mecp2[queryHits(hits)], gr_h2a[subjectHits(hits)]))
mecp2_widths <- width(gr_mecp2[queryHits(hits)])
h2a_widths <- width(gr_h2a[subjectHits(hits)])

overlap_df$overlap_frac_mecp2 <- overlap_widths / mecp2_widths
overlap_df$overlap_frac_h2a <- overlap_widths / h2a_widths

overlap_df <- overlap_df %>%
  filter(overlap_frac_mecp2 >= 0.5 & overlap_frac_h2a >= 0.5)

n_overlap_loci <- nrow(overlap_df)
cat("Overlapping loci after >= 50% reciprocal filter:", n_overlap_loci, "\n")

overlap_df$sig <- ifelse(
  overlap_df$mecp2_FDR < fdr_cutoff & overlap_df$h2a_FDR < fdr_cutoff,
  "Significant",
  "Not significant"
)
overlap_df$sig <- factor(overlap_df$sig, levels = c("Not significant", "Significant"))

overlap_df$quadrant <- case_when(
  overlap_df$mecp2_Fold >= 0 & overlap_df$h2a_Fold >= 0 ~ "Q1 (MeCP2 up / H2A up)",
  overlap_df$mecp2_Fold <  0 & overlap_df$h2a_Fold >= 0 ~ "Q2 (MeCP2 down / H2A up)",
  overlap_df$mecp2_Fold <  0 & overlap_df$h2a_Fold <  0 ~ "Q3 (MeCP2 down / H2A down)",
  overlap_df$mecp2_Fold >= 0 & overlap_df$h2a_Fold <  0 ~ "Q4 (MeCP2 up / H2A down)"
)

overlap_quad_sig <- overlap_df %>%
  filter(sig == "Significant") %>%
  count(quadrant) %>%
  complete(
    quadrant = c(
      "Q1 (MeCP2 up / H2A up)",
      "Q2 (MeCP2 down / H2A up)",
      "Q3 (MeCP2 down / H2A down)",
      "Q4 (MeCP2 up / H2A down)"
    ),
    fill = list(n = 0)
  ) %>%
  mutate(
    percent = round(100 * n / n_overlap_loci, 2),
    label = paste0("n = ", n, "\n", percent, "%")
  )

cat("\nOverlap loci quadrant breakdown (significant):\n")
print(overlap_quad_sig)

# Cross-reference: which genes have physically overlapping peaks?
mecp2_anno_for_overlap <- data.frame(
  mecp2_idx = queryHits(hits),
  SYMBOL = anno_mecp2$SYMBOL[queryHits(hits)],
  stringsAsFactors = FALSE
)

overlap_frac_mecp2_all <- overlap_widths / mecp2_widths
overlap_frac_h2a_all <- overlap_widths / h2a_widths

overlap_pass <- overlap_frac_mecp2_all >= 0.5 & overlap_frac_h2a_all >= 0.5

genes_with_overlap <- unique(mecp2_anno_for_overlap$SYMBOL[overlap_pass])
genes_with_overlap <- genes_with_overlap[!is.na(genes_with_overlap)]

gene_df$has_physical_overlap <- gene_df$SYMBOL %in% genes_with_overlap

cat("Genes with at least one physically overlapping peak pair:",
    sum(gene_df$has_physical_overlap), "of", nrow(gene_df), "\n")

############################################################
# 9. Diagnostic: peak count vs. significance
############################################################

gene_df <- gene_df %>%
  mutate(n_total_peaks = n_mecp2_peaks + n_h2a_peaks)

wilcox_peaks <- wilcox.test(
  n_total_peaks ~ sig,
  data = gene_df
)

cat("\nPeak count vs. significance (Wilcoxon):\n")
print(wilcox_peaks)

peak_summary <- gene_df %>%
  group_by(sig) %>%
  summarise(
    n = n(),
    median_val = median(n_total_peaks),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("n = ", n, "\nmedian = ", median_val))

p_peak_diag <- ggplot(gene_df, aes(x = sig, y = n_total_peaks, fill = sig)) +
  geom_violin(alpha = 0.6, show.legend = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, show.legend = FALSE) +
  geom_text(
    data = peak_summary,
    aes(x = sig, y = max(gene_df$n_total_peaks) * 1.5, label = label),
    size = 4, fontface = "bold", vjust = 1, show.legend = FALSE
  ) +
  scale_fill_manual(values = c("Not significant" = "grey70", "Significant" = "deeppink")) +
  scale_y_log10() +
  labs(
    x = "Significance status",
    y = "Total peaks per gene (MeCP2 + H2AUb, log10)",
    title = "Peak count vs. significance status",
    subtitle = paste0("Wilcoxon p = ", format.pval(wilcox_peaks$p.value, digits = 2))
  ) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

############################################################
# 10. Diagnostic: gene length vs. significance
############################################################

mm10_genes <- genes(txdb)

entrez_to_symbol <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = names(mm10_genes),
  columns = "SYMBOL",
  keytype = "ENTREZID"
)

gene_length_df <- data.frame(
  ENTREZID = names(mm10_genes),
  gene_length = width(mm10_genes),
  stringsAsFactors = FALSE
) %>%
  inner_join(entrez_to_symbol, by = "ENTREZID") %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(SYMBOL) %>%
  summarise(gene_length = max(gene_length), .groups = "drop")

gene_df <- gene_df %>%
  left_join(gene_length_df, by = "SYMBOL")

wilcox_length <- wilcox.test(
  gene_length ~ sig,
  data = gene_df %>% filter(!is.na(gene_length))
)

cat("\nGene length vs. significance (Wilcoxon):\n")
print(wilcox_length)

gene_df_length <- gene_df %>% filter(!is.na(gene_length))

length_summary <- gene_df_length %>%
  group_by(sig) %>%
  summarise(
    n = n(),
    median_kb = round(median(gene_length / 1000), 1),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("n = ", n, "\nmedian = ", median_kb, " kb"))

p_length_diag <- ggplot(gene_df_length, aes(x = sig, y = gene_length / 1000, fill = sig)) +
  geom_violin(alpha = 0.6, show.legend = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, show.legend = FALSE) +
  geom_text(
    data = length_summary,
    aes(x = sig, y = max(gene_df_length$gene_length / 1000) * 1.5, label = label),
    size = 4, fontface = "bold", vjust = 1, show.legend = FALSE
  ) +
  scale_fill_manual(values = c("Not significant" = "grey70", "Significant" = "deeppink")) +
  scale_y_log10() +
  labs(
    x = "Significance status",
    y = "Gene length (kb, log10)",
    title = "Gene length vs. significance status",
    subtitle = paste0("Wilcoxon p = ", format.pval(wilcox_length$p.value, digits = 2))
  ) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

############################################################
# 11. Quadrant scatter plot
############################################################

xmax <- max(abs(gene_df$mecp2_Fold), na.rm = TRUE)
ymax <- max(abs(gene_df$h2a_Fold), na.rm = TRUE)

quadrant_counts$x <- c(0.90 * xmax, -0.30 * xmax, -0.30 * xmax, 0.90 * xmax)
quadrant_counts$y <- c(0.90 * ymax, 0.90 * ymax, -0.90 * ymax, -0.90 * ymax)

top_genes_all <- gene_df %>%
  filter(sig == "Significant") %>%
  mutate(label_strength = abs(mecp2_Fold) + abs(h2a_Fold)) %>%
  arrange(desc(label_strength)) %>%
  slice_head(n = 10)

top_genes_Q3 <- gene_df %>%
  filter(sig == "Significant", quadrant == "Q3 (MeCP2 down / H2A down)") %>%
  mutate(label_strength = abs(mecp2_Fold) + abs(h2a_Fold)) %>%
  arrange(desc(label_strength)) %>%
  slice_head(n = 5)

top_genes_Q2 <- gene_df %>%
  filter(sig == "Significant", quadrant == "Q2 (MeCP2 down / H2A up)") %>%
  mutate(label_strength = abs(mecp2_Fold) + abs(h2a_Fold)) %>%
  arrange(desc(label_strength)) %>%
  slice_head(n = 5)

top_genes_Q4 <- gene_df %>%
  filter(sig == "Significant", quadrant == "Q4 (MeCP2 up / H2A down)") %>%
  mutate(label_strength = abs(mecp2_Fold) + abs(h2a_Fold)) %>%
  arrange(desc(label_strength)) %>%
  slice_head(n = 5)

top_genes <- bind_rows(top_genes_all, top_genes_Q3, top_genes_Q2, top_genes_Q4) %>%
  distinct(SYMBOL, .keep_all = TRUE)

gene_df_plot <- gene_df %>% arrange(sig)

p_quadrant <- ggplot(gene_df_plot, aes(x = mecp2_Fold, y = h2a_Fold, color = sig)) +
  geom_point(size = 1.8, alpha = 0.8) +
  scale_color_manual(values = c(
    "Not significant" = "grey85",
    "Significant" = "deeppink"
  )) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_text(
    data = quadrant_counts,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3,
    fontface = "bold",
    color = "black"
  ) +
  geom_text_repel(
    data = top_genes,
    aes(x = mecp2_Fold, y = h2a_Fold, label = SYMBOL),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    fontface = "bold",
    color = "black"
  ) +
  labs(
    x = "MeCP2 log2(FC, BAP1-KO / ctrl)",
    y = "H2AK119Ub log2(FC, BAP1-KO / ctrl)",
    color = "",
    title = "Gene-level MeCP2 vs H2AK119Ub (v2)",
    subtitle = paste0(
      "n (total) = ", total_genes, " genes | ",
      "n (sig) = ", total_sig_genes, " genes\n",
      rho_label, "\n",
      "Aggregation: median FC of significant peaks | ",
      "GO filter: |FC| ≥ ", fc_threshold
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    plot.subtitle = element_text(size = 10, face = "bold")
  )

print(p_quadrant)

# Overlap scatter plot
overlap_df_plot <- overlap_df %>% arrange(sig)

ovl_xmax <- max(abs(overlap_df$mecp2_Fold), na.rm = TRUE)
ovl_ymax <- max(abs(overlap_df$h2a_Fold), na.rm = TRUE)

ovl_labels <- overlap_quad_sig %>%
  mutate(
    x = c(0.95 * ovl_xmax, -0.95 * ovl_xmax, -0.92 * ovl_xmax, 0.95 * ovl_xmax),
    y = c(0.97 * ovl_ymax, 0.95 * ovl_ymax, -0.95 * ovl_ymax, -0.95 * ovl_ymax)
  )

p_overlap <- ggplot(overlap_df_plot, aes(x = mecp2_Fold, y = h2a_Fold, color = sig)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c(
    "Not significant" = "grey85",
    "Significant" = "deeppink"
  )) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_text(
    data = ovl_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  labs(
    x = "MeCP2 log2(FC, BAP1-KO / ctrl)",
    y = "H2AK119Ub log2(FC, BAP1-KO / ctrl)",
    color = "",
    title = paste0(
      "n = ", n_overlap_loci,
      " physically overlapping loci (≥50% reciprocal overlap)"
    )
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top"
  )

print(p_overlap)

############################################################
# 12. GO enrichment — all four quadrants
############################################################

# Correct background: all genes with peaks in both marks
universe_symbols <- unique(gene_df$SYMBOL)
universe_entrez <- bitr(
  universe_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

cat("\nGO universe size (Entrez IDs):", nrow(universe_entrez), "\n")

quadrant_names <- c(
  "Q1 (MeCP2 up / H2A up)",
  "Q2 (MeCP2 down / H2A up)",
  "Q3 (MeCP2 down / H2A down)",
  "Q4 (MeCP2 up / H2A down)"
)

quadrant_short <- c("Q1", "Q2", "Q3", "Q4")
ontologies <- c("BP", "MF", "CC")

min_go_genes <- 10

go_results <- list()
goseq_results <- list()

# Gene lengths for goseq (all genes in universe)
all_gene_lengths <- gene_length_df %>%
  inner_join(universe_entrez, by = "SYMBOL")

for (i in seq_along(quadrant_names)) {
  q_name <- quadrant_names[i]
  q_short <- quadrant_short[i]

  q_genes <- gene_df %>%
    filter(go_eligible, quadrant == q_name) %>%
    pull(SYMBOL)

  cat("\n", q_short, ": ", length(q_genes), " GO-eligible genes\n", sep = "")

  if (length(q_genes) < min_go_genes) {
    cat("  Skipping GO (< ", min_go_genes, " genes)\n", sep = "")
    next
  }

  q_entrez <- bitr(
    q_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )

  cat("  Converted to Entrez:", nrow(q_entrez), "\n")

  # enrichGO for each ontology
  for (ont in ontologies) {
    key <- paste0(q_short, "_", ont)

    go_res <- enrichGO(
      gene = q_entrez$ENTREZID,
      universe = universe_entrez$ENTREZID,
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      readable = TRUE
    )

    go_results[[key]] <- go_res

    if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
      cat("  enrichGO ", ont, ": ", nrow(as.data.frame(go_res)), " terms\n", sep = "")
    } else {
      cat("  enrichGO ", ont, ": 0 terms\n", sep = "")
    }
  }

  # goseq (gene-length corrected, BP only)
  gene_vector <- as.integer(all_gene_lengths$ENTREZID %in% q_entrez$ENTREZID)
  names(gene_vector) <- all_gene_lengths$ENTREZID
  bias_data <- all_gene_lengths$gene_length
  names(bias_data) <- all_gene_lengths$ENTREZID

  pwf <- nullp(gene_vector, bias.data = bias_data, plot.fit = FALSE)

  goseq_res <- goseq(pwf, "mm10", "knownGene")
  goseq_res$over_represented_padj <- p.adjust(goseq_res$over_represented_pvalue, method = "BH")
  goseq_sig <- goseq_res %>% filter(over_represented_padj < 0.05)

  goseq_results[[q_short]] <- goseq_res

  cat("  goseq: ", nrow(goseq_sig), " significant terms (length-corrected)\n", sep = "")
}

############################################################
# 13. Save outputs
############################################################

cat("\n--- Saving outputs ---\n")

# CSVs
write.csv(gene_df, file.path(out_dir, "gene_level_results_v2.csv"), row.names = FALSE)
write.csv(quadrant_counts, file.path(out_dir, "quadrant_counts_v2.csv"), row.names = FALSE)
write.csv(top_genes, file.path(out_dir, "top_labeled_genes_v2.csv"), row.names = FALSE)
write.csv(overlap_df, file.path(out_dir, "overlap_loci_results.csv"), row.names = FALSE)

# Plots via multi-format utility
save_multiformat_ggplot(p_quadrant, file.path(out_dir, "quadrant_plot_v2"), width = 10, height = 10)
save_multiformat_ggplot(p_overlap, file.path(out_dir, "overlap_plot_v2"), width = 7, height = 7)
save_multiformat_ggplot(p_peak_diag, file.path(out_dir, "diagnostic_peak_count_vs_sig"), width = 7, height = 6)
save_multiformat_ggplot(p_length_diag, file.path(out_dir, "diagnostic_gene_length_vs_sig"), width = 7, height = 6)

# GO enrichGO results and dotplots
for (key in names(go_results)) {
  go_res <- go_results[[key]]
  go_df <- as.data.frame(go_res)

  write.csv(go_df, file.path(out_dir, paste0("GO_enrichGO_", key, ".csv")), row.names = FALSE)

  if (nrow(go_df) > 0) {
    p_go <- dotplot(go_res, showCategory = 20, title = paste0("enrichGO ", key))
    save_multiformat_ggplot(p_go, file.path(out_dir, paste0("GO_enrichGO_", key, "_dotplot")), width = 8, height = 10)
  }
}

# goseq results
for (q_short in names(goseq_results)) {
  goseq_res <- goseq_results[[q_short]]
  write.csv(goseq_res, file.path(out_dir, paste0("GO_goseq_", q_short, ".csv")), row.names = FALSE)
}

cat("\nDone.\n")
