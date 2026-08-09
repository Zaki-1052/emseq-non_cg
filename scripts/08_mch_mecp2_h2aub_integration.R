# scripts/08_mch_mecp2_h2aub_integration.R
#
# Three-way integration: MeCP2 + H2AK119Ub + non-CG methylation (mCH)
#
# Reads corrected gene-level MeCP2/H2AK119Ub results from v2 script,
# merges with mCH edgeR differential output, produces scatter plots,
# Spearman correlations, three-way GO enrichment, and diagnostics.
#
# Upstream dependency: scripts/07_mecp2_h2aub_quadrant.R
#   must be run first to produce gene_level_results_v2.csv
#
# Input:
#   results/07_quadrant/gene_level_results_v2.csv
#   results/03_differential/mch_differential_results.tsv
#   data/gencode.vM25.mouse.genes.annotation.bed

############################################################
# 1. Setup
############################################################

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

out_dir <- here("results", "08_three_way")

############################################################
# 2. Load upstream v2 output
############################################################

gene_df <- read_csv(
  here("results", "07_quadrant", "gene_level_results_v2.csv"),
  show_col_types = FALSE
)

cat("Upstream v2 gene-level results loaded:", nrow(gene_df), "genes\n")
cat("  Significant (both marks):", sum(gene_df$sig == "Significant"), "\n")

############################################################
# 3. Load and deduplicate mCH edgeR output
############################################################

mch_raw <- read_tsv(
  here("results", "03_differential", "mch_differential_results.tsv"),
  show_col_types = FALSE
)

cat("\nmCH edgeR results loaded:", nrow(mch_raw), "genes\n")

dup_genes <- mch_raw %>%
  count(gene_name) %>%
  filter(n > 1)

if (nrow(dup_genes) > 0) {
  cat("  Duplicate gene_names:", nrow(dup_genes), "\n")
  cat("  Genes:", paste(dup_genes$gene_name, collapse = ", "), "\n")
}

mch <- mch_raw %>%
  group_by(gene_name) %>%
  slice_max(abs(edger_logFC), n = 1, with_ties = FALSE) %>%
  ungroup()

cat("  After dedup:", nrow(mch), "genes\n")

############################################################
# 4. Load gencode gene body coordinates
############################################################

gencode_bed <- read.table(
  here("data", "gencode.vM25.mouse.genes.annotation.bed"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

gencode_coords <- gencode_bed %>%
  mutate(
    gencode_chr = paste0("chr", Chromosome),
    gencode_start = Start,
    gencode_end = End
  ) %>%
  select(Name, gencode_chr, gencode_start, gencode_end) %>%
  group_by(Name) %>%
  slice_max(gencode_end - gencode_start, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Gencode gene bodies loaded:", nrow(gencode_coords), "genes\n")

############################################################
# 5. Merge
############################################################

merged <- inner_join(
  gene_df,
  mch %>% rename(
    mch_gene_id = gene_id,
    mch_chr = chr,
    mch_start = start,
    mch_end = end,
    mch_strand = strand,
    mch_gene_type = gene_type,
    mch_gene_length = gene_length
  ),
  by = c("SYMBOL" = "gene_name")
)

merged <- left_join(
  merged,
  gencode_coords,
  by = c("SYMBOL" = "Name")
)

cat("\nMerged genes (all three modalities):", nrow(merged), "\n")
cat("  Of which v2-significant (both ChIP marks):",
    sum(merged$sig == "Significant"), "\n")
cat("  Of which mCH-significant (FDR < 0.05):",
    sum(merged$sig_fdr005, na.rm = TRUE), "\n")

############################################################
# 6. Significance flags
############################################################

merged <- merged %>%
  mutate(
    mch_sig = sig_fdr005 == TRUE,
    mch_sig_vs_mecp2 = factor(
      ifelse(mch_sig & mecp2_has_sig, "Significant", "Not significant"),
      levels = c("Not significant", "Significant")
    ),
    mch_sig_vs_h2a = factor(
      ifelse(mch_sig & h2a_has_sig, "Significant", "Not significant"),
      levels = c("Not significant", "Significant")
    ),
    mch_direction = factor(
      case_when(
        !mch_sig ~ "mCH not sig",
        edger_logFC < 0 ~ "mCH down",
        edger_logFC >= 0 ~ "mCH up"
      ),
      levels = c("mCH not sig", "mCH down", "mCH up")
    )
  )

cat("\nmCH direction breakdown:\n")
print(table(merged$mch_direction))

############################################################
# 7. Quadrant assignment
############################################################

quad_labels <- c(
  "Q1 (up / up)",
  "Q2 (down / up)",
  "Q3 (down / down)",
  "Q4 (up / down)"
)

assign_quadrant <- function(x_fold, y_fold) {
  case_when(
    x_fold >= 0 & y_fold >= 0 ~ quad_labels[1],
    x_fold <  0 & y_fold >= 0 ~ quad_labels[2],
    x_fold <  0 & y_fold <  0 ~ quad_labels[3],
    x_fold >= 0 & y_fold <  0 ~ quad_labels[4]
  )
}

merged$quadrant_mecp2 <- assign_quadrant(merged$mecp2_Fold, merged$edger_logFC)
merged$quadrant_h2a   <- assign_quadrant(merged$h2a_Fold, merged$edger_logFC)

############################################################
# 8. Correlation tests
############################################################

cor_mecp2 <- cor.test(
  merged$mecp2_Fold,
  merged$edger_logFC,
  method = "spearman"
)

cor_h2a <- cor.test(
  merged$h2a_Fold,
  merged$edger_logFC,
  method = "spearman"
)

cat("\nSpearman correlation: mCH vs MeCP2\n")
print(cor_mecp2)

cat("\nSpearman correlation: mCH vs H2AK119Ub\n")
print(cor_h2a)

rho_label_mecp2 <- paste0(
  "Spearman rho = ", round(cor_mecp2$estimate, 2),
  " | p = ", format.pval(cor_mecp2$p.value, digits = 2, eps = 1e-10)
)

rho_label_h2a <- paste0(
  "Spearman rho = ", round(cor_h2a$estimate, 2),
  " | p = ", format.pval(cor_h2a$p.value, digits = 2, eps = 1e-10)
)

############################################################
# 9. Quadrant statistics helper
############################################################

build_quadrant_counts <- function(df, quadrant_col, sig_col, x_fold_col, y_fold_col) {

  total_genes <- nrow(df)
  total_sig_genes <- sum(df[[sig_col]] == "Significant")

  quadrant_all <- df %>%
    count(.data[[quadrant_col]]) %>%
    complete(
      !!sym(quadrant_col) := quad_labels,
      fill = list(n = 0)
    ) %>%
    rename(n_total_genes = n) %>%
    mutate(
      percent_total_genes = round(100 * n_total_genes / total_genes, 1)
    )

  quadrant_sig <- df %>%
    filter(.data[[sig_col]] == "Significant") %>%
    count(.data[[quadrant_col]]) %>%
    complete(
      !!sym(quadrant_col) := quad_labels,
      fill = list(n = 0)
    ) %>%
    rename(n_sig_genes = n) %>%
    mutate(
      percent_sig_genes = if (total_sig_genes > 0) {
        round(100 * n_sig_genes / total_sig_genes, 1)
      } else {
        0
      }
    )

  quadrant_counts <- left_join(quadrant_all, quadrant_sig, by = quadrant_col) %>%
    mutate(
      label = paste0(
        "Total: n=", n_total_genes, " (", percent_total_genes, "%)\n",
        "Sig: n=", n_sig_genes, " (", percent_sig_genes, "%)"
      )
    )

  xmax <- max(abs(df[[x_fold_col]]), na.rm = TRUE)
  ymax <- max(abs(df[[y_fold_col]]), na.rm = TRUE)

  quadrant_counts$x <- c(0.90 * xmax, -0.30 * xmax, -0.30 * xmax, 0.90 * xmax)
  quadrant_counts$y <- c(0.90 * ymax, 0.90 * ymax, -0.90 * ymax, -0.90 * ymax)

  quadrant_counts
}

qc_mecp2 <- build_quadrant_counts(
  merged,
  quadrant_col = "quadrant_mecp2",
  sig_col = "mch_sig_vs_mecp2",
  x_fold_col = "mecp2_Fold",
  y_fold_col = "edger_logFC"
)

qc_h2a <- build_quadrant_counts(
  merged,
  quadrant_col = "quadrant_h2a",
  sig_col = "mch_sig_vs_h2a",
  x_fold_col = "h2a_Fold",
  y_fold_col = "edger_logFC"
)

cat("\nmCH vs MeCP2 quadrant counts:\n")
print(qc_mecp2 %>% select(-x, -y))

cat("\nmCH vs H2AK119Ub quadrant counts:\n")
print(qc_h2a %>% select(-x, -y))

############################################################
# 10. Scatter plot: mCH vs MeCP2
############################################################

top_mecp2_q1 <- merged %>%
  filter(
    mch_sig_vs_mecp2 == "Significant",
    quadrant_mecp2 == "Q1 (up / up)"
  ) %>%
  mutate(label_strength = abs(mecp2_Fold) + abs(edger_logFC)) %>%
  arrange(desc(label_strength)) %>%
  slice_head(n = 10)

top_mecp2_q3 <- merged %>%
  filter(
    mch_sig_vs_mecp2 == "Significant",
    quadrant_mecp2 == "Q3 (down / down)"
  ) %>%
  mutate(label_strength = abs(mecp2_Fold) + abs(edger_logFC)) %>%
  arrange(desc(label_strength)) %>%
  slice_head(n = 10)

top_genes_mecp2 <- bind_rows(top_mecp2_q1, top_mecp2_q3) %>%
  distinct(SYMBOL, .keep_all = TRUE)

merged_plot <- merged %>% arrange(mch_sig_vs_mecp2)

p_mecp2 <- ggplot(
  merged_plot,
  aes(x = mecp2_Fold, y = edger_logFC, color = mch_sig_vs_mecp2)
) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = c(
    "Not significant" = "grey85",
    "Significant" = "deeppink"
  )) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_text(
    data = qc_mecp2,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.5,
    fontface = "bold",
    color = "black"
  ) +
  geom_text_repel(
    data = top_genes_mecp2,
    aes(x = mecp2_Fold, y = edger_logFC, label = SYMBOL),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    fontface = "bold",
    color = "black"
  ) +
  labs(
    x = "MeCP2 log2(FC, BAP1-KO / ctrl)",
    y = "mCH log2FC (edgeR, BAP1-KO / ctrl)",
    color = "",
    title = "Non-CG methylation vs MeCP2 (gene-level)",
    subtitle = paste0(
      "n (total) = ", nrow(merged), " genes | ",
      "n (sig in both) = ", sum(merged$mch_sig_vs_mecp2 == "Significant"),
      " genes\n", rho_label_mecp2
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    plot.subtitle = element_text(size = 10, face = "bold")
  )

print(p_mecp2)

############################################################
# 11. Scatter plot: mCH vs H2AK119Ub
############################################################

top_h2a_q1 <- merged %>%
  filter(
    mch_sig_vs_h2a == "Significant",
    quadrant_h2a == "Q1 (up / up)"
  ) %>%
  mutate(label_strength = abs(h2a_Fold) + abs(edger_logFC)) %>%
  arrange(desc(label_strength)) %>%
  slice_head(n = 10)

top_h2a_q3 <- merged %>%
  filter(
    mch_sig_vs_h2a == "Significant",
    quadrant_h2a == "Q3 (down / down)"
  ) %>%
  mutate(label_strength = abs(h2a_Fold) + abs(edger_logFC)) %>%
  arrange(desc(label_strength)) %>%
  slice_head(n = 10)

top_genes_h2a <- bind_rows(top_h2a_q1, top_h2a_q3) %>%
  distinct(SYMBOL, .keep_all = TRUE)

merged_plot_h2a <- merged %>% arrange(mch_sig_vs_h2a)

p_h2a <- ggplot(
  merged_plot_h2a,
  aes(x = h2a_Fold, y = edger_logFC, color = mch_sig_vs_h2a)
) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = c(
    "Not significant" = "grey85",
    "Significant" = "deeppink"
  )) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_text(
    data = qc_h2a,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.5,
    fontface = "bold",
    color = "black"
  ) +
  geom_text_repel(
    data = top_genes_h2a,
    aes(x = h2a_Fold, y = edger_logFC, label = SYMBOL),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    fontface = "bold",
    color = "black"
  ) +
  labs(
    x = "H2AK119Ub log2(FC, BAP1-KO / ctrl)",
    y = "mCH log2FC (edgeR, BAP1-KO / ctrl)",
    color = "",
    title = "Non-CG methylation vs H2AK119Ub (gene-level)",
    subtitle = paste0(
      "n (total) = ", nrow(merged), " genes | ",
      "n (sig in both) = ", sum(merged$mch_sig_vs_h2a == "Significant"),
      " genes\n", rho_label_h2a
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    plot.subtitle = element_text(size = 10, face = "bold")
  )

print(p_h2a)

############################################################
# 12. Three-way scatter: MeCP2 vs H2AK119Ub colored by mCH
############################################################

n_mch_sig <- sum(merged$mch_sig)
n_mch_down <- sum(merged$mch_direction == "mCH down")
n_mch_up <- sum(merged$mch_direction == "mCH up")

three_way_xmax <- max(abs(merged$mecp2_Fold), na.rm = TRUE)
three_way_ymax <- max(abs(merged$h2a_Fold), na.rm = TRUE)

three_way_quad_summary <- merged %>%
  group_by(quadrant) %>%
  summarise(
    n_total = n(),
    n_mch_sig = sum(mch_sig),
    n_mch_down = sum(mch_direction == "mCH down"),
    n_mch_up = sum(mch_direction == "mCH up"),
    .groups = "drop"
  ) %>%
  complete(
    quadrant = c(
      "Q1 (MeCP2 up / H2A up)",
      "Q2 (MeCP2 down / H2A up)",
      "Q3 (MeCP2 down / H2A down)",
      "Q4 (MeCP2 up / H2A down)"
    ),
    fill = list(n_total = 0, n_mch_sig = 0, n_mch_down = 0, n_mch_up = 0)
  ) %>%
  mutate(
    label = paste0(
      "n=", n_total,
      "\nmCH sig: ", n_mch_sig,
      " (down:", n_mch_down,
      " up:", n_mch_up, ")"
    )
  )

three_way_quad_summary$x <- c(
  0.90 * three_way_xmax, -0.30 * three_way_xmax,
  -0.30 * three_way_xmax, 0.90 * three_way_xmax
)
three_way_quad_summary$y <- c(
  0.90 * three_way_ymax, 0.90 * three_way_ymax,
  -0.90 * three_way_ymax, -0.90 * three_way_ymax
)

cat("\nThree-way quadrant summary:\n")
print(three_way_quad_summary %>% select(-x, -y, -label))

merged_plot_3way <- merged %>% arrange(mch_direction)

p_three_way <- ggplot(
  merged_plot_3way,
  aes(x = mecp2_Fold, y = h2a_Fold, color = mch_direction)
) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_manual(values = c(
    "mCH not sig" = "grey85",
    "mCH down" = "steelblue",
    "mCH up" = "firebrick"
  )) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_text(
    data = three_way_quad_summary,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.5,
    fontface = "bold",
    color = "black"
  ) +
  labs(
    x = "MeCP2 log2(FC, BAP1-KO / ctrl)",
    y = "H2AK119Ub log2(FC, BAP1-KO / ctrl)",
    color = "",
    title = "MeCP2 vs H2AK119Ub colored by mCH status",
    subtitle = paste0(
      "n (total) = ", nrow(merged), " genes | ",
      "mCH sig: ", n_mch_sig,
      " (down: ", n_mch_down,
      ", up: ", n_mch_up, ")"
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    plot.subtitle = element_text(size = 10, face = "bold")
  )

print(p_three_way)

############################################################
# 13. Diagnostic: gene length vs mCH significance
############################################################

merged_length <- merged %>%
  filter(!is.na(mch_gene_length)) %>%
  mutate(
    mch_sig_label = factor(
      ifelse(mch_sig, "mCH sig", "mCH not sig"),
      levels = c("mCH not sig", "mCH sig")
    )
  )

wilcox_mch_length <- wilcox.test(
  mch_gene_length ~ mch_sig_label,
  data = merged_length
)

cat("\nGene length vs mCH significance (Wilcoxon):\n")
print(wilcox_mch_length)

length_summary <- merged_length %>%
  group_by(mch_sig_label) %>%
  summarise(
    n = n(),
    median_kb = round(median(mch_gene_length / 1000), 1),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("n = ", n, "\nmedian = ", median_kb, " kb"))

p_mch_length <- ggplot(
  merged_length,
  aes(x = mch_sig_label, y = mch_gene_length / 1000, fill = mch_sig_label)
) +
  geom_violin(alpha = 0.6, show.legend = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, show.legend = FALSE) +
  geom_text(
    data = length_summary,
    aes(
      x = mch_sig_label,
      y = max(merged_length$mch_gene_length / 1000) * 1.5,
      label = label
    ),
    size = 4, fontface = "bold", vjust = 1, show.legend = FALSE
  ) +
  scale_fill_manual(values = c("mCH not sig" = "grey70", "mCH sig" = "steelblue")) +
  scale_y_log10() +
  labs(
    x = "mCH significance status",
    y = "Gene length (kb, log10)",
    title = "Gene length vs mCH significance",
    subtitle = paste0("Wilcoxon p = ", format.pval(wilcox_mch_length$p.value, digits = 2))
  ) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

print(p_mch_length)

############################################################
# 14. Diagnostic: coverage vs mCH significance
############################################################

merged_cov <- merged %>%
  filter(!is.na(total_coverage_mean)) %>%
  mutate(
    mch_sig_label = factor(
      ifelse(mch_sig, "mCH sig", "mCH not sig"),
      levels = c("mCH not sig", "mCH sig")
    )
  )

wilcox_mch_cov <- wilcox.test(
  total_coverage_mean ~ mch_sig_label,
  data = merged_cov
)

cat("\nCoverage vs mCH significance (Wilcoxon):\n")
print(wilcox_mch_cov)

cov_summary <- merged_cov %>%
  group_by(mch_sig_label) %>%
  summarise(
    n = n(),
    median_cov = round(median(total_coverage_mean / 1e6), 2),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("n = ", n, "\nmedian = ", median_cov, "M"))

p_mch_cov <- ggplot(
  merged_cov,
  aes(x = mch_sig_label, y = total_coverage_mean / 1e6, fill = mch_sig_label)
) +
  geom_violin(alpha = 0.6, show.legend = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, show.legend = FALSE) +
  geom_text(
    data = cov_summary,
    aes(
      x = mch_sig_label,
      y = max(merged_cov$total_coverage_mean / 1e6) * 1.5,
      label = label
    ),
    size = 4, fontface = "bold", vjust = 1, show.legend = FALSE
  ) +
  scale_fill_manual(values = c("mCH not sig" = "grey70", "mCH sig" = "steelblue")) +
  scale_y_log10() +
  labs(
    x = "mCH significance status",
    y = "Total coverage (millions, log10)",
    title = "Total coverage vs mCH significance",
    subtitle = paste0("Wilcoxon p = ", format.pval(wilcox_mch_cov$p.value, digits = 2))
  ) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

print(p_mch_cov)

############################################################
# 15. Three-way GO enrichment
############################################################

universe_symbols <- unique(merged$SYMBOL)
universe_entrez <- bitr(
  universe_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

cat("\nGO universe size (merged genes with Entrez IDs):", nrow(universe_entrez), "\n")

v2_quadrant_names <- c(
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

gene_length_for_goseq <- merged %>%
  select(SYMBOL, mch_gene_length) %>%
  inner_join(universe_entrez, by = "SYMBOL") %>%
  filter(!is.na(mch_gene_length)) %>%
  group_by(ENTREZID) %>%
  slice_head(n = 1) %>%
  ungroup()

for (i in seq_along(v2_quadrant_names)) {
  q_name <- v2_quadrant_names[i]
  q_short <- quadrant_short[i]

  q_genes <- merged %>%
    filter(
      sig == "Significant",
      mch_sig,
      quadrant == q_name
    ) %>%
    pull(SYMBOL)

  cat("\n", q_short, ": ", length(q_genes),
      " three-way significant genes (ChIP sig + mCH sig)\n", sep = "")

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

    n_terms <- if (!is.null(go_res)) nrow(as.data.frame(go_res)) else 0
    cat("  enrichGO ", ont, ": ", n_terms, " terms\n", sep = "")
  }

  gene_vector <- as.integer(gene_length_for_goseq$ENTREZID %in% q_entrez$ENTREZID)
  names(gene_vector) <- gene_length_for_goseq$ENTREZID
  bias_data <- gene_length_for_goseq$mch_gene_length
  names(bias_data) <- gene_length_for_goseq$ENTREZID

  pwf <- nullp(gene_vector, bias.data = bias_data, plot.fit = FALSE)

  goseq_res <- goseq(pwf, "mm10", "knownGene")
  goseq_res$over_represented_padj <- p.adjust(
    goseq_res$over_represented_pvalue, method = "BH"
  )
  goseq_sig <- goseq_res %>% filter(over_represented_padj < 0.05)

  goseq_results[[q_short]] <- goseq_res

  cat("  goseq: ", nrow(goseq_sig), " significant terms (length-corrected)\n", sep = "")
}

############################################################
# 16. Save outputs
############################################################

cat("\n--- Saving outputs ---\n")

write.csv(
  merged,
  file.path(out_dir, "mch_mecp2_h2aub_merged_v2.csv"),
  row.names = FALSE
)

write.csv(
  qc_mecp2 %>% select(-x, -y),
  file.path(out_dir, "quadrant_counts_mecp2_v2.csv"),
  row.names = FALSE
)

write.csv(
  qc_h2a %>% select(-x, -y),
  file.path(out_dir, "quadrant_counts_h2a_v2.csv"),
  row.names = FALSE
)

write.csv(
  three_way_quad_summary %>% select(-x, -y, -label),
  file.path(out_dir, "three_way_quadrant_summary.csv"),
  row.names = FALSE
)

save_multiformat_ggplot(
  p_mecp2,
  file.path(out_dir, "mch_vs_mecp2_v2"),
  width = 10, height = 10
)

save_multiformat_ggplot(
  p_h2a,
  file.path(out_dir, "mch_vs_h2aub_v2"),
  width = 10, height = 10
)

save_multiformat_ggplot(
  p_three_way,
  file.path(out_dir, "three_way_scatter_v2"),
  width = 10, height = 10
)

save_multiformat_ggplot(
  p_mch_length,
  file.path(out_dir, "diagnostic_mch_genelength"),
  width = 7, height = 6
)

save_multiformat_ggplot(
  p_mch_cov,
  file.path(out_dir, "diagnostic_mch_coverage"),
  width = 7, height = 6
)

for (key in names(go_results)) {
  go_res <- go_results[[key]]
  go_df <- as.data.frame(go_res)

  write.csv(
    go_df,
    file.path(out_dir, paste0("GO_three_way_", key, ".csv")),
    row.names = FALSE
  )

  if (nrow(go_df) > 0) {
    p_go <- dotplot(go_res, showCategory = 20, title = paste0("Three-way GO ", key))
    save_multiformat_ggplot(
      p_go,
      file.path(out_dir, paste0("GO_three_way_", key, "_dotplot")),
      width = 8, height = 10
    )
  }
}

for (q_short in names(goseq_results)) {
  goseq_res <- goseq_results[[q_short]]
  write.csv(
    goseq_res,
    file.path(out_dir, paste0("GO_goseq_three_way_", q_short, ".csv")),
    row.names = FALSE
  )
}

cat("\nDone.\n")
