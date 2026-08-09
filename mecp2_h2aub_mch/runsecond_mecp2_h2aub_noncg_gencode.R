############################################################
# Non-CG methylation (mCH) vs MeCP2 / H2AK119Ub gene body analysis
#
# Flow:
# gene_df (MeCP2 + H2AK119Ub gene-level results, already saved as csv)
#      |
#      v
# mch_differential_results.tsv (non-CG edgeR results, gene-level)
#      |
#      v
# Merge on gene name
#      |
#      v
# Correlation + quadrant plots: mCH vs MeCP2, mCH vs H2AK119Ub
############################################################

#libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

############################################################
# Load inputs
############################################################

# gene-level MeCP2 + H2AK119Ub results saved by the earlier script
gene_df <- read_csv(
  "/Users/janviramchandra/Desktop/gencode_gene_bodies/MeCP2_H2AK119Ub_gene_level_results.csv"
)

# non-CG methylation differential results (edgeR, gene-level)
mch <- read_tsv(
  "/Users/janviramchandra/Desktop/ca/mch_differential_results.tsv"
)

############################################################
# Sanity check gene name overlap before trusting the join
############################################################

cat("gene_df SYMBOL examples:\n")
print(head(sort(gene_df$SYMBOL)))

cat("\nmch gene_name examples:\n")
print(head(sort(mch$gene_name)))

n_overlap <- sum(gene_df$SYMBOL %in% mch$gene_name)
cat(
  "\nGenes in gene_df also found in mch$gene_name:",
  n_overlap, "/", nrow(gene_df), "\n"
)

# If n_overlap is low relative to nrow(gene_df), check casing or
# whether mch$gene_id (ENSMUSG...) is a better join key than gene_name.

############################################################
# Merge
############################################################

merged <- inner_join(
  gene_df,
  mch,
  by = c("SYMBOL" = "gene_name")
)

cat(
  "\nGenes with MeCP2, H2AK119Ub, and mCH data:",
  nrow(merged), "\n"
)

############################################################
# Correlation tests
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
  "Spearman \u03C1 = ", round(cor_mecp2$estimate, 2),
  " | p = ", format.pval(cor_mecp2$p.value, digits = 2, eps = 1e-10)
)

rho_label_h2a <- paste0(
  "Spearman \u03C1 = ", round(cor_h2a$estimate, 2),
  " | p = ", format.pval(cor_h2a$p.value, digits = 2, eps = 1e-10)
)

############################################################
# Significance flags for coloring points
#
# Pink = significant in BOTH datasets being plotted, not just mCH alone.
# mecp2_sig / h2a_sig already come from gene_df (mutant/control peak FDR).
# sig_fdr005 is the mCH edgeR FDR.
############################################################

merged$mch_sig_vs_mecp2 <- factor(
  ifelse(merged$sig_fdr005 & merged$mecp2_sig, "Significant", "Not significant"),
  levels = c("Not significant", "Significant")
)

merged$mch_sig_vs_h2a <- factor(
  ifelse(merged$sig_fdr005 & merged$h2a_sig, "Significant", "Not significant"),
  levels = c("Not significant", "Significant")
)

############################################################
# Quadrant assignment
# (mCH log2FC on y-axis for both plots; ChIP fold change on x-axis)
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
# Quadrant count/label helper
#
# Builds a data frame with, per quadrant:
#   n_total_genes / percent_total_genes  (% of all genes in the plot)
#   n_sig_genes   / percent_sig_genes    (% of all significant genes in the plot)
# and a combined text label plus x/y placement coordinates.
############################################################

build_quadrant_counts <- function(df, quadrant_col, sig_col, x_fold_col, y_fold_col) {
  
  total_genes <- nrow(df)
  total_sig_genes <- sum(df[[sig_col]] == "Significant")
  
  quadrant_all <- df %>%
    count(.data[[quadrant_col]]) %>%
    complete(
      !!quadrant_col := quad_labels,
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
      !!quadrant_col := quad_labels,
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
  
  quadrant_counts$x <- c(0.95 * xmax, -0.30 * xmax, -0.30 * xmax, 0.90 * xmax)
  quadrant_counts$y <- c(0.10 * ymax, 0.95 * ymax, -0.90 * ymax, -0.90 * ymax)
  
  quadrant_counts
}

############################################################
# Calculate quadrant counts for plotting
############################################################

quadrant_counts_mecp2 <- build_quadrant_counts(
  merged,
  quadrant_col = "quadrant_mecp2",
  sig_col = "mch_sig_vs_mecp2",
  x_fold_col = "mecp2_Fold",
  y_fold_col = "edger_logFC"
)

quadrant_counts_h2a <- build_quadrant_counts(
  merged,
  quadrant_col = "quadrant_h2a",
  sig_col = "mch_sig_vs_h2a",
  x_fold_col = "h2a_Fold",
  y_fold_col = "edger_logFC"
)

############################################################
# Plot: mCH vs MeCP2
############################################################

# Top 10 significant genes in Q1 and Q3 — MeCP2
top_genes_mecp2_q1 <- merged %>%
  filter(
    mch_sig_vs_mecp2 == "Significant",
    quadrant_mecp2 == "Q1 (up / up)"
  ) %>%
  mutate(
    label_strength = abs(mecp2_Fold) + abs(edger_logFC)
  ) %>%
  arrange(desc(label_strength)) %>%
  slice_head(n = 10)

top_genes_mecp2_q3 <- merged %>%
  filter(
    mch_sig_vs_mecp2 == "Significant",
    quadrant_mecp2 == "Q3 (down / down)"
  ) %>%
  mutate(
    label_strength = abs(mecp2_Fold) + abs(edger_logFC)
  ) %>%
  arrange(desc(label_strength)) %>%
  slice_head(n = 10)

top_genes_mecp2 <- bind_rows(
  top_genes_mecp2_q1,
  top_genes_mecp2_q3
)

p_mecp2 <- ggplot(
  merged,
  aes(x = mecp2_Fold, y = edger_logFC, color = mch_sig_vs_mecp2)
) +
  
  geom_point(size = 1, alpha = 0.8) +
  
  scale_color_manual(
    values = c(
      "Not significant" = "grey85",
      "Significant" = "deeppink"
    )
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_text(
    data = quadrant_counts_mecp2,
    aes(
      x = x,
      y = y,
      label = paste0(
        quadrant_mecp2, "\n",
        "Total: n=", n_total_genes,
        " (", percent_total_genes, "%)\n",
        "Sig: n=", n_sig_genes,
        " (", percent_sig_genes, "%)"
      )
    ),
    inherit.aes = FALSE,
    size = 2,
    fontface = "bold",
    color = "black"
  ) +
  
  geom_text_repel(
    data = top_genes_mecp2,
    aes(
      x = mecp2_Fold,
      y = edger_logFC,
      label = SYMBOL
    ),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    fontface = "bold",
    color = "black"
  ) +
  labs(
    x = "MeCP2 log2(Fold Change, mutant/control)",
    y = "Non-CG methylation log2FC (edgeR, mutant/control)",
    color = "",
    title = "Gene Bodies Non-CG Methylation vs MeCP2",
    subtitle = paste0(
      "n (total) = ", nrow(merged),
      " genes | n (sig in both) = ", sum(merged$mch_sig_vs_mecp2 == "Significant"),
      " genes\n", rho_label_mecp2
    )
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    plot.subtitle = element_text(size = 12, face = "bold")
  )

print(p_mecp2)

############################################################
# Plot: mCH vs H2AK119Ub
############################################################

# Top 10 significant genes in Q1 and Q3 — H2AK119Ub

top_genes_h2a_q1 <- merged %>%
  filter(
    mch_sig_vs_h2a == "Significant",
    quadrant_h2a == "Q1 (up / up)"
  ) %>%
  mutate(
    label_strength = abs(h2a_Fold) + abs(edger_logFC)
  ) %>%
  arrange(desc(label_strength)) %>%
  slice_head(n = 10)

top_genes_h2a_q3 <- merged %>%
  filter(
    mch_sig_vs_h2a == "Significant",
    quadrant_h2a == "Q3 (down / down)"
  ) %>%
  mutate(
    label_strength = abs(h2a_Fold) + abs(edger_logFC)
  ) %>%
  arrange(desc(label_strength)) %>%
  slice_head(n = 10)

top_genes_h2a <- bind_rows(
  top_genes_h2a_q1,
  top_genes_h2a_q3
)

p_h2a <- ggplot(
  merged,
  aes(
    x = h2a_Fold,
    y = edger_logFC,
    color = mch_sig_vs_h2a
  )
) +
  
  geom_point(
    size = 1,
    alpha = 0.8
  ) +
  
  scale_color_manual(
    values = c(
      "Not significant" = "grey85",
      "Significant" = "deeppink"
    )
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  
  # ---- Quadrant counts ----
geom_text(
  data = quadrant_counts_h2a,
  aes(
    x = x,
    y = y,
    label = paste0(
      quadrant_h2a, "\n",
      "Total: n=", n_total_genes,
      " (", percent_total_genes, "%)\n",
      "Sig: n=", n_sig_genes,
      " (", percent_sig_genes, "%)"
    )
  ),
  inherit.aes = FALSE,
  size = 2,
  fontface = "bold",
  color = "black"
) +
  
  # ---- Top 10 Q1 + Top 10 Q3 genes ----
geom_text_repel(
  data = top_genes_h2a,
  aes(
    x = h2a_Fold,
    y = edger_logFC,
    label = SYMBOL
  ),
  size = 3,
  max.overlaps = Inf,
  box.padding = 0.5,
  point.padding = 0.3,
  fontface = "bold",
  color = "black"
) +
  
  labs(
    x = "H2AK119Ub log2(Fold Change, mutant/control)",
    y = "Non-CG methylation log2FC (edgeR, mutant/control)",
    color = "",
    title = "Gene Bodies Non-CG Methylation vs H2AK119Ub",
    subtitle = paste0(
      "n (total) = ", nrow(merged),
      " genes | n (sig in both) = ",
      sum(merged$mch_sig_vs_h2a == "Significant"),
      " genes\n",
      rho_label_h2a
    )
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    plot.subtitle = element_text(
      size = 12,
      face = "bold"
    )
  )

print(p_h2a)
############################################################
# Save outputs
############################################################

setwd("/Users/janviramchandra/Desktop/ca")

write.csv(
  merged,
  "mCH_MeCP2_H2AK119Ub_merged_gene_level.csv",
  row.names = FALSE
)

ggsave(
  "mCH_vs_MeCP2_gene_level_plot.png",
  p_mecp2,
  width = 10,
  height = 10,
  dpi = 300
)
ggsave(
  "mCH_vs_MeCP2_gene_level_plot.svg",
  p_mecp2,
  width = 10,
  height = 10
)

ggsave(
  "mCH_vs_H2AK119Ub_gene_level_plot.png",
  p_h2a,
  width = 10,
  height = 10,
  dpi = 300
)
ggsave(
  "mCH_vs_H2AK119Ub_gene_level_plot.svg",
  p_h2a,
  width = 10,
  height = 10
)