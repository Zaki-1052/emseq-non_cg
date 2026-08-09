# scripts/plot_mca_volcano.R
# Volcano plots for gene-body mCA differential methylation (mca_differential_results.tsv).
#
# Produces four figures, each written as PDF + SVG + PNG:
#   {prefix}_direction   x = mCA difference, y = -log10(p), colour = direction of change
#   {prefix}_robustness  same coordinates, colour = number of leave-one-out refits retaining significance
#   {prefix}_genelength  same coordinates, colour = log10(gene length)
#   {prefix}_panel       the three above combined
#
# Usage:
#   Rscript scripts/plot_mca_volcano.R --results mca_results/mca_differential_results.tsv

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(svglite)
})

# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

option_list <- list(
  make_option("--results", type = "character", default = NULL,
              help = "Path to mca_differential_results.tsv [required]"),
  make_option("--outdir", type = "character", default = "mca_results/plots",
              help = "Directory for output figures [default: %default]"),
  make_option("--prefix", type = "character", default = "mca_volcano",
              help = "Output filename prefix [default: %default]"),
  make_option("--alpha", type = "double", default = 0.05,
              help = "Significance level for both Bonferroni and BH thresholds [default: %default]"),
  make_option("--sig-basis", type = "character", default = "bonferroni", dest = "sig_basis",
              help = "Which correction defines the coloured/labelled gene set: bonferroni or fdr [default: %default]"),
  make_option("--label-top", type = "integer", default = 20, dest = "label_top",
              help = "Number of most-significant genes to label [default: %default]"),
  make_option("--gene-list", type = "character", default = NULL, dest = "gene_list",
              help = "Optional file of gene names (one per line); significant members are additionally labelled"),
  make_option("--set-name", type = "character", default = "Gene set", dest = "set_name",
              help = "Name of the --gene-list set, used in the caption [default: %default]"),
  make_option("--width", type = "double", default = 10,
              help = "Width of each single-panel figure, inches [default: %default]"),
  make_option("--height", type = "double", default = 8,
              help = "Height of each single-panel figure, inches [default: %default]"),
  make_option("--dpi", type = "integer", default = 300,
              help = "Raster resolution for PNG output [default: %default]"),
  make_option("--subfolders", action = "store_true", default = FALSE,
              help = "Write each figure into its own subfolder (matches the biomodal save_multiformat_ggplot layout)")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  description = "Volcano plots for gene-body mCA differential methylation results."
))

if (is.null(opt$results)) stop("--results is required")
if (!file.exists(opt$results)) stop("--results file not found: ", opt$results)
if (!opt$sig_basis %in% c("bonferroni", "fdr")) {
  stop("--sig-basis must be 'bonferroni' or 'fdr', got: ", opt$sig_basis)
}

# =============================================================================
# STYLE
# =============================================================================

# Mirrors theme_biomodal() from the biomodal downstream pipeline
# (scripts/viz_sections/_shared_config.R) so figures from the two pipelines match.
theme_mca <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, size = base_size),
      plot.caption = element_text(hjust = 0, size = base_size - 3, colour = "grey30"),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold")
    )
}

# Diverging pair carried over from COLORS$direction in the biomodal pipeline.
COLOR_DIRECTION <- c(
  "Higher in mutant" = "#D7191C",
  "Lower in mutant"  = "#2C7BB6",
  "Not significant"  = "grey70"
)

# Single-hue sequential ramps: light = low, dark = high.
RAMP_ROBUSTNESS <- c("#807DBA", "#3F007D")
RAMP_GENELENGTH <- c("#74C476", "#00441B")

save_multiformat <- function(plot, base_path, width, height, dpi, subfolders) {
  figure_name <- basename(base_path)
  if (subfolders) {
    output_dir <- file.path(dirname(base_path), figure_name)
    file_prefix <- file.path(output_dir, figure_name)
  } else {
    output_dir <- dirname(base_path)
    file_prefix <- base_path
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ggsave(paste0(file_prefix, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(file_prefix, ".svg"), plot, width = width, height = height, device = svglite::svglite)
  ggsave(paste0(file_prefix, ".png"), plot, width = width, height = height, dpi = dpi, device = "png")

  cat(sprintf("  Saved: %s.{pdf,svg,png}\n", figure_name))
  invisible(plot)
}

# =============================================================================
# DATA
# =============================================================================

REQUIRED_COLS <- c("gene_name", "gene_length", "mca_diff", "dss_pval",
                   "p_bonferroni", "p_bh", "sig_bonferroni", "sig_fdr", "loo_n_sig")

load_results <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  missing <- setdiff(REQUIRED_COLS, names(df))
  if (length(missing) > 0) {
    stop("Results file is missing required columns: ", paste(missing, collapse = ", "))
  }
  if (any(!is.finite(df$dss_pval)) || any(df$dss_pval <= 0)) {
    stop("dss_pval contains non-finite or non-positive values; -log10 transform would be undefined")
  }
  df
}

# The number of leave-one-out refits equals the number of samples, read off the
# per-sample mCA columns rather than assumed.
count_samples <- function(df) {
  n <- sum(grepl("_mca$", names(df)) & !grepl("^mca_", names(df)))
  if (n == 0) stop("No per-sample '*_mca' columns found; cannot determine leave-one-out refit count")
  n
}

results <- load_results(opt$results)
n_samples <- count_samples(results)
n_tested <- nrow(results)

sig_col <- if (opt$sig_basis == "bonferroni") "sig_bonferroni" else "sig_fdr"

results <- results %>%
  mutate(
    mca_diff_pct = mca_diff * 100,
    neg_log10_p  = -log10(dss_pval),
    log10_length = log10(gene_length),
    is_sig       = .data[[sig_col]],
    loo_robust   = is_sig & loo_n_sig == n_samples,
    direction = case_when(
      !is_sig          ~ "Not significant",
      mca_diff > 0     ~ "Higher in mutant",
      TRUE             ~ "Lower in mutant"
    ),
    # Identity scales: the significant points are drawn larger and more opaque
    # than the background cloud without spending a legend on it.
    pt_size  = ifelse(is_sig, 1.9, 0.6),
    pt_alpha = ifelse(is_sig, 0.9, 0.30)
  ) %>%
  # Row order controls z-order within a layer: the grey cloud is drawn first.
  arrange(is_sig, loo_robust)

# Raw-p cutoffs corresponding to each correction, so both can be drawn as
# horizontal lines on a raw-p axis. BH is a step-up procedure, so its cutoff is
# data-dependent and read back off the genes that pass.
bonferroni_cut <- opt$alpha / n_tested
bh_cut <- if (any(results$sig_fdr)) max(results$dss_pval[results$sig_fdr]) else NA_real_

# =============================================================================
# LABEL SET
# =============================================================================

label_genes <- results$gene_name[order(results$dss_pval)][seq_len(min(opt$label_top, n_tested))]
label_genes <- union(label_genes, results$gene_name[results$loo_robust])

gene_list_hits <- character(0)
if (!is.null(opt$gene_list)) {
  if (!file.exists(opt$gene_list)) stop("--gene-list file not found: ", opt$gene_list)
  wanted <- readLines(opt$gene_list)
  wanted <- trimws(wanted[nzchar(trimws(wanted))])
  gene_list_hits <- intersect(wanted, results$gene_name[results$is_sig])
  label_genes <- union(label_genes, gene_list_hits)
  cat(sprintf("Gene list '%s': %d of %d names matched a significant gene\n",
              opt$gene_list, length(gene_list_hits), length(wanted)))
}

results$label <- ifelse(results$gene_name %in% label_genes, results$gene_name, "")

# =============================================================================
# SHARED VOLCANO GEOMETRY
# =============================================================================

n_up     <- sum(results$is_sig & results$mca_diff > 0)
n_down   <- sum(results$is_sig & results$mca_diff < 0)
n_robust <- sum(results$loo_robust)

basis_label <- if (opt$sig_basis == "bonferroni") "Bonferroni" else "BH FDR"

caption_lines <- c(
  sprintf("%d genes tested. Significance basis: %s p < %g (%d genes). Raw p on the y-axis; horizontal lines mark the corrected cutoffs.",
          n_tested, basis_label, opt$alpha, sum(results$is_sig)),
  sprintf("Leave-one-out: each gene was refit %d times, dropping one sample each time. Ringed points retain significance in all %d refits.",
          n_samples, n_samples)
)
if (length(gene_list_hits) > 0) {
  caption_lines <- c(caption_lines,
                     sprintf("%s: %d significant members labelled.", opt$set_name, length(gene_list_hits)))
}
volcano_caption <- paste(caption_lines, collapse = "\n")

# Threshold lines, direct-labelled at the left edge so they are not line-style-only.
x_min <- min(results$mca_diff_pct)
threshold_layers <- function() {
  layers <- list(
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4),
    geom_hline(yintercept = -log10(bonferroni_cut), linetype = "dashed",
               colour = "grey40", linewidth = 0.4),
    annotate("text", x = x_min, y = -log10(bonferroni_cut), hjust = 0, vjust = -0.6,
             label = sprintf("Bonferroni %g", opt$alpha), size = 3, colour = "grey30")
  )
  if (!is.na(bh_cut)) {
    layers <- c(layers, list(
      geom_hline(yintercept = -log10(bh_cut), linetype = "dotted",
                 colour = "grey40", linewidth = 0.4),
      annotate("text", x = x_min, y = -log10(bh_cut), hjust = 0, vjust = -0.6,
               label = sprintf("BH FDR %g", opt$alpha), size = 3, colour = "grey30")
    ))
  }
  layers
}

# Ring separates the leave-one-out-robust genes from the crowd without spending a
# colour channel on them. Both of these layers set inherit.aes = FALSE and carry
# their own data and mapping: panels B and C build from a base ggplot() that has
# no data, and panel A's base mapping would otherwise collide with the fixed
# colour/size params set here.
robust_ring_layer <- function() {
  geom_point(data = filter(results, loo_robust),
             aes(x = mca_diff_pct, y = neg_log10_p), inherit.aes = FALSE,
             shape = 21, fill = NA, colour = "grey15", size = 2.6, stroke = 0.6)
}

gene_label_layer <- function() {
  geom_text_repel(data = filter(results, nzchar(label)),
                  aes(x = mca_diff_pct, y = neg_log10_p, label = label),
                  inherit.aes = FALSE,
                  size = 3, colour = "grey15",
                  max.overlaps = 30, box.padding = 0.45, min.segment.length = 0.2,
                  segment.colour = "grey55", segment.size = 0.3, seed = 1)
}

volcano_scaffold <- function(p, title, subtitle) {
  p +
    threshold_layers() +
    robust_ring_layer() +
    gene_label_layer() +
    labs(
      title = title,
      subtitle = subtitle,
      x = "mCA difference (mutant − control, percentage points)",
      y = expression(-log[10] * (italic(p))),
      caption = volcano_caption
    ) +
    theme_mca()
}

# =============================================================================
# PANEL A: DIRECTION
# =============================================================================

p_direction <- ggplot(results, aes(x = mca_diff_pct, y = neg_log10_p)) +
  geom_point(aes(colour = direction, size = pt_size, alpha = pt_alpha)) +
  scale_colour_manual(values = COLOR_DIRECTION, name = "Direction",
                      breaks = names(COLOR_DIRECTION)) +
  scale_size_identity() +
  scale_alpha_identity() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

p_direction <- volcano_scaffold(
  p_direction,
  "Gene-body mCA differential methylation",
  sprintf("BAP1-KO vs control  |  %d higher, %d lower in mutant (%s p < %g)  |  %d robust to leave-one-out",
          n_up, n_down, basis_label, opt$alpha, n_robust)
)

# =============================================================================
# PANEL B: LEAVE-ONE-OUT ROBUSTNESS
# =============================================================================

p_robust <- ggplot(mapping = aes(x = mca_diff_pct, y = neg_log10_p)) +
  geom_point(data = filter(results, !is_sig), colour = "grey70", size = 0.6, alpha = 0.30) +
  geom_point(data = filter(results, is_sig), aes(colour = loo_n_sig), size = 1.9, alpha = 0.9) +
  scale_colour_gradient(low = RAMP_ROBUSTNESS[1], high = RAMP_ROBUSTNESS[2],
                        limits = c(0, n_samples), breaks = seq(0, n_samples, 2),
                        name = "Refits still\nsignificant")

p_robust <- volcano_scaffold(
  p_robust,
  "mCA volcano coloured by leave-one-out robustness",
  sprintf("Only %d of %d significant genes survive all %d refits; grey points are not significant",
          n_robust, sum(results$is_sig), n_samples)
)

# =============================================================================
# PANEL C: GENE LENGTH
# =============================================================================

median_len_sig <- median(results$gene_length[results$is_sig])
median_len_all <- median(results$gene_length)

p_length <- ggplot(mapping = aes(x = mca_diff_pct, y = neg_log10_p)) +
  geom_point(data = filter(results, !is_sig), colour = "grey70", size = 0.6, alpha = 0.30) +
  geom_point(data = filter(results, is_sig), aes(colour = log10_length), size = 1.9, alpha = 0.9) +
  scale_colour_gradient(low = RAMP_GENELENGTH[1], high = RAMP_GENELENGTH[2],
                        name = expression(log[10] * "(gene length, bp)"))

p_length <- volcano_scaffold(
  p_length,
  "mCA volcano coloured by gene length",
  sprintf("Median length %s kb among significant genes vs %s kb across all tested genes",
          format(round(median_len_sig / 1000, 1), nsmall = 1),
          format(round(median_len_all / 1000, 1), nsmall = 1))
)

# =============================================================================
# WRITE
# =============================================================================

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
out <- function(suffix) file.path(opt$outdir, paste0(opt$prefix, "_", suffix))

cat("\nWriting figures to ", normalizePath(opt$outdir), "\n", sep = "")
save_multiformat(p_direction, out("direction"), opt$width, opt$height, opt$dpi, opt$subfolders)
save_multiformat(p_robust,    out("robustness"), opt$width, opt$height, opt$dpi, opt$subfolders)
save_multiformat(p_length,    out("genelength"), opt$width, opt$height, opt$dpi, opt$subfolders)

p_panel <- p_direction / (p_robust | p_length) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold")))
save_multiformat(p_panel, out("panel"), opt$width * 1.6, opt$height * 1.9, opt$dpi, opt$subfolders)

# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

cat("\n")
cat("Genes tested:            ", n_tested, "\n", sep = "")
cat("Significance basis:      ", basis_label, " p < ", opt$alpha, "\n", sep = "")
cat("Significant:             ", sum(results$is_sig),
    "  (", n_up, " higher, ", n_down, " lower in mutant)\n", sep = "")
cat("Robust to leave-one-out: ", n_robust, "\n", sep = "")
cat("Bonferroni raw-p cutoff: ", format(bonferroni_cut, digits = 4), "\n", sep = "")
cat("BH raw-p cutoff:         ",
    if (is.na(bh_cut)) "none (no gene passes)" else format(bh_cut, digits = 4), "\n", sep = "")
cat("Largest -log10(p):       ", format(max(results$neg_log10_p), digits = 4), "\n", sep = "")
cat("Genes labelled:          ", sum(nzchar(results$label)), "\n", sep = "")
