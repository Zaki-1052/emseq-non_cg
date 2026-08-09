# scripts/05_mch_volcano.R
# Volcano plots for gene-body mCH differential methylation (mch_differential_results.tsv).
#
# Produces four figures per context, each written as PDF + SVG + PNG:
#   {prefix}_direction       x = logFC, y = -log10(p), colour = direction of change
#   {prefix}_fdr_stringency  same coordinates, colour = -log10(FDR)
#   {prefix}_genelength      same coordinates, colour = log10(gene length)
#   {prefix}_panel           the three above combined
#
# Supports running both all-CH and CA-only contexts in a single invocation
# via --results (required) and --results-ca (optional).
#
# Usage:
#   Rscript scripts/05_mch_volcano.R \
#       --results results/03_differential/mch_differential_results.tsv
#
#   Rscript scripts/05_mch_volcano.R \
#       --results results/03_differential/mch_differential_results.tsv \
#       --results-ca results/ca/03_differential/mch_differential_results.tsv

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(svglite)
})

script_dir <- if (interactive()) "scripts" else dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "utils", "multi_format_output.R"))

# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

option_list <- list(
  make_option("--results", type = "character", default = NULL,
              help = "Path to all-CH results/03_differential/mch_differential_results.tsv [required]"),
  make_option("--results-ca", type = "character", default = NULL, dest = "results_ca",
              help = "Path to CA-only results/ca/03_differential/mch_differential_results.tsv [optional; produces a second set of figures]"),
  make_option("--outdir", type = "character", default = "results/plots",
              help = "Output directory for all-CH figures [default: %default]"),
  make_option("--outdir-ca", type = "character", default = "results/ca/plots", dest = "outdir_ca",
              help = "Output directory for CA-only figures [default: %default]"),
  make_option("--prefix", type = "character", default = "05_mch_volcano",
              help = "Output filename prefix for all-CH [default: %default]"),
  make_option("--prefix-ca", type = "character", default = "05_mca_volcano", dest = "prefix_ca",
              help = "Output filename prefix for CA-only [default: %default]"),
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
              help = "Raster resolution for PNG output [default: %default]")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  description = "Volcano plots for gene-body mCH differential methylation results."
))

if (is.null(opt$results)) stop("--results is required")
if (!file.exists(opt$results)) stop("--results file not found: ", opt$results)
if (!is.null(opt$results_ca) && !file.exists(opt$results_ca)) {
  stop("--results-ca file not found: ", opt$results_ca)
}
if (!opt$sig_basis %in% c("bonferroni", "fdr")) {
  stop("--sig-basis must be 'bonferroni' or 'fdr', got: ", opt$sig_basis)
}

# =============================================================================
# STYLE
# =============================================================================

theme_volcano <- function(base_size = 12) {
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

COLOR_DIRECTION <- c(
  "Higher in mutant" = "#D7191C",
  "Lower in mutant"  = "#2C7BB6",
  "Not significant"  = "grey70"
)

RAMP_FDR_STRINGENCY <- c("#807DBA", "#3F007D")
RAMP_GENELENGTH     <- c("#74C476", "#00441B")

# =============================================================================
# BUILD VOLCANO SET
# =============================================================================

build_volcano_set <- function(results_path, outdir, prefix, context_label, opt) {

  cat(sprintf("\n========== %s context ==========\n", context_label))
  cat("Reading:", results_path, "\n")

  # ---------------------------------------------------------------------------
  # Data
  # ---------------------------------------------------------------------------

  REQUIRED_COLS <- c("gene_name", "gene_length", "mch_diff",
                     "edger_logFC", "edger_pval", "edger_fdr",
                     "p_bonferroni", "sig_bonferroni", "sig_fdr005")

  df <- read_tsv(results_path, show_col_types = FALSE, progress = FALSE)
  missing <- setdiff(REQUIRED_COLS, names(df))
  if (length(missing) > 0) {
    stop("Results file is missing required columns: ", paste(missing, collapse = ", "))
  }
  if (any(!is.finite(df$edger_pval)) || any(df$edger_pval <= 0)) {
    stop("edger_pval contains non-finite or non-positive values; -log10 transform would be undefined")
  }

  n_tested <- nrow(df)
  sig_col <- if (opt$sig_basis == "bonferroni") "sig_bonferroni" else "sig_fdr005"

  results <- df %>%
    mutate(
      neg_log10_p  = -log10(edger_pval),
      log10_length = log10(gene_length),
      fdr_stringency = -log10(edger_fdr),
      is_sig       = .data[[sig_col]],
      direction = case_when(
        !is_sig            ~ "Not significant",
        edger_logFC > 0    ~ "Higher in mutant",
        TRUE               ~ "Lower in mutant"
      ),
      pt_size  = ifelse(is_sig, 1.9, 0.6),
      pt_alpha = ifelse(is_sig, 0.9, 0.30)
    ) %>%
    arrange(is_sig)

  bonferroni_cut <- opt$alpha / n_tested
  bh_cut <- if (any(results$sig_fdr005)) {
    max(results$edger_pval[results$sig_fdr005])
  } else {
    NA_real_
  }

  # ---------------------------------------------------------------------------
  # Label set
  # ---------------------------------------------------------------------------

  label_genes <- results$gene_name[order(results$edger_pval)][seq_len(min(opt$label_top, n_tested))]

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

  # ---------------------------------------------------------------------------
  # Shared volcano geometry
  # ---------------------------------------------------------------------------

  n_up     <- sum(results$is_sig & results$edger_logFC > 0)
  n_down   <- sum(results$is_sig & results$edger_logFC < 0)

  basis_label <- if (opt$sig_basis == "bonferroni") "Bonferroni" else "BH FDR"

  caption_lines <- c(
    sprintf("%d genes tested. Significance basis: %s p < %g (%d genes). Raw p on the y-axis; horizontal lines mark the corrected cutoffs.",
            n_tested, basis_label, opt$alpha, sum(results$is_sig))
  )
  if (length(gene_list_hits) > 0) {
    caption_lines <- c(caption_lines,
                       sprintf("%s: %d significant members labelled.", opt$set_name, length(gene_list_hits)))
  }
  volcano_caption <- paste(caption_lines, collapse = "\n")

  x_min <- min(results$edger_logFC)

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

  gene_label_layer <- function() {
    geom_text_repel(data = filter(results, nzchar(label)),
                    aes(x = edger_logFC, y = neg_log10_p, label = label),
                    inherit.aes = FALSE,
                    size = 3, colour = "grey15",
                    max.overlaps = 30, box.padding = 0.45, min.segment.length = 0.2,
                    segment.colour = "grey55", segment.size = 0.3, seed = 1)
  }

  x_lab <- sprintf("edgeR log fold change (mutant / control, %s context)", context_label)

  volcano_scaffold <- function(p, title, subtitle) {
    p +
      threshold_layers() +
      gene_label_layer() +
      labs(
        title = title,
        subtitle = subtitle,
        x = x_lab,
        y = expression(-log[10] * (italic(p))),
        caption = volcano_caption
      ) +
      theme_volcano()
  }

  # ---------------------------------------------------------------------------
  # Panel A: Direction
  # ---------------------------------------------------------------------------

  p_direction <- ggplot(results, aes(x = edger_logFC, y = neg_log10_p)) +
    geom_point(aes(colour = direction, size = pt_size, alpha = pt_alpha)) +
    scale_colour_manual(values = COLOR_DIRECTION, name = "Direction",
                        breaks = names(COLOR_DIRECTION)) +
    scale_size_identity() +
    scale_alpha_identity() +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

  p_direction <- volcano_scaffold(
    p_direction,
    sprintf("Gene-body %s differential methylation", context_label),
    sprintf("BAP1-KO vs control  |  %d higher, %d lower in mutant (%s p < %g)",
            n_up, n_down, basis_label, opt$alpha)
  )

  # ---------------------------------------------------------------------------
  # Panel B: FDR stringency
  # ---------------------------------------------------------------------------

  p_fdr <- ggplot(mapping = aes(x = edger_logFC, y = neg_log10_p)) +
    geom_point(data = filter(results, !is_sig), colour = "grey70", size = 0.6, alpha = 0.30) +
    geom_point(data = filter(results, is_sig), aes(colour = fdr_stringency), size = 1.9, alpha = 0.9) +
    scale_colour_gradient(low = RAMP_FDR_STRINGENCY[1], high = RAMP_FDR_STRINGENCY[2],
                          name = expression(-log[10] * "(FDR)"))

  p_fdr <- volcano_scaffold(
    p_fdr,
    sprintf("%s volcano coloured by FDR stringency", context_label),
    sprintf("%d genes at %s < %g; colour depth = distance from significance boundary",
            sum(results$is_sig), basis_label, opt$alpha)
  )

  # ---------------------------------------------------------------------------
  # Panel C: Gene length
  # ---------------------------------------------------------------------------

  median_len_sig <- median(results$gene_length[results$is_sig])
  median_len_all <- median(results$gene_length)

  p_length <- ggplot(mapping = aes(x = edger_logFC, y = neg_log10_p)) +
    geom_point(data = filter(results, !is_sig), colour = "grey70", size = 0.6, alpha = 0.30) +
    geom_point(data = filter(results, is_sig), aes(colour = log10_length), size = 1.9, alpha = 0.9) +
    scale_colour_gradient(low = RAMP_GENELENGTH[1], high = RAMP_GENELENGTH[2],
                          name = expression(log[10] * "(gene length, bp)"))

  p_length <- volcano_scaffold(
    p_length,
    sprintf("%s volcano coloured by gene length", context_label),
    sprintf("Median length %s kb among significant genes vs %s kb across all tested genes",
            format(round(median_len_sig / 1000, 1), nsmall = 1),
            format(round(median_len_all / 1000, 1), nsmall = 1))
  )

  # ---------------------------------------------------------------------------
  # Write
  # ---------------------------------------------------------------------------

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  out <- function(suffix) file.path(outdir, paste0(prefix, "_", suffix))

  cat("\nWriting figures to ", normalizePath(outdir, mustWork = FALSE), "\n", sep = "")
  save_multiformat_ggplot(p_direction, out("direction"), opt$width, opt$height, opt$dpi)
  save_multiformat_ggplot(p_fdr,       out("fdr_stringency"), opt$width, opt$height, opt$dpi)
  save_multiformat_ggplot(p_length,    out("genelength"), opt$width, opt$height, opt$dpi)

  p_panel <- p_direction / (p_fdr | p_length) +
    plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold")))
  save_multiformat_ggplot(p_panel, out("panel"), opt$width * 1.6, opt$height * 1.9, opt$dpi)

  # ---------------------------------------------------------------------------
  # Console summary
  # ---------------------------------------------------------------------------

  n_fdr05 <- sum(results$sig_fdr005)

  cat("\n")
  cat("Context:                 ", context_label, "\n", sep = "")
  cat("Genes tested:            ", n_tested, "\n", sep = "")
  cat("Significance basis:      ", basis_label, " p < ", opt$alpha, "\n", sep = "")
  cat("Significant:             ", sum(results$is_sig),
      "  (", n_up, " higher, ", n_down, " lower in mutant)\n", sep = "")
  cat("FDR < 0.05:              ", n_fdr05, "\n", sep = "")
  cat("Bonferroni raw-p cutoff: ", format(bonferroni_cut, digits = 4), "\n", sep = "")
  cat("BH raw-p cutoff:         ",
      if (is.na(bh_cut)) "none (no gene passes)" else format(bh_cut, digits = 4), "\n", sep = "")
  cat("Largest -log10(p):       ", format(max(results$neg_log10_p), digits = 4), "\n", sep = "")
  cat("Genes labelled:          ", sum(nzchar(results$label)), "\n", sep = "")
  cat("Median gene length (sig):", format(round(median_len_sig / 1000, 1), nsmall = 1), "kb\n", sep = " ")
  cat("Median gene length (all):", format(round(median_len_all / 1000, 1), nsmall = 1), "kb\n", sep = " ")

  invisible(list(
    p_direction = p_direction,
    p_fdr = p_fdr,
    p_length = p_length,
    p_panel = p_panel,
    results = results
  ))
}

# =============================================================================
# MAIN
# =============================================================================

cat("=== mCH Volcano Plots ===\n")

build_volcano_set(opt$results, opt$outdir, opt$prefix, "mCH", opt)

if (!is.null(opt$results_ca)) {
  build_volcano_set(opt$results_ca, opt$outdir_ca, opt$prefix_ca, "mCA", opt)
}

cat("\n=== Complete ===\n")
