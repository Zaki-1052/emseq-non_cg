# scripts/sections/80_cross_modality/80_01_cross_modality.R
#
# Section 80_01: CG and non-CG methylation cross-modality comparison.
#
# What this tests
#   The EM-seq pipeline measures non-CG methylation (mCH). The companion
#   Biomodal evoC pipeline measures CG 5mC and CG 5hmC. If H2AK119ub
#   instructs DNMT3A and TET activity broadly, CG and non-CG methylation
#   changes should correlate. This section puts both data sets on comparable
#   log-fold-change scales and tests the correlation directly.
#
# Why the prior comparison was inadequate
#   Step 04d compared CG mod_difference (arithmetic: mutant_fraction minus
#   control_fraction) against mCH edger_logFC (model-adjusted log2 fold
#   change). Those are different scales. Here we use log2(mod_fold_change)
#   for CG, which gives a log2 fold change comparable in kind to edger_logFC.
#
# Analyses
#   1. Three Spearman correlations (mCH vs CG 5mC, mCH vs CG 5hmC,
#      CG 5mC vs CG 5hmC) on log-scale fold changes.
#   2. Scatter plots with Spearman rho, coloured by mCH significance.
#   3. Rho comparison with 95% CIs (Fisher z) and Steiger's test.
#   4. Direction agreement for genes significant in both modalities.
#
# Reads
#   mch_results  (shared config)
#   data/biomodal/cg_mc_dmr_genes.bed
#   data/biomodal/cg_hmc_dmr_genes.bed
#
# Writes
#   Figures and tables into OUTPUT_PATHS$cross_modality.
#   Fisher test shards into HANDOFF_PATHS$fisher_registry.

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)

SECTION_ID <- "80_01"

# =============================================================================
# OPTIONS
# =============================================================================

parse_options <- function() {
  option_list <- list(
    make_option("--output-dir", type = "character", default = OUTPUT_PATHS$cross_modality,
                dest = "output_dir",
                help = "Directory for figures and tables [default: %default]"),
    make_option("--fdr-threshold", type = "double", default = Q_THRESHOLD,
                dest = "fdr_threshold",
                help = "FDR cutoff for significance calls [default: %default]")
  )
  parse_args(OptionParser(option_list = option_list))
}

# =============================================================================
# DATA LOADING
# =============================================================================

#' Read a Biomodal DMR BED file and compute log2 fold change.
#'
#' The Biomodal modality XPLR pipeline reports mod_fold_change as the
#' arithmetic ratio mean_mod_group_2 / mean_mod_group_1. log2 of that ratio
#' gives a log2 fold change comparable to edgeR logFC in kind.
#'
#' Genes with mod_fold_change <= 0 or NA are dropped (these occur when
#' one group's mean is exactly zero, which is rare for CG methylation).
load_biomodal_dmr <- function(filepath, modality_name) {
  if (!file.exists(filepath)) {
    stop("Biomodal DMR file not found: ", filepath)
  }
  dmr <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                    quote = "", comment.char = "")

  required <- c("Name", "mod_fold_change", "mod_difference", "dmr_qvalue",
                 "mean_mod_group_1", "mean_mod_group_2")
  missing <- setdiff(required, colnames(dmr))
  if (length(missing) > 0) {
    stop("Biomodal DMR file is missing columns: ", paste(missing, collapse = ", "))
  }

  dmr <- dmr[!is.na(dmr$mod_fold_change) & dmr$mod_fold_change > 0, , drop = FALSE]
  dmr$log2fc <- log2(dmr$mod_fold_change)

  cat(sprintf("  %s: %s genes, %s significant at FDR < %.3f\n",
              modality_name,
              format(nrow(dmr), big.mark = ","),
              format(sum(dmr$dmr_qvalue < Q_THRESHOLD), big.mark = ","),
              Q_THRESHOLD))
  dmr
}

# =============================================================================
# CORRELATION HELPERS
# =============================================================================

#' Fisher z-transform 95% CI for a Spearman rho.
rho_ci <- function(rho, n, level = 0.95) {
  z <- atanh(rho)
  se <- 1 / sqrt(n - 3)
  crit <- qnorm(1 - (1 - level) / 2)
  lo <- tanh(z - crit * se)
  hi <- tanh(z + crit * se)
  c(lo = lo, hi = hi)
}

#' Steiger's test for comparing two dependent correlations that share a
#' variable. Tests r(XY) vs r(XZ) where r(YZ) is also observed.
steiger_test <- function(rho_xy, rho_xz, rho_yz, n) {
  r_det <- 1 - rho_xy^2 - rho_xz^2 - rho_yz^2 + 2 * rho_xy * rho_xz * rho_yz
  r_bar <- (rho_xy + rho_xz) / 2
  denom <- sqrt((2 * (1 - rho_yz)) *
                  ((1 - r_bar^2)^2) /
                  ((1 - r_det) * (n - 3) / (n - 1)))
  if (denom == 0) return(list(z = NA, p = NA))
  z_stat <- (atanh(rho_xy) - atanh(rho_xz)) * sqrt((n - 3) / 2) /
            sqrt(1 - rho_yz)
  p_val <- 2 * pnorm(-abs(z_stat))
  list(z = z_stat, p = p_val)
}

# =============================================================================
# PLOTS
# =============================================================================

#' Cross-modality scatter with Spearman rho, coloured by mCH significance.
cross_scatter <- function(df, x_col, y_col, x_lab, y_lab, rho, rho_ci_lo,
                          rho_ci_hi, n, key_genes, mch_sig_col = "mch_sig") {
  df$sig_label <- ifelse(df[[mch_sig_col]], "mCH Significant", "Not Significant")
  df$sig_label <- factor(df$sig_label, levels = c("mCH Significant", "Not Significant"))

  subtitle <- sprintf("Spearman rho = %.3f [%.3f, %.3f], n = %s",
                       rho, rho_ci_lo, rho_ci_hi,
                       format(n, big.mark = ","))

  key_df <- df[df$gene_name %in% key_genes, , drop = FALSE]

  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(aes(colour = sig_label), size = 0.6, alpha = 0.4) +
    scale_colour_manual(values = c("mCH Significant" = COLORS$direction[["Hypomethylated"]],
                                   "Not Significant" = "grey70"),
                        name = NULL) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    labs(x = x_lab, y = y_lab, subtitle = subtitle) +
    theme_emseq() +
    theme(legend.position = "bottom")

  if (nrow(key_df) > 0) {
    p <- p +
      geom_point(data = key_df, colour = "black", size = 1.5, shape = 21,
                 fill = COLORS$direction[["Hypermethylated"]]) +
      geom_text_repel(data = key_df, aes(label = gene_name),
                      size = 2.8, max.overlaps = 20, seed = 42,
                      min.segment.length = 0)
  }
  p
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opts <- parse_options()
  out_dir <- opts$output_dir
  fdr <- opts$fdr_threshold
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("=== Section 80_01: Cross-modality CG/hmC/mCH comparison ===\n\n")

  # --- Load Biomodal data ---
  cat("Loading Biomodal CG DMR data...\n")
  mc_path  <- file.path(CODE_DIR, "data", "biomodal", "cg_mc_dmr_genes.bed")
  hmc_path <- file.path(CODE_DIR, "data", "biomodal", "cg_hmc_dmr_genes.bed")
  stopifnot(file.exists(mc_path))
  stopifnot(file.exists(hmc_path))

  mc_dmr  <- load_biomodal_dmr(mc_path, "CG 5mC")
  hmc_dmr <- load_biomodal_dmr(hmc_path, "CG 5hmC")

  # --- Merge with mCH data ---
  cat("\nMerging with mCH differential results...\n")

  mc_sub <- data.frame(
    gene_name = mc_dmr$Name,
    cg_mc_log2fc = mc_dmr$log2fc,
    cg_mc_diff = mc_dmr$mod_difference,
    cg_mc_qvalue = mc_dmr$dmr_qvalue,
    cg_mc_ctrl = mc_dmr$mean_mod_group_1,
    cg_mc_mut = mc_dmr$mean_mod_group_2,
    stringsAsFactors = FALSE
  )

  hmc_sub <- data.frame(
    gene_name = hmc_dmr$Name,
    cg_hmc_log2fc = hmc_dmr$log2fc,
    cg_hmc_diff = hmc_dmr$mod_difference,
    cg_hmc_qvalue = hmc_dmr$dmr_qvalue,
    cg_hmc_ctrl = hmc_dmr$mean_mod_group_1,
    cg_hmc_mut = hmc_dmr$mean_mod_group_2,
    stringsAsFactors = FALSE
  )

  mch_sub <- data.frame(
    gene_name = mch_results$gene_name,
    mch_logfc = mch_results$edger_logFC,
    mch_fdr = mch_results$edger_fdr,
    mch_sig = mch_results$mch_sig,
    mch_direction = mch_results$mch_direction,
    mch_diff = mch_results$mch_diff,
    stringsAsFactors = FALSE
  )

  merged <- merge(mch_sub, mc_sub, by = "gene_name", all = FALSE)
  merged <- merge(merged, hmc_sub, by = "gene_name", all = FALSE)

  n_complete <- sum(complete.cases(merged[, c("mch_logfc", "cg_mc_log2fc", "cg_hmc_log2fc")]))
  merged <- merged[complete.cases(merged[, c("mch_logfc", "cg_mc_log2fc", "cg_hmc_log2fc")]), ]

  cat(sprintf("  Genes with all three modalities: %s\n",
              format(nrow(merged), big.mark = ",")))
  cat(sprintf("  mCH significant: %s\n",
              format(sum(merged$mch_sig), big.mark = ",")))

  # --- Spearman correlations ---
  cat("\nComputing Spearman correlations on log2 fold changes...\n")

  cor_mch_mc <- cor.test(merged$mch_logfc, merged$cg_mc_log2fc, method = "spearman",
                         exact = FALSE)
  cor_mch_hmc <- cor.test(merged$mch_logfc, merged$cg_hmc_log2fc, method = "spearman",
                          exact = FALSE)
  cor_mc_hmc <- cor.test(merged$cg_mc_log2fc, merged$cg_hmc_log2fc, method = "spearman",
                         exact = FALSE)

  n <- nrow(merged)
  rho_mch_mc  <- cor_mch_mc$estimate
  rho_mch_hmc <- cor_mch_hmc$estimate
  rho_mc_hmc  <- cor_mc_hmc$estimate

  ci_mch_mc  <- rho_ci(rho_mch_mc, n)
  ci_mch_hmc <- rho_ci(rho_mch_hmc, n)
  ci_mc_hmc  <- rho_ci(rho_mc_hmc, n)

  cat(sprintf("  mCH logFC vs CG 5mC log2FC:  rho = %.4f [%.4f, %.4f], p = %.2e\n",
              rho_mch_mc, ci_mch_mc["lo"], ci_mch_mc["hi"], cor_mch_mc$p.value))
  cat(sprintf("  mCH logFC vs CG 5hmC log2FC: rho = %.4f [%.4f, %.4f], p = %.2e\n",
              rho_mch_hmc, ci_mch_hmc["lo"], ci_mch_hmc["hi"], cor_mch_hmc$p.value))
  cat(sprintf("  CG 5mC vs CG 5hmC log2FC:    rho = %.4f [%.4f, %.4f], p = %.2e\n",
              rho_mc_hmc, ci_mc_hmc["lo"], ci_mc_hmc["hi"], cor_mc_hmc$p.value))

  # Steiger's test: is rho(mCH, 5mC) different from rho(mCH, 5hmC)?
  steiger <- steiger_test(rho_mch_mc, rho_mch_hmc, rho_mc_hmc, n)
  cat(sprintf("  Steiger's test (mCH~5mC vs mCH~5hmC): z = %.3f, p = %.2e\n",
              steiger$z, steiger$p))

  # --- Correlation summary table ---
  cor_table <- data.frame(
    comparison = c("mCH_logFC vs CG_5mC_log2FC",
                   "mCH_logFC vs CG_5hmC_log2FC",
                   "CG_5mC_log2FC vs CG_5hmC_log2FC"),
    n = n,
    rho = c(rho_mch_mc, rho_mch_hmc, rho_mc_hmc),
    ci_lo = c(ci_mch_mc["lo"], ci_mch_hmc["lo"], ci_mc_hmc["lo"]),
    ci_hi = c(ci_mch_mc["hi"], ci_mch_hmc["hi"], ci_mc_hmc["hi"]),
    p_value = c(cor_mch_mc$p.value, cor_mch_hmc$p.value, cor_mc_hmc$p.value),
    steiger_z = c(steiger$z, steiger$z, NA),
    steiger_p = c(steiger$p, steiger$p, NA),
    stringsAsFactors = FALSE
  )
  write_section_table(cor_table, file.path(out_dir, "80_01_correlation_comparison.tsv"))

  # --- Scatter plots ---
  cat("\nGenerating scatter plots...\n")

  p_mch_mc <- cross_scatter(
    merged, "cg_mc_log2fc", "mch_logfc",
    x_lab = expression("CG 5mC " * log[2] * "(mut/ctrl)"),
    y_lab = expression("mCH edgeR " * log[2] * "FC"),
    rho = rho_mch_mc, rho_ci_lo = ci_mch_mc["lo"], rho_ci_hi = ci_mch_mc["hi"],
    n = n, key_genes = KEY_GENES
  ) + ggtitle("mCH vs CG 5-methylcytosine")

  p_mch_hmc <- cross_scatter(
    merged, "cg_hmc_log2fc", "mch_logfc",
    x_lab = expression("CG 5hmC " * log[2] * "(mut/ctrl)"),
    y_lab = expression("mCH edgeR " * log[2] * "FC"),
    rho = rho_mch_hmc, rho_ci_lo = ci_mch_hmc["lo"], rho_ci_hi = ci_mch_hmc["hi"],
    n = n, key_genes = KEY_GENES
  ) + ggtitle("mCH vs CG 5-hydroxymethylcytosine")

  p_mc_hmc <- cross_scatter(
    merged, "cg_mc_log2fc", "cg_hmc_log2fc",
    x_lab = expression("CG 5mC " * log[2] * "(mut/ctrl)"),
    y_lab = expression("CG 5hmC " * log[2] * "(mut/ctrl)"),
    rho = rho_mc_hmc, rho_ci_lo = ci_mc_hmc["lo"], rho_ci_hi = ci_mc_hmc["hi"],
    n = n, key_genes = KEY_GENES
  ) + ggtitle("CG 5mC vs CG 5hmC (reciprocal pattern)")

  save_multiformat_ggplot(p_mch_mc, file.path(out_dir, "80_01a_mch_vs_cg_mc"),
                          width = 8, height = 7)
  save_multiformat_ggplot(p_mch_hmc, file.path(out_dir, "80_01b_mch_vs_cg_hmc"),
                          width = 8, height = 7)
  save_multiformat_ggplot(p_mc_hmc, file.path(out_dir, "80_01c_cg_mc_vs_hmc"),
                          width = 8, height = 7)

  # Composite
  p_composite <- (p_mch_mc | p_mch_hmc | p_mc_hmc) +
    plot_annotation(title = "Cross-modality correlation on log2 fold change scales",
                    theme = theme(plot.title = element_text(face = "bold", size = 14)))
  save_multiformat_ggplot(p_composite, file.path(out_dir, "80_01d_composite"),
                          width = 22, height = 7)

  # --- Direction agreement ---
  cat("\nDirection agreement analysis...\n")

  # mCH vs CG 5mC
  both_mc <- merged[merged$mch_sig & merged$cg_mc_qvalue < fdr, , drop = FALSE]
  if (nrow(both_mc) > 0) {
    both_mc$mch_up <- both_mc$mch_logfc > 0
    both_mc$mc_up  <- both_mc$cg_mc_log2fc > 0
    agree_mc <- sum(both_mc$mch_up == both_mc$mc_up)
    disagree_mc <- sum(both_mc$mch_up != both_mc$mc_up)

    ct_mc <- table(
      mch_direction = ifelse(both_mc$mch_up, "mCH_up", "mCH_down"),
      mc_direction  = ifelse(both_mc$mc_up, "5mC_up", "5mC_down")
    )
    fisher_mc <- fisher.test(ct_mc)

    cat(sprintf("  mCH vs CG 5mC: %s genes sig in both. Agreement: %s (%.1f%%)\n",
                format(nrow(both_mc), big.mark = ","),
                format(agree_mc, big.mark = ","),
                100 * agree_mc / nrow(both_mc)))
    cat(sprintf("    Fisher OR = %.2f, p = %.2e\n", fisher_mc$estimate, fisher_mc$p.value))

    gene_df_mc <- data.frame(
      gene_name = both_mc$gene_name,
      chr = mch_results$chr[match(both_mc$gene_name, mch_results$gene_name)],
      mch_hyper = both_mc$mch_up,
      cg_mc_up = both_mc$mc_up,
      stringsAsFactors = FALSE
    )
    register_fisher_test(SECTION_ID, "mch_direction_vs_cg_mc_direction",
                         "mCH direction vs CG 5mC direction (both significant)",
                         gene_df_mc, "mch_hyper", "cg_mc_up", out_dir)
  }

  # mCH vs CG 5hmC
  both_hmc <- merged[merged$mch_sig & merged$cg_hmc_qvalue < fdr, , drop = FALSE]
  if (nrow(both_hmc) > 0) {
    both_hmc$mch_up <- both_hmc$mch_logfc > 0
    both_hmc$hmc_up <- both_hmc$cg_hmc_log2fc > 0
    agree_hmc <- sum(both_hmc$mch_up == both_hmc$hmc_up)
    disagree_hmc <- sum(both_hmc$mch_up != both_hmc$hmc_up)

    ct_hmc <- table(
      mch_direction = ifelse(both_hmc$mch_up, "mCH_up", "mCH_down"),
      hmc_direction = ifelse(both_hmc$hmc_up, "5hmC_up", "5hmC_down")
    )
    fisher_hmc <- fisher.test(ct_hmc)

    cat(sprintf("  mCH vs CG 5hmC: %s genes sig in both. Agreement: %s (%.1f%%)\n",
                format(nrow(both_hmc), big.mark = ","),
                format(agree_hmc, big.mark = ","),
                100 * agree_hmc / nrow(both_hmc)))
    cat(sprintf("    Fisher OR = %.2f, p = %.2e\n", fisher_hmc$estimate, fisher_hmc$p.value))

    gene_df_hmc <- data.frame(
      gene_name = both_hmc$gene_name,
      chr = mch_results$chr[match(both_hmc$gene_name, mch_results$gene_name)],
      mch_hyper = both_hmc$mch_up,
      cg_hmc_up = both_hmc$hmc_up,
      stringsAsFactors = FALSE
    )
    register_fisher_test(SECTION_ID, "mch_direction_vs_cg_hmc_direction",
                         "mCH direction vs CG 5hmC direction (both significant)",
                         gene_df_hmc, "mch_hyper", "cg_hmc_up", out_dir)
  }

  # Direction agreement summary table
  dir_summary <- data.frame(
    comparison = c("mCH vs CG 5mC", "mCH vs CG 5hmC"),
    n_both_sig = c(nrow(both_mc), nrow(both_hmc)),
    n_agree = c(agree_mc, agree_hmc),
    n_disagree = c(disagree_mc, disagree_hmc),
    pct_agreement = c(100 * agree_mc / nrow(both_mc),
                      100 * agree_hmc / nrow(both_hmc)),
    fisher_or = c(fisher_mc$estimate, fisher_hmc$estimate),
    fisher_p = c(fisher_mc$p.value, fisher_hmc$p.value),
    stringsAsFactors = FALSE
  )
  write_section_table(dir_summary, file.path(out_dir, "80_01_direction_agreement.tsv"))

  # --- Full merged table ---
  write_section_table(merged, file.path(out_dir, "80_01_cross_modality_merged.tsv"))

  cat(sprintf("\nSection 80_01 complete.\n"))
}

main()
