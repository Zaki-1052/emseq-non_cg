# scripts/03_mch_differential.R
#
# Differential gene-body mCH testing between BAP1-KO and wildtype cerebellum.
#
# Reads per-sample aggregated gene-body mCH counts (from mch_aggregate_sample.R),
# filters genes, and runs edgeR quasi-likelihood F-test with
# model ~ genotype + sex + lambda_ch_rate and log(total_coverage) as offset.
# Per-sample lambda conversion noise rates enter as a covariate rather than
# being pre-subtracted from counts: all-CH rates (~1-3%) sit near the lambda
# noise floor (~1%), so pre-subtraction zeroes a large fraction of genes and
# creates a systematic bias when lambda correlates with genotype.
#
# Usage:
#   Rscript 03_mch_differential.R \
#       --agg-dir /path/to/results/02_aggregate/aggregated \
#       --out-dir /path/to/results/03_differential \
#       --plot-dir /path/to/results/plots

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(ggrepel)
})

# Source multi-format output utility
script_dir <- if (interactive()) "scripts" else dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "utils", "multi_format_output.R"))

# ---------------------------------------------------------------------------
# Sample metadata
# ---------------------------------------------------------------------------

SAMPLE_META <- data.frame(
  sample_id = c("ctrl_M1", "ctrl_M2", "ctrl_F1", "ctrl_F2",
                "mut_M1", "mut_M2", "mut_F1", "mut_F2"),
  genotype = factor(c(rep("ctrl", 4), rep("mut", 4)),
                    levels = c("ctrl", "mut")),
  sex = factor(c("M", "M", "F", "F", "M", "M", "F", "F")),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

parse_cli_args <- function() {
  option_list <- list(
    make_option("--agg-dir", type = "character", dest = "agg_dir",
                default = NULL,
                help = "Directory with *_genebody_mch.tsv and *_lambda.tsv files"),
    make_option("--out-dir", type = "character", dest = "out_dir",
                default = NULL,
                help = "Output directory for result tables"),
    make_option("--plot-dir", type = "character", dest = "plot_dir",
                default = NULL,
                help = "Output directory for plots"),
    make_option("--min-ch-sites", type = "integer", dest = "min_ch_sites",
                default = 100L,
                help = "Min CH sites per gene in every sample [default %default]"),
    make_option("--min-total-coverage", type = "integer", dest = "min_total_coverage",
                default = 2500L,
                help = "Min total coverage per gene in every sample [default %default]"),
    make_option("--alpha", type = "double", default = 0.05,
                help = "Significance level [default %default]"),
    make_option("--threads", type = "integer", default = getDTthreads(),
                help = "Threads for data.table [default %default]")
  )

  opt <- parse_args(OptionParser(option_list = option_list))
  for (req in c("agg_dir", "out_dir", "plot_dir")) {
    if (is.null(opt[[req]])) stop("Missing required argument: --", gsub("_", "-", req))
  }
  opt
}

# ---------------------------------------------------------------------------
# Load aggregated data
# ---------------------------------------------------------------------------

load_aggregated <- function(agg_dir) {
  sample_ids <- SAMPLE_META$sample_id
  tables <- list()

  for (sid in sample_ids) {
    path <- file.path(agg_dir, paste0(sid, "_genebody_mch.tsv"))
    if (!file.exists(path)) stop("Missing aggregated file: ", path)
    tables[[sid]] <- fread(path)
    cat("  Loaded:", sid, "(", nrow(tables[[sid]]), "genes )\n")
  }

  ref_genes <- tables[[sample_ids[1]]]$gene_name
  for (sid in sample_ids[-1]) {
    if (!identical(tables[[sid]]$gene_name, ref_genes)) {
      stop("Gene order mismatch between ", sample_ids[1], " and ", sid)
    }
  }

  gene_info <- tables[[sample_ids[1]]][, .(gene_name, gene_id, chr, start, end,
                                           strand, gene_type, gene_length)]
  n_genes <- nrow(gene_info)
  n_samples <- length(sample_ids)

  M <- matrix(0, nrow = n_genes, ncol = n_samples,
              dimnames = list(NULL, sample_ids))
  N <- matrix(0, nrow = n_genes, ncol = n_samples,
              dimnames = list(NULL, sample_ids))
  sites <- matrix(0L, nrow = n_genes, ncol = n_samples,
                  dimnames = list(NULL, sample_ids))

  for (i in seq_along(sample_ids)) {
    sid <- sample_ids[i]
    M[, i] <- tables[[sid]]$meth_reads
    N[, i] <- tables[[sid]]$total_coverage
    sites[, i] <- tables[[sid]]$n_ch_sites
  }

  list(gene_info = gene_info, M = M, N = N, sites = sites)
}

load_lambda <- function(agg_dir) {
  sample_ids <- SAMPLE_META$sample_id
  lambda <- data.table()

  for (sid in sample_ids) {
    path <- file.path(agg_dir, paste0(sid, "_lambda.tsv"))
    if (!file.exists(path)) stop("Missing lambda file: ", path)
    lambda <- rbind(lambda, fread(path))
  }

  lambda
}

# ---------------------------------------------------------------------------
# Gene filtering
# ---------------------------------------------------------------------------

filter_genes <- function(sites, N, min_sites, min_cov) {
  pass_sites <- apply(sites, 1, function(x) all(x >= min_sites))
  pass_cov <- apply(N, 1, function(x) all(x >= min_cov))
  pass <- pass_sites & pass_cov

  cat("\n--- Gene filtering ---\n")
  cat("Total genes:", nrow(sites), "\n")
  cat("Failed n_ch_sites >=", min_sites, "in all samples:",
      sum(!pass_sites), "\n")
  cat("Failed total_coverage >=", min_cov, "in all samples:",
      sum(pass_sites & !pass_cov), "\n")
  cat("Passed both filters:", sum(pass), "\n")

  pass
}

# ---------------------------------------------------------------------------
# QQ plot
# ---------------------------------------------------------------------------

make_qq_plot <- function(pvals, plot_dir) {
  n <- length(pvals)
  expected <- -log10(ppoints(n))
  observed <- -log10(sort(pvals))
  lambda_gc <- median(qchisq(1 - pvals, df = 1)) / qchisq(0.5, df = 1)

  qq_df <- data.frame(expected = expected, observed = observed)

  p <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(size = 0.3, alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.5) +
    labs(x = expression(Expected ~ -log[10](p)),
         y = expression(Observed ~ -log[10](p)),
         title = paste0("QQ Plot — mCH Differential"),
         subtitle = bquote(lambda[GC] == .(sprintf("%.3f", lambda_gc)))) +
    theme_bw(base_size = 12) +
    coord_equal()

  save_multiformat_ggplot(p, file.path(plot_dir, "03_mch_qq_plot"), width = 6, height = 6)
  cat("QQ plot saved (lambda_GC =", sprintf("%.3f", lambda_gc), ")\n")

  lambda_gc
}

# ---------------------------------------------------------------------------
# BCV plot
# ---------------------------------------------------------------------------

make_bcv_plot <- function(y, plot_dir) {
  save_multiformat_base(
    quote(plotBCV(y, main = "BCV Plot — mCH Gene Bodies")),
    file.path(plot_dir, "03_mch_bcv_plot"),
    width = 7, height = 6
  )

  cat("BCV plot saved\n")
  cat("  Common dispersion:", sprintf("%.6f", y$common.dispersion), "\n")
  cat("  Common BCV:", sprintf("%.4f", sqrt(y$common.dispersion)), "\n")
  cat("  Tagwise dispersion range:",
      sprintf("%.6f - %.6f",
              min(y$tagwise.dispersion), max(y$tagwise.dispersion)), "\n")
}

# ---------------------------------------------------------------------------
# Summary report
# ---------------------------------------------------------------------------

write_summary <- function(out_dir, n_total, n_tested, lambda_df, y,
                          top, pass, alpha, lambda_gc, gene_info,
                          sample_meta) {
  pvals <- top$table$edger_pval
  n <- n_tested
  bonf <- pmin(pvals * n, 1)
  fdr <- top$table$edger_fdr
  logfc <- top$table$edger_logFC

  sink(file.path(out_dir, "mch_analysis_summary.txt"))
  cat("=== mCH Gene-Body Differential Analysis Summary ===\n\n")

  cat("--- Model ---\n")
  cat("  Formula: ~ genotype + sex + lambda_ch_rate\n")
  cat("  Lambda enters as covariate (no pre-subtraction of counts)\n")
  cat("  Per-sample rates are raw (uncorrected); genotype effect is\n")
  cat("  lambda-adjusted via the model.\n\n")

  cat("--- Gene counts ---\n")
  cat("Total genes:", n_total, "\n")
  cat("Testable genes:", n_tested, "\n\n")

  cat("--- Lambda CH rates per sample ---\n")
  for (i in seq_len(nrow(lambda_df))) {
    cat(sprintf("  %s: %.4f (%.2f%%)\n",
                lambda_df$sample_id[i],
                lambda_df$lambda_ch_rate[i],
                lambda_df$lambda_ch_rate[i] * 100))
  }
  cat("\n")

  cat("--- Covariate correlations ---\n")
  numeric_design <- data.frame(
    genotype_mut = as.integer(sample_meta$genotype == "mut"),
    sex_F = as.integer(sample_meta$sex == "F"),
    lambda_ch_rate = sample_meta$lambda_ch_rate
  )
  cormat <- cor(numeric_design)
  for (rn in rownames(cormat)) {
    cat(sprintf("  %-16s %s\n", rn,
                paste(sprintf("%+.3f", cormat[rn, ]), collapse = "  ")))
  }
  genotype_lambda_r <- cormat["genotype_mut", "lambda_ch_rate"]
  cat(sprintf("\n  genotype vs lambda_ch_rate r = %+.3f\n", genotype_lambda_r))
  if (abs(genotype_lambda_r) > 0.7) {
    cat("  WARNING: high genotype-lambda correlation; genotype SE inflated\n")
  }
  cat("\n")

  cat("--- edgeR dispersion estimates ---\n")
  cat(sprintf("  Common dispersion: %.6f (BCV: %.4f)\n",
              y$common.dispersion, sqrt(y$common.dispersion)))
  cat(sprintf("  Tagwise range: %.6f - %.6f\n",
              min(y$tagwise.dispersion), max(y$tagwise.dispersion)))
  cat(sprintf("  Tagwise median: %.6f\n", median(y$tagwise.dispersion)))
  prior_df <- y$prior.df
  cat(sprintf("  Prior df (NB): median %.1f, range %.1f - %.1f\n",
              median(prior_df), min(prior_df), max(prior_df)))
  cat("\n")

  n_bonf <- sum(bonf < alpha)
  n_fdr05 <- sum(fdr < 0.05)
  n_fdr10 <- sum(fdr < 0.10)

  cat("--- Significant genes ---\n")
  cat(sprintf("  Bonferroni < %.2f: %d\n", alpha, n_bonf))
  cat(sprintf("  FDR < 0.05: %d\n", n_fdr05))
  cat(sprintf("  FDR < 0.10: %d\n", n_fdr10))
  cat("\n")

  if (n_fdr05 > 0) {
    sig_up <- sum(fdr < 0.05 & logfc > 0)
    sig_down <- sum(fdr < 0.05 & logfc < 0)
    cat("--- Direction (FDR < 0.05) ---\n")
    cat(sprintf("  mCH higher in mutant: %d\n", sig_up))
    cat(sprintf("  mCH lower in mutant: %d\n", sig_down))
    cat("\n")
  }

  cat("--- Effect sizes (all testable genes) ---\n")
  cat(sprintf("  logFC range: %.4f to %.4f\n", min(logfc), max(logfc)))
  cat(sprintf("  logFC median: %.4f\n", median(logfc)))
  cat(sprintf("  logFC IQR: %.4f to %.4f\n",
              quantile(logfc, 0.25), quantile(logfc, 0.75)))
  cat("\n")

  cat("--- Genomic inflation ---\n")
  cat(sprintf("  lambda_GC: %.3f\n", lambda_gc))
  cat("\n")

  cat("--- Top 20 genes by p-value ---\n")
  top20 <- head(top$table, 20)
  for (i in seq_len(nrow(top20))) {
    cat(sprintf("  %2d. %-15s logFC=%+.4f  F=%.1f  p=%.2e  FDR=%.2e\n",
                i, top20$gene_name[i],
                top20$edger_logFC[i], top20$edger_F[i],
                top20$edger_pval[i], top20$edger_fdr[i]))
  }

  sink()
  cat("Summary written to:", file.path(out_dir, "mch_analysis_summary.txt"), "\n")
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  opt <- parse_cli_args()
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(opt$plot_dir, recursive = TRUE, showWarnings = FALSE)

  cat("=== mCH Gene-Body Differential Testing ===\n\n")
  setDTthreads(opt$threads)
  cat("data.table threads:", getDTthreads(), "\n")

  # --- Load data ---
  cat("Loading aggregated data:\n")
  dat <- load_aggregated(opt$agg_dir)
  cat("\nLoading lambda rates:\n")
  lambda_df <- load_lambda(opt$agg_dir)
  cat("  ", paste(sprintf("%s=%.4f", lambda_df$sample_id, lambda_df$lambda_ch_rate),
                  collapse = ", "), "\n")

  # --- Merge lambda into sample metadata ---
  lambda_rates <- setNames(lambda_df$lambda_ch_rate, lambda_df$sample_id)
  SAMPLE_META$lambda_ch_rate <- lambda_rates[SAMPLE_META$sample_id]
  cat("\nLambda rates merged into sample metadata\n")
  cat("  ctrl mean:", sprintf("%.4f", mean(SAMPLE_META$lambda_ch_rate[SAMPLE_META$genotype == "ctrl"])), "\n")
  cat("  mut mean: ", sprintf("%.4f", mean(SAMPLE_META$lambda_ch_rate[SAMPLE_META$genotype == "mut"])), "\n")

  # --- Gene filtering ---
  pass <- filter_genes(dat$sites, dat$N, opt$min_ch_sites, opt$min_total_coverage)
  n_total <- nrow(dat$gene_info)
  n_tested <- sum(pass)

  # --- Compute per-sample raw rates (no lambda subtraction) ---
  ctrl_cols <- grep("^ctrl_", SAMPLE_META$sample_id, value = TRUE)
  mut_cols <- grep("^mut_", SAMPLE_META$sample_id, value = TRUE)
  rate_matrix <- dat$M[pass, ] / dat$N[pass, ]

  gene_annot <- dat$gene_info[pass, ]
  gene_annot$mch_ctrl <- rowMeans(rate_matrix[, ctrl_cols])
  gene_annot$mch_mut <- rowMeans(rate_matrix[, mut_cols])
  gene_annot$mch_diff <- gene_annot$mch_mut - gene_annot$mch_ctrl
  for (sid in SAMPLE_META$sample_id) {
    gene_annot[[paste0(sid, "_mch")]] <- rate_matrix[, sid]
  }
  gene_annot$n_ch_sites_mean <- rowMeans(dat$sites[pass, ])
  gene_annot$total_coverage_mean <- rowMeans(dat$N[pass, ])

  # --- edgeR (raw counts, lambda as covariate) ---
  cat("\n--- edgeR quasi-likelihood F-test ---\n")
  y <- DGEList(counts = dat$M[pass, ])
  y$offset <- log(dat$N[pass, ])
  y$genes <- gene_annot

  design <- model.matrix(~ genotype + sex + lambda_ch_rate, data = SAMPLE_META)
  cat("Design matrix columns:", paste(colnames(design), collapse = ", "), "\n")
  cat("Residual df:", nrow(design) - ncol(design), "\n")
  cat("Estimating dispersions (robust=TRUE)...\n")
  y <- estimateDisp(y, design, robust = TRUE)

  cat("Fitting QL model...\n")
  fit <- glmQLFit(y, design)

  cat("Testing genotype effect...\n")
  res <- glmQLFTest(fit, coef = "genotypemut")
  top <- topTags(res, n = Inf, sort.by = "PValue")

  # --- Multiple testing ---
  pvals <- top$table$PValue
  bonf <- pmin(pvals * n_tested, 1)
  top$table$p_bonferroni <- bonf
  top$table$sig_bonferroni <- bonf < opt$alpha
  top$table$sig_fdr005 <- top$table$FDR < 0.05
  top$table$sig_fdr010 <- top$table$FDR < 0.10

  cat(sprintf("\nSignificant genes:\n"))
  cat(sprintf("  Bonferroni < %.2f: %d\n", opt$alpha, sum(bonf < opt$alpha)))
  cat(sprintf("  FDR < 0.05: %d\n", sum(top$table$FDR < 0.05)))
  cat(sprintf("  FDR < 0.10: %d\n", sum(top$table$FDR < 0.10)))

  # --- Rename edgeR columns for output ---
  names(top$table)[names(top$table) == "logFC"] <- "edger_logFC"
  names(top$table)[names(top$table) == "logCPM"] <- "edger_logCPM"
  names(top$table)[names(top$table) == "F"] <- "edger_F"
  names(top$table)[names(top$table) == "PValue"] <- "edger_pval"
  names(top$table)[names(top$table) == "FDR"] <- "edger_fdr"

  # --- Write results ---
  results_file <- file.path(opt$out_dir, "mch_differential_results.tsv")
  fwrite(top$table, results_file, sep = "\t")
  cat("\nResults written to:", results_file, "(", nrow(top$table), "genes )\n")

  # --- Sample matrix ---
  sample_matrix <- as.data.table(rate_matrix)
  sample_matrix <- cbind(dat$gene_info[pass, .(gene_name, gene_id)], sample_matrix)
  matrix_file <- file.path(opt$out_dir, "mch_sample_matrix.tsv")
  fwrite(sample_matrix, matrix_file, sep = "\t")
  cat("Sample matrix written to:", matrix_file, "\n")

  # --- Plots ---
  cat("\nGenerating plots...\n")
  lambda_gc <- make_qq_plot(top$table$edger_pval, opt$plot_dir)
  make_bcv_plot(y, opt$plot_dir)

  # --- Summary report ---
  write_summary(opt$out_dir, n_total, n_tested, lambda_df, y,
                top, pass, opt$alpha, lambda_gc, dat$gene_info,
                SAMPLE_META)

  cat("\n=== Complete ===\n")
}

main()
