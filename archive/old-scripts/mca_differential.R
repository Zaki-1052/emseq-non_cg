# scripts/mca_differential.R
#
# Differential gene-body mCA testing between BAP1-KO and wildtype cerebellum.
#
# Each gene body is one locus: methylated and total read counts are summed over
# every CA site in the gene. DSS fits a beta-binomial model per locus and shrinks
# the dispersion estimate across loci, which is what makes n=4 per group testable.
# The per-sample lambda conversion noise floor enters as a covariate.
#
# Leave-one-sample-out refits show whether a result depends on a single library.
#
# Usage:
#   Rscript mca_differential.R \
#       --genebody-dir /path/to/genebody_mca \
#       --spike-in /path/to/spike_in_rates.tsv \
#       --output-dir /path/to/mca_results

suppressPackageStartupMessages({
  library(optparse)
  library(DSS)
  library(bsseq)
  library(GenomicRanges)
})

GENE_METADATA_COLUMNS <- c("gene_name", "gene_id", "chr", "start", "end",
                           "strand", "gene_length", "gene_type")

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

parse_cli_args <- function() {
  option_list <- list(
    make_option("--genebody-dir", type = "character", dest = "genebody_dir",
                default = NULL,
                help = "Directory of *.genebody_mca.tsv files from Step 3"),
    make_option("--spike-in", type = "character", dest = "spike_in",
                default = NULL,
                help = "spike_in_rates.tsv from Step 2"),
    make_option("--output-dir", type = "character", dest = "output_dir",
                default = NULL,
                help = "Directory for result tables and the summary report"),
    make_option("--min-ca-sites", type = "integer", dest = "min_ca_sites",
                default = 100L,
                help = "Minimum CA sites per gene in every sample [default %default]"),
    make_option("--min-total-coverage", type = "integer", dest = "min_total_coverage",
                default = 2500L,
                help = "Minimum total coverage per gene in every sample [default %default]"),
    make_option("--alpha", type = "double", default = 0.05,
                help = "Significance level for Bonferroni and FDR [default %default]"),
    make_option("--neuronal-gene-set", type = "character", dest = "neuronal_gene_set",
                default = NULL,
                help = "Optional one-column TSV of neuronal gene symbols for enrichment"),
    make_option("--skip-loo", action = "store_true", dest = "skip_loo",
                default = FALSE,
                help = "Skip leave-one-sample-out influence refits")
  )

  opt <- parse_args(OptionParser(option_list = option_list))

  for (required in c("genebody_dir", "spike_in", "output_dir")) {
    if (is.null(opt[[required]])) {
      stop("Missing required argument: --", gsub("_", "-", required))
    }
  }
  opt
}

# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

load_genebody_tables <- function(genebody_dir) {
  paths <- sort(list.files(genebody_dir, pattern = "\\.genebody_mca\\.tsv$",
                           full.names = TRUE))
  if (length(paths) == 0) {
    stop("No *.genebody_mca.tsv files in ", genebody_dir)
  }

  tables <- lapply(paths, function(path) {
    tab <- read.delim(path, stringsAsFactors = FALSE)
    if (length(unique(tab$sample_id)) != 1L) {
      stop(path, " contains more than one sample_id")
    }
    tab
  })
  names(tables) <- vapply(tables, function(tab) tab$sample_id[1], character(1))

  gene_keys <- lapply(tables, function(tab) tab$gene_name)
  if (!all(vapply(gene_keys, identical, logical(1), gene_keys[[1]]))) {
    stop("Per-sample gene-body tables do not share an identical gene ordering")
  }

  message(sprintf("Loaded %d samples x %d genes", length(tables), nrow(tables[[1]])))
  tables
}

load_spike_in <- function(spike_in_path, sample_ids) {
  spike_in <- read.delim(spike_in_path, stringsAsFactors = FALSE)

  missing <- setdiff(sample_ids, spike_in$sample_id)
  if (length(missing) > 0) {
    stop("spike_in_rates.tsv is missing samples: ", paste(missing, collapse = ", "))
  }

  spike_in <- spike_in[match(sample_ids, spike_in$sample_id), , drop = FALSE]
  rownames(spike_in) <- NULL
  spike_in
}

# ---------------------------------------------------------------------------
# Matrix assembly and filtering
# ---------------------------------------------------------------------------

build_count_matrices <- function(tables) {
  # DSS needs integer read counts; meth_reads is fractional because methylKit
  # stores a percentage per site.
  meth <- vapply(tables, function(tab) as.numeric(tab$meth_reads), numeric(nrow(tables[[1]])))
  cov <- vapply(tables, function(tab) as.numeric(tab$total_coverage), numeric(nrow(tables[[1]])))

  list(M = round(meth), Cov = round(cov))
}

select_testable_genes <- function(tables, min_ca_sites, min_total_coverage) {
  ca_sites <- vapply(tables, function(tab) as.numeric(tab$n_ca_sites),
                     numeric(nrow(tables[[1]])))
  coverage <- vapply(tables, function(tab) as.numeric(tab$total_coverage),
                     numeric(nrow(tables[[1]])))

  passes_sites <- apply(ca_sites >= min_ca_sites, 1, all)
  passes_coverage <- apply(coverage >= min_total_coverage, 1, all)
  keep <- passes_sites & passes_coverage

  list(
    keep = keep,
    n_total = length(keep),
    n_keep = sum(keep),
    n_fail_sites = sum(!passes_sites),
    n_fail_coverage = sum(!passes_coverage),
    n_fail_both = sum(!passes_sites & !passes_coverage)
  )
}

# ---------------------------------------------------------------------------
# Design
# ---------------------------------------------------------------------------

build_design <- function(spike_in) {
  data.frame(
    genotype = factor(spike_in$genotype, levels = c("ctrl", "mut")),
    sex = factor(spike_in$sex, levels = c("M", "F")),
    lambda_ca_rate = as.numeric(spike_in$lambda_ca_rate),
    stringsAsFactors = FALSE
  )
}

check_design_rank <- function(design, formula, label) {
  mm <- model.matrix(formula, data = design)
  rank <- qr(mm)$rank

  message(sprintf("  %s: %d samples, %d coefficients, rank %d",
                  label, nrow(mm), ncol(mm), rank))
  message(sprintf("    coefficients: %s", paste(colnames(mm), collapse = ", ")))

  if (rank < ncol(mm)) {
    stop(label, " design matrix is rank deficient (rank ", rank, " < ",
         ncol(mm), " coefficients). Covariates are collinear; ",
         "the genotype effect is not identifiable.")
  }

  residual_df <- nrow(mm) - ncol(mm)
  message(sprintf("    residual df: %d", residual_df))
  invisible(mm)
}

report_covariate_collinearity <- function(design) {
  numeric_design <- data.frame(
    genotype_mut = as.integer(design$genotype == "mut"),
    sex_F = as.integer(design$sex == "F"),
    lambda_ca_rate = design$lambda_ca_rate
  )

  correlations <- cor(numeric_design)
  message("Covariate correlation matrix:")
  print(round(correlations, 3))

  genotype_lambda <- correlations["genotype_mut", "lambda_ca_rate"]
  message(sprintf("\ngenotype vs lambda_ca_rate correlation: %+.3f", genotype_lambda))
  if (abs(genotype_lambda) > 0.7) {
    message("  WARNING: the conversion noise floor tracks genotype closely. ",
            "The genotype coefficient and the noise covariate compete for the ",
            "same variance, which inflates the genotype standard error.")
  }

  correlations
}

# ---------------------------------------------------------------------------
# DSS fitting
# ---------------------------------------------------------------------------

build_bsseq <- function(M, Cov, genes, sample_ids) {
  # DSS identifies loci by chromosome and position. Gene bodies are not points,
  # so each gene gets a unique synthetic position and is mapped back afterwards.
  # Nothing downstream of the fit uses these coordinates as genomic distances.
  BSseq(chr = genes$chr,
        pos = seq_len(nrow(genes)),
        M = M,
        Cov = Cov,
        sampleNames = sample_ids)
}

fit_and_test <- function(bsobj, design, formula, coef_name, label) {
  message(sprintf("Fitting %s: %s", label, deparse(formula)))
  fit <- DMLfit.multiFactor(bsobj, design = design, formula = formula)

  available <- colnames(fit$X)
  if (!coef_name %in% available) {
    stop(label, ": coefficient '", coef_name, "' not in design matrix. ",
         "Available: ", paste(available, collapse = ", "))
  }

  test <- DMLtest.multiFactor(fit, coef = coef_name)

  # DSS returns rows sorted by chr/pos, not in input order. The synthetic pos
  # is the gene's row index, so ordering by it restores the input order.
  order_idx <- order(test$pos)
  test <- test[order_idx, , drop = FALSE]
  phi <- fit$fit$phi[order_idx]
  stopifnot(identical(test$pos, seq_len(nrow(test))))

  list(test = test, phi = phi)
}

run_leave_one_out <- function(M, Cov, genes, design, formula, coef_name,
                              sample_ids, alpha, n_tested) {
  message("Leave-one-sample-out refits:")
  bonferroni_cutoff <- alpha / n_tested
  significance <- matrix(NA, nrow = nrow(genes), ncol = length(sample_ids),
                         dimnames = list(NULL, sample_ids))

  for (i in seq_along(sample_ids)) {
    dropped <- sample_ids[i]
    keep_samples <- setdiff(seq_along(sample_ids), i)
    sub_design <- droplevels(design[keep_samples, , drop = FALSE])

    if (nlevels(sub_design$genotype) < 2L) {
      stop("Dropping ", dropped, " removes an entire genotype group")
    }

    sub_bsobj <- build_bsseq(M[, keep_samples, drop = FALSE],
                             Cov[, keep_samples, drop = FALSE],
                             genes, sample_ids[keep_samples])

    sub_result <- fit_and_test(sub_bsobj, sub_design, formula, coef_name,
                               sprintf("  without %s", dropped))
    significance[, i] <- sub_result$test$pvals < bonferroni_cutoff
  }

  rowSums(significance, na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# Annotation
# ---------------------------------------------------------------------------

flag_overlapping_significant <- function(genes, is_significant) {
  ranges <- GRanges(seqnames = genes$chr,
                    ranges = IRanges(start = genes$start + 1L, end = genes$end))

  significant_index <- which(is_significant)
  overlaps_flag <- rep(FALSE, nrow(genes))
  if (length(significant_index) < 2L) {
    return(overlaps_flag)
  }

  hits <- findOverlaps(ranges[significant_index], ranges[significant_index],
                       ignore.strand = TRUE)
  hits <- hits[queryHits(hits) != subjectHits(hits)]
  overlaps_flag[significant_index[unique(queryHits(hits))]] <- TRUE
  overlaps_flag
}

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

write_results <- function(results, output_dir) {
  path <- file.path(output_dir, "mca_differential_results.tsv")
  write.table(results, path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote ", path)
}

write_sample_matrix <- function(genes, rate_matrix, output_dir) {
  matrix_out <- cbind(genes[, c("gene_name", "gene_id", "chr", "start", "end",
                                "gene_length")],
                      as.data.frame(rate_matrix))
  path <- file.path(output_dir, "mca_sample_matrix.tsv")
  write.table(matrix_out, path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote ", path)
}

summarise_run <- function(con, opt, filter_stats, spike_in, correlations,
                          results, neuronal_genes) {
  emit <- function(...) cat(..., "\n", sep = "", file = con)

  emit("Gene-body mCA differential methylation")
  emit("======================================")
  emit("")
  emit("Genes considered:        ", filter_stats$n_total)
  emit("Genes tested:            ", filter_stats$n_keep)
  emit("Failed CA-site minimum:  ", filter_stats$n_fail_sites,
       "  (--min-ca-sites ", opt$min_ca_sites, ")")
  emit("Failed coverage minimum: ", filter_stats$n_fail_coverage,
       "  (--min-total-coverage ", opt$min_total_coverage, ")")
  emit("Failed both:             ", filter_stats$n_fail_both)
  emit("")

  emit("Conversion noise floor (lambda CA rate)")
  emit("---------------------------------------")
  for (i in seq_len(nrow(spike_in))) {
    emit(sprintf("  %-10s %s %s  %.4f%%",
                 spike_in$sample_id[i], spike_in$genotype[i], spike_in$sex[i],
                 spike_in$lambda_ca_rate[i] * 100))
  }
  emit("")

  emit("Covariate correlations")
  emit("----------------------")
  correlation_lines <- capture.output(print(round(correlations, 3)))
  for (line in correlation_lines) emit("  ", line)
  emit("")

  emit("Model: ~ genotype + sex + lambda_ca_rate")
  emit(sprintf("  Median dispersion (phi): %.6f", median(results$dispersion)))
  emit(sprintf("  Dispersion range:        [%.6f, %.6f]",
               min(results$dispersion), max(results$dispersion)))
  emit("")

  emit("Significance (reported model)")
  emit("-----------------------------")
  emit("  Bonferroni p < ", opt$alpha, ":  ", sum(results$sig_bonferroni))
  emit("  BH FDR < 0.05:      ", sum(results$p_bh < 0.05))
  emit("  BH FDR < 0.10:      ", sum(results$p_bh < 0.10))
  emit("")

  significant <- results[results$sig_bonferroni, , drop = FALSE]
  if (nrow(significant) > 0) {
    emit("  Direction among Bonferroni-significant genes:")
    emit("    mCA higher in mutant: ", sum(significant$mca_diff > 0))
    emit("    mCA lower in mutant:  ", sum(significant$mca_diff < 0))
    emit(sprintf("    median |mca_diff|:    %.5f%%",
                 median(abs(significant$mca_diff)) * 100))
    emit("    genes overlapping another significant gene: ",
         sum(significant$overlaps_sig_gene))
    emit("")
  }

  if (!is.null(neuronal_genes)) {
    tested_neuronal <- results$gene_name %in% neuronal_genes
    emit("Neuronal gene set")
    emit("-----------------")
    emit("  Gene set size:              ", length(neuronal_genes))
    emit("  Tested genes in set:        ", sum(tested_neuronal))
    emit("  Significant genes in set:   ", sum(tested_neuronal & results$sig_bonferroni))
    if (sum(results$sig_bonferroni) > 0 && sum(tested_neuronal) > 0) {
      contingency <- table(neuronal = tested_neuronal,
                           significant = results$sig_bonferroni)
      if (all(dim(contingency) == c(2, 2))) {
        fisher_result <- fisher.test(contingency)
        emit(sprintf("  Fisher odds ratio:          %.3f", fisher_result$estimate))
        emit(sprintf("  Fisher p-value:             %.3g", fisher_result$p.value))
      }
    }
    emit("")
  }

  if (!"loo_n_sig" %in% names(results)) {
    emit("Leave-one-out influence: skipped (--skip-loo)")
    emit("")
  } else if (nrow(significant) > 0) {
    emit("Leave-one-out influence (Bonferroni-significant genes)")
    emit("------------------------------------------------------")
    emit("  Robust in all 8 refits: ", sum(significant$loo_n_sig == 8))
    emit("  Lost in >= 1 refit:     ", sum(significant$loo_n_sig < 8))
    emit("  Lost in >= 3 refits:    ", sum(significant$loo_n_sig <= 5))
    emit("")
  }

  emit("Top 20 genes by p-value")
  emit("-----------------------")
  top <- head(results[order(results$dss_pval), ], 20)
  emit(sprintf("  %-16s %10s %10s %11s %11s %s",
               "gene", "length", "mca_ctrl", "mca_mut", "p", "bonf"))
  for (i in seq_len(nrow(top))) {
    emit(sprintf("  %-16s %10d %9.4f%% %10.4f%% %11.3g %s",
                 top$gene_name[i], top$gene_length[i],
                 top$mca_ctrl[i] * 100, top$mca_mut[i] * 100,
                 top$dss_pval[i], ifelse(top$sig_bonferroni[i], "*", "")))
  }
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  opt <- parse_cli_args()
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

  tables <- load_genebody_tables(opt$genebody_dir)
  sample_ids <- names(tables)
  spike_in <- load_spike_in(opt$spike_in, sample_ids)

  gene_metadata <- tables[[1]][, GENE_METADATA_COLUMNS]
  counts <- build_count_matrices(tables)

  filter_stats <- select_testable_genes(tables, opt$min_ca_sites,
                                        opt$min_total_coverage)
  message(sprintf("Testable genes: %d of %d", filter_stats$n_keep,
                  filter_stats$n_total))
  if (filter_stats$n_keep == 0L) {
    stop("No genes pass the coverage filters")
  }

  genes <- gene_metadata[filter_stats$keep, , drop = FALSE]
  rownames(genes) <- NULL
  M <- counts$M[filter_stats$keep, , drop = FALSE]
  Cov <- counts$Cov[filter_stats$keep, , drop = FALSE]
  colnames(M) <- sample_ids
  colnames(Cov) <- sample_ids

  design <- build_design(spike_in)
  message("")
  correlations <- report_covariate_collinearity(design)
  message("")

  model_formula <- ~ genotype + sex + lambda_ca_rate
  coef_name <- "genotypemut"

  message("Design check:")
  check_design_rank(design, model_formula, "Model")
  message("")

  bsobj <- build_bsseq(M, Cov, genes, sample_ids)
  fit_result <- fit_and_test(bsobj, design, model_formula, coef_name, "Model")
  reported <- fit_result$test
  phi <- fit_result$phi
  message("")

  rate_matrix <- M / Cov
  colnames(rate_matrix) <- paste0(sample_ids, "_mca")
  is_ctrl <- spike_in$genotype == "ctrl"
  mca_ctrl <- rowMeans(rate_matrix[, is_ctrl, drop = FALSE])
  mca_mut <- rowMeans(rate_matrix[, !is_ctrl, drop = FALSE])

  p_bonferroni <- pmin(reported$pvals * filter_stats$n_keep, 1)
  p_bh <- p.adjust(reported$pvals, method = "BH")
  sig_bonferroni <- p_bonferroni < opt$alpha

  ca_site_matrix <- vapply(tables, function(tab) as.numeric(tab$n_ca_sites),
                           numeric(nrow(tables[[1]])))[filter_stats$keep, , drop = FALSE]

  results <- data.frame(
    genes,
    n_ca_sites_mean = rowMeans(ca_site_matrix),
    total_coverage_mean = rowMeans(Cov),
    mca_ctrl = mca_ctrl,
    mca_mut = mca_mut,
    mca_diff = mca_mut - mca_ctrl,
    log2_fc = log2(mca_mut / mca_ctrl),
    as.data.frame(rate_matrix),
    dss_stat = reported$stat,
    dss_pval = reported$pvals,
    dss_fdr = reported$fdrs,
    p_bonferroni = p_bonferroni,
    p_bh = p_bh,
    sig_bonferroni = sig_bonferroni,
    sig_fdr = p_bh < opt$alpha,
    dispersion = phi,
    overlaps_sig_gene = flag_overlapping_significant(genes, sig_bonferroni),
    stringsAsFactors = FALSE
  )

  if (!opt$skip_loo) {
    message("")
    results$loo_n_sig <- run_leave_one_out(M, Cov, genes, design,
                                           model_formula, coef_name,
                                           sample_ids, opt$alpha,
                                           filter_stats$n_keep)
  }

  results <- results[order(results$dss_pval), , drop = FALSE]
  rownames(results) <- NULL

  neuronal_genes <- NULL
  if (!is.null(opt$neuronal_gene_set)) {
    neuronal_table <- read.delim(opt$neuronal_gene_set, stringsAsFactors = FALSE)
    neuronal_genes <- unique(neuronal_table[[1]])
    message(sprintf("Loaded %d neuronal genes from %s",
                    length(neuronal_genes), opt$neuronal_gene_set))
  }

  write_results(results, opt$output_dir)
  write_sample_matrix(genes, rate_matrix, opt$output_dir)

  summary_path <- file.path(opt$output_dir, "mca_analysis_summary.txt")
  con <- file(summary_path, open = "wt")
  on.exit(close(con), add = TRUE)
  summarise_run(con, opt, filter_stats, spike_in, correlations,
                results, neuronal_genes)
  message("Wrote ", summary_path)
}

main()
