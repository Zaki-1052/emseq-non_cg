# scripts/sections/20_chip_integration/20_04_mch_mecp2_by_mark.R
#
# Section 20_04: does the mCH / MeCP2 relationship depend on histone mark context?
#
# What this tests
#   MeCP2 reads methylated cytosines. Section 20_01 measured one gene-level
#   relationship between mCH change and MeCP2 binding change across all genes.
#   This section asks whether that relationship is the same in every chromatin
#   context. Genes are split into four mutually exclusive categories from the
#   promoter_state and body_state columns of section 10_01, and the mCH / MeCP2
#   relationship is measured inside each category.
#
# Analyses
#   1. Four-category stratification of the gene universe, with group sizes.
#   2. Per-category scatter of mCH change against MeCP2 fold change, carrying
#      Spearman rho, quadrant counts from assign_quadrant(), and labelled key
#      genes. The four panels share one x scale and one y scale.
#   3. Quadrant composition by category, with a chi-square test of the 4x4
#      contingency table.
#   4. Spearman rho per category with a Fisher z confidence interval, a
#      chi-square heterogeneity test across the four categories, and pairwise
#      Fisher z tests between categories.
#   5. Violin of mCH change and violin of MeCP2 fold change across the four
#      categories, each carrying n and the median on the figure, with
#      Kruskal-Wallis and pairwise Wilcoxon tests.
#   6. One registered gene-level Fisher test per category: mCH gain against
#      MeCP2 gain.
#
# Units
#   Every axis and every violin shows the mCH change in percentage points, the
#   column mch_diff_pct = 100 * mch_diff. The per-gene table carries both
#   columns. MeCP2 change is the DiffBind log2 fold change, mutant over control.
#
# Reads
#   HANDOFF_PATHS$gene_level_all_marks   gene-level table written by 20_02,
#                                        carrying promoter_state and body_state
#   mecp2_diffbind                       shared config
#
# Writes
#   Figures 20_04a to 20_04e and one TSV per figure into OUTPUT_PATHS$chip
#   (override with --output-dir), plus the per-gene table behind every figure.
#   One row per registered Fisher test into HANDOFF_PATHS$fisher_registry.
#
# Adapted from Biomodal section 68 (modC vs MeCP2 scatter by histone mark
# category), with mCH change on the x axis in place of total modified C.

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)

# =============================================================================
# SECTION CONSTANTS
# =============================================================================

SECTION_ID <- "20_04"

# Mark category rule. The categories come from the promoter_state and
# body_state columns that section 10_01 writes and section 20_02 carries into
# the handoff table. Assignment runs in this order, so the four categories stay
# mutually exclusive:
#   1. K27ac only    body_state == "Active"
#   2. Bivalent      promoter_state == "Bivalent_Promoter" OR body_state == "Mixed"
#   3. K27me3 only   body_state == "Polycomb" AND promoter_state != "Bivalent_Promoter"
#   4. Neither       every remaining gene
MARK_CATEGORY_ORDER <- c("K27ac only", "Bivalent", "K27me3 only", "Neither")

MARK_CATEGORY_SUBTITLE <- paste(
  "Categories from the chromatin state of section 10_01: K27ac only (body state",
  "Active), Bivalent (bivalent promoter or Mixed body), K27me3 only (Polycomb",
  "body, promoter not bivalent), Neither"
)

MARK_CATEGORY_COLORS <- c(
  "K27ac only"  = unname(COLORS$h3k27ac["H3K27ac Gained"]),
  "Bivalent"    = unname(PROMOTER_STATE_COLORS["Bivalent_Promoter"]),
  "K27me3 only" = unname(BODY_STATE_COLORS["Polycomb"]),
  "Neither"     = "grey55"
)

MARK_CATEGORY_SLUGS <- c(
  "K27ac only"  = "k27ac_only",
  "Bivalent"    = "bivalent",
  "K27me3 only" = "k27me3_only",
  "Neither"     = "neither"
)

MCH_DIRECTION_LEVELS <- c("Hypermethylated", "Hypomethylated")

QUADRANT_LEVELS <- c("Q1", "Q2", "Q3", "Q4")

QUADRANT_DESCRIPTIONS <- c(
  "Q1" = "mCH up, MeCP2 up",
  "Q2" = "mCH down, MeCP2 up",
  "Q3" = "mCH down, MeCP2 down",
  "Q4" = "mCH up, MeCP2 down"
)

# The 27 columns section 20_02 writes into gene_level_all_marks.tsv.
REQUIRED_HANDOFF_COLUMNS <- c(
  "gene_name", "gene_id", "chr", "start", "end", "gene_length",
  "mch_ctrl", "mch_mut", "mch_diff", "edger_logFC", "edger_fdr",
  "mch_sig", "mch_direction",
  "promoter_state", "body_state",
  "atac_fold", "atac_fdr", "atac_n_peaks",
  "k27ac_fold", "k27ac_fdr", "k27ac_n_peaks",
  "k27me3_fold", "k27me3_fdr", "k27me3_n_peaks",
  "k119ub_fold", "k119ub_fdr", "k119ub_n_peaks"
)

# Smallest category the section will analyse. Correlations, Fisher z tests and
# the per-category Fisher test all need more genes than this.
MIN_CATEGORY_SIZE <- 30

# Variance of the Fisher z transform of Spearman rho, as var(z) = FACTOR/(n-3).
# FACTOR is 1 for Pearson r; Fieller, Hartley and Pearson (1957) give 1.06 for
# Spearman rho. Every z statistic and z confidence interval below uses it.
SPEARMAN_Z_VARIANCE_FACTOR <- 1.06

# Axis clipping for the shared-scale scatter panels and the violins.
SCATTER_CLIP_PROBABILITY <- 0.995
VIOLIN_CLIP_PROBABILITIES <- c(0.005, 0.995)

# =============================================================================
# COMMAND LINE
# =============================================================================

parse_section_args <- function() {
  option_list <- list(
    make_option("--output-dir", dest = "output_dir", type = "character",
                default = OUTPUT_PATHS$chip,
                help = "Directory for figures and tables [default: %default]"),
    make_option("--fdr-threshold", dest = "fdr_threshold", type = "double",
                default = Q_THRESHOLD,
                help = "FDR cutoff for MeCP2 peak significance [default: %default]")
  )
  opt <- parse_args(OptionParser(option_list = option_list))

  if (opt$fdr_threshold <= 0 || opt$fdr_threshold >= 1) {
    stop("--fdr-threshold must be between 0 and 1, got ", opt$fdr_threshold)
  }
  opt
}

# =============================================================================
# SMALL UTILITIES
# =============================================================================

#' Format one p-value for a figure label.
format_p_value <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

#' Symmetric display limits from a quantile of the absolute values.
symmetric_display_limits <- function(x, probability = SCATTER_CLIP_PROBABILITY) {
  limit <- as.numeric(quantile(abs(x), probability, na.rm = TRUE))
  if (!is.finite(limit) || limit <= 0) {
    stop("Cannot build display limits: the ", probability,
         " quantile of the absolute values is ", limit)
  }
  c(-limit, limit)
}

#' Display limits from a lower and an upper quantile.
quantile_display_limits <- function(x, probabilities = VIOLIN_CLIP_PROBABILITIES) {
  limits <- as.numeric(quantile(x, probabilities, na.rm = TRUE))
  if (!all(is.finite(limits)) || limits[2] <= limits[1]) {
    stop("Cannot build display limits from quantiles: got ",
         paste(limits, collapse = " to "))
  }
  limits
}

#' Axis labels that carry the gene count of each category.
category_axis_labels <- function(df) {
  counts <- table(df$mark_category)
  setNames(
    sprintf("%s\n(n = %s)", names(counts),
            format(as.integer(counts), big.mark = ",", trim = TRUE)),
    names(counts)
  )
}

# =============================================================================
# INPUT LOADING
# =============================================================================

#' Read the gene-level master table written by section 20_02.
load_gene_level_all_marks <- function(path) {
  if (!file.exists(path)) {
    stop("Gene-level master table not found: ", path,
         "\nRun section 20_02 (20_02_multi_mark_diffbind.R) first. ",
         "It writes gene_level_all_marks.tsv.")
  }

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")

  missing <- setdiff(REQUIRED_HANDOFF_COLUMNS, colnames(df))
  if (length(missing) > 0) {
    stop("Gene-level master table is missing columns: ",
         paste(missing, collapse = ", "), " (", path, ")",
         "\nSection 20_02 writes this file and its column contract.")
  }

  df$mch_sig <- as.logical(df$mch_sig)
  if (any(is.na(df$mch_sig))) {
    stop(sum(is.na(df$mch_sig)), " rows of ", path,
         " have a mch_sig value that is not TRUE or FALSE.")
  }
  for (col in c("start", "end", "mch_diff", "edger_logFC", "edger_fdr")) {
    df[[col]] <- as.numeric(df[[col]])
    if (any(is.na(df[[col]]))) {
      stop(sum(is.na(df[[col]])), " rows of ", path,
           " have a missing value in column ", col)
    }
  }

  cat(sprintf("  Gene-level master table: %s genes (%s mCH significant)\n",
              format(nrow(df), big.mark = ","),
              format(sum(df$mch_sig), big.mark = ",")))
  df
}

#' Keep one row per gene name, the row with the largest absolute edgeR logFC.
#'
#' MeCP2 fold change joins on gene name, so two rows that share a name would
#' carry the same MeCP2 value into the correlations and the Fisher tables.
deduplicate_by_gene_name <- function(df) {
  ord <- order(df$gene_name, -abs(df$edger_logFC))
  keep <- sort(ord[!duplicated(df$gene_name[ord])])
  out <- df[keep, , drop = FALSE]

  cat(sprintf("  Deduplicated %s rows to %s gene names (%s rows dropped)\n",
              format(nrow(df), big.mark = ","),
              format(nrow(out), big.mark = ","),
              format(nrow(df) - nrow(out), big.mark = ",")))
  out
}

#' Restrict the MeCP2 DiffBind peaks to canonical chromosomes and complete rows.
prepare_mecp2_peaks <- function(db, fdr_threshold) {
  n_all <- nrow(db)
  db <- db[db$Chr %in% CANONICAL_CHRS, , drop = FALSE]
  db <- db[!is.na(db$Fold) & !is.na(db$FDR), , drop = FALSE]

  if (nrow(db) == 0) {
    stop("No usable MeCP2 peaks after filtering to canonical chromosomes and ",
         "complete Fold and FDR values.")
  }

  cat(sprintf("  MeCP2 peaks: %s in, %s kept (%s gained, %s lost at FDR<%.3g)\n",
              format(n_all, big.mark = ","), format(nrow(db), big.mark = ","),
              format(sum(db$FDR < fdr_threshold & db$Fold > 0), big.mark = ","),
              format(sum(db$FDR < fdr_threshold & db$Fold < 0), big.mark = ","),
              fdr_threshold))
  db
}

#' Collapse the MeCP2 peaks to one fold change per gene by the nearest-TSS rule.
aggregate_mecp2_to_genes <- function(db, fdr_threshold) {
  annotated <- annotate_peaks_to_genes(db, "MeCP2")

  usable <- !is.na(annotated$SYMBOL) & nzchar(annotated$SYMBOL) &
    !is.na(annotated$Fold) & !is.na(annotated$FDR)
  cat(sprintf("  %s of %s annotated MeCP2 peaks carry a gene symbol\n",
              format(sum(usable), big.mark = ","),
              format(nrow(annotated), big.mark = ",")))
  annotated <- annotated[usable, , drop = FALSE]

  gene_table <- aggregate_diffbind_by_gene(annotated, method = "nearest_tss",
                                           fdr_threshold = fdr_threshold,
                                           prefix = "mecp2")
  if (anyDuplicated(gene_table$gene_name) > 0) {
    stop("aggregate_diffbind_by_gene() returned duplicate gene names for MeCP2.")
  }

  cat(sprintf("  Genes with a MeCP2 fold change: %s\n",
              format(nrow(gene_table), big.mark = ",")))
  gene_table
}

# =============================================================================
# MARK CATEGORY STRATIFICATION
# =============================================================================

#' Add the four-way mark category from promoter_state and body_state.
#'
#' The rule, applied in this order so the four categories are mutually
#' exclusive:
#'   K27ac only    body_state == "Active"
#'   Bivalent      promoter_state == "Bivalent_Promoter" OR body_state == "Mixed"
#'   K27me3 only   body_state == "Polycomb" AND promoter_state != "Bivalent_Promoter"
#'   Neither       every remaining gene
#'
#' @param df Gene-level table carrying promoter_state and body_state.
#' @return the same table with a mark_category factor added
classify_mark_category <- function(df) {
  unknown_promoter <- setdiff(unique(df$promoter_state), PROMOTER_STATE_ORDER)
  if (length(unknown_promoter) > 0) {
    stop("Unknown promoter_state value(s) in the handoff table: ",
         paste(unknown_promoter, collapse = ", "),
         "\nAllowed: ", paste(PROMOTER_STATE_ORDER, collapse = ", "),
         "\nMissing values mean section 10_01 did not cover every gene.")
  }
  unknown_body <- setdiff(unique(df$body_state), BODY_STATE_ORDER)
  if (length(unknown_body) > 0) {
    stop("Unknown body_state value(s) in the handoff table: ",
         paste(unknown_body, collapse = ", "),
         "\nAllowed: ", paste(BODY_STATE_ORDER, collapse = ", "),
         "\nMissing values mean section 10_01 did not cover every gene.")
  }

  category <- rep("Neither", nrow(df))

  is_k27ac <- df$body_state == "Active"
  category[is_k27ac] <- "K27ac only"

  is_bivalent <- !is_k27ac &
    (df$promoter_state == "Bivalent_Promoter" | df$body_state == "Mixed")
  category[is_bivalent] <- "Bivalent"

  is_k27me3 <- !is_k27ac & !is_bivalent &
    df$body_state == "Polycomb" & df$promoter_state != "Bivalent_Promoter"
  category[is_k27me3] <- "K27me3 only"

  df$mark_category <- factor(category, levels = MARK_CATEGORY_ORDER)

  if (any(is.na(df$mark_category))) {
    stop(sum(is.na(df$mark_category)), " genes fell outside the four mark ",
         "categories. The classification rule must cover every gene.")
  }
  df
}

#' Build the analysis table: deduplicated genes with a mark category and a
#' MeCP2 fold change.
build_analysis_table <- function(genes, mecp2_gene, fdr_threshold) {
  merged <- dplyr::left_join(genes, mecp2_gene, by = "gene_name")
  if (nrow(merged) != nrow(genes)) {
    stop("The MeCP2 join changed the row count: ", nrow(genes), " in, ",
         nrow(merged), " out.")
  }

  merged$has_mecp2 <- !is.na(merged$mecp2_fold)
  df <- merged[merged$has_mecp2 & is.finite(merged$mch_diff), , drop = FALSE]

  df$mch_diff_pct <- 100 * df$mch_diff
  df$mch_direction <- factor(df$mch_direction, levels = MCH_DIRECTION_LEVELS)
  if (any(is.na(df$mch_direction))) {
    stop(sum(is.na(df$mch_direction)), " genes carry an mch_direction value ",
         "outside ", paste(MCH_DIRECTION_LEVELS, collapse = " and "), ".")
  }

  df$quadrant <- factor(assign_quadrant(df$mch_diff_pct, df$mecp2_fold),
                        levels = QUADRANT_LEVELS)
  df$mch_gain <- df$mch_diff_pct > 0
  df$mecp2_gain <- df$mecp2_fold > 0
  df$mecp2_sig_gain <- !is.na(df$mecp2_fdr) & df$mecp2_fdr < fdr_threshold &
    df$mecp2_fold > 0
  df$mecp2_sig_loss <- !is.na(df$mecp2_fdr) & df$mecp2_fdr < fdr_threshold &
    df$mecp2_fold < 0

  cat(sprintf("  Analysis table: %s of %s genes carry both mCH and MeCP2 data\n",
              format(nrow(df), big.mark = ","),
              format(nrow(merged), big.mark = ",")))
  df
}

#' Report the group sizes of the four mark categories and check they are usable.
report_category_sizes <- function(all_genes, analysis_df, out_dir) {
  totals <- as.data.frame(table(all_genes$mark_category),
                          stringsAsFactors = FALSE)
  colnames(totals) <- c("mark_category", "n_genes_total")

  with_mecp2 <- as.data.frame(table(analysis_df$mark_category),
                              stringsAsFactors = FALSE)
  colnames(with_mecp2) <- c("mark_category", "n_genes_with_mecp2")

  sizes <- dplyr::left_join(totals, with_mecp2, by = "mark_category")
  sizes$pct_of_universe <- 100 * sizes$n_genes_total / sum(sizes$n_genes_total)
  sizes$pct_with_mecp2 <- 100 * sizes$n_genes_with_mecp2 / sizes$n_genes_total
  sizes$n_bivalent_promoter <- as.integer(tapply(
    all_genes$promoter_state == "Bivalent_Promoter",
    all_genes$mark_category, sum)[sizes$mark_category])
  for (state in BODY_STATE_ORDER) {
    sizes[[paste0("n_body_", tolower(state))]] <- as.integer(tapply(
      all_genes$body_state == state,
      all_genes$mark_category, sum)[sizes$mark_category])
  }

  for (i in seq_len(nrow(sizes))) {
    cat(sprintf("  %-12s %s genes (%.1f%% of universe), %s with MeCP2 (%.1f%%)\n",
                sizes$mark_category[i],
                format(sizes$n_genes_total[i], big.mark = ","),
                sizes$pct_of_universe[i],
                format(sizes$n_genes_with_mecp2[i], big.mark = ","),
                sizes$pct_with_mecp2[i]))
  }

  too_small <- sizes$mark_category[sizes$n_genes_with_mecp2 < MIN_CATEGORY_SIZE]
  if (length(too_small) > 0) {
    stop("These mark categories carry fewer than ", MIN_CATEGORY_SIZE,
         " genes with MeCP2 data: ", paste(too_small, collapse = ", "),
         ". The per-category correlations and Fisher tests need more genes.")
  }

  write_section_table(sizes, file.path(out_dir, "20_04_mark_category_sizes.tsv"))
  sizes
}

#' Cross-tabulate the mark category against both chromatin state columns.
#'
#' The categories are built from these two columns, so the crosstabs show
#' exactly which state combinations landed in each category.
write_category_state_crosstabs <- function(analysis_df, out_dir) {
  promoter_tab <- table(
    mark_category = analysis_df$mark_category,
    promoter_state = factor(analysis_df$promoter_state,
                            levels = PROMOTER_STATE_ORDER)
  )
  promoter_long <- as.data.frame(promoter_tab, stringsAsFactors = FALSE)
  colnames(promoter_long) <- c("mark_category", "promoter_state", "n_genes")
  write_section_table(promoter_long,
                      file.path(out_dir, "20_04_mark_category_vs_promoter_state.tsv"))

  body_tab <- table(
    mark_category = analysis_df$mark_category,
    body_state = factor(analysis_df$body_state, levels = BODY_STATE_ORDER)
  )
  body_long <- as.data.frame(body_tab, stringsAsFactors = FALSE)
  colnames(body_long) <- c("mark_category", "body_state", "n_genes")
  write_section_table(body_long,
                      file.path(out_dir, "20_04_mark_category_vs_body_state.tsv"))

  list(promoter = promoter_long, body = body_long)
}

# =============================================================================
# PER-CATEGORY CORRELATIONS
# =============================================================================

#' Spearman rho of mCH change against MeCP2 fold change inside each category,
#' with a Fisher z confidence interval.
correlation_by_category <- function(df) {
  rows <- lapply(MARK_CATEGORY_ORDER, function(category) {
    sub <- df[df$mark_category == category, , drop = FALSE]
    if (nrow(sub) < MIN_CATEGORY_SIZE) {
      stop("Category ", category, " has ", nrow(sub),
           " genes with both values, fewer than ", MIN_CATEGORY_SIZE, ".")
    }

    test <- suppressWarnings(cor.test(sub$mch_diff_pct, sub$mecp2_fold,
                                      method = "spearman", exact = FALSE))
    pearson <- cor.test(sub$mch_diff_pct, sub$mecp2_fold, method = "pearson")

    rho <- unname(test$estimate)
    z <- atanh(rho)
    se_z <- sqrt(SPEARMAN_Z_VARIANCE_FACTOR / (nrow(sub) - 3))

    data.frame(
      mark_category = category,
      n_genes = nrow(sub),
      spearman_rho = rho,
      spearman_p = test$p.value,
      fisher_z = z,
      fisher_z_se = se_z,
      rho_ci_lower = tanh(z - qnorm(0.975) * se_z),
      rho_ci_upper = tanh(z + qnorm(0.975) * se_z),
      pearson_r = unname(pearson$estimate),
      pearson_p = pearson$p.value,
      median_mch_diff_pct = median(sub$mch_diff_pct),
      median_mecp2_fold = median(sub$mecp2_fold),
      z_variance_factor = SPEARMAN_Z_VARIANCE_FACTOR,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  for (i in seq_len(nrow(out))) {
    cat(sprintf("  %-12s n = %s, rho = %+.3f [%+.3f, %+.3f], p = %.3g\n",
                out$mark_category[i], format(out$n_genes[i], big.mark = ","),
                out$spearman_rho[i], out$rho_ci_lower[i], out$rho_ci_upper[i],
                out$spearman_p[i]))
  }
  out
}

#' Compare the four correlations with the Fisher z transformation.
#'
#' Returns the six pairwise z tests with Holm-adjusted p-values and one
#' chi-square heterogeneity test over all four categories.
compare_correlations_fisher_z <- function(cor_tbl) {
  pairs <- combn(cor_tbl$mark_category, 2, simplify = FALSE)

  pairwise <- do.call(rbind, lapply(pairs, function(pair) {
    a <- cor_tbl[cor_tbl$mark_category == pair[1], ]
    b <- cor_tbl[cor_tbl$mark_category == pair[2], ]
    se_diff <- sqrt(a$fisher_z_se^2 + b$fisher_z_se^2)
    z_stat <- (a$fisher_z - b$fisher_z) / se_diff

    data.frame(
      category_1 = pair[1], category_2 = pair[2],
      n_1 = a$n_genes, n_2 = b$n_genes,
      rho_1 = a$spearman_rho, rho_2 = b$spearman_rho,
      rho_difference = a$spearman_rho - b$spearman_rho,
      z_1 = a$fisher_z, z_2 = b$fisher_z,
      z_difference = a$fisher_z - b$fisher_z,
      se_difference = se_diff,
      z_statistic = z_stat,
      p_value = 2 * pnorm(-abs(z_stat)),
      stringsAsFactors = FALSE
    )
  }))
  pairwise$p_holm <- p.adjust(pairwise$p_value, method = "holm")
  pairwise$pair_label <- sprintf("%s vs %s", pairwise$category_1,
                                 pairwise$category_2)

  weights <- 1 / cor_tbl$fisher_z_se^2
  z_bar <- sum(weights * cor_tbl$fisher_z) / sum(weights)
  q_statistic <- sum(weights * (cor_tbl$fisher_z - z_bar)^2)
  q_df <- nrow(cor_tbl) - 1L

  overall <- data.frame(
    test = "chi_square_heterogeneity_of_fisher_z",
    n_categories = nrow(cor_tbl),
    n_genes = sum(cor_tbl$n_genes),
    pooled_fisher_z = z_bar,
    pooled_rho = tanh(z_bar),
    chi_squared = q_statistic,
    df = q_df,
    p_value = pchisq(q_statistic, df = q_df, lower.tail = FALSE),
    z_variance_factor = SPEARMAN_Z_VARIANCE_FACTOR,
    stringsAsFactors = FALSE
  )

  cat(sprintf("  Heterogeneity of the four rho values: chi-squared = %.2f, df = %d, p = %.3g\n",
              overall$chi_squared, overall$df, overall$p_value))
  for (i in seq_len(nrow(pairwise))) {
    cat(sprintf("    %-28s z = %+.2f, p = %.3g (Holm %.3g)\n",
                pairwise$pair_label[i], pairwise$z_statistic[i],
                pairwise$p_value[i], pairwise$p_holm[i]))
  }

  list(pairwise = pairwise, overall = overall)
}

# =============================================================================
# QUADRANT COUNTS
# =============================================================================

#' Quadrant counts and percentages inside each mark category.
quadrant_counts_by_category <- function(df) {
  tab <- table(mark_category = df$mark_category, quadrant = df$quadrant)
  counts <- as.data.frame(tab, stringsAsFactors = FALSE)
  colnames(counts) <- c("mark_category", "quadrant", "n_genes")

  totals <- tapply(counts$n_genes, counts$mark_category, sum)
  counts$n_category <- as.integer(totals[counts$mark_category])
  counts$percentage <- 100 * counts$n_genes / counts$n_category
  counts$quadrant_meaning <- unname(QUADRANT_DESCRIPTIONS[counts$quadrant])
  counts$mark_category <- factor(counts$mark_category,
                                 levels = MARK_CATEGORY_ORDER)
  counts$quadrant <- factor(counts$quadrant, levels = QUADRANT_LEVELS)

  chisq <- chisq.test(tab)
  cat(sprintf("  Quadrant composition across categories: chi-squared = %.1f, df = %d, p = %.3g\n",
              unname(chisq$statistic), unname(chisq$parameter), chisq$p.value))

  list(counts = counts[order(counts$mark_category, counts$quadrant), ],
       chisq = chisq, table = tab)
}

# =============================================================================
# FIGURE 20_04a: PER-CATEGORY SCATTER WITH A SHARED SCALE
# =============================================================================

#' Corner positions for the four quadrant count labels.
quadrant_label_positions <- function(x_lim, y_lim) {
  pad_x <- diff(x_lim) * 0.03
  pad_y <- diff(y_lim) * 0.03
  data.frame(
    quadrant = QUADRANT_LEVELS,
    x = c(x_lim[2] - pad_x, x_lim[1] + pad_x, x_lim[1] + pad_x, x_lim[2] - pad_x),
    y = c(y_lim[2] - pad_y, y_lim[2] - pad_y, y_lim[1] + pad_y, y_lim[1] + pad_y),
    hjust = c(1, 0, 0, 1),
    vjust = c(1, 1, 0, 0),
    stringsAsFactors = FALSE
  )
}

#' One scatter panel: mCH change against MeCP2 fold change in one category.
build_scatter_panel <- function(df, category, cor_row, quad_counts,
                                x_lim, y_lim) {
  sub <- df[df$mark_category == category, , drop = FALSE]

  labels <- quadrant_label_positions(x_lim, y_lim)
  quad_sub <- quad_counts[quad_counts$mark_category == category,
                          c("quadrant", "n_genes", "percentage"), drop = FALSE]
  quad_sub$quadrant <- as.character(quad_sub$quadrant)
  labels <- dplyr::left_join(labels, quad_sub, by = "quadrant")
  if (any(is.na(labels$n_genes))) {
    stop("Quadrant counts are missing for category ", category, ": ",
         paste(labels$quadrant[is.na(labels$n_genes)], collapse = ", "))
  }
  labels$label <- sprintf("%s\nn = %s (%.1f%%)", labels$quadrant,
                          format(labels$n_genes, big.mark = ",", trim = TRUE),
                          labels$percentage)

  panel <- ggplot(sub, aes(x = mch_diff_pct, y = mecp2_fold)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.3) +
    geom_point(aes(color = mch_direction), alpha = 0.35, size = 0.8) +
    geom_text(data = labels,
              aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
              inherit.aes = FALSE, size = 2.7, color = "grey25",
              lineheight = 1.05) +
    scale_color_manual(values = COLORS$direction, name = "mCH direction",
                       drop = FALSE) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    labs(
      title = category,
      subtitle = sprintf("Spearman rho = %.3f [%.3f, %.3f], %s | n = %s genes",
                         cor_row$spearman_rho, cor_row$rho_ci_lower,
                         cor_row$rho_ci_upper, format_p_value(cor_row$spearman_p),
                         format(cor_row$n_genes, big.mark = ",")),
      x = "mCH change (mutant - control, percentage points)",
      y = "MeCP2 log2 fold change (mutant / control)"
    ) +
    theme_emseq() +
    theme(plot.title = element_text(color = MARK_CATEGORY_COLORS[[category]],
                                    face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 9, color = "grey35"))

  key_df <- sub[sub$gene_name %in% KEY_GENES, , drop = FALSE]
  if (nrow(key_df) == 0) return(panel)

  panel +
    geom_point(data = key_df, size = 1.7, color = "black") +
    ggrepel::geom_text_repel(data = key_df, aes(label = gene_name),
                             size = 2.6, fontface = "italic", color = "black",
                             segment.color = "grey50", segment.size = 0.3,
                             min.segment.length = 0, max.overlaps = Inf)
}

plot_scatter_panels <- function(df, cor_tbl, quad_counts, out_dir) {
  cat("--- Figure 20_04a: mCH change against MeCP2 fold change by mark category ---\n")

  x_lim <- symmetric_display_limits(df$mch_diff_pct)
  y_lim <- symmetric_display_limits(df$mecp2_fold)
  cat(sprintf("  Shared axes: x %.3f to %.3f, y %.3f to %.3f\n",
              x_lim[1], x_lim[2], y_lim[1], y_lim[2]))

  panels <- lapply(MARK_CATEGORY_ORDER, function(category) {
    cor_row <- cor_tbl[cor_tbl$mark_category == category, ]
    panel <- build_scatter_panel(df, category, cor_row, quad_counts,
                                 x_lim, y_lim)
    save_multiformat_ggplot(
      panel,
      file.path(out_dir, sprintf("20_04a_scatter_%s",
                                 MARK_CATEGORY_SLUGS[[category]])),
      width = 8, height = 7
    )
    panel
  })
  names(panels) <- MARK_CATEGORY_ORDER

  caption <- sprintf(paste("All four panels share one x scale and one y scale,",
                           "clipped to the %.1fth percentile of the absolute",
                           "value. Every gene enters the correlation and the",
                           "quadrant counts."),
                     100 * SCATTER_CLIP_PROBABILITY)

  figure <- wrap_plots(panels, nrow = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "mCH change against MeCP2 binding change by histone mark context",
      subtitle = MARK_CATEGORY_SUBTITLE,
      caption = caption,
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
        plot.caption = element_text(hjust = 0.5, size = 8, color = "grey45")
      )
    ) &
    theme(legend.position = "bottom")

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "20_04a_mch_mecp2_by_mark"),
                          width = 15, height = 13)

  scatter_stats <- cor_tbl[, c("mark_category", "n_genes", "spearman_rho",
                               "spearman_p", "rho_ci_lower", "rho_ci_upper",
                               "pearson_r", "pearson_p")]
  scatter_stats$x_limit_lower <- x_lim[1]
  scatter_stats$x_limit_upper <- x_lim[2]
  scatter_stats$y_limit_lower <- y_lim[1]
  scatter_stats$y_limit_upper <- y_lim[2]
  write_section_table(scatter_stats,
                      file.path(out_dir, "20_04a_scatter_statistics.tsv"))
  cat("\n")
  invisible(panels)
}

# =============================================================================
# FIGURE 20_04b: QUADRANT COMPOSITION
# =============================================================================

plot_quadrant_composition <- function(df, quad_result, out_dir) {
  cat("--- Figure 20_04b: quadrant composition by mark category ---\n")

  counts <- quad_result$counts
  chisq <- quad_result$chisq
  axis_labels <- category_axis_labels(df)

  quadrant_labels <- sprintf("%s (%s)", QUADRANT_LEVELS,
                             unname(QUADRANT_DESCRIPTIONS[QUADRANT_LEVELS]))
  names(quadrant_labels) <- QUADRANT_LEVELS

  figure <- ggplot(counts, aes(x = mark_category, y = percentage,
                               fill = quadrant)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.1f%%\n(%s)", percentage,
                                  format(n_genes, big.mark = ",", trim = TRUE))),
              position = position_stack(vjust = 0.5), size = 3,
              lineheight = 1.0) +
    scale_fill_manual(values = COLORS$quadrant, name = "Quadrant",
                      labels = quadrant_labels) +
    scale_x_discrete(labels = axis_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = "Quadrant composition of the mCH / MeCP2 relationship by mark context",
      subtitle = sprintf(paste("%s | chi-square of the 4x4 table:",
                               "chi-squared = %.1f, df = %d, p = %.3g"),
                         MARK_CATEGORY_SUBTITLE, unname(chisq$statistic),
                         unname(chisq$parameter), chisq$p.value),
      x = "", y = "Percent of genes in the category"
    ) +
    theme_emseq() +
    theme(legend.position = "right",
          plot.subtitle = element_text(size = 9, color = "grey40"))

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "20_04b_quadrant_composition"),
                          width = 12, height = 8)

  counts_out <- counts
  counts_out$chi_squared <- unname(chisq$statistic)
  counts_out$chi_squared_df <- unname(chisq$parameter)
  counts_out$chi_squared_p <- chisq$p.value
  write_section_table(counts_out,
                      file.path(out_dir, "20_04b_quadrant_counts.tsv"))

  expected <- as.data.frame(as.table(chisq$expected), stringsAsFactors = FALSE)
  colnames(expected) <- c("mark_category", "quadrant", "expected_n_genes")
  write_section_table(expected,
                      file.path(out_dir, "20_04b_quadrant_expected_counts.tsv"))
  cat("\n")
  invisible(counts_out)
}

# =============================================================================
# FIGURE 20_04c: DO THE CORRELATIONS DIFFER BETWEEN CATEGORIES?
# =============================================================================

plot_correlation_comparison <- function(cor_tbl, comparison, out_dir) {
  cat("--- Figure 20_04c: correlation comparison across mark categories ---\n")

  forest_df <- cor_tbl
  forest_df$mark_category <- factor(forest_df$mark_category,
                                    levels = rev(MARK_CATEGORY_ORDER))
  forest_df$label <- sprintf("rho = %.3f [%.3f, %.3f], n = %s",
                             forest_df$spearman_rho, forest_df$rho_ci_lower,
                             forest_df$rho_ci_upper,
                             format(forest_df$n_genes, big.mark = ","))

  p_forest <- ggplot(forest_df, aes(x = spearman_rho, y = mark_category,
                                    color = mark_category)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(xmin = rho_ci_lower, xmax = rho_ci_upper),
                  orientation = "y", width = 0.18, linewidth = 0.8) +
    geom_point(size = 3.5) +
    geom_text(aes(label = label), vjust = -1.3, size = 3.1, show.legend = FALSE) +
    scale_color_manual(values = MARK_CATEGORY_COLORS, guide = "none") +
    scale_y_discrete(expand = expansion(add = 0.7)) +
    labs(
      title = "Spearman rho per mark category",
      subtitle = sprintf("95%% confidence interval from the Fisher z transform, var(z) = %.2f/(n-3)",
                         SPEARMAN_Z_VARIANCE_FACTOR),
      x = "Spearman rho (mCH change against MeCP2 fold change)", y = ""
    ) +
    theme_emseq() +
    theme(axis.text.y = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"))

  pairwise <- comparison$pairwise
  pairwise$abs_z <- abs(pairwise$z_statistic)
  pairwise <- pairwise[order(pairwise$abs_z), ]
  pairwise$pair_label <- factor(pairwise$pair_label, levels = pairwise$pair_label)
  pairwise$significance <- ifelse(pairwise$p_holm < 0.05, "Holm p < 0.05",
                                  "Holm p >= 0.05")
  critical_value <- qnorm(0.975)

  p_pairs <- ggplot(pairwise, aes(x = abs_z, y = pair_label,
                                  color = significance)) +
    geom_vline(xintercept = critical_value, linetype = "dashed",
               color = "grey40") +
    geom_segment(aes(x = 0, xend = abs_z, y = pair_label, yend = pair_label),
                 linewidth = 0.7) +
    geom_point(size = 3.2) +
    geom_text(aes(label = sprintf("z = %+.2f, p = %.3g (Holm %.3g)",
                                  z_statistic, p_value, p_holm)),
              hjust = -0.08, vjust = -1.0, size = 3, show.legend = FALSE) +
    scale_color_manual(values = c("Holm p < 0.05" = "#E41A1C",
                                  "Holm p >= 0.05" = "grey45"),
                       name = "") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.45))) +
    scale_y_discrete(expand = expansion(add = 0.7)) +
    labs(
      title = "Pairwise Fisher z tests between categories",
      subtitle = sprintf("Dashed line is the two-sided 5%% critical value, |z| = %.2f",
                         critical_value),
      x = "|z| for the difference of two Fisher z values", y = ""
    ) +
    theme_emseq() +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(size = 9, color = "grey40"))

  overall <- comparison$overall
  figure <- p_forest + p_pairs +
    plot_annotation(
      title = "Does the mCH / MeCP2 correlation depend on histone mark context?",
      subtitle = sprintf(paste("Heterogeneity of the four correlations:",
                               "chi-squared = %.2f, df = %d, p = %.3g |",
                               "pooled rho = %.3f over %s genes"),
                         overall$chi_squared, overall$df, overall$p_value,
                         overall$pooled_rho,
                         format(overall$n_genes, big.mark = ",")),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40")
      )
    )

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "20_04c_correlation_comparison"),
                          width = 16, height = 7)

  write_section_table(cor_tbl,
                      file.path(out_dir, "20_04c_correlation_by_category.tsv"))
  write_section_table(comparison$pairwise,
                      file.path(out_dir, "20_04c_correlation_pairwise_fisher_z.tsv"))
  write_section_table(comparison$overall,
                      file.path(out_dir, "20_04c_correlation_heterogeneity.tsv"))
  cat("\n")
  invisible(figure)
}

# =============================================================================
# FIGURES 20_04d AND 20_04e: DISTRIBUTIONS ACROSS THE FOUR CATEGORIES
# =============================================================================

#' Kruskal-Wallis test plus every pairwise Wilcoxon test across the categories.
distribution_tests <- function(df, value_col) {
  tmp <- data.frame(value = df[[value_col]], group = df$mark_category)
  kruskal <- kruskal.test(value ~ group, data = tmp)

  pairs <- combn(MARK_CATEGORY_ORDER, 2, simplify = FALSE)
  pairwise <- do.call(rbind, lapply(pairs, function(pair) {
    a <- tmp$value[tmp$group == pair[1]]
    b <- tmp$value[tmp$group == pair[2]]
    wt <- suppressWarnings(wilcox.test(a, b))
    data.frame(
      value_column = value_col,
      group_1 = pair[1], group_2 = pair[2],
      n_1 = length(a), n_2 = length(b),
      median_1 = median(a), median_2 = median(b),
      median_difference = median(a) - median(b),
      W = unname(wt$statistic), p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  }))
  pairwise$p_holm <- p.adjust(pairwise$p_value, method = "holm")

  list(kruskal = kruskal, pairwise = pairwise)
}

#' Violin of one value across the four mark categories.
#'
#' Every group carries its gene count and its median on the figure:
#' summarise_groups() computes them and group_label() writes the text. The y
#' axis is clipped to the quantiles in VIOLIN_CLIP_PROBABILITIES; the violin
#' shape and every statistic use all genes.
plot_violin_by_category <- function(df, value_col, value_label, figure_id,
                                    plot_title, out_dir, digits) {
  cat(sprintf("--- Figure %s: %s by mark category ---\n", figure_id, value_label))

  group_stats <- summarise_groups(df, "mark_category", value_col)
  group_stats$mark_category <- factor(group_stats$mark_category,
                                      levels = MARK_CATEGORY_ORDER)

  tests <- distribution_tests(df, value_col)
  kruskal <- tests$kruskal

  for (i in seq_len(nrow(group_stats))) {
    cat(sprintf("  %-12s n = %s, median = %.5f, mean = %.5f\n",
                as.character(group_stats$mark_category[i]),
                format(group_stats$n[i], big.mark = ","),
                group_stats$median[i], group_stats$mean[i]))
  }
  cat(sprintf("  Kruskal-Wallis chi-squared = %.1f, df = %d, p = %.3g\n",
              unname(kruskal$statistic), unname(kruskal$parameter),
              kruskal$p.value))

  y_limits <- quantile_display_limits(df[[value_col]])
  y_span <- diff(y_limits)
  y_upper <- y_limits[2] + 0.22 * y_span
  y_lower <- y_limits[1] - 0.05 * y_span
  group_stats$label_y <- y_limits[2] + 0.10 * y_span
  group_stats$label <- group_label(group_stats, digits = digits)

  axis_labels <- category_axis_labels(df)

  figure <- ggplot(df, aes(x = mark_category, y = .data[[value_col]],
                           fill = mark_category)) +
    geom_violin(alpha = 0.6, scale = "width",
                draw_quantiles = c(0.25, 0.5, 0.75), color = "grey25") +
    geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white", alpha = 0.85,
                 linewidth = 0.35) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black",
               linewidth = 0.4) +
    geom_text(data = group_stats,
              aes(x = mark_category, y = label_y, label = label),
              inherit.aes = FALSE, size = 3.2, lineheight = 1.1) +
    scale_fill_manual(values = MARK_CATEGORY_COLORS, guide = "none") +
    scale_x_discrete(labels = axis_labels) +
    coord_cartesian(ylim = c(y_lower, y_upper)) +
    labs(
      title = plot_title,
      subtitle = sprintf(paste("%s | Kruskal-Wallis chi-squared = %.1f,",
                               "df = %d, p = %.3g"),
                         MARK_CATEGORY_SUBTITLE, unname(kruskal$statistic),
                         unname(kruskal$parameter), kruskal$p.value),
      caption = sprintf(paste("y axis clipped to the %.1fth and %.1fth",
                              "percentiles; all %s genes enter the violin",
                              "shape, the medians and the tests."),
                        100 * VIOLIN_CLIP_PROBABILITIES[1],
                        100 * VIOLIN_CLIP_PROBABILITIES[2],
                        format(nrow(df), big.mark = ",")),
      x = "", y = value_label
    ) +
    theme_emseq() +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"),
          plot.caption = element_text(size = 8, color = "grey45"))

  save_multiformat_ggplot(figure,
                          file.path(out_dir, sprintf("%s_%s_by_mark", figure_id,
                                                     value_col)),
                          width = 11, height = 8)

  summary_out <- group_stats[, c("mark_category", "n", "median", "mean")]
  colnames(summary_out) <- c("mark_category", "n_genes", "median", "mean")
  summary_out$value_column <- value_col
  summary_out$kruskal_chisq <- unname(kruskal$statistic)
  summary_out$kruskal_df <- unname(kruskal$parameter)
  summary_out$kruskal_p <- kruskal$p.value
  write_section_table(summary_out,
                      file.path(out_dir, sprintf("%s_%s_summary.tsv", figure_id,
                                                 value_col)))
  write_section_table(tests$pairwise,
                      file.path(out_dir, sprintf("%s_%s_wilcoxon.tsv", figure_id,
                                                 value_col)))
  cat("\n")
  invisible(summary_out)
}

# =============================================================================
# REGISTERED FISHER TESTS
# =============================================================================

#' One registered gene-level Fisher test per mark category.
#'
#' The 2x2 table is the quadrant structure of the scatter: mCH gain against
#' MeCP2 gain. Both TRUE is Q1, both FALSE is Q3.
run_category_fisher_tests <- function(df, out_dir) {
  cat("--- Registered Fisher tests: mCH gain against MeCP2 gain per category ---\n")

  rows <- lapply(MARK_CATEGORY_ORDER, function(category) {
    sub <- df[df$mark_category == category, , drop = FALSE]
    gene_df <- data.frame(
      gene_name = sub$gene_name,
      chr = sub$chr,
      mch_gain = sub$mch_gain,
      mecp2_gain = sub$mecp2_gain,
      stringsAsFactors = FALSE
    )

    ft <- register_fisher_test(
      section = SECTION_ID,
      test_id = paste0("mch_gain_vs_mecp2_gain_",
                       MARK_CATEGORY_SLUGS[[category]]),
      description = sprintf(paste("Among genes in the %s mark category, gene-body",
                                  "mCH gain against MeCP2 binding gain."),
                            category),
      gene_df = gene_df, row_var = "mch_gain", col_var = "mecp2_gain",
      output_dir = out_dir
    )

    data.frame(
      mark_category = category,
      test_id = paste0("mch_gain_vs_mecp2_gain_",
                       MARK_CATEGORY_SLUGS[[category]]),
      n_genes = nrow(gene_df),
      n_mch_gain = sum(gene_df$mch_gain),
      n_mecp2_gain = sum(gene_df$mecp2_gain),
      n_both_gain = sum(gene_df$mch_gain & gene_df$mecp2_gain),
      n_both_loss = sum(!gene_df$mch_gain & !gene_df$mecp2_gain),
      odds_ratio = unname(ft$estimate),
      ci_lower = ft$conf.int[1],
      ci_upper = ft$conf.int[2],
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  })

  summary_tbl <- do.call(rbind, rows)
  write_section_table(summary_tbl,
                      file.path(out_dir, "20_04_fisher_by_mark_category.tsv"))
  cat("\n")
  summary_tbl
}

# =============================================================================
# PER-GENE OUTPUT
# =============================================================================

write_per_gene_table <- function(df, out_dir) {
  cols <- c("gene_name", "gene_id", "chr", "start", "end", "gene_length",
            "promoter_state", "body_state",
            "mark_category", "mch_ctrl", "mch_mut", "mch_diff", "mch_diff_pct",
            "edger_logFC", "edger_fdr", "mch_sig", "mch_direction",
            "mecp2_fold", "mecp2_fdr", "mecp2_min_fdr", "mecp2_n_peaks",
            "mecp2_n_sig", "mecp2_sig_gain", "mecp2_sig_loss",
            "mch_gain", "mecp2_gain", "quadrant")
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop("The analysis table is missing output columns: ",
         paste(missing, collapse = ", "))
  }
  write_section_table(df[, cols],
                      file.path(out_dir, "20_04_per_gene_mch_mecp2_mark.tsv"))
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opt <- parse_section_args()
  OUT_DIR <- opt$output_dir
  fdr_threshold <- opt$fdr_threshold

  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

  cat("================================================================================\n")
  cat("SECTION 20_04: mCH vs MeCP2 BY HISTONE MARK CONTEXT\n")
  cat("================================================================================\n")
  cat("Output dir:    ", OUT_DIR, "\n", sep = "")
  cat("FDR threshold: ", fdr_threshold, "\n\n", sep = "")

  cat("--- Loading the gene-level master table from section 20_02 ---\n")
  handoff <- load_gene_level_all_marks(HANDOFF_PATHS$gene_level_all_marks)
  genes <- deduplicate_by_gene_name(handoff)

  cat("\n--- Classifying genes by promoter state and gene-body state ---\n")
  genes <- classify_mark_category(genes)

  cat("\n--- Aggregating MeCP2 peaks to genes (nearest TSS) ---\n")
  mecp2_peaks <- prepare_mecp2_peaks(mecp2_diffbind, fdr_threshold)
  mecp2_gene <- aggregate_mecp2_to_genes(mecp2_peaks, fdr_threshold)

  cat("\n--- Building the analysis table ---\n")
  analysis_df <- build_analysis_table(genes, mecp2_gene, fdr_threshold)

  cat("\n--- Mark category group sizes ---\n")
  sizes <- report_category_sizes(genes, analysis_df, OUT_DIR)
  write_category_state_crosstabs(analysis_df, OUT_DIR)

  cat("\n--- Per-category Spearman correlations ---\n")
  cor_tbl <- correlation_by_category(analysis_df)

  cat("\n--- Comparing the correlations with the Fisher z transformation ---\n")
  comparison <- compare_correlations_fisher_z(cor_tbl)

  cat("\n--- Quadrant counts ---\n")
  quad_result <- quadrant_counts_by_category(analysis_df)

  cat("\n")
  plot_scatter_panels(analysis_df, cor_tbl, quad_result$counts, OUT_DIR)
  plot_quadrant_composition(analysis_df, quad_result, OUT_DIR)
  plot_correlation_comparison(cor_tbl, comparison, OUT_DIR)

  plot_violin_by_category(
    analysis_df, "mch_diff_pct",
    value_label = "mCH change (mutant - control, percentage points)",
    figure_id = "20_04d",
    plot_title = "Gene-body mCH change by histone mark context",
    out_dir = OUT_DIR, digits = 4
  )
  plot_violin_by_category(
    analysis_df, "mecp2_fold",
    value_label = "MeCP2 log2 fold change (mutant / control)",
    figure_id = "20_04e",
    plot_title = "MeCP2 binding change by histone mark context",
    out_dir = OUT_DIR, digits = 3
  )

  fisher_tbl <- run_category_fisher_tests(analysis_df, OUT_DIR)
  write_per_gene_table(analysis_df, OUT_DIR)

  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 20_04 SUMMARY\n")
  cat("================================================================================\n")
  cat(sprintf("Genes with mCH and MeCP2 data:  %s\n",
              format(nrow(analysis_df), big.mark = ",")))
  cat("Mark category sizes (genes with MeCP2 data):\n")
  for (i in seq_len(nrow(sizes))) {
    cat(sprintf("  %-12s %s of %s genes in the category\n",
                sizes$mark_category[i],
                format(sizes$n_genes_with_mecp2[i], big.mark = ","),
                format(sizes$n_genes_total[i], big.mark = ",")))
  }
  cat("Spearman rho (mCH change against MeCP2 fold):\n")
  for (i in seq_len(nrow(cor_tbl))) {
    cat(sprintf("  %-12s rho = %+.3f [%+.3f, %+.3f], p = %.3g (n = %s)\n",
                cor_tbl$mark_category[i], cor_tbl$spearman_rho[i],
                cor_tbl$rho_ci_lower[i], cor_tbl$rho_ci_upper[i],
                cor_tbl$spearman_p[i],
                format(cor_tbl$n_genes[i], big.mark = ",")))
  }
  cat(sprintf("Fisher z heterogeneity:         chi-squared = %.2f, df = %d, p = %.3g\n",
              comparison$overall$chi_squared, comparison$overall$df,
              comparison$overall$p_value))
  cat(sprintf("Quadrant composition chi-square: chi-squared = %.1f, df = %d, p = %.3g\n",
              unname(quad_result$chisq$statistic),
              unname(quad_result$chisq$parameter),
              quad_result$chisq$p.value))
  cat("Registered Fisher tests (mCH gain against MeCP2 gain):\n")
  for (i in seq_len(nrow(fisher_tbl))) {
    cat(sprintf("  %-12s OR = %.3f [%.3f, %.3f], p = %.3g (n = %s)\n",
                fisher_tbl$mark_category[i], fisher_tbl$odds_ratio[i],
                fisher_tbl$ci_lower[i], fisher_tbl$ci_upper[i],
                fisher_tbl$p_value[i],
                format(fisher_tbl$n_genes[i], big.mark = ",")))
  }
  cat("\nSection 20_04 complete.\n\n")
}

main()
