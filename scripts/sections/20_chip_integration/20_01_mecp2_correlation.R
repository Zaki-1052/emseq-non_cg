# scripts/sections/20_chip_integration/20_01_mecp2_correlation.R
#
# Section 20_01: mCH direction against MeCP2 binding change.
#
# What this tests
#   MeCP2 reads methylated cytosines. If BAP1 loss changes gene-body mCH, and
#   MeCP2 occupancy tracks that change, then genes that gain mCH must gain
#   MeCP2 and genes that lose mCH must lose MeCP2.
#
# Analyses
#   1. Gene-body overlap of significant mCH genes with MeCP2 up and down peaks.
#   2. Two registered gene-level Fisher tests: mCH hyper against MeCP2 gained,
#      mCH hypo against MeCP2 lost.
#   3. Spearman correlation of mch_diff against gene-level MeCP2 fold.
#   4. Wilcoxon tests of MeCP2 fold across the three mCH direction groups.
#   5. Linear model of MeCP2 fold on mch_diff with log10 gene length.
#   6. Logistic model of MeCP2 gain on mch_diff with log10 gene length, and a
#      predicted-probability curve.
#   7. Cross-check of the two gene-level aggregation routes for MeCP2 fold.
#
# Reads
#   mch_results, gene_bodies, mecp2_diffbind  (shared config)
#   MECP2_PATHS$annotated, MECP2_PATHS$up, MECP2_PATHS$down
#
# Writes
#   Figures and TSV tables into OUTPUT_PATHS$chip (override with --output-dir).
#   One row per registered Fisher test into HANDOFF_PATHS$fisher_registry.
#
# Adapted from Biomodal section 11, MeCP2 correlation at DMRs.

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)
library(broom)

SECTION_ID <- "20_01"

MCH_CATEGORY_LEVELS <- c("Hypermethylated", "Hypomethylated", "Not Significant")

MCH_CATEGORY_COLORS <- c(
  "Hypermethylated" = unname(COLORS$direction["Hypermethylated"]),
  "Hypomethylated"  = unname(COLORS$direction["Hypomethylated"]),
  "Not Significant" = "grey70"
)

# =============================================================================
# OPTIONS
# =============================================================================

parse_options <- function() {
  option_list <- list(
    make_option("--output-dir", type = "character", default = OUTPUT_PATHS$chip,
                dest = "output_dir",
                help = "Directory for figures and tables [default: %default]"),
    make_option("--fdr-threshold", type = "double", default = Q_THRESHOLD,
                dest = "fdr_threshold",
                help = "FDR cutoff for MeCP2 peak significance [default: %default]")
  )
  parse_args(OptionParser(option_list = option_list))
}

# =============================================================================
# INPUT LOADING
# =============================================================================

#' Read the MeCP2 peak-level annotation table.
#'
#' Stops when a required column is absent or when Fold, FDR or distanceToTSS
#' holds a missing value. Peaks without a gene symbol are dropped and counted.
#'
#' @param filepath Path to MeCP2_annotated.txt.
#' @return data.frame of annotated peaks carrying a gene symbol
load_mecp2_annotated <- function(filepath) {
  if (!file.exists(filepath)) {
    stop("MeCP2 annotated peak file not found: ", filepath)
  }

  df <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "", fill = TRUE)

  required <- c("seqnames", "start", "end", "Fold", "p.value", "FDR",
                "geneStart", "geneEnd", "geneLength", "distanceToTSS", "SYMBOL")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("MeCP2 annotated file is missing columns: ",
         paste(missing, collapse = ", "), " (", filepath, ")")
  }

  numeric_cols <- c("start", "end", "Fold", "p.value", "FDR",
                    "geneStart", "geneEnd", "geneLength", "distanceToTSS")
  for (col in numeric_cols) df[[col]] <- as.numeric(df[[col]])

  for (col in c("Fold", "FDR", "distanceToTSS")) {
    n_na <- sum(is.na(df[[col]]))
    if (n_na > 0) {
      stop("MeCP2 annotated file has ", n_na, " missing values in column ", col)
    }
  }

  n_all <- nrow(df)
  df <- df[!is.na(df$SYMBOL) & nzchar(df$SYMBOL), , drop = FALSE]

  cat(sprintf("  MeCP2 annotated peaks: %s total, %s with a gene symbol\n",
              format(n_all, big.mark = ","), format(nrow(df), big.mark = ",")))
  cat(sprintf("  Peak-level significance at FDR<%.2f: %s up, %s down\n",
              Q_THRESHOLD,
              format(sum(df$FDR < Q_THRESHOLD & df$Fold > 0), big.mark = ","),
              format(sum(df$FDR < Q_THRESHOLD & df$Fold < 0), big.mark = ",")))
  df
}

# =============================================================================
# GENE-LEVEL AGGREGATION
# =============================================================================

#' Row indices of mch_results that keep one row per gene symbol.
#'
#' Some gene names carry more than one ENSMUSG identifier. The row with the
#' largest absolute edgeR log fold change is kept for each name.
#'
#' @param mch mch_results data.frame.
#' @return integer vector of row indices in ascending order
deduplicate_mch_row_indices <- function(mch) {
  ord <- order(mch$gene_name, -abs(mch$edger_logFC))
  keep <- ord[!duplicated(mch$gene_name[ord])]
  sort(keep)
}

#' Collapse annotated MeCP2 peaks to one fold change per gene.
#'
#' Uses the shared nearest-TSS collapse rule and adds the distance of that
#' peak's gene to its TSS plus the annotated gene length.
#'
#' @param annotated data.frame from load_mecp2_annotated().
#' @param fdr_threshold FDR cutoff for the per-gene significant peak count.
#' @return data.frame with one row per gene
aggregate_annotated_by_gene <- function(annotated, fdr_threshold) {
  gene_tbl <- aggregate_diffbind_by_gene(annotated, method = "nearest_tss",
                                         fdr_threshold = fdr_threshold,
                                         prefix = "mecp2")

  extras <- annotated %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(
      mecp2_min_abs_dist_tss = min(abs(distanceToTSS)),
      mecp2_annot_gene_length = geneLength[which.min(abs(distanceToTSS))],
      .groups = "drop"
    ) %>%
    dplyr::rename(gene_name = SYMBOL) %>%
    as.data.frame()

  out <- dplyr::left_join(gene_tbl, extras, by = "gene_name")
  cat(sprintf("  Genes with a MeCP2 annotated peak: %s\n",
              format(nrow(out), big.mark = ",")))
  out
}

#' Aggregate the filtered MeCP2 DiffBind peakset to genes with ChIPseeker.
#'
#' @param diffbind mecp2_diffbind from the shared config.
#' @param fdr_threshold FDR cutoff passed to the collapse rule.
#' @return data.frame with one row per gene, prefix mecp2db
aggregate_diffbind_route <- function(diffbind, fdr_threshold) {
  annotated <- annotate_peaks_to_genes(diffbind, "MeCP2")
  out <- aggregate_diffbind_by_gene(annotated, method = "nearest_tss",
                                    fdr_threshold = fdr_threshold,
                                    prefix = "mecp2db")
  cat(sprintf("  Genes from the DiffBind route: %s\n",
              format(nrow(out), big.mark = ",")))
  out
}

#' Compare the two gene-level MeCP2 fold estimates.
#'
#' Writes the per-gene comparison and a one-row summary, then draws a scatter
#' of one route against the other.
#'
#' @param annot_gene data.frame from aggregate_annotated_by_gene().
#' @param db_gene data.frame from aggregate_diffbind_route().
#' @param out_dir Section output directory.
#' @return data.frame of the shared genes
crosscheck_aggregation_routes <- function(annot_gene, db_gene, out_dir) {
  shared <- dplyr::inner_join(
    annot_gene[, c("gene_name", "mecp2_fold", "mecp2_fdr", "mecp2_n_peaks")],
    db_gene[, c("gene_name", "mecp2db_fold", "mecp2db_fdr", "mecp2db_n_peaks")],
    by = "gene_name"
  )

  if (nrow(shared) < 100) {
    stop("Only ", nrow(shared), " genes are shared by the two MeCP2 ",
         "aggregation routes. Expected the annotated and DiffBind peaksets ",
         "to cover a common gene set.")
  }

  rho_test <- cor.test(shared$mecp2_fold, shared$mecp2db_fold, method = "spearman")
  pearson_test <- cor.test(shared$mecp2_fold, shared$mecp2db_fold, method = "pearson")
  same_sign <- sign(shared$mecp2_fold) == sign(shared$mecp2db_fold)
  pct_same_sign <- 100 * mean(same_sign)

  cat(sprintf("  Cross-check on %s shared genes: Spearman rho = %.3f (p = %.3g)\n",
              format(nrow(shared), big.mark = ","),
              unname(rho_test$estimate), rho_test$p.value))
  cat(sprintf("  Sign agreement between routes: %.1f%%\n", pct_same_sign))

  write_section_table(shared,
                      file.path(out_dir, "mecp2_aggregation_crosscheck_genes.tsv"))

  summary_tbl <- data.frame(
    n_genes_annotated_route = nrow(annot_gene),
    n_genes_diffbind_route = nrow(db_gene),
    n_genes_shared = nrow(shared),
    spearman_rho = unname(rho_test$estimate),
    spearman_p = rho_test$p.value,
    pearson_r = unname(pearson_test$estimate),
    pearson_p = pearson_test$p.value,
    pct_same_sign = pct_same_sign,
    median_fold_annotated = median(shared$mecp2_fold),
    median_fold_diffbind = median(shared$mecp2db_fold),
    stringsAsFactors = FALSE
  )
  write_section_table(summary_tbl,
                      file.path(out_dir, "mecp2_aggregation_crosscheck_summary.tsv"))

  p <- ggplot(shared, aes(x = mecp2_fold, y = mecp2db_fold)) +
    geom_point(alpha = 0.25, size = 1.2, color = "grey30") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "#D95F02", linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
    labs(
      title = "MeCP2 Gene-Level Fold Change: Two Aggregation Routes",
      subtitle = sprintf(
        "Spearman ρ = %.3f, p = %.2e | sign agreement %.1f%% | n = %s genes",
        unname(rho_test$estimate), rho_test$p.value, pct_same_sign,
        format(nrow(shared), big.mark = ",")),
      x = "MeCP2 log2 fold change (annotated peakset, nearest TSS)",
      y = "MeCP2 log2 fold change (DiffBind peakset, nearest TSS)"
    ) +
    theme_emseq()

  save_multiformat_ggplot(p, file.path(out_dir, "20_01f_aggregation_crosscheck"),
                          width = 9, height = 8)
  shared
}

# =============================================================================
# GENE TABLE
# =============================================================================

#' Build the one-row-per-gene analysis table.
#'
#' Joins the deduplicated mCH results, the gene-body overlap with MeCP2 up and
#' down peaks, and the gene-level MeCP2 fold change.
#'
#' @param mch Deduplicated mch_results rows.
#' @param bodies GRanges of the same rows, in the same order.
#' @param annot_gene Gene-level MeCP2 table from the annotated peakset.
#' @param mecp2_up_gr GRanges of MeCP2 up peaks.
#' @param mecp2_down_gr GRanges of MeCP2 down peaks.
#' @param fdr_threshold FDR cutoff defining a significant MeCP2 gain or loss.
#' @return data.frame with one row per gene
build_gene_table <- function(mch, bodies, annot_gene,
                             mecp2_up_gr, mecp2_down_gr, fdr_threshold) {
  stopifnot(length(bodies) == nrow(mch))

  tbl <- data.frame(
    gene_name = mch$gene_name,
    gene_id = mch$gene_id,
    chr = mch$chr,
    start = mch$start,
    end = mch$end,
    gene_length = mch$gene_length,
    mch_ctrl = mch$mch_ctrl,
    mch_mut = mch$mch_mut,
    mch_diff = mch$mch_diff,
    edger_logFC = mch$edger_logFC,
    edger_fdr = mch$edger_fdr,
    mch_sig = mch$mch_sig,
    mch_hyper = mch$mch_hyper,
    mch_hypo = mch$mch_hypo,
    stringsAsFactors = FALSE
  )

  tbl$mch_diff_pct <- 100 * tbl$mch_diff
  tbl$log10_gene_length <- log10(tbl$gene_length)

  tbl$mecp2_up_overlap <- countOverlaps(bodies, mecp2_up_gr) > 0
  tbl$mecp2_down_overlap <- countOverlaps(bodies, mecp2_down_gr) > 0

  tbl$mch_category <- factor(
    ifelse(tbl$mch_hyper, "Hypermethylated",
      ifelse(tbl$mch_hypo, "Hypomethylated", "Not Significant")),
    levels = MCH_CATEGORY_LEVELS
  )

  tbl <- dplyr::left_join(tbl, annot_gene, by = "gene_name")

  tbl$mecp2_gain <- tbl$mecp2_fdr < fdr_threshold & tbl$mecp2_fold > 0
  tbl$mecp2_loss <- tbl$mecp2_fdr < fdr_threshold & tbl$mecp2_fold < 0
  tbl$has_mecp2 <- !is.na(tbl$mecp2_fold)

  cat(sprintf("  Gene table: %s genes, %s with a gene-level MeCP2 fold change\n",
              format(nrow(tbl), big.mark = ","),
              format(sum(tbl$has_mecp2), big.mark = ",")))
  tbl
}

# =============================================================================
# ANALYSIS 1: GENE-BODY OVERLAP WITH MeCP2 UP AND DOWN PEAKS
# =============================================================================

#' Count and plot how often each mCH group overlaps MeCP2 up and down peaks.
#'
#' @param gene_tbl Gene table from build_gene_table().
#' @param out_dir Section output directory.
#' @return data.frame of counts and percentages
run_overlap_analysis <- function(gene_tbl, out_dir) {
  overlap_df <- gene_tbl %>%
    dplyr::group_by(mch_category) %>%
    dplyr::summarise(
      n_genes = dplyr::n(),
      n_mecp2_up = sum(mecp2_up_overlap),
      n_mecp2_down = sum(mecp2_down_overlap),
      .groups = "drop"
    ) %>%
    as.data.frame()

  long_df <- overlap_df %>%
    tidyr::pivot_longer(cols = c("n_mecp2_up", "n_mecp2_down"),
                        names_to = "mecp2_direction", values_to = "count") %>%
    dplyr::mutate(
      mecp2_direction = ifelse(mecp2_direction == "n_mecp2_up",
                               "MeCP2 Up", "MeCP2 Down"),
      mecp2_direction = factor(mecp2_direction,
                               levels = c("MeCP2 Up", "MeCP2 Down")),
      percentage = 100 * count / n_genes
    ) %>%
    as.data.frame()

  for (i in seq_len(nrow(long_df))) {
    cat(sprintf("  %-16s x %-11s: %s/%s genes (%.1f%%)\n",
                as.character(long_df$mch_category[i]),
                as.character(long_df$mecp2_direction[i]),
                format(long_df$count[i], big.mark = ","),
                format(long_df$n_genes[i], big.mark = ","),
                long_df$percentage[i]))
  }

  write_section_table(long_df,
                      file.path(out_dir, "mecp2_overlap_by_mch_direction.tsv"))

  p <- ggplot(long_df, aes(x = mch_category, y = percentage, fill = mecp2_direction)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.75),
             width = 0.65, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%\n(%s/%s)", percentage,
                                  format(count, big.mark = ",", trim = TRUE),
                                  format(n_genes, big.mark = ",", trim = TRUE))),
              position = position_dodge(width = 0.75), vjust = -0.25, size = 3) +
    scale_fill_manual(values = COLORS$mecp2[c("MeCP2 Up", "MeCP2 Down")],
                      name = "MeCP2 peak set") +
    scale_y_continuous(limits = c(0, max(long_df$percentage) * 1.35),
                       expand = c(0, 0)) +
    labs(
      title = "MeCP2 Peak Overlap at mCH Gene Bodies",
      subtitle = "Percentage of gene bodies in each mCH group that overlap a differential MeCP2 peak",
      x = "mCH direction (edgeR FDR < 0.05)",
      y = "% of gene bodies overlapping a MeCP2 peak"
    ) +
    theme_emseq() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p, file.path(out_dir, "20_01a_mecp2_overlap_by_mch_direction"),
                          width = 9, height = 7)
  long_df
}

#' Register the two directional Fisher tests and write their summary.
#'
#' @param gene_tbl Gene table from build_gene_table().
#' @param out_dir Section output directory.
#' @return data.frame with one row per test
run_direction_fisher_tests <- function(gene_tbl, out_dir) {
  ft_hyper <- register_fisher_test(
    section = SECTION_ID, test_id = "hyper_vs_mecp2_up",
    description = paste("Do gene bodies that gain mCH in the mutant overlap",
                        "MeCP2 up peaks more often than other genes?"),
    gene_df = gene_tbl, row_var = "mch_hyper", col_var = "mecp2_up_overlap",
    output_dir = out_dir)

  ft_hypo <- register_fisher_test(
    section = SECTION_ID, test_id = "hypo_vs_mecp2_down",
    description = paste("Do gene bodies that lose mCH in the mutant overlap",
                        "MeCP2 down peaks more often than other genes?"),
    gene_df = gene_tbl, row_var = "mch_hypo", col_var = "mecp2_down_overlap",
    output_dir = out_dir)

  summary_tbl <- data.frame(
    test_id = c("hyper_vs_mecp2_up", "hypo_vs_mecp2_down"),
    row_var = c("mch_hyper", "mch_hypo"),
    col_var = c("mecp2_up_overlap", "mecp2_down_overlap"),
    n_genes = nrow(gene_tbl),
    n_row_true = c(sum(gene_tbl$mch_hyper), sum(gene_tbl$mch_hypo)),
    n_col_true = c(sum(gene_tbl$mecp2_up_overlap), sum(gene_tbl$mecp2_down_overlap)),
    n_both_true = c(sum(gene_tbl$mch_hyper & gene_tbl$mecp2_up_overlap),
                    sum(gene_tbl$mch_hypo & gene_tbl$mecp2_down_overlap)),
    odds_ratio = c(unname(ft_hyper$estimate), unname(ft_hypo$estimate)),
    ci_low = c(ft_hyper$conf.int[1], ft_hypo$conf.int[1]),
    ci_high = c(ft_hyper$conf.int[2], ft_hypo$conf.int[2]),
    p_value = c(ft_hyper$p.value, ft_hypo$p.value),
    stringsAsFactors = FALSE
  )
  write_section_table(summary_tbl,
                      file.path(out_dir, "mecp2_direction_fisher_summary.tsv"))
  summary_tbl
}

# =============================================================================
# ANALYSIS 2: MeCP2 FOLD CHANGE BY mCH DIRECTION
# =============================================================================

#' Violin of MeCP2 fold change split by mCH direction, with Wilcoxon tests.
#'
#' Each group carries its gene count and median on the figure.
#'
#' @param gene_tbl Gene table from build_gene_table().
#' @param out_dir Section output directory.
#' @return data.frame of the pairwise Wilcoxon results
run_fold_by_direction <- function(gene_tbl, out_dir) {
  df <- gene_tbl[gene_tbl$has_mecp2, , drop = FALSE]

  grp_summary <- summarise_groups(df, "mch_category", "mecp2_fold")
  write_section_table(grp_summary,
                      file.path(out_dir, "mecp2_fold_by_mch_direction_summary.tsv"))

  for (i in seq_len(nrow(grp_summary))) {
    cat(sprintf("  %-16s n = %s, median MeCP2 fold = %.4f, mean = %.4f\n",
                as.character(grp_summary$mch_category[i]),
                format(grp_summary$n[i], big.mark = ","),
                grp_summary$median[i], grp_summary$mean[i]))
  }

  fold_of <- function(level) df$mecp2_fold[df$mch_category == level]
  pairs <- list(
    c("Hypermethylated", "Hypomethylated"),
    c("Hypermethylated", "Not Significant"),
    c("Hypomethylated", "Not Significant")
  )

  wilcox_tbl <- do.call(rbind, lapply(pairs, function(pr) {
    a <- fold_of(pr[1])
    b <- fold_of(pr[2])
    wt <- wilcox.test(a, b)
    data.frame(
      group_1 = pr[1], group_2 = pr[2],
      n_1 = length(a), n_2 = length(b),
      median_1 = median(a), median_2 = median(b),
      median_difference = median(a) - median(b),
      W = unname(wt$statistic), p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  }))
  write_section_table(wilcox_tbl,
                      file.path(out_dir, "mecp2_fold_wilcoxon_tests.tsv"))

  for (i in seq_len(nrow(wilcox_tbl))) {
    cat(sprintf("  Wilcoxon %s vs %s: W = %.0f, p = %.3g\n",
                wilcox_tbl$group_1[i], wilcox_tbl$group_2[i],
                wilcox_tbl$W[i], wilcox_tbl$p_value[i]))
  }

  y_min <- min(df$mecp2_fold)
  y_max <- max(df$mecp2_fold)
  y_span <- y_max - y_min
  y_label <- y_max + 0.10 * y_span
  y_stats <- y_max + 0.30 * y_span

  wilcox_label <- paste(
    sprintf("Wilcoxon %s vs %s: p = %.2e",
            wilcox_tbl$group_1, wilcox_tbl$group_2, wilcox_tbl$p_value),
    collapse = "\n"
  )

  grp_summary$label <- group_label(grp_summary)

  p <- ggplot(df, aes(x = mch_category, y = mecp2_fold, fill = mch_category)) +
    geom_violin(alpha = 0.6, scale = "width",
                draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black",
               linewidth = 0.4) +
    geom_text(data = grp_summary,
              aes(x = mch_category, y = y_label, label = label),
              inherit.aes = FALSE, size = 3.2, lineheight = 1.1) +
    annotate("text", x = 2, y = y_stats, label = wilcox_label,
             size = 3.2, fontface = "italic", lineheight = 1.2) +
    scale_fill_manual(values = MCH_CATEGORY_COLORS) +
    scale_y_continuous(limits = c(y_min - 0.05 * y_span, y_stats + 0.12 * y_span)) +
    labs(
      title = "MeCP2 Fold Change by mCH Direction",
      subtitle = "Gene-level MeCP2 fold change from the peak nearest the TSS",
      x = "mCH direction (edgeR FDR < 0.05)",
      y = "MeCP2 log2 fold change (mutant / control)"
    ) +
    theme_emseq() +
    theme(legend.position = "none")

  save_multiformat_ggplot(p, file.path(out_dir, "20_01b_mecp2_fold_by_mch_direction"),
                          width = 9, height = 8)
  wilcox_tbl
}

# =============================================================================
# ANALYSIS 3: mCH CHANGE AGAINST MeCP2 FOLD CHANGE
# =============================================================================

#' Spearman correlation and scatter of mch_diff against MeCP2 fold change.
#'
#' Computes the correlation over all genes with MeCP2 data and again over the
#' mCH-significant subset.
#'
#' @param gene_tbl Gene table from build_gene_table().
#' @param out_dir Section output directory.
#' @return data.frame with one row per gene set
run_scatter_correlation <- function(gene_tbl, out_dir) {
  df <- gene_tbl[gene_tbl$has_mecp2, , drop = FALSE]
  sig_df <- df[df$mch_sig, , drop = FALSE]

  correlate <- function(d, label) {
    st <- suppressWarnings(cor.test(d$mch_diff_pct, d$mecp2_fold, method = "spearman"))
    pt <- cor.test(d$mch_diff_pct, d$mecp2_fold, method = "pearson")
    data.frame(
      gene_set = label,
      n_genes = nrow(d),
      spearman_rho = unname(st$estimate), spearman_p = st$p.value,
      pearson_r = unname(pt$estimate), pearson_p = pt$p.value,
      stringsAsFactors = FALSE
    )
  }

  cor_tbl <- rbind(
    correlate(df, "all_genes_with_mecp2"),
    correlate(sig_df, "mch_significant")
  )
  write_section_table(cor_tbl,
                      file.path(out_dir, "mecp2_mch_spearman_correlation.tsv"))

  for (i in seq_len(nrow(cor_tbl))) {
    cat(sprintf("  %s (n = %s): Spearman rho = %.3f, p = %.3g\n",
                cor_tbl$gene_set[i], format(cor_tbl$n_genes[i], big.mark = ","),
                cor_tbl$spearman_rho[i], cor_tbl$spearman_p[i]))
  }

  df$quadrant <- assign_quadrant(df$mch_diff_pct, df$mecp2_fold)
  quad_tbl <- df %>%
    dplyr::count(quadrant, name = "n_genes") %>%
    dplyr::mutate(percentage = 100 * n_genes / sum(n_genes)) %>%
    as.data.frame()
  write_section_table(quad_tbl,
                      file.path(out_dir, "mecp2_mch_quadrant_counts.tsv"))

  df$label_gene <- ifelse(df$gene_name %in% KEY_GENES, df$gene_name, "")

  x_max <- max(df$mch_diff_pct)
  x_min <- min(df$mch_diff_pct)
  y_max <- max(df$mecp2_fold)
  y_min <- min(df$mecp2_fold)
  quad_n <- setNames(quad_tbl$n_genes, quad_tbl$quadrant)
  quad_count <- function(q) if (q %in% names(quad_n)) quad_n[[q]] else 0L

  p <- ggplot(df, aes(x = mch_diff_pct, y = mecp2_fold)) +
    geom_point(aes(color = mch_category), alpha = 0.35, size = 1.3) +
    geom_smooth(method = "lm", formula = y ~ x, color = "black",
                linewidth = 0.8, se = TRUE, alpha = 0.15) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    geom_text_repel(aes(label = label_gene), size = 3, max.overlaps = 20,
                    fontface = "italic", color = "grey15",
                    segment.color = "grey60", segment.size = 0.3,
                    min.segment.length = 0) +
    annotate("text", x = 0.75 * x_max, y = 0.92 * y_max,
             label = sprintf("mCH↑ MeCP2↑\nn = %s",
                             format(quad_count("Q1"), big.mark = ",")),
             size = 3.2, color = "#D95F02", fontface = "bold") +
    annotate("text", x = 0.75 * x_min, y = 0.92 * y_max,
             label = sprintf("mCH↓ MeCP2↑\nn = %s",
                             format(quad_count("Q2"), big.mark = ",")),
             size = 3.2, color = "grey40") +
    annotate("text", x = 0.75 * x_min, y = 0.92 * y_min,
             label = sprintf("mCH↓ MeCP2↓\nn = %s",
                             format(quad_count("Q3"), big.mark = ",")),
             size = 3.2, color = "#7570B3", fontface = "bold") +
    annotate("text", x = 0.75 * x_max, y = 0.92 * y_min,
             label = sprintf("mCH↑ MeCP2↓\nn = %s",
                             format(quad_count("Q4"), big.mark = ",")),
             size = 3.2, color = "grey40") +
    scale_color_manual(values = MCH_CATEGORY_COLORS, name = "mCH direction") +
    labs(
      title = "mCH Change vs MeCP2 Fold Change",
      subtitle = sprintf(
        "All genes with MeCP2 data: Spearman ρ = %.3f, p = %.2e, n = %s | mCH-significant subset: ρ = %.3f, p = %.2e, n = %s",
        cor_tbl$spearman_rho[1], cor_tbl$spearman_p[1],
        format(cor_tbl$n_genes[1], big.mark = ","),
        cor_tbl$spearman_rho[2], cor_tbl$spearman_p[2],
        format(cor_tbl$n_genes[2], big.mark = ",")),
      x = "mCH change (mutant - control, percentage points)",
      y = "MeCP2 log2 fold change (mutant / control)"
    ) +
    theme_emseq() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p, file.path(out_dir, "20_01c_mch_diff_vs_mecp2_fold"),
                          width = 11, height = 9)
  cor_tbl
}

# =============================================================================
# ANALYSIS 4: REGRESSION MODELS
# =============================================================================

#' Fit MeCP2 fold change on mCH change with gene length as a covariate.
#'
#' Fits the model twice: over every gene with MeCP2 data, and over the
#' mCH-significant subset. Gene length enters as log10 base pairs.
#'
#' @param gene_tbl Gene table from build_gene_table().
#' @param out_dir Section output directory.
#' @return list with the all-gene fit and the tidy coefficient table
fit_linear_models <- function(gene_tbl, out_dir) {
  df <- gene_tbl[gene_tbl$has_mecp2, , drop = FALSE]
  sig_df <- df[df$mch_sig, , drop = FALSE]

  if (nrow(sig_df) < 50) {
    stop("Only ", nrow(sig_df), " mCH-significant genes carry MeCP2 data. ",
         "The linear model needs at least 50.")
  }

  fit_one <- function(d, scope) {
    fit <- lm(mecp2_fold ~ mch_diff_pct + log10_gene_length, data = d)
    tidy_tbl <- broom::tidy(fit, conf.int = TRUE)
    glance_tbl <- broom::glance(fit)
    tidy_tbl$model_scope <- scope
    tidy_tbl$n_genes <- nrow(d)
    tidy_tbl$r_squared <- glance_tbl$r.squared
    tidy_tbl$adj_r_squared <- glance_tbl$adj.r.squared
    tidy_tbl$model_f_statistic <- glance_tbl$statistic
    tidy_tbl$model_p_value <- glance_tbl$p.value
    list(fit = fit, tidy = as.data.frame(tidy_tbl), glance = glance_tbl)
  }

  all_fit <- fit_one(df, "all_genes_with_mecp2")
  sig_fit <- fit_one(sig_df, "mch_significant")

  coef_tbl <- rbind(all_fit$tidy, sig_fit$tidy)
  write_section_table(coef_tbl,
                      file.path(out_dir, "mecp2_linear_model_coefficients.tsv"))

  cat(sprintf("  Linear model (all genes, n = %s): R2 = %.4f, F p = %.3g\n",
              format(nrow(df), big.mark = ","),
              all_fit$glance$r.squared, all_fit$glance$p.value))
  for (i in seq_len(nrow(all_fit$tidy))) {
    cat(sprintf("    %-20s beta = %+.4f (SE = %.4f, p = %.3g)\n",
                all_fit$tidy$term[i], all_fit$tidy$estimate[i],
                all_fit$tidy$std.error[i], all_fit$tidy$p.value[i]))
  }
  cat(sprintf("  Linear model (mCH significant, n = %s): R2 = %.4f, F p = %.3g\n",
              format(nrow(sig_df), big.mark = ","),
              sig_fit$glance$r.squared, sig_fit$glance$p.value))

  list(fit = all_fit$fit, coefficients = coef_tbl, glance = all_fit$glance,
       n_all = nrow(df), n_sig = nrow(sig_df))
}

#' Fit the probability of a significant MeCP2 gain on mCH change.
#'
#' Fits over every gene with MeCP2 data and again over the mCH-significant
#' subset. Gene length enters as log10 base pairs.
#'
#' @param gene_tbl Gene table from build_gene_table().
#' @param out_dir Section output directory.
#' @return list with the all-gene fit and the odds-ratio table
fit_logistic_models <- function(gene_tbl, out_dir) {
  df <- gene_tbl[gene_tbl$has_mecp2, , drop = FALSE]
  df$mecp2_gain_int <- as.integer(df$mecp2_gain)

  n_gain <- sum(df$mecp2_gain_int == 1)
  n_other <- sum(df$mecp2_gain_int == 0)
  cat(sprintf("  MeCP2 significant gain: %s genes, other: %s genes\n",
              format(n_gain, big.mark = ","), format(n_other, big.mark = ",")))
  if (n_gain < 20 || n_other < 20) {
    stop("The logistic model needs at least 20 genes in each class. Found ",
         n_gain, " with a MeCP2 gain and ", n_other, " without.")
  }

  sig_df <- df[df$mch_sig, , drop = FALSE]

  fit_one <- function(d, scope) {
    fit <- glm(mecp2_gain_int ~ mch_diff_pct + log10_gene_length,
               data = d, family = binomial)
    tidy_tbl <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
    tidy_raw <- broom::tidy(fit, conf.int = FALSE, exponentiate = FALSE)
    glance_tbl <- broom::glance(fit)
    tidy_tbl$log_odds_estimate <- tidy_raw$estimate
    tidy_tbl$model_scope <- scope
    tidy_tbl$n_genes <- nrow(d)
    tidy_tbl$n_events <- sum(d$mecp2_gain_int)
    tidy_tbl$aic <- glance_tbl$AIC
    tidy_tbl$null_deviance <- glance_tbl$null.deviance
    tidy_tbl$residual_deviance <- glance_tbl$deviance
    list(fit = fit, tidy = as.data.frame(tidy_tbl), glance = glance_tbl)
  }

  all_fit <- fit_one(df, "all_genes_with_mecp2")
  sig_fit <- fit_one(sig_df, "mch_significant")

  or_tbl <- rbind(all_fit$tidy, sig_fit$tidy)
  write_section_table(or_tbl,
                      file.path(out_dir, "mecp2_logistic_model_odds_ratios.tsv"))

  cat(sprintf("  Logistic model (all genes, n = %s): AIC = %.1f\n",
              format(nrow(df), big.mark = ","), all_fit$glance$AIC))
  for (i in seq_len(nrow(all_fit$tidy))) {
    cat(sprintf("    %-20s OR = %.4f (95%% CI %.4f to %.4f, p = %.3g)\n",
                all_fit$tidy$term[i], all_fit$tidy$estimate[i],
                all_fit$tidy$conf.low[i], all_fit$tidy$conf.high[i],
                all_fit$tidy$p.value[i]))
  }

  list(fit = all_fit$fit, odds_ratios = or_tbl, glance = all_fit$glance,
       data = df)
}

#' Forest plot of the linear and logistic model terms.
#'
#' @param lm_res Result of fit_linear_models().
#' @param glm_res Result of fit_logistic_models().
#' @param out_dir Section output directory.
#' @return the patchwork object
plot_model_coefficients <- function(lm_res, glm_res, out_dir) {
  term_labels <- c(
    "mch_diff_pct" = "mCH change\n(percentage points)",
    "log10_gene_length" = "Gene length\n(log10 bp)"
  )

  lm_df <- lm_res$coefficients %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
      term_label = unname(term_labels[term]),
      significance = ifelse(p.value < 0.05, "p < 0.05", "p >= 0.05")
    )

  glm_df <- glm_res$odds_ratios %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
      term_label = unname(term_labels[term]),
      significance = ifelse(p.value < 0.05, "p < 0.05", "p >= 0.05")
    )

  sig_colors <- c("p < 0.05" = "#E41A1C", "p >= 0.05" = "grey50")

  p_lm <- ggplot(lm_df, aes(x = estimate, y = term_label, color = significance)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.5) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high, shape = model_scope,
                        group = model_scope),
                    size = 0.7, linewidth = 0.8,
                    position = position_dodge(width = 0.5)) +
    scale_color_manual(values = sig_colors, name = "Significance") +
    scale_shape_manual(values = c("all_genes_with_mecp2" = 16,
                                  "mch_significant" = 17),
                       name = "Gene set") +
    labs(
      title = "Linear Model: MeCP2 Fold Change",
      subtitle = sprintf("mecp2_fold ~ mch_diff_pct + log10_gene_length | R² = %.4f (all genes, n = %s)",
                         lm_res$glance$r.squared,
                         format(lm_res$n_all, big.mark = ",")),
      x = "Coefficient estimate (95% CI)", y = ""
    ) +
    theme_emseq() +
    theme(legend.position = "bottom", legend.box = "vertical")

  p_glm <- ggplot(glm_df, aes(x = estimate, y = term_label, color = significance)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50",
               linewidth = 0.5) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high, shape = model_scope,
                        group = model_scope),
                    size = 0.7, linewidth = 0.8,
                    position = position_dodge(width = 0.5)) +
    scale_color_manual(values = sig_colors, name = "Significance") +
    scale_shape_manual(values = c("all_genes_with_mecp2" = 16,
                                  "mch_significant" = 17),
                       name = "Gene set") +
    scale_x_log10() +
    labs(
      title = "Logistic Model: Probability of MeCP2 Gain",
      subtitle = sprintf("mecp2_gain ~ mch_diff_pct + log10_gene_length | AIC = %.1f (all genes)",
                         glm_res$glance$AIC),
      x = "Odds ratio (95% CI, log scale)", y = ""
    ) +
    theme_emseq() +
    theme(legend.position = "bottom", legend.box = "vertical")

  combined <- p_lm + p_glm +
    plot_annotation(
      title = "Does mCH Change Predict MeCP2 Binding Change?",
      subtitle = "Gene length enters both models as log10 base pairs",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(combined,
                          file.path(out_dir, "20_01d_mecp2_model_coefficients"),
                          width = 14, height = 7)
  combined
}

#' Predicted probability of a MeCP2 gain across the range of mCH change.
#'
#' Holds log10 gene length at its median and draws a 95% confidence band from
#' the standard error on the link scale.
#'
#' @param glm_res Result of fit_logistic_models().
#' @param out_dir Section output directory.
#' @return data.frame of the prediction grid
plot_predicted_probability <- function(glm_res, out_dir) {
  df <- glm_res$data
  median_length <- median(df$log10_gene_length)

  grid <- data.frame(
    mch_diff_pct = seq(min(df$mch_diff_pct), max(df$mch_diff_pct),
                       length.out = 400),
    log10_gene_length = median_length
  )

  link <- predict(glm_res$fit, newdata = grid, type = "link", se.fit = TRUE)
  grid$probability <- plogis(link$fit)
  grid$ci_low <- plogis(link$fit - 1.96 * link$se.fit)
  grid$ci_high <- plogis(link$fit + 1.96 * link$se.fit)
  grid$gene_length_bp <- 10^median_length

  write_section_table(grid,
                      file.path(out_dir, "mecp2_gain_predicted_probability.tsv"))

  or_tbl <- glm_res$odds_ratios
  mch_row <- or_tbl[or_tbl$term == "mch_diff_pct" &
                      or_tbl$model_scope == "all_genes_with_mecp2", ]

  rug_df <- df
  rug_df$gain_label <- ifelse(rug_df$mecp2_gain, "MeCP2 gain", "No MeCP2 gain")

  p <- ggplot() +
    geom_ribbon(data = grid,
                aes(x = mch_diff_pct, ymin = ci_low, ymax = ci_high),
                fill = "#D95F02", alpha = 0.20) +
    geom_line(data = grid, aes(x = mch_diff_pct, y = probability),
              color = "#D95F02", linewidth = 1) +
    geom_rug(data = rug_df,
             aes(x = mch_diff_pct, color = gain_label),
             sides = "b", alpha = 0.25, length = unit(0.02, "npc")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    scale_color_manual(values = c("MeCP2 gain" = "#D95F02",
                                  "No MeCP2 gain" = "grey60"),
                       name = "Observed genes") +
    labs(
      title = "Predicted Probability of a Significant MeCP2 Gain",
      subtitle = sprintf(
        "Gene length held at the median (%s bp) | OR per percentage point of mCH = %.3f (95%% CI %.3f to %.3f, p = %.2e) | n = %s genes",
        format(round(10^median_length), big.mark = ","),
        mch_row$estimate[1], mch_row$conf.low[1], mch_row$conf.high[1],
        mch_row$p.value[1], format(nrow(df), big.mark = ",")),
      x = "mCH change (mutant - control, percentage points)",
      y = "P(MeCP2 significantly gained)"
    ) +
    theme_emseq() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p,
                          file.path(out_dir, "20_01e_mecp2_gain_predicted_probability"),
                          width = 10, height = 7)
  grid
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opt <- parse_options()
  OUT_DIR <- opt$output_dir
  fdr_threshold <- opt$fdr_threshold

  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

  cat("================================================================================\n")
  cat("SECTION 20_01: mCH DIRECTION vs MeCP2 BINDING CHANGE\n")
  cat("================================================================================\n")
  cat("Output dir:    ", OUT_DIR, "\n", sep = "")
  cat("FDR threshold: ", fdr_threshold, "\n\n", sep = "")

  stopifnot(
    "MeCP2 annotated peak file not found" = file.exists(MECP2_PATHS$annotated),
    "MeCP2 up peak BED not found" = file.exists(MECP2_PATHS$up),
    "MeCP2 down peak BED not found" = file.exists(MECP2_PATHS$down)
  )

  cat("--- Loading MeCP2 peak data ---\n")
  mecp2_annotated <- load_mecp2_annotated(MECP2_PATHS$annotated)
  mecp2_up_gr <- load_chip_peaks(MECP2_PATHS$up, "MeCP2 up")
  mecp2_down_gr <- load_chip_peaks(MECP2_PATHS$down, "MeCP2 down")

  cat("\n--- Aggregating MeCP2 peaks to genes ---\n")
  annot_gene <- aggregate_annotated_by_gene(mecp2_annotated, fdr_threshold)
  db_gene <- aggregate_diffbind_route(mecp2_diffbind, fdr_threshold)

  cat("\n--- Cross-checking the two aggregation routes ---\n")
  crosscheck_aggregation_routes(annot_gene, db_gene, OUT_DIR)

  cat("\n--- Building the gene table ---\n")
  keep_idx <- deduplicate_mch_row_indices(mch_results)
  cat(sprintf("  Deduplicated %s mCH rows to %s gene names\n",
              format(nrow(mch_results), big.mark = ","),
              format(length(keep_idx), big.mark = ",")))
  gene_tbl <- build_gene_table(mch_results[keep_idx, , drop = FALSE],
                               gene_bodies[keep_idx],
                               annot_gene, mecp2_up_gr, mecp2_down_gr,
                               fdr_threshold)

  gene_tbl <- dplyr::left_join(
    gene_tbl,
    db_gene[, c("gene_name", "mecp2db_fold", "mecp2db_fdr", "mecp2db_n_peaks")],
    by = "gene_name"
  )
  write_section_table(gene_tbl, file.path(OUT_DIR, "mecp2_mch_gene_level.tsv"))

  cat("\n--- Analysis 1: gene-body overlap with MeCP2 peaks ---\n")
  run_overlap_analysis(gene_tbl, OUT_DIR)

  cat("\n--- Analysis 2: directional Fisher tests ---\n")
  fisher_tbl <- run_direction_fisher_tests(gene_tbl, OUT_DIR)

  cat("\n--- Analysis 3: MeCP2 fold change by mCH direction ---\n")
  wilcox_tbl <- run_fold_by_direction(gene_tbl, OUT_DIR)

  cat("\n--- Analysis 4: mCH change vs MeCP2 fold change ---\n")
  cor_tbl <- run_scatter_correlation(gene_tbl, OUT_DIR)

  cat("\n--- Analysis 5: linear model ---\n")
  lm_res <- fit_linear_models(gene_tbl, OUT_DIR)

  cat("\n--- Analysis 6: logistic model ---\n")
  glm_res <- fit_logistic_models(gene_tbl, OUT_DIR)

  cat("\n--- Model figures ---\n")
  plot_model_coefficients(lm_res, glm_res, OUT_DIR)
  plot_predicted_probability(glm_res, OUT_DIR)

  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 20_01 SUMMARY\n")
  cat("================================================================================\n")
  cat(sprintf("Genes analysed:                 %s\n",
              format(nrow(gene_tbl), big.mark = ",")))
  cat(sprintf("Genes with MeCP2 fold change:   %s\n",
              format(sum(gene_tbl$has_mecp2), big.mark = ",")))
  cat(sprintf("mCH hypermethylated:            %s\n",
              format(sum(gene_tbl$mch_hyper), big.mark = ",")))
  cat(sprintf("mCH hypomethylated:             %s\n",
              format(sum(gene_tbl$mch_hypo), big.mark = ",")))
  for (i in seq_len(nrow(fisher_tbl))) {
    cat(sprintf("Fisher %-20s OR = %.3f (95%% CI %.3f to %.3f), p = %.3g\n",
                fisher_tbl$test_id[i], fisher_tbl$odds_ratio[i],
                fisher_tbl$ci_low[i], fisher_tbl$ci_high[i],
                fisher_tbl$p_value[i]))
  }
  cat(sprintf("Spearman (all genes):           rho = %.3f, p = %.3g\n",
              cor_tbl$spearman_rho[1], cor_tbl$spearman_p[1]))
  cat(sprintf("Spearman (mCH significant):     rho = %.3f, p = %.3g\n",
              cor_tbl$spearman_rho[2], cor_tbl$spearman_p[2]))
  cat(sprintf("Wilcoxon hyper vs hypo:         p = %.3g\n",
              wilcox_tbl$p_value[wilcox_tbl$group_1 == "Hypermethylated" &
                                   wilcox_tbl$group_2 == "Hypomethylated"]))
  lm_mch <- lm_res$coefficients[
    lm_res$coefficients$term == "mch_diff_pct" &
      lm_res$coefficients$model_scope == "all_genes_with_mecp2", ]
  cat(sprintf("Linear model mch_diff_pct:      beta = %+.4f, p = %.3g\n",
              lm_mch$estimate[1], lm_mch$p.value[1]))
  glm_mch <- glm_res$odds_ratios[
    glm_res$odds_ratios$term == "mch_diff_pct" &
      glm_res$odds_ratios$model_scope == "all_genes_with_mecp2", ]
  cat(sprintf("Logistic model mch_diff_pct:    OR = %.4f, p = %.3g\n",
              glm_mch$estimate[1], glm_mch$p.value[1]))
  cat("\nSection 20_01 complete.\n\n")
}

main()
