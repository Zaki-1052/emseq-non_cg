# scripts/sections/10_chromatin/10_02_ab_compartment.R
#
# Section 10_02: A/B compartment asymmetry of gene-body mCH change.
#
# What this tests
#   A DNMT3A-redistribution model predicts an asymmetric mCH response to BAP1
#   loss: gain of mCH in the A compartment (euchromatin) and loss of mCH in the
#   B compartment (heterochromatin). This section assigns every tested gene to
#   a Hi-C compartment bin and asks whether the mCH direction follows that
#   prediction.
#
# Reads
#   HIC_PATHS$compartments   HOMER getDiffExpression output. 25 kb bins, six
#                            per-sample PC1 columns (ctrl_M1..mut_M3), the
#                            "ctrl vs. mut Difference" column, and the
#                            "ctrl vs. mut adj. p-value" column.
#   mch_results              gene-level mCH differential table, already loaded
#                            by _shared_config.R.
#
# Writes (into OUT_DIR, default OUTPUT_PATHS$chromatin)
#   gene_compartment_assignment.tsv        one row per gene: bin, compartment,
#                                          shift, PC1 values, mCH statistics
#   compartment_bin_summary.tsv            bin counts per compartment and shift
#   compartment_gene_summary.tsv           mCH summary per compartment and shift
#   compartment_direction_proportions.tsv  numbers behind the stacked bar
#   compartment_fisher_tests.tsv           the four registered Fisher tests
#   compartment_wilcoxon_tests.tsv         compartment and shift Wilcoxon tests
#   compartment_pc1_spearman.tsv           PC1 against mch_diff correlation
#   violin_stats_compartment.tsv           n and median behind each violin
#   violin_stats_shift.tsv                 n and median behind each violin
#   shift_direction_observed_expected.tsv  shift x direction O/E ratios
#   figures: 10_02a_mch_violin_by_compartment, 10_02b_mch_violin_by_shift,
#            10_02c_mch_direction_stacked_bar, 10_02d_pc1_vs_mch_scatter,
#            10_02e_composite_compartment_summary
#
# Adapted from Biomodal section 29 (A/B compartment methylation mapping),
# written for the single mCH modality that EM-seq measures.
#
# Run:
#   Rscript scripts/sections/10_chromatin/10_02_ab_compartment.R 2>&1 | tee logs/10_02.txt

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)

# =============================================================================
# SECTION CONSTANTS
# =============================================================================

SECTION_ID <- "10_02"

# Order of the six per-sample PC1 columns in the HOMER header.
PC1_SAMPLE_ORDER <- c("ctrl_M1", "ctrl_M2", "ctrl_M3",
                      "mut_M1", "mut_M2", "mut_M3")

SHIFT_LEVELS <- c("B to A", "A to B", "Stable")

SHIFT_COLORS <- c("B to A" = "#D7191C",
                  "A to B" = "#2C7BB6",
                  "Stable" = "grey70")

COMPARTMENT_LEVELS <- c("A", "B")

# =============================================================================
# OPTIONS
# =============================================================================

parse_options <- function() {
  option_list <- list(
    make_option("--output-dir", type = "character", dest = "output_dir",
                default = OUTPUT_PATHS$chromatin,
                help = "Directory for figures and tables [default: %default]"),
    make_option("--shift-fdr", type = "double", dest = "shift_fdr",
                default = 0.05,
                help = "Adjusted p cutoff for a compartment shift [default: %default]"),
    make_option("--shift-diff", type = "double", dest = "shift_diff",
                default = 0.30,
                help = "Minimum |PC1 difference| for a compartment shift [default: %default]")
  )
  parse_args(OptionParser(option_list = option_list))
}

# =============================================================================
# SMALL HELPERS
# =============================================================================

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

write_tsv_table <- function(df, out_dir, filename) {
  write_section_table(df, file.path(out_dir, filename))
}

# =============================================================================
# STEP 1: LOAD HOMER COMPARTMENT BINS
# =============================================================================

#' Read the HOMER compartment table and classify every bin.
#'
#' Column names hold full command lines and file paths, so the reader keeps
#' them verbatim (check.names = FALSE) and finds the columns by pattern.
#' Compartment A is mean control PC1 > 0. The shift category comes from the
#' mutant-minus-control PC1 difference and its adjusted p-value.
load_compartment_bins <- function(path, shift_fdr, shift_diff) {
  if (!file.exists(path)) stop("Compartment file not found: ", path)

  raw <- read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE, comment.char = "", quote = "")
  cat(sprintf("  Loaded %s bins from %s\n",
              format(nrow(raw), big.mark = ","), basename(path)))

  pc1_cols <- grep("bedGraph avg over given bp", names(raw), value = TRUE)
  if (length(pc1_cols) != 6) {
    stop("Expected 6 per-sample PC1 columns matching ",
         "'bedGraph avg over given bp', found ", length(pc1_cols))
  }
  for (i in seq_along(PC1_SAMPLE_ORDER)) {
    if (!grepl(PC1_SAMPLE_ORDER[i], pc1_cols[i], fixed = TRUE)) {
      stop("PC1 column ", i, " does not name sample ", PC1_SAMPLE_ORDER[i],
           ". Column header: ", pc1_cols[i])
    }
  }

  difference_col <- grep("ctrl vs\\. mut Difference", names(raw), value = TRUE)
  if (length(difference_col) != 1) {
    stop("Expected exactly 1 'ctrl vs. mut Difference' column, found ",
         length(difference_col))
  }
  adj_pvalue_col <- grep("ctrl vs\\. mut adj\\. p-value", names(raw), value = TRUE)
  if (length(adj_pvalue_col) != 1) {
    stop("Expected exactly 1 'ctrl vs. mut adj. p-value' column, found ",
         length(adj_pvalue_col))
  }

  missing_coords <- setdiff(c("Chr", "Start", "End"), names(raw))
  if (length(missing_coords) > 0) {
    stop("Compartment file is missing coordinate columns: ",
         paste(missing_coords, collapse = ", "))
  }

  cat(sprintf("  Control PC1 columns: %s\n", paste(PC1_SAMPLE_ORDER[1:3], collapse = ", ")))
  cat(sprintf("  Mutant PC1 columns:  %s\n", paste(PC1_SAMPLE_ORDER[4:6], collapse = ", ")))
  cat(sprintf("  Difference column:   %s\n", difference_col))
  cat(sprintf("  Adjusted p column:   %s\n", adj_pvalue_col))

  ctrl_pc1 <- vapply(pc1_cols[1:3], function(cl) as.numeric(raw[[cl]]),
                     numeric(nrow(raw)))
  mut_pc1 <- vapply(pc1_cols[4:6], function(cl) as.numeric(raw[[cl]]),
                    numeric(nrow(raw)))

  bins <- data.frame(
    bin_chr = as.character(raw[["Chr"]]),
    bin_start = as.integer(raw[["Start"]]),
    bin_end = as.integer(raw[["End"]]),
    mean_ctrl_pc1 = rowMeans(ctrl_pc1),
    mean_mut_pc1 = rowMeans(mut_pc1),
    pc1_difference = as.numeric(raw[[difference_col]]),
    pc1_adj_pvalue = as.numeric(raw[[adj_pvalue_col]]),
    stringsAsFactors = FALSE
  )

  n_na <- sum(!complete.cases(bins))
  if (n_na > 0) {
    stop(n_na, " compartment bins hold NA in a coordinate, PC1, difference, ",
         "or adjusted p-value field.")
  }

  bins$compartment <- factor(ifelse(bins$mean_ctrl_pc1 > 0, "A", "B"),
                             levels = COMPARTMENT_LEVELS)

  shift <- rep("Stable", nrow(bins))
  passes <- bins$pc1_adj_pvalue < shift_fdr
  shift[passes & bins$pc1_difference > shift_diff] <- "B to A"
  shift[passes & bins$pc1_difference < -shift_diff] <- "A to B"
  bins$shift <- factor(shift, levels = SHIFT_LEVELS)

  n_before <- nrow(bins)
  bins <- bins[bins$bin_chr %in% CANONICAL_CHRS, , drop = FALSE]
  cat(sprintf("  Bins on canonical chromosomes: %s of %s\n",
              format(nrow(bins), big.mark = ","),
              format(n_before, big.mark = ",")))

  cat(sprintf("  Compartment A: %s bins (%.1f%%), B: %s bins (%.1f%%)\n",
              format(sum(bins$compartment == "A"), big.mark = ","),
              100 * mean(bins$compartment == "A"),
              format(sum(bins$compartment == "B"), big.mark = ","),
              100 * mean(bins$compartment == "B")))
  cat(sprintf("  Shift (adj p < %.2f, |difference| > %.2f): B to A %s, A to B %s, Stable %s\n",
              shift_fdr, shift_diff,
              format(sum(bins$shift == "B to A"), big.mark = ","),
              format(sum(bins$shift == "A to B"), big.mark = ","),
              format(sum(bins$shift == "Stable"), big.mark = ",")))

  bins
}

#' Convert compartment bins to GRanges.
#'
#' Each bin End equals the next bin Start, so the HOMER intervals are
#' half-open. Adding 1 to Start makes neighbouring bins disjoint. The
#' countOverlaps check below confirms that on the loaded data.
compartment_bins_to_granges <- function(bins) {
  gr <- GRanges(
    seqnames = factor(bins$bin_chr, levels = CANONICAL_CHRS),
    ranges = IRanges(start = bins$bin_start + 1L, end = bins$bin_end),
    bin_chr = bins$bin_chr,
    bin_start = bins$bin_start,
    bin_end = bins$bin_end,
    mean_ctrl_pc1 = bins$mean_ctrl_pc1,
    mean_mut_pc1 = bins$mean_mut_pc1,
    pc1_difference = bins$pc1_difference,
    pc1_adj_pvalue = bins$pc1_adj_pvalue,
    compartment = bins$compartment,
    shift = bins$shift
  )

  max_self_overlap <- max(countOverlaps(gr, gr))
  if (max_self_overlap != 1) {
    stop("Compartment bins overlap each other (max self-overlap = ",
         max_self_overlap, "). Gene assignment would be ambiguous.")
  }
  gr
}

# =============================================================================
# STEP 2: MAP GENES TO BINS
# =============================================================================

#' Keep one row per gene symbol, the one with the largest |edger_logFC|.
dedup_mch_by_gene <- function(mch) {
  n_dup <- sum(duplicated(mch$gene_name))
  out <- mch %>%
    dplyr::group_by(gene_name) %>%
    dplyr::slice_max(abs(edger_logFC), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    as.data.frame()
  cat(sprintf("  mCH genes: %s tested, %d duplicate symbols collapsed, %s kept\n",
              format(nrow(mch), big.mark = ","), n_dup,
              format(nrow(out), big.mark = ",")))
  out
}

#' Assign each gene to the compartment bin holding its midpoint.
#'
#' mch_results coordinates are 0-based half-open, so the 1-based midpoint is
#' floor((start + 1 + end) / 2). Genes whose midpoint falls in no bin drop out.
assign_genes_to_bins <- function(mch, bins_gr) {
  n_before <- nrow(mch)
  mch <- mch[mch$chr %in% CANONICAL_CHRS, , drop = FALSE]
  cat(sprintf("  Genes on canonical chromosomes: %s of %s\n",
              format(nrow(mch), big.mark = ","),
              format(n_before, big.mark = ",")))

  midpoint <- floor((mch$start + 1L + mch$end) / 2)
  gene_mid_gr <- GRanges(
    seqnames = factor(mch$chr, levels = CANONICAL_CHRS),
    ranges = IRanges(start = midpoint, width = 1)
  )

  hits <- findOverlaps(gene_mid_gr, bins_gr)
  if (any(duplicated(queryHits(hits)))) {
    stop("A gene midpoint overlaps more than one compartment bin.")
  }

  bin_meta <- as.data.frame(mcols(bins_gr))[subjectHits(hits), , drop = FALSE]
  merged <- cbind(mch[queryHits(hits), , drop = FALSE], bin_meta)
  rownames(merged) <- NULL

  match_rate <- nrow(merged) / nrow(mch)
  cat(sprintf("  Genes assigned to a bin: %s of %s (%.1f%%)\n",
              format(nrow(merged), big.mark = ","),
              format(nrow(mch), big.mark = ","), 100 * match_rate))
  if (match_rate <= 0.50) {
    stop("Only ", sprintf("%.1f%%", 100 * match_rate),
         " of genes fall inside a compartment bin. Check chromosome naming ",
         "and bin coverage.")
  }

  na_counts <- colSums(is.na(merged[, c("mch_diff", "mch_sig", "mch_hyper",
                                        "mch_hypo", "edger_logFC",
                                        "edger_fdr"), drop = FALSE]))
  if (any(na_counts > 0)) {
    stop("mCH columns hold NA after bin assignment: ",
         paste(sprintf("%s = %d", names(na_counts)[na_counts > 0],
                       na_counts[na_counts > 0]), collapse = ", "))
  }

  merged$compartment <- factor(as.character(merged$compartment),
                               levels = COMPARTMENT_LEVELS)
  merged$shift <- factor(as.character(merged$shift), levels = SHIFT_LEVELS)
  merged$in_compartment_a <- merged$compartment == "A"
  merged$in_compartment_b <- merged$compartment == "B"
  merged$shift_b_to_a <- merged$shift == "B to A"
  merged$shift_a_to_b <- merged$shift == "A to B"
  merged$sig_label <- ifelse(merged$mch_sig, "Significant", "Not Significant")

  cat(sprintf("  Assigned genes by compartment: A %s, B %s\n",
              format(sum(merged$compartment == "A"), big.mark = ","),
              format(sum(merged$compartment == "B"), big.mark = ",")))
  cat(sprintf("  Assigned genes by shift: B to A %s, A to B %s, Stable %s\n",
              format(sum(merged$shift == "B to A"), big.mark = ","),
              format(sum(merged$shift == "A to B"), big.mark = ","),
              format(sum(merged$shift == "Stable"), big.mark = ",")))
  merged
}

# =============================================================================
# STEP 3: FISHER TESTS
# =============================================================================

#' Extract the 2x2 counts and effect size of one registered Fisher test.
fisher_summary_row <- function(ft, test_id, description, gene_df,
                               row_var, col_var) {
  row_true <- gene_df[[row_var]]
  col_true <- gene_df[[col_var]]

  n_both <- sum(row_true & col_true)
  n_row_only <- sum(row_true & !col_true)
  n_col_only <- sum(!row_true & col_true)
  n_neither <- sum(!row_true & !col_true)

  data.frame(
    section = SECTION_ID,
    test_id = test_id,
    description = description,
    row_var = row_var,
    col_var = col_var,
    n_genes = nrow(gene_df),
    n_row_true_col_true = n_both,
    n_row_true_col_false = n_row_only,
    n_row_false_col_true = n_col_only,
    n_row_false_col_false = n_neither,
    pct_row_true_when_col_true = 100 * n_both / (n_both + n_col_only),
    pct_row_true_when_col_false = 100 * n_row_only / (n_row_only + n_neither),
    odds_ratio = unname(ft$estimate),
    ci_lower = ft$conf.int[1],
    ci_upper = ft$conf.int[2],
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

#' Run the four gene-level compartment Fisher tests through the registry.
run_compartment_fisher_tests <- function(gene_comp, out_dir) {
  specs <- list(
    list(test_id = "hyper_in_compartment_a",
         description = paste("mCH hypermethylated genes are enriched in the A",
                             "compartment (control PC1 > 0)."),
         row_var = "mch_hyper", col_var = "in_compartment_a"),
    list(test_id = "hypo_in_compartment_b",
         description = paste("mCH hypomethylated genes are enriched in the B",
                             "compartment (control PC1 <= 0)."),
         row_var = "mch_hypo", col_var = "in_compartment_b"),
    list(test_id = "hyper_in_shift_b_to_a",
         description = paste("mCH hypermethylated genes are enriched in bins",
                             "that shift from B to A in the mutant."),
         row_var = "mch_hyper", col_var = "shift_b_to_a"),
    list(test_id = "hypo_in_shift_a_to_b",
         description = paste("mCH hypomethylated genes are enriched in bins",
                             "that shift from A to B in the mutant."),
         row_var = "mch_hypo", col_var = "shift_a_to_b")
  )

  rows <- lapply(specs, function(spec) {
    gene_df <- gene_comp[, c("gene_name", "chr", spec$row_var, spec$col_var)]
    ft <- register_fisher_test(
      section = SECTION_ID, test_id = spec$test_id,
      description = spec$description,
      gene_df = gene_df, row_var = spec$row_var, col_var = spec$col_var,
      output_dir = out_dir)
    fisher_summary_row(ft, spec$test_id, spec$description, gene_df,
                       spec$row_var, spec$col_var)
  })

  do.call(rbind, rows)
}

# =============================================================================
# STEP 4: WILCOXON AND SPEARMAN TESTS
# =============================================================================

#' Wilcoxon of mch_diff between the A and B compartments.
wilcoxon_by_compartment <- function(gene_comp) {
  wt <- wilcox.test(mch_diff ~ compartment, data = gene_comp)
  data.frame(
    comparison = "A vs B (compartment)",
    group1 = "A", group2 = "B",
    n_group1 = sum(gene_comp$compartment == "A"),
    n_group2 = sum(gene_comp$compartment == "B"),
    median_group1 = median(gene_comp$mch_diff[gene_comp$compartment == "A"]),
    median_group2 = median(gene_comp$mch_diff[gene_comp$compartment == "B"]),
    statistic = unname(wt$statistic),
    p_value = wt$p.value,
    p_adjusted = wt$p.value,
    stringsAsFactors = FALSE
  )
}

#' Pairwise Wilcoxon of mch_diff between the three shift categories.
wilcoxon_by_shift <- function(gene_comp) {
  group_sizes <- table(gene_comp$shift)
  if (any(group_sizes < 3)) {
    stop("Shift categories with fewer than 3 genes: ",
         paste(names(group_sizes)[group_sizes < 3], collapse = ", "))
  }

  pw <- pairwise.wilcox.test(gene_comp$mch_diff, gene_comp$shift,
                             p.adjust.method = "BH")
  pw_long <- as.data.frame(as.table(pw$p.value), stringsAsFactors = FALSE)
  names(pw_long) <- c("group1", "group2", "p_adjusted")
  pw_long <- pw_long[!is.na(pw_long$p_adjusted), , drop = FALSE]

  medians <- tapply(gene_comp$mch_diff, gene_comp$shift, median)
  counts <- as.integer(group_sizes)
  names(counts) <- names(group_sizes)

  data.frame(
    comparison = paste0(pw_long$group1, " vs ", pw_long$group2, " (shift)"),
    group1 = pw_long$group1,
    group2 = pw_long$group2,
    n_group1 = counts[pw_long$group1],
    n_group2 = counts[pw_long$group2],
    median_group1 = medians[pw_long$group1],
    median_group2 = medians[pw_long$group2],
    statistic = NA_real_,
    p_value = NA_real_,
    p_adjusted = pw_long$p_adjusted,
    stringsAsFactors = FALSE
  )
}

#' Spearman correlation of control PC1 against mch_diff.
spearman_pc1_vs_mch <- function(gene_comp) {
  ct <- cor.test(gene_comp$mean_ctrl_pc1, gene_comp$mch_diff,
                 method = "spearman", exact = FALSE)
  data.frame(
    x = "mean_ctrl_pc1",
    y = "mch_diff",
    n_genes = nrow(gene_comp),
    rho = unname(ct$estimate),
    statistic_S = unname(ct$statistic),
    p_value = ct$p.value,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# STEP 5: DIRECTION PROPORTIONS AND O/E
# =============================================================================

#' Count mCH direction within each compartment and each shift category.
#'
#' Significant genes only, so the proportions describe the genes that changed.
direction_proportions <- function(gene_comp) {
  sig <- gene_comp[gene_comp$mch_sig, , drop = FALSE]

  by_compartment <- sig %>%
    dplyr::count(group_label = compartment, mch_direction) %>%
    dplyr::group_by(group_label) %>%
    dplyr::mutate(pct = 100 * n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(grouping = "By compartment")

  by_shift <- sig %>%
    dplyr::count(group_label = shift, mch_direction) %>%
    dplyr::group_by(group_label) %>%
    dplyr::mutate(pct = 100 * n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(grouping = "By compartment shift")

  out <- dplyr::bind_rows(
    dplyr::mutate(by_compartment, group_label = as.character(group_label)),
    dplyr::mutate(by_shift, group_label = as.character(group_label))
  )
  out$group_label <- factor(out$group_label,
                            levels = c(COMPARTMENT_LEVELS, SHIFT_LEVELS))
  out$grouping <- factor(out$grouping,
                         levels = c("By compartment", "By compartment shift"))
  as.data.frame(out)
}

#' Observed over expected counts for shift category by mCH direction.
shift_direction_observed_expected <- function(gene_comp) {
  sig <- gene_comp[gene_comp$mch_sig, , drop = FALSE]
  observed <- table(sig$shift, sig$mch_direction)
  expected <- outer(rowSums(observed), colSums(observed)) / sum(observed)

  out <- as.data.frame(as.table(observed), stringsAsFactors = FALSE)
  names(out) <- c("shift", "mch_direction", "observed")
  out$expected <- as.vector(expected)
  out$oe_ratio <- out$observed / out$expected
  out
}

#' mCH summary per compartment and per shift category.
compartment_gene_summary <- function(gene_comp) {
  summarise_one <- function(df, group_col, grouping_name) {
    df %>%
      dplyr::group_by(group_label = .data[[group_col]]) %>%
      dplyr::summarise(
        n_genes = dplyr::n(),
        mean_mch_diff = mean(mch_diff),
        median_mch_diff = median(mch_diff),
        mean_edger_logFC = mean(edger_logFC),
        n_sig = sum(mch_sig),
        n_hyper = sum(mch_hyper),
        n_hypo = sum(mch_hypo),
        pct_sig = 100 * sum(mch_sig) / dplyr::n(),
        pct_hyper_of_sig = 100 * sum(mch_hyper) / sum(mch_sig),
        median_mean_ctrl_pc1 = median(mean_ctrl_pc1),
        .groups = "drop"
      ) %>%
      dplyr::mutate(grouping = grouping_name,
                    group_label = as.character(group_label))
  }

  out <- dplyr::bind_rows(
    summarise_one(gene_comp, "compartment", "compartment"),
    summarise_one(gene_comp, "shift", "shift")
  )
  dplyr::select(as.data.frame(out), grouping, group_label, dplyr::everything())
}

#' Bin counts per compartment and shift category.
compartment_bin_summary <- function(bins) {
  counts <- as.data.frame(table(compartment = bins$compartment,
                                shift = bins$shift),
                          stringsAsFactors = FALSE)
  names(counts)[3] <- "n_bins"
  counts$pct_of_bins <- 100 * counts$n_bins / nrow(bins)
  counts
}

# =============================================================================
# STEP 6: FIGURES
# =============================================================================

#' n and median of mch_diff per group, with the grouping column as an ordered
#' factor so the violin panels keep the level order of the section.
violin_stats <- function(df, group_col, levels_order) {
  stats <- summarise_groups(df, group_col, "mch_diff")
  stats[[group_col]] <- factor(as.character(stats[[group_col]]),
                               levels = levels_order)
  stats
}

#' Violin of mch_diff across one grouping, annotated with n and median.
plot_mch_violin <- function(df, group_col, stats, fill_colors,
                            title, subtitle, x_lab) {
  y_range <- unname(quantile(df$mch_diff, c(0.005, 0.995), na.rm = TRUE))
  pad <- diff(y_range) * 0.12
  stats$label_y <- y_range[2] + pad
  stats$label <- group_label(stats, digits = 4)

  ggplot(df, aes(x = .data[[group_col]], y = mch_diff,
                 fill = .data[[group_col]])) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_text(data = stats,
              aes(x = .data[[group_col]], y = label_y, label = label),
              inherit.aes = FALSE, size = 3.4, lineheight = 0.95) +
    coord_cartesian(ylim = c(y_range[1] - pad, y_range[2] + 3 * pad)) +
    scale_fill_manual(values = fill_colors) +
    labs(title = title, subtitle = subtitle, x = x_lab,
         y = "mCH difference (mutant - control)") +
    theme_emseq() +
    theme(legend.position = "none")
}

#' Stacked proportion bar of mCH direction per compartment and per shift.
plot_direction_bar <- function(proportions, n_sig) {
  ggplot(proportions, aes(x = group_label, y = pct, fill = mch_direction)) +
    geom_col(position = "stack", width = 0.7) +
    geom_text(aes(label = sprintf("n = %s\n%.1f%%",
                                  format(n, big.mark = ","), pct)),
              position = position_stack(vjust = 0.5), size = 3,
              lineheight = 0.95) +
    facet_wrap(~grouping, scales = "free_x") +
    scale_fill_manual(values = COLORS$direction) +
    labs(
      title = "mCH direction by compartment and compartment shift",
      subtitle = sprintf("Significant genes only (FDR < %.2f), n = %s",
                         Q_THRESHOLD, format(n_sig, big.mark = ",")),
      x = NULL, y = "Percent of significant genes", fill = "mCH direction"
    ) +
    theme_emseq() +
    theme(legend.position = "top")
}

#' Scatter of control PC1 against mch_diff with a loess fit and KEY_GENES.
plot_pc1_scatter <- function(gene_comp, spearman) {
  key_genes <- gene_comp[gene_comp$gene_name %in% KEY_GENES, , drop = FALSE]
  y_range <- unname(quantile(gene_comp$mch_diff, c(0.005, 0.995), na.rm = TRUE))
  x_range <- range(gene_comp$mean_ctrl_pc1)

  ggplot(gene_comp, aes(x = mean_ctrl_pc1, y = mch_diff)) +
    geom_point(aes(color = sig_label), alpha = 0.20, size = 0.5) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                color = "#333333", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggrepel::geom_text_repel(
      data = key_genes, aes(label = gene_name),
      size = 3, min.segment.length = 0, max.overlaps = Inf,
      box.padding = 0.4, segment.color = "grey30") +
    annotate("text", x = x_range[2] * 0.55, y = y_range[2] * 0.90,
             label = sprintf("Spearman rho = %.3f\n%s",
                             spearman$rho, fmt_p(spearman$p_value)),
             hjust = 0, size = 4, fontface = "italic") +
    annotate("text", x = x_range[2] * 0.75, y = y_range[1] * 0.85,
             label = "A compartment", color = COLORS$compartment[["A"]],
             fontface = "bold", size = 4) +
    annotate("text", x = x_range[1] * 0.75, y = y_range[1] * 0.85,
             label = "B compartment", color = COLORS$compartment[["B"]],
             fontface = "bold", size = 4) +
    coord_cartesian(ylim = y_range) +
    scale_color_manual(values = COLORS$significant) +
    labs(
      title = "Control PC1 against mCH change",
      subtitle = sprintf("Spearman rho = %.3f, %s | n = %s genes",
                         spearman$rho, fmt_p(spearman$p_value),
                         format(spearman$n_genes, big.mark = ",")),
      x = "Mean control PC1 (A/B compartment score)",
      y = "mCH difference (mutant - control)",
      color = NULL
    ) +
    theme_emseq() +
    theme(legend.position = "top")
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opt <- parse_options()
  OUT_DIR <- opt$output_dir
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 10_02: A/B COMPARTMENT ASYMMETRY OF mCH CHANGE\n")
  cat("================================================================================\n")
  cat("Output dir:  ", OUT_DIR, "\n", sep = "")
  cat("Shift FDR:   ", opt$shift_fdr, "\n", sep = "")
  cat("Shift |diff|:", opt$shift_diff, "\n", sep = "")
  cat("\n")

  cat("STEP 1: Load HOMER compartment bins\n")
  bins <- load_compartment_bins(HIC_PATHS$compartments,
                                opt$shift_fdr, opt$shift_diff)
  bins_gr <- compartment_bins_to_granges(bins)
  cat("\n")

  cat("STEP 2: Assign genes to compartment bins\n")
  mch_dedup <- dedup_mch_by_gene(mch_results)
  gene_comp <- assign_genes_to_bins(mch_dedup, bins_gr)
  cat("\n")

  cat("STEP 3: Gene-level Fisher tests\n")
  fisher_table <- run_compartment_fisher_tests(gene_comp, OUT_DIR)
  cat("\n")

  cat("STEP 4: Wilcoxon and Spearman tests\n")
  wilcoxon_table <- rbind(wilcoxon_by_compartment(gene_comp),
                          wilcoxon_by_shift(gene_comp))
  rownames(wilcoxon_table) <- NULL
  spearman_table <- spearman_pc1_vs_mch(gene_comp)
  for (i in seq_len(nrow(wilcoxon_table))) {
    cat(sprintf("  Wilcoxon %s: median %s = %.5f, median %s = %.5f, adjusted %s\n",
                wilcoxon_table$comparison[i],
                wilcoxon_table$group1[i], wilcoxon_table$median_group1[i],
                wilcoxon_table$group2[i], wilcoxon_table$median_group2[i],
                fmt_p(wilcoxon_table$p_adjusted[i])))
  }
  cat(sprintf("  Spearman PC1 vs mCH: rho = %.3f, %s\n",
              spearman_table$rho, fmt_p(spearman_table$p_value)))
  cat("\n")

  cat("STEP 5: Direction proportions\n")
  proportions <- direction_proportions(gene_comp)
  oe_table <- shift_direction_observed_expected(gene_comp)
  gene_summary <- compartment_gene_summary(gene_comp)
  bin_summary <- compartment_bin_summary(bins)
  for (i in seq_len(nrow(gene_summary))) {
    cat(sprintf("  %s %s: n = %s, %s significant, %.1f%% hypermethylated of significant\n",
                gene_summary$grouping[i], gene_summary$group_label[i],
                format(gene_summary$n_genes[i], big.mark = ","),
                format(gene_summary$n_sig[i], big.mark = ","),
                gene_summary$pct_hyper_of_sig[i]))
  }
  cat("\n")

  cat("STEP 6: Figures\n")
  stats_compartment <- violin_stats(gene_comp, "compartment", COMPARTMENT_LEVELS)
  stats_shift <- violin_stats(gene_comp, "shift", SHIFT_LEVELS)
  wilcox_compartment_p <- wilcoxon_table$p_adjusted[1]

  p_compartment <- plot_mch_violin(
    gene_comp, "compartment", stats_compartment, COLORS$compartment,
    title = "mCH change by A/B compartment",
    subtitle = sprintf("Wilcoxon A vs B: %s", fmt_p(wilcox_compartment_p)),
    x_lab = "Compartment (control PC1)")
  save_multiformat_ggplot(p_compartment,
                          file.path(OUT_DIR, "10_02a_mch_violin_by_compartment"),
                          width = 8, height = 7)

  shift_subtitle <- paste(
    sprintf("%s vs %s: %s",
            wilcoxon_table$group1[-1], wilcoxon_table$group2[-1],
            vapply(wilcoxon_table$p_adjusted[-1], fmt_p, character(1))),
    collapse = " | ")
  p_shift <- plot_mch_violin(
    gene_comp, "shift", stats_shift, SHIFT_COLORS,
    title = "mCH change by compartment shift",
    subtitle = sprintf("Pairwise Wilcoxon, BH adjusted: %s", shift_subtitle),
    x_lab = sprintf("Compartment shift (adj p < %.2f, |difference| > %.2f)",
                    opt$shift_fdr, opt$shift_diff))
  save_multiformat_ggplot(p_shift,
                          file.path(OUT_DIR, "10_02b_mch_violin_by_shift"),
                          width = 10, height = 7)

  p_bar <- plot_direction_bar(proportions, sum(gene_comp$mch_sig))
  save_multiformat_ggplot(p_bar,
                          file.path(OUT_DIR, "10_02c_mch_direction_stacked_bar"),
                          width = 11, height = 7)

  p_scatter <- plot_pc1_scatter(gene_comp, spearman_table)
  save_multiformat_ggplot(p_scatter,
                          file.path(OUT_DIR, "10_02d_pc1_vs_mch_scatter"),
                          width = 10, height = 9)

  p_composite <- (p_compartment + p_shift + p_bar) +
    plot_layout(widths = c(1, 1.3, 1.6)) +
    plot_annotation(
      title = "A/B compartment asymmetry of gene-body mCH change in BAP1-KO",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = 14))
    )
  save_multiformat_ggplot(p_composite,
                          file.path(OUT_DIR, "10_02e_composite_compartment_summary"),
                          width = 20, height = 7)
  cat("\n")

  cat("STEP 7: Tables\n")
  gene_export <- gene_comp[, c(
    "gene_name", "gene_id", "chr", "start", "end", "strand", "gene_length",
    "bin_chr", "bin_start", "bin_end",
    "mean_ctrl_pc1", "mean_mut_pc1", "pc1_difference", "pc1_adj_pvalue",
    "compartment", "shift",
    "mch_ctrl", "mch_mut", "mch_diff", "edger_logFC", "edger_fdr",
    "mch_sig", "mch_direction", "mch_hyper", "mch_hypo")]
  write_tsv_table(gene_export, OUT_DIR, "gene_compartment_assignment.tsv")
  write_tsv_table(bin_summary, OUT_DIR, "compartment_bin_summary.tsv")
  write_tsv_table(gene_summary, OUT_DIR, "compartment_gene_summary.tsv")
  write_tsv_table(proportions, OUT_DIR, "compartment_direction_proportions.tsv")
  write_tsv_table(oe_table, OUT_DIR, "shift_direction_observed_expected.tsv")
  write_tsv_table(fisher_table, OUT_DIR, "compartment_fisher_tests.tsv")
  write_tsv_table(wilcoxon_table, OUT_DIR, "compartment_wilcoxon_tests.tsv")
  write_tsv_table(spearman_table, OUT_DIR, "compartment_pc1_spearman.tsv")
  write_tsv_table(stats_compartment, OUT_DIR, "violin_stats_compartment.tsv")
  write_tsv_table(stats_shift, OUT_DIR, "violin_stats_shift.tsv")
  cat("\n")

  cat("================================================================================\n")
  cat("SECTION 10_02 SUMMARY\n")
  cat("================================================================================\n")
  cat(sprintf("Genes assigned to compartment bins: %s\n",
              format(nrow(gene_comp), big.mark = ",")))
  for (i in seq_len(nrow(fisher_table))) {
    cat(sprintf("  Fisher %s: OR = %.3f (95%% CI %.3f-%.3f), %s\n",
                fisher_table$test_id[i], fisher_table$odds_ratio[i],
                fisher_table$ci_lower[i], fisher_table$ci_upper[i],
                fmt_p(fisher_table$p_value[i])))
  }
  cat(sprintf("  Median mCH difference: A = %.5f, B = %.5f\n",
              stats_compartment$median[stats_compartment$compartment == "A"],
              stats_compartment$median[stats_compartment$compartment == "B"]))
  cat(sprintf("  Spearman control PC1 against mCH difference: rho = %.3f, %s\n",
              spearman_table$rho, fmt_p(spearman_table$p_value)))
  cat(sprintf("\nFigures and tables written to: %s\n", OUT_DIR))
  cat("Section 10_02 complete.\n\n")
}

main()
