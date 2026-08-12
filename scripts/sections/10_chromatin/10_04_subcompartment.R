# scripts/sections/10_chromatin/10_04_subcompartment.R
#
# Section 10_04: gene-body mCH change stratified by CALDER2 Hi-C subcompartment.
#
# What this tests
#   Whether differential mCH concentrates in particular chromatin
#   subcompartments. CALDER2 splits the genome into A.1 (strong active),
#   A.2 (weak active), B.1 (facultative heterochromatin) and B.2 (constitutive
#   heterochromatin) at 100 kb resolution. Each tested gene is placed in one bin
#   by its midpoint. The section then measures, per subcompartment:
#     - the fraction of genes with a significant mCH change,
#     - the hyper/hypo direction split among significant genes,
#     - the control mCH level and the mCH difference distributions,
#     - the H3K27me3 and H3K27ac peak overlap rate.
#   It also compares mCH difference between bins whose subcompartment label
#   changed in the mutant and bins whose label stayed the same.
#
# Reads
#   HIC_PATHS$subcompartments   CALDER2 labels, 100 kb bins, 1-based bin_start
#   CHIP_PATHS$h3k27me3         H3K27me3 consensus peaks (BED)
#   CHIP_PATHS$h3k27ac          H3K27ac consensus peaks (BED)
#   mch_results                 loaded by _shared_config.R
#
# Writes (into OUTPUT_PATHS$chromatin, override with --output-dir)
#   Figures  10_04a_significance_rate_by_subcompartment
#            10_04b_direction_split_by_subcompartment
#            10_04c_mch_level_by_subcompartment
#            10_04d_mch_diff_by_subcompartment
#            10_04e_histone_overlap_by_subcompartment
#            10_04_subcompartment_composite
#   Tables   10_04_*.tsv, one per computed statistic
#   Fisher   fisher_tables/10_04_*_genes.tsv plus rows in the shared registry
#
# Adapted from Biomodal section 66 (section_66_subcompartment_methylation.R).
# EM-seq measures mCH as a single modality, so the 5hmC violin of Biomodal
# panel 66c is dropped and the 5mC violin becomes the mCH level violin.

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)

# =============================================================================
# SECTION CONSTANTS
# =============================================================================

SECTION_ID <- "10_04"

# Axis and facet text for each CALDER2 label.
SUBCMPT_LABELS <- c(
  "A.1" = "A.1\n(strong active)",
  "A.2" = "A.2\n(weak active)",
  "B.1" = "B.1\n(facultative het)",
  "B.2" = "B.2\n(constitutive het)"
)

# Histone bar colours reuse the "gained" entry of each mark palette.
HISTONE_COLORS <- c(
  "H3K27me3" = COLORS$h3k27me3[["H3K27me3 Gained"]],
  "H3K27ac"  = COLORS$h3k27ac[["H3K27ac Gained"]]
)

# Violin fills for the control and mutant mCH levels.
CONDITION_COLORS <- c(
  "Control" = COLORS$genotype[["ctrl"]],
  "Mutant"  = COLORS$genotype[["mut"]]
)

# Bins that switched subcompartment label between control and mutant.
CHANGED_LEVELS <- c("Label stable", "Label changed")
CHANGED_COLORS <- c(
  "Label stable"  = COLORS$significant[["Not Significant"]],
  "Label changed" = COLORS$significant[["Significant"]]
)

BIN_SIZE <- 100000L

# =============================================================================
# OPTIONS
# =============================================================================

parse_options <- function() {
  option_list <- list(
    make_option("--output-dir", type = "character", dest = "output_dir",
                default = OUTPUT_PATHS$chromatin,
                help = "Directory for figures and tables [default: %default]")
  )
  parse_args(OptionParser(option_list = option_list))
}

# =============================================================================
# FORMATTING HELPERS
# =============================================================================

#' Format a p-value for a plot annotation.
format_p <- function(p) {
  ifelse(p < 2.2e-16, "< 2.2e-16", sprintf("%.2e", p))
}

#' Format an integer with thousands separators and no padding.
format_n <- function(x) {
  format(x, big.mark = ",", trim = TRUE)
}

#' Write a data.frame as a TSV in the section output directory.
write_tsv_table <- function(df, out_dir, filename) {
  write_section_table(df, file.path(out_dir, filename))
}

# =============================================================================
# SUBCOMPARTMENT LOADING AND GENE ASSIGNMENT
# =============================================================================

#' Read the CALDER2 subcompartment label table and keep the labelled bins.
#'
#' The file carries the literal string "NA" in the label columns for bins
#' CALDER2 could not call. Those bins are dropped.
#'
#' @param path TSV with chr, bin_start, bin_end, ctrl_label, mut_label,
#'   continous_rank_ctrl, continous_rank_mut, label_changed.
#' @return data.frame of labelled bins with a logical label_changed column
load_subcompartment_bins <- function(path) {
  if (!file.exists(path)) {
    stop("CALDER2 subcompartment file not found: ", path)
  }

  bins <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                     quote = "")

  required <- c("chr", "bin_start", "bin_end", "ctrl_label", "mut_label",
                "continous_rank_ctrl", "continous_rank_mut", "label_changed")
  missing <- setdiff(required, colnames(bins))
  if (length(missing) > 0) {
    stop("Subcompartment file is missing columns: ",
         paste(missing, collapse = ", "), " in ", path)
  }

  n_all <- nrow(bins)
  bins <- bins[!is.na(bins$ctrl_label) & bins$ctrl_label != "NA", , drop = FALSE]

  unexpected <- setdiff(unique(bins$ctrl_label), SUBCOMPARTMENT_ORDER)
  if (length(unexpected) > 0) {
    stop("Unexpected ctrl_label values in ", path, ": ",
         paste(unexpected, collapse = ", "))
  }

  bins$label_changed <- as.logical(bins$label_changed)
  if (any(is.na(bins$label_changed))) {
    stop("label_changed has values that are neither TRUE nor FALSE in ", path)
  }

  bins$ctrl_label <- factor(bins$ctrl_label, levels = SUBCOMPARTMENT_ORDER)

  cat(sprintf("  %s of %s bins carry a control label\n",
              format_n(nrow(bins)), format_n(n_all)))
  for (lvl in SUBCOMPARTMENT_ORDER) {
    cat(sprintf("    %s: %s bins\n", lvl, format_n(sum(bins$ctrl_label == lvl))))
  }
  cat(sprintf("    label changed in mutant: %s bins, stable: %s bins\n",
              format_n(sum(bins$label_changed)),
              format_n(sum(!bins$label_changed))))
  bins
}

#' Build a GRanges from the labelled bins. bin_start is already 1-based.
bins_to_granges <- function(bins) {
  GRanges(
    seqnames = bins$chr,
    ranges = IRanges(start = bins$bin_start, end = bins$bin_end)
  )
}

#' Place each tested gene in the bin that contains its midpoint.
#'
#' mch_results uses BED-style 0-based starts, so the 1-based midpoint is
#' (start + 1 + end) %/% 2. Genes whose midpoint falls outside every labelled
#' bin (for example on chrX or chrY, which CALDER2 does not label) are dropped.
#'
#' @param mch_df mch_results
#' @param bins data.frame from load_subcompartment_bins()
#' @return data.frame of assigned genes with subcompartment columns
assign_genes_to_bins <- function(mch_df, bins) {
  midpoint <- (mch_df$start + 1L + mch_df$end) %/% 2L
  gene_mid_gr <- GRanges(
    seqnames = mch_df$chr,
    ranges = IRanges(start = midpoint, width = 1L)
  )

  bin_idx <- findOverlaps(gene_mid_gr, bins_to_granges(bins), select = "first")

  genes <- mch_df
  genes$gene_midpoint <- midpoint
  genes$bin_start  <- bins$bin_start[bin_idx]
  genes$bin_end    <- bins$bin_end[bin_idx]
  genes$subcompartment <- bins$ctrl_label[bin_idx]
  genes$mut_subcompartment <- bins$mut_label[bin_idx]
  genes$label_changed <- bins$label_changed[bin_idx]
  genes$rank_ctrl <- bins$continous_rank_ctrl[bin_idx]
  genes$rank_mut  <- bins$continous_rank_mut[bin_idx]

  n_total <- nrow(genes)
  genes <- genes[!is.na(bin_idx), , drop = FALSE]
  cat(sprintf("  %s of %s genes fall in a labelled %d kb bin (%.1f%%)\n",
              format_n(nrow(genes)), format_n(n_total), BIN_SIZE %/% 1000L,
              100 * nrow(genes) / n_total))

  n_dup <- sum(duplicated(genes$gene_name))
  cat(sprintf("  %s assigned gene bodies share a gene symbol with another row\n",
              format_n(n_dup)))

  genes$subcompartment <- factor(as.character(genes$subcompartment),
                                 levels = SUBCOMPARTMENT_ORDER)
  genes$changed_group <- factor(
    ifelse(genes$label_changed, "Label changed", "Label stable"),
    levels = CHANGED_LEVELS
  )
  genes
}

# =============================================================================
# COUNT SUMMARIES AND CATEGORICAL TESTS
# =============================================================================

#' Count genes, significant genes, and direction per subcompartment.
summarise_subcompartments <- function(genes) {
  genes %>%
    dplyr::group_by(subcompartment, .drop = FALSE) %>%
    dplyr::summarise(
      n_genes  = dplyr::n(),
      n_sig    = sum(mch_sig),
      n_hyper  = sum(mch_hyper),
      n_hypo   = sum(mch_hypo),
      median_mch_ctrl = median(mch_ctrl),
      median_mch_mut  = median(mch_mut),
      median_mch_diff = median(mch_diff),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      pct_sig   = 100 * n_sig / n_genes,
      pct_hyper_of_sig = 100 * n_hyper / (n_hyper + n_hypo),
      pct_hypo_of_sig  = 100 * n_hypo / (n_hyper + n_hypo)
    ) %>%
    as.data.frame()
}

#' Chi-squared test of subcompartment against mCH significance.
#'
#' @return list with the htest and a one-row summary data.frame
test_subcompartment_significance <- function(summary_df) {
  counts <- matrix(
    c(summary_df$n_sig, summary_df$n_genes - summary_df$n_sig),
    ncol = 2,
    dimnames = list(as.character(summary_df$subcompartment),
                    c("Significant", "NotSignificant"))
  )
  ct <- chisq.test(counts)

  cat(sprintf("  Chi-squared (subcompartment x mCH significance): X2 = %.1f, df = %d, p = %s\n",
              unname(ct$statistic), unname(ct$parameter), format_p(ct$p.value)))

  list(
    test = ct,
    table = data.frame(
      test = "subcompartment x mch_sig",
      chi_squared = unname(ct$statistic),
      df = unname(ct$parameter),
      p_value = ct$p.value,
      stringsAsFactors = FALSE
    ),
    counts = counts
  )
}

#' Chi-squared test of subcompartment against direction, significant genes only.
test_subcompartment_direction <- function(summary_df) {
  counts <- matrix(
    c(summary_df$n_hyper, summary_df$n_hypo),
    ncol = 2,
    dimnames = list(as.character(summary_df$subcompartment),
                    c("Hypermethylated", "Hypomethylated"))
  )
  ct <- chisq.test(counts)

  cat(sprintf("  Chi-squared (subcompartment x direction): X2 = %.1f, df = %d, p = %s\n",
              unname(ct$statistic), unname(ct$parameter), format_p(ct$p.value)))

  list(
    test = ct,
    table = data.frame(
      test = "subcompartment x direction (significant genes)",
      chi_squared = unname(ct$statistic),
      df = unname(ct$parameter),
      p_value = ct$p.value,
      stringsAsFactors = FALSE
    )
  )
}

#' Run the gene-level 2x2 Fisher tests of this section through the registry.
#'
#' Every test compares a subcompartment or histone membership against an mCH
#' outcome, so section 40_04 can validate each odds ratio by permutation.
register_subcompartment_fisher_tests <- function(genes, out_dir) {
  genes$in_a_compartment <- genes$subcompartment %in% c("A.1", "A.2")
  genes$in_b2            <- genes$subcompartment == "B.2"

  register_fisher_test(
    section = SECTION_ID, test_id = "a_compartment_sig",
    description = paste("Genes in A subcompartments (A.1 or A.2) versus B",
                        "subcompartments, against significant mCH change."),
    gene_df = genes, row_var = "in_a_compartment", col_var = "mch_sig",
    output_dir = out_dir)

  register_fisher_test(
    section = SECTION_ID, test_id = "b2_sig",
    description = paste("Genes in the B.2 constitutive heterochromatin",
                        "subcompartment versus all other subcompartments,",
                        "against significant mCH change."),
    gene_df = genes, row_var = "in_b2", col_var = "mch_sig",
    output_dir = out_dir)

  register_fisher_test(
    section = SECTION_ID, test_id = "label_changed_sig",
    description = paste("Genes in bins whose subcompartment label changed in",
                        "the mutant versus stable bins, against significant",
                        "mCH change."),
    gene_df = genes, row_var = "label_changed", col_var = "mch_sig",
    output_dir = out_dir)

  register_fisher_test(
    section = SECTION_ID, test_id = "k27me3_sig",
    description = paste("Genes whose body overlaps an H3K27me3 consensus peak,",
                        "against significant mCH change."),
    gene_df = genes, row_var = "h3k27me3_overlap", col_var = "mch_sig",
    output_dir = out_dir)

  register_fisher_test(
    section = SECTION_ID, test_id = "k27ac_sig",
    description = paste("Genes whose body overlaps an H3K27ac consensus peak,",
                        "against significant mCH change."),
    gene_df = genes, row_var = "h3k27ac_overlap", col_var = "mch_sig",
    output_dir = out_dir)

  sig_genes <- genes[genes$mch_sig, , drop = FALSE]
  register_fisher_test(
    section = SECTION_ID, test_id = "a_compartment_hyper",
    description = paste("Among genes with significant mCH change, whether A",
                        "subcompartment membership associates with gaining",
                        "rather than losing mCH."),
    gene_df = sig_genes, row_var = "in_a_compartment", col_var = "mch_hyper",
    output_dir = out_dir)

  invisible(NULL)
}

# =============================================================================
# DISTRIBUTION TESTS
# =============================================================================

#' Paired Wilcoxon of control against mutant mCH within each subcompartment.
paired_wilcoxon_by_subcompartment <- function(genes) {
  rows <- lapply(SUBCOMPARTMENT_ORDER, function(lvl) {
    sub <- genes[genes$subcompartment == lvl, , drop = FALSE]
    wt <- wilcox.test(sub$mch_ctrl, sub$mch_mut, paired = TRUE)
    data.frame(
      subcompartment = lvl,
      n_genes = nrow(sub),
      median_mch_ctrl = median(sub$mch_ctrl),
      median_mch_mut = median(sub$mch_mut),
      median_mch_diff = median(sub$mch_diff),
      statistic = unname(wt$statistic),
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)

  for (i in seq_len(nrow(out))) {
    cat(sprintf("    %s: n = %s, ctrl median = %.4f, mut median = %.4f, paired p = %s\n",
                out$subcompartment[i], format_n(out$n_genes[i]),
                out$median_mch_ctrl[i], out$median_mch_mut[i],
                format_p(out$p_value[i])))
  }
  out
}

#' Kruskal-Wallis across groups plus BH-adjusted pairwise Wilcoxon tests.
#'
#' @param values Numeric vector.
#' @param groups Factor of the same length.
#' @param measure Name of the measured quantity, copied into the output.
#' @return list with an omnibus one-row data.frame and a pairwise data.frame
across_group_wilcoxon <- function(values, groups, measure) {
  keep <- !is.na(values) & !is.na(groups)
  values <- values[keep]
  groups <- droplevels(factor(groups[keep]))

  kw <- kruskal.test(values, groups)
  omnibus <- data.frame(
    measure = measure,
    test = "Kruskal-Wallis",
    n = length(values),
    n_groups = nlevels(groups),
    statistic = unname(kw$statistic),
    df = unname(kw$parameter),
    p_value = kw$p.value,
    stringsAsFactors = FALSE
  )

  pw <- pairwise.wilcox.test(values, groups, p.adjust.method = "BH")
  pw_long <- as.data.frame(as.table(pw$p.value), stringsAsFactors = FALSE)
  colnames(pw_long) <- c("group_1", "group_2", "p_adjusted")
  pw_long <- pw_long[!is.na(pw_long$p_adjusted), , drop = FALSE]

  medians <- tapply(values, groups, median)
  counts  <- tapply(values, groups, length)
  pw_long$measure     <- measure
  pw_long$n_group_1   <- as.integer(counts[pw_long$group_1])
  pw_long$n_group_2   <- as.integer(counts[pw_long$group_2])
  pw_long$median_group_1 <- as.numeric(medians[pw_long$group_1])
  pw_long$median_group_2 <- as.numeric(medians[pw_long$group_2])
  pw_long$median_difference <- pw_long$median_group_1 - pw_long$median_group_2
  pw_long <- pw_long[, c("measure", "group_1", "group_2", "n_group_1",
                         "n_group_2", "median_group_1", "median_group_2",
                         "median_difference", "p_adjusted")]

  cat(sprintf("  %s: Kruskal-Wallis chi2 = %.1f, df = %d, p = %s\n",
              measure, unname(kw$statistic), unname(kw$parameter),
              format_p(kw$p.value)))
  for (i in seq_len(nrow(pw_long))) {
    cat(sprintf("    %s vs %s: median difference = %.5f, BH p = %s\n",
                pw_long$group_1[i], pw_long$group_2[i],
                pw_long$median_difference[i], format_p(pw_long$p_adjusted[i])))
  }

  list(omnibus = omnibus, pairwise = pw_long)
}

#' Wilcoxon rank-sum test of mCH difference, changed bins against stable bins.
test_label_changed_effect <- function(genes) {
  changed <- genes$mch_diff[genes$label_changed]
  stable  <- genes$mch_diff[!genes$label_changed]
  wt <- wilcox.test(changed, stable)

  out <- data.frame(
    measure = "mch_diff",
    test = "Wilcoxon rank-sum, changed vs stable subcompartment label",
    n_changed = length(changed),
    n_stable = length(stable),
    median_changed = median(changed),
    median_stable = median(stable),
    median_difference = median(changed) - median(stable),
    statistic = unname(wt$statistic),
    p_value = wt$p.value,
    stringsAsFactors = FALSE
  )

  cat(sprintf("  Label changed (n = %s, median = %.5f) vs stable (n = %s, median = %.5f): p = %s\n",
              format_n(out$n_changed), out$median_changed,
              format_n(out$n_stable), out$median_stable,
              format_p(out$p_value)))
  out
}

# =============================================================================
# HISTONE MARK OVERLAP
# =============================================================================

#' Add H3K27me3 and H3K27ac gene-body overlap flags and a mark category.
add_histone_overlaps <- function(genes) {
  if (!file.exists(CHIP_PATHS$h3k27me3)) {
    stop("H3K27me3 consensus peaks not found: ", CHIP_PATHS$h3k27me3)
  }
  if (!file.exists(CHIP_PATHS$h3k27ac)) {
    stop("H3K27ac consensus peaks not found: ", CHIP_PATHS$h3k27ac)
  }

  k27me3 <- load_chip_peaks(CHIP_PATHS$h3k27me3, "H3K27me3")
  k27ac  <- load_chip_peaks(CHIP_PATHS$h3k27ac,  "H3K27ac")

  gene_gr <- GRanges(
    seqnames = genes$chr,
    ranges = IRanges(start = genes$start + 1L, end = genes$end)
  )

  overlaps <- compute_chip_overlaps(gene_gr,
                                    list(h3k27me3 = k27me3, h3k27ac = k27ac))
  genes$h3k27me3_overlap <- overlaps$h3k27me3_overlap
  genes$h3k27ac_overlap  <- overlaps$h3k27ac_overlap

  genes$mark_category <- dplyr::case_when(
    genes$h3k27me3_overlap & genes$h3k27ac_overlap  ~ "H3K27me3 + H3K27ac",
    genes$h3k27me3_overlap                          ~ "H3K27me3 only",
    genes$h3k27ac_overlap                           ~ "H3K27ac only",
    TRUE                                            ~ "Neither"
  )
  genes$mark_category <- factor(
    genes$mark_category,
    levels = c("H3K27ac only", "H3K27me3 + H3K27ac", "H3K27me3 only", "Neither")
  )

  cat(sprintf("  Gene bodies overlapping H3K27me3: %s, H3K27ac: %s\n",
              format_n(sum(genes$h3k27me3_overlap)),
              format_n(sum(genes$h3k27ac_overlap))))
  genes
}

#' Overlap rate for each mark within each subcompartment, in long form.
summarise_histone_by_subcompartment <- function(genes) {
  long <- rbind(
    data.frame(subcompartment = genes$subcompartment, mark = "H3K27me3",
               overlap = genes$h3k27me3_overlap, stringsAsFactors = FALSE),
    data.frame(subcompartment = genes$subcompartment, mark = "H3K27ac",
               overlap = genes$h3k27ac_overlap, stringsAsFactors = FALSE)
  )

  long %>%
    dplyr::group_by(subcompartment, mark) %>%
    dplyr::summarise(
      n_genes = dplyr::n(),
      n_overlap = sum(overlap),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      pct_overlap = 100 * n_overlap / n_genes,
      subcompartment = factor(as.character(subcompartment),
                              levels = SUBCOMPARTMENT_ORDER),
      mark = factor(mark, levels = c("H3K27me3", "H3K27ac"))
    ) %>%
    as.data.frame()
}

#' Gene counts and mCH outcome for each H3K27me3 / H3K27ac category.
summarise_mark_category <- function(genes) {
  genes %>%
    dplyr::group_by(mark_category, .drop = FALSE) %>%
    dplyr::summarise(
      n_genes = dplyr::n(),
      n_sig   = sum(mch_sig),
      n_hyper = sum(mch_hyper),
      n_hypo  = sum(mch_hypo),
      median_mch_ctrl = median(mch_ctrl),
      median_mch_diff = median(mch_diff),
      .groups = "drop"
    ) %>%
    dplyr::mutate(pct_sig = 100 * n_sig / n_genes) %>%
    as.data.frame()
}

# =============================================================================
# FIGURES
# =============================================================================

#' Bar of the significant-gene percentage in each subcompartment.
plot_significance_rate <- function(summary_df, chisq_p) {
  ggplot(summary_df, aes(x = subcompartment, y = pct_sig,
                         fill = subcompartment)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%s / %s\n(%.1f%%)",
                                  format_n(n_sig), format_n(n_genes), pct_sig)),
              vjust = -0.25, size = 3.1, lineheight = 0.9) +
    scale_fill_manual(values = SUBCOMPARTMENT_COLORS, guide = "none") +
    scale_x_discrete(labels = SUBCMPT_LABELS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
    labs(
      title = "Significant mCH genes by CALDER2 subcompartment",
      subtitle = sprintf("Chi-squared subcompartment x significance: p = %s",
                         format_p(chisq_p)),
      x = "CALDER2 subcompartment (control)",
      y = sprintf("%% of genes with significant mCH change (FDR < %.2f)",
                  Q_THRESHOLD)
    ) +
    theme_emseq()
}

#' Stacked bar of the hyper/hypo split among significant genes.
plot_direction_split <- function(summary_df, chisq_p) {
  dir_df <- summary_df %>%
    dplyr::select(subcompartment, n_hyper, n_hypo) %>%
    tidyr::pivot_longer(cols = c(n_hyper, n_hypo),
                        names_to = "direction", values_to = "n") %>%
    dplyr::mutate(
      direction = factor(ifelse(direction == "n_hyper",
                                "Hypermethylated", "Hypomethylated"),
                         levels = c("Hypermethylated", "Hypomethylated"))
    ) %>%
    dplyr::group_by(subcompartment) %>%
    dplyr::mutate(pct = 100 * n / sum(n)) %>%
    dplyr::ungroup() %>%
    as.data.frame()

  ggplot(dir_df, aes(x = subcompartment, y = pct, fill = direction)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%\n(%s)", pct, format_n(n))),
              position = position_stack(vjust = 0.5), size = 3.1,
              colour = "white", fontface = "bold", lineheight = 0.9) +
    scale_fill_manual(values = COLORS$direction, name = "Direction") +
    scale_x_discrete(labels = SUBCMPT_LABELS) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = "mCH change direction by subcompartment",
      subtitle = sprintf("Significant genes only | chi-squared subcompartment x direction: p = %s",
                         format_p(chisq_p)),
      x = "CALDER2 subcompartment (control)",
      y = "% of significant genes"
    ) +
    theme_emseq() +
    theme(legend.position = "top")
}

#' Violin of control and mutant mCH level, one facet per subcompartment.
#'
#' Facet strips carry the paired Wilcoxon p-value. Each violin is annotated
#' with its n and median.
plot_mch_level <- function(level_df, level_summary, paired_tests, kruskal_p) {
  level_summary$label <- group_label(level_summary, digits = 3)

  facet_labels <- setNames(
    sprintf("%s\npaired p = %s",
            SUBCMPT_LABELS[paired_tests$subcompartment],
            format_p(paired_tests$p_value)),
    paired_tests$subcompartment
  )

  ggplot(level_df, aes(x = condition, y = mch_pct, fill = condition)) +
    geom_violin(alpha = 0.75, scale = "width", trim = TRUE,
                linewidth = 0.3) +
    geom_boxplot(width = 0.16, outlier.size = 0.2, alpha = 0.85,
                 linewidth = 0.3) +
    geom_text(data = level_summary,
              aes(x = condition, y = Inf, label = label),
              inherit.aes = FALSE, vjust = 1.05, size = 2.8,
              lineheight = 0.9) +
    facet_wrap(~ subcompartment, nrow = 1,
               labeller = labeller(subcompartment = facet_labels)) +
    scale_fill_manual(values = CONDITION_COLORS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    labs(
      title = "Gene-body mCH level by subcompartment",
      subtitle = sprintf("Kruskal-Wallis of control mCH across subcompartments: p = %s",
                         format_p(kruskal_p)),
      x = NULL, y = "mCH level (%)"
    ) +
    theme_emseq()
}

#' Violin of the mCH difference in each subcompartment.
plot_mch_diff <- function(genes, diff_summary, kruskal_p) {
  diff_summary$label <- group_label(diff_summary, digits = 3)

  ggplot(genes, aes(x = subcompartment, y = mch_diff_pct,
                    fill = subcompartment)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_violin(alpha = 0.75, scale = "width", trim = TRUE, linewidth = 0.3) +
    geom_boxplot(width = 0.14, outlier.size = 0.2, alpha = 0.85,
                 linewidth = 0.3) +
    geom_text(data = diff_summary,
              aes(x = subcompartment, y = Inf, label = label),
              inherit.aes = FALSE, vjust = 1.05, size = 2.8,
              lineheight = 0.9) +
    scale_fill_manual(values = SUBCOMPARTMENT_COLORS, guide = "none") +
    scale_x_discrete(labels = SUBCMPT_LABELS) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    labs(
      title = "mCH difference by subcompartment",
      subtitle = sprintf("Mutant minus control | Kruskal-Wallis across subcompartments: p = %s",
                         format_p(kruskal_p)),
      x = "CALDER2 subcompartment (control)",
      y = "mCH difference (percentage points)"
    ) +
    theme_emseq()
}

#' Violin of the mCH difference in bins that changed label against stable bins.
plot_label_changed <- function(genes, changed_summary, wilcox_p) {
  changed_summary$label <- group_label(changed_summary, digits = 3)

  ggplot(genes, aes(x = changed_group, y = mch_diff_pct, fill = changed_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_violin(alpha = 0.75, scale = "width", trim = TRUE, linewidth = 0.3) +
    geom_boxplot(width = 0.14, outlier.size = 0.2, alpha = 0.85,
                 linewidth = 0.3) +
    geom_text(data = changed_summary,
              aes(x = changed_group, y = Inf, label = label),
              inherit.aes = FALSE, vjust = 1.05, size = 2.8,
              lineheight = 0.9) +
    scale_fill_manual(values = CHANGED_COLORS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    labs(
      title = "mCH difference by subcompartment label change",
      subtitle = sprintf("Wilcoxon rank-sum: p = %s", format_p(wilcox_p)),
      x = "Bin label between control and mutant",
      y = "mCH difference (percentage points)"
    ) +
    theme_emseq()
}

#' Grouped bar of the H3K27me3 and H3K27ac overlap rate per subcompartment.
plot_histone_overlay <- function(histone_df) {
  ggplot(histone_df, aes(x = subcompartment, y = pct_overlap, fill = mark)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65,
             colour = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%s\n(%.1f%%)",
                                  format_n(n_overlap), pct_overlap)),
              position = position_dodge(width = 0.75), vjust = -0.25,
              size = 2.9, lineheight = 0.9) +
    scale_fill_manual(values = HISTONE_COLORS, name = "Consensus peak") +
    scale_x_discrete(labels = SUBCMPT_LABELS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
    labs(
      title = "Histone mark overlap by subcompartment",
      subtitle = "Fraction of gene bodies overlapping a consensus peak",
      x = "CALDER2 subcompartment (control)",
      y = "% of gene bodies overlapping the mark"
    ) +
    theme_emseq() +
    theme(legend.position = "top")
}

# =============================================================================
# LONG-FORM DATA FOR THE VIOLINS
# =============================================================================

#' Stack control and mutant mCH levels into one long data.frame.
build_level_long <- function(genes) {
  long <- rbind(
    data.frame(gene_name = genes$gene_name,
               subcompartment = as.character(genes$subcompartment),
               condition = "Control", mch_pct = genes$mch_ctrl * 100,
               stringsAsFactors = FALSE),
    data.frame(gene_name = genes$gene_name,
               subcompartment = as.character(genes$subcompartment),
               condition = "Mutant", mch_pct = genes$mch_mut * 100,
               stringsAsFactors = FALSE)
  )
  long$subcompartment <- factor(long$subcompartment,
                                levels = SUBCOMPARTMENT_ORDER)
  long$condition <- factor(long$condition, levels = c("Control", "Mutant"))
  long$group_key <- paste(as.character(long$subcompartment),
                          as.character(long$condition), sep = "|")
  long
}

#' Split the combined "subcompartment|condition" key of a summary table.
split_group_key <- function(summary_df) {
  parts <- strsplit(summary_df$group_key, "|", fixed = TRUE)
  summary_df$subcompartment <- factor(
    vapply(parts, `[`, character(1), 1), levels = SUBCOMPARTMENT_ORDER)
  summary_df$condition <- factor(
    vapply(parts, `[`, character(1), 2), levels = c("Control", "Mutant"))
  summary_df$group_key <- NULL
  summary_df[, c("subcompartment", "condition", "n", "median", "mean",
                 "q25", "q75")]
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opt <- parse_options()
  out_dir <- opt$output_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("================================================================================\n")
  cat("SECTION 10_04: mCH CHANGE BY CALDER2 SUBCOMPARTMENT\n")
  cat("================================================================================\n")
  cat("Output directory: ", out_dir, "\n\n", sep = "")

  # ---- Load and assign ------------------------------------------------------

  cat("Loading CALDER2 subcompartment labels...\n")
  bins <- load_subcompartment_bins(HIC_PATHS$subcompartments)

  cat("\nAssigning genes to subcompartments by gene midpoint...\n")
  genes <- assign_genes_to_bins(mch_results, bins)

  cat("\nComputing histone consensus peak overlap...\n")
  genes <- add_histone_overlaps(genes)

  genes$mch_diff_pct <- genes$mch_diff * 100

  # ---- Counts and categorical tests -----------------------------------------

  cat("\nSummarising genes per subcompartment...\n")
  summary_df <- summarise_subcompartments(genes)
  for (i in seq_len(nrow(summary_df))) {
    cat(sprintf("    %s: %s genes, %s significant (%.1f%%), %s hyper, %s hypo\n",
                as.character(summary_df$subcompartment[i]),
                format_n(summary_df$n_genes[i]),
                format_n(summary_df$n_sig[i]), summary_df$pct_sig[i],
                format_n(summary_df$n_hyper[i]), format_n(summary_df$n_hypo[i])))
  }

  cat("\nCategorical tests...\n")
  sig_chisq <- test_subcompartment_significance(summary_df)
  dir_chisq <- test_subcompartment_direction(summary_df)

  cat("\nGene-level Fisher tests...\n")
  register_subcompartment_fisher_tests(genes, out_dir)

  # ---- Distribution tests ---------------------------------------------------

  cat("\nPaired Wilcoxon of control vs mutant mCH per subcompartment...\n")
  paired_tests <- paired_wilcoxon_by_subcompartment(genes)

  cat("\nControl mCH level across subcompartments...\n")
  ctrl_tests <- across_group_wilcoxon(genes$mch_ctrl, genes$subcompartment,
                                      "mch_ctrl")

  cat("\nmCH difference across subcompartments...\n")
  diff_tests <- across_group_wilcoxon(genes$mch_diff, genes$subcompartment,
                                      "mch_diff")

  cat("\nmCH difference in changed versus stable bins...\n")
  changed_test <- test_label_changed_effect(genes)

  # ---- Group summaries backing the violins ----------------------------------

  level_df <- build_level_long(genes)
  level_summary <- split_group_key(
    summarise_groups(level_df, "group_key", "mch_pct"))
  diff_summary <- summarise_groups(genes, "subcompartment", "mch_diff_pct")
  changed_summary <- summarise_groups(genes, "changed_group", "mch_diff_pct")

  histone_df <- summarise_histone_by_subcompartment(genes)
  mark_summary <- summarise_mark_category(genes)

  # ---- Figures --------------------------------------------------------------

  cat("\nBuilding figures...\n")

  p_sig <- plot_significance_rate(summary_df, sig_chisq$table$p_value)
  save_multiformat_ggplot(
    p_sig, file.path(out_dir, "10_04a_significance_rate_by_subcompartment"),
    width = 9, height = 7)

  p_dir <- plot_direction_split(summary_df, dir_chisq$table$p_value)
  save_multiformat_ggplot(
    p_dir, file.path(out_dir, "10_04b_direction_split_by_subcompartment"),
    width = 9, height = 7)

  p_level <- plot_mch_level(level_df, level_summary, paired_tests,
                            ctrl_tests$omnibus$p_value)
  save_multiformat_ggplot(
    p_level, file.path(out_dir, "10_04c_mch_level_by_subcompartment"),
    width = 14, height = 7)

  p_diff <- plot_mch_diff(genes, diff_summary, diff_tests$omnibus$p_value)
  p_changed <- plot_label_changed(genes, changed_summary,
                                  changed_test$p_value)
  save_multiformat_ggplot(
    p_diff + p_changed + patchwork::plot_layout(widths = c(2, 1)),
    file.path(out_dir, "10_04d_mch_diff_by_subcompartment"),
    width = 15, height = 7)

  p_histone <- plot_histone_overlay(histone_df)
  save_multiformat_ggplot(
    p_histone, file.path(out_dir, "10_04e_histone_overlap_by_subcompartment"),
    width = 10, height = 7)

  composite <- ((p_sig | p_dir) / p_level / (p_diff | p_changed | p_histone)) +
    patchwork::plot_annotation(
      title = "mCH change stratified by CALDER2 subcompartment",
      subtitle = sprintf("%s genes assigned to %s labelled 100 kb bins",
                         format_n(nrow(genes)), format_n(nrow(bins))),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 11)
      )
    )
  save_multiformat_ggplot(
    composite, file.path(out_dir, "10_04_subcompartment_composite"),
    width = 22, height = 22)

  # ---- Tables ---------------------------------------------------------------

  cat("\nWriting tables...\n")

  gene_table <- genes[, c("gene_name", "gene_id", "chr", "start", "end",
                          "gene_length", "gene_midpoint", "bin_start", "bin_end",
                          "subcompartment", "mut_subcompartment",
                          "label_changed", "rank_ctrl", "rank_mut",
                          "h3k27me3_overlap", "h3k27ac_overlap", "mark_category",
                          "mch_ctrl", "mch_mut", "mch_diff", "edger_logFC",
                          "edger_fdr", "mch_sig", "mch_direction",
                          "mch_hyper", "mch_hypo")]
  write_tsv_table(gene_table, out_dir, "10_04_gene_subcompartment_assignment.tsv")

  write_tsv_table(summary_df, out_dir, "10_04_subcompartment_summary.tsv")
  write_tsv_table(rbind(sig_chisq$table, dir_chisq$table), out_dir,
                  "10_04_subcompartment_chisq.tsv")
  write_tsv_table(paired_tests, out_dir,
                  "10_04_mch_level_paired_wilcoxon.tsv")
  write_tsv_table(rbind(ctrl_tests$omnibus, diff_tests$omnibus), out_dir,
                  "10_04_across_subcompartment_kruskal.tsv")
  write_tsv_table(rbind(ctrl_tests$pairwise, diff_tests$pairwise), out_dir,
                  "10_04_across_subcompartment_pairwise_wilcoxon.tsv")
  write_tsv_table(changed_test, out_dir,
                  "10_04_label_changed_wilcoxon.tsv")
  write_tsv_table(histone_df, out_dir,
                  "10_04_histone_overlap_by_subcompartment.tsv")
  write_tsv_table(mark_summary, out_dir,
                  "10_04_histone_mark_category_summary.tsv")

  level_summary$measure <- "mch_level_pct"
  diff_summary$measure <- "mch_diff_pct"
  changed_summary$measure <- "mch_diff_pct"
  write_tsv_table(level_summary, out_dir,
                  "10_04_violin_group_summary_mch_level.tsv")
  write_tsv_table(diff_summary, out_dir,
                  "10_04_violin_group_summary_mch_diff.tsv")
  write_tsv_table(changed_summary, out_dir,
                  "10_04_violin_group_summary_label_changed.tsv")

  cat("\n================================================================================\n")
  cat("SECTION 10_04 COMPLETE\n")
  cat("================================================================================\n")
}

main()
