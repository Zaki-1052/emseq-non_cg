# scripts/sections/70_neuronal/70_02_chromatin_remodeling.R
#
# Section 70_02: chromatin remodeling at neuronal genes after BAP1 loss.
#
# Neuronal genes carry a constitutively high H2AK119ub level (section 70_01).
# This section asks whether they also remodel their chromatin more strongly
# than other genes when BAP1 is lost. The model predicts four directions in the
# mutant: ATAC down, H3K27ac down, H3K27me3 up, H2AK119ub up. Every test here
# compares the neuronal gene class against every other gene of the same
# universe, on the DiffBind log2 fold change of each mark.
#
# Two properties travel with neuronal identity and could produce the same
# result on their own, so both are held fixed:
#   K119ub signal deciles   separate neuronal identity from the K119ub level
#   gene-length quintiles   separate neuronal identity from gene length
# A linear model per mark holds both at once and reports the adjusted neuronal
# coefficient beside the unadjusted one.
#
# No methylation data enters this section. The question is about chromatin, so
# the methylation columns of the input table are dropped on load.
#
# Reads:
#   HANDOFF_PATHS$neuronal_gene_set      neuronal gene set, written by 70_01
#   HANDOFF_PATHS$gene_level_all_marks   gene-level mark folds, written by 20_02
#   DIFFBIND_PATHS$k119ub_gene_signal    per-gene K119ub level in ctrl and mut
#
# Writes into --output-dir (default results/sections/70_neuronal/):
#   70_02a  mark fold distributions, neuronal against other
#   70_02b  effect sizes across marks, unadjusted and adjusted
#   70_02c  the comparison inside K119ub signal deciles
#   70_02d  the comparison inside gene-length quintiles
#   70_02e  count of marks moving in the predicted direction per gene
#   70_02f  K119ub level against mark fold change, by gene class
#   70_02_universe_per_gene.tsv   per-gene table behind every figure
#   fisher_tables/                gene tables behind the registered Fisher tests
#
# Adapted from Biomodal section 73 (neuronal chromatin remodeling).

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)

# =============================================================================
# SECTION CONSTANTS
# =============================================================================

SECTION_ID <- "70_02"

MARK_ORDER <- c("atac", "k27ac", "k27me3", "k119ub")

# display    name used in figures and tables
# predicted  direction of the log2FC that the model predicts in the mutant
MARK_META <- list(
  atac   = list(display = "ATAC-seq",  predicted = "down"),
  k27ac  = list(display = "H3K27ac",   predicted = "down"),
  k27me3 = list(display = "H3K27me3",  predicted = "up"),
  k119ub = list(display = "H2AK119ub", predicted = "up")
)

GROUP_ORDER <- c("Neuronal", "Other")

GROUP_COLORS <- c(
  "Neuronal" = COLORS$k119ub[["K119ub Gained"]],
  "Other"    = "grey70"
)

# Point colour for an effect whose sign agrees with the predicted direction.
PREDICTION_COLORS <- c(
  "Consistent" = unname(COLORS$quadrant[["Q1"]]),
  "Opposite"   = unname(COLORS$quadrant[["Q3"]])
)

N_DECILES <- 10
N_LENGTH_QUINTILES <- 5

# Column contracts of the three input tables.
MARK_TABLE_COLUMNS <- c("gene_name", "gene_id", "chr", "start", "end",
                        "gene_length",
                        paste0(MARK_ORDER, "_fold"),
                        paste0(MARK_ORDER, "_fdr"),
                        paste0(MARK_ORDER, "_n_peaks"))

NEURONAL_SET_COLUMNS <- c("gene", "source", "is_derived", "is_external",
                          "k119ub_ctrl_signal", "k119ub_log2fc",
                          "k119ub_decile")

SIGNAL_COLUMNS <- c("symbol", "chr", "width", "gb_ctrl_signal", "gb_mut_signal",
                    "gb_log2fc", "gb_signal_class")

# Genes of this K119ub signal class carry a usable level in both conditions.
QUANTIFIABLE_CLASS <- "quantifiable"

# The handoff and the signal table must report the same control K119ub level
# for a shared gene. write.table rounds, so the check allows this much drift.
SIGNAL_TOLERANCE <- 1e-6

# =============================================================================
# COMMAND LINE
# =============================================================================

parse_section_args <- function() {
  option_list <- list(
    make_option("--output-dir", dest = "output_dir", type = "character",
                default = OUTPUT_PATHS$neuronal,
                help = "Directory for figures and tables [default: %default]"),
    make_option("--fdr-threshold", dest = "fdr_threshold", type = "double",
                default = Q_THRESHOLD,
                help = "FDR cutoff for peak significance [default: %default]"),
    make_option("--min-group-size", dest = "min_group_size", type = "integer",
                default = 10L,
                help = paste("Smallest group allowed inside a decile or",
                             "quintile cell [default: %default]"))
  )
  opt <- parse_args(OptionParser(option_list = option_list))

  if (opt$fdr_threshold <= 0 || opt$fdr_threshold >= 1) {
    stop("--fdr-threshold must be between 0 and 1, got ", opt$fdr_threshold)
  }
  if (opt$min_group_size < 3) {
    stop("--min-group-size must be at least 3, got ", opt$min_group_size)
  }
  opt
}

# =============================================================================
# SMALL OUTPUT HELPERS
# =============================================================================


format_p <- function(p) {
  ifelse(is.na(p), "p = NA",
    ifelse(p < 2.2e-16, "p < 2.2e-16", sprintf("p = %.2e", p)))
}

mark_display <- function(mark) MARK_META[[mark]]$display

mark_displays <- function() {
  vapply(MARK_META[MARK_ORDER], `[[`, character(1), "display")
}

#' Sign of the log2 fold change that the model predicts for a mark.
predicted_sign <- function(mark) {
  if (MARK_META[[mark]]$predicted == "down") return(-1)
  1
}

#' Summarise n, median, mean and quartiles for every combination of the key
#' columns. This generalises summarise_groups() from the shared config to more
#' than one grouping column. It returns data only. The on-plot text comes from
#' group_label() at the plot site.
summarise_cells <- function(df, keys, value_col) {
  df %>%
    dplyr::filter(!is.na(.data[[value_col]])) %>%
    dplyr::group_by(!!!rlang::syms(keys)) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median = median(.data[[value_col]]),
      mean = mean(.data[[value_col]]),
      q25 = unname(quantile(.data[[value_col]], 0.25)),
      q75 = unname(quantile(.data[[value_col]], 0.75)),
      .groups = "drop"
    ) %>%
    as.data.frame()
}

#' Cut a numeric vector into equal-count bins, stopping when the quantile
#' breaks are not strictly increasing.
cut_into_bins <- function(values, n_bins, label) {
  breaks <- quantile(values, probs = seq(0, 1, length.out = n_bins + 1))
  if (any(duplicated(breaks))) {
    stop("Cannot split ", label, " into ", n_bins,
         " equal-count bins: the quantile breaks repeat a value (",
         paste(sprintf("%.6g", breaks), collapse = ", "), ").")
  }
  bins <- as.integer(cut(values, breaks = breaks, include.lowest = TRUE,
                         labels = FALSE))
  if (anyNA(bins)) {
    stop(sum(is.na(bins)), " values fell outside the ", label, " bin breaks.")
  }
  bins
}

# =============================================================================
# WILCOXON COMPARISON OF THE TWO GENE CLASSES
# =============================================================================

#' Compare one numeric vector between neuronal and other genes.
#'
#' Reports the median of each class, the shift between them, the rank-biserial
#' correlation as the effect size, and the Wilcoxon p value. With with_ci the
#' Hodges-Lehmann location shift and its 95% confidence interval come from the
#' same test.
#'
#' @param values Numeric vector, one entry per gene. NA entries are dropped.
#' @param is_neuronal Logical vector of the same length.
#' @param with_ci Ask wilcox.test for the location shift and its interval.
#' @return one-row data.frame
compare_gene_classes <- function(values, is_neuronal, with_ci = FALSE) {
  neuronal <- values[is_neuronal & !is.na(values)]
  other <- values[!is_neuronal & !is.na(values)]
  if (length(neuronal) == 0 || length(other) == 0) {
    stop("compare_gene_classes(): one class holds no value (neuronal ",
         length(neuronal), ", other ", length(other), ").")
  }

  test <- wilcox.test(neuronal, other, exact = FALSE, conf.int = with_ci)
  w <- unname(test$statistic)
  rank_biserial <- 2 * w / (length(neuronal) * length(other)) - 1

  data.frame(
    n_neuronal = length(neuronal),
    n_other = length(other),
    median_neuronal = median(neuronal),
    median_other = median(other),
    mean_neuronal = mean(neuronal),
    mean_other = mean(other),
    delta_median = median(neuronal) - median(other),
    rank_biserial = rank_biserial,
    location_shift = if (with_ci) unname(test$estimate) else NA_real_,
    shift_ci_lower = if (with_ci) test$conf.int[1] else NA_real_,
    shift_ci_upper = if (with_ci) test$conf.int[2] else NA_real_,
    wilcoxon_w = w,
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

#' Stop when a stratum cell holds too few genes of either class to compare.
require_group_sizes <- function(n_neuronal, n_other, min_group_size, context) {
  if (n_neuronal >= min_group_size && n_other >= min_group_size) return(invisible(NULL))
  stop(context, " holds ", n_neuronal, " neuronal and ", n_other,
       " other genes with a measured fold change. The comparison needs at ",
       "least ", min_group_size, " of each. Lower --min-group-size or use ",
       "fewer strata.")
}

#' Add the predicted-direction columns that every effect table carries.
annotate_prediction <- function(stats_df, mark) {
  stats_df$mark <- mark
  stats_df$mark_display <- mark_display(mark)
  stats_df$predicted_direction <- MARK_META[[mark]]$predicted
  stats_df$observed_direction <- ifelse(stats_df$delta_median > 0, "up", "down")
  stats_df$matches_prediction <- sign(stats_df$delta_median) == predicted_sign(mark)
  stats_df
}

# =============================================================================
# STEP 1: INPUT TABLES
# =============================================================================

load_mark_table <- function() {
  path <- HANDOFF_PATHS$gene_level_all_marks
  if (!file.exists(path)) {
    stop("Gene-level all-marks table not found: ", path,
         "\nRun section 20_02 (20_02_multi_mark_diffbind.R) first. ",
         "It writes this file.")
  }

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")
  missing <- setdiff(MARK_TABLE_COLUMNS, colnames(df))
  if (length(missing) > 0) {
    stop("The gene-level all-marks table is missing columns: ",
         paste(missing, collapse = ", "), " in ", path,
         "\nSection 20_02 writes the full column contract.")
  }

  # This section tests chromatin only, so the methylation columns and the two
  # chromatin-state columns of the handoff (promoter_state and body_state) are
  # dropped here rather than carried unused.
  df <- df[, MARK_TABLE_COLUMNS, drop = FALSE]

  cat(sprintf("  Gene-level marks: %s rows\n", format(nrow(df), big.mark = ",")))
  for (mark in MARK_ORDER) {
    n_with <- sum(!is.na(df[[paste0(mark, "_fold")]]))
    cat(sprintf("    %-10s %s genes with a peak (%.1f%%)\n",
                mark_display(mark), format(n_with, big.mark = ","),
                100 * n_with / nrow(df)))
  }
  df
}

#' Collapse rows that repeat a gene symbol.
#'
#' Section 20_02 joins every mark onto the mCH table by gene symbol, so two
#' rows sharing a symbol carry the same mark values but different gene bodies.
#' The longest gene body is kept. The function stops when two rows sharing a
#' symbol disagree on a fold change, because the length rule cannot then decide
#' which chromatin measurement belongs to the symbol.
collapse_duplicate_symbols <- function(marks) {
  duplicated_names <- unique(marks$gene_name[duplicated(marks$gene_name)])
  if (length(duplicated_names) == 0) {
    cat("  No repeated gene symbol in the marks table\n")
    return(marks)
  }

  repeated <- marks[marks$gene_name %in% duplicated_names, , drop = FALSE]
  for (mark in MARK_ORDER) {
    fold_col <- paste0(mark, "_fold")
    n_distinct <- tapply(repeated[[fold_col]], repeated$gene_name,
                         function(v) length(unique(v)))
    conflicting <- names(n_distinct)[n_distinct > 1]
    if (length(conflicting) > 0) {
      stop(length(conflicting), " gene symbols carry more than one ",
           mark_display(mark), " fold change in the all-marks table. ",
           "First conflicts: ", paste(head(conflicting, 5), collapse = ", "),
           ". The symbol-level collapse cannot pick between them.")
    }
  }

  ordered <- marks[order(marks$gene_name, -marks$gene_length), , drop = FALSE]
  kept <- ordered[!duplicated(ordered$gene_name), , drop = FALSE]

  cat(sprintf("  %s gene symbols repeat; kept the longest gene body of each, ",
              format(length(duplicated_names), big.mark = ",")))
  cat(sprintf("%s rows dropped\n",
              format(nrow(marks) - nrow(kept), big.mark = ",")))
  kept
}

load_neuronal_set <- function() {
  path <- HANDOFF_PATHS$neuronal_gene_set
  if (!file.exists(path)) {
    stop("Neuronal gene set not found: ", path,
         "\nRun section 70_01 (70_01_k119ub_neuronal.R) first. ",
         "It writes this file.")
  }

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")
  missing <- setdiff(NEURONAL_SET_COLUMNS, colnames(df))
  if (length(missing) > 0) {
    stop("The neuronal gene set is missing columns: ",
         paste(missing, collapse = ", "), " in ", path,
         "\nSection 70_01 writes the full column contract.")
  }

  df$is_derived <- as.logical(df$is_derived)
  df$is_external <- as.logical(df$is_external)
  if (anyNA(df$is_derived) || anyNA(df$is_external)) {
    stop("is_derived or is_external is not readable as TRUE/FALSE in ", path)
  }
  if (anyDuplicated(df$gene) > 0) {
    stop("The neuronal gene set repeats a gene symbol in ", path)
  }

  cat(sprintf("  Neuronal gene set: %s genes (%s derived, %s external, %s both)\n",
              format(nrow(df), big.mark = ","),
              format(sum(df$source == "derived"), big.mark = ","),
              format(sum(df$source == "external"), big.mark = ","),
              format(sum(df$source == "both"), big.mark = ",")))
  df
}

#' Load the per-gene K119ub level and keep the genes quantifiable in both
#' conditions. The decile control needs this level for every gene, not only for
#' the neuronal ones the handoff carries.
load_k119ub_signal <- function() {
  path <- DIFFBIND_PATHS$k119ub_gene_signal
  if (!file.exists(path)) {
    stop("K119ub gene signal table not found: ", path,
         "\nRun scripts/utils/copy_reference_data.sh to place it under data/.")
  }

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")
  missing <- setdiff(SIGNAL_COLUMNS, colnames(df))
  if (length(missing) > 0) {
    stop("The K119ub gene signal table is missing columns: ",
         paste(missing, collapse = ", "), " in ", path)
  }
  if (anyDuplicated(df$symbol) > 0) {
    stop("The K119ub gene signal table repeats a gene symbol in ", path)
  }

  quantifiable <- df[df$gb_signal_class == QUANTIFIABLE_CLASS, , drop = FALSE]
  if (nrow(quantifiable) < 1000) {
    stop("Only ", nrow(quantifiable), " genes carry a quantifiable K119ub ",
         "gene-body signal. Too few to stratify the comparison by K119ub level.")
  }
  if (any(!is.finite(quantifiable$gb_ctrl_signal)) ||
      any(quantifiable$gb_ctrl_signal <= 0)) {
    stop("A quantifiable gene carries a control K119ub signal that is not a ",
         "positive number. The decile figure uses a log10 axis of it.")
  }

  cat(sprintf("  K119ub signal: %s genes, %s quantifiable\n",
              format(nrow(df), big.mark = ","),
              format(nrow(quantifiable), big.mark = ",")))
  quantifiable
}

# =============================================================================
# STEP 2: THE ANALYSIS UNIVERSE
# =============================================================================

#' Confirm that the handoff and the signal table report the same control K119ub
#' level for the genes they share.
check_signal_agreement <- function(universe, neuronal) {
  lookup <- neuronal[!is.na(neuronal$k119ub_ctrl_signal), , drop = FALSE]
  shared <- match(universe$gene_name, lookup$gene)
  has_both <- !is.na(shared)
  if (sum(has_both) == 0) {
    stop("No gene of the universe carries a K119ub level in both the neuronal ",
         "handoff and the K119ub signal table. Sections 70_01 and 70_02 are ",
         "reading different data.")
  }

  from_handoff <- lookup$k119ub_ctrl_signal[shared[has_both]]
  from_signal <- universe$k119ub_ctrl_signal[has_both]
  worst <- max(abs(from_handoff - from_signal))
  if (worst > SIGNAL_TOLERANCE) {
    stop("The neuronal handoff and the K119ub signal table disagree on the ",
         "control K119ub level by up to ", sprintf("%.3g", worst),
         " for the ", sum(has_both), " genes they share. Section 70_01 ran ",
         "against a different version of ", DIFFBIND_PATHS$k119ub_gene_signal)
  }
  cat(sprintf("  Control K119ub level agrees with the 70_01 handoff for %s ",
              format(sum(has_both), big.mark = ",")))
  cat(sprintf("shared genes (largest difference %.3g)\n", worst))
  invisible(worst)
}

#' Flag the marks that move in the predicted direction for every gene.
add_prediction_flags <- function(universe, fdr_threshold) {
  for (mark in MARK_ORDER) {
    fold <- universe[[paste0(mark, "_fold")]]
    fdr <- universe[[paste0(mark, "_fdr")]]
    moved <- if (MARK_META[[mark]]$predicted == "down") {
      !is.na(fold) & fold < 0
    } else {
      !is.na(fold) & fold > 0
    }
    universe[[paste0(mark, "_predicted")]] <- moved
    universe[[paste0(mark, "_predicted_sig")]] <- moved & !is.na(fdr) &
      fdr < fdr_threshold
  }

  fold_cols <- paste0(MARK_ORDER, "_fold")
  universe$n_marks_available <- as.integer(rowSums(!is.na(universe[, fold_cols])))
  universe$n_predicted <- as.integer(
    rowSums(universe[, paste0(MARK_ORDER, "_predicted")]))
  universe$n_predicted_sig <- as.integer(
    rowSums(universe[, paste0(MARK_ORDER, "_predicted_sig")]))
  universe
}

build_universe <- function(marks, neuronal, signal, fdr_threshold, out_dir) {
  cat("--- Step 2: building the analysis universe ---\n")

  signal_cols <- signal[, c("symbol", "gb_ctrl_signal", "gb_log2fc"), drop = FALSE]
  colnames(signal_cols) <- c("gene_name", "k119ub_ctrl_signal",
                             "k119ub_signal_log2fc")

  universe <- dplyr::inner_join(marks, signal_cols, by = "gene_name")
  if (anyDuplicated(universe$gene_name) > 0) {
    stop("The join with the K119ub signal table produced repeated gene ",
         "symbols. Both inputs must hold one row per symbol.")
  }
  if (nrow(universe) < 1000) {
    stop("Only ", nrow(universe), " genes carry both a mark fold change and a ",
         "quantifiable K119ub level. Too few for this section.")
  }

  cat(sprintf("  %s of %s marks-table genes also carry a K119ub level (%.1f%%)\n",
              format(nrow(universe), big.mark = ","),
              format(nrow(marks), big.mark = ","),
              100 * nrow(universe) / nrow(marks)))

  check_signal_agreement(universe, neuronal)

  derived_genes <- neuronal$gene[neuronal$is_derived]
  universe$is_neuronal <- universe$gene_name %in% neuronal$gene
  universe$is_neuronal_derived <- universe$gene_name %in% derived_genes
  universe$group <- factor(ifelse(universe$is_neuronal, "Neuronal", "Other"),
                           levels = GROUP_ORDER)
  universe$is_key_gene <- universe$gene_name %in% KEY_GENES

  if (sum(universe$is_neuronal) < 100 || sum(!universe$is_neuronal) < 100) {
    stop("The universe holds ", sum(universe$is_neuronal), " neuronal and ",
         sum(!universe$is_neuronal), " other genes. One class is too small to ",
         "compare.")
  }

  universe$k119ub_signal_decile <- cut_into_bins(universe$k119ub_ctrl_signal,
                                                 N_DECILES,
                                                 "control K119ub signal")
  universe$length_quintile <- cut_into_bins(universe$gene_length,
                                            N_LENGTH_QUINTILES, "gene length")
  universe$log10_k119ub_ctrl_signal <- log10(universe$k119ub_ctrl_signal)
  universe$log10_gene_length <- log10(universe$gene_length)

  universe <- add_prediction_flags(universe, fdr_threshold)

  n_complete <- sum(universe$n_marks_available == length(MARK_ORDER))
  cat(sprintf("  Neuronal: %s of %s genes (%.1f%%); %s of them are in the GO-derived set\n",
              format(sum(universe$is_neuronal), big.mark = ","),
              format(nrow(universe), big.mark = ","),
              100 * mean(universe$is_neuronal),
              format(sum(universe$is_neuronal_derived), big.mark = ",")))
  cat(sprintf("  Genes with all %d marks: %s (%.1f%%)\n",
              length(MARK_ORDER), format(n_complete, big.mark = ","),
              100 * n_complete / nrow(universe)))
  cat(sprintf("  Key genes present: %d of %d (%s)\n",
              sum(universe$is_key_gene), length(KEY_GENES),
              paste(universe$gene_name[universe$is_key_gene], collapse = ", ")))

  write_bin_bounds(universe, out_dir)
  cat("\n")
  universe
}

write_bin_bounds <- function(universe, out_dir) {
  decile_bounds <- universe %>%
    dplyr::group_by(k119ub_signal_decile) %>%
    dplyr::summarise(n_genes = dplyr::n(),
                     n_neuronal = sum(is_neuronal),
                     frac_neuronal = mean(is_neuronal),
                     signal_min = min(k119ub_ctrl_signal),
                     signal_max = max(k119ub_ctrl_signal),
                     .groups = "drop") %>%
    as.data.frame()
  write_section_table(decile_bounds, file.path(out_dir, "70_02c_decile_bounds.tsv"))

  quintile_bounds <- universe %>%
    dplyr::group_by(length_quintile) %>%
    dplyr::summarise(n_genes = dplyr::n(),
                     n_neuronal = sum(is_neuronal),
                     frac_neuronal = mean(is_neuronal),
                     length_min = min(gene_length),
                     length_max = max(gene_length),
                     .groups = "drop") %>%
    as.data.frame()
  write_section_table(quintile_bounds,
                  file.path(out_dir, "70_02d_quintile_bounds.tsv"))
}

write_universe_table <- function(universe, out_dir) {
  cols <- c("gene_name", "gene_id", "chr", "start", "end", "gene_length",
            "is_neuronal", "is_neuronal_derived", "is_key_gene",
            "k119ub_ctrl_signal", "k119ub_signal_log2fc",
            "k119ub_signal_decile", "length_quintile",
            unlist(lapply(MARK_ORDER, function(mark) {
              paste0(mark, c("_fold", "_fdr", "_n_peaks", "_predicted",
                             "_predicted_sig"))
            })),
            "n_marks_available", "n_predicted", "n_predicted_sig")
  out <- universe[order(universe$gene_name), cols]
  write_section_table(out, file.path(out_dir, "70_02_universe_per_gene.tsv"))
}

#' Long-format fold changes, one row per gene and mark.
build_fold_long <- function(universe) {
  displays <- mark_displays()
  long <- do.call(rbind, lapply(MARK_ORDER, function(mark) {
    data.frame(
      gene_name = universe$gene_name,
      mark = mark,
      mark_display = mark_display(mark),
      fold = universe[[paste0(mark, "_fold")]],
      group = universe$group,
      is_neuronal = universe$is_neuronal,
      k119ub_signal_decile = universe$k119ub_signal_decile,
      length_quintile = universe$length_quintile,
      k119ub_ctrl_signal = universe$k119ub_ctrl_signal,
      log10_k119ub_ctrl_signal = universe$log10_k119ub_ctrl_signal,
      is_key_gene = universe$is_key_gene,
      stringsAsFactors = FALSE
    )
  }))
  long <- long[!is.na(long$fold), , drop = FALSE]
  long$mark_display <- factor(long$mark_display, levels = unname(displays))
  long$group <- factor(long$group, levels = GROUP_ORDER)
  long
}

# =============================================================================
# FIGURE 70_02a: MARK FOLD CHANGE, NEURONAL AGAINST OTHER
# =============================================================================

#' Wilcoxon per mark for one definition of the neuronal class.
test_marks_for_gene_set <- function(universe, membership_col, gene_set_label) {
  do.call(rbind, lapply(MARK_ORDER, function(mark) {
    stats <- compare_gene_classes(universe[[paste0(mark, "_fold")]],
                                  universe[[membership_col]], with_ci = TRUE)
    stats$gene_set <- gene_set_label
    annotate_prediction(stats, mark)
  }))
}

run_mark_tests <- function(universe, out_dir) {
  cat("--- Figure 70_02a: mark fold change, neuronal against other ---\n")

  union_stats <- test_marks_for_gene_set(universe, "is_neuronal", "union")
  derived_stats <- test_marks_for_gene_set(universe, "is_neuronal_derived",
                                           "derived")
  stats <- rbind(union_stats, derived_stats)
  stats$q_value <- p.adjust(stats$p_value, method = "BH")

  ordered_cols <- c("gene_set", "mark", "mark_display", "predicted_direction",
                    "observed_direction", "matches_prediction",
                    "n_neuronal", "n_other",
                    "median_neuronal", "median_other", "delta_median",
                    "mean_neuronal", "mean_other",
                    "location_shift", "shift_ci_lower", "shift_ci_upper",
                    "rank_biserial", "wilcoxon_w", "p_value", "q_value")
  stats <- stats[, ordered_cols]
  write_section_table(stats, file.path(out_dir, "70_02a_mark_wilcoxon.tsv"))

  for (i in which(stats$gene_set == "union")) {
    r <- stats[i, ]
    cat(sprintf("  %-10s predicted %-4s | neuronal med %+.4f (n=%s) vs other %+.4f (n=%s) | ",
                r$mark_display, r$predicted_direction,
                r$median_neuronal, format(r$n_neuronal, big.mark = ","),
                r$median_other, format(r$n_other, big.mark = ",")))
    cat(sprintf("r = %+.3f, %s, %s\n", r$rank_biserial, format_p(r$p_value),
                ifelse(r$matches_prediction, "consistent", "opposite")))
  }
  stats
}

plot_mark_folds <- function(universe, mark_stats, out_dir) {
  fold_long <- build_fold_long(universe)

  union_stats <- mark_stats[mark_stats$gene_set == "union", , drop = FALSE]
  strip_labels <- setNames(
    sprintf("%s | predicted %s\nshift = %+.3f, r = %+.3f\n%s",
            union_stats$mark_display, union_stats$predicted_direction,
            union_stats$location_shift, union_stats$rank_biserial,
            format_p(union_stats$p_value)),
    union_stats$mark_display
  )
  fold_long$panel <- factor(strip_labels[as.character(fold_long$mark_display)],
                            levels = unname(strip_labels[levels(fold_long$mark_display)]))

  # mark_display is a one-to-one function of panel, so grouping by both gives
  # the same cells and carries the plain mark name into the summary table.
  labels <- summarise_cells(fold_long, c("mark_display", "panel", "group"),
                            "fold")
  spans <- fold_long %>%
    dplyr::group_by(panel) %>%
    dplyr::summarise(top = max(fold), span = diff(range(fold)), .groups = "drop") %>%
    as.data.frame()
  labels$label_y <- spans$top[match(labels$panel, spans$panel)] +
    0.10 * spans$span[match(labels$panel, spans$panel)]
  labels$label <- group_label(labels)

  figure <- ggplot(fold_long, aes(x = group, y = fold, fill = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_violin(scale = "width", alpha = 0.65, linewidth = 0.3,
                show.legend = FALSE) +
    geom_boxplot(width = 0.14, outlier.size = 0.2, outlier.alpha = 0.15,
                 show.legend = FALSE) +
    geom_text(data = labels, aes(x = group, y = label_y, label = label),
              inherit.aes = FALSE, size = 3.1, lineheight = 1.1) +
    facet_wrap(~ panel, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = GROUP_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(
      title = "Chromatin mark change in BAP1-KO: neuronal against other genes",
      subtitle = paste("DiffBind log2FC (mut / ctrl). The model predicts ATAC",
                       "down, H3K27ac down, H3K27me3 up, H2AK119ub up.",
                       "\nShift is the Hodges-Lehmann location shift,",
                       "neuronal minus other; r is the rank-biserial correlation."),
      x = "", y = "DiffBind log2FC (mut / ctrl)"
    ) +
    theme_emseq() +
    theme(strip.text = element_text(size = 8.5),
          plot.subtitle = element_text(size = 9, color = "grey40"))

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "70_02a_mark_fold_neuronal_vs_other"),
                          width = 16, height = 8)

  # panel holds the three-line facet strip text and label holds the on-plot
  # annotation. Both carry newlines, which write_section_table refuses, so the
  # written table names each mark with mark_display instead.
  panel_stats <- labels[, c("mark_display", "group", "n", "median",
                            "mean", "q25", "q75")]
  write_section_table(panel_stats,
                      file.path(out_dir, "70_02a_mark_fold_group_summary.tsv"))
  cat("\n")
  invisible(figure)
}

# =============================================================================
# FIGURE 70_02b: EFFECT SIZES ACROSS MARKS
# =============================================================================

#' Fit the adjusted models for one mark.
#'
#' The additive model asks whether the neuronal class still shifts the fold
#' change once the constitutive K119ub level and the gene length are held
#' fixed. The interaction model asks whether the neuronal class responds more
#' per unit of K119ub level.
fit_mark_models <- function(universe, mark) {
  fold_col <- paste0(mark, "_fold")
  df <- universe[!is.na(universe[[fold_col]]), , drop = FALSE]
  df$fold <- df[[fold_col]]

  additive_model <- lm(
    fold ~ log10_k119ub_ctrl_signal + log10_gene_length + is_neuronal,
    data = df)
  interaction_model <- lm(
    fold ~ log10_k119ub_ctrl_signal * is_neuronal + log10_gene_length,
    data = df)

  additive_summary <- summary(additive_model)
  interaction_summary <- summary(interaction_model)
  model_comparison <- anova(additive_model, interaction_model)

  additive_ci <- confint(additive_model)
  neuronal_term <- "is_neuronalTRUE"
  interaction_term <- "log10_k119ub_ctrl_signal:is_neuronalTRUE"
  if (!neuronal_term %in% rownames(additive_summary$coefficients)) {
    stop("The additive model for ", mark_display(mark), " has no ",
         neuronal_term, " coefficient. The neuronal class has one level only.")
  }
  if (!interaction_term %in% rownames(interaction_summary$coefficients)) {
    stop("The interaction model for ", mark_display(mark), " has no ",
         interaction_term, " coefficient.")
  }

  out <- data.frame(
    n_genes = nrow(df),
    adjusted_neuronal_beta = unname(additive_summary$coefficients[neuronal_term, 1]),
    adjusted_neuronal_se = unname(additive_summary$coefficients[neuronal_term, 2]),
    adjusted_neuronal_ci_lower = unname(additive_ci[neuronal_term, 1]),
    adjusted_neuronal_ci_upper = unname(additive_ci[neuronal_term, 2]),
    adjusted_neuronal_p = unname(additive_summary$coefficients[neuronal_term, 4]),
    k119ub_signal_beta = unname(additive_summary$coefficients["log10_k119ub_ctrl_signal", 1]),
    k119ub_signal_p = unname(additive_summary$coefficients["log10_k119ub_ctrl_signal", 4]),
    gene_length_beta = unname(additive_summary$coefficients["log10_gene_length", 1]),
    gene_length_p = unname(additive_summary$coefficients["log10_gene_length", 4]),
    interaction_beta = unname(interaction_summary$coefficients[interaction_term, 1]),
    interaction_p = unname(interaction_summary$coefficients[interaction_term, 4]),
    r2_additive = additive_summary$r.squared,
    r2_interaction = interaction_summary$r.squared,
    anova_p = model_comparison$`Pr(>F)`[2],
    stringsAsFactors = FALSE
  )
  out$delta_median <- out$adjusted_neuronal_beta
  out <- annotate_prediction(out, mark)
  out$delta_median <- NULL
  out
}

run_adjusted_models <- function(universe, out_dir) {
  cat("--- Adjusted models: mark fold ~ K119ub level + gene length + neuronal ---\n")

  models <- do.call(rbind, lapply(MARK_ORDER, function(mark) {
    fit_mark_models(universe, mark)
  }))
  models$adjusted_neuronal_q <- p.adjust(models$adjusted_neuronal_p,
                                         method = "BH")
  models$interaction_q <- p.adjust(models$interaction_p, method = "BH")

  ordered_cols <- c("mark", "mark_display", "predicted_direction",
                    "observed_direction", "matches_prediction", "n_genes",
                    "adjusted_neuronal_beta", "adjusted_neuronal_se",
                    "adjusted_neuronal_ci_lower", "adjusted_neuronal_ci_upper",
                    "adjusted_neuronal_p", "adjusted_neuronal_q",
                    "k119ub_signal_beta", "k119ub_signal_p",
                    "gene_length_beta", "gene_length_p",
                    "interaction_beta", "interaction_p", "interaction_q",
                    "r2_additive", "r2_interaction", "anova_p")
  models <- models[, ordered_cols]
  write_section_table(models, file.path(out_dir, "70_02b_adjusted_models.tsv"))

  for (i in seq_len(nrow(models))) {
    r <- models[i, ]
    cat(sprintf("  %-10s adjusted neuronal beta = %+.4f [%+.4f, %+.4f], %s (%s)\n",
                r$mark_display, r$adjusted_neuronal_beta,
                r$adjusted_neuronal_ci_lower, r$adjusted_neuronal_ci_upper,
                format_p(r$adjusted_neuronal_p),
                ifelse(r$matches_prediction, "consistent", "opposite")))
    cat(sprintf("             K119ub level beta = %+.4f, gene length beta = %+.4f, ",
                r$k119ub_signal_beta, r$gene_length_beta))
    cat(sprintf("interaction beta = %+.4f (%s)\n", r$interaction_beta,
                format_p(r$interaction_p)))
  }
  cat("\n")
  models
}

#' One forest panel: an estimate with its interval, one row per mark.
build_forest_panel <- function(df, estimate_col, lower_col, upper_col, p_col,
                               panel_title, x_title) {
  displays <- unname(mark_displays())
  df$mark_display <- factor(df$mark_display, levels = rev(displays))
  df$prediction <- factor(ifelse(df$matches_prediction, "Consistent", "Opposite"),
                          levels = names(PREDICTION_COLORS))
  df$point_label <- sprintf("%+.4f [%+.4f, %+.4f]  %s",
                            df[[estimate_col]], df[[lower_col]],
                            df[[upper_col]], format_p(df[[p_col]]))

  ggplot(df, aes(x = .data[[estimate_col]], y = mark_display)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_errorbar(aes(xmin = .data[[lower_col]], xmax = .data[[upper_col]]),
                  orientation = "y", width = 0.2, linewidth = 0.8) +
    geom_point(aes(color = prediction), size = 4.5) +
    geom_text(aes(label = point_label), vjust = -1.2, size = 3.1,
              color = "grey20") +
    scale_color_manual(values = PREDICTION_COLORS,
                       name = "Sign against prediction", drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.25, 0.25))) +
    scale_y_discrete(expand = expansion(add = 0.7)) +
    labs(title = panel_title, x = x_title, y = "") +
    theme_emseq() +
    theme(axis.text.y = element_text(size = 11, face = "bold"),
          legend.position = "top")
}

plot_effect_forest <- function(mark_stats, models, out_dir) {
  cat("--- Figure 70_02b: effect sizes across marks ---\n")

  unadjusted <- mark_stats[mark_stats$gene_set == "union", , drop = FALSE]

  p_unadjusted <- build_forest_panel(
    unadjusted, "location_shift", "shift_ci_lower", "shift_ci_upper", "p_value",
    "Unadjusted: Hodges-Lehmann shift, neuronal minus other",
    "Shift in DiffBind log2FC"
  )

  p_adjusted <- build_forest_panel(
    models, "adjusted_neuronal_beta", "adjusted_neuronal_ci_lower",
    "adjusted_neuronal_ci_upper", "adjusted_neuronal_p",
    "Adjusted: neuronal coefficient with K119ub level and gene length held fixed",
    "Coefficient in DiffBind log2FC"
  )

  figure <- p_unadjusted / p_adjusted +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Effect of neuronal identity on chromatin remodeling, by mark",
      subtitle = paste("Consistent means the sign agrees with the predicted",
                       "direction: ATAC down, H3K27ac down, H3K27me3 up,",
                       "H2AK119ub up"),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = 16),
                    plot.subtitle = element_text(hjust = 0.5, size = 10,
                                                 color = "grey40"),
                    legend.position = "top")
    )

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "70_02b_effect_size_forest"),
                          width = 13, height = 10)

  forest_numbers <- rbind(
    data.frame(panel = "unadjusted", mark = unadjusted$mark,
               mark_display = unadjusted$mark_display,
               predicted_direction = unadjusted$predicted_direction,
               estimate = unadjusted$location_shift,
               ci_lower = unadjusted$shift_ci_lower,
               ci_upper = unadjusted$shift_ci_upper,
               p_value = unadjusted$p_value,
               matches_prediction = unadjusted$matches_prediction,
               stringsAsFactors = FALSE),
    data.frame(panel = "adjusted", mark = models$mark,
               mark_display = models$mark_display,
               predicted_direction = models$predicted_direction,
               estimate = models$adjusted_neuronal_beta,
               ci_lower = models$adjusted_neuronal_ci_lower,
               ci_upper = models$adjusted_neuronal_ci_upper,
               p_value = models$adjusted_neuronal_p,
               matches_prediction = models$matches_prediction,
               stringsAsFactors = FALSE)
  )
  write_section_table(forest_numbers,
                  file.path(out_dir, "70_02b_effect_size_forest.tsv"))
  cat("\n")
  invisible(figure)
}

# =============================================================================
# STRATIFIED COMPARISONS
# =============================================================================

#' Run the neuronal against other comparison inside every level of a stratum.
#'
#' @param universe The analysis universe.
#' @param stratum_col Column holding the integer stratum of each gene.
#' @param n_strata Number of strata.
#' @param min_group_size Smallest class allowed in a cell.
#' @param stratum_name Word used in the stop message and the output table.
#' @return data.frame with one row per mark and stratum
test_within_strata <- function(universe, stratum_col, n_strata, min_group_size,
                               stratum_name) {
  rows <- list()
  for (mark in MARK_ORDER) {
    fold_col <- paste0(mark, "_fold")
    for (level in seq_len(n_strata)) {
      in_stratum <- universe[[stratum_col]] == level
      values <- universe[[fold_col]][in_stratum]
      neuronal <- universe$is_neuronal[in_stratum]

      require_group_sizes(
        sum(neuronal & !is.na(values)), sum(!neuronal & !is.na(values)),
        min_group_size,
        sprintf("%s %d for %s", stratum_name, level, mark_display(mark))
      )

      stats <- compare_gene_classes(values, neuronal, with_ci = FALSE)
      stats$stratum <- level
      rows[[length(rows) + 1]] <- annotate_prediction(stats, mark)
    }
  }
  out <- do.call(rbind, rows)
  out$q_value <- p.adjust(out$p_value, method = "BH")
  out$significance <- ifelse(out$q_value < Q_THRESHOLD, "Significant",
                             "Not Significant")
  out
}

#' Per-cell n, median and quartiles for a stratified figure.
summarise_within_strata <- function(fold_long, stratum_col) {
  summary_df <- summarise_cells(fold_long,
                                c("mark_display", stratum_col, "group"), "fold")
  summary_df
}

# =============================================================================
# FIGURE 70_02c: THE COMPARISON INSIDE K119ub SIGNAL DECILES
# =============================================================================

plot_decile_stratified <- function(universe, min_group_size, out_dir) {
  cat("--- Figure 70_02c: the comparison inside K119ub signal deciles ---\n")

  tests <- test_within_strata(universe, "k119ub_signal_decile", N_DECILES,
                              min_group_size, "K119ub signal decile")
  write_section_table(tests, file.path(out_dir, "70_02c_decile_tests.tsv"))

  fold_long <- build_fold_long(universe)
  summary_df <- summarise_within_strata(fold_long, "k119ub_signal_decile")
  # summary_df holds data only. The three panels below draw the median, the
  # interquartile range and n from its columns, so no label text is needed.
  write_section_table(summary_df,
                      file.path(out_dir, "70_02c_decile_group_summary.tsv"))

  # A numeric x keeps the ribbon, the bars and the counts on one shared axis
  # across the three stacked panels.
  decile_labels <- sprintf("D%d", seq_len(N_DECILES))
  decile_axis <- scale_x_continuous(breaks = seq_len(N_DECILES),
                                    labels = decile_labels)
  summary_df$decile <- summary_df$k119ub_signal_decile
  tests$decile <- tests$stratum
  tests$mark_display <- factor(tests$mark_display,
                               levels = unname(mark_displays()))
  summary_df$group <- factor(summary_df$group, levels = GROUP_ORDER)

  p_median <- ggplot(summary_df, aes(x = decile, y = median,
                                     color = group, group = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = group), alpha = 0.15,
                color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(~ mark_display, nrow = 1, scales = "free_y") +
    decile_axis +
    scale_color_manual(values = GROUP_COLORS, name = "Gene class") +
    scale_fill_manual(values = GROUP_COLORS, guide = "none") +
    labs(
      title = "Median mark fold change across control K119ub signal deciles",
      subtitle = "Band shows the interquartile range of each class in each decile",
      x = "", y = "Median log2FC (mut / ctrl)"
    ) +
    theme_emseq(base_size = 11) +
    theme(axis.text.x = element_text(size = 7),
          plot.subtitle = element_text(size = 9, color = "grey40"))

  p_delta <- ggplot(tests, aes(x = decile, y = delta_median,
                               fill = significance)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%+.2f", rank_biserial),
                  vjust = ifelse(delta_median >= 0, -0.5, 1.3)),
              size = 2.2, color = "grey25") +
    facet_wrap(~ mark_display, nrow = 1, scales = "free_y") +
    decile_axis +
    scale_fill_manual(values = COLORS$significant, name = "BH q < 0.05") +
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
    labs(
      title = "Median difference, neuronal minus other, inside each decile",
      subtitle = "Number above each bar is the rank-biserial correlation of that cell",
      x = "", y = "Difference in median log2FC"
    ) +
    theme_emseq(base_size = 11) +
    theme(axis.text.x = element_text(size = 7),
          plot.subtitle = element_text(size = 9, color = "grey40"))

  p_counts <- ggplot(summary_df, aes(x = decile, y = n, fill = group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_text(aes(label = format(n, big.mark = ",", trim = TRUE)),
              position = position_dodge(width = 0.8), vjust = -0.4, size = 2.1) +
    facet_wrap(~ mark_display, nrow = 1) +
    decile_axis +
    scale_fill_manual(values = GROUP_COLORS, name = "Gene class") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "Genes with a measured fold change in each decile",
         x = "Control K119ub signal decile (D1 lowest, D10 highest)",
         y = "Genes") +
    theme_emseq(base_size = 11) +
    theme(axis.text.x = element_text(size = 7))

  figure <- p_median / p_delta / p_counts +
    plot_layout(heights = c(1, 0.9, 0.7)) +
    plot_annotation(
      title = "Chromatin remodeling at neuronal genes inside K119ub level strata",
      subtitle = paste("Equal-count deciles of the control K119ub gene-body",
                       "signal, so the K119ub level cannot drive a",
                       "within-decile difference"),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = 16),
                    plot.subtitle = element_text(hjust = 0.5, size = 10,
                                                 color = "grey40"))
    )

  save_multiformat_ggplot(figure, file.path(out_dir, "70_02c_decile_stratified"),
                          width = 18, height = 13)

  for (mark in MARK_ORDER) {
    rows <- tests[tests$mark == mark, ]
    n_consistent <- sum(rows$matches_prediction & rows$q_value < Q_THRESHOLD)
    cat(sprintf("  %-10s %d of %d deciles consistent with the prediction at q < %.2f",
                mark_display(mark), n_consistent, N_DECILES, Q_THRESHOLD))
    cat(sprintf(" | D1 delta %+.4f, D10 delta %+.4f\n",
                rows$delta_median[rows$stratum == 1],
                rows$delta_median[rows$stratum == N_DECILES]))
  }
  cat("\n")
  tests
}

# =============================================================================
# FIGURE 70_02d: THE COMPARISON INSIDE GENE-LENGTH QUINTILES
# =============================================================================

#' Axis labels naming the gene-length range of each quintile.
build_quintile_labels <- function(universe) {
  bounds <- universe %>%
    dplyr::group_by(length_quintile) %>%
    dplyr::summarise(length_min = min(gene_length), length_max = max(gene_length),
                     .groups = "drop") %>%
    as.data.frame()
  labels <- sprintf("L%d\n%s-%s kb", bounds$length_quintile,
                    format(round(bounds$length_min / 1000, 1), trim = TRUE),
                    format(round(bounds$length_max / 1000, 1), trim = TRUE))
  names(labels) <- as.character(bounds$length_quintile)
  labels
}

plot_quintile_stratified <- function(universe, min_group_size, out_dir) {
  cat("--- Figure 70_02d: the comparison inside gene-length quintiles ---\n")

  tests <- test_within_strata(universe, "length_quintile", N_LENGTH_QUINTILES,
                              min_group_size, "gene-length quintile")
  write_section_table(tests, file.path(out_dir, "70_02d_quintile_tests.tsv"))

  survival <- do.call(rbind, lapply(MARK_ORDER, function(mark) {
    rows <- tests[tests$mark == mark, , drop = FALSE]
    consistent <- rows$matches_prediction & rows$q_value < Q_THRESHOLD
    data.frame(
      mark = mark,
      mark_display = mark_display(mark),
      predicted_direction = MARK_META[[mark]]$predicted,
      n_quintiles = nrow(rows),
      n_consistent_and_significant = sum(consistent),
      n_consistent_sign = sum(rows$matches_prediction),
      n_significant = sum(rows$q_value < Q_THRESHOLD),
      min_delta_median = min(rows$delta_median),
      max_delta_median = max(rows$delta_median),
      effect_survives_length_control = all(consistent),
      stringsAsFactors = FALSE
    )
  }))
  write_section_table(survival, file.path(out_dir, "70_02d_quintile_survival.tsv"))

  fold_long <- build_fold_long(universe)
  summary_df <- summarise_within_strata(fold_long, "length_quintile")
  # The written table holds data only. The violin annotation below adds the
  # label column after the write.
  write_section_table(summary_df,
                      file.path(out_dir, "70_02d_quintile_group_summary.tsv"))
  summary_df$label <- group_label(summary_df)

  quintile_labels <- build_quintile_labels(universe)
  fold_long$quintile_label <- factor(
    quintile_labels[as.character(fold_long$length_quintile)],
    levels = unname(quintile_labels)
  )
  summary_df$quintile_label <- factor(
    quintile_labels[as.character(summary_df$length_quintile)],
    levels = unname(quintile_labels)
  )
  tests$quintile_label <- factor(
    quintile_labels[as.character(tests$stratum)],
    levels = unname(quintile_labels)
  )
  tests$mark_display <- factor(tests$mark_display,
                               levels = unname(mark_displays()))

  spans <- fold_long %>%
    dplyr::group_by(mark_display) %>%
    dplyr::summarise(top = max(fold), span = diff(range(fold)), .groups = "drop") %>%
    as.data.frame()
  summary_df$label_y <- spans$top[match(summary_df$mark_display,
                                        spans$mark_display)] +
    0.12 * spans$span[match(summary_df$mark_display, spans$mark_display)]

  p_violin <- ggplot(fold_long, aes(x = group, y = fold, fill = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_violin(scale = "width", alpha = 0.65, linewidth = 0.25,
                show.legend = FALSE) +
    geom_boxplot(width = 0.14, outlier.size = 0.1, outlier.alpha = 0.1,
                 show.legend = FALSE) +
    geom_text(data = summary_df,
              aes(x = group, y = label_y, label = label),
              inherit.aes = FALSE, size = 2.2, lineheight = 1.05) +
    facet_grid(mark_display ~ quintile_label, scales = "free_y") +
    scale_fill_manual(values = GROUP_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    labs(
      title = "Mark fold change inside gene-length quintiles",
      subtitle = "Equal-count quintiles of gene length, so length cannot drive a within-panel difference",
      x = "", y = "DiffBind log2FC (mut / ctrl)"
    ) +
    theme_emseq(base_size = 11) +
    theme(axis.text.x = element_text(size = 8),
          strip.text.x = element_text(size = 8),
          strip.text.y = element_text(size = 9),
          plot.subtitle = element_text(size = 9, color = "grey40"))

  p_effect <- ggplot(tests, aes(x = quintile_label, y = rank_biserial,
                                group = 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_line(color = "grey45", linewidth = 0.6) +
    geom_point(aes(color = significance), size = 3) +
    geom_text(aes(label = sprintf("%+.3f\n%s", rank_biserial,
                                  format_p(p_value))),
              vjust = -0.6, size = 2.3, lineheight = 1.05, color = "grey25") +
    facet_wrap(~ mark_display, nrow = 1) +
    scale_color_manual(values = COLORS$significant, name = "BH q < 0.05") +
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.4))) +
    labs(
      title = "Rank-biserial correlation, neuronal against other, inside each quintile",
      subtitle = "Positive means neuronal genes carry the higher log2FC in that quintile",
      x = "Gene-length quintile", y = "Rank-biserial correlation"
    ) +
    theme_emseq(base_size = 11) +
    theme(axis.text.x = element_text(size = 7),
          legend.position = "top",
          plot.subtitle = element_text(size = 9, color = "grey40"))

  figure <- p_violin / p_effect +
    plot_layout(heights = c(1, 0.45)) +
    plot_annotation(
      title = "Gene length does not explain the neuronal chromatin response",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = 16))
    )

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "70_02d_length_quintile_stratified"),
                          width = 18, height = 15)

  for (i in seq_len(nrow(survival))) {
    r <- survival[i, ]
    cat(sprintf("  %-10s %d of %d quintiles consistent and significant | survives: %s\n",
                r$mark_display, r$n_consistent_and_significant, r$n_quintiles,
                ifelse(r$effect_survives_length_control, "yes", "no")))
  }
  cat("\n")
  survival
}

# =============================================================================
# GENE-LEVEL FISHER TESTS
# =============================================================================

#' Register one Fisher test of neuronal membership against a logical gene flag.
#'
#' Figure 70_02e uses it for the concordance thresholds and the per-mark tests
#' use it for the direction of a single mark.
#'
#' @param df Genes to test. Needs gene_name, chr, is_neuronal and flag_col.
#' @param flag_col Name of the logical column forming the second variable.
#' @param test_id Identifier of this test inside section 70_02.
#' @param description What the test asks, in one sentence.
#' @param out_dir Section output directory.
#' @return one-row data.frame with the counts and the odds ratio
test_neuronal_flag <- function(df, flag_col, test_id, description, out_dir) {
  gene_df <- data.frame(
    gene_name = df$gene_name,
    chr = df$chr,
    is_neuronal = df$is_neuronal,
    stringsAsFactors = FALSE
  )
  gene_df[[flag_col]] <- df[[flag_col]]

  ft <- register_fisher_test(
    section = SECTION_ID, test_id = test_id, description = description,
    gene_df = gene_df, row_var = "is_neuronal", col_var = flag_col,
    output_dir = out_dir
  )

  data.frame(
    test_id = test_id,
    description = description,
    n_genes = nrow(gene_df),
    n_neuronal = sum(gene_df$is_neuronal),
    n_flagged = sum(gene_df[[flag_col]]),
    n_neuronal_flagged = sum(gene_df$is_neuronal & gene_df[[flag_col]]),
    frac_flagged_neuronal = mean(gene_df[[flag_col]][gene_df$is_neuronal]),
    frac_flagged_other = mean(gene_df[[flag_col]][!gene_df$is_neuronal]),
    odds_ratio = unname(ft$estimate),
    ci_lower = ft$conf.int[1],
    ci_upper = ft$conf.int[2],
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# FIGURE 70_02e: COUNT OF MARKS MOVING IN THE PREDICTED DIRECTION
# =============================================================================

#' Distribution of a per-gene count across the two gene classes.
#'
#' Every count from 0 to the number of marks appears for both classes, with a
#' zero where no gene carries that count, so the two bars of a count always
#' stand side by side.
summarise_count_distribution <- function(df, count_col, count_label) {
  observed <- df %>%
    dplyr::mutate(n_predicted = .data[[count_col]]) %>%
    dplyr::count(group, n_predicted, name = "n_genes") %>%
    as.data.frame()

  totals <- df %>%
    dplyr::count(group, name = "n_group") %>%
    as.data.frame()

  grid <- expand.grid(
    group = factor(GROUP_ORDER, levels = GROUP_ORDER),
    n_predicted = seq(0, length(MARK_ORDER)),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )

  distribution <- dplyr::left_join(grid, observed, by = c("group", "n_predicted"))
  distribution$n_genes[is.na(distribution$n_genes)] <- 0L
  distribution <- dplyr::left_join(distribution, totals, by = "group")
  distribution$pct <- 100 * distribution$n_genes / distribution$n_group
  distribution$count_type <- count_label
  distribution
}

plot_concordance <- function(universe, fdr_threshold, out_dir) {
  cat("--- Figure 70_02e: marks moving in the predicted direction per gene ---\n")

  complete <- universe[universe$n_marks_available == length(MARK_ORDER), ,
                       drop = FALSE]
  if (nrow(complete) < 100) {
    stop("Only ", nrow(complete), " genes carry a fold change for all four ",
         "marks. The concordance count needs the same four marks in every gene.")
  }
  cat(sprintf("  %s genes carry all four marks (%s neuronal)\n",
              format(nrow(complete), big.mark = ","),
              format(sum(complete$is_neuronal), big.mark = ",")))

  sign_test <- compare_gene_classes(complete$n_predicted, complete$is_neuronal,
                                    with_ci = FALSE)
  sig_test <- compare_gene_classes(complete$n_predicted_sig,
                                   complete$is_neuronal, with_ci = FALSE)
  sign_test$count_type <- "Direction of the fold change"
  sig_test$count_type <- sprintf("Direction and FDR < %.3g", fdr_threshold)
  count_tests <- rbind(sign_test, sig_test)
  write_section_table(count_tests,
                  file.path(out_dir, "70_02e_concordance_wilcoxon.tsv"))

  cat(sprintf("  Direction only: neuronal median %.2f (mean %.3f) vs other %.2f (mean %.3f), %s\n",
              sign_test$median_neuronal, sign_test$mean_neuronal,
              sign_test$median_other, sign_test$mean_other,
              format_p(sign_test$p_value)))
  cat(sprintf("  Direction and FDR: neuronal median %.2f (mean %.3f) vs other %.2f (mean %.3f), %s\n",
              sig_test$median_neuronal, sig_test$mean_neuronal,
              sig_test$median_other, sig_test$mean_other,
              format_p(sig_test$p_value)))

  complete$all_four_predicted <- complete$n_predicted == length(MARK_ORDER)
  complete$three_plus_predicted <- complete$n_predicted >= 3
  complete$two_plus_predicted_sig <- complete$n_predicted_sig >= 2

  fisher_rows <- rbind(
    test_neuronal_flag(
      complete, "all_four_predicted", "all_four_predicted",
      paste("Among genes with all four marks measured, neuronal genes carry",
            "all four marks moving in the predicted direction."),
      out_dir),
    test_neuronal_flag(
      complete, "three_plus_predicted", "three_plus_predicted",
      paste("Among genes with all four marks measured, neuronal genes carry",
            "three or more marks moving in the predicted direction."),
      out_dir),
    test_neuronal_flag(
      complete, "two_plus_predicted_sig", "two_plus_predicted_sig",
      paste("Among genes with all four marks measured, neuronal genes carry",
            "two or more marks moving in the predicted direction at a",
            "significant FDR."),
      out_dir)
  )
  write_section_table(fisher_rows,
                  file.path(out_dir, "70_02e_concordance_fisher.tsv"))

  distribution <- rbind(
    summarise_count_distribution(complete, "n_predicted",
                                 "Direction of the fold change"),
    summarise_count_distribution(complete, "n_predicted_sig",
                                 sprintf("Direction and FDR < %.3g",
                                         fdr_threshold))
  )
  write_section_table(distribution,
                  file.path(out_dir, "70_02e_concordance_distribution.tsv"))

  axis_labels <- setNames(
    sprintf("%s\n(n = %s)", GROUP_ORDER,
            format(c(sum(complete$is_neuronal), sum(!complete$is_neuronal)),
                   big.mark = ",", trim = TRUE)),
    GROUP_ORDER
  )

  build_count_panel <- function(count_type, test_row, panel_title) {
    df <- distribution[distribution$count_type == count_type, , drop = FALSE]
    df$n_predicted <- factor(df$n_predicted,
                             levels = sort(unique(df$n_predicted)))
    df$group <- factor(df$group, levels = GROUP_ORDER)

    ggplot(df, aes(x = n_predicted, y = pct, fill = group)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.75) +
      geom_text(aes(label = sprintf("%.1f%%\nn = %s", pct,
                                    format(n_genes, big.mark = ",",
                                           trim = TRUE))),
                position = position_dodge(width = 0.8), vjust = -0.3,
                size = 2.7, lineheight = 1.05) +
      scale_fill_manual(values = GROUP_COLORS, name = "Gene class",
                        labels = axis_labels) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
      labs(
        title = panel_title,
        subtitle = sprintf(paste("Median: neuronal %.2f, other %.2f |",
                                 "mean: neuronal %.3f, other %.3f | %s"),
                           test_row$median_neuronal, test_row$median_other,
                           test_row$mean_neuronal, test_row$mean_other,
                           format_p(test_row$p_value)),
        x = "Marks moving in the predicted direction", y = "Percent of genes"
      ) +
      theme_emseq(base_size = 11) +
      theme(plot.subtitle = element_text(size = 9, color = "grey40"))
  }

  p_sign <- build_count_panel(sign_test$count_type, sign_test,
                              "Counted by the direction of the fold change")
  p_sig <- build_count_panel(sig_test$count_type, sig_test,
                             sprintf("Counted by direction and FDR < %.3g",
                                     fdr_threshold))

  all_four <- fisher_rows[fisher_rows$test_id == "all_four_predicted", ]
  figure <- p_sign / p_sig +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Marks moving in the predicted direction per gene",
      subtitle = sprintf(paste("Genes with all four marks measured (n = %s).",
                               "All four predicted: %.1f%% of neuronal genes",
                               "against %.1f%% of other genes, OR = %.2f",
                               "[%.2f, %.2f], %s"),
                         format(nrow(complete), big.mark = ","),
                         100 * all_four$frac_flagged_neuronal,
                         100 * all_four$frac_flagged_other,
                         all_four$odds_ratio, all_four$ci_lower,
                         all_four$ci_upper, format_p(all_four$p_value)),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = 16),
                    plot.subtitle = element_text(hjust = 0.5, size = 10,
                                                 color = "grey40"))
    )

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "70_02e_concordant_mark_counts"),
                          width = 14, height = 11)
  cat("\n")
  list(count_tests = count_tests, fisher = fisher_rows)
}

# =============================================================================
# PER-MARK DIRECTION FISHER TESTS
# =============================================================================

run_mark_direction_fishers <- function(universe, out_dir) {
  cat("--- Fisher tests: mark direction against neuronal identity ---\n")

  rows <- lapply(MARK_ORDER, function(mark) {
    fold_col <- paste0(mark, "_fold")
    df <- universe[!is.na(universe[[fold_col]]), , drop = FALSE]
    flag_col <- paste0(mark, "_predicted")

    row <- test_neuronal_flag(
      df, flag_col, paste0(mark, "_predicted"),
      sprintf(paste("Among genes with a %s peak, the predicted %s log2 fold",
                    "change is more common at neuronal genes."),
              mark_display(mark), MARK_META[[mark]]$predicted),
      out_dir
    )
    row$mark <- mark
    row$mark_display <- mark_display(mark)
    row$predicted_direction <- MARK_META[[mark]]$predicted
    row
  })

  results <- do.call(rbind, rows)
  results$q_value <- p.adjust(results$p_value, method = "BH")
  results <- results[, c("mark", "mark_display", "predicted_direction",
                         setdiff(colnames(results),
                                 c("mark", "mark_display",
                                   "predicted_direction")))]
  write_section_table(results,
                  file.path(out_dir, "70_02_mark_direction_fisher.tsv"))
  cat("\n")
  results
}

# =============================================================================
# FIGURE 70_02f: K119ub LEVEL AGAINST MARK FOLD CHANGE
# =============================================================================

plot_signal_scatter <- function(universe, out_dir) {
  cat("--- Figure 70_02f: K119ub level against mark fold change ---\n")

  fold_long <- build_fold_long(universe)

  annotations <- do.call(rbind, lapply(MARK_ORDER, function(mark) {
    df <- fold_long[fold_long$mark == mark, , drop = FALSE]
    rows <- lapply(GROUP_ORDER, function(grp) {
      sub <- df[df$group == grp, , drop = FALSE]
      test <- cor.test(sub$log10_k119ub_ctrl_signal, sub$fold,
                       method = "spearman", exact = FALSE)
      data.frame(mark = mark, mark_display = mark_display(mark), group = grp,
                 n_genes = nrow(sub), spearman_rho = unname(test$estimate),
                 p_value = test$p.value, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }))
  write_section_table(annotations,
                  file.path(out_dir, "70_02f_signal_vs_fold_correlation.tsv"))

  panel_text <- annotations %>%
    dplyr::group_by(mark_display) %>%
    dplyr::summarise(label = paste(sprintf("%s: rho = %+.3f (n = %s)",
                                           group, spearman_rho,
                                           format(n_genes, big.mark = ",",
                                                  trim = TRUE)),
                                   collapse = "\n"),
                     .groups = "drop") %>%
    as.data.frame()
  panel_text$mark_display <- factor(panel_text$mark_display,
                                    levels = unname(mark_displays()))
  # A log10 x axis cannot place text at -Inf, so the annotation sits at the
  # smallest signal of the universe.
  panel_text$x_pos <- min(fold_long$k119ub_ctrl_signal)

  key_df <- fold_long[fold_long$is_key_gene, , drop = FALSE]

  figure <- ggplot(fold_long, aes(x = k119ub_ctrl_signal, y = fold)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55",
               linewidth = 0.3) +
    geom_point(aes(color = group), size = 0.35, alpha = 0.25) +
    geom_smooth(aes(color = group), method = "loess", formula = y ~ x,
                se = FALSE, linewidth = 0.9) +
    geom_text(data = panel_text, aes(x = x_pos, y = Inf, label = label),
              inherit.aes = FALSE, hjust = 0, vjust = 1.2, size = 2.9,
              lineheight = 1.1, color = "grey20") +
    facet_wrap(~ mark_display, nrow = 1, scales = "free_y") +
    scale_color_manual(values = GROUP_COLORS, name = "Gene class",
                       guide = guide_legend(override.aes = list(size = 3,
                                                                alpha = 1))) +
    scale_x_log10() +
    labs(
      title = "Constitutive K119ub level against the chromatin response",
      subtitle = paste("Loess fit per gene class. The decile panel tests the",
                       "same relationship inside equal-count strata."),
      x = "Control K119ub gene-body signal (log10 scale)",
      y = "DiffBind log2FC (mut / ctrl)"
    ) +
    theme_emseq(base_size = 11) +
    theme(legend.position = "top",
          plot.subtitle = element_text(size = 9, color = "grey40"))

  if (nrow(key_df) > 0) {
    figure <- figure +
      geom_point(data = key_df, size = 1.4, color = "black") +
      ggrepel::geom_text_repel(data = key_df,
                               aes(label = gene_name), size = 2.4,
                               min.segment.length = 0, max.overlaps = Inf,
                               color = "black")
  }

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "70_02f_k119ub_signal_vs_fold"),
                          width = 18, height = 6.5)

  key_out <- key_df[, c("gene_name", "mark_display", "fold", "group",
                        "k119ub_ctrl_signal", "k119ub_signal_decile",
                        "length_quintile")]
  write_section_table(key_out, file.path(out_dir, "70_02f_key_genes.tsv"))
  cat("\n")
  invisible(annotations)
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opt <- parse_section_args()
  out_dir <- opt$output_dir
  fdr_threshold <- opt$fdr_threshold
  min_group_size <- opt$min_group_size

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("===========================================================================\n")
  cat("  SECTION 70_02: CHROMATIN REMODELING AT NEURONAL GENES\n")
  cat("===========================================================================\n")
  cat("  Output directory:   ", out_dir, "\n", sep = "")
  cat("  Peak FDR threshold: ", fdr_threshold, "\n", sep = "")
  cat("  Minimum group size: ", min_group_size, "\n", sep = "")
  cat("  Predicted directions: ")
  cat(paste(sprintf("%s %s", unname(mark_displays()),
                    vapply(MARK_META[MARK_ORDER], `[[`, character(1),
                           "predicted")),
            collapse = ", "), "\n\n", sep = "")

  cat("--- Step 1: reading the input tables ---\n")
  marks <- load_mark_table()
  marks <- collapse_duplicate_symbols(marks)
  neuronal <- load_neuronal_set()
  signal <- load_k119ub_signal()
  cat("\n")

  universe <- build_universe(marks, neuronal, signal, fdr_threshold, out_dir)
  write_universe_table(universe, out_dir)
  cat("\n")

  mark_stats <- run_mark_tests(universe, out_dir)
  plot_mark_folds(universe, mark_stats, out_dir)

  models <- run_adjusted_models(universe, out_dir)
  plot_effect_forest(mark_stats, models, out_dir)

  decile_tests <- plot_decile_stratified(universe, min_group_size, out_dir)
  quintile_survival <- plot_quintile_stratified(universe, min_group_size,
                                                out_dir)
  concordance <- plot_concordance(universe, fdr_threshold, out_dir)
  direction_fishers <- run_mark_direction_fishers(universe, out_dir)
  plot_signal_scatter(universe, out_dir)

  cat("---------------------------------------------------------------------------\n")
  cat("  SUMMARY\n")
  cat("---------------------------------------------------------------------------\n")
  cat(sprintf("  Universe: %s genes with mark folds and a K119ub level, %s neuronal (%.1f%%)\n",
              format(nrow(universe), big.mark = ","),
              format(sum(universe$is_neuronal), big.mark = ","),
              100 * mean(universe$is_neuronal)))

  cat("  Unadjusted Wilcoxon per mark (union neuronal set):\n")
  union_stats <- mark_stats[mark_stats$gene_set == "union", ]
  for (i in seq_len(nrow(union_stats))) {
    r <- union_stats[i, ]
    cat(sprintf("    %-10s predicted %-4s | shift %+.4f [%+.4f, %+.4f] | r %+.3f | %s | %s\n",
                r$mark_display, r$predicted_direction, r$location_shift,
                r$shift_ci_lower, r$shift_ci_upper, r$rank_biserial,
                format_p(r$q_value),
                ifelse(r$matches_prediction, "consistent", "opposite")))
  }

  cat("  Adjusted for K119ub level and gene length:\n")
  for (i in seq_len(nrow(models))) {
    r <- models[i, ]
    cat(sprintf("    %-10s beta %+.4f [%+.4f, %+.4f] | %s | %s\n",
                r$mark_display, r$adjusted_neuronal_beta,
                r$adjusted_neuronal_ci_lower, r$adjusted_neuronal_ci_upper,
                format_p(r$adjusted_neuronal_q),
                ifelse(r$matches_prediction, "consistent", "opposite")))
  }

  cat("  K119ub decile control (deciles consistent and significant, of 10):\n")
  for (mark in MARK_ORDER) {
    rows <- decile_tests[decile_tests$mark == mark, ]
    cat(sprintf("    %-10s %d\n", mark_display(mark),
                sum(rows$matches_prediction & rows$q_value < Q_THRESHOLD)))
  }

  cat("  Gene-length control (quintiles consistent and significant, of 5):\n")
  for (i in seq_len(nrow(quintile_survival))) {
    r <- quintile_survival[i, ]
    cat(sprintf("    %-10s %d | survives: %s\n", r$mark_display,
                r$n_consistent_and_significant,
                ifelse(r$effect_survives_length_control, "yes", "no")))
  }

  cat("  Fisher odds ratios, neuronal against the predicted mark direction:\n")
  for (i in seq_len(nrow(direction_fishers))) {
    r <- direction_fishers[i, ]
    cat(sprintf("    %-10s OR = %.3f [%.3f, %.3f], %s\n", r$mark_display,
                r$odds_ratio, r$ci_lower, r$ci_upper, format_p(r$q_value)))
  }

  cat("  Concordance count Fisher odds ratios:\n")
  for (i in seq_len(nrow(concordance$fisher))) {
    r <- concordance$fisher[i, ]
    cat(sprintf("    %-24s OR = %.3f [%.3f, %.3f], %s (%.1f%% neuronal vs %.1f%% other)\n",
                r$test_id, r$odds_ratio, r$ci_lower, r$ci_upper,
                format_p(r$p_value), 100 * r$frac_flagged_neuronal,
                100 * r$frac_flagged_other))
  }

  cat(sprintf("\n  Figures and tables: %s\n", out_dir))
  cat(sprintf("  Fisher registry:    %s\n", HANDOFF_PATHS$fisher_registry))
  cat("\n=== Section 70_02 complete ===\n\n")
}

main()
