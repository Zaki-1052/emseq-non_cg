# scripts/sections/60_mecp2/60_03_reconciliation.R
#
# Section 60_03: MeCP2 peak-level gains against gene-level MeCP2 loss.
#
# What this tests
#   DiffBind calls far more MeCP2 peaks gained than lost in the mutant, yet
#   gene-level MeCP2 aggregates fall. Three counting rules can produce that
#   pattern and this section measures all three.
#     1. Peak counting gives every peak the same weight, so a few wide or
#        strong lost peaks can carry more signal than many narrow gained peaks.
#     2. The rule that collapses many peaks into one number per gene decides
#        which peaks reach the gene level at all.
#     3. The peakset itself. This project holds two MeCP2 DiffBind tables and
#        they are not the same peaks (see "Two peaksets" below).
#
# Two peaksets
#   MECP2_PATHS$annotated and DIFFBIND_PATHS$mecp2 are separate DiffBind
#   outputs with different peak counts, different Conc ranges, and different
#   significant-peak counts. This script never merges them. Every peak-level
#   table and figure carries a peakset column and reports both. Gene-level work
#   uses the annotated peakset, which is the only one carrying SYMBOL and
#   distanceToTSS.
#
# Analyses
#   1. Peak counts and fold distribution by direction, in both peaksets.
#   2. Peak width and concentration distribution by direction, with Wilcoxon
#      tests of gained against lost in each peakset.
#   3. The gained-against-lost comparison under four weighting schemes: peak
#      count, total width, summed linear concentration, and summed signal mass
#      (width x linear concentration).
#   4. Signal mass balance: total mutant minus control signal mass per
#      direction class and over each whole peakset.
#   5. Overlap of each peakset with MeCP2_up.bed, MeCP2_down.bed, and the
#      MeCP2 consensus peakset.
#   6. Gene-level aggregation under both collapse rules of
#      aggregate_diffbind_by_gene(), nearest_tss and median_significant,
#      compared gene by gene.
#   7. Distribution of annotated MeCP2 peaks per gene.
#   8. GSEA of genes ranked by MeCP2 fold (both collapse rules) and by
#      H2AK119ub gene-body log2 fold change, compared by how many significant
#      neuronal GO BP terms each ranking recovers.
#
# Reads
#   mecp2_diffbind, mecp2_consensus            (shared config)
#   MECP2_PATHS$annotated, MECP2_PATHS$up, MECP2_PATHS$down
#   DIFFBIND_PATHS$k119ub_gene_signal
#
# Writes
#   Figures and TSV tables into OUTPUT_PATHS$mecp2 (override with --output-dir).
#   Two rows into HANDOFF_PATHS$fisher_registry.
#
# Adapted from Biomodal section 75, MeCP2 signal direction reconciliation.

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)

SECTION_ID <- "60_03"

DIRECTION_LEVELS <- c("Gained", "Lost", "Unchanged")

DIRECTION_COLORS <- c(
  "Gained"    = unname(COLORS$mecp2["MeCP2 Up"]),
  "Lost"      = unname(COLORS$mecp2["MeCP2 Down"]),
  "Unchanged" = "grey70"
)

PEAKSET_ANNOTATED <- "Annotated peakset"
PEAKSET_DIFFBIND <- "DiffBind Conc4 peakset"
PEAKSET_LEVELS <- c(PEAKSET_ANNOTATED, PEAKSET_DIFFBIND)

# Regular expression that marks a GO term description as neuronal.
NEURONAL_PATTERN <- "synap|neuron|axon|dendrit|nervous|brain|cerebel"

# Columns each input table must carry.
PEAK_REQUIRED <- c("width", "Conc", "Conc_ctrl", "Conc_mut", "Fold", "FDR")
ANNOTATED_REQUIRED <- c("seqnames", "start", "end", "width", "Conc",
                        "Conc_ctrl", "Conc_mut", "Fold", "p.value", "FDR",
                        "annotation", "geneLength", "distanceToTSS", "SYMBOL")
K119UB_SIGNAL_REQUIRED <- c("symbol", "chr", "gb_log2fc", "gb_signal_class")

# =============================================================================
# OPTIONS
# =============================================================================

parse_options <- function() {
  option_list <- list(
    make_option("--output-dir", type = "character", default = OUTPUT_PATHS$mecp2,
                dest = "output_dir",
                help = "Directory for figures and tables [default: %default]"),
    make_option("--fdr-threshold", type = "double", default = Q_THRESHOLD,
                dest = "fdr_threshold",
                help = "FDR cutoff for peak and term significance [default: %default]"),
    make_option("--gsea-min-size", type = "integer", default = 15L,
                dest = "gsea_min_size",
                help = "Smallest GO gene set tested by GSEA [default: %default]"),
    make_option("--gsea-max-size", type = "integer", default = 500L,
                dest = "gsea_max_size",
                help = "Largest GO gene set tested by GSEA [default: %default]"),
    make_option("--seed", type = "integer", default = 42L, dest = "seed",
                help = "Random seed for the GSEA permutations [default: %default]")
  )
  parse_args(OptionParser(option_list = option_list))
}

# =============================================================================
# FORMATTING HELPERS
# =============================================================================

#' Format a number with thousands separators and four significant digits.
#'
#' @param x Numeric vector.
#' @return character vector
fmt_big <- function(x) {
  formatC(x, format = "fg", digits = 4, big.mark = ",")
}

#' Format a count with thousands separators.
#'
#' @param x Numeric vector.
#' @return character vector
fmt_count <- function(x) {
  format(x, big.mark = ",", trim = TRUE)
}

# =============================================================================
# INPUT LOADING
# =============================================================================

#' Check that a table carries every required column.
#'
#' @param df data.frame to check.
#' @param required Character vector of column names.
#' @param source_label Name used in the error message.
#' @return the data.frame, unchanged
require_columns <- function(df, required, source_label) {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(source_label, " is missing columns: ", paste(missing, collapse = ", "))
  }
  df
}

#' Classify one peak table at a single FDR cutoff and add the weight columns.
#'
#' Recomputes direction so the --fdr-threshold option controls every peak count
#' in this section. linear_conc is 2^Conc. mass_ctrl and mass_mut multiply peak
#' width by the linear-scale group concentration, and mass_delta is their
#' difference.
#'
#' @param df Peak table carrying PEAK_REQUIRED and either Chr or seqnames.
#' @param fdr_threshold FDR cutoff separating gained and lost from unchanged.
#' @param peakset_label Name recorded in the peakset column.
#' @return data.frame with Chr, Start, End, direction, and the weight columns
prepare_peaks <- function(df, fdr_threshold, peakset_label) {
  df <- require_columns(df, PEAK_REQUIRED, peakset_label)

  if ("Chr" %in% colnames(df)) {
    chr <- as.character(df$Chr)
    start <- df$Start
    end <- df$End
  } else if ("seqnames" %in% colnames(df)) {
    chr <- as.character(df$seqnames)
    start <- df$start
    end <- df$end
  } else {
    stop(peakset_label, " has neither a Chr nor a seqnames column.")
  }

  needs_prefix <- !grepl("^chr", chr)
  chr[needs_prefix] <- paste0("chr", chr[needs_prefix])

  peaks <- data.frame(
    peakset = factor(peakset_label, levels = PEAKSET_LEVELS),
    Chr = chr,
    Start = as.numeric(start),
    End = as.numeric(end),
    width = as.numeric(df$width),
    Conc = as.numeric(df$Conc),
    Conc_ctrl = as.numeric(df$Conc_ctrl),
    Conc_mut = as.numeric(df$Conc_mut),
    Fold = as.numeric(df$Fold),
    FDR = as.numeric(df$FDR),
    stringsAsFactors = FALSE
  )

  for (col in c("Start", "End", "width", "Conc", "Conc_ctrl", "Conc_mut",
                "Fold", "FDR")) {
    n_bad <- sum(!is.finite(peaks[[col]]))
    if (n_bad > 0) {
      stop(peakset_label, " has ", n_bad, " non-finite values in column ", col)
    }
  }

  # The width violin and the width weighting both need a positive width.
  n_zero_width <- sum(peaks$width <= 0)
  if (n_zero_width > 0) {
    stop(peakset_label, " has ", n_zero_width, " peaks of width zero or less.")
  }

  peaks$direction <- "Unchanged"
  peaks$direction[peaks$FDR < fdr_threshold & peaks$Fold > 0] <- "Gained"
  peaks$direction[peaks$FDR < fdr_threshold & peaks$Fold < 0] <- "Lost"
  peaks$direction <- factor(peaks$direction, levels = DIRECTION_LEVELS)

  peaks$linear_conc <- 2^peaks$Conc
  peaks$mass_ctrl <- peaks$width * 2^peaks$Conc_ctrl
  peaks$mass_mut <- peaks$width * 2^peaks$Conc_mut
  peaks$mass_delta <- peaks$mass_mut - peaks$mass_ctrl
  peaks$signal_mass <- peaks$width * peaks$linear_conc

  cat(sprintf("  %-24s %s peaks: %s gained, %s lost, %s unchanged at FDR<%.3f\n",
              peakset_label, fmt_count(nrow(peaks)),
              fmt_count(sum(peaks$direction == "Gained")),
              fmt_count(sum(peaks$direction == "Lost")),
              fmt_count(sum(peaks$direction == "Unchanged")),
              fdr_threshold))
  cat(sprintf("  %-24s Conc range %.3f to %.3f, width median %.0f bp, %s peaks off the canonical chromosomes\n",
              peakset_label, min(peaks$Conc), max(peaks$Conc),
              median(peaks$width),
              fmt_count(sum(!peaks$Chr %in% CANONICAL_CHRS))))
  peaks
}

#' Read the ChIPseeker-annotated MeCP2 peak table.
#'
#' @param filepath Path to MeCP2_annotated.txt.
#' @return data.frame of every annotated peak
load_mecp2_annotated <- function(filepath) {
  df <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "", fill = TRUE)
  df <- require_columns(df, ANNOTATED_REQUIRED, basename(filepath))

  numeric_cols <- c("start", "end", "width", "Conc", "Conc_ctrl", "Conc_mut",
                    "Fold", "p.value", "FDR", "geneLength", "distanceToTSS")
  for (col in numeric_cols) df[[col]] <- as.numeric(df[[col]])

  for (col in c("Fold", "FDR", "distanceToTSS", "width")) {
    n_bad <- sum(!is.finite(df[[col]]))
    if (n_bad > 0) {
      stop(basename(filepath), " has ", n_bad,
           " non-finite values in column ", col)
    }
  }

  cat(sprintf("  Annotated peak file: %s peaks, %s carry a gene symbol\n",
              fmt_count(nrow(df)),
              fmt_count(sum(!is.na(df$SYMBOL) & nzchar(df$SYMBOL)))))
  df
}

#' Keep the annotated peaks that carry a gene symbol.
#'
#' @param annotated data.frame from load_mecp2_annotated().
#' @return data.frame of peaks with a non-empty SYMBOL
filter_to_symbol <- function(annotated) {
  annotated[!is.na(annotated$SYMBOL) & nzchar(annotated$SYMBOL), , drop = FALSE]
}

#' Read the H2AK119ub gene-body signal table and keep quantifiable genes.
#'
#' @param filepath Path to k119ub_gene_signal.tsv.
#' @return data.frame with symbol, chr, and gb_log2fc
load_k119ub_gene_signal <- function(filepath) {
  df <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")
  df <- require_columns(df, K119UB_SIGNAL_REQUIRED, basename(filepath))

  df$gb_log2fc <- as.numeric(df$gb_log2fc)
  keep <- df$gb_signal_class == "quantifiable" & is.finite(df$gb_log2fc) &
    !is.na(df$symbol) & nzchar(df$symbol)
  out <- df[keep, , drop = FALSE]

  cat(sprintf("  K119ub gene signal: %s rows, %s quantifiable genes\n",
              fmt_count(nrow(df)), fmt_count(nrow(out))))
  out
}

# =============================================================================
# ANALYSIS 1: PEAK COUNTS BY DIRECTION
# =============================================================================

#' Count peaks per direction in each peakset and draw the paired bar.
#'
#' @param peaks Combined peak table carrying the peakset column.
#' @param out_dir Section output directory.
#' @return list with the counts table and the ggplot
plot_direction_counts <- function(peaks, out_dir) {
  counts <- peaks %>%
    dplyr::count(peakset, direction, name = "n_peaks", .drop = FALSE) %>%
    dplyr::group_by(peakset) %>%
    dplyr::mutate(percentage = 100 * n_peaks / sum(n_peaks)) %>%
    dplyr::ungroup() %>%
    as.data.frame()

  write_section_table(counts,
                      file.path(out_dir, "mecp2_peak_direction_counts.tsv"))

  for (i in seq_len(nrow(counts))) {
    cat(sprintf("  %-24s %-10s %s peaks (%.1f%%)\n",
                as.character(counts$peakset[i]),
                as.character(counts$direction[i]),
                fmt_count(counts$n_peaks[i]), counts$percentage[i]))
  }

  ratio_tbl <- counts %>%
    dplyr::filter(direction %in% c("Gained", "Lost")) %>%
    dplyr::select(peakset, direction, n_peaks) %>%
    tidyr::pivot_wider(names_from = direction, values_from = n_peaks) %>%
    dplyr::mutate(gained_over_lost = Gained / Lost) %>%
    as.data.frame()

  for (i in seq_len(nrow(ratio_tbl))) {
    cat(sprintf("  %-24s gained-to-lost ratio %.2f\n",
                as.character(ratio_tbl$peakset[i]),
                ratio_tbl$gained_over_lost[i]))
  }

  paired <- counts[counts$direction %in% c("Gained", "Lost"), , drop = FALSE]
  paired$direction <- droplevels(paired$direction)

  p <- ggplot(paired, aes(x = direction, y = n_peaks, fill = direction)) +
    geom_col(width = 0.6, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("n = %s\n(%.1f%% of peakset)",
                                  fmt_count(n_peaks), percentage)),
              vjust = -0.25, size = 3.6, lineheight = 1.1) +
    facet_wrap(~ peakset, scales = "free_y") +
    scale_fill_manual(values = DIRECTION_COLORS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.30)),
                       labels = scales::comma) +
    labs(
      title = "MeCP2 Differential Peaks by Direction",
      subtitle = paste0(
        "Gained-to-lost ratio: ",
        paste(sprintf("%s = %.2f", as.character(ratio_tbl$peakset),
                      ratio_tbl$gained_over_lost), collapse = " | ")),
      x = "MeCP2 peak direction (mutant vs control)",
      y = "Number of peaks"
    ) +
    theme_emseq()

  save_multiformat_ggplot(p, file.path(out_dir, "60_03a_peak_direction_counts"),
                          width = 11, height = 7)
  list(counts = counts, ratios = ratio_tbl, plot = p)
}

# =============================================================================
# ANALYSIS 2: PEAK PROPERTY DISTRIBUTIONS
# =============================================================================

#' Per-direction summary of one property, computed inside each peakset, with the
#' figure text for every row.
#'
#' The result goes to a geom_text layer only. No table written by this section
#' comes from it.
#'
#' @param peaks Combined peak table.
#' @param value_col Numeric column to summarise.
#' @param digits Digits for the median shown on the figure.
#' @return data.frame with direction, n, median, mean, q25, q75, peakset, label
label_summary_by_peakset <- function(peaks, value_col, digits) {
  do.call(rbind, lapply(PEAKSET_LEVELS, function(ps) {
    sub <- peaks[peaks$peakset == ps, , drop = FALSE]
    out <- summarise_groups(sub, "direction", value_col)
    out$peakset <- factor(ps, levels = PEAKSET_LEVELS)
    out$label <- group_label(out, digits = digits)
    out
  }))
}

#' Violin of one peak property split by direction and faceted by peakset.
#'
#' Every group carries its peak count and median on the figure.
#'
#' @param peaks Combined peak table.
#' @param value_col Numeric column to plot.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param y_lab Y axis label.
#' @param digits Digits used for the median shown on the figure.
#' @param log_scale TRUE to draw the y axis on a log10 scale.
#' @return list with the ggplot and the per-group figure text
plot_property_violin <- function(peaks, value_col, title, subtitle, y_lab,
                                 digits = 3, log_scale = FALSE) {
  grp <- label_summary_by_peakset(peaks, value_col, digits)
  vals <- peaks[[value_col]]

  if (log_scale) {
    y_label_pos <- max(vals) * 2.5
    y_limits <- c(min(vals), max(vals) * 7)
  } else {
    span <- max(vals) - min(vals)
    y_label_pos <- max(vals) + 0.12 * span
    y_limits <- c(min(vals) - 0.05 * span, max(vals) + 0.28 * span)
  }

  p <- ggplot(peaks, aes(x = direction, y = .data[[value_col]],
                         fill = direction)) +
    geom_violin(alpha = 0.6, scale = "width",
                draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white", alpha = 0.85) +
    geom_text(data = grp, aes(x = direction, y = y_label_pos, label = label),
              inherit.aes = FALSE, size = 3, lineheight = 1.1) +
    facet_wrap(~ peakset) +
    scale_fill_manual(values = DIRECTION_COLORS, guide = "none") +
    labs(title = title, subtitle = subtitle,
         x = "MeCP2 peak direction (mutant vs control)", y = y_lab) +
    theme_emseq()

  if (log_scale) {
    p <- p + scale_y_log10(limits = y_limits, labels = scales::comma)
  } else {
    p <- p + scale_y_continuous(limits = y_limits)
  }

  list(plot = p, labels = grp)
}

#' Full distribution summary of one property per peakset and direction.
#'
#' @param peaks Combined peak table.
#' @param value_col Numeric column to summarise.
#' @return data.frame with one row per peakset and direction
summarise_property <- function(peaks, value_col) {
  peaks %>%
    dplyr::group_by(peakset, direction) %>%
    dplyr::summarise(
      property = value_col,
      n = dplyr::n(),
      total = sum(.data[[value_col]]),
      mean = mean(.data[[value_col]]),
      sd = sd(.data[[value_col]]),
      min = min(.data[[value_col]]),
      q25 = unname(quantile(.data[[value_col]], 0.25)),
      median = median(.data[[value_col]]),
      q75 = unname(quantile(.data[[value_col]], 0.75)),
      max = max(.data[[value_col]]),
      .groups = "drop"
    ) %>%
    dplyr::select(property, peakset, direction, dplyr::everything()) %>%
    as.data.frame()
}

#' Wilcoxon rank-sum test of one property, gained against lost, per peakset.
#'
#' @param peaks Combined peak table.
#' @param value_col Numeric column to test.
#' @return data.frame with one row per peakset
wilcoxon_gained_vs_lost <- function(peaks, value_col) {
  do.call(rbind, lapply(PEAKSET_LEVELS, function(ps) {
    sub <- peaks[peaks$peakset == ps, , drop = FALSE]
    a <- sub[[value_col]][sub$direction == "Gained"]
    b <- sub[[value_col]][sub$direction == "Lost"]
    wt <- wilcox.test(a, b)
    data.frame(
      peakset = ps, property = value_col,
      n_gained = length(a), n_lost = length(b),
      median_gained = median(a), median_lost = median(b),
      median_difference = median(a) - median(b),
      mean_gained = mean(a), mean_lost = mean(b),
      W = unname(wt$statistic), p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  }))
}

#' Draw the three property violins and write their summary and test tables.
#'
#' @param peaks Combined peak table.
#' @param out_dir Section output directory.
#' @return list with the three ggplots and the Wilcoxon table
run_property_analysis <- function(peaks, out_dir) {
  fold_res <- plot_property_violin(
    peaks, "Fold",
    title = "MeCP2 Peak Fold Change by Direction",
    subtitle = "DiffBind log2 fold change, mutant over control, one point per peak",
    y_lab = "MeCP2 log2 fold change (mutant / control)",
    digits = 3)

  width_res <- plot_property_violin(
    peaks, "width",
    title = "MeCP2 Peak Width by Direction",
    subtitle = "Peak width sets how much sequence each peak contributes to a gene-level aggregate",
    y_lab = "Peak width (bp, log10 scale)",
    digits = 0, log_scale = TRUE)

  conc_res <- plot_property_violin(
    peaks, "Conc",
    title = "MeCP2 Peak Concentration by Direction",
    subtitle = "DiffBind mean normalised read concentration across the libraries (log2)",
    y_lab = "Peak concentration (log2 normalised reads)",
    digits = 3)

  save_multiformat_ggplot(fold_res$plot,
                          file.path(out_dir, "60_03b_peak_fold_by_direction"),
                          width = 11, height = 7)
  save_multiformat_ggplot(width_res$plot,
                          file.path(out_dir, "60_03c_peak_width_by_direction"),
                          width = 11, height = 7)
  save_multiformat_ggplot(conc_res$plot,
                          file.path(out_dir, "60_03d_peak_concentration_by_direction"),
                          width = 11, height = 7)

  summary_tbl <- rbind(
    summarise_property(peaks, "Fold"),
    summarise_property(peaks, "width"),
    summarise_property(peaks, "Conc"),
    summarise_property(peaks, "linear_conc"),
    summarise_property(peaks, "signal_mass"),
    summarise_property(peaks, "mass_delta")
  )
  write_section_table(summary_tbl,
                      file.path(out_dir, "mecp2_peak_property_summary.tsv"))

  wilcox_tbl <- rbind(
    wilcoxon_gained_vs_lost(peaks, "Fold"),
    wilcoxon_gained_vs_lost(peaks, "width"),
    wilcoxon_gained_vs_lost(peaks, "Conc"),
    wilcoxon_gained_vs_lost(peaks, "signal_mass")
  )
  write_section_table(wilcox_tbl,
                      file.path(out_dir, "mecp2_peak_property_wilcoxon.tsv"))

  for (i in seq_len(nrow(wilcox_tbl))) {
    cat(sprintf("  %-24s Wilcoxon %-12s gained median = %s, lost median = %s, W = %.0f, p = %.3g\n",
                wilcox_tbl$peakset[i], wilcox_tbl$property[i],
                fmt_big(wilcox_tbl$median_gained[i]),
                fmt_big(wilcox_tbl$median_lost[i]),
                wilcox_tbl$W[i], wilcox_tbl$p_value[i]))
  }

  list(fold_plot = fold_res$plot, width_plot = width_res$plot,
       conc_plot = conc_res$plot, wilcoxon = wilcox_tbl)
}

#' Peak-level Fisher test of wide peaks against the gained direction.
#'
#' Wide means wider than the median width of the peakset the peak belongs to.
#' The test runs over peaks that pass the FDR cutoff in either direction, once
#' per peakset.
#'
#' @param peaks Combined peak table.
#' @param out_dir Section output directory.
#' @return data.frame with one row per peakset
run_peak_width_fisher <- function(peaks, out_dir) {
  out <- do.call(rbind, lapply(PEAKSET_LEVELS, function(ps) {
    sub <- peaks[peaks$peakset == ps, , drop = FALSE]
    median_width <- median(sub$width)
    sig <- sub[sub$direction %in% c("Gained", "Lost"), , drop = FALSE]
    tab <- table(factor(sig$width > median_width, levels = c(TRUE, FALSE)),
                 factor(sig$direction == "Gained", levels = c(TRUE, FALSE)))
    ft <- fisher.test(tab)
    data.frame(
      peakset = ps,
      median_width_all_peaks = median_width,
      n_significant_peaks = nrow(sig),
      n_wide_gained = tab["TRUE", "TRUE"],
      n_wide_lost = tab["TRUE", "FALSE"],
      n_narrow_gained = tab["FALSE", "TRUE"],
      n_narrow_lost = tab["FALSE", "FALSE"],
      odds_ratio = unname(ft$estimate),
      ci_low = ft$conf.int[1], ci_high = ft$conf.int[2],
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  }))
  write_section_table(out,
                      file.path(out_dir, "mecp2_peak_width_fisher.tsv"))

  for (i in seq_len(nrow(out))) {
    cat(sprintf("  %-24s Fisher (wide peak x gained): OR = %.3f (95%% CI %.3f to %.3f), p = %.3g\n",
                out$peakset[i], out$odds_ratio[i], out$ci_low[i],
                out$ci_high[i], out$p_value[i]))
  }
  out
}

# =============================================================================
# ANALYSIS 3: WEIGHTED GAINED-AGAINST-LOST COMPARISON
# =============================================================================

#' Compare gained and lost peaks under four weighting schemes.
#'
#' peak_count gives every peak the same weight. total_width_bp weights each
#' peak by its width. sum_linear_concentration weights each peak by 2^Conc.
#' sum_signal_mass multiplies width by 2^Conc.
#'
#' @param peaks Combined peak table.
#' @param out_dir Section output directory.
#' @return list with the long table, the wide table, and the ggplot
run_weighted_comparison <- function(peaks, out_dir) {
  sig <- peaks[peaks$direction %in% c("Gained", "Lost"), , drop = FALSE]
  sig$direction <- droplevels(sig$direction)

  scheme_levels <- c("peak_count", "total_width_bp",
                     "sum_linear_concentration", "sum_signal_mass")
  scheme_labels <- c(
    peak_count = "Unweighted peak count",
    total_width_bp = "Weighted by peak width (total bp)",
    sum_linear_concentration = "Weighted by concentration (sum of 2^Conc)",
    sum_signal_mass = "Weighted by signal mass (width x 2^Conc)"
  )

  long <- sig %>%
    dplyr::group_by(peakset, direction) %>%
    dplyr::summarise(
      peak_count = as.numeric(dplyr::n()),
      total_width_bp = sum(width),
      sum_linear_concentration = sum(linear_conc),
      sum_signal_mass = sum(signal_mass),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(scheme_levels),
                        names_to = "scheme", values_to = "value") %>%
    as.data.frame()

  long$scheme <- factor(long$scheme, levels = scheme_levels)
  long$scheme_label <- factor(unname(scheme_labels[as.character(long$scheme)]),
                              levels = unname(scheme_labels[scheme_levels]))

  wide <- long %>%
    dplyr::select(peakset, scheme, direction, value) %>%
    tidyr::pivot_wider(names_from = direction, values_from = value) %>%
    dplyr::mutate(
      gained_over_lost = Gained / Lost,
      gained_share_pct = 100 * Gained / (Gained + Lost)
    ) %>%
    as.data.frame()

  write_section_table(
    long, file.path(out_dir, "mecp2_weighted_gained_vs_lost_long.tsv"))
  write_section_table(
    wide, file.path(out_dir, "mecp2_weighted_gained_vs_lost.tsv"))

  for (i in seq_len(nrow(wide))) {
    cat(sprintf("  %-24s %-26s gained = %s, lost = %s, ratio = %.2f, gained share = %.1f%%\n",
                as.character(wide$peakset[i]), as.character(wide$scheme[i]),
                fmt_big(wide$Gained[i]), fmt_big(wide$Lost[i]),
                wide$gained_over_lost[i], wide$gained_share_pct[i]))
  }

  p <- ggplot(long, aes(x = direction, y = value, fill = direction)) +
    geom_col(width = 0.6, color = "black", linewidth = 0.3) +
    geom_text(aes(label = fmt_big(value)), vjust = -0.35, size = 3.1) +
    facet_wrap(~ peakset + scheme_label, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = DIRECTION_COLORS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25)),
                       labels = scales::comma) +
    labs(
      title = "Gained Against Lost MeCP2 Peaks Under Four Weighting Schemes",
      subtitle = "Each panel reweights the same peaks. The gained-to-lost balance moves with the weight.",
      x = "MeCP2 peak direction (mutant vs control)",
      y = "Total weight"
    ) +
    theme_emseq() +
    theme(strip.text = element_text(size = 8))

  save_multiformat_ggplot(p, file.path(out_dir, "60_03e_weighted_gained_vs_lost"),
                          width = 18, height = 10)
  list(long = long, wide = wide, plot = p)
}

#' Total mutant-minus-control signal mass per direction class and per peakset.
#'
#' Signal mass is peak width multiplied by the linear-scale group
#' concentration. Each bar sums mass_mut minus mass_ctrl over the peaks of one
#' class, and the Net bar sums it over the whole peakset.
#'
#' @param peaks Combined peak table.
#' @param out_dir Section output directory.
#' @return list with the balance table and the ggplot
run_signal_mass_balance <- function(peaks, out_dir) {
  by_class <- peaks %>%
    dplyr::group_by(peakset, direction) %>%
    dplyr::summarise(
      n_peaks = dplyr::n(),
      total_mass_ctrl = sum(mass_ctrl),
      total_mass_mut = sum(mass_mut),
      total_mass_delta = sum(mass_delta),
      .groups = "drop"
    ) %>%
    dplyr::mutate(class = as.character(direction)) %>%
    dplyr::select(-direction) %>%
    as.data.frame()

  net <- peaks %>%
    dplyr::group_by(peakset) %>%
    dplyr::summarise(
      n_peaks = dplyr::n(),
      total_mass_ctrl = sum(mass_ctrl),
      total_mass_mut = sum(mass_mut),
      total_mass_delta = sum(mass_delta),
      .groups = "drop"
    ) %>%
    dplyr::mutate(class = "Net (all peaks)") %>%
    as.data.frame()

  balance <- rbind(by_class, net)
  balance$pct_change <- 100 * balance$total_mass_delta / balance$total_mass_ctrl
  balance$class <- factor(balance$class,
                          levels = c(DIRECTION_LEVELS, "Net (all peaks)"))
  balance <- balance[order(balance$peakset, balance$class), ]

  write_section_table(balance,
                      file.path(out_dir, "mecp2_signal_mass_balance.tsv"))

  for (i in seq_len(nrow(balance))) {
    cat(sprintf("  %-24s %-16s n = %s, mass change = %s (%.2f%% of control mass)\n",
                as.character(balance$peakset[i]),
                as.character(balance$class[i]), fmt_count(balance$n_peaks[i]),
                fmt_big(balance$total_mass_delta[i]), balance$pct_change[i]))
  }

  balance$delta_millions <- balance$total_mass_delta / 1e6
  bar_colors <- c(DIRECTION_COLORS, "Net (all peaks)" = "grey25")

  p <- ggplot(balance, aes(x = class, y = delta_millions, fill = class)) +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
    geom_col(width = 0.62, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%s\n(%.2f%%)", fmt_big(delta_millions),
                                  pct_change),
                  vjust = ifelse(delta_millions >= 0, -0.35, 1.25)),
              size = 3, lineheight = 1.1) +
    facet_wrap(~ peakset, scales = "free_y") +
    scale_fill_manual(values = bar_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.28, 0.28))) +
    labs(
      title = "MeCP2 Signal Mass Change, Mutant Minus Control",
      subtitle = paste("Signal mass is peak width multiplied by the linear-scale group concentration.",
                       "Percentages are of the control mass of the same class."),
      x = "Peak class",
      y = "Signal mass change (millions of bp x normalised reads)"
    ) +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  save_multiformat_ggplot(p, file.path(out_dir, "60_03f_signal_mass_balance"),
                          width = 12, height = 7)
  list(balance = balance, plot = p)
}

# =============================================================================
# ANALYSIS 4: PEAKSET OVERLAP
# =============================================================================

#' Overlap of each peakset with the MeCP2 up, down, and consensus interval sets.
#'
#' For each peakset the table reports how many gained peaks overlap
#' MeCP2_up.bed, how many lost peaks overlap MeCP2_down.bed, and how many peaks
#' of each direction fall inside the MeCP2 consensus peakset.
#'
#' @param peaks Combined peak table.
#' @param up_gr GRanges from MeCP2_up.bed.
#' @param down_gr GRanges from MeCP2_down.bed.
#' @param consensus_gr mecp2_consensus from the shared config.
#' @param out_dir Section output directory.
#' @return data.frame with one row per peakset and direction
run_interval_overlap <- function(peaks, up_gr, down_gr, consensus_gr, out_dir) {
  peaks_gr <- GRanges(seqnames = peaks$Chr,
                      ranges = IRanges(start = peaks$Start, end = peaks$End))

  peaks$overlaps_up_bed <- countOverlaps(peaks_gr, up_gr) > 0
  peaks$overlaps_down_bed <- countOverlaps(peaks_gr, down_gr) > 0
  peaks$in_consensus <- countOverlaps(peaks_gr, consensus_gr) > 0

  out <- peaks %>%
    dplyr::group_by(peakset, direction) %>%
    dplyr::summarise(
      n_peaks = dplyr::n(),
      n_overlap_up_bed = sum(overlaps_up_bed),
      n_overlap_down_bed = sum(overlaps_down_bed),
      n_in_consensus = sum(in_consensus),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      pct_overlap_up_bed = 100 * n_overlap_up_bed / n_peaks,
      pct_overlap_down_bed = 100 * n_overlap_down_bed / n_peaks,
      pct_in_consensus = 100 * n_in_consensus / n_peaks
    ) %>%
    as.data.frame()

  bed_sizes <- data.frame(
    interval_set = c("MeCP2_up.bed", "MeCP2_down.bed", "MeCP2 consensus"),
    n_intervals = c(length(up_gr), length(down_gr), length(consensus_gr)),
    stringsAsFactors = FALSE
  )

  write_section_table(out, file.path(out_dir, "mecp2_interval_overlap.tsv"))
  write_section_table(bed_sizes,
                      file.path(out_dir, "mecp2_interval_set_sizes.tsv"))

  for (i in seq_len(nrow(bed_sizes))) {
    cat(sprintf("  %-16s %s intervals\n", bed_sizes$interval_set[i],
                fmt_count(bed_sizes$n_intervals[i])))
  }
  for (i in seq_len(nrow(out))) {
    cat(sprintf("  %-24s %-10s up.bed %.1f%%, down.bed %.1f%%, consensus %.1f%% of %s peaks\n",
                as.character(out$peakset[i]), as.character(out$direction[i]),
                out$pct_overlap_up_bed[i], out$pct_overlap_down_bed[i],
                out$pct_in_consensus[i], fmt_count(out$n_peaks[i])))
  }
  out
}

# =============================================================================
# ANALYSIS 5: GENE-LEVEL COLLAPSE RULES
# =============================================================================

#' Aggregate the annotated peaks to genes under both collapse rules.
#'
#' nearest_tss keeps the peak closest to the TSS and calls the gene from that
#' peak's own FDR. median_significant takes the median fold of the peaks that
#' pass the FDR cutoff and calls the gene from whether any peak passes. Both
#' rules come from aggregate_diffbind_by_gene() in the shared config.
#'
#' @param annotated_sym Annotated peaks carrying a gene symbol.
#' @param fdr_threshold FDR cutoff passed to both rules.
#' @return data.frame with one row per gene
build_collapse_table <- function(annotated_sym, fdr_threshold) {
  tss_tbl <- aggregate_diffbind_by_gene(annotated_sym, method = "nearest_tss",
                                        fdr_threshold = fdr_threshold,
                                        prefix = "mecp2_tss")
  med_tbl <- aggregate_diffbind_by_gene(annotated_sym,
                                        method = "median_significant",
                                        fdr_threshold = fdr_threshold,
                                        prefix = "mecp2_med")

  if (!setequal(tss_tbl$gene_name, med_tbl$gene_name)) {
    stop("The two collapse rules returned different gene sets: ",
         nrow(tss_tbl), " and ", nrow(med_tbl), " genes.")
  }

  gene_extras <- annotated_sym %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(
      chr = as.character(seqnames[which.min(abs(distanceToTSS))]),
      min_abs_dist_tss = min(abs(distanceToTSS)),
      total_peak_width = sum(width),
      annot_gene_length = geneLength[which.min(abs(distanceToTSS))],
      .groups = "drop"
    ) %>%
    dplyr::rename(gene_name = SYMBOL) %>%
    as.data.frame()

  tbl <- tss_tbl %>%
    dplyr::inner_join(med_tbl, by = "gene_name") %>%
    dplyr::inner_join(gene_extras, by = "gene_name") %>%
    as.data.frame()

  needs_prefix <- !grepl("^chr", tbl$chr)
  tbl$chr[needs_prefix] <- paste0("chr", tbl$chr[needs_prefix])

  tbl$gained_nearest_tss <- tbl$mecp2_tss_fdr < fdr_threshold &
    tbl$mecp2_tss_fold > 0
  tbl$lost_nearest_tss <- tbl$mecp2_tss_fdr < fdr_threshold &
    tbl$mecp2_tss_fold < 0
  tbl$gained_median_sig <- tbl$mecp2_med_has_sig & tbl$mecp2_med_fold > 0
  tbl$lost_median_sig <- tbl$mecp2_med_has_sig & tbl$mecp2_med_fold < 0
  tbl$sign_flip <- sign(tbl$mecp2_tss_fold) != sign(tbl$mecp2_med_fold)
  tbl$many_peaks <- tbl$mecp2_tss_n_peaks > median(tbl$mecp2_tss_n_peaks)

  cat(sprintf("  Genes with an annotated MeCP2 peak: %s\n",
              fmt_count(nrow(tbl))))
  cat(sprintf("  Peaks per gene: median %.0f, mean %.2f, max %.0f\n",
              median(tbl$mecp2_tss_n_peaks), mean(tbl$mecp2_tss_n_peaks),
              max(tbl$mecp2_tss_n_peaks)))
  tbl
}

#' Compare the two collapse rules gene by gene and draw the scatter.
#'
#' @param collapse_tbl data.frame from build_collapse_table().
#' @param out_dir Section output directory.
#' @return list with the summary table and the ggplot
run_collapse_scatter <- function(collapse_tbl, out_dir) {
  rho_test <- suppressWarnings(cor.test(collapse_tbl$mecp2_tss_fold,
                                        collapse_tbl$mecp2_med_fold,
                                        method = "spearman"))
  pearson_test <- cor.test(collapse_tbl$mecp2_tss_fold,
                           collapse_tbl$mecp2_med_fold, method = "pearson")
  pct_flip <- 100 * mean(collapse_tbl$sign_flip)

  summary_tbl <- data.frame(
    n_genes = nrow(collapse_tbl),
    spearman_rho = unname(rho_test$estimate), spearman_p = rho_test$p.value,
    pearson_r = unname(pearson_test$estimate), pearson_p = pearson_test$p.value,
    pct_sign_flip = pct_flip,
    n_sign_flip = sum(collapse_tbl$sign_flip),
    median_fold_nearest_tss = median(collapse_tbl$mecp2_tss_fold),
    median_fold_median_sig = median(collapse_tbl$mecp2_med_fold),
    mean_fold_nearest_tss = mean(collapse_tbl$mecp2_tss_fold),
    mean_fold_median_sig = mean(collapse_tbl$mecp2_med_fold),
    stringsAsFactors = FALSE
  )
  write_section_table(summary_tbl,
                      file.path(out_dir, "mecp2_collapse_rule_summary.tsv"))

  cat(sprintf("  Collapse rules on %s genes: Spearman rho = %.3f (p = %.3g), sign differs for %.1f%%\n",
              fmt_count(nrow(collapse_tbl)), summary_tbl$spearman_rho,
              summary_tbl$spearman_p, pct_flip))
  cat(sprintf("  Median gene-level fold: nearest_tss = %.4f, median_significant = %.4f\n",
              summary_tbl$median_fold_nearest_tss,
              summary_tbl$median_fold_median_sig))

  plot_df <- collapse_tbl
  plot_df$flip_label <- ifelse(plot_df$sign_flip,
                               "Sign differs between rules",
                               "Same sign under both rules")
  plot_df$label_gene <- ifelse(plot_df$gene_name %in% KEY_GENES,
                               plot_df$gene_name, "")

  p <- ggplot(plot_df, aes(x = mecp2_tss_fold, y = mecp2_med_fold)) +
    geom_point(aes(color = flip_label), alpha = 0.3, size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "black", linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
    geom_text_repel(aes(label = label_gene), size = 3, max.overlaps = 20,
                    fontface = "italic", color = "grey15",
                    segment.color = "grey60", segment.size = 0.3,
                    min.segment.length = 0) +
    scale_color_manual(values = c("Sign differs between rules" = "#D95F02",
                                  "Same sign under both rules" = "grey45"),
                       name = NULL) +
    labs(
      title = "Gene-Level MeCP2 Fold Change Under Two Collapse Rules",
      subtitle = sprintf(
        "Spearman rho = %.3f, p = %.2e | sign differs for %s of %s genes (%.1f%%)",
        summary_tbl$spearman_rho, summary_tbl$spearman_p,
        fmt_count(summary_tbl$n_sign_flip), fmt_count(nrow(plot_df)), pct_flip),
      x = "MeCP2 log2 fold change (nearest TSS peak)",
      y = "MeCP2 log2 fold change (median of significant peaks)"
    ) +
    theme_emseq() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p, file.path(out_dir, "60_03g_collapse_rule_scatter"),
                          width = 9, height = 8)
  list(summary = summary_tbl, plot = p)
}

#' Gene calls under each collapse rule, side by side.
#'
#' @param collapse_tbl data.frame from build_collapse_table().
#' @param out_dir Section output directory.
#' @return list with the counts table and the ggplot
run_collapse_gene_calls <- function(collapse_tbl, out_dir) {
  call_of <- function(gained, lost) {
    ifelse(gained, "Gained", ifelse(lost, "Lost", "Unchanged"))
  }

  calls <- rbind(
    data.frame(rule = "nearest_tss",
               call = call_of(collapse_tbl$gained_nearest_tss,
                              collapse_tbl$lost_nearest_tss),
               stringsAsFactors = FALSE),
    data.frame(rule = "median_significant",
               call = call_of(collapse_tbl$gained_median_sig,
                              collapse_tbl$lost_median_sig),
               stringsAsFactors = FALSE)
  )
  calls$rule <- factor(calls$rule,
                       levels = c("nearest_tss", "median_significant"))
  calls$call <- factor(calls$call, levels = DIRECTION_LEVELS)

  counts <- calls %>%
    dplyr::count(rule, call, name = "n_genes", .drop = FALSE) %>%
    dplyr::group_by(rule) %>%
    dplyr::mutate(percentage = 100 * n_genes / sum(n_genes)) %>%
    dplyr::ungroup() %>%
    as.data.frame()

  write_section_table(counts,
                      file.path(out_dir, "mecp2_collapse_rule_gene_calls.tsv"))

  for (i in seq_len(nrow(counts))) {
    cat(sprintf("  %-19s %-10s %s genes (%.1f%%)\n",
                as.character(counts$rule[i]), as.character(counts$call[i]),
                fmt_count(counts$n_genes[i]), counts$percentage[i]))
  }

  p <- ggplot(counts, aes(x = call, y = n_genes, fill = call)) +
    geom_col(width = 0.65, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%s\n(%.1f%%)", fmt_count(n_genes),
                                  percentage)),
              vjust = -0.25, size = 3.2, lineheight = 1.1) +
    facet_wrap(~ rule) +
    scale_fill_manual(values = DIRECTION_COLORS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25)),
                       labels = scales::comma) +
    labs(
      title = "Gene-Level MeCP2 Calls Under Each Collapse Rule",
      subtitle = sprintf("Same %s genes, same peaks, two aggregation rules",
                         fmt_count(nrow(collapse_tbl))),
      x = "Gene-level MeCP2 call", y = "Number of genes"
    ) +
    theme_emseq()

  save_multiformat_ggplot(p, file.path(out_dir, "60_03h_collapse_rule_gene_calls"),
                          width = 10, height = 7)
  list(counts = counts, plot = p)
}

#' Register the two gene-level Fisher tests of this section.
#'
#' @param collapse_tbl data.frame from build_collapse_table().
#' @param out_dir Section output directory.
#' @return data.frame with one row per test
run_gene_fisher_tests <- function(collapse_tbl, out_dir) {
  ft_rules <- register_fisher_test(
    section = SECTION_ID, test_id = "collapse_rule_gain_agreement",
    description = paste("Do the nearest-TSS and median-significant collapse",
                        "rules call the same genes MeCP2 gained?"),
    gene_df = collapse_tbl, row_var = "gained_nearest_tss",
    col_var = "gained_median_sig", output_dir = out_dir)

  ft_peaks <- register_fisher_test(
    section = SECTION_ID, test_id = "many_peaks_vs_gain",
    description = paste("Are genes carrying more MeCP2 peaks than the median",
                        "gene more often called MeCP2 gained by the",
                        "nearest-TSS rule?"),
    gene_df = collapse_tbl, row_var = "many_peaks",
    col_var = "gained_nearest_tss", output_dir = out_dir)

  out <- data.frame(
    test_id = c("collapse_rule_gain_agreement", "many_peaks_vs_gain"),
    row_var = c("gained_nearest_tss", "many_peaks"),
    col_var = c("gained_median_sig", "gained_nearest_tss"),
    n_genes = nrow(collapse_tbl),
    odds_ratio = c(unname(ft_rules$estimate), unname(ft_peaks$estimate)),
    ci_low = c(ft_rules$conf.int[1], ft_peaks$conf.int[1]),
    ci_high = c(ft_rules$conf.int[2], ft_peaks$conf.int[2]),
    p_value = c(ft_rules$p.value, ft_peaks$p.value),
    stringsAsFactors = FALSE
  )
  write_section_table(out,
                      file.path(out_dir, "mecp2_gene_fisher_summary.tsv"))
  out
}

#' Histogram of annotated MeCP2 peaks per gene.
#'
#' @param collapse_tbl data.frame from build_collapse_table().
#' @param out_dir Section output directory.
#' @return list with the count table, the summary, and the ggplot
run_peaks_per_gene <- function(collapse_tbl, out_dir) {
  n_peaks <- collapse_tbl$mecp2_tss_n_peaks

  dist_tbl <- as.data.frame(table(n_peaks), stringsAsFactors = FALSE)
  colnames(dist_tbl) <- c("n_peaks_per_gene", "n_genes")
  dist_tbl$n_peaks_per_gene <- as.integer(dist_tbl$n_peaks_per_gene)
  dist_tbl <- dist_tbl[order(dist_tbl$n_peaks_per_gene), ]
  dist_tbl$percentage_of_genes <- 100 * dist_tbl$n_genes / nrow(collapse_tbl)
  write_section_table(
    dist_tbl, file.path(out_dir, "mecp2_peaks_per_gene_distribution.tsv"))

  stats_tbl <- data.frame(
    n_genes = length(n_peaks),
    total_peaks = sum(n_peaks),
    min = min(n_peaks),
    q25 = unname(quantile(n_peaks, 0.25)),
    median = median(n_peaks),
    mean = mean(n_peaks),
    q75 = unname(quantile(n_peaks, 0.75)),
    q99 = unname(quantile(n_peaks, 0.99)),
    max = max(n_peaks),
    pct_single_peak = 100 * mean(n_peaks == 1),
    stringsAsFactors = FALSE
  )
  write_section_table(stats_tbl,
                      file.path(out_dir, "mecp2_peaks_per_gene_summary.tsv"))

  cat(sprintf("  Peaks per gene: n = %s genes, median = %.0f, mean = %.2f, max = %.0f\n",
              fmt_count(stats_tbl$n_genes), stats_tbl$median, stats_tbl$mean,
              stats_tbl$max))
  cat(sprintf("  Genes with exactly one peak: %.1f%%\n",
              stats_tbl$pct_single_peak))

  annotation_label <- sprintf(
    "n = %s genes\nmedian = %.0f\nmean = %.2f\nmax = %.0f\none peak only = %.1f%%",
    fmt_count(stats_tbl$n_genes), stats_tbl$median, stats_tbl$mean,
    stats_tbl$max, stats_tbl$pct_single_peak)

  p <- ggplot(data.frame(n_peaks = n_peaks), aes(x = n_peaks)) +
    geom_histogram(binwidth = 1, fill = unname(COLORS$mecp2["MeCP2 Up"]),
                   color = "black", linewidth = 0.2) +
    geom_vline(xintercept = stats_tbl$median, linetype = "dashed",
               color = "black", linewidth = 0.6) +
    annotate("text", x = stats_tbl$max * 0.55, y = stats_tbl$n_genes * 0.4,
             label = annotation_label, hjust = 0, size = 3.6, lineheight = 1.2) +
    scale_y_log10(labels = scales::comma) +
    labs(
      title = "MeCP2 Peaks per Gene",
      subtitle = paste("Annotated peakset, one count per gene symbol.",
                       "Dashed line marks the median. Counts on a log10 axis."),
      x = "Number of MeCP2 peaks assigned to the gene",
      y = "Number of genes (log10 scale)"
    ) +
    theme_emseq()

  save_multiformat_ggplot(p, file.path(out_dir, "60_03i_peaks_per_gene"),
                          width = 10, height = 7)
  list(distribution = dist_tbl, stats = stats_tbl, plot = p)
}

# =============================================================================
# ANALYSIS 6: GSEA TERM RECOVERY
# =============================================================================

#' Build a decreasing ranked gene list keyed by Entrez identifier.
#'
#' Symbols are mapped to Entrez with bitr. When several symbols map to the
#' same identifier, the one with the largest absolute value is kept.
#'
#' @param symbols Character vector of gene symbols.
#' @param values Numeric vector of the same length.
#' @param label Name used in the progress message.
#' @return named numeric vector sorted in decreasing order
build_ranked_list <- function(symbols, values, label) {
  valid <- !is.na(symbols) & nzchar(symbols) & is.finite(values)
  df <- data.frame(SYMBOL = symbols[valid], value = values[valid],
                   stringsAsFactors = FALSE)

  mapping <- clusterProfiler::bitr(unique(df$SYMBOL), fromType = "SYMBOL",
                                   toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  df <- merge(df, mapping, by = "SYMBOL")
  df <- df[order(-abs(df$value)), ]
  df <- df[!duplicated(df$ENTREZID), ]

  ranked <- sort(setNames(df$value, df$ENTREZID), decreasing = TRUE)
  if (length(ranked) < 500) {
    stop("Ranked list '", label, "' has only ", length(ranked),
         " genes after Entrez mapping. GSEA needs a genome-scale ranking.")
  }

  cat(sprintf("  %-32s %s genes ranked (%.3f to %.3f)\n", label,
              fmt_count(length(ranked)), max(ranked), min(ranked)))
  ranked
}

#' Run GO BP GSEA on one ranked list and mark neuronal terms.
#'
#' pvalueCutoff is 1 so the returned table holds every tested term. The
#' is_significant column applies the section FDR cutoff.
#'
#' @param ranked Named numeric vector from build_ranked_list().
#' @param label Ranking name recorded in the output.
#' @param opt Parsed options carrying the GSEA sizes, seed, and FDR cutoff.
#' @return data.frame of GSEA results with ranking, is_significant, is_neuronal
run_gsea <- function(ranked, label, opt) {
  cat(sprintf("\n  GSEA on %s (%s genes)...\n", label,
              fmt_count(length(ranked))))
  set.seed(opt$seed)

  gsea <- clusterProfiler::gseGO(
    geneList = ranked,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    minGSSize = opt$gsea_min_size,
    maxGSSize = opt$gsea_max_size,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE,
    seed = TRUE
  )

  res <- as.data.frame(gsea@result)
  if (nrow(res) == 0) {
    stop("gseGO returned no tested GO BP terms for the ranking '", label,
         "'. Check the ranked list and the gene set size limits.")
  }

  res$ranking <- label
  res$is_significant <- !is.na(res$p.adjust) & res$p.adjust < opt$fdr_threshold
  res$is_neuronal <- grepl(NEURONAL_PATTERN, res$Description, ignore.case = TRUE)
  res$n_genes_ranked <- length(ranked)

  cat(sprintf("    %s terms tested, %s significant at q < %.3f, %s of those neuronal\n",
              fmt_count(nrow(res)), fmt_count(sum(res$is_significant)),
              opt$fdr_threshold,
              fmt_count(sum(res$is_significant & res$is_neuronal))))
  res
}

#' Fisher tests of the neuronal term share between every pair of rankings.
#'
#' These are term-level tests, so they call fisher.test() directly rather than
#' the gene-level registry.
#'
#' @param comparison data.frame from run_gsea_comparison().
#' @param out_dir Section output directory.
#' @return data.frame with one row per ranking pair
run_gsea_pairwise_tests <- function(comparison, out_dir) {
  pairs <- utils::combn(seq_len(nrow(comparison)), 2, simplify = FALSE)

  out <- do.call(rbind, lapply(pairs, function(idx) {
    a <- comparison[idx[1], ]
    b <- comparison[idx[2], ]
    tab <- matrix(c(a$n_significant_neuronal, b$n_significant_neuronal,
                    a$n_significant_other, b$n_significant_other),
                  nrow = 2)
    ft <- fisher.test(tab)
    data.frame(
      ranking_1 = a$ranking, ranking_2 = b$ranking,
      neuronal_1 = a$n_significant_neuronal, other_1 = a$n_significant_other,
      neuronal_2 = b$n_significant_neuronal, other_2 = b$n_significant_other,
      odds_ratio = unname(ft$estimate),
      ci_low = ft$conf.int[1], ci_high = ft$conf.int[2],
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  }))
  write_section_table(out,
                      file.path(out_dir, "gsea_neuronal_share_fisher.tsv"))

  for (i in seq_len(nrow(out))) {
    cat(sprintf("  Fisher neuronal share, %s vs %s: OR = %.3f, p = %.3g\n",
                out$ranking_1[i], out$ranking_2[i], out$odds_ratio[i],
                out$p_value[i]))
  }
  out
}

#' Run every GSEA ranking, compare neuronal term recovery, and draw the bar.
#'
#' @param collapse_tbl data.frame from build_collapse_table().
#' @param k119ub_signal data.frame from load_k119ub_gene_signal().
#' @param opt Parsed options.
#' @param out_dir Section output directory.
#' @return list with the comparison table, the pairwise tests, and the ggplot
run_gsea_comparison <- function(collapse_tbl, k119ub_signal, opt, out_dir) {
  ranked_lists <- list(
    "MeCP2 fold (nearest TSS)" = build_ranked_list(
      collapse_tbl$gene_name, collapse_tbl$mecp2_tss_fold,
      "MeCP2 fold (nearest TSS)"),
    "MeCP2 fold (median significant)" = build_ranked_list(
      collapse_tbl$gene_name, collapse_tbl$mecp2_med_fold,
      "MeCP2 fold (median significant)"),
    "K119ub gene-body log2FC" = build_ranked_list(
      k119ub_signal$symbol, k119ub_signal$gb_log2fc,
      "K119ub gene-body log2FC")
  )

  file_stub <- c(
    "MeCP2 fold (nearest TSS)" = "gsea_mecp2_nearest_tss_go_bp.tsv",
    "MeCP2 fold (median significant)" = "gsea_mecp2_median_significant_go_bp.tsv",
    "K119ub gene-body log2FC" = "gsea_k119ub_gene_signal_go_bp.tsv"
  )

  results <- lapply(names(ranked_lists), function(label) {
    res <- run_gsea(ranked_lists[[label]], label, opt)
    write_section_table(res, file.path(out_dir, file_stub[[label]]))
    res
  })
  names(results) <- names(ranked_lists)

  comparison <- do.call(rbind, lapply(names(results), function(label) {
    res <- results[[label]]
    sig <- res[res$is_significant, , drop = FALSE]
    data.frame(
      ranking = label,
      n_genes_ranked = length(ranked_lists[[label]]),
      n_terms_tested = nrow(res),
      n_significant = nrow(sig),
      n_significant_neuronal = sum(sig$is_neuronal),
      n_significant_other = sum(!sig$is_neuronal),
      n_neuronal_positive_nes = sum(sig$is_neuronal & sig$NES > 0),
      n_neuronal_negative_nes = sum(sig$is_neuronal & sig$NES < 0),
      pct_significant_neuronal = 100 * sum(sig$is_neuronal) / max(nrow(sig), 1),
      stringsAsFactors = FALSE
    )
  }))
  write_section_table(comparison,
                      file.path(out_dir, "gsea_term_comparison.tsv"))

  for (i in seq_len(nrow(comparison))) {
    cat(sprintf("  %-32s %s significant terms, %s neuronal (%.1f%%)\n",
                comparison$ranking[i], fmt_count(comparison$n_significant[i]),
                fmt_count(comparison$n_significant_neuronal[i]),
                comparison$pct_significant_neuronal[i]))
  }

  neuronal_terms <- do.call(rbind, lapply(results, function(res) {
    keep <- res$is_significant & res$is_neuronal
    res[keep, c("ranking", "ID", "Description", "setSize", "NES",
                "pvalue", "p.adjust"), drop = FALSE]
  }))
  write_section_table(neuronal_terms,
                      file.path(out_dir, "gsea_neuronal_terms.tsv"))
  cat(sprintf("  Significant neuronal terms across all rankings: %s\n",
              fmt_count(nrow(neuronal_terms))))

  pairwise <- run_gsea_pairwise_tests(comparison, out_dir)

  long <- comparison %>%
    dplyr::select(ranking, n_significant_neuronal, n_significant_other) %>%
    tidyr::pivot_longer(cols = -ranking, names_to = "category",
                        values_to = "count") %>%
    dplyr::mutate(
      category = ifelse(category == "n_significant_neuronal",
                        "Neuronal terms", "Other terms"),
      category = factor(category, levels = c("Neuronal terms", "Other terms")),
      ranking = factor(ranking, levels = comparison$ranking)
    ) %>%
    as.data.frame()

  totals <- comparison
  totals$ranking <- factor(totals$ranking, levels = comparison$ranking)

  p <- ggplot(long, aes(x = ranking, y = count, fill = category)) +
    geom_col(position = "stack", width = 0.62, color = "black",
             linewidth = 0.3) +
    geom_text(aes(label = fmt_count(count)),
              position = position_stack(vjust = 0.5), size = 4,
              color = "white", fontface = "bold") +
    geom_text(data = totals,
              aes(x = ranking, y = n_significant,
                  label = sprintf("%s significant\n%.1f%% neuronal",
                                  fmt_count(n_significant),
                                  pct_significant_neuronal)),
              inherit.aes = FALSE, vjust = -0.3, size = 3.4, lineheight = 1.1) +
    scale_fill_manual(values = c("Neuronal terms" = "#756BB1",
                                 "Other terms" = "grey70"), name = NULL) +
    scale_x_discrete(labels = function(x) gsub(" \\(", "\n(", x)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.20)),
                       labels = scales::comma) +
    labs(
      title = "Significant GSEA GO BP Terms by Ranking Metric",
      subtitle = sprintf(
        "Neuronal terms match '%s' in the term description | q < %.3f",
        NEURONAL_PATTERN, opt$fdr_threshold),
      x = "Gene ranking metric", y = "Number of significant GO BP terms"
    ) +
    theme_emseq() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p, file.path(out_dir, "60_03j_gsea_term_comparison"),
                          width = 11, height = 8)
  list(comparison = comparison, pairwise = pairwise,
       neuronal_terms = neuronal_terms, plot = p)
}

# =============================================================================
# COMPOSITE FIGURE
# =============================================================================

#' Assemble the section panels into one composite figure.
#'
#' @param plots Named list of ggplot objects.
#' @param out_dir Section output directory.
#' @return the patchwork object
build_composite <- function(plots, out_dir) {
  composite <- (plots$counts | plots$fold | plots$width) /
    (plots$weighted) /
    (plots$collapse | plots$peaks_per_gene | plots$gsea) +
    plot_annotation(
      title = "Section 60_03: MeCP2 Peak Gains Against Gene-Level MeCP2 Loss",
      subtitle = paste("Peak counting, peak weighting, the gene-level collapse",
                       "rule, and GO term recovery"),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40")
      )
    )

  save_multiformat_ggplot(composite, file.path(out_dir, "60_03k_composite"),
                          width = 26, height = 24)
  composite
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
  cat("SECTION 60_03: MeCP2 PEAK GAINS vs GENE-LEVEL MeCP2 LOSS\n")
  cat("================================================================================\n")
  cat("Output dir:    ", OUT_DIR, "\n", sep = "")
  cat("FDR threshold: ", fdr_threshold, "\n", sep = "")
  cat("GSEA set size: ", opt$gsea_min_size, " to ", opt$gsea_max_size, "\n", sep = "")
  cat("Seed:          ", opt$seed, "\n\n", sep = "")

  stopifnot(
    "MeCP2 annotated peak file not found" = file.exists(MECP2_PATHS$annotated),
    "MeCP2 up peak BED not found" = file.exists(MECP2_PATHS$up),
    "MeCP2 down peak BED not found" = file.exists(MECP2_PATHS$down),
    "K119ub gene signal table not found" =
      file.exists(DIFFBIND_PATHS$k119ub_gene_signal)
  )

  cat("--- Loading inputs ---\n")
  annotated <- load_mecp2_annotated(MECP2_PATHS$annotated)
  annotated_peaks <- prepare_peaks(annotated, fdr_threshold, PEAKSET_ANNOTATED)
  diffbind_peaks <- prepare_peaks(mecp2_diffbind, fdr_threshold,
                                  PEAKSET_DIFFBIND)
  peaks <- rbind(annotated_peaks, diffbind_peaks)

  mecp2_up_gr <- load_chip_peaks(MECP2_PATHS$up, "MeCP2 up")
  mecp2_down_gr <- load_chip_peaks(MECP2_PATHS$down, "MeCP2 down")
  k119ub_signal <- load_k119ub_gene_signal(DIFFBIND_PATHS$k119ub_gene_signal)

  write_section_table(peaks,
                      file.path(OUT_DIR, "mecp2_peak_level_classified.tsv"))

  cat("\n--- Analysis 1: peak counts by direction ---\n")
  counts_res <- plot_direction_counts(peaks, OUT_DIR)

  cat("\n--- Analysis 2: peak fold, width, and concentration by direction ---\n")
  property_res <- run_property_analysis(peaks, OUT_DIR)
  width_fisher <- run_peak_width_fisher(peaks, OUT_DIR)

  cat("\n--- Analysis 3: weighted gained against lost ---\n")
  weighted_res <- run_weighted_comparison(peaks, OUT_DIR)

  cat("\n--- Analysis 4: signal mass balance ---\n")
  mass_res <- run_signal_mass_balance(peaks, OUT_DIR)

  cat("\n--- Analysis 5: peakset interval overlap ---\n")
  run_interval_overlap(peaks, mecp2_up_gr, mecp2_down_gr, mecp2_consensus,
                       OUT_DIR)

  cat("\n--- Analysis 6: gene-level collapse rules ---\n")
  annotated_sym <- filter_to_symbol(annotated)
  collapse_tbl <- build_collapse_table(annotated_sym, fdr_threshold)
  write_section_table(
    collapse_tbl, file.path(OUT_DIR, "mecp2_gene_level_collapse_rules.tsv"))

  scatter_res <- run_collapse_scatter(collapse_tbl, OUT_DIR)
  calls_res <- run_collapse_gene_calls(collapse_tbl, OUT_DIR)

  cat("\n--- Analysis 7: gene-level Fisher tests ---\n")
  fisher_tbl <- run_gene_fisher_tests(collapse_tbl, OUT_DIR)

  cat("\n--- Analysis 8: peaks per gene ---\n")
  peaks_per_gene_res <- run_peaks_per_gene(collapse_tbl, OUT_DIR)

  cat("\n--- Analysis 9: GSEA term recovery ---\n")
  gsea_res <- run_gsea_comparison(collapse_tbl, k119ub_signal, opt, OUT_DIR)

  cat("\n--- Composite figure ---\n")
  build_composite(list(
    counts = counts_res$plot,
    fold = property_res$fold_plot,
    width = property_res$width_plot,
    weighted = weighted_res$plot,
    collapse = scatter_res$plot,
    peaks_per_gene = peaks_per_gene_res$plot,
    gsea = gsea_res$plot
  ), OUT_DIR)

  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 60_03 SUMMARY\n")
  cat("================================================================================\n")
  for (i in seq_len(nrow(counts_res$ratios))) {
    cat(sprintf("%-24s gained %s / lost %s (ratio %.2f)\n",
                as.character(counts_res$ratios$peakset[i]),
                fmt_count(counts_res$ratios$Gained[i]),
                fmt_count(counts_res$ratios$Lost[i]),
                counts_res$ratios$gained_over_lost[i]))
  }
  for (i in seq_len(nrow(weighted_res$wide))) {
    cat(sprintf("%-24s gained share, %-26s %.1f%%\n",
                as.character(weighted_res$wide$peakset[i]),
                as.character(weighted_res$wide$scheme[i]),
                weighted_res$wide$gained_share_pct[i]))
  }
  net_rows <- mass_res$balance[mass_res$balance$class == "Net (all peaks)", ]
  for (i in seq_len(nrow(net_rows))) {
    cat(sprintf("%-24s net signal mass change %s (%.2f%% of control mass)\n",
                as.character(net_rows$peakset[i]),
                fmt_big(net_rows$total_mass_delta[i]), net_rows$pct_change[i]))
  }
  width_rows <- property_res$wilcoxon[property_res$wilcoxon$property == "width", ]
  for (i in seq_len(nrow(width_rows))) {
    cat(sprintf("%-24s peak width gained vs lost: median %s vs %s, p = %.3g\n",
                width_rows$peakset[i], fmt_big(width_rows$median_gained[i]),
                fmt_big(width_rows$median_lost[i]), width_rows$p_value[i]))
  }
  for (i in seq_len(nrow(width_fisher))) {
    cat(sprintf("%-24s wide-peak Fisher: OR = %.3f, p = %.3g\n",
                width_fisher$peakset[i], width_fisher$odds_ratio[i],
                width_fisher$p_value[i]))
  }
  cat(sprintf("Genes with MeCP2 peaks:           %s\n",
              fmt_count(nrow(collapse_tbl))))
  cat(sprintf("Collapse-rule Spearman rho:       %.3f (sign differs for %.1f%% of genes)\n",
              scatter_res$summary$spearman_rho,
              scatter_res$summary$pct_sign_flip))
  for (i in seq_len(nrow(calls_res$counts))) {
    cat(sprintf("Gene calls, %-19s %-10s %s\n",
                as.character(calls_res$counts$rule[i]),
                as.character(calls_res$counts$call[i]),
                fmt_count(calls_res$counts$n_genes[i])))
  }
  for (i in seq_len(nrow(fisher_tbl))) {
    cat(sprintf("Fisher %-30s OR = %.3f (95%% CI %.3f to %.3f), p = %.3g\n",
                fisher_tbl$test_id[i], fisher_tbl$odds_ratio[i],
                fisher_tbl$ci_low[i], fisher_tbl$ci_high[i],
                fisher_tbl$p_value[i]))
  }
  for (i in seq_len(nrow(gsea_res$comparison))) {
    cat(sprintf("GSEA %-32s %s significant, %s neuronal\n",
                gsea_res$comparison$ranking[i],
                fmt_count(gsea_res$comparison$n_significant[i]),
                fmt_count(gsea_res$comparison$n_significant_neuronal[i])))
  }
  cat("\nSection 60_03 complete.\n\n")
}

main()
