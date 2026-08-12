# scripts/sections/70_neuronal/70_04_synapse_chromatin.R
#
# Section 70_04: chromatin remodeling at neuronal and synapse genes that also
# gain MeCP2 binding.
#
# Asks whether genes that belong to two sets at once carry larger chromatin
# changes than genes that belong to neither. The two sets are a gene-identity
# set (neuronal, then synapse and axon) and the set of genes that gain MeCP2
# binding in BAP1-KO. Chromatin change is the gene-level DiffBind log2 fold
# change of ATAC, H3K27ac, H3K27me3 and H2AK119ub, taken from the section 20_02
# handoff.
#
# The Biomodal original built a TRIPLE overlap: neuronal genes, MeCP2-up genes,
# and genes with a coordinated 5mC-up / 5hmC-down change. That third criterion
# needs two separate methylation modalities. EM-seq measures non-CG methylation
# as one modality, so the criterion does not exist here and the intersection in
# this section is DOUBLE: gene set AND MeCP2-up.
#
# Four groups partition every tested gene, once against the neuronal set and
# again against the synapse set:
#   Both            in the gene set and MeCP2-up
#   Gene set only   in the gene set, not MeCP2-up
#   MeCP2-up only   MeCP2-up, not in the gene set
#   Neither         in neither set; the reference group of every comparison
#
# MeCP2-up means the peak nearest the TSS of the gene gains binding at the
# section FDR threshold. A gene with no MeCP2 peak is not MeCP2-up.
#
# Long genes could produce a group difference without any gene-identity effect,
# so the whole four-group comparison repeats inside equal-count gene-length
# quintiles. A quintile cell holding fewer than MIN_STRATUM_SIZE genes in either
# arm of a comparison is written to the output table with test columns NA and a
# reason string, rather than tested.
#
# Reads:
#   HANDOFF_PATHS$neuronal_gene_set     neuronal gene set, written by 70_01
#   GENESET_PATHS$synapse               synapse and axon gene set, placed under
#                                       data/neuronal/ by copy_reference_data.sh
#   HANDOFF_PATHS$gene_level_all_marks  gene-level mark folds, written by 20_02
#   mecp2_diffbind                      pre-loaded by _shared_config.R
#
# Writes into --output-dir (default results/sections/70_neuronal/):
#   70_04a .. 70_04h              figures, each beside the TSV of its numbers
#   70_04_per_gene_groups.tsv     per-gene membership, folds and strata
#   fisher_tables/                gene tables behind the registered Fisher tests
#
# Adapted from Biomodal section 76 (triple overlap and synapse chromatin).

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)

# =============================================================================
# SECTION CONSTANTS
# =============================================================================

SECTION_ID <- "70_04"

MARK_ORDER <- c("atac", "k27ac", "k27me3", "k119ub")

# display  name used in figures and tables
# color    accent colour, taken from the shared COLORS palettes
MARK_META <- list(
  atac   = list(display = "ATAC-seq",  color = COLORS$atac[["ATAC Up"]]),
  k27ac  = list(display = "H3K27ac",   color = COLORS$h3k27ac[["H3K27ac Gained"]]),
  k27me3 = list(display = "H3K27me3",  color = COLORS$h3k27me3[["H3K27me3 Gained"]]),
  k119ub = list(display = "H2AK119ub", color = COLORS$k119ub[["K119ub Gained"]])
)

MARK_DISPLAY <- vapply(MARK_META[MARK_ORDER], `[[`, character(1), "display")

# The four groups every comparison uses. "Both" is the double overlap that
# replaces the Biomodal triple overlap.
GROUP_ORDER <- c("Both", "Gene set only", "MeCP2-up only", "Neither")
REFERENCE_GROUP <- "Neither"
COMPARISON_GROUPS <- setdiff(GROUP_ORDER, REFERENCE_GROUP)

GROUP_COLORS <- c(
  "Both"          = COLORS$mecp2[["MeCP2 Up"]],
  "Gene set only" = COLORS$k119ub[["K119ub Gained"]],
  "MeCP2-up only" = COLORS$atac[["ATAC Up"]],
  "Neither"       = "grey75"
)

# Two-line axis labels, so the group names fit under a narrow violin.
GROUP_AXIS_LABELS <- c(
  "Both"          = "Both",
  "Gene set only" = "Gene set\nonly",
  "MeCP2-up only" = "MeCP2-up\nonly",
  "Neither"       = "Neither"
)

# The two gene-identity sets, reported side by side everywhere.
SET_SLUGS <- c("neuronal", "synapse")
SET_LABELS <- c(neuronal = "Neuronal gene set", synapse = "Synapse gene set")

# Columns the section 20_02 handoff must carry for this section to run.
REQUIRED_MASTER_COLUMNS <- c(
  "gene_name", "gene_id", "chr", "start", "end", "gene_length",
  "mch_diff", "edger_logFC", "edger_fdr", "mch_sig", "mch_direction",
  paste0(MARK_ORDER, "_fold"), paste0(MARK_ORDER, "_fdr")
)

N_DECILES <- 10
N_LENGTH_QUINTILES <- 5

# A group holding fewer genes than this for a mark stops the main comparison.
MIN_GROUP_SIZE <- 10

# A gene-length quintile cell holding fewer genes than this is recorded as
# untested instead of tested.
MIN_STRATUM_SIZE <- 10

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
                help = paste("FDR cutoff for MeCP2-up calling and for the",
                             "per-mark remodelled flag [default: %default]"))
  )
  opt <- parse_args(OptionParser(option_list = option_list))

  if (opt$fdr_threshold <= 0 || opt$fdr_threshold >= 1) {
    stop("--fdr-threshold must be between 0 and 1, got ", opt$fdr_threshold)
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

significance_stars <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", ""))))
}

#' Summarise n, median and mean of one value column for every cell of an
#' arbitrary set of grouping columns, with the on-plot label text.
summarise_cells <- function(df, group_cols, value_col, digits = 3) {
  df %>%
    dplyr::filter(!is.na(.data[[value_col]])) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median = median(.data[[value_col]]),
      mean = mean(.data[[value_col]]),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = sprintf("n = %s\nmed = %s",
                      format(n, big.mark = ",", trim = TRUE),
                      format(round(median, digits), nsmall = digits, trim = TRUE))
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

#' Per-row y positions for on-plot annotation, one value per facet row.
#'
#' facet_grid(mark ~ ., scales = "free_y") gives every mark its own y range,
#' shared by the columns of that row. The label sits above the row maximum and
#' the test annotation below the row minimum.
mark_label_positions <- function(long) {
  long %>%
    dplyr::filter(!is.na(fold)) %>%
    dplyr::group_by(mark) %>%
    dplyr::summarise(row_min = min(fold), row_max = max(fold), .groups = "drop") %>%
    dplyr::mutate(
      span = row_max - row_min,
      label_y = row_max + 0.10 * span,
      test_y = row_min - 0.07 * span
    ) %>%
    as.data.frame()
}

# =============================================================================
# STEP 1: INPUTS
# =============================================================================

#' Read the gene-level all-marks table written by section 20_02.
load_gene_level_marks <- function() {
  path <- HANDOFF_PATHS$gene_level_all_marks
  if (!file.exists(path)) {
    stop("Gene-level all-marks table not found: ", path,
         "\nRun section 20_02 (20_chip_integration/20_02_multi_mark_diffbind.R) ",
         "first. It writes this file.")
  }

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")
  missing <- setdiff(REQUIRED_MASTER_COLUMNS, colnames(df))
  if (length(missing) > 0) {
    stop("The section 20_02 handoff is missing columns: ",
         paste(missing, collapse = ", "), " in ", path)
  }
  if (nrow(df) == 0) {
    stop("The section 20_02 handoff holds no rows: ", path)
  }

  df$mch_sig <- as.logical(df$mch_sig)
  cat(sprintf("  Gene-level all-marks table: %s genes\n",
              format(nrow(df), big.mark = ",")))
  for (mark in MARK_ORDER) {
    n_with <- sum(!is.na(df[[paste0(mark, "_fold")]]))
    cat(sprintf("    %-10s %s genes carry a peak (%.1f%%)\n",
                MARK_META[[mark]]$display, format(n_with, big.mark = ","),
                100 * n_with / nrow(df)))
  }
  df
}

#' Read the neuronal gene set handoff written by section 70_01.
load_neuronal_gene_set <- function() {
  path <- HANDOFF_PATHS$neuronal_gene_set
  if (!file.exists(path)) {
    stop("Neuronal gene set handoff not found: ", path,
         "\nRun section 70_01 (70_neuronal/70_01_k119ub_neuronal.R) first. ",
         "It writes this file.")
  }

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")
  if (!"gene" %in% colnames(df)) {
    stop("The section 70_01 handoff has no 'gene' column: ", path,
         "\nSection 70_01 writes gene, source, is_derived, is_external, ",
         "k119ub_ctrl_signal, k119ub_log2fc, k119ub_decile.")
  }

  genes <- unique(df$gene[!is.na(df$gene) & nzchar(df$gene)])
  if (length(genes) == 0) {
    stop("The section 70_01 handoff holds no gene symbol: ", path)
  }
  cat(sprintf("  Neuronal gene set: %s genes\n",
              format(length(genes), big.mark = ",")))
  genes
}

#' Read the synapse and axon gene set placed under data/neuronal/.
load_synapse_gene_set <- function() {
  path <- GENESET_PATHS$synapse
  if (!file.exists(path)) {
    stop("Synapse and axon gene set not found: ", path,
         "\nRun scripts/utils/copy_reference_data.sh to place it under ",
         "data/neuronal/.")
  }
  load_gene_set(path, "Synapse and axon gene set")
}

#' Collapse the MeCP2 DiffBind peakset to one fold change per gene.
#'
#' Peaks are restricted to canonical chromosomes, annotated to genes with
#' ChIPseeker, then collapsed by the nearest-TSS rule: the peak closest to the
#' TSS supplies the gene-level fold change and FDR.
derive_mecp2_gene_table <- function(fdr_threshold) {
  cat("--- Step 1d: MeCP2 peaks to gene-level fold changes ---\n")

  db <- mecp2_diffbind
  n_all <- nrow(db)
  db <- db[db$Chr %in% CANONICAL_CHRS, , drop = FALSE]
  db <- db[!is.na(db$Fold) & !is.na(db$FDR), , drop = FALSE]
  if (nrow(db) == 0) {
    stop("No usable MeCP2 peaks after filtering to canonical chromosomes and ",
         "complete Fold and FDR values.")
  }
  cat(sprintf("  %s MeCP2 peaks in, %s kept on canonical chromosomes\n",
              format(n_all, big.mark = ","), format(nrow(db), big.mark = ",")))

  annotated <- annotate_peaks_to_genes(db, "MeCP2")
  usable <- !is.na(annotated$SYMBOL) & nzchar(annotated$SYMBOL) &
    !is.na(annotated$Fold) & !is.na(annotated$FDR)
  cat(sprintf("  %s of %s annotated peaks carry a gene symbol\n",
              format(sum(usable), big.mark = ","),
              format(nrow(annotated), big.mark = ",")))
  annotated <- annotated[usable, , drop = FALSE]

  gene_table <- aggregate_diffbind_by_gene(annotated, method = "nearest_tss",
                                           fdr_threshold = fdr_threshold,
                                           prefix = "mecp2")
  if (anyDuplicated(gene_table$gene_name) > 0) {
    stop("aggregate_diffbind_by_gene() returned duplicate genes for MeCP2.")
  }

  n_up <- sum(gene_table$mecp2_fdr < fdr_threshold & gene_table$mecp2_fold > 0,
              na.rm = TRUE)
  n_down <- sum(gene_table$mecp2_fdr < fdr_threshold & gene_table$mecp2_fold < 0,
                na.rm = TRUE)
  cat(sprintf("  %s genes carry a MeCP2 peak: %s up, %s down at FDR < %.3g\n\n",
              format(nrow(gene_table), big.mark = ","),
              format(n_up, big.mark = ","), format(n_down, big.mark = ","),
              fdr_threshold))
  gene_table
}

# =============================================================================
# STEP 2: THE GENE UNIVERSE
# =============================================================================

#' Keep one row per gene symbol, the row with the largest absolute edgeR log
#' fold change. Set membership and MeCP2 aggregation both key on the symbol.
deduplicate_by_gene_name <- function(master) {
  if (anyNA(master$edger_logFC)) {
    stop(sum(is.na(master$edger_logFC)), " genes carry no edger_logFC, so the ",
         "duplicate gene symbols cannot be resolved by effect size.")
  }
  n_dup <- sum(duplicated(master$gene_name))
  ord <- order(master$gene_name, -abs(master$edger_logFC))
  keep <- sort(ord[!duplicated(master$gene_name[ord])])
  cat(sprintf("  %s genes in, %s kept after collapsing %s duplicate symbols\n",
              format(nrow(master), big.mark = ","),
              format(length(keep), big.mark = ","),
              format(n_dup, big.mark = ",")))
  master[keep, , drop = FALSE]
}

#' Join MeCP2 to the master table and add every derived flag and stratum.
build_universe <- function(master, mecp2_gene_table, neuronal_genes,
                           synapse_genes, fdr_threshold) {
  cat("--- Step 2: building the gene universe ---\n")

  universe <- deduplicate_by_gene_name(master)
  n_before <- nrow(universe)
  universe <- dplyr::left_join(universe, mecp2_gene_table, by = "gene_name")
  if (nrow(universe) != n_before) {
    stop("The MeCP2 join changed the row count: ", n_before, " genes in, ",
         nrow(universe), " out.")
  }

  universe$is_mecp2_up <- !is.na(universe$mecp2_fold) &
    !is.na(universe$mecp2_fdr) &
    universe$mecp2_fdr < fdr_threshold & universe$mecp2_fold > 0
  universe$is_neuronal <- universe$gene_name %in% neuronal_genes
  universe$is_synapse <- universe$gene_name %in% synapse_genes

  for (mark in MARK_ORDER) {
    fold <- universe[[paste0(mark, "_fold")]]
    fdr <- universe[[paste0(mark, "_fdr")]]
    remodeled <- fdr < fdr_threshold
    remodeled[is.na(fold) | is.na(fdr)] <- NA
    universe[[paste0(mark, "_remodeled")]] <- remodeled
  }

  if (any(!is.finite(universe$gene_length)) || any(universe$gene_length <= 0)) {
    stop("A gene carries a non-finite or non-positive gene_length, so the ",
         "gene-length quintiles cannot be built.")
  }
  universe$length_quintile <- cut_into_bins(universe$gene_length,
                                            N_LENGTH_QUINTILES, "gene length")

  # K119ub deciles rank the genes that carry a K119ub peak by their gene-level
  # log2 fold change. D1 is the strongest loss, D10 the strongest gain.
  has_k119ub <- !is.na(universe$k119ub_fold)
  if (sum(has_k119ub) < N_DECILES * MIN_GROUP_SIZE) {
    stop("Only ", sum(has_k119ub), " genes carry a K119ub fold change. Too few ",
         "for ", N_DECILES, " deciles.")
  }
  universe$k119ub_decile <- NA_integer_
  universe$k119ub_decile[has_k119ub] <- cut_into_bins(
    universe$k119ub_fold[has_k119ub], N_DECILES, "K119ub log2 fold change")

  cat(sprintf("  Universe: %s genes\n", format(nrow(universe), big.mark = ",")))
  cat(sprintf("  MeCP2-up:    %s genes (%.1f%%)\n",
              format(sum(universe$is_mecp2_up), big.mark = ","),
              100 * mean(universe$is_mecp2_up)))
  cat(sprintf("  Neuronal:    %s genes (%.1f%%)\n",
              format(sum(universe$is_neuronal), big.mark = ","),
              100 * mean(universe$is_neuronal)))
  cat(sprintf("  Synapse:     %s genes (%.1f%%)\n",
              format(sum(universe$is_synapse), big.mark = ","),
              100 * mean(universe$is_synapse)))
  cat(sprintf("  With a K119ub decile: %s genes\n\n",
              format(sum(has_k119ub), big.mark = ",")))
  universe
}

#' Add the four-group assignment for one gene-identity set.
build_set_table <- function(universe, slug, set_genes) {
  df <- universe
  df$set_slug <- slug
  df$set_label <- factor(SET_LABELS[[slug]], levels = unname(SET_LABELS))
  df$in_set <- df$gene_name %in% set_genes
  df$group <- factor(
    dplyr::case_when(
      df$in_set & df$is_mecp2_up ~ "Both",
      df$in_set                  ~ "Gene set only",
      df$is_mecp2_up             ~ "MeCP2-up only",
      TRUE                       ~ "Neither"
    ),
    levels = GROUP_ORDER
  )
  df$is_both <- df$group == "Both"
  df
}

#' Long form of the four mark folds, one row per gene and mark.
build_mark_long <- function(combined) {
  long <- do.call(rbind, lapply(MARK_ORDER, function(mark) {
    data.frame(
      gene_name = combined$gene_name,
      chr = combined$chr,
      set_slug = combined$set_slug,
      set_label = combined$set_label,
      group = combined$group,
      length_quintile = combined$length_quintile,
      mark = MARK_META[[mark]]$display,
      fold = combined[[paste0(mark, "_fold")]],
      stringsAsFactors = FALSE
    )
  }))
  long$mark <- factor(long$mark, levels = unname(MARK_DISPLAY))
  long$group <- factor(long$group, levels = GROUP_ORDER)
  long[!is.na(long$fold), , drop = FALSE]
}

# =============================================================================
# FIGURE 70_04a: GROUP SIZES
# =============================================================================

summarise_group_sizes <- function(combined) {
  sizes <- combined %>%
    dplyr::group_by(set_label, group) %>%
    dplyr::summarise(n_genes = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(set_label) %>%
    dplyr::mutate(frac_of_universe = n_genes / sum(n_genes)) %>%
    dplyr::ungroup() %>%
    as.data.frame()

  for (mark in MARK_ORDER) {
    col <- paste0(mark, "_fold")
    counts <- combined %>%
      dplyr::group_by(set_label, group) %>%
      dplyr::summarise(n = sum(!is.na(.data[[col]])), .groups = "drop") %>%
      as.data.frame()
    sizes[[paste0("n_with_", mark)]] <- counts$n[
      match(paste(sizes$set_label, sizes$group),
            paste(counts$set_label, counts$group))]
  }
  sizes
}

plot_group_sizes <- function(sizes, out_dir) {
  cat("--- Figure 70_04a: group sizes ---\n")

  for (i in seq_len(nrow(sizes))) {
    cat(sprintf("  %-18s %-14s %s genes (%.1f%%)\n",
                as.character(sizes$set_label[i]), as.character(sizes$group[i]),
                format(sizes$n_genes[i], big.mark = ","),
                100 * sizes$frac_of_universe[i]))
  }

  plot_df <- sizes
  plot_df$bar_label <- sprintf("%s\n%.1f%%",
                               format(plot_df$n_genes, big.mark = ",",
                                      trim = TRUE),
                               100 * plot_df$frac_of_universe)

  figure <- ggplot(plot_df, aes(x = group, y = n_genes, fill = group)) +
    geom_col(width = 0.72, show.legend = FALSE) +
    geom_text(aes(label = bar_label), vjust = -0.25, size = 3.2,
              lineheight = 1.05) +
    facet_wrap(~ set_label) +
    scale_fill_manual(values = GROUP_COLORS) +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0, 0.18))) +
    scale_x_discrete(labels = GROUP_AXIS_LABELS) +
    labs(
      title = "Genes by gene-set and MeCP2 membership",
      subtitle = paste("Both = in the gene set and gaining MeCP2 binding;",
                       "Neither is the reference group of every comparison"),
      x = "", y = "Genes"
    ) +
    theme_emseq() +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"),
          axis.text.x = element_text(size = 9))

  save_multiformat_ggplot(figure, file.path(out_dir, "70_04a_group_sizes"),
                          width = 11, height = 6)
  write_section_table(sizes, file.path(out_dir, "70_04a_group_sizes.tsv"))
  cat("\n")
  figure
}

# =============================================================================
# STEP 3: MARK FOLD CHANGES ACROSS THE FOUR GROUPS
# =============================================================================

#' Stop when any group holds too few genes with a fold change for a mark.
check_group_sizes <- function(long) {
  cells <- long %>%
    dplyr::group_by(set_label, mark, group) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    as.data.frame()

  expected <- nrow(expand.grid(set = unname(SET_LABELS),
                               mark = unname(MARK_DISPLAY),
                               group = GROUP_ORDER))
  if (nrow(cells) != expected) {
    stop("Expected ", expected, " set-by-mark-by-group cells with data, found ",
         nrow(cells), ". A group holds no gene with a fold change for a mark.")
  }

  small <- cells[cells$n < MIN_GROUP_SIZE, , drop = FALSE]
  if (nrow(small) > 0) {
    stop(nrow(small), " set-by-mark-by-group cells hold fewer than ",
         MIN_GROUP_SIZE, " genes, so the four-group comparison has no power. ",
         "First cell: ", small$set_label[1], " / ", small$mark[1], " / ",
         small$group[1], " with n = ", small$n[1], ".")
  }
  invisible(cells)
}

#' Kruskal-Wallis across the four groups, per gene set and mark.
test_kruskal_by_mark <- function(long) {
  rows <- lapply(unname(SET_LABELS), function(set_label) {
    do.call(rbind, lapply(unname(MARK_DISPLAY), function(mark) {
      sub <- long[long$set_label == set_label & long$mark == mark, , drop = FALSE]
      test <- kruskal.test(fold ~ group, data = sub)
      data.frame(
        set_label = set_label,
        mark = mark,
        n_genes = nrow(sub),
        chi_squared = unname(test$statistic),
        df = unname(test$parameter),
        p_value = test$p.value,
        stringsAsFactors = FALSE
      )
    }))
  })
  out <- do.call(rbind, rows)
  out$q_value <- p.adjust(out$p_value, method = "BH")
  out
}

#' Wilcoxon rank-sum of one group against the reference group.
#'
#' The effect size is the Hodges-Lehmann location shift with its confidence
#' interval, which is what wilcox.test() reports with conf.int = TRUE.
compare_group_to_reference <- function(values, reference_values) {
  test <- wilcox.test(values, reference_values, conf.int = TRUE, exact = FALSE)
  list(
    location_shift = unname(test$estimate),
    ci_lower = test$conf.int[1],
    ci_upper = test$conf.int[2],
    p_value = test$p.value,
    statistic = unname(test$statistic)
  )
}

#' Pairwise Wilcoxon of every non-reference group against Neither, per gene set
#' and mark.
test_groups_against_reference <- function(long) {
  rows <- list()
  for (set_label in unname(SET_LABELS)) {
    for (mark in unname(MARK_DISPLAY)) {
      sub <- long[long$set_label == set_label & long$mark == mark, , drop = FALSE]
      reference <- sub$fold[sub$group == REFERENCE_GROUP]
      for (grp in COMPARISON_GROUPS) {
        values <- sub$fold[sub$group == grp]
        test <- compare_group_to_reference(values, reference)
        rows[[length(rows) + 1]] <- data.frame(
          set_label = set_label,
          mark = mark,
          group = grp,
          reference_group = REFERENCE_GROUP,
          n_group = length(values),
          n_reference = length(reference),
          median_group = median(values),
          median_reference = median(reference),
          median_difference = median(values) - median(reference),
          location_shift = test$location_shift,
          ci_lower = test$ci_lower,
          ci_upper = test$ci_upper,
          wilcoxon_w = test$statistic,
          p_value = test$p_value,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  out$q_value <- p.adjust(out$p_value, method = "BH")
  out
}

plot_mark_violins <- function(long, kruskal, out_dir) {
  cat("--- Figure 70_04b: mark fold change across the four groups ---\n")

  labels <- summarise_cells(long, c("set_label", "mark", "group"), "fold")
  positions <- mark_label_positions(long)
  labels <- dplyr::left_join(labels, positions[, c("mark", "label_y")],
                             by = "mark")

  kruskal_df <- kruskal
  kruskal_df$mark <- factor(kruskal_df$mark, levels = unname(MARK_DISPLAY))
  kruskal_df$set_label <- factor(kruskal_df$set_label,
                                 levels = unname(SET_LABELS))
  kruskal_df <- dplyr::left_join(kruskal_df, positions[, c("mark", "test_y")],
                                 by = "mark")
  kruskal_df$text <- sprintf("Kruskal-Wallis %s", format_p(kruskal_df$p_value))

  figure <- ggplot(long, aes(x = group, y = fold, fill = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.3) +
    geom_violin(scale = "width", alpha = 0.65, linewidth = 0.3,
                show.legend = FALSE) +
    geom_boxplot(width = 0.14, outlier.size = 0.15, outlier.alpha = 0.15,
                 show.legend = FALSE) +
    geom_text(data = labels, aes(x = group, y = label_y, label = label),
              inherit.aes = FALSE, size = 2.5, lineheight = 1.05) +
    geom_text(data = kruskal_df, aes(x = 2.5, y = test_y, label = text),
              inherit.aes = FALSE, size = 2.8, color = "grey25") +
    facet_grid(mark ~ set_label, scales = "free_y") +
    scale_fill_manual(values = GROUP_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0.16, 0.20))) +
    scale_x_discrete(labels = GROUP_AXIS_LABELS) +
    labs(
      title = "Chromatin mark change by gene-set and MeCP2 membership",
      subtitle = paste("DiffBind log2 fold change (mut / ctrl) per gene;",
                       "n and median printed above every group"),
      x = "", y = "DiffBind log2 fold change (mut / ctrl)"
    ) +
    theme_emseq(base_size = 11) +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"),
          axis.text.x = element_text(size = 8))

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "70_04b_mark_folds_by_group"),
                          width = 14, height = 15)

  summary_out <- labels[, c("set_label", "mark", "group", "n", "median", "mean")]
  colnames(summary_out) <- c("set_label", "mark", "group", "n_genes",
                             "median_fold", "mean_fold")
  write_section_table(summary_out,
                  file.path(out_dir, "70_04b_mark_fold_summary.tsv"))
  write_section_table(kruskal, file.path(out_dir, "70_04b_kruskal_wallis.tsv"))

  for (i in seq_len(nrow(kruskal))) {
    cat(sprintf("  %-18s %-10s Kruskal-Wallis chi-squared = %.1f, df = %d, %s\n",
                kruskal$set_label[i], kruskal$mark[i], kruskal$chi_squared[i],
                kruskal$df[i], format_p(kruskal$p_value[i])))
  }
  cat("\n")
  figure
}

plot_effect_forest <- function(pairwise, out_dir) {
  cat("--- Figure 70_04c: effect size against the Neither group ---\n")

  for (i in seq_len(nrow(pairwise))) {
    r <- pairwise[i, ]
    cat(sprintf("  %-18s %-10s %-14s shift = %+.4f [%+.4f, %+.4f], %s\n",
                r$set_label, r$mark, r$group, r$location_shift,
                r$ci_lower, r$ci_upper, format_p(r$p_value)))
  }

  plot_df <- pairwise
  plot_df$mark <- factor(plot_df$mark, levels = unname(MARK_DISPLAY))
  plot_df$set_label <- factor(plot_df$set_label, levels = unname(SET_LABELS))
  plot_df$group <- factor(plot_df$group, levels = rev(COMPARISON_GROUPS))
  plot_df$point_label <- sprintf("%+.3f  %s%s",
                                 plot_df$location_shift,
                                 format_p(plot_df$p_value),
                                 significance_stars(plot_df$q_value))

  figure <- ggplot(plot_df, aes(x = location_shift, y = group, color = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y",
                  width = 0.2, linewidth = 0.7, show.legend = FALSE) +
    geom_point(size = 3, show.legend = FALSE) +
    geom_text(aes(label = point_label), vjust = -1.1, size = 2.7,
              color = "grey20", show.legend = FALSE) +
    facet_grid(mark ~ set_label, scales = "free_x") +
    scale_color_manual(values = GROUP_COLORS) +
    scale_x_continuous(expand = expansion(mult = c(0.20, 0.20))) +
    scale_y_discrete(expand = expansion(add = 0.75),
                     labels = GROUP_AXIS_LABELS) +
    labs(
      title = "Chromatin change against the Neither group",
      subtitle = paste("Hodges-Lehmann location shift with 95% confidence",
                       "interval; stars mark BH q < 0.05"),
      x = "Location shift of DiffBind log2 fold change against Neither",
      y = ""
    ) +
    theme_emseq(base_size = 11) +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"),
          axis.text.y = element_text(size = 8))

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "70_04c_effect_size_forest"),
                          width = 14, height = 13)
  write_section_table(pairwise,
                  file.path(out_dir, "70_04c_wilcoxon_vs_neither.tsv"))
  cat("\n")
  figure
}

# =============================================================================
# STEP 4: K119ub DECILE ANALYSIS
# =============================================================================

#' Per-decile counts and fractions of each group, per gene set.
summarise_deciles <- function(combined) {
  df <- combined[!is.na(combined$k119ub_decile), , drop = FALSE]

  bounds <- df %>%
    dplyr::filter(set_slug == SET_SLUGS[1]) %>%
    dplyr::group_by(k119ub_decile) %>%
    dplyr::summarise(fold_min = min(k119ub_fold), fold_max = max(k119ub_fold),
                     .groups = "drop") %>%
    as.data.frame()

  summary_df <- df %>%
    dplyr::group_by(set_label, k119ub_decile) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_both = sum(is_both),
      n_in_set = sum(in_set),
      n_mecp2_up = sum(is_mecp2_up),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      frac_both = n_both / n_total,
      frac_in_set = n_in_set / n_total,
      frac_mecp2_up = n_mecp2_up / n_total
    ) %>%
    as.data.frame()

  ci <- mapply(function(x, n) binom.test(x, n)$conf.int,
               summary_df$n_both, summary_df$n_total)
  summary_df$both_ci_lower <- ci[1, ]
  summary_df$both_ci_upper <- ci[2, ]

  summary_df <- dplyr::left_join(summary_df, bounds, by = "k119ub_decile")
  list(summary = summary_df, universe = df)
}

#' Group composition of every decile, per gene set.
summarise_decile_composition <- function(decile_universe) {
  decile_universe %>%
    dplyr::count(set_label, k119ub_decile, group, .drop = FALSE) %>%
    dplyr::group_by(set_label, k119ub_decile) %>%
    dplyr::mutate(frac = n / sum(n)) %>%
    dplyr::ungroup() %>%
    as.data.frame()
}

#' Register the Fisher tests asking whether the double-overlap genes sit in the
#' extreme K119ub deciles.
register_decile_fisher <- function(decile_universe, out_dir) {
  rows <- list()
  for (slug in SET_SLUGS) {
    df <- decile_universe[decile_universe$set_slug == slug, , drop = FALSE]
    both_col <- paste0("is_", slug, "_and_mecp2_up")

    specs <- list(
      list(col = "in_k119ub_top_decile", decile = N_DECILES,
           id = sprintf("%s_both_x_k119ub_top_decile", slug),
           label = "Top decile (D10, strongest K119ub gain)"),
      list(col = "in_k119ub_bottom_decile", decile = 1L,
           id = sprintf("%s_both_x_k119ub_bottom_decile", slug),
           label = "Bottom decile (D1, strongest K119ub loss)")
    )

    for (spec in specs) {
      gene_df <- data.frame(gene_name = df$gene_name, chr = df$chr,
                            stringsAsFactors = FALSE)
      gene_df[[both_col]] <- df$is_both
      gene_df[[spec$col]] <- df$k119ub_decile == spec$decile

      ft <- register_fisher_test(
        section = SECTION_ID, test_id = spec$id,
        description = sprintf(paste("Genes in both the %s and the MeCP2-up set",
                                    "are enriched in decile %d of the gene-level",
                                    "H2AK119ub log2 fold change."),
                              tolower(SET_LABELS[[slug]]), spec$decile),
        gene_df = gene_df, row_var = both_col, col_var = spec$col,
        output_dir = out_dir
      )

      in_decile <- gene_df[[spec$col]]
      rows[[length(rows) + 1]] <- data.frame(
        set_label = SET_LABELS[[slug]],
        test_id = spec$id,
        decile = spec$decile,
        decile_label = spec$label,
        n_genes = nrow(gene_df),
        n_in_decile = sum(in_decile),
        n_both = sum(gene_df[[both_col]]),
        n_both_in_decile = sum(in_decile & gene_df[[both_col]]),
        frac_both_in_decile = sum(in_decile & gene_df[[both_col]]) /
          sum(in_decile),
        frac_both_overall = mean(gene_df[[both_col]]),
        odds_ratio = unname(ft$estimate),
        ci_lower = ft$conf.int[1],
        ci_upper = ft$conf.int[2],
        p_value = ft$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out$q_value <- p.adjust(out$p_value, method = "BH")
  out
}

plot_decile_panel <- function(deciles, composition, fisher_rows, out_dir) {
  cat("--- Figure 70_04d: K119ub decile analysis ---\n")

  summary_df <- deciles$summary

  overall <- deciles$universe %>%
    dplyr::group_by(set_label) %>%
    dplyr::summarise(frac_both = mean(is_both), .groups = "drop") %>%
    as.data.frame()

  for (i in seq_len(nrow(summary_df))) {
    r <- summary_df[i, ]
    cat(sprintf("  %-18s D%-2d n = %-6s Both = %-5s (%.1f%%)  fold %+.3f to %+.3f\n",
                as.character(r$set_label), r$k119ub_decile,
                format(r$n_total, big.mark = ","),
                format(r$n_both, big.mark = ","), 100 * r$frac_both,
                r$fold_min, r$fold_max))
  }

  summary_df$decile_label <- factor(sprintf("D%d", summary_df$k119ub_decile),
                                    levels = sprintf("D%d", seq_len(N_DECILES)))
  summary_df$bar_label <- sprintf("%.1f%%\nn = %s",
                                  100 * summary_df$frac_both,
                                  format(summary_df$n_both, big.mark = ",",
                                         trim = TRUE))

  p_fraction <- ggplot(summary_df, aes(x = decile_label, y = frac_both)) +
    geom_col(fill = GROUP_COLORS[["Both"]], width = 0.75) +
    geom_errorbar(aes(ymin = both_ci_lower, ymax = both_ci_upper), width = 0.25,
                  linewidth = 0.5) +
    geom_text(aes(label = bar_label, y = both_ci_upper), vjust = -0.3,
              size = 2.6, lineheight = 1.05) +
    geom_hline(data = overall, aes(yintercept = frac_both), linetype = "dashed",
               color = "grey35", linewidth = 0.7) +
    facet_wrap(~ set_label) +
    scale_y_continuous(labels = scales::percent_format(),
                       expand = expansion(mult = c(0, 0.22))) +
    labs(
      title = "Double-overlap genes across H2AK119ub fold-change deciles",
      subtitle = paste("D1 is the strongest K119ub loss, D10 the strongest",
                       "gain; dashed line is the overall Both fraction"),
      x = "H2AK119ub log2 fold-change decile", y = "Fraction in the Both group"
    ) +
    theme_emseq(base_size = 11) +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"))

  composition$decile_label <- factor(sprintf("D%d", composition$k119ub_decile),
                                     levels = sprintf("D%d", seq_len(N_DECILES)))

  p_composition <- ggplot(composition,
                          aes(x = decile_label, y = frac, fill = group)) +
    geom_col(width = 0.8) +
    facet_wrap(~ set_label) +
    scale_fill_manual(values = GROUP_COLORS, name = "") +
    scale_y_continuous(labels = scales::percent_format(),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = "Group composition of every H2AK119ub decile",
      x = "H2AK119ub log2 fold-change decile", y = "Fraction of the decile"
    ) +
    theme_emseq(base_size = 11) +
    theme(legend.position = "top")

  figure <- p_fraction / p_composition +
    plot_layout(heights = c(1, 0.85)) +
    plot_annotation(
      title = "H2AK119ub decile position of the double-overlap genes",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = 15))
    )

  save_multiformat_ggplot(figure, file.path(out_dir, "70_04d_k119ub_decile"),
                          width = 14, height = 12)

  write_section_table(summary_df[, c("set_label", "k119ub_decile", "fold_min",
                                 "fold_max", "n_total", "n_both", "n_in_set",
                                 "n_mecp2_up", "frac_both", "frac_in_set",
                                 "frac_mecp2_up", "both_ci_lower",
                                 "both_ci_upper")],
                  file.path(out_dir, "70_04d_k119ub_decile_summary.tsv"))
  write_section_table(composition[, c("set_label", "k119ub_decile", "group", "n",
                                  "frac")],
                  file.path(out_dir, "70_04d_k119ub_decile_composition.tsv"))
  write_section_table(fisher_rows,
                  file.path(out_dir, "70_04d_extreme_decile_fisher.tsv"))
  cat("\n")

  list(figure = figure, fraction_panel = p_fraction, summary = summary_df)
}

# =============================================================================
# STEP 5: GENE-LENGTH CONTROL
# =============================================================================

#' Bounds of every equal-count gene-length quintile.
summarise_quintile_bounds <- function(universe) {
  universe %>%
    dplyr::group_by(length_quintile) %>%
    dplyr::summarise(
      n_genes = dplyr::n(),
      length_min = min(gene_length),
      length_max = max(gene_length),
      length_median = median(gene_length),
      .groups = "drop"
    ) %>%
    as.data.frame()
}

quintile_labels <- function(bounds) {
  labels <- sprintf("Q%d\n%s-%s kb", bounds$length_quintile,
                    format(round(bounds$length_min / 1000, 1), trim = TRUE),
                    format(round(bounds$length_max / 1000, 1), trim = TRUE))
  names(labels) <- as.character(bounds$length_quintile)
  labels
}

#' Repeat the against-Neither comparison inside every gene-length quintile.
#'
#' A cell with fewer than MIN_STRATUM_SIZE genes in either arm is recorded with
#' NA test columns and a reason string instead of being tested.
test_length_quintiles <- function(long) {
  rows <- list()
  for (set_label in unname(SET_LABELS)) {
    for (mark in unname(MARK_DISPLAY)) {
      for (quintile in seq_len(N_LENGTH_QUINTILES)) {
        sub <- long[long$set_label == set_label & long$mark == mark &
                      long$length_quintile == quintile, , drop = FALSE]
        reference <- sub$fold[sub$group == REFERENCE_GROUP]

        for (grp in COMPARISON_GROUPS) {
          values <- sub$fold[sub$group == grp]
          row <- data.frame(
            set_label = set_label,
            mark = mark,
            length_quintile = quintile,
            group = grp,
            reference_group = REFERENCE_GROUP,
            n_group = length(values),
            n_reference = length(reference),
            median_group = ifelse(length(values) > 0, median(values), NA_real_),
            median_reference = ifelse(length(reference) > 0, median(reference),
                                      NA_real_),
            location_shift = NA_real_,
            ci_lower = NA_real_,
            ci_upper = NA_real_,
            p_value = NA_real_,
            untested_reason = NA_character_,
            stringsAsFactors = FALSE
          )

          if (length(values) < MIN_STRATUM_SIZE ||
              length(reference) < MIN_STRATUM_SIZE) {
            row$untested_reason <- sprintf(
              "fewer than %d genes in an arm (group n = %d, reference n = %d)",
              MIN_STRATUM_SIZE, length(values), length(reference))
          } else {
            test <- compare_group_to_reference(values, reference)
            row$location_shift <- test$location_shift
            row$ci_lower <- test$ci_lower
            row$ci_upper <- test$ci_upper
            row$p_value <- test$p_value
          }
          rows[[length(rows) + 1]] <- row
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  out$q_value <- NA_real_
  tested <- !is.na(out$p_value)
  out$q_value[tested] <- p.adjust(out$p_value[tested], method = "BH")
  out
}

plot_quintile_violins <- function(long, slug, bounds, out_dir) {
  set_label <- SET_LABELS[[slug]]
  cat(sprintf("--- Figure 70_04e: gene-length quintiles, %s ---\n", set_label))

  labels_map <- quintile_labels(bounds)
  plot_df <- long[long$set_label == set_label, , drop = FALSE]
  plot_df$quintile_label <- factor(
    labels_map[as.character(plot_df$length_quintile)],
    levels = unname(labels_map))

  cell_labels <- summarise_cells(plot_df, c("mark", "quintile_label", "group"),
                                 "fold")
  positions <- mark_label_positions(plot_df)
  cell_labels <- dplyr::left_join(cell_labels, positions[, c("mark", "label_y")],
                                  by = "mark")

  figure <- ggplot(plot_df, aes(x = group, y = fold, fill = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.3) +
    geom_violin(scale = "width", alpha = 0.65, linewidth = 0.25,
                show.legend = FALSE) +
    geom_boxplot(width = 0.14, outlier.size = 0.1, outlier.alpha = 0.12,
                 show.legend = FALSE) +
    geom_text(data = cell_labels, aes(x = group, y = label_y, label = label),
              inherit.aes = FALSE, size = 1.9, lineheight = 1.02) +
    facet_grid(mark ~ quintile_label, scales = "free_y") +
    scale_fill_manual(values = GROUP_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.22))) +
    scale_x_discrete(labels = GROUP_AXIS_LABELS) +
    labs(
      title = sprintf("Chromatin mark change inside gene-length quintiles (%s)",
                      set_label),
      subtitle = paste("Equal-count quintiles of gene length, so gene length",
                       "cannot drive a within-panel difference"),
      x = "", y = "DiffBind log2 fold change (mut / ctrl)"
    ) +
    theme_emseq(base_size = 10) +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"),
          axis.text.x = element_text(size = 6.5),
          strip.text.x = element_text(size = 8))

  save_multiformat_ggplot(
    figure,
    file.path(out_dir, sprintf("70_04e_length_quintile_violins_%s", slug)),
    width = 20, height = 15
  )
  cat("\n")
  figure
}

plot_quintile_effects <- function(quintile_stats, bounds, out_dir) {
  cat("--- Figure 70_04f: effect size inside gene-length quintiles ---\n")

  n_untested <- sum(!is.na(quintile_stats$untested_reason))
  cat(sprintf("  %d of %d quintile cells were tested; %d were too small\n",
              nrow(quintile_stats) - n_untested, nrow(quintile_stats),
              n_untested))
  if (n_untested > 0) {
    untested <- quintile_stats[!is.na(quintile_stats$untested_reason), ]
    for (i in seq_len(nrow(untested))) {
      r <- untested[i, ]
      cat(sprintf("    untested: %s / %s / Q%d / %s: %s\n",
                  r$set_label, r$mark, r$length_quintile, r$group,
                  r$untested_reason))
    }
  }

  labels_map <- quintile_labels(bounds)
  plot_df <- quintile_stats[is.na(quintile_stats$untested_reason), , drop = FALSE]
  plot_df$mark <- factor(plot_df$mark, levels = unname(MARK_DISPLAY))
  plot_df$set_label <- factor(plot_df$set_label, levels = unname(SET_LABELS))
  plot_df$group <- factor(plot_df$group, levels = COMPARISON_GROUPS)
  plot_df$quintile_label <- factor(
    labels_map[as.character(plot_df$length_quintile)],
    levels = unname(labels_map))

  figure <- ggplot(plot_df, aes(x = quintile_label, y = location_shift,
                                color = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.25,
                  linewidth = 0.6,
                  position = position_dodge(width = 0.6)) +
    geom_point(size = 2.4, position = position_dodge(width = 0.6)) +
    geom_text(aes(label = sprintf("n = %s", format(n_group, big.mark = ",",
                                                   trim = TRUE))),
              position = position_dodge(width = 0.6), vjust = -1.1, size = 2.2,
              show.legend = FALSE) +
    facet_grid(mark ~ set_label, scales = "free_y") +
    scale_color_manual(values = GROUP_COLORS, name = "") +
    scale_y_continuous(expand = expansion(mult = c(0.18, 0.22))) +
    labs(
      title = "Chromatin change against Neither inside gene-length quintiles",
      subtitle = paste("Hodges-Lehmann location shift with 95% confidence",
                       "interval; Q1 is the shortest quintile"),
      x = "Gene-length quintile",
      y = "Location shift against Neither"
    ) +
    theme_emseq(base_size = 11) +
    theme(legend.position = "top",
          plot.subtitle = element_text(size = 9, color = "grey40"),
          axis.text.x = element_text(size = 7.5))

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "70_04f_length_quintile_effects"),
                          width = 15, height = 13)
  write_section_table(quintile_stats,
                  file.path(out_dir, "70_04f_length_quintile_wilcoxon.tsv"))
  cat("\n")
  figure
}

# =============================================================================
# FIGURE 70_04g: MeCP2 AGAINST H2AK119ub PER GENE
# =============================================================================

plot_mecp2_against_k119ub <- function(combined, out_dir) {
  cat("--- Figure 70_04g: MeCP2 against H2AK119ub fold change ---\n")

  df <- combined[!is.na(combined$mecp2_fold) & !is.na(combined$k119ub_fold), ,
                 drop = FALSE]
  if (nrow(df) < 100) {
    stop("Only ", nrow(df), " genes carry both a MeCP2 and a K119ub fold ",
         "change. Too few for the per-gene scatter.")
  }

  one_set <- df[df$set_slug == SET_SLUGS[1], , drop = FALSE]
  test <- cor.test(one_set$mecp2_fold, one_set$k119ub_fold, method = "spearman",
                   exact = FALSE)
  cat(sprintf("  %s genes with both folds, Spearman rho = %.3f, %s\n",
              format(nrow(one_set), big.mark = ","), unname(test$estimate),
              format_p(test$p.value)))

  df <- df[order(df$group, decreasing = TRUE), , drop = FALSE]
  key_df <- df[df$gene_name %in% KEY_GENES, , drop = FALSE]
  cat(sprintf("  Key genes on the scatter: %d of %d\n",
              length(unique(key_df$gene_name)), length(KEY_GENES)))

  figure <- ggplot(df, aes(x = mecp2_fold, y = k119ub_fold, color = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60",
               linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60",
               linewidth = 0.3) +
    geom_point(size = 0.5, alpha = 0.4) +
    facet_wrap(~ set_label) +
    scale_color_manual(values = GROUP_COLORS, name = "",
                       guide = guide_legend(override.aes = list(size = 3,
                                                                alpha = 1))) +
    labs(
      title = "MeCP2 against H2AK119ub change, coloured by group",
      subtitle = sprintf(paste("Spearman rho = %.3f, %s, n = %s genes with",
                               "both fold changes"),
                         unname(test$estimate), format_p(test$p.value),
                         format(nrow(one_set), big.mark = ",")),
      x = "MeCP2 DiffBind log2 fold change (mut / ctrl)",
      y = "H2AK119ub DiffBind log2 fold change (mut / ctrl)"
    ) +
    theme_emseq() +
    theme(legend.position = "top",
          plot.subtitle = element_text(size = 9, color = "grey40"))

  if (nrow(key_df) > 0) {
    figure <- figure +
      geom_point(data = key_df, size = 1.8, color = "black") +
      ggrepel::geom_text_repel(data = key_df, aes(label = gene_name),
                               size = 3, min.segment.length = 0,
                               max.overlaps = Inf, color = "black",
                               show.legend = FALSE)
  }

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "70_04g_mecp2_vs_k119ub"),
                          width = 13, height = 7)

  key_out <- key_df[, c("gene_name", "chr", "set_label", "group", "gene_length",
                        "mecp2_fold", "mecp2_fdr", "k119ub_fold", "k119ub_fdr",
                        "mch_diff", "edger_fdr")]
  key_out$spearman_rho <- unname(test$estimate)
  key_out$spearman_p <- test$p.value
  write_section_table(key_out,
                  file.path(out_dir, "70_04g_key_genes.tsv"))
  cat("\n")
  figure
}

# =============================================================================
# REGISTERED FISHER TESTS
# =============================================================================

#' Fisher test of gene-set membership against MeCP2-up membership.
register_overlap_fisher <- function(set_tables, out_dir) {
  cat("--- Fisher: gene set against MeCP2-up ---\n")

  rows <- lapply(SET_SLUGS, function(slug) {
    df <- set_tables[[slug]]
    set_col <- paste0("is_", slug)

    gene_df <- data.frame(gene_name = df$gene_name, chr = df$chr,
                          stringsAsFactors = FALSE)
    gene_df[[set_col]] <- df$in_set
    gene_df$is_mecp2_up <- df$is_mecp2_up

    ft <- register_fisher_test(
      section = SECTION_ID,
      test_id = sprintf("%s_x_mecp2_up", slug),
      description = sprintf(paste("Genes of the %s are enriched among genes",
                                  "that gain MeCP2 binding in BAP1-KO."),
                            tolower(SET_LABELS[[slug]])),
      gene_df = gene_df, row_var = set_col, col_var = "is_mecp2_up",
      output_dir = out_dir
    )

    data.frame(
      set_label = SET_LABELS[[slug]],
      test_id = sprintf("%s_x_mecp2_up", slug),
      n_genes = nrow(gene_df),
      n_in_set = sum(df$in_set),
      n_mecp2_up = sum(df$is_mecp2_up),
      n_both = sum(df$is_both),
      odds_ratio = unname(ft$estimate),
      ci_lower = ft$conf.int[1],
      ci_upper = ft$conf.int[2],
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out$q_value <- p.adjust(out$p_value, method = "BH")
  write_section_table(out, file.path(out_dir, "70_04_overlap_fisher.tsv"))
  cat("\n")
  out
}

#' Fisher test of double-overlap membership against a significantly changed
#' peak, one test per gene set and mark.
register_remodeling_fisher <- function(set_tables, fdr_threshold, out_dir) {
  cat("--- Fisher: double overlap against a significantly changed peak ---\n")

  rows <- list()
  for (slug in SET_SLUGS) {
    df <- set_tables[[slug]]
    both_col <- paste0("is_", slug, "_and_mecp2_up")

    for (mark in MARK_ORDER) {
      remodeled_col <- paste0(mark, "_remodeled")
      gene_df <- data.frame(gene_name = df$gene_name, chr = df$chr,
                            stringsAsFactors = FALSE)
      gene_df[[both_col]] <- df$is_both
      gene_df[[remodeled_col]] <- df[[remodeled_col]]

      ft <- register_fisher_test(
        section = SECTION_ID,
        test_id = sprintf("%s_both_x_%s_remodeled", slug, mark),
        description = sprintf(paste("Genes in both the %s and the MeCP2-up set",
                                    "are enriched among genes whose %s peak",
                                    "changes at FDR < %.3g."),
                              tolower(SET_LABELS[[slug]]),
                              MARK_META[[mark]]$display, fdr_threshold),
        gene_df = gene_df, row_var = both_col, col_var = remodeled_col,
        output_dir = out_dir
      )

      tested <- !is.na(df[[remodeled_col]])
      rows[[length(rows) + 1]] <- data.frame(
        set_label = SET_LABELS[[slug]],
        mark = MARK_META[[mark]]$display,
        test_id = sprintf("%s_both_x_%s_remodeled", slug, mark),
        n_genes_with_peak = sum(tested),
        n_both = sum(df$is_both & tested),
        n_remodeled = sum(df[[remodeled_col]], na.rm = TRUE),
        n_both_remodeled = sum(df$is_both & df[[remodeled_col]], na.rm = TRUE),
        frac_remodeled_both = sum(df$is_both & df[[remodeled_col]],
                                  na.rm = TRUE) / sum(df$is_both & tested),
        frac_remodeled_other = sum(!df$is_both & df[[remodeled_col]],
                                   na.rm = TRUE) / sum(!df$is_both & tested),
        odds_ratio = unname(ft$estimate),
        ci_lower = ft$conf.int[1],
        ci_upper = ft$conf.int[2],
        p_value = ft$p.value,
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, rows)
  out$q_value <- p.adjust(out$p_value, method = "BH")
  write_section_table(out, file.path(out_dir, "70_04_remodeling_fisher.tsv"))
  cat("\n")
  out
}

# =============================================================================
# COMPOSITE AND PER-GENE TABLE
# =============================================================================

plot_composite <- function(panels, sizes, out_dir) {
  cat("--- Figure 70_04h: composite ---\n")

  n_both <- sizes$n_genes[sizes$group == "Both"]
  set_names <- as.character(sizes$set_label[sizes$group == "Both"])

  figure <- (panels$sizes | panels$decile_fraction) /
    panels$violins / panels$forest +
    plot_layout(heights = c(1, 2.2, 1.8)) +
    plot_annotation(
      title = paste("Section 70_04: chromatin remodeling at gene-set and",
                    "MeCP2-up genes"),
      subtitle = sprintf("Double overlap (gene set and MeCP2-up): %s",
                         paste(sprintf("%s %s genes", set_names,
                                       format(n_both, big.mark = ",",
                                              trim = TRUE)),
                               collapse = " | ")),
      theme = theme(
        plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey35")
      )
    )

  save_multiformat_ggplot(figure, file.path(out_dir, "70_04h_composite"),
                          width = 20, height = 24)
  cat("\n")
  invisible(figure)
}

write_per_gene_table <- function(universe, set_tables, out_dir) {
  out <- universe[, c("gene_name", "gene_id", "chr", "start", "end",
                      "gene_length", "length_quintile",
                      "mch_diff", "edger_logFC", "edger_fdr", "mch_sig",
                      "mch_direction",
                      "mecp2_fold", "mecp2_fdr", "mecp2_n_peaks",
                      "is_mecp2_up", "is_neuronal", "is_synapse",
                      "k119ub_decile",
                      paste0(MARK_ORDER, "_fold"),
                      paste0(MARK_ORDER, "_fdr"),
                      paste0(MARK_ORDER, "_remodeled"))]

  for (slug in SET_SLUGS) {
    out[[paste0(slug, "_group")]] <- as.character(set_tables[[slug]]$group)
  }

  out <- out[order(out$gene_name), ]
  write_section_table(out, file.path(out_dir, "70_04_per_gene_groups.tsv"))
  invisible(out)
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opt <- parse_section_args()
  out_dir <- opt$output_dir
  fdr_threshold <- opt$fdr_threshold

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("===========================================================================\n")
  cat("  SECTION 70_04: SYNAPSE AND NEURONAL CHROMATIN AT MeCP2-UP GENES\n")
  cat("===========================================================================\n")
  cat("  Output directory: ", out_dir, "\n", sep = "")
  cat("  FDR threshold:    ", fdr_threshold, "\n", sep = "")
  cat("  Overlap:          double (gene set AND MeCP2-up)\n\n")

  cat("--- Step 1: inputs ---\n")
  master <- load_gene_level_marks()
  neuronal_genes <- load_neuronal_gene_set()
  synapse_genes <- load_synapse_gene_set()
  cat("\n")
  mecp2_gene_table <- derive_mecp2_gene_table(fdr_threshold)

  universe <- build_universe(master, mecp2_gene_table, neuronal_genes,
                             synapse_genes, fdr_threshold)

  set_tables <- lapply(SET_SLUGS, function(slug) {
    build_set_table(universe, slug,
                    if (slug == "neuronal") neuronal_genes else synapse_genes)
  })
  names(set_tables) <- SET_SLUGS
  combined <- do.call(rbind, set_tables)

  long <- build_mark_long(combined)
  check_group_sizes(long)

  sizes <- summarise_group_sizes(combined)
  p_sizes <- plot_group_sizes(sizes, out_dir)

  kruskal <- test_kruskal_by_mark(long)
  p_violins <- plot_mark_violins(long, kruskal, out_dir)

  pairwise <- test_groups_against_reference(long)
  p_forest <- plot_effect_forest(pairwise, out_dir)

  deciles <- summarise_deciles(combined)
  composition <- summarise_decile_composition(deciles$universe)
  decile_fisher <- register_decile_fisher(deciles$universe, out_dir)
  decile <- plot_decile_panel(deciles, composition, decile_fisher, out_dir)

  bounds <- summarise_quintile_bounds(universe)
  write_section_table(bounds, file.path(out_dir, "70_04e_length_quintile_bounds.tsv"))
  for (slug in SET_SLUGS) {
    plot_quintile_violins(long, slug, bounds, out_dir)
  }
  quintile_stats <- test_length_quintiles(long)
  plot_quintile_effects(quintile_stats, bounds, out_dir)

  plot_mecp2_against_k119ub(combined, out_dir)

  overlap_fisher <- register_overlap_fisher(set_tables, out_dir)
  remodeling_fisher <- register_remodeling_fisher(set_tables, fdr_threshold,
                                                  out_dir)

  plot_composite(list(sizes = p_sizes, decile_fraction = decile$fraction_panel,
                      violins = p_violins, forest = p_forest),
                 sizes, out_dir)

  write_per_gene_table(universe, set_tables, out_dir)

  cat("\n---------------------------------------------------------------------------\n")
  cat("  SUMMARY\n")
  cat("---------------------------------------------------------------------------\n")
  cat(sprintf("  Universe: %s genes, %s MeCP2-up\n",
              format(nrow(universe), big.mark = ","),
              format(sum(universe$is_mecp2_up), big.mark = ",")))
  for (slug in SET_SLUGS) {
    df <- set_tables[[slug]]
    cat(sprintf("  %-18s in set %s, double overlap %s\n",
                SET_LABELS[[slug]], format(sum(df$in_set), big.mark = ","),
                format(sum(df$is_both), big.mark = ",")))
  }
  cat("  Gene set against MeCP2-up:\n")
  for (i in seq_len(nrow(overlap_fisher))) {
    r <- overlap_fisher[i, ]
    cat(sprintf("    %-18s OR = %.2f [%.2f, %.2f], %s\n", r$set_label,
                r$odds_ratio, r$ci_lower, r$ci_upper, format_p(r$p_value)))
  }
  cat("  Double overlap against Neither, location shift:\n")
  both_rows <- pairwise[pairwise$group == "Both", ]
  for (i in seq_len(nrow(both_rows))) {
    r <- both_rows[i, ]
    cat(sprintf("    %-18s %-10s %+.4f [%+.4f, %+.4f], %s (q = %.3g)\n",
                r$set_label, r$mark, r$location_shift, r$ci_lower, r$ci_upper,
                format_p(r$p_value), r$q_value))
  }
  cat("  Double overlap against a significantly changed peak:\n")
  for (i in seq_len(nrow(remodeling_fisher))) {
    r <- remodeling_fisher[i, ]
    cat(sprintf("    %-18s %-10s OR = %.2f [%.2f, %.2f], %s\n", r$set_label,
                r$mark, r$odds_ratio, r$ci_lower, r$ci_upper,
                format_p(r$p_value)))
  }
  cat("  K119ub extreme deciles:\n")
  for (i in seq_len(nrow(decile_fisher))) {
    r <- decile_fisher[i, ]
    cat(sprintf("    %-18s %-42s OR = %.2f, %s\n", r$set_label, r$decile_label,
                r$odds_ratio, format_p(r$p_value)))
  }
  n_untested <- sum(!is.na(quintile_stats$untested_reason))
  cat(sprintf("  Gene-length quintiles: %d of %d cells tested, %d too small\n",
              nrow(quintile_stats) - n_untested, nrow(quintile_stats),
              n_untested))
  cat("\n=== Section 70_04 complete ===\n\n")
}

main()
