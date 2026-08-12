# scripts/sections/20_chip_integration/20_03_quadrant_scatters.R
#
# Section 20_03: pairwise effect-size scatters relating the chromatin marks and
# gene-body mCH to each other at gene level.
#
# What this tests
#   BAP1 loss moves H2AK119ub, MeCP2, H3K27ac, H3K27me3 and gene-body mCH. If
#   those changes belong to one regulatory axis, the gene-level effect sizes
#   must agree in sign. Each scatter puts two effect sizes on one panel, counts
#   the genes in each quadrant, reports the Spearman correlation, and registers
#   a gene-level Fisher test of the two signs.
#
# Five scatters
#   20_03a  H2AK119ub gene-body log2FC  against  MeCP2 fold change
#   20_03b  MeCP2 fold change           against  H3K27ac fold change   (euchromatin)
#   20_03c  MeCP2 fold change           against  H3K27me3 fold change  (heterochromatin)
#   20_03d  H2AK119ub gene-body log2FC  against  mCH difference
#   20_03e  MeCP2 fold change           against  mCH difference
#   20_03f  composite of all five
#
# Reads
#   HANDOFF_PATHS$gene_level_all_marks   gene-level mark table written by 20_02
#   DIFFBIND_PATHS$k119ub_gene_signal    H2AK119ub bigwig signal per gene body
#   mecp2_diffbind                       shared config, aggregated to genes here
#                                        with the nearest-TSS rule of 20_01
#
# Writes into --output-dir (default results/sections/20_chip_integration/)
#   20_03a .. 20_03f            figures in every output format
#   20_03_scatter_statistics.tsv       one row per pair
#   20_03_quadrant_counts.tsv          one row per pair per quadrant
#   20_03_quadrant_master.tsv          one row per gene, coordinates and
#                                      quadrant for each of the five pairs
#   20_03_chromatin_group_counts.tsv   genes per gene-body state, promoter state
#                                      and chromatin group
#   20_03_key_gene_coordinates.tsv     the labelled genes on every pair
#   20_03_fisher_summary.tsv           the registered Fisher tests
#   fisher_tables/                     gene tables behind those tests
#
# Adapted from Biomodal section 59 (log2 quadrant scatters). That section drew
# seven scatters, two of which needed the second methylation modality of the
# dual-modality assay. This pipeline measures non-CG methylation as one
# modality, so those two are dropped and five remain.

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)

# =============================================================================
# SECTION CONSTANTS
# =============================================================================

SECTION_ID <- "20_03"

# Axes are clipped symmetrically at this quantile of the absolute value, so a
# few extreme genes cannot flatten the cloud that carries the correlation.
AXIS_CLIP_PERCENTILE <- 0.995

# A scatter with fewer genes than this stops the script rather than producing a
# correlation the reader cannot rely on.
MIN_GENES_PER_SCATTER <- 50

QUADRANT_LEVELS <- c("Q1", "Q2", "Q3", "Q4")

# Geometry and meaning of each quadrant under assign_quadrant() in
# _shared_config.R: Q1 both up, Q2 x down and y up, Q3 both down, Q4 x up and
# y down. A value of exactly zero counts as "down" under that rule.
# x_side / y_side name the panel corner the label sits in; hjust / vjust keep
# the text inside the panel.
QUADRANT_GEOMETRY <- data.frame(
  quadrant    = QUADRANT_LEVELS,
  x_direction = c("up", "down", "down", "up"),
  y_direction = c("up", "up", "down", "down"),
  x_side      = c("right", "left", "left", "right"),
  y_side      = c("top", "top", "bottom", "bottom"),
  hjust       = c(1, 0, 0, 1),
  vjust       = c(1, 1, 0, 0),
  stringsAsFactors = FALSE
)

# Chromatin groups for the two restricted scatters, defined on body_state from
# classify_body_state() in _shared_config.R. The two marks in these scatters are
# gene-body marks, so the gene-body state decides the restriction.
#
#   Euchromatin      body_state is Active or Enhancer_Marked
#                    H3K27ac, or H3K4me1 without H3K27ac, and no H3K27me3
#   Heterochromatin  body_state is Polycomb or Mixed
#                    H3K27me3, with or without H3K27ac
#
# Unmarked gene bodies join neither group and are absent from both restricted
# scatters.
EUCHROMATIN_BODY_STATES <- c("Active", "Enhancer_Marked")
HETEROCHROMATIN_BODY_STATES <- c("Polycomb", "Mixed")

CHROMATIN_GROUP_LEVELS <- c("Euchromatin", "Heterochromatin", "Neither")

# The 27 columns section 20_02 writes into gene_level_all_marks.tsv. The script
# indexes these by name.
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

# Columns the H2AK119ub bigwig signal table must carry.
REQUIRED_K119UB_SIGNAL_COLUMNS <- c(
  "entrez_id", "symbol", "start", "end",
  "gb_ctrl_signal", "gb_mut_signal", "gb_log2fc", "gb_signal_class"
)

MCH_SIGNIFICANCE_LEVELS <- names(COLORS$significant)

# One entry per scatter. The master table, the figures, the statistics table
# and the Fisher tests all read this list, so a pair is defined once.
SCATTER_SPECS <- list(
  list(
    key = "a", pair_id = "k119ub_vs_mecp2",
    file_stem = "20_03a_k119ub_vs_mecp2",
    x_col = "k119ub_gb_log2fc", y_col = "mecp2_fold",
    x_short = "K119ub", y_short = "MeCP2",
    x_lab = "H2AK119ub gene-body log2FC (mutant / control)",
    y_lab = "MeCP2 log2 fold change (mutant / control)",
    title = "H2AK119ub against MeCP2 at gene bodies",
    subset_id = "k119ub_quantifiable",
    subset_label = "genes with quantifiable H2AK119ub bigwig signal",
    color_col = NULL,
    fisher_test_id = "k119ub_up_vs_mecp2_up",
    fisher_description = paste("Among genes with quantifiable H2AK119ub signal",
                               "and a MeCP2 peak, does an H2AK119ub gain occur",
                               "with a MeCP2 gain?"),
    width = 9.5, height = 8.5
  ),
  list(
    key = "b", pair_id = "mecp2_vs_k27ac",
    file_stem = "20_03b_mecp2_vs_k27ac_euchromatin",
    x_col = "mecp2_fold", y_col = "k27ac_fold",
    x_short = "MeCP2", y_short = "K27ac",
    x_lab = "MeCP2 log2 fold change (mutant / control)",
    y_lab = "H3K27ac log2 fold change (mutant / control)",
    title = "MeCP2 against H3K27ac in euchromatin",
    subset_id = "euchromatin",
    subset_label = "euchromatic genes (gene-body state Active or Enhancer_Marked)",
    color_col = NULL,
    fisher_test_id = "mecp2_up_vs_k27ac_up_euchromatin",
    fisher_description = paste("Among euchromatic genes with both a MeCP2 and",
                               "an H3K27ac peak, does a MeCP2 gain occur with",
                               "an H3K27ac gain?"),
    width = 9.5, height = 8.5
  ),
  list(
    key = "c", pair_id = "mecp2_vs_k27me3",
    file_stem = "20_03c_mecp2_vs_k27me3_heterochromatin",
    x_col = "mecp2_fold", y_col = "k27me3_fold",
    x_short = "MeCP2", y_short = "K27me3",
    x_lab = "MeCP2 log2 fold change (mutant / control)",
    y_lab = "H3K27me3 log2 fold change (mutant / control)",
    title = "MeCP2 against H3K27me3 in heterochromatin",
    subset_id = "heterochromatin",
    subset_label = "heterochromatic genes (gene-body state Polycomb or Mixed)",
    color_col = NULL,
    fisher_test_id = "mecp2_up_vs_k27me3_up_heterochromatin",
    fisher_description = paste("Among heterochromatic genes with both a MeCP2",
                               "and an H3K27me3 peak, does a MeCP2 gain occur",
                               "with an H3K27me3 gain?"),
    width = 9.5, height = 8.5
  ),
  list(
    key = "d", pair_id = "k119ub_vs_mch",
    file_stem = "20_03d_k119ub_vs_mch",
    x_col = "k119ub_gb_log2fc", y_col = "mch_diff",
    x_short = "K119ub", y_short = "mCH",
    x_lab = "H2AK119ub gene-body log2FC (mutant / control)",
    y_lab = "mCH difference (mutant - control)",
    title = "H2AK119ub against gene-body mCH change",
    subset_id = "k119ub_quantifiable",
    subset_label = "genes with quantifiable H2AK119ub bigwig signal",
    color_col = "mch_significance",
    fisher_test_id = "k119ub_up_vs_mch_hyper",
    fisher_description = paste("Among genes with quantifiable H2AK119ub signal,",
                               "does an H2AK119ub gain occur with a gene-body",
                               "mCH gain?"),
    width = 9.5, height = 8.5
  ),
  list(
    key = "e", pair_id = "mecp2_vs_mch",
    file_stem = "20_03e_mecp2_vs_mch",
    x_col = "mecp2_fold", y_col = "mch_diff",
    x_short = "MeCP2", y_short = "mCH",
    x_lab = "MeCP2 log2 fold change (mutant / control)",
    y_lab = "mCH difference (mutant - control)",
    title = "MeCP2 against gene-body mCH change",
    subset_id = "all_genes",
    subset_label = "every gene with a MeCP2 peak",
    color_col = "mch_significance",
    fisher_test_id = "mecp2_up_vs_mch_hyper",
    fisher_description = paste("Among genes with a MeCP2 peak, does a MeCP2",
                               "gain occur with a gene-body mCH gain?"),
    width = 9.5, height = 8.5
  )
)

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
                help = "FDR cutoff for MeCP2 peak significance [default: %default]"),
    make_option("--axis-clip", dest = "axis_clip", type = "double",
                default = AXIS_CLIP_PERCENTILE,
                help = "Quantile of |value| that sets the symmetric axis limits [default: %default]")
  )
  opt <- parse_args(OptionParser(option_list = option_list))

  if (opt$fdr_threshold <= 0 || opt$fdr_threshold >= 1) {
    stop("--fdr-threshold must be between 0 and 1, got ", opt$fdr_threshold)
  }
  if (opt$axis_clip <= 0.5 || opt$axis_clip > 1) {
    stop("--axis-clip must be above 0.5 and at most 1, got ", opt$axis_clip)
  }
  opt
}

# =============================================================================
# OUTPUT HELPERS
# =============================================================================

count_label <- function(n) format(n, big.mark = ",", trim = TRUE)

# =============================================================================
# INPUT LOADING
# =============================================================================

#' Read the gene-level mark table written by section 20_02.
#'
#' Stops when the file is absent or when a column the scatters index by name is
#' missing.
#'
#' @param path HANDOFF_PATHS$gene_level_all_marks.
#' @return data.frame, one row per row of mch_results
load_gene_level_all_marks <- function(path) {
  if (!file.exists(path)) {
    stop("Gene-level mark table not found: ", path,
         "\nRun section 20_02 (20_02_multi_mark_diffbind.R) first. ",
         "It writes gene_level_all_marks.tsv.")
  }

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")
  missing <- setdiff(REQUIRED_HANDOFF_COLUMNS, colnames(df))
  if (length(missing) > 0) {
    stop("gene_level_all_marks.tsv is missing columns: ",
         paste(missing, collapse = ", "),
         "\nSection 20_02 writes the column contract this section reads.")
  }

  if (!is.logical(df$mch_sig)) {
    stop("gene_level_all_marks.tsv column mch_sig is ", class(df$mch_sig)[1],
         ", not logical. Section 20_02 must write TRUE or FALSE.")
  }

  cat(sprintf("  Gene-level mark table: %s rows\n", count_label(nrow(df))))
  cat(sprintf("    mCH significant: %s | H3K27ac fold: %s | H3K27me3 fold: %s\n",
              count_label(sum(df$mch_sig)),
              count_label(sum(!is.na(df$k27ac_fold))),
              count_label(sum(!is.na(df$k27me3_fold)))))
  df
}

#' Keep one row per gene name.
#'
#' Some gene names carry more than one ENSMUSG identifier. The row with the
#' largest absolute edgeR log fold change is kept for each name, so the
#' gene-keyed joins below cannot change the row count.
#'
#' @param df Gene-level table with gene_name and edger_logFC.
#' @return data.frame with unique gene_name, in the input order
deduplicate_by_gene_name <- function(df) {
  ord <- order(df$gene_name, -abs(df$edger_logFC))
  keep <- sort(ord[!duplicated(df$gene_name[ord])])
  out <- df[keep, , drop = FALSE]
  cat(sprintf("  Collapsed %s rows to %s gene names (%s rows dropped)\n",
              count_label(nrow(df)), count_label(nrow(out)),
              count_label(nrow(df) - nrow(out))))
  out
}

#' Read the H2AK119ub bigwig signal table and key it by gene symbol.
#'
#' The file holds the absolute mean H2AK119ub signal over each gene body in
#' control and mutant, the log2 ratio of the two, and a signal class naming the
#' genes where both conditions clear the pseudocount. It is keyed on Entrez
#' identifier, so symbols that appear twice are collapsed to the widest gene
#' body.
#'
#' @param path DIFFBIND_PATHS$k119ub_gene_signal.
#' @return data.frame with gene_name and four k119ub_ columns
load_k119ub_gene_signal <- function(path) {
  if (!file.exists(path)) {
    stop("H2AK119ub gene-body signal table not found: ", path,
         "\nIt holds the per-gene H2AK119ub bigwig signal and must exist ",
         "before section 20_03 runs.")
  }

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")
  missing <- setdiff(REQUIRED_K119UB_SIGNAL_COLUMNS, colnames(df))
  if (length(missing) > 0) {
    stop("H2AK119ub gene signal table is missing columns: ",
         paste(missing, collapse = ", "), " (", path, ")")
  }

  n_all <- nrow(df)
  df <- df[!is.na(df$symbol) & nzchar(df$symbol), , drop = FALSE]

  df$gene_span <- df$end - df$start
  ord <- order(df$symbol, -df$gene_span, df$entrez_id)
  keep <- sort(ord[!duplicated(df$symbol[ord])])
  n_collapsed <- nrow(df) - length(keep)
  df <- df[keep, , drop = FALSE]

  out <- data.frame(
    gene_name              = df$symbol,
    k119ub_gb_ctrl_signal  = as.numeric(df$gb_ctrl_signal),
    k119ub_gb_mut_signal   = as.numeric(df$gb_mut_signal),
    k119ub_gb_log2fc       = as.numeric(df$gb_log2fc),
    k119ub_gb_signal_class = as.character(df$gb_signal_class),
    stringsAsFactors = FALSE
  )

  n_quantifiable <- sum(out$k119ub_gb_signal_class == "quantifiable")
  cat(sprintf("  H2AK119ub signal: %s rows in, %s gene symbols kept (%s collapsed)\n",
              count_label(n_all), count_label(nrow(out)),
              count_label(n_collapsed)))
  cat(sprintf("    quantifiable in both conditions: %s (%.1f%%)\n",
              count_label(n_quantifiable), 100 * n_quantifiable / nrow(out)))
  out
}

#' Collapse the MeCP2 DiffBind peakset to one fold change per gene.
#'
#' Uses the same route as section 20_01: ChIPseeker annotation of the DiffBind
#' peaks, then the nearest-TSS collapse rule.
#'
#' @param diffbind mecp2_diffbind from the shared config.
#' @param fdr_threshold FDR cutoff passed to the collapse rule.
#' @return data.frame with gene_name and the mecp2_ summary columns
aggregate_mecp2_to_genes <- function(diffbind, fdr_threshold) {
  annotated <- annotate_peaks_to_genes(diffbind, "MeCP2")

  usable <- !is.na(annotated$SYMBOL) & nzchar(annotated$SYMBOL) &
    !is.na(annotated$Fold) & !is.na(annotated$FDR)
  cat(sprintf("    %s of %s annotated MeCP2 peaks carry a gene symbol\n",
              count_label(sum(usable)), count_label(nrow(annotated))))
  annotated <- annotated[usable, , drop = FALSE]

  out <- aggregate_diffbind_by_gene(annotated, method = "nearest_tss",
                                    fdr_threshold = fdr_threshold,
                                    prefix = "mecp2")
  if (anyDuplicated(out$gene_name) > 0) {
    stop("aggregate_diffbind_by_gene() returned duplicate genes for MeCP2.")
  }
  cat(sprintf("    %s genes carry a MeCP2 peak\n", count_label(nrow(out))))
  out
}

# =============================================================================
# MASTER TABLE
# =============================================================================

#' Assign each gene to Euchromatin, Heterochromatin, or Neither.
#'
#' Reads body_state, not promoter_state: both restricted scatters plot gene-body
#' marks. A gene whose body_state is missing joins Neither.
#'
#' @param body_state Character vector of gene-body chromatin states.
#' @return factor with levels CHROMATIN_GROUP_LEVELS
assign_chromatin_group <- function(body_state) {
  group <- rep("Neither", length(body_state))
  group[body_state %in% EUCHROMATIN_BODY_STATES] <- "Euchromatin"
  group[body_state %in% HETEROCHROMATIN_BODY_STATES] <- "Heterochromatin"
  factor(group, levels = CHROMATIN_GROUP_LEVELS)
}

#' Join the H2AK119ub signal and MeCP2 fold onto the gene-level mark table.
#'
#' @param marks Deduplicated gene-level mark table.
#' @param k119ub_signal Table from load_k119ub_gene_signal().
#' @param mecp2_genes Table from aggregate_mecp2_to_genes().
#' @return data.frame with one row per gene name
build_master_table <- function(marks, k119ub_signal, mecp2_genes) {
  master <- dplyr::left_join(marks, k119ub_signal, by = "gene_name")
  master <- dplyr::left_join(
    master,
    mecp2_genes[, c("gene_name", "mecp2_fold", "mecp2_fdr", "mecp2_min_fdr",
                    "mecp2_n_peaks", "mecp2_n_sig", "mecp2_has_sig")],
    by = "gene_name"
  )

  if (nrow(master) != nrow(marks)) {
    stop("The joins changed the row count: ", nrow(marks), " genes in, ",
         nrow(master), " out. A joined table has duplicate gene names.")
  }

  master$chromatin_group <- assign_chromatin_group(master$body_state)
  master$mch_significance <- factor(
    ifelse(master$mch_sig, "Significant", "Not Significant"),
    levels = MCH_SIGNIFICANCE_LEVELS
  )

  cat(sprintf("  Master table: %s genes\n", count_label(nrow(master))))
  cat(sprintf("    H2AK119ub gene-body log2FC: %s | quantifiable: %s\n",
              count_label(sum(!is.na(master$k119ub_gb_log2fc))),
              count_label(sum(!is.na(master$k119ub_gb_signal_class) &
                                master$k119ub_gb_signal_class == "quantifiable"))))
  cat(sprintf("    MeCP2 fold: %s | H3K27ac fold: %s | H3K27me3 fold: %s\n",
              count_label(sum(!is.na(master$mecp2_fold))),
              count_label(sum(!is.na(master$k27ac_fold))),
              count_label(sum(!is.na(master$k27me3_fold)))))
  # Both distributions are printed, so a restriction that lands on zero genes is
  # visible here rather than only in the MIN_GENES_PER_SCATTER error below.
  for (state in BODY_STATE_ORDER) {
    cat(sprintf("    body state      %-16s %s genes\n", state,
                count_label(sum(master$body_state == state, na.rm = TRUE))))
  }
  for (group in CHROMATIN_GROUP_LEVELS) {
    cat(sprintf("    chromatin group %-16s %s genes\n", group,
                count_label(sum(master$chromatin_group == group, na.rm = TRUE))))
  }
  master
}

#' Count genes per gene-body chromatin state and per chromatin group.
write_chromatin_group_counts <- function(master, out_dir) {
  counts <- master %>%
    dplyr::count(body_state, promoter_state, chromatin_group, name = "n_genes") %>%
    dplyr::arrange(chromatin_group, dplyr::desc(n_genes)) %>%
    as.data.frame()
  counts$percent_of_genes <- 100 * counts$n_genes / nrow(master)
  write_section_table(counts, file.path(out_dir, "20_03_chromatin_group_counts.tsv"))
  invisible(counts)
}

#' Restrict the master table to the genes one scatter plots.
#'
#' @param master Master table from build_master_table().
#' @param subset_id Identifier named by a SCATTER_SPECS entry.
#' @return data.frame
apply_scatter_subset <- function(master, subset_id) {
  if (subset_id == "all_genes") return(master)

  if (subset_id == "k119ub_quantifiable") {
    keep <- !is.na(master$k119ub_gb_signal_class) &
      master$k119ub_gb_signal_class == "quantifiable"
    return(master[keep, , drop = FALSE])
  }
  if (subset_id == "euchromatin") {
    keep <- !is.na(master$chromatin_group) &
      master$chromatin_group == "Euchromatin"
    return(master[keep, , drop = FALSE])
  }
  if (subset_id == "heterochromatin") {
    keep <- !is.na(master$chromatin_group) &
      master$chromatin_group == "Heterochromatin"
    return(master[keep, , drop = FALSE])
  }
  stop("apply_scatter_subset(): unknown subset id '", subset_id, "'")
}

# =============================================================================
# SCATTER HELPERS
# =============================================================================

#' Symmetric axis limits from a quantile of the absolute values.
#'
#' @param values Numeric vector.
#' @param percentile Quantile of |values| that sets the limit.
#' @return numeric vector of length 2, centred on zero
clip_symmetric <- function(values, percentile = AXIS_CLIP_PERCENTILE) {
  limit <- unname(quantile(abs(values), percentile, na.rm = TRUE))
  if (!is.finite(limit) || limit <= 0) {
    stop("clip_symmetric(): the ", percentile, " quantile of |values| is ",
         limit, ". The axis cannot be scaled.")
  }
  c(-limit, limit)
}

#' Count genes per quadrant and place one label in each panel corner.
#'
#' @param quadrants Character vector from assign_quadrant().
#' @param x_range Two-element x limits.
#' @param y_range Two-element y limits.
#' @param x_short Short name of the x variable, used in the label text.
#' @param y_short Short name of the y variable.
#' @return data.frame with one row per quadrant: counts, text, and position
quadrant_corner_labels <- function(quadrants, x_range, y_range,
                                   x_short, y_short) {
  counts <- table(factor(quadrants, levels = QUADRANT_LEVELS))
  total <- sum(counts)
  if (total == 0) stop("quadrant_corner_labels(): no genes to count.")

  pad_x <- diff(x_range) * 0.04
  pad_y <- diff(y_range) * 0.04

  out <- QUADRANT_GEOMETRY
  out$n_genes <- as.integer(counts[out$quadrant])
  out$percent <- 100 * out$n_genes / total
  out$x <- ifelse(out$x_side == "right", x_range[2] - pad_x, x_range[1] + pad_x)
  out$y <- ifelse(out$y_side == "top", y_range[2] - pad_y, y_range[1] + pad_y)
  out$meaning <- sprintf("%s %s, %s %s", x_short, out$x_direction,
                         y_short, out$y_direction)
  out$label <- sprintf("%s: %s\nn = %s (%.1f%%)", out$quadrant, out$meaning,
                       count_label(out$n_genes), out$percent)
  out
}

#' Format a Spearman correlation for a plot subtitle.
format_rho_label <- function(rho, p_value) {
  if (p_value < 2.2e-16) return(sprintf("rho = %.3f, p < 2.2e-16", rho))
  sprintf("rho = %.3f, p = %.2e", rho, p_value)
}

#' Build one quadrant scatter.
#'
#' Drops genes missing either coordinate, clips both axes symmetrically,
#' assigns quadrants with assign_quadrant() from the shared config, writes the
#' four corner counts onto the panel, and labels KEY_GENES with ggrepel.
#'
#' @param df Genes to plot, already restricted to the scatter's subset.
#' @param spec One entry of SCATTER_SPECS.
#' @param axis_clip Quantile of |value| that sets the axis limits.
#' @return list with plot, stats, corners, plotted genes, and key gene rows
build_quadrant_scatter <- function(df, spec, axis_clip) {
  x_col <- spec$x_col
  y_col <- spec$y_col

  for (col in c(x_col, y_col)) {
    if (!col %in% colnames(df)) {
      stop("build_quadrant_scatter(): pair ", spec$pair_id,
           " needs column ", col, ", which the master table does not have.")
    }
  }

  plotted <- df[!is.na(df[[x_col]]) & !is.na(df[[y_col]]), , drop = FALSE]
  if (nrow(plotted) < MIN_GENES_PER_SCATTER) {
    stop("Pair ", spec$pair_id, " has ", nrow(plotted), " genes with both ",
         x_col, " and ", y_col, " (subset: ", spec$subset_id,
         "). The section needs at least ", MIN_GENES_PER_SCATTER, ".")
  }

  x_vals <- plotted[[x_col]]
  y_vals <- plotted[[y_col]]
  plotted$quadrant <- assign_quadrant(x_vals, y_vals)

  spearman <- cor.test(x_vals, y_vals, method = "spearman", exact = FALSE)
  pearson <- cor.test(x_vals, y_vals, method = "pearson")

  x_lim <- clip_symmetric(x_vals, axis_clip)
  y_lim <- clip_symmetric(y_vals, axis_clip)
  inside_panel <- x_vals >= x_lim[1] & x_vals <= x_lim[2] &
    y_vals >= y_lim[1] & y_vals <= y_lim[2]
  n_outside <- sum(!inside_panel)

  corners <- quadrant_corner_labels(plotted$quadrant, x_lim, y_lim,
                                    spec$x_short, spec$y_short)

  # Grey points first, coloured points on top of them.
  if (!is.null(spec$color_col)) {
    plotted <- plotted[order(plotted[[spec$color_col]], decreasing = TRUE), ,
                       drop = FALSE]
  }

  key_df <- plotted[plotted$gene_name %in% KEY_GENES, , drop = FALSE]
  key_df <- key_df[key_df[[x_col]] >= x_lim[1] & key_df[[x_col]] <= x_lim[2] &
                     key_df[[y_col]] >= y_lim[1] & key_df[[y_col]] <= y_lim[2], ,
                   drop = FALSE]

  point_layer <- if (is.null(spec$color_col)) {
    geom_point(alpha = 0.30, size = 0.7, color = "#404040")
  } else {
    geom_point(aes(color = .data[[spec$color_col]]), alpha = 0.35, size = 0.8)
  }

  plot <- ggplot(plotted, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.35) +
    point_layer +
    geom_text(data = corners,
              aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
              inherit.aes = FALSE, size = 2.8, color = "grey25",
              lineheight = 1.05) +
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    labs(
      title = spec$title,
      subtitle = sprintf("Spearman %s | n = %s %s",
                         format_rho_label(unname(spearman$estimate),
                                          spearman$p.value),
                         count_label(nrow(plotted)), spec$subset_label),
      caption = sprintf(paste("Axes clipped at the %.1fth percentile of the",
                              "absolute value; %s genes fall outside the panel",
                              "and still count in the quadrant totals."),
                        100 * axis_clip, count_label(n_outside)),
      x = spec$x_lab, y = spec$y_lab
    ) +
    theme_emseq() +
    theme(plot.caption = element_text(size = 8, color = "grey40", hjust = 0))

  if (!is.null(spec$color_col)) {
    plot <- plot +
      scale_color_manual(values = COLORS$significant, name = "mCH",
                         drop = FALSE) +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
      theme(legend.position = "top")
  }

  if (nrow(key_df) > 0) {
    plot <- plot +
      geom_point(data = key_df, aes(x = .data[[x_col]], y = .data[[y_col]]),
                 inherit.aes = FALSE, size = 1.7, color = "black") +
      ggrepel::geom_text_repel(
        data = key_df,
        aes(x = .data[[x_col]], y = .data[[y_col]], label = gene_name),
        inherit.aes = FALSE, size = 2.7, fontface = "italic", color = "black",
        segment.color = "grey50", segment.size = 0.3,
        min.segment.length = 0, max.overlaps = Inf
      )
  }

  quadrant_n <- setNames(corners$n_genes, corners$quadrant)
  stats <- data.frame(
    pair_id = spec$pair_id,
    figure = spec$file_stem,
    x_variable = x_col,
    y_variable = y_col,
    subset_id = spec$subset_id,
    subset_label = spec$subset_label,
    n_genes = nrow(plotted),
    spearman_rho = unname(spearman$estimate),
    spearman_p = spearman$p.value,
    pearson_r = unname(pearson$estimate),
    pearson_p = pearson$p.value,
    x_clip_low = x_lim[1], x_clip_high = x_lim[2],
    y_clip_low = y_lim[1], y_clip_high = y_lim[2],
    axis_clip_percentile = axis_clip,
    n_outside_panel = n_outside,
    n_q1 = quadrant_n[["Q1"]], n_q2 = quadrant_n[["Q2"]],
    n_q3 = quadrant_n[["Q3"]], n_q4 = quadrant_n[["Q4"]],
    n_key_genes_labelled = nrow(key_df),
    stringsAsFactors = FALSE
  )

  list(plot = plot, stats = stats, corners = corners, plotted = plotted,
       key_genes = key_df)
}

# =============================================================================
# FIGURES AND TABLES
# =============================================================================

#' Build every scatter, save each figure, and report its correlation.
run_all_scatters <- function(master, axis_clip, out_dir) {
  results <- list()

  for (spec in SCATTER_SPECS) {
    cat(sprintf("--- Figure %s: %s ---\n", spec$file_stem, spec$title))
    subset_df <- apply_scatter_subset(master, spec$subset_id)
    cat(sprintf("  Subset '%s': %s genes\n", spec$subset_id,
                count_label(nrow(subset_df))))

    result <- build_quadrant_scatter(subset_df, spec, axis_clip)
    save_multiformat_ggplot(result$plot, file.path(out_dir, spec$file_stem),
                            width = spec$width, height = spec$height)

    cat(sprintf("  n = %s | Spearman rho = %+.3f (p = %.3g) | Pearson r = %+.3f\n",
                count_label(result$stats$n_genes), result$stats$spearman_rho,
                result$stats$spearman_p, result$stats$pearson_r))
    for (i in seq_len(nrow(result$corners))) {
      cat(sprintf("    %-3s %-28s %s genes (%.1f%%)\n",
                  result$corners$quadrant[i], result$corners$meaning[i],
                  count_label(result$corners$n_genes[i]),
                  result$corners$percent[i]))
    }
    cat("\n")

    result$spec <- spec
    results[[spec$key]] <- result
  }
  results
}

#' Lay all five scatters into one composite panel.
save_composite_figure <- function(results, out_dir) {
  cat("--- Figure 20_03f: composite of all five scatters ---\n")

  composite <- (results$a$plot + results$b$plot) /
    (results$c$plot + results$d$plot) /
    (results$e$plot + patchwork::plot_spacer()) +
    plot_annotation(
      title = "Gene-level effect sizes: H2AK119ub, MeCP2, H3K27ac, H3K27me3 and mCH",
      subtitle = paste("Every panel counts genes per quadrant and reports the",
                       "Spearman correlation over the genes it plots"),
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 17),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(composite, file.path(out_dir, "20_03f_quadrant_composite"),
                          width = 19, height = 25)
  cat("\n")
  invisible(composite)
}

#' Write the per-pair statistics table.
write_scatter_statistics <- function(results, out_dir) {
  stats <- do.call(rbind, lapply(results, `[[`, "stats"))
  rownames(stats) <- NULL
  write_section_table(stats, file.path(out_dir, "20_03_scatter_statistics.tsv"))
  invisible(stats)
}

#' Write one row per pair per quadrant, matching the four corner labels.
write_quadrant_counts <- function(results, out_dir) {
  rows <- lapply(results, function(result) {
    corners <- result$corners
    data.frame(
      pair_id = result$spec$pair_id,
      figure = result$spec$file_stem,
      subset_id = result$spec$subset_id,
      quadrant = corners$quadrant,
      x_variable = result$spec$x_col,
      y_variable = result$spec$y_col,
      x_direction = corners$x_direction,
      y_direction = corners$y_direction,
      meaning = corners$meaning,
      n_genes = corners$n_genes,
      percent = corners$percent,
      stringsAsFactors = FALSE
    )
  })
  counts <- do.call(rbind, rows)
  rownames(counts) <- NULL
  write_section_table(counts, file.path(out_dir, "20_03_quadrant_counts.tsv"))
  invisible(counts)
}

#' Write the coordinates of the labelled key genes on every pair.
write_key_gene_coordinates <- function(results, out_dir) {
  rows <- lapply(results, function(result) {
    key_df <- result$key_genes
    if (nrow(key_df) == 0) return(NULL)
    data.frame(
      pair_id = result$spec$pair_id,
      figure = result$spec$file_stem,
      gene_name = key_df$gene_name,
      gene_id = key_df$gene_id,
      chr = key_df$chr,
      gene_length = key_df$gene_length,
      x_variable = result$spec$x_col,
      x_value = key_df[[result$spec$x_col]],
      y_variable = result$spec$y_col,
      y_value = key_df[[result$spec$y_col]],
      quadrant = key_df$quadrant,
      mch_diff = key_df$mch_diff,
      mch_sig = key_df$mch_sig,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    stop("No key gene from KEY_GENES falls inside any panel. ",
         "Check the gene names against the master table before interpreting ",
         "the figures.")
  }
  key_tbl <- do.call(rbind, rows)
  rownames(key_tbl) <- NULL
  write_section_table(key_tbl, file.path(out_dir, "20_03_key_gene_coordinates.tsv"))
  invisible(key_tbl)
}

#' Write the quadrant master table.
#'
#' Carries every gene of the master table with, for each of the five pairs, the
#' two coordinates and the quadrant. The three per-pair columns hold NA when
#' the gene is outside that pair's subset or missing a coordinate, so counting
#' the quadrant column reproduces the corner counts on the figure.
write_quadrant_master <- function(master, results, out_dir) {
  base_cols <- c("gene_name", "gene_id", "chr", "start", "end", "gene_length",
                 "promoter_state", "body_state", "chromatin_group",
                 "mch_ctrl", "mch_mut", "mch_diff", "edger_logFC", "edger_fdr",
                 "mch_sig", "mch_direction",
                 "k119ub_gb_ctrl_signal", "k119ub_gb_mut_signal",
                 "k119ub_gb_log2fc", "k119ub_gb_signal_class",
                 "k119ub_fold", "k119ub_fdr",
                 "mecp2_fold", "mecp2_fdr", "mecp2_n_peaks",
                 "k27ac_fold", "k27ac_fdr", "k27me3_fold", "k27me3_fdr")
  missing <- setdiff(base_cols, colnames(master))
  if (length(missing) > 0) {
    stop("write_quadrant_master(): the master table is missing columns: ",
         paste(missing, collapse = ", "))
  }

  out <- master[, base_cols, drop = FALSE]
  out$chromatin_group <- as.character(out$chromatin_group)

  for (result in results) {
    spec <- result$spec
    plotted <- result$plotted
    idx <- match(out$gene_name, plotted$gene_name)

    out[[paste0(spec$pair_id, "_x")]] <- plotted[[spec$x_col]][idx]
    out[[paste0(spec$pair_id, "_y")]] <- plotted[[spec$y_col]][idx]
    out[[paste0(spec$pair_id, "_quadrant")]] <- plotted$quadrant[idx]
  }

  write_section_table(out, file.path(out_dir, "20_03_quadrant_master.tsv"))
  invisible(out)
}

# =============================================================================
# GENE-LEVEL FISHER TESTS
# =============================================================================

#' Stop when a 2x2 table has an empty row or column.
#'
#' fisher.test() returns a zero, infinite or undefined odds ratio for such a
#' table. This check turns that into a named error instead.
#'
#' @param gene_df data.frame holding the two logical columns.
#' @param row_var Name of the row column.
#' @param col_var Name of the column column.
#' @param test_id Identifier used in the error message.
stop_on_degenerate_2x2 <- function(gene_df, row_var, col_var, test_id) {
  tab <- table(factor(gene_df[[row_var]], levels = c(TRUE, FALSE)),
               factor(gene_df[[col_var]], levels = c(TRUE, FALSE)))
  if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) {
    stop("Test ", test_id, " has an empty row or column in its 2x2 table: ",
         paste(as.vector(tab), collapse = ", "),
         ". Every gene shares one sign of ", row_var, " or ", col_var, ".")
  }
  invisible(tab)
}

#' Register the sign-concordance Fisher test behind one scatter.
#'
#' The two logical columns are the sign of the two plotted coordinates, so the
#' 2x2 table is the quadrant table of the figure: both TRUE is Q1, both FALSE
#' is Q3.
#'
#' @param result One element of run_all_scatters().
#' @param out_dir Section output directory.
#' @return one-row data.frame summarising the test
register_scatter_fisher <- function(result, out_dir) {
  spec <- result$spec
  plotted <- result$plotted

  x_flag <- paste0(spec$x_short, "_up")
  y_flag <- paste0(spec$y_short, "_up")

  gene_df <- data.frame(
    gene_name = plotted$gene_name,
    chr = plotted$chr,
    stringsAsFactors = FALSE
  )
  gene_df[[x_flag]] <- plotted[[spec$x_col]] > 0
  gene_df[[y_flag]] <- plotted[[spec$y_col]] > 0
  stop_on_degenerate_2x2(gene_df, x_flag, y_flag, spec$fisher_test_id)

  ft <- register_fisher_test(
    section = SECTION_ID,
    test_id = spec$fisher_test_id,
    description = spec$fisher_description,
    gene_df = gene_df, row_var = x_flag, col_var = y_flag,
    output_dir = out_dir
  )

  data.frame(
    test_id = spec$fisher_test_id,
    pair_id = spec$pair_id,
    gene_set = spec$subset_label,
    row_var = x_flag,
    col_var = y_flag,
    n_genes = nrow(gene_df),
    n_row_true = sum(gene_df[[x_flag]]),
    n_col_true = sum(gene_df[[y_flag]]),
    n_both_true = sum(gene_df[[x_flag]] & gene_df[[y_flag]]),
    odds_ratio = unname(ft$estimate),
    ci_low = ft$conf.int[1],
    ci_high = ft$conf.int[2],
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

#' Repeat the two mCH concordance tests over mCH-significant genes only.
#'
#' The scatter versions of these tests use every plotted gene, where the sign
#' of a non-significant mCH difference carries no genotype claim. These two
#' restrict the same 2x2 table to genes that pass the edgeR FDR cutoff.
#'
#' @param results Output of run_all_scatters().
#' @param out_dir Section output directory.
#' @return data.frame with one row per test
register_mch_significant_fisher <- function(results, out_dir) {
  specs <- list(
    list(key = "d", mark_flag = "K119ub_up", value_col = "k119ub_gb_log2fc",
         test_id = "k119ub_up_vs_mch_hyper_sig",
         description = paste("Among mCH-significant genes with quantifiable",
                             "H2AK119ub signal, does an H2AK119ub gain occur",
                             "with mCH hypermethylation?")),
    list(key = "e", mark_flag = "MeCP2_up", value_col = "mecp2_fold",
         test_id = "mecp2_up_vs_mch_hyper_sig",
         description = paste("Among mCH-significant genes with a MeCP2 peak,",
                             "does a MeCP2 gain occur with mCH",
                             "hypermethylation?"))
  )

  rows <- lapply(specs, function(spec) {
    plotted <- results[[spec$key]]$plotted
    sig <- plotted[plotted$mch_sig, , drop = FALSE]
    if (nrow(sig) < MIN_GENES_PER_SCATTER) {
      stop("Only ", nrow(sig), " mCH-significant genes remain for test ",
           spec$test_id, ". The section needs at least ",
           MIN_GENES_PER_SCATTER, ".")
    }

    gene_df <- data.frame(
      gene_name = sig$gene_name,
      chr = sig$chr,
      mch_hyper = sig$mch_diff > 0,
      stringsAsFactors = FALSE
    )
    gene_df[[spec$mark_flag]] <- sig[[spec$value_col]] > 0
    stop_on_degenerate_2x2(gene_df, spec$mark_flag, "mch_hyper", spec$test_id)

    ft <- register_fisher_test(
      section = SECTION_ID,
      test_id = spec$test_id,
      description = spec$description,
      gene_df = gene_df, row_var = spec$mark_flag, col_var = "mch_hyper",
      output_dir = out_dir
    )

    data.frame(
      test_id = spec$test_id,
      pair_id = results[[spec$key]]$spec$pair_id,
      gene_set = "mCH-significant genes only",
      row_var = spec$mark_flag,
      col_var = "mch_hyper",
      n_genes = nrow(gene_df),
      n_row_true = sum(gene_df[[spec$mark_flag]]),
      n_col_true = sum(gene_df$mch_hyper),
      n_both_true = sum(gene_df[[spec$mark_flag]] & gene_df$mch_hyper),
      odds_ratio = unname(ft$estimate),
      ci_low = ft$conf.int[1],
      ci_high = ft$conf.int[2],
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

#' Run and write every gene-level Fisher test of this section.
run_fisher_tests <- function(results, out_dir) {
  cat("--- Gene-level Fisher tests ---\n")
  scatter_rows <- do.call(rbind, lapply(results, register_scatter_fisher,
                                        out_dir = out_dir))
  sig_rows <- register_mch_significant_fisher(results, out_dir)

  summary_tbl <- rbind(scatter_rows, sig_rows)
  rownames(summary_tbl) <- NULL
  write_section_table(summary_tbl, file.path(out_dir, "20_03_fisher_summary.tsv"))
  cat("\n")
  summary_tbl
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opt <- parse_section_args()
  out_dir <- opt$output_dir
  fdr_threshold <- opt$fdr_threshold
  axis_clip <- opt$axis_clip

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("===========================================================================\n")
  cat("  SECTION 20_03: QUADRANT SCATTERS OF CHROMATIN MARKS AND mCH\n")
  cat("===========================================================================\n")
  cat("  Output directory:     ", out_dir, "\n", sep = "")
  cat("  MeCP2 FDR threshold:  ", fdr_threshold, "\n", sep = "")
  cat("  Axis clip percentile: ", axis_clip, "\n\n", sep = "")

  cat("--- Loading inputs ---\n")
  marks <- load_gene_level_all_marks(HANDOFF_PATHS$gene_level_all_marks)
  marks <- deduplicate_by_gene_name(marks)
  k119ub_signal <- load_k119ub_gene_signal(DIFFBIND_PATHS$k119ub_gene_signal)

  cat("  Aggregating MeCP2 DiffBind peaks to genes...\n")
  mecp2_genes <- aggregate_mecp2_to_genes(mecp2_diffbind, fdr_threshold)
  cat("\n")

  cat("--- Building the master table ---\n")
  master <- build_master_table(marks, k119ub_signal, mecp2_genes)
  write_chromatin_group_counts(master, out_dir)
  cat("\n")

  results <- run_all_scatters(master, axis_clip, out_dir)
  save_composite_figure(results, out_dir)

  cat("--- Writing tables ---\n")
  stats <- write_scatter_statistics(results, out_dir)
  write_quadrant_counts(results, out_dir)
  write_key_gene_coordinates(results, out_dir)
  write_quadrant_master(master, results, out_dir)
  cat("\n")

  fisher_tbl <- run_fisher_tests(results, out_dir)

  cat("---------------------------------------------------------------------------\n")
  cat("  SUMMARY\n")
  cat("---------------------------------------------------------------------------\n")
  for (i in seq_len(nrow(stats))) {
    cat(sprintf("  %-18s n = %7s  rho = %+.3f (p = %.3g)\n",
                stats$pair_id[i], count_label(stats$n_genes[i]),
                stats$spearman_rho[i], stats$spearman_p[i]))
  }
  cat("  Fisher tests (sign concordance):\n")
  for (i in seq_len(nrow(fisher_tbl))) {
    cat(sprintf("    %-40s OR = %.3f [%.3f, %.3f], p = %.3g (n = %s)\n",
                fisher_tbl$test_id[i], fisher_tbl$odds_ratio[i],
                fisher_tbl$ci_low[i], fisher_tbl$ci_high[i],
                fisher_tbl$p_value[i], count_label(fisher_tbl$n_genes[i])))
  }
  cat("\n=== Section 20_03 complete ===\n\n")
}

main()
