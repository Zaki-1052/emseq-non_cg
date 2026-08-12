# scripts/sections/60_mecp2/60_02_k119ub_unmethylated.R
#
# Section 60_02: H2AK119ub at genes where MeCP2 binding moves without mCH.
#
# What this tests
#   MeCP2 reads methylated cytosines, so a change in MeCP2 binding is usually
#   explained by a change in DNA methylation. Section 60_01 isolates the genes
#   that break that rule: MeCP2 binding changes at the FDR threshold while
#   gene-body mCH does not. This section asks what does change at those genes.
#
#   The answer is measured against two backgrounds drawn from the same gene
#   universe: genes where MeCP2 and mCH both change, and genes where neither
#   changes. Five chromatin marks enter the comparison, so H2AK119ub is judged
#   against H3K27me3, H3K27ac and ATAC rather than on its own.
#
#   "Unmethylated" throughout this script means the edgeR mCH test is not
#   significant at the section FDR threshold. It does not mean zero mCH.
#
# Marks quantified per gene
#   H2AK119ub gene-body control signal   absolute level, from the signal table
#   H2AK119ub gene-body signal log2FC    mutant over control, from the signal table
#   H2AK119ub peak log2FC                DiffBind peak nearest the TSS
#   H3K27me3, H3K27ac, ATAC peak log2FC  DiffBind peak nearest the TSS
#
# Analyses
#   1. Mann-Whitney U for every mark, comparing the MeCP2-without-mCH set with
#      each background. Two-sided, with Cliff's delta and its confidence
#      interval as the effect size.
#   2. Rank of the two H2AK119ub measures among the five log2FC marks, by the
#      absolute Cliff's delta of each background comparison.
#   3. Two registered gene-level Fisher tests: membership in the
#      MeCP2-without-mCH set against H2AK119ub gain, measured once from the
#      gene-body signal and once from significant peak gain.
#   4. Logistic model of set membership on the four peak-based mark log2FCs,
#      fitted on raw predictors and on standardised predictors so the marks can
#      be compared with each other.
#   5. Spearman correlation of MeCP2 fold change against H2AK119ub gene-body
#      log2FC inside the MeCP2-without-mCH set.
#
# Reads
#   mch_results, mecp2_diffbind, DIFFBIND_TABLES  pre-loaded by _shared_config.R
#   HANDOFF_PATHS$mecp2_no_meth_genes             written by section 60_01
#   MECP2_PATHS$annotated                         per-peak MeCP2 gene annotation
#   DIFFBIND_PATHS$k119ub_gene_signal             gene-body H2AK119ub signal
#
# Writes to results/sections/60_mecp2/ (OUTPUT_PATHS$mecp2, override with
#   --output-dir): seven multi-format figures, twelve TSV tables, and two
#   registered Fisher gene tables under fisher_tables/.
#
# Adapted from Biomodal section 67, MeCP2 binding at unmethylated genes.

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)
library(effsize)
library(pROC)

# =============================================================================
# SECTION CONSTANTS
# =============================================================================

SECTION_ID <- "60_02"

# Smallest number of genes a group may contribute to any test.
MIN_GROUP_N <- 5

# Table-safe group identifiers used in every TSV.
GROUP_ORDER <- c("MeCP2_only", "MeCP2_and_mCH", "Neither")

# Plot-facing labels, one per GROUP_ORDER entry, in the same order. The line
# breaks keep the axis readable and never reach a TSV: tables carry
# GROUP_TABLE_LABELS instead.
GROUP_LABELS <- c(
  MeCP2_only    = "MeCP2 changes\n(no mCH change)",
  MeCP2_and_mCH = "MeCP2 and mCH\n(both change)",
  Neither       = "Neither changes"
)

GROUP_TABLE_LABELS <- gsub("\n", " ", GROUP_LABELS, fixed = TRUE)

GROUP_COLORS <- c(
  "MeCP2 changes\n(no mCH change)" = unname(COLORS$mecp2[["MeCP2 Up"]]),
  "MeCP2 and mCH\n(both change)"   = unname(COLORS$mecp2[["MeCP2 Down"]]),
  "Neither changes"                = "grey70"
)

FOCAL_GROUP <- "MeCP2_only"
BACKGROUND_GROUPS <- c("MeCP2_and_mCH", "Neither")

# One sentence per background, used in figure subtitles and table columns.
BACKGROUND_MEANING <- c(
  MeCP2_and_mCH = "genes where MeCP2 binding and mCH both change",
  Neither       = "genes where neither MeCP2 binding nor mCH changes"
)

# Marks measured on a log2 fold-change scale. These share the violin figure and
# the effect-size ranking.
FOLD_METRICS <- c("k119ub_fold", "k119ub_gb_log2fc",
                  "k27me3_fold", "k27ac_fold", "atac_fold")

# Marks measured on an absolute signal scale. Own units, own figure.
ABSOLUTE_METRICS <- c("k119ub_gb_ctrl_signal")

# Peak-derived log2 fold changes that enter the logistic model together.
MODEL_METRICS <- c("k119ub_fold", "k27me3_fold", "k27ac_fold", "atac_fold")

# display  panel title and forest-plot row label
# short    compact label for log lines and ranking tables
# color    panel accent, taken from the shared COLORS palettes
# axis     y-axis label of the panel showing this metric
METRIC_META <- list(
  k119ub_fold = list(
    display = "H2AK119ub (peak log2FC)", short = "H2AK119ub peak",
    color = unname(COLORS$k119ub[["K119ub Gained"]]),
    axis = "H2AK119ub peak log2FC (mutant / control)"
  ),
  k119ub_gb_log2fc = list(
    display = "H2AK119ub (gene-body signal log2FC)", short = "H2AK119ub signal",
    color = "#54278F",
    axis = "H2AK119ub gene-body log2FC (mutant / control)"
  ),
  k27me3_fold = list(
    display = "H3K27me3 (peak log2FC)", short = "H3K27me3",
    color = unname(COLORS$h3k27me3[["H3K27me3 Gained"]]),
    axis = "H3K27me3 peak log2FC (mutant / control)"
  ),
  k27ac_fold = list(
    display = "H3K27ac (peak log2FC)", short = "H3K27ac",
    color = unname(COLORS$h3k27ac[["H3K27ac Gained"]]),
    axis = "H3K27ac peak log2FC (mutant / control)"
  ),
  atac_fold = list(
    display = "ATAC (peak log2FC)", short = "ATAC",
    color = unname(COLORS$atac[["ATAC Up"]]),
    axis = "ATAC peak log2FC (mutant / control)"
  ),
  k119ub_gb_ctrl_signal = list(
    display = "H2AK119ub (control gene-body signal)", short = "H2AK119ub level",
    color = "#9E9AC8",
    axis = "H2AK119ub gene-body signal, control"
  )
)

# DiffBind marks aggregated to genes here. k119ub also supplies the peak-based
# H2AK119ub fold used by the model.
PEAK_MARKS <- c("k119ub", "k27me3", "k27ac", "atac")

PEAK_MARK_DISPLAY <- c(
  k119ub = "H2AK119ub",
  k27me3 = "H3K27me3",
  k27ac  = "H3K27ac",
  atac   = "ATAC"
)

# Columns section 60_01 writes into mecp2_no_meth_genes.tsv.
HANDOFF_COLUMNS <- c(
  "gene_name", "chr", "start", "end", "gene_length",
  "mch_diff", "edger_fdr", "mch_sig",
  "mecp2_fold", "mecp2_fdr", "mecp2_sig", "mecp2_direction"
)

N_DECILES <- 10

# =============================================================================
# COMMAND LINE
# =============================================================================

parse_section_args <- function() {
  option_list <- list(
    make_option("--output-dir", dest = "output_dir", type = "character",
                default = OUTPUT_PATHS$mecp2,
                help = "Directory for figures and tables [default: %default]"),
    make_option("--fdr-threshold", dest = "fdr_threshold", type = "double",
                default = Q_THRESHOLD,
                help = paste("FDR cutoff for peak significance. Must match the",
                             "value section 60_01 ran with [default: %default]"))
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

fmt_int <- function(x) format(x, big.mark = ",", trim = TRUE)

significance_stars <- function(p) {
  ifelse(is.na(p), "ns",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", "ns"))))
}

metric_display <- function(metric) METRIC_META[[metric]]$display
metric_short <- function(metric) METRIC_META[[metric]]$short
metric_color <- function(metric) METRIC_META[[metric]]$color
metric_axis <- function(metric) METRIC_META[[metric]]$axis

#' Turn table-safe group identifiers into the plot-facing factor.
group_label_factor <- function(group) {
  factor(unname(GROUP_LABELS[group]), levels = unname(GROUP_LABELS[GROUP_ORDER]))
}

# =============================================================================
# INPUT LOADING
# =============================================================================

#' Read the MeCP2-without-mCH gene list written by section 60_01.
#'
#' @param path HANDOFF_PATHS$mecp2_no_meth_genes.
#' @return data.frame with the HANDOFF_COLUMNS
load_handoff_genes <- function(path) {
  if (!file.exists(path)) {
    stop("MeCP2-without-mCH gene list not found: ", path,
         "\nRun section 60_01 (60_01_methylation_scale.R) first. ",
         "It writes mecp2_no_meth_genes.tsv.")
  }

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")
  missing <- setdiff(HANDOFF_COLUMNS, colnames(df))
  if (length(missing) > 0) {
    stop("Handoff table ", path, " is missing columns: ",
         paste(missing, collapse = ", "),
         "\nSection 60_01 must write every column in HANDOFF_COLUMNS.")
  }
  if (nrow(df) == 0) {
    stop("Handoff table ", path, " has no rows. Section 60_02 needs a ",
         "non-empty MeCP2-without-mCH gene set.")
  }
  if (anyDuplicated(df$gene_name) > 0) {
    stop("Handoff table ", path, " repeats gene names. Section 60_01 must ",
         "write one row per gene.")
  }

  cat(sprintf("  MeCP2-without-mCH genes from 60_01: %s\n", fmt_int(nrow(df))))
  cat(sprintf("    MeCP2 gained: %s | MeCP2 lost: %s\n",
              fmt_int(sum(df$mecp2_direction == "Gained")),
              fmt_int(sum(df$mecp2_direction == "Lost"))))
  df
}

#' Read the per-peak MeCP2 gene annotation.
#'
#' @param path MECP2_PATHS$annotated.
#' @return data.frame of peaks carrying a gene symbol
load_mecp2_annotated <- function(path) {
  if (!file.exists(path)) {
    stop("MeCP2 annotated peak file not found: ", path)
  }

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "", fill = TRUE)

  required <- c("SYMBOL", "Fold", "FDR", "distanceToTSS")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("MeCP2 annotated file is missing columns: ",
         paste(missing, collapse = ", "), " (", path, ")")
  }

  for (col in c("Fold", "FDR", "distanceToTSS")) {
    df[[col]] <- as.numeric(df[[col]])
  }

  usable <- !is.na(df$SYMBOL) & nzchar(df$SYMBOL) & df$SYMBOL != "NA" &
    !is.na(df$Fold) & !is.na(df$FDR) & !is.na(df$distanceToTSS)
  if (sum(usable) == 0) {
    stop("No MeCP2 annotated peak carries a gene symbol and complete ",
         "Fold, FDR and distanceToTSS values: ", path)
  }

  cat(sprintf("  MeCP2 annotated peaks: %s total, %s usable\n",
              fmt_int(nrow(df)), fmt_int(sum(usable))))
  df[usable, , drop = FALSE]
}

#' Read the H2AK119ub gene-body signal table.
#'
#' Keeps the gene-body columns and renames them with a k119ub prefix. When a
#' gene symbol appears on more than one row, the row with the largest gene
#' width is kept and the number of collapsed symbols is reported.
#'
#' @param path DIFFBIND_PATHS$k119ub_gene_signal.
#' @return data.frame with one row per gene symbol
load_k119ub_gene_signal <- function(path) {
  if (!file.exists(path)) {
    stop("H2AK119ub gene-body signal table not found: ", path)
  }

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")

  required <- c("symbol", "chr", "start", "end", "width",
                "gb_ctrl_signal", "gb_mut_signal", "gb_log2fc", "gb_signal_class")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("H2AK119ub gene signal table is missing columns: ",
         paste(missing, collapse = ", "), " (", path, ")")
  }

  for (col in c("width", "gb_ctrl_signal", "gb_mut_signal", "gb_log2fc")) {
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  }

  n_all <- nrow(df)
  keep <- !is.na(df$symbol) & nzchar(df$symbol) & df$symbol != "NA"
  df <- df[keep, , drop = FALSE]

  ord <- order(df$symbol, -df$width)
  df <- df[ord, , drop = FALSE]
  n_dup <- sum(duplicated(df$symbol))
  df <- df[!duplicated(df$symbol), , drop = FALSE]

  out <- data.frame(
    gene_name = df$symbol,
    k119ub_gb_ctrl_signal = df$gb_ctrl_signal,
    k119ub_gb_mut_signal = df$gb_mut_signal,
    k119ub_gb_log2fc = df$gb_log2fc,
    k119ub_gb_class = df$gb_signal_class,
    stringsAsFactors = FALSE
  )

  cat(sprintf("  H2AK119ub gene signal: %s rows in, %s genes kept (%s duplicate symbols collapsed)\n",
              fmt_int(n_all), fmt_int(nrow(out)), fmt_int(n_dup)))
  cat(sprintf("    Finite gene-body log2FC: %s genes | signal classes: %s\n",
              fmt_int(sum(is.finite(out$k119ub_gb_log2fc))),
              paste(sprintf("%s=%s", names(table(out$k119ub_gb_class)),
                            as.integer(table(out$k119ub_gb_class))),
                    collapse = ", ")))
  out
}

#' Row indices of mch_results that keep one row per gene symbol.
#'
#' The row with the largest absolute edgeR log fold change is kept for each
#' gene name. Section 60_01 deduplicates the same way, so the two sections
#' describe the same genes.
#'
#' @param mch mch_results data.frame.
#' @return integer vector of row indices in ascending order
deduplicate_mch_row_indices <- function(mch) {
  ord <- order(mch$gene_name, -abs(mch$edger_logFC))
  keep <- ord[!duplicated(mch$gene_name[ord])]
  sort(keep)
}

# =============================================================================
# GENE-LEVEL MARK TABLES
# =============================================================================

#' Collapse a DiffBind table to one fold change per gene.
#'
#' Peaks are annotated with ChIPseeker, then the peak nearest the TSS supplies
#' the gene-level fold change and FDR.
#'
#' @param db data.frame from load_diffbind_flex().
#' @param mark Short mark identifier used as the column prefix.
#' @param fdr_threshold FDR cutoff for the per-gene significant peak count.
#' @return data.frame with one row per gene
annotate_mark_to_genes <- function(db, mark, fdr_threshold) {
  display <- PEAK_MARK_DISPLAY[[mark]]

  db <- db[db$Chr %in% CANONICAL_CHRS & !is.na(db$Fold) & !is.na(db$FDR), ,
           drop = FALSE]
  if (nrow(db) == 0) {
    stop("No usable ", display, " peaks on canonical chromosomes with complete ",
         "Fold and FDR values.")
  }

  annotated <- annotate_peaks_to_genes(db, display)
  usable <- !is.na(annotated$SYMBOL) & nzchar(annotated$SYMBOL) &
    !is.na(annotated$Fold) & !is.na(annotated$FDR)
  annotated <- annotated[usable, , drop = FALSE]
  if (nrow(annotated) == 0) {
    stop("No ", display, " peak was annotated to a gene symbol.")
  }

  gene_table <- aggregate_diffbind_by_gene(
    annotated, method = "nearest_tss",
    fdr_threshold = fdr_threshold, prefix = mark
  )
  if (anyDuplicated(gene_table$gene_name) > 0) {
    stop("aggregate_diffbind_by_gene() returned duplicate genes for ", display)
  }

  cat(sprintf("    %-10s %s peaks annotated, %s genes carry a peak\n",
              display, fmt_int(nrow(annotated)), fmt_int(nrow(gene_table))))
  gene_table
}

#' Aggregate every DiffBind mark in PEAK_MARKS to genes.
#'
#' @param fdr_threshold FDR cutoff passed to the aggregation.
#' @return named list of gene-level data.frames
annotate_all_peak_marks <- function(fdr_threshold) {
  tables <- list()
  for (mark in PEAK_MARKS) {
    tables[[mark]] <- annotate_mark_to_genes(DIFFBIND_TABLES[[mark]], mark,
                                             fdr_threshold)
  }
  tables
}

#' Collapse MeCP2 annotated peaks to one fold change per gene.
#'
#' Uses the same nearest-TSS rule as section 60_01, so the group assignment
#' below reproduces the handoff gene set.
#'
#' @param annotated data.frame from load_mecp2_annotated().
#' @param fdr_threshold FDR cutoff for the per-gene significant peak count.
#' @return data.frame with one row per gene, prefix mecp2
aggregate_mecp2_by_gene <- function(annotated, fdr_threshold) {
  out <- aggregate_diffbind_by_gene(annotated, method = "nearest_tss",
                                    fdr_threshold = fdr_threshold,
                                    prefix = "mecp2")
  if (anyDuplicated(out$gene_name) > 0) {
    stop("aggregate_diffbind_by_gene() returned duplicate genes for MeCP2.")
  }
  cat(sprintf("  Genes with a MeCP2 annotated peak: %s\n", fmt_int(nrow(out))))
  out
}

# =============================================================================
# GENE UNIVERSE AND GROUP ASSIGNMENT
# =============================================================================

#' Build the one-row-per-gene table carrying groups and every mark.
#'
#' The universe is the deduplicated mCH-tested genes that also carry a MeCP2
#' peak, because MeCP2 significance is undefined without a peak. Every gene in
#' that universe falls into one of four cells; the three cells this section
#' compares are named in GROUP_ORDER and the fourth (mCH changes, MeCP2 does
#' not) is carried as "mCH_only" and excluded from the comparisons.
#'
#' @param mch Deduplicated mch_results rows.
#' @param mecp2_gene Gene-level MeCP2 table.
#' @param peak_tables Named list of gene-level DiffBind tables.
#' @param k119ub_signal Gene-level H2AK119ub signal table.
#' @param fdr_threshold FDR cutoff defining a significant MeCP2 change.
#' @return data.frame with one row per gene in the universe
build_gene_table <- function(mch, mecp2_gene, peak_tables, k119ub_signal,
                             fdr_threshold) {
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
    stringsAsFactors = FALSE
  )

  tbl <- dplyr::left_join(tbl, mecp2_gene, by = "gene_name")
  tbl <- tbl[!is.na(tbl$mecp2_fold) & !is.na(tbl$mecp2_fdr), , drop = FALSE]
  if (nrow(tbl) == 0) {
    stop("No mCH-tested gene name matched a MeCP2 annotated peak. Check that ",
         MECP2_PATHS$annotated, " uses the same gene symbols as ",
         DATA_PATHS$mch_results)
  }

  tbl$mecp2_sig <- tbl$mecp2_fdr < fdr_threshold
  tbl$mecp2_direction <- ifelse(!tbl$mecp2_sig, "Unchanged",
                                ifelse(tbl$mecp2_fold > 0, "Gained", "Lost"))

  tbl$group <- ifelse(tbl$mecp2_sig & !tbl$mch_sig, "MeCP2_only",
                 ifelse(tbl$mecp2_sig & tbl$mch_sig, "MeCP2_and_mCH",
                   ifelse(!tbl$mecp2_sig & !tbl$mch_sig, "Neither", "mCH_only")))

  for (mark in PEAK_MARKS) {
    tbl <- dplyr::left_join(tbl, peak_tables[[mark]], by = "gene_name")
  }
  tbl <- dplyr::left_join(tbl, k119ub_signal, by = "gene_name")

  if (anyDuplicated(tbl$gene_name) > 0) {
    stop("The gene table repeats gene names after joining the mark tables. ",
         "A joined table has duplicate keys.")
  }

  cat(sprintf("  Universe: %s genes with an mCH test and a MeCP2 peak\n",
              fmt_int(nrow(tbl))))
  cell_counts <- table(factor(tbl$group,
                              levels = c(GROUP_ORDER, "mCH_only")))
  for (nm in names(cell_counts)) {
    cat(sprintf("    %-14s %8s genes (%5.2f%%)\n", nm,
                fmt_int(as.integer(cell_counts[[nm]])),
                100 * cell_counts[[nm]] / nrow(tbl)))
  }
  for (metric in c(FOLD_METRICS, ABSOLUTE_METRICS)) {
    cat(sprintf("    %-24s measured at %s genes\n", metric_short(metric),
                fmt_int(sum(is.finite(tbl[[metric]])))))
  }
  tbl
}

#' Stop unless the reconstructed MeCP2-without-mCH set equals the handoff set.
#'
#' Section 60_01 built its handoff list with the same nearest-TSS rule and the
#' same deduplication. The two sets differ only when the two sections ran with
#' different FDR thresholds or different input files.
#'
#' @param gene_tbl Gene table from build_gene_table().
#' @param handoff data.frame from load_handoff_genes().
#' @param fdr_threshold FDR cutoff this section is using.
check_handoff_consistency <- function(gene_tbl, handoff, fdr_threshold) {
  reconstructed <- gene_tbl$gene_name[gene_tbl$group == FOCAL_GROUP]
  only_in_handoff <- setdiff(handoff$gene_name, reconstructed)
  only_in_reconstruction <- setdiff(reconstructed, handoff$gene_name)

  if (length(only_in_handoff) > 0 || length(only_in_reconstruction) > 0) {
    stop("The MeCP2-without-mCH set rebuilt here does not match the set ",
         "section 60_01 wrote.\n",
         "  In the 60_01 handoff only: ", length(only_in_handoff),
         " genes (first: ", paste(head(only_in_handoff, 5), collapse = ", "), ")\n",
         "  Rebuilt here only: ", length(only_in_reconstruction),
         " genes (first: ", paste(head(only_in_reconstruction, 5), collapse = ", "), ")\n",
         "  This section ran with --fdr-threshold ", fdr_threshold,
         ". Run 60_01 and 60_02 with the same threshold and the same inputs.")
  }

  cat(sprintf("  Handoff check passed: %s genes match section 60_01 exactly\n",
              fmt_int(length(reconstructed))))
  invisible(TRUE)
}

#' Restrict the gene table to the three compared groups.
#'
#' @param gene_tbl Gene table from build_gene_table().
#' @return data.frame carrying group as a factor and in_mecp2_no_mch
build_analysis_table <- function(gene_tbl) {
  df <- gene_tbl[gene_tbl$group %in% GROUP_ORDER, , drop = FALSE]
  df$group <- factor(df$group, levels = GROUP_ORDER)
  df$in_mecp2_no_mch <- df$group == FOCAL_GROUP

  counts <- table(df$group)
  for (g in GROUP_ORDER) {
    if (counts[[g]] < MIN_GROUP_N) {
      stop("Group ", g, " has ", counts[[g]], " genes, below the minimum of ",
           MIN_GROUP_N, ". The three-group comparison cannot run.")
    }
  }

  cat(sprintf("  Analysis table: %s genes across the three compared groups\n",
              fmt_int(nrow(df))))
  df
}

#' Count the genes each group contributes to each metric.
#'
#' @param df Analysis table.
#' @return data.frame with one row per group and metric
build_group_counts <- function(df) {
  rows <- list()
  for (g in GROUP_ORDER) {
    sub <- df[df$group == g, , drop = FALSE]
    for (metric in c(FOLD_METRICS, ABSOLUTE_METRICS)) {
      values <- sub[[metric]]
      values <- values[is.finite(values)]
      rows[[length(rows) + 1]] <- data.frame(
        group = g,
        group_label = unname(GROUP_TABLE_LABELS[[g]]),
        metric = metric,
        metric_display = metric_display(metric),
        n_genes_in_group = nrow(sub),
        n_with_metric = length(values),
        pct_with_metric = 100 * length(values) / nrow(sub),
        median = if (length(values) > 0) median(values) else NA_real_,
        mean = if (length(values) > 0) mean(values) else NA_real_,
        q25 = if (length(values) > 0) unname(quantile(values, 0.25)) else NA_real_,
        q75 = if (length(values) > 0) unname(quantile(values, 0.75)) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# =============================================================================
# ANALYSIS 1: MANN-WHITNEY U WITH EFFECT SIZES
# =============================================================================

#' Compare one metric between the focal group and one background.
#'
#' The test is a two-sided Mann-Whitney U. Direction is reported by the median
#' difference and by Cliff's delta, which is positive when the focal group
#' carries the larger values.
#'
#' @param df Analysis table.
#' @param metric Column name of the numeric measure.
#' @param background One entry of BACKGROUND_GROUPS.
#' @return data.frame with one row
compare_metric <- function(df, metric, background) {
  focal <- df[[metric]][df$group == FOCAL_GROUP]
  other <- df[[metric]][df$group == background]
  focal <- focal[is.finite(focal)]
  other <- other[is.finite(other)]

  if (length(focal) < MIN_GROUP_N || length(other) < MIN_GROUP_N) {
    stop("Metric ", metric, " has ", length(focal), " focal and ",
         length(other), " ", background, " genes with a finite value. ",
         "Both must reach ", MIN_GROUP_N, ".")
  }

  wt <- wilcox.test(focal, other, alternative = "two.sided", exact = FALSE)
  cd <- effsize::cliff.delta(focal, other)

  data.frame(
    metric = metric,
    metric_display = metric_display(metric),
    metric_short = metric_short(metric),
    comparison = sprintf("%s vs %s", FOCAL_GROUP, background),
    background = background,
    background_meaning = unname(BACKGROUND_MEANING[[background]]),
    n_focal = length(focal),
    n_background = length(other),
    median_focal = median(focal),
    median_background = median(other),
    median_difference = median(focal) - median(other),
    mean_focal = mean(focal),
    mean_background = mean(other),
    U = unname(wt$statistic),
    p_value = wt$p.value,
    cliffs_delta = unname(cd$estimate),
    delta_ci_low = unname(cd$conf.int[1]),
    delta_ci_high = unname(cd$conf.int[2]),
    delta_magnitude = as.character(cd$magnitude),
    stringsAsFactors = FALSE
  )
}

#' Run compare_metric() for every metric against every background.
#'
#' @param df Analysis table.
#' @return data.frame with one row per metric and background
run_all_comparisons <- function(df) {
  rows <- list()
  for (metric in c(FOLD_METRICS, ABSOLUTE_METRICS)) {
    for (background in BACKGROUND_GROUPS) {
      rows[[length(rows) + 1]] <- compare_metric(df, metric, background)
    }
  }
  out <- do.call(rbind, rows)
  out$fdr <- p.adjust(out$p_value, method = "BH")
  out$stars <- significance_stars(out$fdr)

  for (i in seq_len(nrow(out))) {
    cat(sprintf("  %-22s vs %-14s med %+8.4f vs %+8.4f | delta = %+.3f [%+.3f, %+.3f] (%s) | p = %.3g, BH = %.3g\n",
                out$metric_short[i], out$background[i],
                out$median_focal[i], out$median_background[i],
                out$cliffs_delta[i], out$delta_ci_low[i], out$delta_ci_high[i],
                out$delta_magnitude[i], out$p_value[i], out$fdr[i]))
  }
  out
}

# =============================================================================
# ANALYSIS 2: RANK OF H2AK119ub AMONG THE MARKS
# =============================================================================

#' Rank the log2FC marks by absolute Cliff's delta inside each comparison.
#'
#' @param comparisons data.frame from run_all_comparisons().
#' @return data.frame with one row per metric and background, carrying rank
build_effect_size_ranking <- function(comparisons) {
  df <- comparisons[comparisons$metric %in% FOLD_METRICS, , drop = FALSE]
  df$abs_delta <- abs(df$cliffs_delta)

  out <- df %>%
    dplyr::group_by(background) %>%
    dplyr::mutate(
      rank = rank(-abs_delta, ties.method = "min"),
      n_marks_ranked = dplyr::n()
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(background, rank) %>%
    dplyr::select(background, background_meaning, rank, n_marks_ranked,
                  metric, metric_short, metric_display,
                  cliffs_delta, abs_delta, delta_ci_low, delta_ci_high,
                  delta_magnitude, median_difference, p_value, fdr) %>%
    as.data.frame()

  for (background in BACKGROUND_GROUPS) {
    cat(sprintf("  Ranked against %s:\n", background))
    sub <- out[out$background == background, , drop = FALSE]
    for (i in seq_len(nrow(sub))) {
      cat(sprintf("    %d. %-22s |delta| = %.3f (%s)\n",
                  sub$rank[i], sub$metric_short[i], sub$abs_delta[i],
                  sub$delta_magnitude[i]))
    }
  }
  out
}

#' One sentence naming where the H2AK119ub measures rank in each comparison.
#'
#' @param ranking data.frame from build_effect_size_ranking().
#' @return single string
ranking_sentence <- function(ranking) {
  parts <- vapply(BACKGROUND_GROUPS, function(background) {
    sub <- ranking[ranking$background == background, , drop = FALSE]
    peak <- sub[sub$metric == "k119ub_fold", , drop = FALSE]
    signal <- sub[sub$metric == "k119ub_gb_log2fc", , drop = FALSE]
    sprintf("against %s, H2AK119ub ranks %d (signal) and %d (peak) of %d marks",
            background, signal$rank[1], peak$rank[1], sub$n_marks_ranked[1])
  }, character(1))
  paste(parts, collapse = "; ")
}

# =============================================================================
# ANALYSIS 3: REGISTERED FISHER TESTS
# =============================================================================

#' Register one Fisher test of set membership against an H2AK119ub gain flag.
#'
#' @param df Analysis table carrying in_mecp2_no_mch and the gain column.
#' @param gain_col Logical column naming H2AK119ub gain.
#' @param test_id Identifier inside this section.
#' @param description What the test asks.
#' @param out_dir Section output directory.
#' @return data.frame with one row summarising the test
register_gain_fisher <- function(df, gain_col, test_id, description, out_dir) {
  sub <- df[!is.na(df[[gain_col]]), , drop = FALSE]
  if (sum(sub$in_mecp2_no_mch) < MIN_GROUP_N) {
    stop("Only ", sum(sub$in_mecp2_no_mch), " focal genes carry ", gain_col,
         ". The Fisher test ", test_id, " cannot run.")
  }

  ft <- register_fisher_test(
    section = SECTION_ID, test_id = test_id, description = description,
    gene_df = sub, row_var = "in_mecp2_no_mch", col_var = gain_col,
    output_dir = out_dir)

  focal <- sub[sub$in_mecp2_no_mch, , drop = FALSE]
  rest <- sub[!sub$in_mecp2_no_mch, , drop = FALSE]

  out <- data.frame(
    test_id = test_id,
    row_var = "in_mecp2_no_mch",
    col_var = gain_col,
    n_genes = nrow(sub),
    n_focal = nrow(focal),
    n_background = nrow(rest),
    n_focal_gained = sum(focal[[gain_col]]),
    n_background_gained = sum(rest[[gain_col]]),
    pct_focal_gained = 100 * mean(focal[[gain_col]]),
    pct_background_gained = 100 * mean(rest[[gain_col]]),
    odds_ratio = unname(ft$estimate),
    ci_low = ft$conf.int[1],
    ci_high = ft$conf.int[2],
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )

  cat(sprintf("    Focal gained %.1f%% (%s/%s), background gained %.1f%% (%s/%s)\n",
              out$pct_focal_gained, fmt_int(out$n_focal_gained), fmt_int(out$n_focal),
              out$pct_background_gained, fmt_int(out$n_background_gained),
              fmt_int(out$n_background)))
  out
}

#' Add the two H2AK119ub gain flags to the analysis table.
#'
#' Signal gain is a positive gene-body log2FC. Peak gain is a significant
#' positive DiffBind fold change at the section FDR threshold. Both are NA for
#' a gene the matching measure does not cover.
#'
#' @param df Analysis table.
#' @param fdr_threshold FDR cutoff for peak gain.
#' @return the analysis table with k119ub_signal_gain and k119ub_peak_gain
add_k119ub_gain_flags <- function(df, fdr_threshold) {
  df$k119ub_signal_gain <- ifelse(is.finite(df$k119ub_gb_log2fc),
                                  df$k119ub_gb_log2fc > 0, NA)
  df$k119ub_peak_gain <- ifelse(
    !is.na(df$k119ub_fold) & !is.na(df$k119ub_fdr),
    df$k119ub_fdr < fdr_threshold & df$k119ub_fold > 0, NA)

  cat(sprintf("  H2AK119ub signal gain measured at %s genes, peak gain at %s genes\n",
              fmt_int(sum(!is.na(df$k119ub_signal_gain))),
              fmt_int(sum(!is.na(df$k119ub_peak_gain)))))
  df
}

#' Register both Fisher tests of set membership against H2AK119ub gain.
#'
#' @param df Analysis table carrying both gain flags.
#' @param out_dir Section output directory.
#' @return data.frame with one row per test
run_gain_fisher_tests <- function(df, out_dir) {
  signal_row <- register_gain_fisher(
    df, "k119ub_signal_gain", "mecp2_no_mch_vs_k119ub_signal_gain",
    paste("Do genes where MeCP2 binding changes without an mCH change gain",
          "H2AK119ub gene-body signal more often than other genes with a",
          "MeCP2 peak?"),
    out_dir)

  peak_row <- register_gain_fisher(
    df, "k119ub_peak_gain", "mecp2_no_mch_vs_k119ub_peak_gain",
    paste("Do genes where MeCP2 binding changes without an mCH change carry a",
          "significantly gained H2AK119ub peak more often than other genes",
          "with a MeCP2 peak?"),
    out_dir)

  rbind(signal_row, peak_row)
}

# =============================================================================
# ANALYSIS 4: LOGISTIC MODEL OF SET MEMBERSHIP
# =============================================================================

#' Odds ratios with Wald confidence intervals from a fitted glm.
#'
#' @param fit glm object.
#' @return data.frame with one row per coefficient
extract_odds_ratios <- function(fit) {
  coefs <- coef(fit)
  conf <- confint.default(fit)
  smry <- summary(fit)$coefficients
  data.frame(
    term = names(coefs),
    estimate = unname(coefs),
    std_error = unname(smry[, 2]),
    z_value = unname(smry[, 3]),
    p_value = unname(smry[, 4]),
    odds_ratio = exp(unname(coefs)),
    or_lower = exp(unname(conf[, 1])),
    or_upper = exp(unname(conf[, 2])),
    stringsAsFactors = FALSE
  )
}

#' Fit set membership on the four peak-based mark log2FCs.
#'
#' Two fits run on the same genes. The raw fit gives odds ratios per unit
#' log2FC. The standardised fit divides each predictor by its own standard
#' deviation, so the four odds ratios are on one scale and can be ranked.
#'
#' @param df Analysis table.
#' @return list with both fits, the model data and the fit statistics
fit_membership_model <- function(df) {
  complete <- stats::complete.cases(df[, MODEL_METRICS, drop = FALSE])
  d <- df[complete, , drop = FALSE]
  d$outcome <- as.integer(d$in_mecp2_no_mch)

  if (length(unique(d$outcome)) < 2) {
    stop("Set membership has one class only among the ", nrow(d),
         " genes carrying all four peak marks.")
  }
  if (sum(d$outcome) < MIN_GROUP_N) {
    stop("Only ", sum(d$outcome), " focal genes carry all four peak marks. ",
         "The logistic model cannot run.")
  }

  model_formula <- as.formula(paste("outcome ~",
                                    paste(MODEL_METRICS, collapse = " + ")))
  raw_fit <- glm(model_formula, data = d, family = binomial)

  sds <- vapply(MODEL_METRICS, function(col) sd(d[[col]]), numeric(1))
  if (any(!is.finite(sds)) || any(sds == 0)) {
    stop("A predictor has zero or non-finite standard deviation among the ",
         nrow(d), " model genes: ",
         paste(sprintf("%s=%.4g", MODEL_METRICS, sds), collapse = ", "))
  }
  z <- d
  for (col in MODEL_METRICS) z[[col]] <- (z[[col]] - mean(z[[col]])) / sds[[col]]
  std_fit <- glm(model_formula, data = z, family = binomial)

  null_fit <- glm(outcome ~ 1, data = d, family = binomial)
  mcfadden <- 1 - as.numeric(logLik(raw_fit)) / as.numeric(logLik(null_fit))

  probs <- predict(raw_fit, type = "response")
  roc_obj <- pROC::roc(d$outcome, probs, quiet = TRUE)
  auc_ci <- as.numeric(pROC::ci.auc(roc_obj))

  cat(sprintf("  %s genes carry all four peak marks, %s in the focal set (%.2f%%)\n",
              fmt_int(nrow(d)), fmt_int(sum(d$outcome)),
              100 * mean(d$outcome)))
  cat(sprintf("  AUC = %.3f [%.3f, %.3f], McFadden R2 = %.4f\n",
              auc_ci[2], auc_ci[1], auc_ci[3], mcfadden))

  list(raw_fit = raw_fit, std_fit = std_fit, data = d, sds = sds,
       mcfadden = mcfadden, auc = auc_ci[2], auc_lower = auc_ci[1],
       auc_upper = auc_ci[3], n_genes = nrow(d), n_focal = sum(d$outcome))
}

#' Build the coefficient table of one fit, with mark labels attached.
#'
#' @param fit glm object.
#' @param scale_label "raw" or "standardised".
#' @param sds Predictor standard deviations from the model data.
#' @return data.frame with one row per mark
build_coefficient_table <- function(fit, scale_label, sds) {
  or_df <- extract_odds_ratios(fit)
  or_df <- or_df[or_df$term != "(Intercept)", , drop = FALSE]
  or_df$scale <- scale_label
  or_df$metric_short <- unname(vapply(or_df$term, metric_short, character(1)))
  or_df$metric_display <- unname(vapply(or_df$term, metric_display, character(1)))
  or_df$predictor_sd <- unname(sds[or_df$term])
  or_df$stars <- significance_stars(or_df$p_value)
  or_df[, c("scale", "term", "metric_short", "metric_display", "predictor_sd",
            "estimate", "std_error", "z_value", "p_value",
            "odds_ratio", "or_lower", "or_upper", "stars")]
}

#' Print the standardised coefficients ordered by how far they sit from 1.
print_coefficient_ranking <- function(coef_df) {
  ord <- order(-abs(coef_df$estimate))
  for (i in ord) {
    cat(sprintf("  %-22s OR = %.3f [%.3f, %.3f] per 1 SD, p = %.3g %s\n",
                coef_df$metric_short[i], coef_df$odds_ratio[i],
                coef_df$or_lower[i], coef_df$or_upper[i],
                coef_df$p_value[i], coef_df$stars[i]))
  }
}

# =============================================================================
# ANALYSIS 5: MeCP2 FOLD AGAINST H2AK119ub INSIDE THE FOCAL SET
# =============================================================================

#' Spearman correlation of MeCP2 fold change with H2AK119ub gene-body log2FC.
#'
#' @param df Analysis table.
#' @return list with the correlation table and the genes it was computed on
correlate_mecp2_with_k119ub <- function(df) {
  sub <- df[df$group == FOCAL_GROUP &
              is.finite(df$k119ub_gb_log2fc) & is.finite(df$mecp2_fold), ,
            drop = FALSE]

  if (nrow(sub) < MIN_GROUP_N) {
    stop("Only ", nrow(sub), " focal genes carry both a MeCP2 fold change and ",
         "an H2AK119ub gene-body log2FC. The correlation cannot run.")
  }

  ct <- suppressWarnings(cor.test(sub$mecp2_fold, sub$k119ub_gb_log2fc,
                                  method = "spearman"))
  tbl <- data.frame(
    set = FOCAL_GROUP,
    x = "mecp2_fold",
    y = "k119ub_gb_log2fc",
    n_genes = nrow(sub),
    spearman_rho = unname(ct$estimate),
    S = unname(ct$statistic),
    p_value = ct$p.value,
    stringsAsFactors = FALSE
  )

  cat(sprintf("  Spearman rho = %+.3f, p = %.3g over %s focal genes\n",
              tbl$spearman_rho, tbl$p_value, fmt_int(tbl$n_genes)))
  list(summary = tbl, genes = sub)
}

# =============================================================================
# FIGURE 60_02a AND 60_02b: GROUPED VIOLINS
# =============================================================================

#' One violin panel of a metric across the three gene groups.
#'
#' The panel carries the gene count and the median of every group, and the two
#' Mann-Whitney p-values comparing the focal group with each background.
#'
#' @param df Analysis table.
#' @param metric Column name of the numeric measure.
#' @param comparisons data.frame from run_all_comparisons().
#' @param zero_line TRUE draws a dashed line at zero.
#' @return ggplot object
violin_panel <- function(df, metric, comparisons, zero_line = TRUE) {
  d <- df[is.finite(df[[metric]]), c("group", metric), drop = FALSE]
  names(d)[2] <- "value"
  d$group_label <- group_label_factor(as.character(d$group))

  grp <- summarise_groups(d, "group_label", "value")
  grp$label <- group_label(grp)

  # The view is clipped to the central 99% so one extreme gene cannot flatten
  # the violins. Every gene still enters the violin and the summary numbers.
  y_lo <- unname(quantile(d$value, 0.005))
  y_hi <- unname(quantile(d$value, 0.995))
  span <- y_hi - y_lo
  if (!is.finite(span) || span == 0) {
    stop("Metric ", metric, " has no spread across the plotted genes.")
  }
  label_y <- y_lo - 0.07 * span
  stat_y <- y_hi + 0.16 * span

  stats_sub <- comparisons[comparisons$metric == metric, , drop = FALSE]
  stat_label <- paste(
    sprintf("vs %s: p = %.2e %s",
            stats_sub$background, stats_sub$p_value, stats_sub$stars),
    collapse = "\n")

  p <- ggplot(d, aes(x = group_label, y = value, fill = group_label)) +
    geom_violin(alpha = 0.6, scale = "width",
                draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white", alpha = 0.85) +
    geom_text(data = grp, aes(x = group_label, y = label_y, label = label),
              inherit.aes = FALSE, size = 2.9, lineheight = 1.1) +
    annotate("text", x = 2, y = stat_y, label = stat_label,
             size = 2.9, fontface = "italic", lineheight = 1.2) +
    scale_fill_manual(values = GROUP_COLORS, guide = "none") +
    coord_cartesian(ylim = c(y_lo - 0.18 * span, y_hi + 0.30 * span)) +
    labs(title = metric_display(metric), x = "", y = metric_axis(metric)) +
    theme_emseq(base_size = 11) +
    theme(plot.title = element_text(color = metric_color(metric), size = 12),
          axis.text.x = element_text(size = 8.5))

  if (zero_line) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "grey40",
                        linewidth = 0.4)
  }
  p
}

#' Draw the five log2FC violins as one figure.
#'
#' @param df Analysis table.
#' @param comparisons data.frame from run_all_comparisons().
#' @param ranking data.frame from build_effect_size_ranking().
#' @param out_dir Section output directory.
#' @return the combined ggplot object
plot_fold_violins <- function(df, comparisons, ranking, out_dir) {
  panels <- lapply(FOLD_METRICS, function(metric) {
    violin_panel(df, metric, comparisons, zero_line = TRUE)
  })

  n_focal <- sum(df$group == FOCAL_GROUP)
  combined <- patchwork::wrap_plots(panels, ncol = 3) +
    plot_annotation(
      title = "Chromatin Mark Change at Genes Where MeCP2 Moves Without mCH",
      subtitle = sprintf(
        "%s focal genes against %s genes where both change and %s where neither changes. %s.",
        fmt_int(n_focal),
        fmt_int(sum(df$group == "MeCP2_and_mCH")),
        fmt_int(sum(df$group == "Neither")),
        ranking_sentence(ranking)),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 9.5, color = "grey30")
      )
    )

  save_multiformat_ggplot(combined,
                          file.path(out_dir, "60_02a_mark_folds_by_group"),
                          width = 18, height = 12)
  combined
}

#' Draw the absolute H2AK119ub control-signal violin.
#'
#' @param df Analysis table.
#' @param comparisons data.frame from run_all_comparisons().
#' @param out_dir Section output directory.
#' @return the ggplot object
plot_absolute_signal_violin <- function(df, comparisons, out_dir) {
  metric <- "k119ub_gb_ctrl_signal"
  p <- violin_panel(df, metric, comparisons, zero_line = FALSE) +
    labs(
      title = "Control H2AK119ub Level at Genes Where MeCP2 Moves Without mCH",
      subtitle = paste("Absolute gene-body signal in control, before any",
                       "mutant comparison")
    ) +
    theme_emseq() +
    theme(plot.title = element_text(color = metric_color(metric)))

  save_multiformat_ggplot(p,
                          file.path(out_dir, "60_02b_k119ub_control_level"),
                          width = 9, height = 8)
  p
}

# =============================================================================
# FIGURE 60_02c: EFFECT-SIZE FOREST
# =============================================================================

#' Forest plot of Cliff's delta for every metric and background.
#'
#' @param comparisons data.frame from run_all_comparisons().
#' @param ranking data.frame from build_effect_size_ranking().
#' @param out_dir Section output directory.
#' @return the ggplot object
plot_effect_size_forest <- function(comparisons, ranking, out_dir) {
  d <- comparisons
  d$metric_row <- factor(
    d$metric_short,
    levels = rev(unname(vapply(c(FOLD_METRICS, ABSOLUTE_METRICS),
                               metric_short, character(1)))))
  d$background_label <- factor(
    sprintf("vs %s", d$background),
    levels = sprintf("vs %s", BACKGROUND_GROUPS))

  background_colors <- setNames(
    c(unname(GROUP_COLORS[[unname(GROUP_LABELS[["MeCP2_and_mCH"]])]]),
      "grey40"),
    sprintf("vs %s", BACKGROUND_GROUPS))

  # One facet per background, so every metric keeps its own row and no two
  # intervals share a y position.
  p <- ggplot(d, aes(x = cliffs_delta, y = metric_row,
                     color = background_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_errorbar(aes(xmin = delta_ci_low, xmax = delta_ci_high),
                  orientation = "y", width = 0.25, linewidth = 0.7) +
    geom_point(size = 3) +
    geom_text(aes(label = sprintf("%+.3f %s", cliffs_delta, stars)),
              vjust = -1.1, size = 2.9, show.legend = FALSE) +
    facet_wrap(~ background_label, nrow = 1) +
    scale_color_manual(values = background_colors, name = "Background") +
    scale_y_discrete(expand = expansion(add = 0.6)) +
    labs(
      title = "Effect Size of Each Mark at Genes Where MeCP2 Moves Without mCH",
      subtitle = sprintf(
        paste("Cliff's delta with 95%% confidence interval. Positive means the",
              "focal set carries the larger values. Stars mark the BH-adjusted",
              "Mann-Whitney p-value. %s."),
        ranking_sentence(ranking)),
      x = "Cliff's delta (focal set minus background)", y = ""
    ) +
    theme_emseq() +
    theme(axis.text.y = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"))

  save_multiformat_ggplot(p,
                          file.path(out_dir, "60_02c_effect_size_forest"),
                          width = 13, height = 8)
  p
}

# =============================================================================
# FIGURE 60_02d: LOGISTIC COEFFICIENT FOREST
# =============================================================================

#' Forest plot of the standardised logistic odds ratios.
#'
#' @param model list from fit_membership_model().
#' @param coef_df Standardised coefficient table.
#' @param out_dir Section output directory.
#' @return the ggplot object
plot_logistic_forest <- function(model, coef_df, out_dir) {
  d <- coef_df
  d$mark_row <- factor(d$metric_short,
                       levels = rev(unname(vapply(MODEL_METRICS, metric_short,
                                                  character(1)))))
  d$point_color <- unname(vapply(d$term, metric_color, character(1)))

  p <- ggplot(d, aes(x = odds_ratio, y = mark_row)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    geom_errorbar(aes(xmin = or_lower, xmax = or_upper), orientation = "y",
                  width = 0.2, linewidth = 0.8) +
    geom_point(aes(fill = mark_row), size = 4, shape = 21, color = "grey20") +
    geom_text(aes(label = sprintf("OR = %.3f [%.3f, %.3f] %s", odds_ratio,
                                  or_lower, or_upper, stars)),
              vjust = -1.2, size = 3.1) +
    scale_fill_manual(values = setNames(d$point_color, d$mark_row),
                      guide = "none") +
    scale_x_log10() +
    scale_y_discrete(expand = expansion(add = 0.7)) +
    labs(
      title = "Which Mark Predicts MeCP2 Change Without mCH Change",
      subtitle = sprintf(
        paste("Logistic model on standardised peak log2FCs. n = %s genes",
              "carrying all four marks, %s in the focal set. AUC = %.3f",
              "[%.3f, %.3f], McFadden R2 = %.4f."),
        fmt_int(model$n_genes), fmt_int(model$n_focal),
        model$auc, model$auc_lower, model$auc_upper, model$mcfadden),
      x = "Odds ratio per 1 SD of log2FC (log scale)", y = ""
    ) +
    theme_emseq() +
    theme(axis.text.y = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"))

  save_multiformat_ggplot(p,
                          file.path(out_dir, "60_02d_logistic_forest"),
                          width = 11, height = 7)
  p
}

# =============================================================================
# FIGURE 60_02e: H2AK119ub DECILE DISTRIBUTION
# =============================================================================

#' Assign genome-wide deciles of one metric and count each group per decile.
#'
#' Deciles come from dplyr::ntile() over every gene in the analysis table with
#' a finite value, so each decile holds the same number of genes overall.
#'
#' @param df Analysis table.
#' @param metric Column name of the numeric measure.
#' @return data.frame with one row per group and decile
build_decile_table <- function(df, metric) {
  d <- df[is.finite(df[[metric]]), c("group", "gene_name", metric),
          drop = FALSE]
  names(d)[3] <- "value"
  if (nrow(d) < N_DECILES) {
    stop("Metric ", metric, " has ", nrow(d), " finite values, fewer than the ",
         N_DECILES, " deciles.")
  }

  d$decile <- dplyr::ntile(d$value, N_DECILES)

  boundaries <- d %>%
    dplyr::group_by(decile) %>%
    dplyr::summarise(decile_min = min(value), decile_max = max(value),
                     n_universe = dplyr::n(), .groups = "drop") %>%
    as.data.frame()

  counts <- d %>%
    dplyr::group_by(group, decile) %>%
    dplyr::summarise(n_genes = dplyr::n(), .groups = "drop") %>%
    as.data.frame()

  # Every group gets a row for every decile, so a group that reaches no gene in
  # a decile plots as zero rather than dropping out of the panel.
  out <- expand.grid(
    group = factor(GROUP_ORDER, levels = GROUP_ORDER),
    decile = seq_len(N_DECILES),
    stringsAsFactors = FALSE
  )
  out <- dplyr::left_join(out, counts, by = c("group", "decile"))
  out$n_genes[is.na(out$n_genes)] <- 0L

  out <- dplyr::left_join(out, boundaries, by = "decile")
  out$n_group <- ave(out$n_genes, out$group, FUN = sum)
  out$pct_of_group <- 100 * out$n_genes / out$n_group
  out$metric <- metric
  out$metric_display <- metric_display(metric)
  out$group_label <- unname(GROUP_TABLE_LABELS[as.character(out$group)])

  out[order(out$group, out$decile), , drop = FALSE]
}

#' One bar panel of the decile distribution.
#'
#' @param decile_tbl data.frame from build_decile_table().
#' @param metric Column name, used for the panel title.
#' @return ggplot object
decile_panel <- function(decile_tbl, metric) {
  d <- decile_tbl
  d$group_plot <- group_label_factor(as.character(d$group))

  n_by_group <- d %>%
    dplyr::distinct(group_plot, n_group) %>%
    dplyr::arrange(group_plot)
  legend_labels <- setNames(
    sprintf("%s\nn = %s", as.character(n_by_group$group_plot),
            fmt_int(n_by_group$n_group)),
    as.character(n_by_group$group_plot))

  ggplot(d, aes(x = factor(decile), y = pct_of_group, fill = group_plot)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75,
             alpha = 0.9) +
    geom_hline(yintercept = 100 / N_DECILES, linetype = "dashed",
               color = "grey30", linewidth = 0.4) +
    scale_fill_manual(values = GROUP_COLORS, labels = legend_labels,
                      name = "Gene group") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title = metric_display(metric),
      subtitle = sprintf(
        "A group that follows the whole universe sits on the dashed %.0f%% line",
        100 / N_DECILES),
      x = sprintf("%s decile (1 = lowest)", metric_short(metric)),
      y = "Percentage of the group"
    ) +
    theme_emseq(base_size = 11) +
    theme(plot.title = element_text(color = metric_color(metric), size = 12),
          plot.subtitle = element_text(size = 8.5, color = "grey40"))
}

#' Draw the decile distribution for the gene-body log2FC and the control level.
#'
#' @param df Analysis table.
#' @param out_dir Section output directory.
#' @return data.frame holding both decile tables
plot_decile_distribution <- function(df, out_dir) {
  metrics <- c("k119ub_gb_log2fc", "k119ub_gb_ctrl_signal")
  tables <- lapply(metrics, function(metric) build_decile_table(df, metric))
  names(tables) <- metrics

  panels <- lapply(metrics, function(metric) {
    decile_panel(tables[[metric]], metric)
  })

  combined <- (panels[[1]] / panels[[2]]) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "H2AK119ub Decile Distribution by Gene Group",
      subtitle = paste("Deciles are computed over every gene in the universe,",
                       "so each decile holds one tenth of the genes"),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30")
      )
    )

  save_multiformat_ggplot(combined,
                          file.path(out_dir, "60_02e_k119ub_deciles"),
                          width = 13, height = 11)

  for (metric in metrics) {
    focal <- tables[[metric]][tables[[metric]]$group == FOCAL_GROUP, ,
                              drop = FALSE]
    top_pct <- sum(focal$pct_of_group[focal$decile >= N_DECILES - 1])
    cat(sprintf("  %-22s focal genes in the top two deciles: %.1f%%\n",
                metric_short(metric), top_pct))
  }

  do.call(rbind, tables)
}

# =============================================================================
# FIGURE 60_02f: MeCP2 FOLD AGAINST H2AK119ub SCATTER
# =============================================================================

#' Scatter of MeCP2 fold change against H2AK119ub gene-body log2FC.
#'
#' Only the focal set is plotted. Labels go on the KEY_GENES that are present
#' and on the genes with the largest H2AK119ub gain.
#'
#' @param correlation list from correlate_mecp2_with_k119ub().
#' @param out_dir Section output directory.
#' @return the ggplot object
plot_mecp2_k119ub_scatter <- function(correlation, out_dir) {
  d <- correlation$genes
  tbl <- correlation$summary

  top_k119ub <- d$gene_name[order(-d$k119ub_gb_log2fc)][seq_len(min(10, nrow(d)))]
  label_genes <- unique(c(intersect(KEY_GENES, d$gene_name), top_k119ub))
  d$label_gene <- ifelse(d$gene_name %in% label_genes, d$gene_name, "")

  p <- ggplot(d, aes(x = k119ub_gb_log2fc, y = mecp2_fold)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40",
               linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40",
               linewidth = 0.3) +
    geom_point(aes(color = mecp2_direction), alpha = 0.65, size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "grey20",
                linewidth = 0.6) +
    geom_text_repel(aes(label = label_gene), size = 3, max.overlaps = 20,
                    fontface = "italic", color = "grey15",
                    segment.color = "grey60", segment.size = 0.3,
                    min.segment.length = 0) +
    scale_color_manual(
      values = c("Gained" = unname(COLORS$mecp2[["MeCP2 Up"]]),
                 "Lost" = unname(COLORS$mecp2[["MeCP2 Down"]])),
      name = "MeCP2 direction") +
    labs(
      title = "MeCP2 Fold Change against H2AK119ub at Genes Without an mCH Change",
      subtitle = sprintf("Spearman rho = %+.3f, p = %.2e, n = %s genes",
                         tbl$spearman_rho, tbl$p_value, fmt_int(tbl$n_genes)),
      x = "H2AK119ub gene-body log2FC (mutant / control)",
      y = "MeCP2 log2 fold change (mutant / control)"
    ) +
    theme_emseq()

  save_multiformat_ggplot(p,
                          file.path(out_dir, "60_02f_mecp2_vs_k119ub_scatter"),
                          width = 11, height = 9)
  p
}

# =============================================================================
# FIGURE 60_02g: COMBINED PANEL
# =============================================================================

#' Assemble the effect-size forest, the model forest and the scatter.
#'
#' @param p_forest Effect-size forest.
#' @param p_model Logistic coefficient forest.
#' @param p_scatter MeCP2 against H2AK119ub scatter.
#' @param out_dir Section output directory.
plot_combined_panel <- function(p_forest, p_model, p_scatter, out_dir) {
  combined <- (p_forest / (p_model | p_scatter)) +
    plot_annotation(
      title = "MeCP2 Responds to H2AK119ub Where mCH Does Not Change",
      subtitle = paste("Effect size of every mark, the multivariate model of",
                       "set membership, and MeCP2 against H2AK119ub inside",
                       "the focal set"),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 17),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30")
      )
    )

  save_multiformat_ggplot(combined,
                          file.path(out_dir, "60_02g_combined"),
                          width = 17, height = 15)
  invisible(combined)
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
  cat("SECTION 60_02: H2AK119ub AT GENES WHERE MeCP2 MOVES WITHOUT mCH\n")
  cat("================================================================================\n")
  cat("Output dir:    ", OUT_DIR, "\n", sep = "")
  cat("FDR threshold: ", fdr_threshold, "\n\n", sep = "")

  stopifnot(
    "MeCP2 annotated peak file not found" = file.exists(MECP2_PATHS$annotated),
    "H2AK119ub gene signal table not found" =
      file.exists(DIFFBIND_PATHS$k119ub_gene_signal)
  )

  cat("--- Loading section inputs ---\n")
  handoff <- load_handoff_genes(HANDOFF_PATHS$mecp2_no_meth_genes)
  mecp2_annotated <- load_mecp2_annotated(MECP2_PATHS$annotated)
  k119ub_signal <- load_k119ub_gene_signal(DIFFBIND_PATHS$k119ub_gene_signal)

  cat("\n--- Aggregating MeCP2 peaks to genes (nearest TSS) ---\n")
  mecp2_gene <- aggregate_mecp2_by_gene(mecp2_annotated, fdr_threshold)

  cat("\n--- Aggregating DiffBind marks to genes (nearest TSS) ---\n")
  peak_tables <- annotate_all_peak_marks(fdr_threshold)

  cat("\n--- Building the gene universe ---\n")
  keep_idx <- deduplicate_mch_row_indices(mch_results)
  cat(sprintf("  Deduplicated %s mCH rows to %s gene names\n",
              fmt_int(nrow(mch_results)), fmt_int(length(keep_idx))))
  gene_tbl <- build_gene_table(mch_results[keep_idx, , drop = FALSE],
                               mecp2_gene, peak_tables, k119ub_signal,
                               fdr_threshold)
  check_handoff_consistency(gene_tbl, handoff, fdr_threshold)

  analysis_tbl <- build_analysis_table(gene_tbl)
  analysis_tbl <- add_k119ub_gain_flags(analysis_tbl, fdr_threshold)

  write_section_table(gene_tbl, file.path(OUT_DIR, "60_02_gene_universe.tsv"))
  write_section_table(analysis_tbl,
                      file.path(OUT_DIR, "60_02_gene_groups_with_marks.tsv"))

  group_counts <- build_group_counts(analysis_tbl)
  write_section_table(group_counts,
                      file.path(OUT_DIR, "60_02_group_metric_counts.tsv"))

  cat("\n--- Analysis 1: Mann-Whitney U with effect sizes ---\n")
  comparisons <- run_all_comparisons(analysis_tbl)
  write_section_table(comparisons, file.path(OUT_DIR, "60_02_mann_whitney.tsv"))

  cat("\n--- Analysis 2: rank of H2AK119ub among the marks ---\n")
  ranking <- build_effect_size_ranking(comparisons)
  write_section_table(ranking,
                      file.path(OUT_DIR, "60_02_effect_size_ranking.tsv"))

  cat("\n--- Analysis 3: registered Fisher tests ---\n")
  fisher_summary <- run_gain_fisher_tests(analysis_tbl, OUT_DIR)
  write_section_table(fisher_summary,
                      file.path(OUT_DIR, "60_02_fisher_summary.tsv"))

  cat("\n--- Analysis 4: logistic model of set membership ---\n")
  model <- fit_membership_model(analysis_tbl)
  raw_coefs <- build_coefficient_table(model$raw_fit, "raw", model$sds)
  std_coefs <- build_coefficient_table(model$std_fit, "standardised", model$sds)
  print_coefficient_ranking(std_coefs)
  write_section_table(rbind(std_coefs, raw_coefs),
                      file.path(OUT_DIR, "60_02_logistic_coefficients.tsv"))
  write_section_table(
    data.frame(
      n_genes = model$n_genes,
      n_focal = model$n_focal,
      pct_focal = 100 * model$n_focal / model$n_genes,
      auc = model$auc, auc_lower = model$auc_lower, auc_upper = model$auc_upper,
      mcfadden_r2 = model$mcfadden,
      predictors = paste(MODEL_METRICS, collapse = ","),
      stringsAsFactors = FALSE
    ),
    file.path(OUT_DIR, "60_02_logistic_fit.tsv"))

  cat("\n--- Analysis 5: MeCP2 against H2AK119ub inside the focal set ---\n")
  correlation <- correlate_mecp2_with_k119ub(analysis_tbl)
  write_section_table(correlation$summary,
                      file.path(OUT_DIR, "60_02_mecp2_k119ub_correlation.tsv"))

  cat("\n--- Figure 60_02a: mark folds by gene group ---\n")
  plot_fold_violins(analysis_tbl, comparisons, ranking, OUT_DIR)

  cat("\n--- Figure 60_02b: control H2AK119ub level by gene group ---\n")
  plot_absolute_signal_violin(analysis_tbl, comparisons, OUT_DIR)

  cat("\n--- Figure 60_02c: effect-size forest ---\n")
  p_forest <- plot_effect_size_forest(comparisons, ranking, OUT_DIR)

  cat("\n--- Figure 60_02d: logistic coefficient forest ---\n")
  p_model <- plot_logistic_forest(model, std_coefs, OUT_DIR)

  cat("\n--- Figure 60_02e: H2AK119ub decile distribution ---\n")
  decile_tbl <- plot_decile_distribution(analysis_tbl, OUT_DIR)
  write_section_table(decile_tbl, file.path(OUT_DIR, "60_02_k119ub_deciles.tsv"))

  cat("\n--- Figure 60_02f: MeCP2 against H2AK119ub scatter ---\n")
  p_scatter <- plot_mecp2_k119ub_scatter(correlation, OUT_DIR)

  cat("\n--- Figure 60_02g: combined panel ---\n")
  plot_combined_panel(p_forest, p_model, p_scatter, OUT_DIR)

  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 60_02 SUMMARY\n")
  cat("================================================================================\n")
  cat(sprintf("Gene universe (mCH tested and a MeCP2 peak): %s\n",
              fmt_int(nrow(gene_tbl))))
  for (g in GROUP_ORDER) {
    cat(sprintf("  %-14s %8s genes\n", g, fmt_int(sum(analysis_tbl$group == g))))
  }

  for (background in BACKGROUND_GROUPS) {
    sub <- ranking[ranking$background == background, , drop = FALSE]
    cat(sprintf("Effect-size order against %s: %s\n", background,
                paste(sprintf("%d) %s |d|=%.3f", sub$rank, sub$metric_short,
                              sub$abs_delta), collapse = "  ")))
  }

  for (i in seq_len(nrow(fisher_summary))) {
    row <- fisher_summary[i, ]
    cat(sprintf("Fisher %-36s OR = %.3f (95%% CI %.3f to %.3f), p = %.3g\n",
                row$test_id, row$odds_ratio, row$ci_low, row$ci_high,
                row$p_value))
  }

  cat(sprintf("Logistic model: n = %s genes, %s focal, AUC = %.3f, McFadden R2 = %.4f\n",
              fmt_int(model$n_genes), fmt_int(model$n_focal), model$auc,
              model$mcfadden))
  top_coef <- std_coefs[which.max(abs(std_coefs$estimate)), ]
  cat(sprintf("Largest standardised coefficient: %s (OR = %.3f per 1 SD, p = %.3g)\n",
              top_coef$metric_short, top_coef$odds_ratio, top_coef$p_value))
  cat(sprintf("MeCP2 against H2AK119ub in the focal set: rho = %+.3f, p = %.3g\n",
              correlation$summary$spearman_rho, correlation$summary$p_value))
  cat(sprintf("Fisher tests registered in: %s\n", HANDOFF_PATHS$fisher_registry))
  cat("\nSection 60_02 complete.\n\n")
}

main()
