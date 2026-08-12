# scripts/sections/60_mecp2/60_01_methylation_scale.R
#
# Section 60_01: the scale of mCH change against the scale of MeCP2 binding change.
#
# What this tests
#   MeCP2 reads methylated cytosines, and in neurons most of the methylation it
#   reads is gene-body mCH. If mCH change alone drove MeCP2 recruitment, the
#   number of genes that change MeCP2 binding would track the number that change
#   mCH. This section counts both and compares them.
#
#   The cascade runs: genes tested, genes with a significant mCH change, genes
#   that gain mCH, those carrying a MeCP2 peak, and of those the genes whose
#   MeCP2 binding also changes. Every step reports its count, the percentage of
#   the step before it, and the percentage of the tested universe.
#
# Analyses
#   1. Cascade table with per-step counts and percentages.
#   2. Two-panel hierarchy bars: the mCH scale beside the MeCP2 scale.
#   3. Two-set Venn: mCH significant against MeCP2 significantly changed.
#   4. Quadrant scatter of mCH change against MeCP2 fold change, assigned with
#      assign_quadrant() and carrying the gene count in each corner.
#   5. Funnel figure of the cascade.
#   6. Violin of MeCP2 fold change across the three mCH direction groups, with
#      pairwise Wilcoxon tests.
#   7. Registered gene-level Fisher test: mCH significance against MeCP2
#      significance, over the genes that carry a MeCP2 peak.
#
# Gene-level MeCP2 values come from the peak nearest the TSS. A gene counts as
# MeCP2-changed when that nearest peak has FDR below the threshold.
#
# Reads
#   mch_results, mecp2_diffbind   pre-loaded by the shared config
#   MECP2_PATHS$annotated         per-peak MeCP2 gene annotation
#
# Writes to results/sections/60_mecp2/ (OUTPUT_PATHS$mecp2, override with
#   --output-dir): five multi-format figures, nine TSV tables, one registered
#   Fisher gene table under fisher_tables/, and the handoff table
#   mecp2_no_meth_genes.tsv that section 60_02 reads.
#
# Adapted from Biomodal section 65, methylation scale against MeCP2 occupancy.

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)
library(ggVennDiagram)

# =============================================================================
# SECTION CONSTANTS
# =============================================================================

SECTION_ID <- "60_01"

MCH_CATEGORY_LEVELS <- c("Hypermethylated", "Hypomethylated", "Not Significant")

MCH_CATEGORY_COLORS <- c(
  "Hypermethylated" = unname(COLORS$direction["Hypermethylated"]),
  "Hypomethylated"  = unname(COLORS$direction["Hypomethylated"]),
  "Not Significant" = "grey70"
)

# Fill ramps for the two hierarchy panels, light to dark down each cascade.
MCH_RAMP   <- c("#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D")
MECP2_RAMP <- c("#DADAEB", "#BCBDDC", "#807DBA", "#54278F")

MCH_PANEL_COLOR   <- "#CB181D"
MECP2_PANEL_COLOR <- "#54278F"

# Fill ramp for the five funnel bars, widest to narrowest.
FUNNEL_RAMP <- c("#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26", "#A50F15")

# What each assign_quadrant() label means on the mCH against MeCP2 scatter.
QUADRANT_MEANING <- c(
  Q1 = "mCH up, MeCP2 up",
  Q2 = "mCH down, MeCP2 up",
  Q3 = "mCH down, MeCP2 down",
  Q4 = "mCH up, MeCP2 down"
)

# Column contract for mecp2_no_meth_genes.tsv. Section 60_02 reads this file
# and indexes these columns by name.
HANDOFF_COLUMNS <- c(
  "gene_name", "chr", "start", "end", "gene_length",
  "mch_diff", "edger_fdr", "mch_sig",
  "mecp2_fold", "mecp2_fdr", "mecp2_sig", "mecp2_direction"
)

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

fmt_int <- function(x) format(x, big.mark = ",", trim = TRUE)

percent_of <- function(numerator, denominator) {
  if (denominator == 0) stop("percent_of(): denominator is zero.")
  100 * numerator / denominator
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
#' @param fdr_threshold FDR cutoff used for the peak-level significance message.
#' @return data.frame of annotated peaks carrying a gene symbol
load_mecp2_annotated <- function(filepath, fdr_threshold) {
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
  if (nrow(df) == 0) {
    stop("No MeCP2 annotated peak carries a gene symbol: ", filepath)
  }

  cat(sprintf("  MeCP2 annotated peaks: %s total, %s with a gene symbol\n",
              fmt_int(n_all), fmt_int(nrow(df))))
  cat(sprintf("  Peak-level significance at FDR<%.2f: %s up, %s down\n",
              fdr_threshold,
              fmt_int(sum(df$FDR < fdr_threshold & df$Fold > 0)),
              fmt_int(sum(df$FDR < fdr_threshold & df$Fold < 0))))
  df
}

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

# =============================================================================
# GENE TABLE
# =============================================================================

#' Collapse annotated MeCP2 peaks to one fold change per gene.
#'
#' Uses the shared nearest-TSS collapse rule, so the gene carries the fold
#' change and FDR of the peak closest to its TSS.
#'
#' @param annotated data.frame from load_mecp2_annotated().
#' @param fdr_threshold FDR cutoff for the per-gene significant peak count.
#' @return data.frame with one row per gene, prefix mecp2
aggregate_mecp2_by_gene <- function(annotated, fdr_threshold) {
  out <- aggregate_diffbind_by_gene(annotated, method = "nearest_tss",
                                    fdr_threshold = fdr_threshold,
                                    prefix = "mecp2")
  cat(sprintf("  Genes with a MeCP2 annotated peak: %s\n", fmt_int(nrow(out))))
  out
}

#' Build the one-row-per-gene analysis table.
#'
#' Joins the deduplicated mCH results with the gene-level MeCP2 summary and
#' derives the significance and direction flags that every analysis below uses.
#'
#' @param mch Deduplicated mch_results rows.
#' @param mecp2_gene Gene-level MeCP2 table from aggregate_mecp2_by_gene().
#' @param fdr_threshold FDR cutoff defining a significant MeCP2 change.
#' @return data.frame with one row per gene
build_gene_table <- function(mch, mecp2_gene, fdr_threshold) {
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

  tbl$mch_category <- factor(
    ifelse(tbl$mch_hyper, "Hypermethylated",
      ifelse(tbl$mch_hypo, "Hypomethylated", "Not Significant")),
    levels = MCH_CATEGORY_LEVELS
  )

  tbl <- dplyr::left_join(tbl, mecp2_gene, by = "gene_name")

  tbl$has_mecp2 <- !is.na(tbl$mecp2_fold) & !is.na(tbl$mecp2_fdr)

  # MeCP2 significance is undefined for a gene with no peak, so it stays NA.
  tbl$mecp2_sig <- ifelse(tbl$has_mecp2, tbl$mecp2_fdr < fdr_threshold, NA)
  tbl$mecp2_direction <- ifelse(
    !tbl$has_mecp2, NA_character_,
    ifelse(!tbl$mecp2_sig, "Unchanged",
      ifelse(tbl$mecp2_fold > 0, "Gained", "Lost")))

  if (sum(tbl$has_mecp2) == 0) {
    stop("No mCH-tested gene name matched a MeCP2 annotated peak. Check that ",
         MECP2_PATHS$annotated, " uses the same gene symbols as ",
         DATA_PATHS$mch_results)
  }

  cat(sprintf("  Gene table: %s genes, %s carry a MeCP2 peak\n",
              fmt_int(nrow(tbl)), fmt_int(sum(tbl$has_mecp2))))
  cat(sprintf("  MeCP2 gene-level change at FDR<%.2f: %s gained, %s lost\n",
              fdr_threshold,
              fmt_int(sum(tbl$mecp2_direction == "Gained", na.rm = TRUE)),
              fmt_int(sum(tbl$mecp2_direction == "Lost", na.rm = TRUE))))
  tbl
}

# =============================================================================
# ANALYSIS 1: THE CASCADE
# =============================================================================

#' Count the five cascade steps from the tested universe down to MeCP2 change.
#'
#' Each step is a subset of the step before it. The table reports the count, the
#' percentage of the previous step, and the percentage of the tested universe.
#'
#' @param gene_tbl Gene table from build_gene_table().
#' @return data.frame with one row per cascade step
build_cascade_table <- function(gene_tbl) {
  hyper <- gene_tbl$mch_hyper
  hyper_with_peak <- hyper & gene_tbl$has_mecp2
  hyper_with_change <- hyper_with_peak & gene_tbl$mecp2_sig %in% TRUE

  counts <- c(
    nrow(gene_tbl),
    sum(gene_tbl$mch_sig),
    sum(hyper),
    sum(hyper_with_peak),
    sum(hyper_with_change)
  )

  labels <- c(
    "Genes tested",
    "mCH significant",
    "mCH hypermethylated",
    "Hypermethylated with a MeCP2 peak",
    "Hypermethylated with a MeCP2 change"
  )

  short_labels <- c(
    "Genes\ntested",
    "mCH\nsignificant",
    "mCH\nhypermethylated",
    "Hyper +\nMeCP2 peak",
    "Hyper +\nMeCP2 change"
  )

  previous <- c(counts[1], counts[-length(counts)])

  data.frame(
    step = seq_along(counts),
    step_label = labels,
    short_label = short_labels,
    n_genes = counts,
    pct_of_previous = 100 * counts / previous,
    pct_of_tested = 100 * counts / counts[1],
    stringsAsFactors = FALSE
  )
}

#' Print the cascade to the log, one line per step.
print_cascade <- function(cascade) {
  for (i in seq_len(nrow(cascade))) {
    cat(sprintf("  %-38s %8s genes  (%6.2f%% of previous, %6.2f%% of tested)\n",
                cascade$step_label[i], fmt_int(cascade$n_genes[i]),
                cascade$pct_of_previous[i], cascade$pct_of_tested[i]))
  }
}

#' Compare the size of the mCH response with the size of the MeCP2 response.
#'
#' Sufficiency is the share of mCH-significant genes that also change MeCP2.
#' Necessity is the share of MeCP2-changed genes that also change mCH. Both are
#' measured over the genes that carry a MeCP2 peak, where both quantities exist.
#'
#' @param gene_tbl Gene table from build_gene_table().
#' @param mecp2_db mecp2_diffbind from the shared config.
#' @param fdr_threshold FDR cutoff used for the peak-level counts.
#' @return data.frame with one row
build_scale_comparison <- function(gene_tbl, mecp2_db, fdr_threshold) {
  with_peak <- gene_tbl[gene_tbl$has_mecp2, , drop = FALSE]
  n_with_peak <- nrow(with_peak)
  n_mch_sig <- sum(with_peak$mch_sig)
  n_mecp2_sig <- sum(with_peak$mecp2_sig)
  n_both <- sum(with_peak$mch_sig & with_peak$mecp2_sig)

  out <- data.frame(
    n_genes_tested = nrow(gene_tbl),
    n_genes_with_mecp2_peak = n_with_peak,
    n_mch_sig_all_genes = sum(gene_tbl$mch_sig),
    n_mch_sig_with_peak = n_mch_sig,
    n_mecp2_sig_with_peak = n_mecp2_sig,
    n_both_sig = n_both,
    ratio_mch_to_mecp2 = n_mch_sig / n_mecp2_sig,
    pct_mch_sig_that_change_mecp2 = percent_of(n_both, n_mch_sig),
    pct_mecp2_changed_that_change_mch = percent_of(n_both, n_mecp2_sig),
    n_mecp2_peaks_total = nrow(mecp2_db),
    n_mecp2_peaks_gained = sum(mecp2_db$Fold > 0 & mecp2_db$FDR < fdr_threshold),
    n_mecp2_peaks_lost = sum(mecp2_db$Fold < 0 & mecp2_db$FDR < fdr_threshold),
    stringsAsFactors = FALSE
  )

  cat(sprintf("  Genes with a MeCP2 peak: %s\n", fmt_int(n_with_peak)))
  cat(sprintf("  mCH significant: %s | MeCP2 changed: %s | ratio %.1fx\n",
              fmt_int(n_mch_sig), fmt_int(n_mecp2_sig), out$ratio_mch_to_mecp2))
  cat(sprintf("  Of mCH-significant genes, %.2f%% also change MeCP2 (sufficiency)\n",
              out$pct_mch_sig_that_change_mecp2))
  cat(sprintf("  Of MeCP2-changed genes, %.2f%% also change mCH (necessity)\n",
              out$pct_mecp2_changed_that_change_mch))
  cat(sprintf("  MeCP2 peak level: %s peaks, %s gained, %s lost at FDR<%.2f\n",
              fmt_int(out$n_mecp2_peaks_total),
              fmt_int(out$n_mecp2_peaks_gained),
              fmt_int(out$n_mecp2_peaks_lost), fdr_threshold))
  out
}

# =============================================================================
# FIGURE 60_01a: TWO-PANEL HIERARCHY
# =============================================================================

#' One horizontal bar panel of a counting hierarchy.
#'
#' @param labels Character vector of step names, widest step first.
#' @param counts Integer vector of counts in the same order.
#' @param palette Fill colours in the same order.
#' @param denominator Count that the printed percentage divides by.
#' @param panel_title Title above the panel.
#' @param title_color Colour of the panel title.
#' @return ggplot object
hierarchy_panel <- function(labels, counts, palette, denominator,
                            panel_title, title_color) {
  df <- data.frame(
    step_label = factor(labels, levels = rev(labels)),
    n_genes = counts,
    stringsAsFactors = FALSE
  )

  ggplot(df, aes(x = step_label, y = n_genes, fill = step_label)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%s (%.1f%%)", fmt_int(n_genes),
                                  100 * n_genes / denominator)),
              hjust = -0.05, size = 3.4) +
    coord_flip(clip = "off") +
    scale_fill_manual(values = setNames(palette, labels), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(title = panel_title, x = "", y = "Number of genes") +
    theme_emseq() +
    theme(plot.title = element_text(color = title_color))
}

#' Draw the mCH hierarchy beside the MeCP2 hierarchy.
#'
#' @param gene_tbl Gene table from build_gene_table().
#' @param scale_tbl One-row table from build_scale_comparison().
#' @param out_dir Section output directory.
#' @return data.frame of both hierarchies
plot_hierarchy_panels <- function(gene_tbl, scale_tbl, out_dir) {
  n_tested <- nrow(gene_tbl)
  mch_labels <- c("Genes tested", "mCH significant",
                  "mCH hypermethylated", "mCH hypomethylated")
  mch_counts <- c(n_tested, sum(gene_tbl$mch_sig),
                  sum(gene_tbl$mch_hyper), sum(gene_tbl$mch_hypo))

  n_peak <- sum(gene_tbl$has_mecp2)
  mecp2_labels <- c("Genes with a MeCP2 peak", "MeCP2 significant",
                    "MeCP2 gained", "MeCP2 lost")
  mecp2_counts <- c(
    n_peak,
    sum(gene_tbl$mecp2_sig, na.rm = TRUE),
    sum(gene_tbl$mecp2_direction == "Gained", na.rm = TRUE),
    sum(gene_tbl$mecp2_direction == "Lost", na.rm = TRUE)
  )

  p_left <- hierarchy_panel(mch_labels, mch_counts, MCH_RAMP, n_tested,
                            "mCH hierarchy", MCH_PANEL_COLOR)
  p_right <- hierarchy_panel(mecp2_labels, mecp2_counts, MECP2_RAMP, n_peak,
                             "MeCP2 binding hierarchy", MECP2_PANEL_COLOR)

  # Both counts in the subtitle come from the genes that carry a MeCP2 peak, so
  # the printed ratio is the ratio of the two printed numbers.
  subtitle <- sprintf(
    paste("Among the %s genes with a MeCP2 peak: %s change mCH against %s that",
          "change MeCP2 (%.1fx). Of those mCH-significant genes, %.2f%% also change MeCP2."),
    fmt_int(scale_tbl$n_genes_with_mecp2_peak),
    fmt_int(scale_tbl$n_mch_sig_with_peak),
    fmt_int(scale_tbl$n_mecp2_sig_with_peak),
    scale_tbl$ratio_mch_to_mecp2,
    scale_tbl$pct_mch_sig_that_change_mecp2)

  combined <- (p_left | p_right) +
    plot_annotation(
      title = "Scale of mCH Change against Scale of MeCP2 Binding Change",
      subtitle = subtitle,
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30")
      )
    )

  save_multiformat_ggplot(combined,
                          file.path(out_dir, "60_01a_scale_hierarchy"),
                          width = 15, height = 7)

  rbind(
    data.frame(hierarchy = "mCH", step_label = mch_labels,
               n_genes = mch_counts,
               pct_of_hierarchy_top = 100 * mch_counts / n_tested,
               stringsAsFactors = FALSE),
    data.frame(hierarchy = "MeCP2", step_label = mecp2_labels,
               n_genes = mecp2_counts,
               pct_of_hierarchy_top = 100 * mecp2_counts / n_peak,
               stringsAsFactors = FALSE)
  )
}

# =============================================================================
# FIGURE 60_01b: TWO-SET VENN
# =============================================================================

#' Venn of mCH-significant genes against MeCP2-changed genes.
#'
#' Both sets sit inside the tested gene universe. Genes in the MeCP2 set but
#' not the mCH set are the handoff set that section 60_02 explains.
#'
#' @param gene_tbl Gene table from build_gene_table().
#' @param out_dir Section output directory.
#' @return data.frame of the four Venn regions
plot_venn <- function(gene_tbl, out_dir) {
  mch_genes <- gene_tbl$gene_name[gene_tbl$mch_sig]
  mecp2_genes <- gene_tbl$gene_name[which(gene_tbl$mecp2_sig)]

  n_both <- length(intersect(mch_genes, mecp2_genes))
  n_mch_only <- length(setdiff(mch_genes, mecp2_genes))
  n_mecp2_only <- length(setdiff(mecp2_genes, mch_genes))
  n_neither <- nrow(gene_tbl) - n_both - n_mch_only - n_mecp2_only

  venn_counts <- data.frame(
    region = c("mCH significant only", "MeCP2 changed only",
               "Both", "Neither", "Tested universe"),
    n_genes = c(n_mch_only, n_mecp2_only, n_both, n_neither, nrow(gene_tbl)),
    stringsAsFactors = FALSE
  )
  venn_counts$pct_of_tested <- 100 * venn_counts$n_genes / nrow(gene_tbl)

  for (i in seq_len(nrow(venn_counts))) {
    cat(sprintf("  %-22s %8s genes (%5.2f%% of tested)\n",
                venn_counts$region[i], fmt_int(venn_counts$n_genes[i]),
                venn_counts$pct_of_tested[i]))
  }

  venn_list <- list(mch_genes, mecp2_genes)
  names(venn_list) <- c(
    sprintf("mCH significant\n(n = %s)", fmt_int(length(mch_genes))),
    sprintf("MeCP2 changed\n(n = %s)", fmt_int(length(mecp2_genes)))
  )

  p <- ggVennDiagram(venn_list, label = "both", label_alpha = 0,
                     set_size = 4.2) +
    scale_fill_gradient(low = "white", high = MCH_PANEL_COLOR, name = "Genes") +
    scale_color_manual(values = c(MCH_PANEL_COLOR, MECP2_PANEL_COLOR)) +
    labs(
      title = "Gene-Level Overlap: mCH Change against MeCP2 Binding Change",
      subtitle = sprintf(
        paste("Universe: %s tested genes. %s genes change mCH without a MeCP2 change;",
              "%s change MeCP2 without an mCH change."),
        fmt_int(nrow(gene_tbl)), fmt_int(n_mch_only), fmt_int(n_mecp2_only))
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
      legend.position = "right"
    )

  save_multiformat_ggplot(p, file.path(out_dir, "60_01b_mch_mecp2_venn"),
                          width = 11, height = 8)
  venn_counts
}

# =============================================================================
# FIGURE 60_01c: QUADRANT SCATTER
# =============================================================================

#' Scatter of mCH change against MeCP2 fold change, split into four quadrants.
#'
#' Quadrants come from assign_quadrant(): Q1 both up, Q2 mCH down and MeCP2 up,
#' Q3 both down, Q4 mCH up and MeCP2 down. The count and percentage of each
#' quadrant sit in the matching corner.
#'
#' @param gene_tbl Gene table from build_gene_table().
#' @param out_dir Section output directory.
#' @return data.frame of quadrant counts
plot_quadrant_scatter <- function(gene_tbl, out_dir) {
  df <- gene_tbl[gene_tbl$has_mecp2, , drop = FALSE]
  df$quadrant <- assign_quadrant(df$mch_diff_pct, df$mecp2_fold)

  quad_tbl <- df %>%
    dplyr::count(quadrant, name = "n_genes") %>%
    dplyr::mutate(
      pct_of_genes = 100 * n_genes / sum(n_genes),
      meaning = unname(QUADRANT_MEANING[quadrant])
    ) %>%
    as.data.frame()

  for (i in seq_len(nrow(quad_tbl))) {
    cat(sprintf("  %s (%-22s): %8s genes (%5.2f%%)\n",
                quad_tbl$quadrant[i], quad_tbl$meaning[i],
                fmt_int(quad_tbl$n_genes[i]), quad_tbl$pct_of_genes[i]))
  }

  quad_n <- setNames(quad_tbl$n_genes, quad_tbl$quadrant)
  quad_pct <- setNames(quad_tbl$pct_of_genes, quad_tbl$quadrant)
  corner_label <- function(q) {
    n <- if (q %in% names(quad_n)) quad_n[[q]] else 0L
    pct <- if (q %in% names(quad_pct)) quad_pct[[q]] else 0
    sprintf("%s\n%s\nn = %s (%.1f%%)", q, QUADRANT_MEANING[[q]], fmt_int(n), pct)
  }

  x_max <- max(df$mch_diff_pct)
  x_min <- min(df$mch_diff_pct)
  y_max <- max(df$mecp2_fold)
  y_min <- min(df$mecp2_fold)

  df$label_gene <- ifelse(df$gene_name %in% KEY_GENES, df$gene_name, "")

  p <- ggplot(df, aes(x = mch_diff_pct, y = mecp2_fold)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40",
               linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40",
               linewidth = 0.4) +
    geom_point(aes(color = quadrant), alpha = 0.35, size = 1.2) +
    geom_text_repel(aes(label = label_gene), size = 3, max.overlaps = 20,
                    fontface = "italic", color = "grey15",
                    segment.color = "grey60", segment.size = 0.3,
                    min.segment.length = 0) +
    annotate("text", x = 0.72 * x_max, y = 0.90 * y_max,
             label = corner_label("Q1"),
             size = 3.2, fontface = "bold", color = COLORS$quadrant[["Q1"]]) +
    annotate("text", x = 0.72 * x_min, y = 0.90 * y_max,
             label = corner_label("Q2"),
             size = 3.2, fontface = "bold", color = COLORS$quadrant[["Q2"]]) +
    annotate("text", x = 0.72 * x_min, y = 0.90 * y_min,
             label = corner_label("Q3"),
             size = 3.2, fontface = "bold", color = COLORS$quadrant[["Q3"]]) +
    annotate("text", x = 0.72 * x_max, y = 0.90 * y_min,
             label = corner_label("Q4"),
             size = 3.2, fontface = "bold", color = COLORS$quadrant[["Q4"]]) +
    scale_color_manual(values = COLORS$quadrant, name = "Quadrant") +
    labs(
      title = "mCH Change against MeCP2 Fold Change",
      subtitle = sprintf(
        "All %s genes with a MeCP2 peak. MeCP2 fold comes from the peak nearest the TSS.",
        fmt_int(nrow(df))),
      x = "mCH change (mutant - control, percentage points)",
      y = "MeCP2 log2 fold change (mutant / control)"
    ) +
    theme_emseq() +
    theme(legend.position = "right")

  save_multiformat_ggplot(p, file.path(out_dir, "60_01c_mch_mecp2_quadrant"),
                          width = 12, height = 9)
  quad_tbl
}

# =============================================================================
# FIGURE 60_01d: CASCADE FUNNEL
# =============================================================================

#' Bar funnel of the five cascade steps.
#'
#' Every bar carries its count, the percentage of the tested universe, and the
#' percentage of the step before it.
#'
#' @param cascade data.frame from build_cascade_table().
#' @param out_dir Section output directory.
#' @return the ggplot object
plot_cascade_funnel <- function(cascade, out_dir) {
  df <- cascade
  df$short_label <- factor(df$short_label, levels = df$short_label)

  final_pct <- df$pct_of_tested[nrow(df)]

  p <- ggplot(df, aes(x = short_label, y = n_genes, fill = short_label)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%s\n%.2f%% of tested\n%.1f%% of previous",
                                  fmt_int(n_genes), pct_of_tested,
                                  pct_of_previous)),
              vjust = -0.25, size = 3.2, lineheight = 1.1) +
    scale_fill_manual(values = setNames(FUNNEL_RAMP, levels(df$short_label)),
                      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(
      title = "From mCH Change to MeCP2 Binding Change",
      subtitle = sprintf(
        paste("Each step is a subset of the step before it.",
              "%.2f%% of tested genes reach the final step (%s genes)."),
        final_pct, fmt_int(df$n_genes[nrow(df)])),
      x = "", y = "Number of genes"
    ) +
    theme_emseq() +
    theme(axis.text.x = element_text(size = 9))

  save_multiformat_ggplot(p, file.path(out_dir, "60_01d_cascade_funnel"),
                          width = 11, height = 8)
  p
}

# =============================================================================
# FIGURE 60_01e: MeCP2 FOLD BY mCH DIRECTION
# =============================================================================

#' Violin of MeCP2 fold change across the three mCH direction groups.
#'
#' Each group carries its gene count and median on the figure, and the pairwise
#' Wilcoxon results are annotated above the violins.
#'
#' @param gene_tbl Gene table from build_gene_table().
#' @param out_dir Section output directory.
#' @return list with the group summary and the Wilcoxon table
plot_fold_by_mch_direction <- function(gene_tbl, out_dir) {
  df <- gene_tbl[gene_tbl$has_mecp2, , drop = FALSE]

  grp_summary <- summarise_groups(df, "mch_category", "mecp2_fold")
  for (i in seq_len(nrow(grp_summary))) {
    cat(sprintf("  %-16s n = %8s, median MeCP2 fold = %+.4f, mean = %+.4f\n",
                as.character(grp_summary$mch_category[i]),
                fmt_int(grp_summary$n[i]),
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

  # Figure text only. grp_summary itself stays data-only, because main() writes it.
  grp_labels <- grp_summary
  grp_labels$label <- group_label(grp_summary)

  p <- ggplot(df, aes(x = mch_category, y = mecp2_fold, fill = mch_category)) +
    geom_violin(alpha = 0.6, scale = "width",
                draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black",
               linewidth = 0.4) +
    geom_text(data = grp_labels,
              aes(x = mch_category, y = y_label, label = label),
              inherit.aes = FALSE, size = 3.2, lineheight = 1.1) +
    annotate("text", x = 2, y = y_stats, label = wilcox_label,
             size = 3.2, fontface = "italic", lineheight = 1.2) +
    scale_fill_manual(values = MCH_CATEGORY_COLORS, guide = "none") +
    scale_y_continuous(limits = c(y_min - 0.05 * y_span,
                                  y_stats + 0.14 * y_span)) +
    labs(
      title = "MeCP2 Fold Change by mCH Direction",
      subtitle = sprintf(
        "Genes with a MeCP2 peak (n = %s). Fold comes from the peak nearest the TSS.",
        fmt_int(nrow(df))),
      x = "mCH direction (edgeR FDR < 0.05)",
      y = "MeCP2 log2 fold change (mutant / control)"
    ) +
    theme_emseq()

  save_multiformat_ggplot(p,
                          file.path(out_dir, "60_01e_mecp2_fold_by_mch_direction"),
                          width = 10, height = 8)
  list(summary = grp_summary, wilcoxon = wilcox_tbl)
}

# =============================================================================
# ANALYSIS: REGISTERED FISHER TEST
# =============================================================================

#' Register the mCH significance against MeCP2 significance Fisher test.
#'
#' Runs over the genes that carry a MeCP2 peak, where both flags are defined.
#'
#' @param gene_tbl Gene table from build_gene_table().
#' @param out_dir Section output directory.
#' @return data.frame with one row summarising the test
run_significance_fisher <- function(gene_tbl, out_dir) {
  df <- gene_tbl[gene_tbl$has_mecp2, , drop = FALSE]
  df$mecp2_sig <- as.logical(df$mecp2_sig)

  ft <- register_fisher_test(
    section = SECTION_ID, test_id = "mch_sig_vs_mecp2_sig",
    description = paste("Among genes with a MeCP2 peak, do genes with a",
                        "significant mCH change also change MeCP2 binding",
                        "more often than other genes?"),
    gene_df = df, row_var = "mch_sig", col_var = "mecp2_sig",
    output_dir = out_dir)

  data.frame(
    test_id = "mch_sig_vs_mecp2_sig",
    row_var = "mch_sig",
    col_var = "mecp2_sig",
    n_genes = nrow(df),
    n_mch_sig = sum(df$mch_sig),
    n_mecp2_sig = sum(df$mecp2_sig),
    n_both = sum(df$mch_sig & df$mecp2_sig),
    n_mch_only = sum(df$mch_sig & !df$mecp2_sig),
    n_mecp2_only = sum(!df$mch_sig & df$mecp2_sig),
    n_neither = sum(!df$mch_sig & !df$mecp2_sig),
    odds_ratio = unname(ft$estimate),
    ci_low = ft$conf.int[1],
    ci_high = ft$conf.int[2],
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# HANDOFF TABLE
# =============================================================================

#' Write the genes where MeCP2 moves without an mCH change.
#'
#' Keeps genes that carry a MeCP2 peak, whose MeCP2 binding changes at the
#' threshold, and whose mCH does not change. Section 60_02 reads this file and
#' indexes HANDOFF_COLUMNS by name.
#'
#' @param gene_tbl Gene table from build_gene_table().
#' @param out_dir Section output directory.
#' @return data.frame written to disk
write_handoff_table <- function(gene_tbl, out_dir) {
  keep <- gene_tbl$has_mecp2 & gene_tbl$mecp2_sig %in% TRUE & !gene_tbl$mch_sig
  handoff <- gene_tbl[keep, , drop = FALSE]

  missing <- setdiff(HANDOFF_COLUMNS, colnames(handoff))
  if (length(missing) > 0) {
    stop("Gene table is missing handoff columns: ",
         paste(missing, collapse = ", "))
  }
  handoff <- handoff[, HANDOFF_COLUMNS, drop = FALSE]

  if (nrow(handoff) == 0) {
    stop("No gene changes MeCP2 binding without an mCH change. Section 60_02 ",
         "needs a non-empty mecp2_no_meth_genes.tsv.")
  }

  path <- file.path(out_dir, basename(HANDOFF_PATHS$mecp2_no_meth_genes))
  write_section_table(handoff, path)

  n_gained <- sum(handoff$mecp2_direction == "Gained")
  n_lost <- sum(handoff$mecp2_direction == "Lost")
  cat(sprintf("  MeCP2 changes without an mCH change: %s genes (%s gained, %s lost)\n",
              fmt_int(nrow(handoff)), fmt_int(n_gained), fmt_int(n_lost)))

  if (normalizePath(out_dir, mustWork = FALSE) !=
      normalizePath(OUTPUT_PATHS$mecp2, mustWork = FALSE)) {
    cat(sprintf("  NOTE: section 60_02 reads %s\n",
                HANDOFF_PATHS$mecp2_no_meth_genes))
  }
  handoff
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
  cat("SECTION 60_01: mCH SCALE vs MeCP2 BINDING SCALE\n")
  cat("================================================================================\n")
  cat("Output dir:    ", OUT_DIR, "\n", sep = "")
  cat("FDR threshold: ", fdr_threshold, "\n\n", sep = "")

  stopifnot(
    "MeCP2 annotated peak file not found" = file.exists(MECP2_PATHS$annotated)
  )

  cat("--- Loading MeCP2 peak annotation ---\n")
  mecp2_annotated <- load_mecp2_annotated(MECP2_PATHS$annotated, fdr_threshold)

  cat("\n--- Aggregating MeCP2 peaks to genes (nearest TSS) ---\n")
  mecp2_gene <- aggregate_mecp2_by_gene(mecp2_annotated, fdr_threshold)

  cat("\n--- Building the gene table ---\n")
  keep_idx <- deduplicate_mch_row_indices(mch_results)
  cat(sprintf("  Deduplicated %s mCH rows to %s gene names\n",
              fmt_int(nrow(mch_results)), fmt_int(length(keep_idx))))
  gene_tbl <- build_gene_table(mch_results[keep_idx, , drop = FALSE],
                               mecp2_gene, fdr_threshold)
  write_section_table(gene_tbl, file.path(OUT_DIR, "60_01_gene_level_mch_mecp2.tsv"))

  cat("\n--- Analysis 1: the cascade ---\n")
  cascade <- build_cascade_table(gene_tbl)
  print_cascade(cascade)
  write_section_table(cascade, file.path(OUT_DIR, "60_01_cascade.tsv"))

  cat("\n--- Analysis 2: scale comparison ---\n")
  scale_tbl <- build_scale_comparison(gene_tbl, mecp2_diffbind, fdr_threshold)
  write_section_table(scale_tbl, file.path(OUT_DIR, "60_01_scale_comparison.tsv"))

  cat("\n--- Figure 60_01a: two-panel hierarchy ---\n")
  hierarchy_tbl <- plot_hierarchy_panels(gene_tbl, scale_tbl, OUT_DIR)
  write_section_table(hierarchy_tbl, file.path(OUT_DIR, "60_01_hierarchy_counts.tsv"))

  cat("\n--- Figure 60_01b: two-set Venn ---\n")
  venn_tbl <- plot_venn(gene_tbl, OUT_DIR)
  write_section_table(venn_tbl, file.path(OUT_DIR, "60_01_venn_counts.tsv"))

  cat("\n--- Figure 60_01c: quadrant scatter ---\n")
  quad_tbl <- plot_quadrant_scatter(gene_tbl, OUT_DIR)
  write_section_table(quad_tbl, file.path(OUT_DIR, "60_01_quadrant_counts.tsv"))

  cat("\n--- Figure 60_01d: cascade funnel ---\n")
  plot_cascade_funnel(cascade, OUT_DIR)

  cat("\n--- Figure 60_01e: MeCP2 fold by mCH direction ---\n")
  fold_res <- plot_fold_by_mch_direction(gene_tbl, OUT_DIR)
  write_section_table(fold_res$summary,
                      file.path(OUT_DIR,
                                "60_01_mecp2_fold_by_mch_direction_summary.tsv"))
  write_section_table(fold_res$wilcoxon,
                      file.path(OUT_DIR,
                                "60_01_mecp2_fold_by_mch_direction_wilcoxon.tsv"))

  cat("\n--- Analysis 3: registered Fisher test ---\n")
  fisher_tbl <- run_significance_fisher(gene_tbl, OUT_DIR)
  write_section_table(fisher_tbl, file.path(OUT_DIR, "60_01_fisher_summary.tsv"))

  cat("\n--- Handoff: MeCP2 changes without an mCH change ---\n")
  handoff <- write_handoff_table(gene_tbl, OUT_DIR)

  # The handoff set is the MeCP2-only region of the Venn, so the two counts
  # must agree.
  n_mecp2_only <- venn_tbl$n_genes[venn_tbl$region == "MeCP2 changed only"]
  stopifnot(
    "Handoff row count does not match the MeCP2-only Venn region" =
      nrow(handoff) == n_mecp2_only
  )

  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 60_01 SUMMARY\n")
  cat("================================================================================\n")
  for (i in seq_len(nrow(cascade))) {
    cat(sprintf("%-38s %8s  (%6.2f%% of previous, %6.2f%% of tested)\n",
                cascade$step_label[i], fmt_int(cascade$n_genes[i]),
                cascade$pct_of_previous[i], cascade$pct_of_tested[i]))
  }
  cat(sprintf("Genes with a MeCP2 peak:               %s\n",
              fmt_int(scale_tbl$n_genes_with_mecp2_peak)))
  cat(sprintf("mCH-significant to MeCP2-changed ratio: %.1fx\n",
              scale_tbl$ratio_mch_to_mecp2))
  cat(sprintf("mCH-significant genes that change MeCP2: %.2f%%\n",
              scale_tbl$pct_mch_sig_that_change_mecp2))
  cat(sprintf("MeCP2-changed genes that change mCH:     %.2f%%\n",
              scale_tbl$pct_mecp2_changed_that_change_mch))
  cat(sprintf("Fisher mch_sig x mecp2_sig:            OR = %.3f (95%% CI %.3f to %.3f), p = %.3g\n",
              fisher_tbl$odds_ratio[1], fisher_tbl$ci_low[1],
              fisher_tbl$ci_high[1], fisher_tbl$p_value[1]))
  cat(sprintf("Wilcoxon hyper vs hypo MeCP2 fold:     p = %.3g\n",
              fold_res$wilcoxon$p_value[
                fold_res$wilcoxon$group_1 == "Hypermethylated" &
                  fold_res$wilcoxon$group_2 == "Hypomethylated"]))
  cat(sprintf("MeCP2 change without an mCH change:    %s genes\n",
              fmt_int(nrow(handoff))))
  cat(sprintf("Handoff table: %s\n",
              file.path(OUT_DIR, basename(HANDOFF_PATHS$mecp2_no_meth_genes))))
  cat(sprintf("Fisher tests registered in: %s\n", HANDOFF_PATHS$fisher_registry))
  cat("\nSection 60_01 complete.\n\n")
}

main()
