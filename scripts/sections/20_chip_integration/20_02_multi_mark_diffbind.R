# scripts/sections/20_chip_integration/20_02_multi_mark_diffbind.R
#
# Section 20_02: multi-mark differential binding integrated with gene-body mCH.
#
# Tests whether the four chromatin marks that BAP1 loss perturbs (ATAC,
# H3K27ac, H3K27me3, H2AK119ub) move together with the gene-body mCH change.
# Each mark enters alone (correlation matrix, direction enrichment, scatter)
# and all four enter together (multivariate logistic model of mCH
# hypermethylation, per-gene convergence count).
#
# Reads:
#   mch_results, DIFFBIND_TABLES   pre-loaded by _shared_config.R
#   HANDOFF_PATHS$chromatin_state  promoter_state and body_state per gene,
#                                  written by 10_01, joined on gene_id
#
# Writes into --output-dir (default results/sections/20_chip_integration/):
#   gene_level_all_marks.tsv   gene-level master table read by sections 20_03,
#                              20_04, 70_02 and 70_04
#   20_02a .. 20_02h           figures, each with the TSV holding its numbers
#   fisher_tables/             gene tables behind the registered Fisher tests
#
# Adapted from Biomodal section 33 (multi-mark DiffBind integration).

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)
library(pROC)

# =============================================================================
# SECTION CONSTANTS
# =============================================================================

SECTION_ID <- "20_02"

MARK_ORDER <- c("atac", "k27ac", "k27me3", "k119ub")

# display  name used in figures and tables
# color    panel accent, taken from the shared COLORS palettes
# gained   label for peaks with Fold > 0 at the section FDR threshold
# lost     label for peaks with Fold < 0 at the section FDR threshold
# short    label used in the convergence combination strings
# concordant  fold-change sign counted by the convergence score
MARK_META <- list(
  atac = list(
    display = "ATAC-seq", color = COLORS$atac[["ATAC Up"]],
    gained = "ATAC Up", lost = "ATAC Down", short = "ATAC", concordant = "down"
  ),
  k27ac = list(
    display = "H3K27ac", color = COLORS$h3k27ac[["H3K27ac Gained"]],
    gained = "K27ac Gained", lost = "K27ac Lost", short = "K27ac", concordant = "down"
  ),
  k27me3 = list(
    display = "H3K27me3", color = COLORS$h3k27me3[["H3K27me3 Gained"]],
    gained = "K27me3 Gained", lost = "K27me3 Lost", short = "K27me3", concordant = "up"
  ),
  k119ub = list(
    display = "H2AK119ub", color = COLORS$k119ub[["K119ub Gained"]],
    gained = "K119ub Gained", lost = "K119ub Lost", short = "K119ub", concordant = "up"
  )
)

PEAK_DIRECTION_COLORS <- c("Gained" = "#B2182B", "Lost" = "#2166AC",
                           "Unchanged" = "grey70")

MCH_GROUP_ORDER <- c("mCH Hyper", "mCH Hypo", "Not Significant")

MCH_GROUP_COLORS <- c(
  "mCH Hyper"       = COLORS$direction[["Hypermethylated"]],
  "mCH Hypo"        = COLORS$direction[["Hypomethylated"]],
  "Not Significant" = "grey70"
)

MCH_DIRECTION_COLORS <- c(
  "mCH Hyper" = COLORS$direction[["Hypermethylated"]],
  "mCH Hypo"  = COLORS$direction[["Hypomethylated"]]
)

CORRELATION_VARS <- c(
  mch_diff    = "mCH difference",
  atac_fold   = "ATAC-seq log2FC",
  k27ac_fold  = "H3K27ac log2FC",
  k27me3_fold = "H3K27me3 log2FC",
  k119ub_fold = "H2AK119ub log2FC"
)

# Column contract for gene_level_all_marks.tsv. Sections 20_03, 20_04, 70_02
# and 70_04 read this file and index these columns by name.
#
# Chromatin state arrives as two columns, promoter_state over the TSS window and
# body_state over the gene body. The marks this section integrates are gene-body
# marks, so the stratifications below read body_state.
HANDOFF_COLUMNS <- c(
  "gene_name", "gene_id", "chr", "start", "end", "gene_length",
  "mch_ctrl", "mch_mut", "mch_diff", "edger_logFC", "edger_fdr",
  "mch_sig", "mch_direction",
  "promoter_state", "body_state",
  "atac_fold", "atac_fdr", "atac_n_peaks",
  "k27ac_fold", "k27ac_fdr", "k27ac_n_peaks",
  "k27me3_fold", "k27me3_fdr", "k27me3_n_peaks",
  "k119ub_fold", "k119ub_fdr", "k119ub_n_peaks"
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
                help = "FDR cutoff for peak significance [default: %default]")
  )
  opt <- parse_args(OptionParser(option_list = option_list))

  if (opt$fdr_threshold <= 0 || opt$fdr_threshold >= 1) {
    stop("--fdr-threshold must be between 0 and 1, got ", opt$fdr_threshold)
  }
  opt
}

# =============================================================================
# TABLE AND FIGURE OUTPUT
# =============================================================================

significance_stars <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", ""))))
}

# =============================================================================
# STEP 1: DIFFBIND PEAKS TO GENE-LEVEL FOLD CHANGES
# =============================================================================

#' Restrict a DiffBind table to canonical chromosomes and re-classify direction
#' at the section FDR threshold.
prepare_diffbind_table <- function(db, mark, fdr_threshold) {
  meta <- MARK_META[[mark]]

  n_all <- nrow(db)
  db <- db[db$Chr %in% CANONICAL_CHRS, , drop = FALSE]
  db <- db[!is.na(db$Fold) & !is.na(db$FDR), , drop = FALSE]
  if (nrow(db) == 0) {
    stop("No usable ", meta$display, " peaks after filtering to canonical ",
         "chromosomes and complete Fold/FDR values.")
  }

  db$direction <- "Unchanged"
  db$direction[db$FDR < fdr_threshold & db$Fold > 0] <- "Gained"
  db$direction[db$FDR < fdr_threshold & db$Fold < 0] <- "Lost"

  cat(sprintf("  %-10s %s peaks in, %s kept (%s gained, %s lost at FDR<%.3g)\n",
              meta$display,
              format(n_all, big.mark = ","), format(nrow(db), big.mark = ","),
              format(sum(db$direction == "Gained"), big.mark = ","),
              format(sum(db$direction == "Lost"), big.mark = ","),
              fdr_threshold))
  db
}

prepare_all_diffbind_tables <- function(fdr_threshold) {
  cat("--- Step 1: preparing DiffBind tables ---\n")
  tables <- list()
  for (mark in MARK_ORDER) {
    tables[[mark]] <- prepare_diffbind_table(DIFFBIND_TABLES[[mark]], mark,
                                             fdr_threshold)
  }
  cat("\n")
  tables
}

#' Annotate one DiffBind table to genes and collapse to one row per gene.
#'
#' The peak nearest the TSS supplies the gene-level fold change and FDR.
annotate_mark_to_genes <- function(db, mark, fdr_threshold) {
  meta <- MARK_META[[mark]]
  annotated <- annotate_peaks_to_genes(db, meta$display)

  usable <- !is.na(annotated$SYMBOL) & nzchar(annotated$SYMBOL) &
    !is.na(annotated$Fold) & !is.na(annotated$FDR)
  cat(sprintf("    %s of %s annotated peaks carry a gene symbol\n",
              format(sum(usable), big.mark = ","),
              format(nrow(annotated), big.mark = ",")))
  annotated <- annotated[usable, , drop = FALSE]

  gene_table <- aggregate_diffbind_by_gene(
    annotated, method = "nearest_tss",
    fdr_threshold = fdr_threshold, prefix = mark
  )
  if (anyDuplicated(gene_table$gene_name) > 0) {
    stop("aggregate_diffbind_by_gene() returned duplicate genes for ", meta$display)
  }
  cat(sprintf("    %s genes carry a %s peak\n",
              format(nrow(gene_table), big.mark = ","), meta$display))
  gene_table
}

annotate_all_marks <- function(diffbind_tables, fdr_threshold) {
  cat("--- Step 2: annotating peaks to genes ---\n")
  gene_tables <- list()
  for (mark in MARK_ORDER) {
    gene_tables[[mark]] <- annotate_mark_to_genes(diffbind_tables[[mark]], mark,
                                                  fdr_threshold)
  }
  cat("\n")
  gene_tables
}

# =============================================================================
# STEP 3: GENE-LEVEL MASTER TABLE
# =============================================================================

#' Read the chromatin state table written by section 10_01.
#'
#' The join key is gene_id. Gene names are not unique across ENSMUSG
#' identifiers, so a gene-name join would attach one gene's state to another
#' gene's row.
#'
#' Rows that repeat a gene_id with the same pair of states collapse; a gene_id
#' carrying two different pairs stops the script. A state value outside
#' PROMOTER_STATE_ORDER or BODY_STATE_ORDER also stops the script.
load_chromatin_state <- function(path) {
  if (!file.exists(path)) {
    stop("Chromatin state table not found: ", path,
         "\nRun section 10_01 (10_01_chromatin_state.R) first. It writes this file.")
  }

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")

  state_cols <- c("gene_id", "promoter_state", "body_state")
  missing <- setdiff(state_cols, colnames(df))
  if (length(missing) > 0) {
    stop("Chromatin state table is missing columns: ",
         paste(missing, collapse = ", "), " (", path, ")",
         "\nSection 10_01 (10_01_chromatin_state.R) writes gene_id, ",
         "promoter_state and body_state. Rerun it.")
  }

  unknown_promoter <- setdiff(unique(df$promoter_state), PROMOTER_STATE_ORDER)
  if (length(unknown_promoter) > 0) {
    stop("Unknown promoter_state value(s) in ", path, ": ",
         paste(unknown_promoter, collapse = ", "),
         "\nAllowed: ", paste(PROMOTER_STATE_ORDER, collapse = ", "))
  }
  unknown_body <- setdiff(unique(df$body_state), BODY_STATE_ORDER)
  if (length(unknown_body) > 0) {
    stop("Unknown body_state value(s) in ", path, ": ",
         paste(unknown_body, collapse = ", "),
         "\nAllowed: ", paste(BODY_STATE_ORDER, collapse = ", "))
  }

  df <- unique(df[, state_cols])
  conflicting <- unique(df$gene_id[duplicated(df$gene_id)])
  if (length(conflicting) > 0) {
    stop(length(conflicting), " gene identifiers carry more than one state pair in ",
         path, ". First conflicts: ",
         paste(head(conflicting, 5), collapse = ", "),
         "\nSection 10_01 must write one promoter_state and one body_state per gene.")
  }

  cat(sprintf("  Chromatin state: %s genes, joined on gene_id\n",
              format(nrow(df), big.mark = ",")))
  df
}

#' Join the four gene-level mark tables and the chromatin state onto every
#' gene in mch_results.
build_master_table <- function(gene_tables, state_table) {
  cat("--- Step 3: building the gene-level master table ---\n")

  master <- mch_results
  for (mark in MARK_ORDER) {
    master <- dplyr::left_join(master, gene_tables[[mark]], by = "gene_name")
  }

  master <- dplyr::left_join(master, state_table, by = "gene_id")

  if (nrow(master) != nrow(mch_results)) {
    stop("Joins changed the row count: ", nrow(mch_results), " genes in, ",
         nrow(master), " out. A joined table has duplicate keys.")
  }

  n_dup_names <- sum(duplicated(master$gene_name))
  cat(sprintf("  Master table: %s genes (%d share a gene name with another row)\n",
              format(nrow(master), big.mark = ","), n_dup_names))

  for (mark in MARK_ORDER) {
    n_with <- sum(!is.na(master[[paste0(mark, "_fold")]]))
    cat(sprintf("    %-10s %s genes with a peak (%.1f%%)\n",
                MARK_META[[mark]]$display, format(n_with, big.mark = ","),
                100 * n_with / nrow(master)))
  }

  n_promoter <- sum(!is.na(master$promoter_state))
  n_body <- sum(!is.na(master$body_state))
  cat(sprintf("    %-10s %s genes with a promoter state (%.1f%%)\n",
              "state", format(n_promoter, big.mark = ","),
              100 * n_promoter / nrow(master)))
  cat(sprintf("    %-10s %s genes with a body state (%.1f%%)\n",
              "", format(n_body, big.mark = ","),
              100 * n_body / nrow(master)))

  # Print both state distributions. An empty state is a defect in 10_01 that
  # would otherwise only show up as an empty group several sections later.
  for (state in PROMOTER_STATE_ORDER) {
    n <- sum(!is.na(master$promoter_state) & master$promoter_state == state)
    cat(sprintf("      promoter %-19s %s genes (%.1f%%)\n", state,
                format(n, big.mark = ","), 100 * n / nrow(master)))
  }
  for (state in BODY_STATE_ORDER) {
    n <- sum(!is.na(master$body_state) & master$body_state == state)
    cat(sprintf("      body     %-19s %s genes (%.1f%%)\n", state,
                format(n, big.mark = ","), 100 * n / nrow(master)))
  }

  n_complete <- sum(stats::complete.cases(master[, paste0(MARK_ORDER, "_fold")]))
  cat(sprintf("    %s genes carry all four marks (%.1f%%)\n\n",
              format(n_complete, big.mark = ","),
              100 * n_complete / nrow(master)))

  master
}

#' Write the cross-section handoff table with its exact column contract.
write_handoff_table <- function(master, out_dir) {
  missing <- setdiff(HANDOFF_COLUMNS, colnames(master))
  if (length(missing) > 0) {
    stop("Master table is missing handoff columns: ",
         paste(missing, collapse = ", "))
  }

  path <- file.path(out_dir, basename(HANDOFF_PATHS$gene_level_all_marks))
  write_section_table(master[, HANDOFF_COLUMNS], path)

  if (normalizePath(out_dir, mustWork = FALSE) !=
      normalizePath(OUTPUT_PATHS$chip, mustWork = FALSE)) {
    cat(sprintf("  NOTE: sections 20_03, 20_04, 70_02 and 70_04 read %s\n",
                HANDOFF_PATHS$gene_level_all_marks))
  }
  path
}

# =============================================================================
# FIGURE 20_02a: PER-MARK VOLCANO PLOTS
# =============================================================================

build_volcano_panel <- function(db, mark, fdr_threshold) {
  meta <- MARK_META[[mark]]
  df <- db
  df$neg_log10_fdr <- -log10(pmax(df$FDR, 1e-300))
  df$direction <- factor(df$direction, levels = names(PEAK_DIRECTION_COLORS))

  counts <- sprintf("Gained: %s\nLost: %s\nUnchanged: %s",
                    format(sum(df$direction == "Gained"), big.mark = ","),
                    format(sum(df$direction == "Lost"), big.mark = ","),
                    format(sum(df$direction == "Unchanged"), big.mark = ","))

  ggplot(df, aes(x = Fold, y = neg_log10_fdr, color = direction)) +
    geom_point(size = 0.5, alpha = 0.4) +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed",
               color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60",
               linewidth = 0.3) +
    scale_color_manual(values = PEAK_DIRECTION_COLORS, name = "") +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2,
             label = counts, size = 3.2, color = "grey20") +
    labs(title = meta$display, x = "DiffBind log2FC (mut / ctrl)",
         y = "-log10(FDR)") +
    theme_emseq() +
    theme(legend.position = "none",
          plot.title = element_text(color = meta$color, hjust = 0.5,
                                    face = "bold"))
}

plot_mark_volcanoes <- function(diffbind_tables, fdr_threshold, out_dir) {
  cat("--- Figure 20_02a: per-mark volcano plots ---\n")

  panels <- lapply(MARK_ORDER, function(mark) {
    build_volcano_panel(diffbind_tables[[mark]], mark, fdr_threshold)
  })
  names(panels) <- MARK_ORDER

  figure <- (panels$atac + panels$k27ac) / (panels$k27me3 + panels$k119ub) +
    plot_annotation(
      title = "Differential binding of four chromatin marks",
      subtitle = sprintf("Dashed line marks FDR = %.3g", fdr_threshold),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(figure, file.path(out_dir, "20_02a_mark_volcanoes"),
                          width = 14, height = 10)

  summary_df <- do.call(rbind, lapply(MARK_ORDER, function(mark) {
    db <- diffbind_tables[[mark]]
    data.frame(
      mark = MARK_META[[mark]]$display,
      n_peaks = nrow(db),
      n_gained = sum(db$direction == "Gained"),
      n_lost = sum(db$direction == "Lost"),
      n_unchanged = sum(db$direction == "Unchanged"),
      median_fold = median(db$Fold),
      min_fdr = min(db$FDR),
      fdr_threshold = fdr_threshold,
      stringsAsFactors = FALSE
    )
  }))
  write_section_table(summary_df, file.path(out_dir, "20_02a_peak_summary.tsv"))
  cat("\n")
  invisible(summary_df)
}

# =============================================================================
# FIGURE 20_02b: CROSS-MARK SPEARMAN CORRELATION MATRIX
# =============================================================================

#' Pairwise-complete Spearman correlations across mch_diff and the four folds.
compute_correlation_matrices <- function(master) {
  vars <- names(CORRELATION_VARS)
  labels <- unname(CORRELATION_VARS)
  n_vars <- length(vars)

  rho <- matrix(NA_real_, n_vars, n_vars, dimnames = list(labels, labels))
  pval <- matrix(NA_real_, n_vars, n_vars, dimnames = list(labels, labels))
  npair <- matrix(NA_integer_, n_vars, n_vars, dimnames = list(labels, labels))

  for (i in seq_len(n_vars)) {
    for (j in seq_len(n_vars)) {
      x <- master[[vars[i]]]
      y <- master[[vars[j]]]
      keep <- !is.na(x) & !is.na(y)
      npair[i, j] <- sum(keep)

      if (i == j) {
        rho[i, j] <- 1
        next
      }
      if (sum(keep) < 10) {
        stop("Only ", sum(keep), " genes carry both ", vars[i], " and ", vars[j],
             ". The gene-level join produced too little overlap to correlate.")
      }
      test <- cor.test(x[keep], y[keep], method = "spearman", exact = FALSE)
      rho[i, j] <- unname(test$estimate)
      pval[i, j] <- test$p.value
    }
  }
  list(rho = rho, pval = pval, n = npair, labels = labels)
}

plot_correlation_heatmap <- function(cors, out_dir) {
  cat("--- Figure 20_02b: cross-mark Spearman correlation matrix ---\n")

  display <- matrix(
    sprintf("%.2f%s", cors$rho, significance_stars(cors$pval)),
    nrow = nrow(cors$rho), dimnames = dimnames(cors$rho)
  )
  diag(display) <- ""

  hm_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  hm_breaks <- seq(-1, 1, length.out = 101)
  hm_title <- "Spearman correlation: mCH difference and mark log2FC (pairwise complete)"

  # bquote splices the matrices into the call, because
  # save_multiformat_pheatmap() evaluates it outside this function.
  hm_call <- bquote(pheatmap(
    mat = .(cors$rho),
    display_numbers = .(display),
    color = .(hm_colors),
    breaks = .(hm_breaks),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize = 12,
    fontsize_number = 11,
    number_color = "black",
    border_color = "grey80",
    na_col = "grey90",
    main = .(hm_title),
    silent = TRUE
  ))

  save_multiformat_pheatmap(hm_call,
                            file.path(out_dir, "20_02b_cross_mark_correlation"),
                            width = 10, height = 9)

  labels <- cors$labels
  long <- expand.grid(variable_1 = labels, variable_2 = labels,
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  long$spearman_rho <- as.vector(cors$rho)
  long$p_value <- as.vector(cors$pval)
  long$n_genes <- as.vector(cors$n)
  write_section_table(long, file.path(out_dir, "20_02b_cross_mark_correlation.tsv"))

  wide <- as.data.frame(cors$rho)
  wide <- cbind(variable = rownames(cors$rho), wide)
  colnames(wide) <- c("variable", labels)
  write_section_table(wide, file.path(out_dir, "20_02b_correlation_matrix_wide.tsv"))
  cat("\n")
  invisible(long)
}

# =============================================================================
# FIGURE 20_02c: mCH DIFFERENCE VERSUS MARK FOLD CHANGE
# =============================================================================

build_scatter_panel <- function(master, mark) {
  meta <- MARK_META[[mark]]
  fold_col <- paste0(mark, "_fold")

  df <- master[!is.na(master[[fold_col]]) & !is.na(master$mch_diff), , drop = FALSE]
  df$sig_group <- ifelse(df$mch_sig, "Significant", "Not Significant")
  df$sig_group <- factor(df$sig_group, levels = names(COLORS$significant))

  test <- cor.test(df[[fold_col]], df$mch_diff, method = "spearman",
                   exact = FALSE)
  key_df <- df[df$gene_name %in% KEY_GENES, , drop = FALSE]

  panel <- ggplot(df, aes(x = .data[[fold_col]], y = mch_diff)) +
    geom_point(aes(color = sig_group), size = 0.4, alpha = 0.35) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60",
               linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60",
               linewidth = 0.3) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "black",
                linewidth = 0.8) +
    scale_color_manual(values = COLORS$significant, name = "mCH") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.3,
             label = sprintf("rho = %.3f\np = %.2g\nn = %s",
                             unname(test$estimate), test$p.value,
                             format(nrow(df), big.mark = ",")),
             size = 3.2, color = "grey20") +
    labs(title = meta$display, x = "DiffBind log2FC (mut / ctrl)",
         y = "mCH difference (mut - ctrl)") +
    theme_emseq() +
    theme(legend.position = "none",
          plot.title = element_text(color = meta$color, hjust = 0.5,
                                    face = "bold"))

  if (nrow(key_df) > 0) {
    panel <- panel +
      geom_point(data = key_df, size = 1.6, color = "black") +
      ggrepel::geom_text_repel(data = key_df, aes(label = gene_name),
                               size = 2.8, min.segment.length = 0,
                               max.overlaps = Inf, color = "black")
  }

  stats_row <- data.frame(
    mark = meta$display,
    n_genes = nrow(df),
    spearman_rho = unname(test$estimate),
    p_value = test$p.value,
    n_key_genes_labelled = nrow(key_df),
    stringsAsFactors = FALSE
  )
  list(panel = panel, stats = stats_row)
}

plot_mch_vs_mark_scatters <- function(master, out_dir) {
  cat("--- Figure 20_02c: mCH difference versus mark log2FC ---\n")

  built <- lapply(MARK_ORDER, function(mark) build_scatter_panel(master, mark))
  panels <- lapply(built, `[[`, "panel")

  figure <- wrap_plots(panels, nrow = 1) +
    plot_annotation(
      title = "Gene-body mCH change against chromatin mark fold change",
      subtitle = "Loess fit in black, key genes labelled",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "20_02c_mch_vs_mark_scatters"),
                          width = 18, height = 5.5)

  stats_df <- do.call(rbind, lapply(built, `[[`, "stats"))
  write_section_table(stats_df, file.path(out_dir, "20_02c_mch_vs_mark_correlations.tsv"))
  cat("\n")
  invisible(stats_df)
}

# =============================================================================
# FIGURE 20_02d: DIRECTION ENRICHMENT (OBSERVED / EXPECTED)
# =============================================================================

#' Test one mark's direction against mCH direction and return the four cells
#' of the 2x2 table with observed, expected and O/E values.
test_mark_direction <- function(master, mark, fdr_threshold, out_dir) {
  meta <- MARK_META[[mark]]
  fold_col <- paste0(mark, "_fold")
  fdr_col <- paste0(mark, "_fdr")
  gained_col <- paste0(mark, "_gained")

  eligible <- master$mch_sig &
    !is.na(master[[fold_col]]) & !is.na(master[[fdr_col]]) &
    master[[fdr_col]] < fdr_threshold
  df <- master[eligible, , drop = FALSE]

  if (nrow(df) < 10) {
    stop("Only ", nrow(df), " mCH-significant genes carry a significant ",
         meta$display, " peak. Too few for the direction test.")
  }

  gene_df <- data.frame(
    gene_name = df$gene_name,
    chr = df$chr,
    mch_hyper = df$edger_logFC > 0,
    stringsAsFactors = FALSE
  )
  gene_df[[gained_col]] <- df[[fold_col]] > 0

  tab <- table(
    factor(ifelse(gene_df$mch_hyper, "mCH Hyper", "mCH Hypo"),
           levels = c("mCH Hyper", "mCH Hypo")),
    factor(ifelse(gene_df[[gained_col]], meta$gained, meta$lost),
           levels = c(meta$gained, meta$lost))
  )
  if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) {
    stop("The 2x2 table for ", meta$display, " has an empty row or column: ",
         paste(as.vector(tab), collapse = ", "))
  }

  ft <- register_fisher_test(
    section = SECTION_ID,
    test_id = paste0(mark, "_direction"),
    description = sprintf(paste("Among mCH-significant genes with a significant",
                                "%s peak, mCH hypermethylation against %s gain."),
                          meta$display, meta$display),
    gene_df = gene_df, row_var = "mch_hyper", col_var = gained_col,
    output_dir = out_dir
  )

  expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)

  cells <- as.data.frame(as.table(tab), stringsAsFactors = FALSE)
  colnames(cells) <- c("mch_direction", "mark_direction", "observed")
  expected_df <- as.data.frame(as.table(expected), stringsAsFactors = FALSE)
  colnames(expected_df) <- c("mch_direction", "mark_direction", "expected")

  cells <- merge(cells, expected_df, by = c("mch_direction", "mark_direction"))
  cells$oe_ratio <- cells$observed / cells$expected
  cells$mark <- meta$display
  cells$n_genes <- sum(tab)
  cells$fisher_or <- unname(ft$estimate)
  cells$fisher_p <- ft$p.value
  cells$fisher_ci_lower <- ft$conf.int[1]
  cells$fisher_ci_upper <- ft$conf.int[2]
  cells[, c("mark", "mch_direction", "mark_direction", "observed", "expected",
            "oe_ratio", "n_genes", "fisher_or", "fisher_ci_lower",
            "fisher_ci_upper", "fisher_p")]
}

plot_direction_enrichment <- function(master, fdr_threshold, out_dir) {
  cat("--- Figure 20_02d: mCH direction against mark direction ---\n")

  cells <- do.call(rbind, lapply(MARK_ORDER, function(mark) {
    test_mark_direction(master, mark, fdr_threshold, out_dir)
  }))

  mark_levels <- vapply(MARK_META[MARK_ORDER], `[[`, character(1), "display")
  cells$facet_label <- sprintf("%s\nOR = %.2f, p = %.2g (n = %s)",
                               cells$mark, cells$fisher_or, cells$fisher_p,
                               format(cells$n_genes, big.mark = ",", trim = TRUE))
  facet_levels <- unique(cells$facet_label[order(match(cells$mark, mark_levels))])
  cells$facet_label <- factor(cells$facet_label, levels = facet_levels)
  cells$mch_direction <- factor(cells$mch_direction,
                                levels = names(MCH_DIRECTION_COLORS))

  direction_levels <- unlist(lapply(MARK_ORDER, function(mark) {
    c(MARK_META[[mark]]$gained, MARK_META[[mark]]$lost)
  }))
  cells$mark_direction <- factor(cells$mark_direction, levels = direction_levels)

  figure <- ggplot(cells, aes(x = mark_direction, y = oe_ratio,
                              color = mch_direction)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    geom_point(size = 4, position = position_dodge(width = 0.6)) +
    geom_text(aes(label = sprintf("%.2f", oe_ratio)),
              position = position_dodge(width = 0.6), vjust = -1.1, size = 3.2,
              show.legend = FALSE) +
    geom_text(aes(label = sprintf("obs %s / exp %.0f",
                                  format(observed, big.mark = ",", trim = TRUE),
                                  expected)),
              position = position_dodge(width = 0.6), vjust = 2.2, size = 2.6,
              color = "grey35", show.legend = FALSE) +
    scale_color_manual(values = MCH_DIRECTION_COLORS, name = "mCH direction") +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.2))) +
    facet_wrap(~ facet_label, nrow = 1, scales = "free_x") +
    labs(
      title = "Observed / expected enrichment of mark direction by mCH direction",
      subtitle = sprintf(paste("mCH-significant genes with a peak at FDR < %.3g;",
                               "dashed line is the null (O/E = 1)"),
                         fdr_threshold),
      x = "", y = "Observed / expected"
    ) +
    theme_emseq() +
    theme(legend.position = "top",
          panel.spacing = unit(1.2, "lines"),
          strip.text = element_text(size = 9))

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "20_02d_direction_enrichment"),
                          width = 16, height = 7)
  write_section_table(cells[, setdiff(colnames(cells), "facet_label")],
                      file.path(out_dir, "20_02d_direction_enrichment.tsv"))
  cat("\n")
  invisible(cells)
}

# =============================================================================
# FIGURES 20_02e AND 20_02f: MULTIVARIATE LOGISTIC MODEL
# =============================================================================

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

fit_hypermethylation_model <- function(master) {
  fold_cols <- paste0(MARK_ORDER, "_fold")
  df <- master[stats::complete.cases(master[, fold_cols]), , drop = FALSE]
  df$mch_hyper_outcome <- as.integer(df$mch_hyper)

  if (length(unique(df$mch_hyper_outcome)) < 2) {
    stop("The mCH hypermethylation outcome has one class only among the ",
         nrow(df), " genes carrying all four marks.")
  }

  fit <- glm(mch_hyper_outcome ~ atac_fold + k27ac_fold + k27me3_fold + k119ub_fold,
             data = df, family = binomial)
  null_fit <- glm(mch_hyper_outcome ~ 1, data = df, family = binomial)
  mcfadden <- 1 - as.numeric(logLik(fit)) / as.numeric(logLik(null_fit))

  probs <- predict(fit, type = "response")
  roc_obj <- pROC::roc(df$mch_hyper_outcome, probs, quiet = TRUE)
  auc_ci <- as.numeric(pROC::ci.auc(roc_obj))

  cat(sprintf("  %s genes, %s hypermethylated (%.1f%%)\n",
              format(nrow(df), big.mark = ","),
              format(sum(df$mch_hyper_outcome), big.mark = ","),
              100 * mean(df$mch_hyper_outcome)))
  cat(sprintf("  AUC = %.3f [%.3f, %.3f], McFadden R2 = %.4f\n",
              auc_ci[2], auc_ci[1], auc_ci[3], mcfadden))

  list(fit = fit, data = df, roc = roc_obj,
       auc = auc_ci[2], auc_lower = auc_ci[1], auc_upper = auc_ci[3],
       mcfadden = mcfadden, n_genes = nrow(df),
       n_hyper = sum(df$mch_hyper_outcome))
}

plot_logistic_forest <- function(model, out_dir) {
  or_df <- extract_odds_ratios(model$fit)
  or_df <- or_df[or_df$term != "(Intercept)", , drop = FALSE]

  term_display <- vapply(MARK_META[MARK_ORDER], `[[`, character(1), "display")
  names(term_display) <- paste0(MARK_ORDER, "_fold")
  or_df$mark <- unname(term_display[or_df$term])
  or_df$stars <- significance_stars(or_df$p_value)
  or_df$stars[or_df$stars == ""] <- "ns"
  or_df$mark <- factor(or_df$mark, levels = rev(unname(term_display)))

  figure <- ggplot(or_df, aes(x = odds_ratio, y = mark)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    geom_errorbar(aes(xmin = or_lower, xmax = or_upper), orientation = "y",
                  width = 0.2, linewidth = 0.8) +
    geom_point(size = 3.5, color = COLORS$direction[["Hypermethylated"]]) +
    geom_text(aes(label = sprintf("OR = %.2f [%.2f, %.2f] %s", odds_ratio,
                                  or_lower, or_upper, stars)),
              vjust = -1.1, size = 3.2) +
    scale_x_log10() +
    scale_y_discrete(expand = expansion(add = 0.6)) +
    labs(
      title = "Four-mark logistic model of mCH hypermethylation",
      subtitle = sprintf(paste("n = %s genes with all four marks, %s",
                               "hypermethylated | AUC = %.3f [%.3f, %.3f] |",
                               "McFadden R2 = %.4f"),
                         format(model$n_genes, big.mark = ","),
                         format(model$n_hyper, big.mark = ","),
                         model$auc, model$auc_lower, model$auc_upper,
                         model$mcfadden),
      x = "Odds ratio per unit log2FC (log scale)", y = ""
    ) +
    theme_emseq() +
    theme(axis.text.y = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"))

  save_multiformat_ggplot(figure, file.path(out_dir, "20_02e_logistic_forest"),
                          width = 11, height = 6)
  write_section_table(or_df, file.path(out_dir, "20_02e_logistic_coefficients.tsv"))
  invisible(or_df)
}

plot_logistic_roc <- function(model, out_dir) {
  roc_df <- data.frame(
    false_positive_rate = 1 - as.numeric(model$roc$specificities),
    true_positive_rate = as.numeric(model$roc$sensitivities),
    threshold = as.numeric(model$roc$thresholds)
  )
  roc_df <- roc_df[order(roc_df$false_positive_rate, roc_df$true_positive_rate), ]

  figure <- ggplot(roc_df, aes(x = false_positive_rate, y = true_positive_rate)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_line(color = COLORS$direction[["Hypermethylated"]], linewidth = 1) +
    annotate("text", x = 0.65, y = 0.15,
             label = sprintf("AUC = %.3f\n95%% CI [%.3f, %.3f]\nn = %s genes",
                             model$auc, model$auc_lower, model$auc_upper,
                             format(model$n_genes, big.mark = ",")),
             size = 4, color = "grey20") +
    coord_equal() +
    labs(
      title = "ROC curve: four-mark model of mCH hypermethylation",
      subtitle = "Predictors are the ATAC, H3K27ac, H3K27me3 and H2AK119ub log2FC",
      x = "False positive rate", y = "True positive rate"
    ) +
    theme_emseq()

  save_multiformat_ggplot(figure, file.path(out_dir, "20_02f_logistic_roc"),
                          width = 7, height = 7)
  write_section_table(roc_df, file.path(out_dir, "20_02f_roc_curve.tsv"))

  fit_df <- data.frame(
    n_genes = model$n_genes,
    n_hyper = model$n_hyper,
    pct_hyper = 100 * model$n_hyper / model$n_genes,
    auc = model$auc,
    auc_ci_lower = model$auc_lower,
    auc_ci_upper = model$auc_upper,
    mcfadden_r2 = model$mcfadden,
    aic = AIC(model$fit),
    stringsAsFactors = FALSE
  )
  write_section_table(fit_df, file.path(out_dir, "20_02f_logistic_model_fit.tsv"))
  invisible(fit_df)
}

run_logistic_model <- function(master, out_dir) {
  cat("--- Figures 20_02e and 20_02f: multivariate logistic model ---\n")
  model <- fit_hypermethylation_model(master)
  plot_logistic_forest(model, out_dir)
  plot_logistic_roc(model, out_dir)
  cat("\n")
  invisible(model)
}

# =============================================================================
# FIGURES 20_02g AND 20_02h: CONVERGENCE COUNT
# =============================================================================

#' Score each gene for the number of marks that pass the FDR threshold in the
#' direction named by MARK_META$concordant.
add_convergence_scores <- function(master, fdr_threshold) {
  for (mark in MARK_ORDER) {
    fold <- master[[paste0(mark, "_fold")]]
    fdr <- master[[paste0(mark, "_fdr")]]
    passes <- !is.na(fold) & !is.na(fdr) & fdr < fdr_threshold
    if (MARK_META[[mark]]$concordant == "down") {
      master[[paste0(mark, "_concordant")]] <- passes & fold < 0
    } else {
      master[[paste0(mark, "_concordant")]] <- passes & fold > 0
    }
    master[[paste0(mark, "_tested")]] <- !is.na(fold)
  }

  master$n_concordant <- rowSums(master[, paste0(MARK_ORDER, "_concordant")])
  master$n_marks_tested <- rowSums(master[, paste0(MARK_ORDER, "_tested")])

  master$mch_group <- ifelse(master$mch_hyper, "mCH Hyper",
                             ifelse(master$mch_hypo, "mCH Hypo", "Not Significant"))
  master$mch_group <- factor(master$mch_group, levels = MCH_GROUP_ORDER)

  short_names <- vapply(MARK_META[MARK_ORDER], `[[`, character(1), "short")
  flags <- as.matrix(master[, paste0(MARK_ORDER, "_concordant")])
  master$mark_combination <- apply(flags, 1, function(row) {
    paste(short_names[row], collapse = " + ")
  })
  master$mark_combination[master$mark_combination == ""] <- "none"

  master
}

plot_convergence_by_group <- function(master, fdr_threshold, out_dir) {
  cat("--- Figure 20_02g: convergence count by mCH group ---\n")

  group_totals <- master %>%
    dplyr::count(mch_group, name = "n_group")

  stacked <- master %>%
    dplyr::count(mch_group, n_concordant, name = "n") %>%
    dplyr::left_join(group_totals, by = "mch_group") %>%
    dplyr::mutate(pct = 100 * n / n_group,
                  n_concordant = factor(n_concordant,
                                        levels = sort(unique(master$n_concordant))))

  axis_labels <- setNames(
    sprintf("%s\n(n = %s)", as.character(group_totals$mch_group),
            format(group_totals$n_group, big.mark = ",", trim = TRUE)),
    as.character(group_totals$mch_group)
  )

  p_stacked <- ggplot(stacked, aes(x = mch_group, y = pct, fill = n_concordant)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = ifelse(pct >= 3, sprintf("%.1f%%", pct), "")),
              position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_brewer(palette = "YlOrRd", name = "Concordant\nmarks") +
    scale_x_discrete(labels = axis_labels) +
    labs(
      title = "Concordant mark count by mCH group",
      subtitle = sprintf(paste("Concordant: ATAC or H3K27ac log2FC < 0, or",
                               "H3K27me3 or H2AK119ub log2FC > 0, at FDR < %.3g"),
                         fdr_threshold),
      x = "", y = "Percent of genes"
    ) +
    theme_emseq() +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"))

  group_stats <- summarise_groups(master, "mch_group", "n_concordant")
  group_stats$mch_group <- factor(group_stats$mch_group, levels = MCH_GROUP_ORDER)

  max_concordant <- max(master$n_concordant)
  group_stats$label_y <- max_concordant + 0.9
  group_stats$label <- group_label(group_stats, digits = 2)

  kruskal <- kruskal.test(n_concordant ~ mch_group, data = master)

  p_box <- ggplot(master, aes(x = mch_group, y = n_concordant, fill = mch_group)) +
    geom_boxplot(width = 0.55, outlier.size = 0.4, outlier.alpha = 0.3) +
    geom_text(data = group_stats, aes(x = mch_group, y = label_y, label = label),
              inherit.aes = FALSE, size = 3.2, lineheight = 1.1) +
    scale_fill_manual(values = MCH_GROUP_COLORS, name = "") +
    scale_x_discrete(labels = axis_labels) +
    scale_y_continuous(limits = c(-0.4, max_concordant + 1.4),
                       breaks = sort(unique(master$n_concordant))) +
    labs(
      title = "Concordant mark count per gene",
      subtitle = sprintf("Kruskal-Wallis chi-squared = %.1f, df = %d, p = %.3g",
                         unname(kruskal$statistic), unname(kruskal$parameter),
                         kruskal$p.value),
      x = "", y = "Concordant marks per gene"
    ) +
    theme_emseq() +
    theme(legend.position = "none",
          plot.subtitle = element_text(size = 9, color = "grey40"))

  figure <- p_stacked + p_box +
    plot_annotation(
      title = "Epigenomic convergence at BAP1-responsive genes",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )

  save_multiformat_ggplot(figure, file.path(out_dir, "20_02g_convergence_by_group"),
                          width = 15, height = 7)

  stacked_out <- as.data.frame(stacked)
  stacked_out$n_concordant <- as.integer(as.character(stacked_out$n_concordant))
  write_section_table(stacked_out,
                      file.path(out_dir, "20_02g_convergence_distribution.tsv"))

  summary_out <- group_stats[, c("mch_group", "n", "median", "mean")]
  colnames(summary_out) <- c("mch_group", "n_genes", "median_concordant",
                             "mean_concordant")
  summary_out$kruskal_chisq <- unname(kruskal$statistic)
  summary_out$kruskal_df <- unname(kruskal$parameter)
  summary_out$kruskal_p <- kruskal$p.value
  write_section_table(summary_out,
                      file.path(out_dir, "20_02g_convergence_summary.tsv"))
  cat("\n")
  invisible(summary_out)
}

plot_mark_combinations <- function(master, out_dir) {
  cat("--- Figure 20_02h: concordant mark combinations ---\n")

  multi <- master[master$n_concordant >= 2, , drop = FALSE]
  if (nrow(multi) == 0) {
    stop("No gene carries two or more concordant marks. ",
         "Check the DiffBind joins before interpreting the convergence figures.")
  }

  combos <- multi %>%
    dplyr::count(mark_combination, name = "n_genes") %>%
    dplyr::arrange(dplyr::desc(n_genes))
  top_combos <- head(combos, 15)
  top_combos$mark_combination <- factor(top_combos$mark_combination,
                                        levels = rev(top_combos$mark_combination))

  figure <- ggplot(top_combos, aes(x = mark_combination, y = n_genes)) +
    geom_col(fill = COLORS$k119ub[["K119ub Gained"]], width = 0.7) +
    geom_text(aes(label = format(n_genes, big.mark = ",", trim = TRUE)),
              hjust = -0.15, size = 3.4) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Mark combinations at genes with two or more concordant marks",
      subtitle = sprintf("n = %s genes, %d distinct combinations, top %d shown",
                         format(nrow(multi), big.mark = ","), nrow(combos),
                         nrow(top_combos)),
      x = "Mark combination", y = "Genes"
    ) +
    theme_emseq()

  save_multiformat_ggplot(figure, file.path(out_dir, "20_02h_mark_combinations"),
                          width = 11, height = 7)
  write_section_table(as.data.frame(combos),
                      file.path(out_dir, "20_02h_mark_combinations.tsv"))
  cat("\n")
  invisible(combos)
}

#' Fisher test asking whether convergent genes are the hypermethylated ones.
test_convergence_direction <- function(master, out_dir) {
  df <- master[master$mch_sig & master$n_marks_tested > 0, , drop = FALSE]
  if (nrow(df) < 10) {
    stop("Only ", nrow(df), " mCH-significant genes carry any mark peak. ",
         "Too few for the convergence direction test.")
  }

  gene_df <- data.frame(
    gene_name = df$gene_name,
    chr = df$chr,
    mch_hyper = df$edger_logFC > 0,
    convergent_2plus = df$n_concordant >= 2,
    stringsAsFactors = FALSE
  )

  register_fisher_test(
    section = SECTION_ID,
    test_id = "convergence_2plus",
    description = paste("Among mCH-significant genes with at least one mark peak,",
                        "mCH hypermethylation against two or more concordant marks."),
    gene_df = gene_df, row_var = "mch_hyper", col_var = "convergent_2plus",
    output_dir = out_dir
  )
}

write_convergence_table <- function(master, out_dir) {
  cols <- c("gene_name", "gene_id", "chr", "gene_length",
            "promoter_state", "body_state",
            "mch_diff", "edger_logFC", "edger_fdr", "mch_sig", "mch_group",
            unlist(lapply(MARK_ORDER, function(mark) {
              paste0(mark, c("_fold", "_fdr", "_concordant"))
            })),
            "n_marks_tested", "n_concordant", "mark_combination")
  write_section_table(master[, cols],
                      file.path(out_dir, "20_02_convergence_per_gene.tsv"))
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
  cat("  SECTION 20_02: MULTI-MARK DIFFBIND INTEGRATION WITH mCH\n")
  cat("===========================================================================\n")
  cat("  Output directory: ", out_dir, "\n", sep = "")
  cat("  Peak FDR threshold: ", fdr_threshold, "\n\n", sep = "")

  diffbind_tables <- prepare_all_diffbind_tables(fdr_threshold)
  gene_tables <- annotate_all_marks(diffbind_tables, fdr_threshold)
  state_table <- load_chromatin_state(HANDOFF_PATHS$chromatin_state)

  master <- build_master_table(gene_tables, state_table)
  handoff_path <- write_handoff_table(master, out_dir)
  cat("\n")

  master <- add_convergence_scores(master, fdr_threshold)

  plot_mark_volcanoes(diffbind_tables, fdr_threshold, out_dir)
  cors <- compute_correlation_matrices(master)
  plot_correlation_heatmap(cors, out_dir)
  plot_mch_vs_mark_scatters(master, out_dir)
  enrichment <- plot_direction_enrichment(master, fdr_threshold, out_dir)
  model <- run_logistic_model(master, out_dir)
  plot_convergence_by_group(master, fdr_threshold, out_dir)
  plot_mark_combinations(master, out_dir)
  test_convergence_direction(master, out_dir)
  write_convergence_table(master, out_dir)

  cat("\n---------------------------------------------------------------------------\n")
  cat("  SUMMARY\n")
  cat("---------------------------------------------------------------------------\n")
  cat(sprintf("  Handoff table: %s\n", handoff_path))
  cat(sprintf("  Genes: %s | all four marks: %s\n",
              format(nrow(master), big.mark = ","),
              format(model$n_genes, big.mark = ",")))
  cat(sprintf("  Logistic model AUC: %.3f [%.3f, %.3f], McFadden R2 = %.4f\n",
              model$auc, model$auc_lower, model$auc_upper, model$mcfadden))
  cat("  Spearman rho against mCH difference:\n")
  for (mark in MARK_ORDER) {
    label <- CORRELATION_VARS[[paste0(mark, "_fold")]]
    cat(sprintf("    %-18s %+.3f\n", label,
                cors$rho["mCH difference", label]))
  }
  cat("  Fisher odds ratios (mCH hyper against mark gain):\n")
  for (mark in MARK_ORDER) {
    display <- MARK_META[[mark]]$display
    row <- enrichment[enrichment$mark == display, ][1, ]
    cat(sprintf("    %-10s OR = %.3f, p = %.3g (n = %s genes)\n",
                display, row$fisher_or, row$fisher_p,
                format(row$n_genes, big.mark = ",")))
  }
  cat(sprintf("  Genes with 3 or more concordant marks: %s\n",
              format(sum(master$n_concordant >= 3), big.mark = ",")))
  cat("\n=== Section 20_02 complete ===\n\n")
}

main()
