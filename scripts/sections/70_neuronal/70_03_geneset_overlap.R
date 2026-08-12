# scripts/sections/70_neuronal/70_03_geneset_overlap.R
#
# Section 70_03: neuronal gene sets against MeCP2 gain, and mCH inside them.
#
# Asks two questions about the same gene universe.
#   1. Do genes that gain MeCP2 binding in the mutant overlap the neuronal gene
#      class more than chance predicts? The section answers this for a broad
#      neuronal set (section 70_01) and for a narrower synapse and axon set.
#   2. Does the mCH level, and the mCH change between genotypes, differ between
#      those gene sets?
#
# Gene length confounds both mCH and MeCP2 binding, so the section also
# compares gene length between every set.
#
# Reads:
#   HANDOFF_PATHS$neuronal_gene_set   broad neuronal gene set, written by 70_01
#   GENESET_PATHS$synapse             narrower synapse and axon gene set
#   mecp2_diffbind                    MeCP2 differential binding, shared config
#   mch_results                       gene-level mCH results, shared config
#
# Writes into --output-dir (default results/sections/70_neuronal/):
#   70_03a_geneset_venn        two-set Venn of neuronal against MeCP2 up
#   70_03b_mch_level_by_set    mCH level per gene set
#   70_03c_mch_diff_by_set     mCH change per gene set
#   70_03d_fisher_forest       odds ratios of the two registered Fisher tests
#   70_03e_gene_length_by_set  gene length per gene set
#   70_03f_composite           all panels on one page
#   *.tsv beside each figure holding every number the figure shows
#   fisher_tables/             gene tables behind the registered Fisher tests
#
# Adapted from Biomodal section 74 (gene set overlap and methylation levels).
# EM-seq measures non-CG methylation as one modality, so every panel below
# reads mCH alone and the Venn compares two gene sets.

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

SECTION_ID <- "70_03"

# The four sets partition the universe. Every gene falls in exactly one.
SET_ORDER <- c("Neuronal + MeCP2 up", "Neuronal only", "MeCP2 up only", "Other")

SET_COLORS <- c(
  "Neuronal + MeCP2 up" = "#7570B3",
  "Neuronal only"       = "#1B9E77",
  "MeCP2 up only"       = unname(COLORS$mecp2[["MeCP2 Up"]]),
  "Other"               = "grey70"
)

# Columns the 70_01 handoff must carry for this section to read it.
NEURONAL_HANDOFF_COLUMNS <- c("gene", "source")

# Fewest genes a MeCP2 gene-level table may hold before the input is treated
# as broken rather than sparse.
MIN_MECP2_GENES <- 1000

# Fewest values a Wilcoxon group may hold.
MIN_WILCOXON_GROUP <- 3

# Fill colour of the Venn circles, dark end of the count gradient.
VENN_FILL_HIGH <- "#7570B3"

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
                help = "FDR cutoff for a MeCP2 peak [default: %default]")
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

format_count <- function(x) format(x, big.mark = ",", trim = TRUE)

#' Extend summarise_groups() with the range of each group.
#'
#' summarise_groups() gives n, median, mean and the quartiles. The minimum and
#' the maximum let a reader see the full spread behind a violin from the table
#' alone. The result is data only, so it is writable as TSV. The on-plot text
#' comes from group_label() at each plot site.
set_summary <- function(df, group_col, value_col) {
  base <- summarise_groups(df, group_col, value_col)

  spread <- df %>%
    dplyr::filter(!is.na(.data[[value_col]])) %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(
      min = min(.data[[value_col]]),
      max = max(.data[[value_col]]),
      .groups = "drop"
    ) %>%
    as.data.frame()
  names(spread)[1] <- group_col

  dplyr::left_join(base, spread, by = group_col)
}

# =============================================================================
# STEP 1: GENE SETS
# =============================================================================

#' Read the broad neuronal gene set written by section 70_01.
load_neuronal_handoff <- function() {
  cat("--- Step 1a: neuronal gene set from section 70_01 ---\n")

  path <- HANDOFF_PATHS$neuronal_gene_set
  if (!file.exists(path)) {
    stop("Neuronal gene set handoff not found: ", path,
         "\nSection 70_01 (scripts/sections/70_neuronal/70_01_k119ub_neuronal.R) ",
         "writes this file. Run section 70_01 before this section.")
  }

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")
  missing <- setdiff(NEURONAL_HANDOFF_COLUMNS, colnames(df))
  if (length(missing) > 0) {
    stop("The section 70_01 neuronal handoff is missing columns: ",
         paste(missing, collapse = ", "), " in ", path)
  }

  genes <- unique(df$gene[!is.na(df$gene) & nzchar(df$gene)])
  if (length(genes) == 0) {
    stop("The section 70_01 neuronal handoff holds no gene symbol: ", path)
  }

  source_counts <- table(df$source)
  cat(sprintf("  %s genes in the handoff, %s unique symbols\n",
              format_count(nrow(df)), format_count(length(genes))))
  for (src in names(source_counts)) {
    cat(sprintf("    source %-10s %s\n", src,
                format_count(as.integer(source_counts[[src]]))))
  }
  cat("\n")
  genes
}

#' Collapse the MeCP2 DiffBind peakset to one fold change per gene.
#'
#' The nearest-TSS rule keeps the peak closest to each gene's TSS, so one gene
#' carries one MeCP2 fold change and one FDR.
aggregate_mecp2_to_genes <- function(fdr_threshold) {
  cat("--- Step 1c: MeCP2 differential binding aggregated to genes ---\n")

  annotated <- annotate_peaks_to_genes(mecp2_diffbind, "MeCP2")
  required <- c("SYMBOL", "Fold", "FDR", "distanceToTSS")
  missing <- setdiff(required, colnames(annotated))
  if (length(missing) > 0) {
    stop("ChIPseeker returned no ", paste(missing, collapse = ", "),
         " column for the MeCP2 peaks, so the nearest-TSS collapse cannot run.")
  }

  gene_tbl <- aggregate_diffbind_by_gene(annotated, method = "nearest_tss",
                                         fdr_threshold = fdr_threshold,
                                         prefix = "mecp2")
  if (nrow(gene_tbl) < MIN_MECP2_GENES) {
    stop("Only ", nrow(gene_tbl), " genes carry a MeCP2 peak. ",
         "Expected at least ", MIN_MECP2_GENES,
         " from ", nrow(mecp2_diffbind), " peaks.")
  }

  n_up <- sum(gene_tbl$mecp2_fdr < fdr_threshold & gene_tbl$mecp2_fold > 0,
              na.rm = TRUE)
  n_down <- sum(gene_tbl$mecp2_fdr < fdr_threshold & gene_tbl$mecp2_fold < 0,
                na.rm = TRUE)
  cat(sprintf("  %s genes carry a MeCP2 peak: %s up, %s down at FDR < %.2f\n\n",
              format_count(nrow(gene_tbl)), format_count(n_up),
              format_count(n_down), fdr_threshold))
  gene_tbl
}

# =============================================================================
# STEP 2: THE GENE UNIVERSE
# =============================================================================

#' Row indices of mch_results that keep one row per gene symbol.
#'
#' Some gene names carry more than one ENSMUSG identifier. The row with the
#' largest absolute edgeR log fold change is kept for each name, so the
#' universe holds one row per symbol and gene set membership joins cleanly.
deduplicate_mch_row_indices <- function(mch) {
  ord <- order(mch$gene_name, -abs(mch$edger_logFC))
  keep <- ord[!duplicated(mch$gene_name[ord])]
  sort(keep)
}

#' Build the gene universe and flag set membership on it.
#'
#' The universe is every gene in mch_results, one row per gene symbol. That
#' choice fixes what "chance" means for the Fisher tests below: a gene enters
#' the 2x2 table only if it carries an mCH measurement and an edgeR test, so
#' the odds ratio compares neuronal genes against other tested genes and never
#' against genes this experiment could not measure.
#'
#' A universe gene with no MeCP2 peak is not MeCP2 up. is_mecp2_up is FALSE for
#' it rather than NA, which keeps the Fisher universe equal to the whole tested
#' gene set. has_mecp2_peak records how much of the universe the MeCP2 peakset
#' reaches.
build_universe <- function(neuronal_genes, synapse_genes, mecp2_gene_tbl,
                           fdr_threshold) {
  cat("--- Step 2: gene universe ---\n")

  keep <- deduplicate_mch_row_indices(mch_results)
  cat(sprintf("  mch_results: %s rows, %s after one row per gene symbol\n",
              format_count(nrow(mch_results)), format_count(length(keep))))

  universe <- mch_results[keep, c("gene_name", "gene_id", "chr", "start", "end",
                                  "gene_length", "mch_ctrl", "mch_mut",
                                  "mch_diff", "edger_logFC", "edger_fdr",
                                  "mch_sig", "mch_direction")]
  rownames(universe) <- NULL

  universe <- dplyr::left_join(
    universe,
    mecp2_gene_tbl[, c("gene_name", "mecp2_fold", "mecp2_fdr", "mecp2_n_peaks")],
    by = "gene_name"
  )
  if (nrow(universe) != length(keep)) {
    stop("The MeCP2 join changed the universe row count: ", length(keep),
         " genes in, ", nrow(universe), " out. The MeCP2 gene table repeats a ",
         "gene symbol.")
  }

  universe$has_mecp2_peak <- !is.na(universe$mecp2_fdr)
  no_fold <- universe$has_mecp2_peak & is.na(universe$mecp2_fold)
  if (any(no_fold)) {
    stop(sum(no_fold), " genes carry a MeCP2 FDR with no fold change, so the ",
         "direction of their MeCP2 change cannot be assigned.")
  }
  universe$is_mecp2_up <- universe$has_mecp2_peak &
    universe$mecp2_fdr < fdr_threshold & universe$mecp2_fold > 0
  universe$is_neuronal <- universe$gene_name %in% neuronal_genes
  universe$is_synapse  <- universe$gene_name %in% synapse_genes

  universe$gene_set <- factor(
    ifelse(universe$is_neuronal & universe$is_mecp2_up, "Neuronal + MeCP2 up",
      ifelse(universe$is_neuronal, "Neuronal only",
        ifelse(universe$is_mecp2_up, "MeCP2 up only", "Other"))),
    levels = SET_ORDER
  )

  empty_sets <- setdiff(SET_ORDER, unique(as.character(universe$gene_set)))
  if (length(empty_sets) > 0) {
    stop("These gene sets hold no gene, so the per-set figures cannot be ",
         "drawn: ", paste(empty_sets, collapse = ", "))
  }
  if (anyNA(universe$gene_length) || any(universe$gene_length <= 0)) {
    stop("gene_length is missing or not positive for ",
         sum(is.na(universe$gene_length) | universe$gene_length <= 0),
         " universe genes. The gene-length figure uses a log10 axis.")
  }
  for (col in c("mch_ctrl", "mch_mut", "mch_diff")) {
    n_bad <- sum(!is.finite(universe[[col]]))
    if (n_bad > 0) {
      stop(n_bad, " universe genes carry a non-finite ", col,
           " value. Every mCH panel in this section reads it.")
    }
  }

  cat(sprintf("  Universe: %s genes\n", format_count(nrow(universe))))
  cat(sprintf("  With a MeCP2 peak: %s (%.1f%%)\n",
              format_count(sum(universe$has_mecp2_peak)),
              100 * mean(universe$has_mecp2_peak)))
  cat(sprintf("  Neuronal: %s (%.1f%%) | Synapse: %s (%.1f%%) | MeCP2 up: %s (%.1f%%)\n",
              format_count(sum(universe$is_neuronal)),
              100 * mean(universe$is_neuronal),
              format_count(sum(universe$is_synapse)),
              100 * mean(universe$is_synapse),
              format_count(sum(universe$is_mecp2_up)),
              100 * mean(universe$is_mecp2_up)))
  for (lvl in SET_ORDER) {
    n <- sum(universe$gene_set == lvl)
    cat(sprintf("    %-22s %s (%.1f%%)\n", lvl, format_count(n),
                100 * n / nrow(universe)))
  }
  cat("\n")
  universe
}

write_per_gene_table <- function(universe, out_dir) {
  cols <- c("gene_name", "gene_id", "chr", "start", "end", "gene_length",
            "mch_ctrl", "mch_mut", "mch_diff", "edger_logFC", "edger_fdr",
            "mch_sig", "mch_direction", "mecp2_fold", "mecp2_fdr",
            "mecp2_n_peaks", "has_mecp2_peak", "is_mecp2_up", "is_neuronal",
            "is_synapse", "gene_set")
  out <- universe[order(universe$gene_set, universe$gene_name), cols]
  write_section_table(out, file.path(out_dir, "70_03_per_gene_gene_sets.tsv"))
}

#' Write the set membership of the pipeline's key genes.
write_key_gene_table <- function(universe, out_dir) {
  cat("--- Key gene membership ---\n")
  out <- universe[universe$gene_name %in% KEY_GENES,
                  c("gene_name", "gene_set", "is_neuronal", "is_synapse",
                    "is_mecp2_up", "mecp2_fold", "mecp2_fdr", "gene_length",
                    "mch_ctrl", "mch_diff", "edger_fdr", "mch_sig")]
  out <- out[order(out$gene_name), ]
  write_section_table(out, file.path(out_dir, "70_03_key_gene_membership.tsv"))

  for (i in seq_len(nrow(out))) {
    cat(sprintf("    %-10s %-22s mCH ctrl = %.5f, mCH diff = %.5f\n",
                out$gene_name[i], as.character(out$gene_set[i]),
                out$mch_ctrl[i], out$mch_diff[i]))
  }
  missing <- setdiff(KEY_GENES, out$gene_name)
  if (length(missing) > 0) {
    cat(sprintf("    Not in the universe: %s\n",
                paste(missing, collapse = ", ")))
  }
  cat("\n")
}

# =============================================================================
# STEP 3: OVERLAP COUNTS
# =============================================================================

#' Count the four cells of one set pair inside the universe.
count_set_pair <- function(universe, a_col, b_col, a_label, b_label) {
  a <- universe[[a_col]]
  b <- universe[[b_col]]
  n <- length(a)
  data.frame(
    set_a = a_label,
    set_b = b_label,
    n_universe = n,
    n_a = sum(a),
    n_b = sum(b),
    n_both = sum(a & b),
    n_a_only = sum(a & !b),
    n_b_only = sum(!a & b),
    n_neither = sum(!a & !b),
    pct_of_a_in_b = 100 * sum(a & b) / sum(a),
    pct_of_b_in_a = 100 * sum(a & b) / sum(b),
    pct_expected_of_a_in_b = 100 * sum(b) / n,
    stringsAsFactors = FALSE
  )
}

write_overlap_counts <- function(universe, out_dir) {
  cat("--- Step 3: overlap counts ---\n")

  counts <- rbind(
    count_set_pair(universe, "is_neuronal", "is_mecp2_up",
                   "Neuronal", "MeCP2 up"),
    count_set_pair(universe, "is_synapse", "is_mecp2_up",
                   "Synapse and axon", "MeCP2 up"),
    count_set_pair(universe, "is_neuronal", "is_synapse",
                   "Neuronal", "Synapse and axon")
  )

  for (i in seq_len(nrow(counts))) {
    r <- counts[i, ]
    cat(sprintf("  %s x %s: both %s, %s only %s, %s only %s, neither %s\n",
                r$set_a, r$set_b, format_count(r$n_both),
                r$set_a, format_count(r$n_a_only),
                r$set_b, format_count(r$n_b_only),
                format_count(r$n_neither)))
    cat(sprintf("    %.1f%% of %s genes are %s, against %.1f%% expected from the universe\n",
                r$pct_of_a_in_b, r$set_a, r$set_b, r$pct_expected_of_a_in_b))
  }

  write_section_table(counts, file.path(out_dir, "70_03a_overlap_counts.tsv"))
  cat("\n")
  counts
}

# =============================================================================
# STEP 4: FISHER TESTS
# =============================================================================

#' Run one registered gene-level Fisher test and return it as one summary row.
run_one_fisher <- function(universe, membership_col, membership_label,
                           test_id, description, out_dir) {
  gene_df <- universe[, c("gene_name", "chr", membership_col, "is_mecp2_up")]

  ft <- register_fisher_test(
    section = SECTION_ID, test_id = test_id,
    description = description,
    gene_df = gene_df, row_var = membership_col, col_var = "is_mecp2_up",
    output_dir = out_dir
  )

  if (!all(is.finite(ft$conf.int))) {
    stop("The Fisher confidence interval for ", membership_label,
         " against MeCP2 up is not finite. One cell of the 2x2 table is empty, ",
         "so the odds ratio cannot be plotted.")
  }

  a <- universe[[membership_col]]
  b <- universe$is_mecp2_up
  data.frame(
    comparison = sprintf("%s x MeCP2 up", membership_label),
    test_id = test_id,
    set_label = membership_label,
    n_universe = nrow(universe),
    n_set = sum(a),
    n_mecp2_up = sum(b),
    n_both = sum(a & b),
    odds_ratio = unname(ft$estimate),
    ci_lower = ft$conf.int[1],
    ci_upper = ft$conf.int[2],
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

run_fisher_tests <- function(universe, out_dir) {
  cat("--- Step 4: Fisher enrichment of MeCP2 gain in the gene sets ---\n")

  results <- rbind(
    run_one_fisher(
      universe, "is_neuronal", "Neuronal", "neuronal_vs_mecp2_up",
      paste("Neuronal genes are enriched among genes that gain MeCP2 binding",
            "in the BAP1 mutant."),
      out_dir),
    run_one_fisher(
      universe, "is_synapse", "Synapse and axon", "synapse_vs_mecp2_up",
      paste("Synapse and axon genes are enriched among genes that gain MeCP2",
            "binding in the BAP1 mutant."),
      out_dir)
  )
  results$p_adj <- p.adjust(results$p_value, method = "BH")
  results$log2_odds_ratio <- log2(results$odds_ratio)

  for (i in seq_len(nrow(results))) {
    r <- results[i, ]
    cat(sprintf("  %-28s OR = %.3f [%.3f, %.3f], %s (BH %s)\n",
                r$comparison, r$odds_ratio, r$ci_lower, r$ci_upper,
                format_p(r$p_value), format_p(r$p_adj)))
  }

  write_section_table(results, file.path(out_dir, "70_03d_fisher_enrichment.tsv"))
  cat("\n")
  results
}

# =============================================================================
# FIGURE 70_03a: TWO-SET VENN
# =============================================================================

plot_two_set_venn <- function(universe, counts, fisher_tbl, out_dir) {
  cat("--- Figure 70_03a: neuronal against MeCP2 up ---\n")

  venn_list <- list(
    "Neuronal" = universe$gene_name[universe$is_neuronal],
    "MeCP2 up" = universe$gene_name[universe$is_mecp2_up]
  )

  pair <- counts[counts$set_a == "Neuronal" & counts$set_b == "MeCP2 up", ]
  neuronal_fisher <- fisher_tbl[fisher_tbl$set_label == "Neuronal", ]

  caption <- paste(
    sprintf("Neuronal only %s | Both %s | MeCP2 up only %s | Neither %s",
            format_count(pair$n_a_only), format_count(pair$n_both),
            format_count(pair$n_b_only), format_count(pair$n_neither)),
    sprintf("Fisher: OR = %.2f [%.2f, %.2f], %s",
            neuronal_fisher$odds_ratio, neuronal_fisher$ci_lower,
            neuronal_fisher$ci_upper, format_p(neuronal_fisher$p_value)),
    sprintf("%.1f%% of neuronal genes gain MeCP2, against %.1f%% of the universe",
            pair$pct_of_a_in_b, pair$pct_expected_of_a_in_b),
    sep = "\n"
  )

  p <- ggVennDiagram(venn_list, label = "count", label_alpha = 0,
                     edge_size = 0.8) +
    scale_fill_gradient(low = "white", high = VENN_FILL_HIGH, guide = "none") +
    labs(
      title = "Neuronal genes and genes that gain MeCP2 binding",
      subtitle = sprintf("Universe: %s genes tested for mCH",
                         format_count(pair$n_universe)),
      caption = caption
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      plot.caption = element_text(hjust = 0, size = 9, color = "grey30")
    )

  save_multiformat_ggplot(p, file.path(out_dir, "70_03a_geneset_venn"),
                          width = 9, height = 8)
  cat("\n")
  p
}

# =============================================================================
# FIGURE 70_03d: FISHER FOREST
# =============================================================================

plot_fisher_forest <- function(fisher_tbl, out_dir) {
  cat("--- Figure 70_03d: Fisher odds ratios ---\n")

  plot_df <- fisher_tbl
  plot_df$comparison <- factor(plot_df$comparison,
                               levels = rev(plot_df$comparison))
  plot_df$point_label <- sprintf("OR = %.2f [%.2f, %.2f], %s (%s of %s in the set)",
                                 plot_df$odds_ratio, plot_df$ci_lower,
                                 plot_df$ci_upper, format_p(plot_df$p_value),
                                 format_count(plot_df$n_both),
                                 format_count(plot_df$n_set))
  plot_df$significance <- ifelse(plot_df$p_adj < Q_THRESHOLD,
                                 "Significant", "Not Significant")

  p <- ggplot(plot_df, aes(x = log2_odds_ratio, y = comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_errorbar(aes(xmin = log2(ci_lower), xmax = log2(ci_upper)),
                  orientation = "y", width = 0.18, linewidth = 0.8) +
    geom_point(aes(color = significance), size = 4) +
    geom_text(aes(label = point_label), vjust = -1.4, size = 3,
              color = "grey20") +
    scale_color_manual(values = COLORS$significant, name = "") +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.35))) +
    scale_y_discrete(expand = expansion(add = 0.8)) +
    labs(
      title = "Enrichment of MeCP2 gain in the neuronal gene sets",
      subtitle = sprintf("Universe: %s genes tested for mCH | %s gain MeCP2",
                         format_count(plot_df$n_universe[1]),
                         format_count(plot_df$n_mecp2_up[1])),
      x = expression(log[2](Odds~Ratio)), y = ""
    ) +
    theme_emseq() +
    theme(legend.position = "top",
          plot.subtitle = element_text(size = 9, color = "grey40"))

  save_multiformat_ggplot(p, file.path(out_dir, "70_03d_fisher_forest"),
                          width = 12, height = 6)
  cat("\n")
  p
}

# =============================================================================
# WILCOXON COMPARISONS BETWEEN GENE SETS
# =============================================================================

#' Build the set comparisons this section tests.
#'
#' The first two compare a whole gene class against the rest of the universe.
#' The last two compare the neuronal and MeCP2-up intersection against each of
#' the two single sets it sits inside.
build_comparisons <- function(universe) {
  neuronal <- universe$is_neuronal
  synapse  <- universe$is_synapse
  mecp2_up <- universe$is_mecp2_up

  list(
    list(a_label = "Neuronal", b_label = "Other genes",
         a = neuronal, b = !neuronal),
    list(a_label = "Synapse and axon", b_label = "Other genes",
         a = synapse, b = !synapse),
    list(a_label = "Neuronal + MeCP2 up", b_label = "Neuronal only",
         a = neuronal & mecp2_up, b = neuronal & !mecp2_up),
    list(a_label = "Neuronal + MeCP2 up", b_label = "MeCP2 up only",
         a = neuronal & mecp2_up, b = !neuronal & mecp2_up)
  )
}

#' Two-sided Wilcoxon rank sum test for every comparison on one value column.
run_wilcoxon_comparisons <- function(universe, value_col, comparisons) {
  values <- universe[[value_col]]

  rows <- lapply(comparisons, function(cmp) {
    a <- values[cmp$a & !is.na(values)]
    b <- values[cmp$b & !is.na(values)]
    if (length(a) < MIN_WILCOXON_GROUP || length(b) < MIN_WILCOXON_GROUP) {
      stop("Wilcoxon on ", value_col, " for ", cmp$a_label, " against ",
           cmp$b_label, " has ", length(a), " and ", length(b),
           " values. At least ", MIN_WILCOXON_GROUP, " are needed in each.")
    }
    wt <- wilcox.test(a, b)
    data.frame(
      value = value_col,
      group_a = cmp$a_label,
      group_b = cmp$b_label,
      n_a = length(a),
      n_b = length(b),
      median_a = median(a),
      median_b = median(b),
      median_difference = median(a) - median(b),
      mean_a = mean(a),
      mean_b = mean(b),
      W = unname(wt$statistic),
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out$p_adj <- p.adjust(out$p_value, method = "BH")
  out
}

report_wilcoxon <- function(tbl, digits = 5) {
  fmt <- paste0("  %-22s vs %-20s n = %s vs %s, median %.", digits,
                "f vs %.", digits, "f, %s\n")
  for (i in seq_len(nrow(tbl))) {
    cat(sprintf(fmt, tbl$group_a[i], tbl$group_b[i],
                format_count(tbl$n_a[i]), format_count(tbl$n_b[i]),
                tbl$median_a[i], tbl$median_b[i], format_p(tbl$p_value[i])))
  }
}

# =============================================================================
# PER-SET DISTRIBUTION FIGURES
# =============================================================================

#' Violin of one value column across the four gene sets.
#'
#' summary_df needs a label column built by group_label(), so that every set
#' carries its gene count and its median on the figure.
build_set_violin <- function(df, value_col, summary_df, title, subtitle,
                             y_label, caption = NULL, log_y = FALSE,
                             zero_line = FALSE) {
  values <- df[[value_col]]
  values <- values[!is.na(values)]

  if (log_y) {
    if (any(values <= 0)) {
      stop("build_set_violin(): ", value_col, " holds a value of zero or less, ",
           "which a log10 axis cannot show.")
    }
    label_y <- max(values) * 2.4
    y_scale <- scale_y_log10(labels = scales::comma,
                             limits = c(min(values) * 0.7, max(values) * 6))
  } else {
    span <- max(values) - min(values)
    label_y <- max(values) + 0.10 * span
    y_scale <- scale_y_continuous(
      limits = c(min(values) - 0.05 * span, max(values) + 0.24 * span))
  }

  label_df <- summary_df
  label_df$label_y <- label_y

  p <- ggplot(df, aes(x = gene_set, y = .data[[value_col]], fill = gene_set))

  if (zero_line) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "grey40",
                        linewidth = 0.4)
  }

  p +
    geom_violin(alpha = 0.65, scale = "width",
                draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white",
                 alpha = 0.85) +
    geom_text(data = label_df,
              aes(x = gene_set, y = label_y, label = label),
              inherit.aes = FALSE, size = 3.1, lineheight = 1.1) +
    y_scale +
    scale_fill_manual(values = SET_COLORS) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12)) +
    labs(title = title, subtitle = subtitle, x = "Gene set", y = y_label,
         caption = caption) +
    theme_emseq() +
    theme(legend.position = "none",
          plot.subtitle = element_text(size = 9, color = "grey40"),
          plot.caption = element_text(hjust = 0, size = 8, color = "grey30"))
}

# =============================================================================
# ANALYSIS: mCH LEVEL PER GENE SET
# =============================================================================

analyse_mch_level <- function(universe, comparisons, out_dir) {
  cat("--- Figure 70_03b: mCH level per gene set ---\n")

  summary_ctrl <- set_summary(universe, "gene_set", "mch_ctrl")
  summary_mut  <- set_summary(universe, "gene_set", "mch_mut")
  summary_ctrl$condition <- "Control"
  summary_mut$condition <- "Mutant"
  combined <- rbind(summary_ctrl, summary_mut)
  # combined holds data only. The figure below annotates from its own copy.
  write_section_table(combined,
                      file.path(out_dir, "70_03b_mch_level_summary.tsv"))

  for (i in seq_len(nrow(summary_ctrl))) {
    cat(sprintf("  %-22s n = %s, control mCH median = %.5f, mutant = %.5f\n",
                as.character(summary_ctrl$gene_set[i]),
                format_count(summary_ctrl$n[i]), summary_ctrl$median[i],
                summary_mut$median[i]))
  }

  wilcox_tbl <- run_wilcoxon_comparisons(universe, "mch_ctrl", comparisons)
  report_wilcoxon(wilcox_tbl)
  write_section_table(wilcox_tbl,
                  file.path(out_dir, "70_03b_mch_ctrl_wilcoxon.tsv"))

  caption <- paste(
    sprintf("Wilcoxon %s vs %s: %s", wilcox_tbl$group_a, wilcox_tbl$group_b,
            format_p(wilcox_tbl$p_value)),
    collapse = "\n"
  )

  plot_summary <- summary_ctrl
  plot_summary$label <- group_label(plot_summary, digits = 5)

  p <- build_set_violin(
    universe, "mch_ctrl", plot_summary,
    title = "Control mCH level by gene set",
    subtitle = "Gene-body non-CG methylation rate in the control condition",
    y_label = "mCH rate (control)",
    caption = caption
  )
  save_multiformat_ggplot(p, file.path(out_dir, "70_03b_mch_level_by_set"),
                          width = 10, height = 8)
  cat("\n")
  list(plot = p, summary = combined, wilcoxon = wilcox_tbl)
}

# =============================================================================
# ANALYSIS: mCH CHANGE PER GENE SET
# =============================================================================

#' Paired Wilcoxon of mutant against control mCH inside each gene set.
run_paired_condition_tests <- function(universe) {
  rows <- lapply(SET_ORDER, function(lvl) {
    sub <- universe[universe$gene_set == lvl, , drop = FALSE]
    sub <- sub[!is.na(sub$mch_ctrl) & !is.na(sub$mch_mut), , drop = FALSE]
    if (nrow(sub) < MIN_WILCOXON_GROUP) {
      stop("Gene set ", lvl, " has ", nrow(sub),
           " genes with both mCH rates. At least ", MIN_WILCOXON_GROUP,
           " are needed for the paired test.")
    }
    wt <- wilcox.test(sub$mch_mut, sub$mch_ctrl, paired = TRUE)
    data.frame(
      gene_set = lvl,
      n_genes = nrow(sub),
      median_mch_ctrl = median(sub$mch_ctrl),
      median_mch_mut = median(sub$mch_mut),
      median_paired_delta = median(sub$mch_mut - sub$mch_ctrl),
      median_mch_diff = median(sub$mch_diff),
      n_mch_significant = sum(sub$mch_sig),
      n_hypermethylated = sum(sub$mch_sig & sub$mch_diff > 0),
      n_hypomethylated = sum(sub$mch_sig & sub$mch_diff < 0),
      V = unname(wt$statistic),
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$p_adj <- p.adjust(out$p_value, method = "BH")
  out
}

analyse_mch_change <- function(universe, comparisons, out_dir) {
  cat("--- Figure 70_03c: mCH change per gene set ---\n")

  summary_diff <- set_summary(universe, "gene_set", "mch_diff")
  # summary_diff holds data only. The figure below annotates from its own copy.
  write_section_table(summary_diff,
                      file.path(out_dir, "70_03c_mch_diff_summary.tsv"))

  paired <- run_paired_condition_tests(universe)
  for (i in seq_len(nrow(paired))) {
    cat(sprintf("  %-22s n = %s, paired median delta = %.5f, %s | %s of %s genes significant\n",
                paired$gene_set[i], format_count(paired$n_genes[i]),
                paired$median_paired_delta[i], format_p(paired$p_value[i]),
                format_count(paired$n_mch_significant[i]),
                format_count(paired$n_genes[i])))
  }
  write_section_table(paired,
                  file.path(out_dir, "70_03c_mch_paired_condition_tests.tsv"))

  wilcox_tbl <- run_wilcoxon_comparisons(universe, "mch_diff", comparisons)
  report_wilcoxon(wilcox_tbl)
  write_section_table(wilcox_tbl,
                  file.path(out_dir, "70_03c_mch_diff_wilcoxon.tsv"))

  caption <- paste(
    paste(sprintf("Wilcoxon %s vs %s: %s", wilcox_tbl$group_a,
                  wilcox_tbl$group_b, format_p(wilcox_tbl$p_value)),
          collapse = "\n"),
    paste(sprintf("Paired mutant vs control, %s: median delta %.5f, %s",
                  paired$gene_set, paired$median_paired_delta,
                  format_p(paired$p_value)),
          collapse = "\n"),
    sep = "\n"
  )

  plot_summary <- summary_diff
  plot_summary$label <- group_label(plot_summary, digits = 5)

  p <- build_set_violin(
    universe, "mch_diff", plot_summary,
    title = "mCH change by gene set",
    subtitle = "Mutant minus control gene-body mCH rate",
    y_label = "mCH difference (mutant - control)",
    caption = caption,
    zero_line = TRUE
  )
  save_multiformat_ggplot(p, file.path(out_dir, "70_03c_mch_diff_by_set"),
                          width = 10, height = 8)
  cat("\n")
  list(plot = p, summary = summary_diff, wilcoxon = wilcox_tbl, paired = paired)
}

# =============================================================================
# ANALYSIS: GENE LENGTH PER GENE SET
# =============================================================================

analyse_gene_length <- function(universe, comparisons, out_dir) {
  cat("--- Figure 70_03e: gene length per gene set ---\n")

  summary_len <- set_summary(universe, "gene_set", "gene_length")
  # summary_len holds lengths in base pairs, data only. The figure below prints
  # the median in kb from its own copy.
  write_section_table(summary_len,
                      file.path(out_dir, "70_03e_gene_length_summary.tsv"))

  for (i in seq_len(nrow(summary_len))) {
    cat(sprintf("  %-22s n = %s, median length = %s bp\n",
                as.character(summary_len$gene_set[i]),
                format_count(summary_len$n[i]),
                format_count(round(summary_len$median[i]))))
  }

  wilcox_tbl <- run_wilcoxon_comparisons(universe, "gene_length", comparisons)
  report_wilcoxon(wilcox_tbl, digits = 0)
  write_section_table(wilcox_tbl,
                  file.path(out_dir, "70_03e_gene_length_wilcoxon.tsv"))

  caption <- paste(
    "Gene length is associated with both mCH and MeCP2 binding, so it is",
    "shown for every set that the mCH panels compare.",
    paste(sprintf("Wilcoxon %s vs %s: %s", wilcox_tbl$group_a,
                  wilcox_tbl$group_b, format_p(wilcox_tbl$p_value)),
          collapse = "\n"),
    sep = "\n"
  )

  # group_label() reads the median column, so the kb label comes from a copy
  # whose median is divided by 1000. The plotted table keeps base pairs.
  length_kb <- summary_len
  length_kb$median <- length_kb$median / 1000
  plot_summary <- summary_len
  plot_summary$label <- group_label(length_kb, digits = 1, unit = " kb")

  p <- build_set_violin(
    universe, "gene_length", plot_summary,
    title = "Gene length by gene set",
    subtitle = "Annotated gene body length, log10 axis",
    y_label = "Gene length (bp)",
    caption = caption,
    log_y = TRUE
  )
  save_multiformat_ggplot(p, file.path(out_dir, "70_03e_gene_length_by_set"),
                          width = 10, height = 8)
  cat("\n")
  list(plot = p, summary = summary_len, wilcoxon = wilcox_tbl)
}

# =============================================================================
# FIGURE 70_03f: COMPOSITE
# =============================================================================

plot_composite <- function(p_venn, p_forest, p_level, p_diff, p_length,
                           universe, out_dir) {
  cat("--- Figure 70_03f: composite ---\n")

  composite <- (p_venn | p_forest) / (p_level | p_diff | p_length) +
    plot_annotation(
      title = "Section 70_03: neuronal gene sets, MeCP2 gain and mCH",
      subtitle = sprintf(
        "Universe %s genes | Neuronal %s | Synapse and axon %s | MeCP2 up %s | Both %s",
        format_count(nrow(universe)),
        format_count(sum(universe$is_neuronal)),
        format_count(sum(universe$is_synapse)),
        format_count(sum(universe$is_mecp2_up)),
        format_count(sum(universe$is_neuronal & universe$is_mecp2_up))),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5)
      )
    )

  save_multiformat_ggplot(composite, file.path(out_dir, "70_03f_composite"),
                          width = 22, height = 16)
  cat("\n")
  invisible(composite)
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
  cat("  SECTION 70_03: GENE SET OVERLAP WITH MeCP2 GAIN\n")
  cat("===========================================================================\n")
  cat("  Output directory: ", out_dir, "\n", sep = "")
  cat("  MeCP2 FDR threshold: ", fdr_threshold, "\n", sep = "")
  cat("  Universe: every gene tested for mCH, one row per gene symbol\n\n")

  neuronal_genes <- load_neuronal_handoff()

  cat("--- Step 1b: synapse and axon gene set ---\n")
  synapse_genes <- load_gene_set(GENESET_PATHS$synapse,
                                 "Synapse and axon gene set")
  cat("\n")

  mecp2_gene_tbl <- aggregate_mecp2_to_genes(fdr_threshold)

  universe <- build_universe(neuronal_genes, synapse_genes, mecp2_gene_tbl,
                             fdr_threshold)
  write_per_gene_table(universe, out_dir)
  write_key_gene_table(universe, out_dir)

  counts <- write_overlap_counts(universe, out_dir)
  fisher_tbl <- run_fisher_tests(universe, out_dir)

  p_venn <- plot_two_set_venn(universe, counts, fisher_tbl, out_dir)
  p_forest <- plot_fisher_forest(fisher_tbl, out_dir)

  comparisons <- build_comparisons(universe)
  level <- analyse_mch_level(universe, comparisons, out_dir)
  change <- analyse_mch_change(universe, comparisons, out_dir)
  length_result <- analyse_gene_length(universe, comparisons, out_dir)

  plot_composite(p_venn, p_forest, level$plot, change$plot,
                 length_result$plot, universe, out_dir)

  cat("---------------------------------------------------------------------------\n")
  cat("  SUMMARY\n")
  cat("---------------------------------------------------------------------------\n")
  cat(sprintf("  Universe: %s genes, %s neuronal, %s synapse and axon, %s MeCP2 up\n",
              format_count(nrow(universe)),
              format_count(sum(universe$is_neuronal)),
              format_count(sum(universe$is_synapse)),
              format_count(sum(universe$is_mecp2_up))))
  cat(sprintf("  Neuronal and MeCP2 up: %s genes\n",
              format_count(sum(universe$is_neuronal & universe$is_mecp2_up))))
  for (i in seq_len(nrow(fisher_tbl))) {
    r <- fisher_tbl[i, ]
    cat(sprintf("  Fisher %-28s OR = %.3f [%.3f, %.3f], %s\n",
                r$comparison, r$odds_ratio, r$ci_lower, r$ci_upper,
                format_p(r$p_value)))
  }
  ctrl_rows <- level$summary[level$summary$condition == "Control", ]
  for (i in seq_len(nrow(ctrl_rows))) {
    cat(sprintf("  %-22s control mCH median = %.5f\n",
                as.character(ctrl_rows$gene_set[i]), ctrl_rows$median[i]))
  }
  for (i in seq_len(nrow(change$summary))) {
    cat(sprintf("  %-22s mCH diff median = %.5f\n",
                as.character(change$summary$gene_set[i]),
                change$summary$median[i]))
  }
  for (i in seq_len(nrow(length_result$summary))) {
    cat(sprintf("  %-22s median gene length = %s bp\n",
                as.character(length_result$summary$gene_set[i]),
                format_count(round(length_result$summary$median[i]))))
  }
  cat("\n=== Section 70_03 complete ===\n\n")
}

main()
