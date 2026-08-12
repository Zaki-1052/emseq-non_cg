# scripts/sections/70_neuronal/70_01_k119ub_neuronal.R
#
# Section 70_01: H2AK119ub is constitutively high at neuronal genes.
#
# Asks whether neuronal, synaptic and axon-guidance genes carry a high
# H2AK119ub level in the control condition. The measure is the absolute
# gene-body ChIP signal, not a fold change, so the answer does not depend on
# the mutant condition. No methylation data and no MeCP2 data enter this
# section: the question is whether the neuronal gene class is K119ub-high on
# its own terms.
#
# Long genes could produce the same association without any neuronal identity
# effect, so the section compares gene length between neuronal and other genes
# and repeats the K119ub comparison inside gene-length strata.
#
# Reads:
#   DIFFBIND_PATHS$k119ub_gene_signal  per-gene H2AK119ub gene-body and
#                                      promoter signal in ctrl and mut
#   GENESET_PATHS$neuronal             external neuronal gene set
#   org.Mm.eg.db and GO.db             GO biological process annotations, used
#                                      to derive a second neuronal set fresh
#
# Writes into --output-dir (default results/sections/70_neuronal/):
#   neuronal_gene_set.tsv    handoff read by sections 70_02, 70_03 and 70_04.
#                            Columns: gene, source, is_derived, is_external,
#                            k119ub_ctrl_signal, k119ub_log2fc, k119ub_decile.
#                            The three k119ub columns hold values for genes in
#                            the quantifiable K119ub universe and NA for every
#                            other gene of the union set.
#   70_01a .. 70_01j         figures, each beside the TSV holding its numbers
#   70_01_per_gene_k119ub_neuronal.tsv   per-gene table for the whole universe
#   fisher_tables/           gene tables behind the registered Fisher tests
#
# Adapted from Biomodal section 72 (K119ub neuronal characterisation).

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)
library(GO.db)

# =============================================================================
# SECTION CONSTANTS
# =============================================================================

SECTION_ID <- "70_01"

# GO biological process term names matching this pattern define the derived
# neuronal gene set, and mark a GO or GSEA result term as neuronal.
NEURONAL_PATTERN <- "synap|neuron|axon|dendrit|nervous"

# Columns the K119ub gene signal table must carry.
SIGNAL_COLUMNS <- c("entrez_id", "symbol", "chr", "start", "end", "width",
                    "gb_ctrl_signal", "gb_mut_signal", "gb_log2fc",
                    "gb_signal_class", "pr_ctrl_signal", "pr_mut_signal",
                    "pr_log2fc", "pr_signal_class")

# Genes of this signal class carry a usable level in both conditions. The
# other classes ("one_condition", "no_signal") are dropped from every analysis.
QUANTIFIABLE_CLASS <- "quantifiable"

# Exact column contract of the cross-section handoff file.
HANDOFF_COLUMNS <- c("gene", "source", "is_derived", "is_external",
                     "k119ub_ctrl_signal", "k119ub_log2fc", "k119ub_decile")

GROUP_ORDER <- c("Neuronal", "Other")

GROUP_COLORS <- c(
  "Neuronal" = COLORS$k119ub[["K119ub Gained"]],
  "Other"    = "grey70"
)

CONDITION_ORDER <- c("Control", "Mutant")

N_DECILES <- 10
N_LENGTH_STRATA <- 5

# Number of top-ranked significant GSEA terms counted as the leading terms.
LEADING_TERM_COUNT <- 50

# Terms shown in the GO and GSEA dot plots.
DOTPLOT_CATEGORIES <- 25

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
                help = "FDR cutoff for GO and GSEA terms [default: %default]")
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

#' Summarise n, median and mean for every combination of two grouping columns.
#'
#' Two-column version of summarise_groups(). Returns data only. Call
#' group_label() on the result at the plot site to build the annotation text.
summarise_groups_by <- function(df, facet_col, group_col, value_col) {
  out <- df %>%
    dplyr::filter(!is.na(.data[[value_col]])) %>%
    dplyr::group_by(.data[[facet_col]], .data[[group_col]]) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median = median(.data[[value_col]]),
      mean = mean(.data[[value_col]]),
      q25 = unname(quantile(.data[[value_col]], 0.25)),
      q75 = unname(quantile(.data[[value_col]], 0.75)),
      .groups = "drop"
    ) %>%
    as.data.frame()

  colnames(out)[1] <- facet_col
  colnames(out)[2] <- group_col
  out
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
# STEP 1: K119ub GENE SIGNAL
# =============================================================================

load_k119ub_signal <- function() {
  cat("--- Step 1: K119ub gene signal ---\n")

  path <- DIFFBIND_PATHS$k119ub_gene_signal
  if (!file.exists(path)) {
    stop("K119ub gene signal table not found: ", path,
         "\nRun scripts/utils/copy_reference_data.sh to place it under data/.")
  }

  signal <- read.table(path, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "")

  missing <- setdiff(SIGNAL_COLUMNS, colnames(signal))
  if (length(missing) > 0) {
    stop("K119ub gene signal table is missing columns: ",
         paste(missing, collapse = ", "), " in ", path)
  }

  signal$entrez_id <- as.character(signal$entrez_id)
  if (anyDuplicated(signal$entrez_id) > 0) {
    stop("K119ub gene signal table repeats an Entrez ID. ",
         "Every GO and GSEA call in this section keys on it.")
  }
  if (anyDuplicated(signal$symbol) > 0) {
    stop("K119ub gene signal table repeats a gene symbol. ",
         "Gene-set membership in this section keys on the symbol.")
  }

  class_counts <- table(signal$gb_signal_class)
  cat(sprintf("  %s genes in the table\n", format(nrow(signal), big.mark = ",")))
  for (cls in names(class_counts)) {
    cat(sprintf("    %-14s %s\n", cls,
                format(as.integer(class_counts[[cls]]), big.mark = ",")))
  }

  universe <- signal[signal$gb_signal_class == QUANTIFIABLE_CLASS, , drop = FALSE]
  if (nrow(universe) < 1000) {
    stop("Only ", nrow(universe), " genes carry a quantifiable K119ub gene-body ",
         "signal. Too few to characterise the neuronal gene class.")
  }

  incomplete <- !is.finite(universe$gb_ctrl_signal) |
    !is.finite(universe$gb_mut_signal) |
    !is.finite(universe$gb_log2fc) |
    !is.finite(universe$width)
  if (any(incomplete)) {
    stop(sum(incomplete), " quantifiable genes carry a non-finite signal, ",
         "log2FC or width value. The table is inconsistent with its own class ",
         "labels.")
  }
  if (any(universe$gb_ctrl_signal <= 0) || any(universe$gb_mut_signal <= 0)) {
    stop("A quantifiable gene carries a K119ub signal of zero or less. ",
         "The signal figures use a log10 axis, which cannot show it.")
  }

  universe$gene_name <- universe$symbol
  cat(sprintf("  Working universe: %s quantifiable genes\n",
              format(nrow(universe), big.mark = ",")))
  cat(sprintf("  Ctrl gene-body signal range: %.3f to %.3f\n\n",
              min(universe$gb_ctrl_signal), max(universe$gb_ctrl_signal)))

  universe
}

# =============================================================================
# STEP 2: NEURONAL GENE SETS
# =============================================================================

#' Derive a neuronal gene set from every GO biological process term whose name
#' matches NEURONAL_PATTERN, using the GOALL mapping so that genes annotated to
#' a descendant term are included.
derive_neuronal_gene_set <- function(out_dir) {
  cat("--- Step 2a: deriving the neuronal gene set from GO ---\n")
  cat(sprintf("  Term name pattern: %s\n", NEURONAL_PATTERN))

  go_ids <- AnnotationDbi::keys(GO.db, keytype = "GOID")
  go_terms <- AnnotationDbi::select(GO.db, keys = go_ids, keytype = "GOID",
                                    columns = c("GOID", "TERM", "ONTOLOGY"))

  bp_terms <- go_terms[go_terms$ONTOLOGY == "BP" & !is.na(go_terms$TERM), ,
                       drop = FALSE]
  neuro_terms <- bp_terms[grepl(NEURONAL_PATTERN, bp_terms$TERM,
                                ignore.case = TRUE), , drop = FALSE]
  if (nrow(neuro_terms) == 0) {
    stop("No GO biological process term name matches '", NEURONAL_PATTERN,
         "'. The GO.db install cannot support this section.")
  }

  cat(sprintf("  %s of %s GO BP terms match\n",
              format(nrow(neuro_terms), big.mark = ","),
              format(nrow(bp_terms), big.mark = ",")))
  cat(sprintf("  Examples: %s\n",
              paste(head(neuro_terms$TERM, 4), collapse = "; ")))

  mapping <- AnnotationDbi::select(org.Mm.eg.db, keys = neuro_terms$GOID,
                                   keytype = "GOALL",
                                   columns = c("ENTREZID", "SYMBOL"))
  mapping <- mapping[!is.na(mapping$SYMBOL) & nzchar(mapping$SYMBOL), ,
                     drop = FALSE]
  genes <- sort(unique(mapping$SYMBOL))
  if (length(genes) == 0) {
    stop("The matching GO terms map to no mouse gene symbol. ",
         "Check the org.Mm.eg.db install.")
  }

  if (!"GOALL" %in% colnames(mapping)) {
    stop("The GOALL mapping has no GOALL column, so the derived gene set ",
         "cannot be traced back to its GO terms. Columns returned: ",
         paste(colnames(mapping), collapse = ", "))
  }
  term_gene_pairs <- unique(mapping[, c("GOALL", "SYMBOL")])
  gene_counts <- table(term_gene_pairs$GOALL)

  term_table <- neuro_terms[, c("GOID", "TERM")]
  term_table$n_mouse_genes <- as.integer(gene_counts[term_table$GOID])
  term_table$n_mouse_genes[is.na(term_table$n_mouse_genes)] <- 0L
  term_table <- term_table[order(-term_table$n_mouse_genes), ]
  write_section_table(term_table,
                  file.path(out_dir, "70_01_neuronal_go_terms.tsv"))

  cat(sprintf("  Derived neuronal gene set: %s symbols\n\n",
              format(length(genes), big.mark = ",")))
  genes
}

#' Compare the derived set with the external set and write the overlap figure.
report_gene_set_overlap <- function(derived, external, universe, out_dir) {
  cat("--- Step 2b: derived against external neuronal gene set ---\n")

  union_genes <- sort(union(derived, external))
  both <- intersect(derived, external)
  derived_only <- setdiff(derived, external)
  external_only <- setdiff(external, derived)
  jaccard <- length(both) / length(union_genes)

  cat(sprintf("  Derived only: %s | External only: %s | Both: %s\n",
              format(length(derived_only), big.mark = ","),
              format(length(external_only), big.mark = ","),
              format(length(both), big.mark = ",")))
  cat(sprintf("  Union: %s genes, Jaccard index = %.3f\n",
              format(length(union_genes), big.mark = ","), jaccard))
  cat(sprintf("  In the K119ub universe: derived %s, external %s, both %s\n\n",
              format(sum(universe$gene_name %in% derived), big.mark = ","),
              format(sum(universe$gene_name %in% external), big.mark = ","),
              format(sum(universe$gene_name %in% both), big.mark = ",")))

  overlap <- data.frame(
    category = c("Derived only", "Both", "External only"),
    n_genes = c(length(derived_only), length(both), length(external_only)),
    n_in_k119ub_universe = c(
      sum(universe$gene_name %in% derived_only),
      sum(universe$gene_name %in% both),
      sum(universe$gene_name %in% external_only)
    ),
    stringsAsFactors = FALSE
  )
  overlap$n_derived_total <- length(derived)
  overlap$n_external_total <- length(external)
  overlap$n_union_total <- length(union_genes)
  overlap$jaccard_index <- jaccard
  write_section_table(overlap, file.path(out_dir, "70_01a_gene_set_overlap.tsv"))

  plot_df <- overlap[, c("category", "n_genes", "n_in_k119ub_universe")]
  plot_df <- tidyr::pivot_longer(plot_df,
                                 cols = c("n_genes", "n_in_k119ub_universe"),
                                 names_to = "scope", values_to = "n")
  plot_df$scope <- factor(
    ifelse(plot_df$scope == "n_genes", "All annotated genes",
           "In the K119ub universe"),
    levels = c("All annotated genes", "In the K119ub universe")
  )
  plot_df$category <- factor(plot_df$category,
                             levels = c("Derived only", "Both", "External only"))

  figure <- ggplot(plot_df, aes(x = category, y = n, fill = scope)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = format(n, big.mark = ",", trim = TRUE)),
              position = position_dodge(width = 0.75), vjust = -0.4, size = 3.4) +
    scale_fill_manual(values = c("All annotated genes" = "#4D4D4D",
                                 "In the K119ub universe" = GROUP_COLORS[["Neuronal"]]),
                      name = "") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Derived and external neuronal gene sets",
      subtitle = sprintf(paste("Derived from GO BP terms matching '%s' |",
                               "Jaccard index = %.3f"),
                         NEURONAL_PATTERN, jaccard),
      x = "", y = "Genes"
    ) +
    theme_emseq() +
    theme(legend.position = "top",
          plot.subtitle = element_text(size = 9, color = "grey40"))

  save_multiformat_ggplot(figure, file.path(out_dir, "70_01a_gene_set_overlap"),
                          width = 9, height = 6)
  cat("\n")

  list(n_union = length(union_genes), n_both = length(both),
       n_derived_only = length(derived_only),
       n_external_only = length(external_only), jaccard = jaccard)
}

# =============================================================================
# STEP 3: ANNOTATE THE UNIVERSE
# =============================================================================

#' Add neuronal membership, K119ub signal deciles, high-signal flags and
#' gene-length strata to every gene of the quantifiable universe.
annotate_universe <- function(universe, derived, external, out_dir) {
  cat("--- Step 3: annotating the K119ub universe ---\n")

  universe$is_neuronal <- universe$gene_name %in% derived
  universe$is_neuronal_external <- universe$gene_name %in% external
  universe$group <- factor(ifelse(universe$is_neuronal, "Neuronal", "Other"),
                           levels = GROUP_ORDER)

  universe$k119ub_decile <- cut_into_bins(universe$gb_ctrl_signal, N_DECILES,
                                          "ctrl K119ub signal")

  ctrl_q3 <- quantile(universe$gb_ctrl_signal, 0.75)
  ctrl_d9 <- quantile(universe$gb_ctrl_signal, 0.90)
  mut_q3  <- quantile(universe$gb_mut_signal, 0.75)
  mut_d9  <- quantile(universe$gb_mut_signal, 0.90)

  universe$ctrl_top_quartile <- universe$gb_ctrl_signal >= ctrl_q3
  universe$ctrl_top_decile   <- universe$gb_ctrl_signal >= ctrl_d9
  universe$mut_top_quartile  <- universe$gb_mut_signal >= mut_q3
  universe$mut_top_decile    <- universe$gb_mut_signal >= mut_d9

  universe$length_stratum <- cut_into_bins(universe$width, N_LENGTH_STRATA,
                                           "gene length")

  cat(sprintf("  Neuronal (derived):  %s of %s genes (%.1f%%)\n",
              format(sum(universe$is_neuronal), big.mark = ","),
              format(nrow(universe), big.mark = ","),
              100 * mean(universe$is_neuronal)))
  cat(sprintf("  Neuronal (external): %s of %s genes (%.1f%%)\n",
              format(sum(universe$is_neuronal_external), big.mark = ","),
              format(nrow(universe), big.mark = ","),
              100 * mean(universe$is_neuronal_external)))
  cat(sprintf("  Ctrl signal thresholds: Q3 = %.4f, D9 = %.4f\n",
              ctrl_q3, ctrl_d9))
  cat(sprintf("  Mut  signal thresholds: Q3 = %.4f, D9 = %.4f\n",
              mut_q3, mut_d9))

  thresholds <- data.frame(
    flag = c("ctrl_top_quartile", "ctrl_top_decile",
             "mut_top_quartile", "mut_top_decile"),
    signal_column = c("gb_ctrl_signal", "gb_ctrl_signal",
                      "gb_mut_signal", "gb_mut_signal"),
    quantile = c(0.75, 0.90, 0.75, 0.90),
    threshold = c(unname(ctrl_q3), unname(ctrl_d9),
                  unname(mut_q3), unname(mut_d9)),
    n_genes_at_or_above = c(sum(universe$ctrl_top_quartile),
                            sum(universe$ctrl_top_decile),
                            sum(universe$mut_top_quartile),
                            sum(universe$mut_top_decile)),
    stringsAsFactors = FALSE
  )
  write_section_table(thresholds,
                  file.path(out_dir, "70_01_signal_thresholds.tsv"))

  decile_bounds <- universe %>%
    dplyr::group_by(k119ub_decile) %>%
    dplyr::summarise(n_genes = dplyr::n(),
                     signal_min = min(gb_ctrl_signal),
                     signal_max = max(gb_ctrl_signal),
                     .groups = "drop") %>%
    as.data.frame()
  write_section_table(decile_bounds,
                  file.path(out_dir, "70_01_decile_bounds.tsv"))
  cat("\n")

  universe
}

# =============================================================================
# FIGURE 70_01b: K119ub SIGNAL, NEURONAL AGAINST OTHER
# =============================================================================

plot_signal_distribution <- function(universe, out_dir) {
  cat("--- Figure 70_01b: K119ub signal, neuronal against other ---\n")

  neuronal <- universe$is_neuronal
  ctrl_test <- wilcox.test(universe$gb_ctrl_signal[neuronal],
                           universe$gb_ctrl_signal[!neuronal],
                           alternative = "greater")
  mut_test <- wilcox.test(universe$gb_mut_signal[neuronal],
                          universe$gb_mut_signal[!neuronal],
                          alternative = "greater")

  long <- rbind(
    data.frame(condition = "Control", group = universe$group,
               signal = universe$gb_ctrl_signal, stringsAsFactors = FALSE),
    data.frame(condition = "Mutant", group = universe$group,
               signal = universe$gb_mut_signal, stringsAsFactors = FALSE)
  )
  long$condition <- factor(long$condition, levels = CONDITION_ORDER)
  long$group <- factor(long$group, levels = GROUP_ORDER)

  labels <- summarise_groups_by(long, "condition", "group", "signal")
  labels$label <- group_label(labels, digits = 3)
  labels$label_y <- max(long$signal) * 1.9

  figure <- ggplot(long, aes(x = group, y = signal, fill = group)) +
    geom_violin(scale = "width", alpha = 0.65, linewidth = 0.3,
                show.legend = FALSE) +
    geom_boxplot(width = 0.14, outlier.size = 0.2, outlier.alpha = 0.2,
                 show.legend = FALSE) +
    geom_text(data = labels, aes(x = group, y = label_y, label = label),
              inherit.aes = FALSE, size = 3.2, lineheight = 1.1) +
    facet_wrap(~ condition) +
    scale_fill_manual(values = GROUP_COLORS) +
    scale_y_log10(expand = expansion(mult = c(0.05, 0.18))) +
    labs(
      title = "H2AK119ub gene-body signal at neuronal and other genes",
      subtitle = sprintf(paste("One-sided Wilcoxon, neuronal greater |",
                               "Control %s | Mutant %s"),
                         format_p(ctrl_test$p.value), format_p(mut_test$p.value)),
      x = "", y = "K119ub gene-body signal (log10 scale)"
    ) +
    theme_emseq() +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"))

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "70_01b_k119ub_signal_neuronal_vs_other"),
                          width = 11, height = 7)

  stats_df <- labels[, c("condition", "group", "n", "median", "mean")]
  colnames(stats_df) <- c("condition", "group", "n_genes", "median_signal",
                          "mean_signal")
  stats_df$wilcoxon_p <- ifelse(stats_df$condition == "Control",
                                ctrl_test$p.value, mut_test$p.value)
  stats_df$wilcoxon_alternative <- "neuronal greater"
  write_section_table(stats_df,
                  file.path(out_dir, "70_01b_signal_distribution.tsv"))
  cat("\n")

  list(ctrl_p = ctrl_test$p.value, mut_p = mut_test$p.value, stats = stats_df)
}

# =============================================================================
# FIGURE 70_01c: NEURONAL ENRICHMENT AMONG K119ub-HIGH GENES
# =============================================================================

#' Register one Fisher test of neuronal membership against a high-signal flag
#' and return the odds ratio row behind it.
test_high_signal_enrichment <- function(universe, membership_col, high_col,
                                        test_id, label, description, out_dir) {
  gene_df <- data.frame(
    gene_name = universe$gene_name,
    chr = universe$chr,
    stringsAsFactors = FALSE
  )
  gene_df[[membership_col]] <- universe[[membership_col]]
  gene_df[[high_col]] <- universe[[high_col]]

  ft <- register_fisher_test(
    section = SECTION_ID, test_id = test_id, description = description,
    gene_df = gene_df, row_var = membership_col, col_var = high_col,
    output_dir = out_dir
  )

  high <- universe[[high_col]]
  member <- universe[[membership_col]]
  data.frame(
    test = label,
    test_id = test_id,
    gene_set = membership_col,
    n_genes = nrow(universe),
    n_high = sum(high),
    n_member_high = sum(high & member),
    frac_member_high = sum(high & member) / sum(high),
    frac_member_universe = mean(member),
    odds_ratio = unname(ft$estimate),
    ci_lower = ft$conf.int[1],
    ci_upper = ft$conf.int[2],
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

plot_high_signal_enrichment <- function(universe, out_dir) {
  cat("--- Figure 70_01c: neuronal enrichment among K119ub-high genes ---\n")

  specs <- list(
    list(membership = "is_neuronal", high = "ctrl_top_quartile",
         id = "neuronal_ctrl_top_quartile", label = "Ctrl top quartile",
         desc = paste("Derived neuronal genes are enriched among genes in the",
                      "top quartile of control H2AK119ub gene-body signal.")),
    list(membership = "is_neuronal", high = "ctrl_top_decile",
         id = "neuronal_ctrl_top_decile", label = "Ctrl top decile",
         desc = paste("Derived neuronal genes are enriched among genes in the",
                      "top decile of control H2AK119ub gene-body signal.")),
    list(membership = "is_neuronal", high = "mut_top_quartile",
         id = "neuronal_mut_top_quartile", label = "Mut top quartile",
         desc = paste("Derived neuronal genes are enriched among genes in the",
                      "top quartile of mutant H2AK119ub gene-body signal.")),
    list(membership = "is_neuronal", high = "mut_top_decile",
         id = "neuronal_mut_top_decile", label = "Mut top decile",
         desc = paste("Derived neuronal genes are enriched among genes in the",
                      "top decile of mutant H2AK119ub gene-body signal.")),
    list(membership = "is_neuronal_external", high = "ctrl_top_quartile",
         id = "external_ctrl_top_quartile", label = "Ctrl top quartile (external set)",
         desc = paste("External neuronal genes are enriched among genes in the",
                      "top quartile of control H2AK119ub gene-body signal."))
  )

  rows <- lapply(specs, function(s) {
    test_high_signal_enrichment(universe, s$membership, s$high, s$id, s$label,
                                s$desc, out_dir)
  })
  results <- do.call(rbind, rows)
  results$q_value <- p.adjust(results$p_value, method = "BH")
  results$log2_odds_ratio <- log2(results$odds_ratio)
  write_section_table(results, file.path(out_dir, "70_01c_threshold_fisher.tsv"))

  plot_df <- results
  plot_df$test <- factor(plot_df$test, levels = rev(plot_df$test))
  plot_df$point_label <- sprintf("OR = %.2f, %s (%s of %s neuronal)",
                                 plot_df$odds_ratio, format_p(plot_df$p_value),
                                 format(plot_df$n_member_high, big.mark = ",",
                                        trim = TRUE),
                                 format(plot_df$n_high, big.mark = ",",
                                        trim = TRUE))
  plot_df$significance <- ifelse(plot_df$q_value < Q_THRESHOLD,
                                 "Significant", "Not Significant")

  figure <- ggplot(plot_df, aes(x = log2_odds_ratio, y = test)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_errorbar(aes(xmin = log2(ci_lower), xmax = log2(ci_upper)),
                  orientation = "y", width = 0.2, linewidth = 0.8) +
    geom_point(aes(color = significance), size = 4) +
    geom_text(aes(label = point_label), vjust = -1.2, size = 3, color = "grey20") +
    scale_color_manual(values = COLORS$significant, name = "") +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.35))) +
    scale_y_discrete(expand = expansion(add = 0.7)) +
    labs(
      title = "Neuronal gene enrichment among K119ub-high genes",
      subtitle = sprintf(paste("Universe: %s quantifiable genes,",
                               "%.1f%% neuronal in the derived set"),
                         format(plot_df$n_genes[1], big.mark = ","),
                         100 * plot_df$frac_member_universe[1]),
      x = expression(log[2](Odds~Ratio)), y = ""
    ) +
    theme_emseq() +
    theme(legend.position = "top",
          plot.subtitle = element_text(size = 9, color = "grey40"))

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "70_01c_neuronal_enrichment_fisher"),
                          width = 12, height = 7)
  cat("\n")
  invisible(results)
}

# =============================================================================
# FIGURE 70_01d: DOSE-RESPONSE ACROSS K119ub SIGNAL DECILES
# =============================================================================

#' Register one Fisher test comparing one decile against the other nine.
test_one_decile <- function(universe, decile, out_dir) {
  gene_df <- data.frame(
    gene_name = universe$gene_name,
    chr = universe$chr,
    is_neuronal = universe$is_neuronal,
    in_decile = universe$k119ub_decile == decile,
    stringsAsFactors = FALSE
  )

  ft <- register_fisher_test(
    section = SECTION_ID,
    test_id = sprintf("decile_%02d", decile),
    description = sprintf(paste("Derived neuronal genes are enriched in decile",
                                "%d of control H2AK119ub signal against the",
                                "other nine deciles."), decile),
    gene_df = gene_df, row_var = "is_neuronal", col_var = "in_decile",
    output_dir = out_dir
  )

  in_decile <- gene_df$in_decile
  n_total <- sum(in_decile)
  n_neuronal <- sum(in_decile & gene_df$is_neuronal)
  binom_ci <- binom.test(n_neuronal, n_total)$conf.int

  data.frame(
    decile = decile,
    signal_min = min(universe$gb_ctrl_signal[in_decile]),
    signal_max = max(universe$gb_ctrl_signal[in_decile]),
    n_total = n_total,
    n_neuronal = n_neuronal,
    frac_neuronal = n_neuronal / n_total,
    ci_lower = binom_ci[1],
    ci_upper = binom_ci[2],
    odds_ratio = unname(ft$estimate),
    or_ci_lower = ft$conf.int[1],
    or_ci_upper = ft$conf.int[2],
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

plot_decile_dose_response <- function(universe, out_dir) {
  cat("--- Figure 70_01d: neuronal fraction across K119ub deciles ---\n")

  deciles <- do.call(rbind, lapply(seq_len(N_DECILES), function(d) {
    test_one_decile(universe, d, out_dir)
  }))
  deciles$q_value <- p.adjust(deciles$p_value, method = "BH")
  deciles$log2_odds_ratio <- log2(deciles$odds_ratio)

  universe_fraction <- mean(universe$is_neuronal)

  spearman <- cor.test(deciles$decile, deciles$frac_neuronal, method = "spearman",
                       exact = FALSE)
  trend <- prop.trend.test(deciles$n_neuronal, deciles$n_total,
                           score = deciles$decile)

  cat(sprintf("  Universe neuronal fraction: %.1f%%\n", 100 * universe_fraction))
  cat(sprintf("  D1 = %.1f%% neuronal, D10 = %.1f%% neuronal\n",
              100 * deciles$frac_neuronal[1], 100 * deciles$frac_neuronal[N_DECILES]))
  cat(sprintf("  Spearman across deciles: rho = %.3f, %s\n",
              unname(spearman$estimate), format_p(spearman$p.value)))
  cat(sprintf("  Cochran-Armitage trend: chi-squared = %.1f, df = %d, %s\n",
              unname(trend$statistic), unname(trend$parameter),
              format_p(trend$p.value)))

  deciles$spearman_rho <- unname(spearman$estimate)
  deciles$spearman_p <- spearman$p.value
  deciles$trend_chisq <- unname(trend$statistic)
  deciles$trend_df <- unname(trend$parameter)
  deciles$trend_p <- trend$p.value
  write_section_table(deciles, file.path(out_dir, "70_01d_decile_summary.tsv"))

  deciles$decile_label <- factor(sprintf("D%d", deciles$decile),
                                 levels = sprintf("D%d", seq_len(N_DECILES)))
  deciles$bar_label <- sprintf("%.1f%%\nn = %s / %s",
                               100 * deciles$frac_neuronal,
                               format(deciles$n_neuronal, big.mark = ",",
                                      trim = TRUE),
                               format(deciles$n_total, big.mark = ",",
                                      trim = TRUE))

  p_fraction <- ggplot(deciles, aes(x = decile_label, y = frac_neuronal)) +
    geom_col(aes(fill = frac_neuronal), width = 0.75, show.legend = FALSE) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.25,
                  linewidth = 0.5) +
    geom_text(aes(label = bar_label, y = ci_upper), vjust = -0.35, size = 2.9,
              lineheight = 1.05) +
    geom_hline(yintercept = universe_fraction, linetype = "dashed",
               color = COLORS$direction[["Hypermethylated"]], linewidth = 0.8) +
    scale_fill_gradient(low = "grey85", high = GROUP_COLORS[["Neuronal"]]) +
    scale_y_continuous(labels = scales::percent_format(),
                       expand = expansion(mult = c(0, 0.22))) +
    labs(
      title = "Neuronal gene fraction across control K119ub signal deciles",
      subtitle = sprintf(paste("Spearman rho = %.3f (%s) |",
                               "trend chi-squared = %.1f (%s) |",
                               "dashed line = universe fraction %.1f%%"),
                         unname(spearman$estimate), format_p(spearman$p.value),
                         unname(trend$statistic), format_p(trend$p.value),
                         100 * universe_fraction),
      x = "Control K119ub signal decile (D1 lowest, D10 highest)",
      y = "Fraction neuronal"
    ) +
    theme_emseq() +
    theme(plot.subtitle = element_text(size = 8.5, color = "grey40"))

  deciles$significance <- ifelse(deciles$q_value < Q_THRESHOLD,
                                 "Significant", "Not Significant")

  p_odds <- ggplot(deciles, aes(x = decile_label, y = log2_odds_ratio)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_segment(aes(xend = decile_label, y = 0, yend = log2_odds_ratio),
                 color = "grey55") +
    geom_errorbar(aes(ymin = log2(or_ci_lower), ymax = log2(or_ci_upper)),
                  width = 0.2, linewidth = 0.5, color = "grey45") +
    geom_point(aes(color = significance), size = 3) +
    geom_text(aes(label = sprintf("%.2f", odds_ratio)), vjust = -1.3,
              size = 2.9, color = "grey20") +
    scale_color_manual(values = COLORS$significant, name = "") +
    scale_y_continuous(expand = expansion(mult = c(0.18, 0.22))) +
    labs(
      title = "Neuronal enrichment odds ratio per decile against the other nine",
      x = "Control K119ub signal decile", y = expression(log[2](Odds~Ratio))
    ) +
    theme_emseq() +
    theme(legend.position = "top")

  figure <- p_fraction / p_odds +
    plot_layout(heights = c(1, 0.85)) +
    plot_annotation(
      title = "Dose-response between K119ub level and neuronal gene identity",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = 16))
    )

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "70_01d_neuronal_fraction_by_decile"),
                          width = 13, height = 11)
  cat("\n")

  list(deciles = deciles, spearman = spearman, trend = trend,
       universe_fraction = universe_fraction)
}

# =============================================================================
# FIGURE 70_01e: GO ENRICHMENT OF THE K119ub TOP QUARTILE
# =============================================================================

run_go_enrichment <- function(universe, fdr_threshold, out_dir) {
  cat("--- Figure 70_01e: GO BP enrichment of the K119ub top quartile ---\n")

  high <- universe$ctrl_top_quartile
  gene_ids <- as.character(universe$entrez_id[high])
  universe_ids <- as.character(universe$entrez_id)

  cat(sprintf("  Query: %s top-quartile genes | Background: %s genes\n",
              format(length(gene_ids), big.mark = ","),
              format(length(universe_ids), big.mark = ",")))

  ego <- enrichGO(
    gene = gene_ids,
    universe = universe_ids,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    keyType = "ENTREZID",
    pAdjustMethod = "BH",
    pvalueCutoff = fdr_threshold,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  if (is.null(ego)) {
    stop("enrichGO() returned NULL for the ", length(gene_ids),
         " top-quartile K119ub genes. None of their Entrez IDs mapped to a GO ",
         "biological process term in org.Mm.eg.db.")
  }

  full_result <- ego@result
  full_result$is_neuronal_term <- grepl(NEURONAL_PATTERN, full_result$Description,
                                        ignore.case = TRUE)
  write_section_table(full_result,
                  file.path(out_dir, "70_01e_go_bp_top_quartile.tsv"))

  significant <- full_result[!is.na(full_result$p.adjust) &
                               full_result$p.adjust < fdr_threshold, , drop = FALSE]
  if (nrow(significant) == 0) {
    stop("GO biological process enrichment of the ", length(gene_ids),
         " top-quartile K119ub genes returned no term at FDR < ", fdr_threshold,
         ". The full result is in 70_01e_go_bp_top_quartile.tsv.")
  }

  significant <- significant[order(significant$p.adjust), , drop = FALSE]
  n_neuronal <- sum(significant$is_neuronal_term)
  top_terms <- head(significant, DOTPLOT_CATEGORIES)

  cat(sprintf("  %s significant GO BP terms, %s neuronal (%.1f%%)\n",
              format(nrow(significant), big.mark = ","),
              format(n_neuronal, big.mark = ","),
              100 * n_neuronal / nrow(significant)))
  cat(sprintf("  Neuronal among the top %d: %d\n",
              nrow(top_terms), sum(top_terms$is_neuronal_term)))
  cat("  Top 10 terms:\n")
  for (i in seq_len(min(10, nrow(significant)))) {
    cat(sprintf("    %-58s q = %.2e%s\n",
                substr(significant$Description[i], 1, 58),
                significant$p.adjust[i],
                ifelse(significant$is_neuronal_term[i], "  [neuronal]", "")))
  }

  summary_df <- data.frame(
    n_query_genes = length(gene_ids),
    n_background_genes = length(universe_ids),
    n_terms_reported = nrow(full_result),
    n_terms_significant = nrow(significant),
    n_terms_neuronal = n_neuronal,
    frac_terms_neuronal = n_neuronal / nrow(significant),
    n_top_shown = nrow(top_terms),
    n_top_neuronal = sum(top_terms$is_neuronal_term),
    fdr_threshold = fdr_threshold,
    stringsAsFactors = FALSE
  )
  write_section_table(summary_df,
                  file.path(out_dir, "70_01e_go_neuronal_term_summary.tsv"))

  figure <- dotplot(ego, showCategory = min(DOTPLOT_CATEGORIES,
                                            nrow(significant))) +
    labs(
      title = "GO biological process: K119ub top-quartile genes",
      subtitle = sprintf(paste("%s query genes against %s background |",
                               "%s of %s significant terms are neuronal"),
                         format(length(gene_ids), big.mark = ","),
                         format(length(universe_ids), big.mark = ","),
                         format(n_neuronal, big.mark = ","),
                         format(nrow(significant), big.mark = ","))
    ) +
    theme_emseq(base_size = 10) +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"))

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "70_01e_go_bp_top_quartile_dotplot"),
                          width = 13, height = 11)
  cat("\n")

  list(ego = ego, significant = significant,
       n_significant = nrow(significant), n_neuronal = n_neuronal)
}

# =============================================================================
# FIGURES 70_01f AND 70_01g: GSEA ON ABSOLUTE K119ub SIGNAL
# =============================================================================

#' Build a named ranked vector of absolute K119ub signal keyed by Entrez ID.
build_signal_ranking <- function(universe, signal_col) {
  ranked <- setNames(universe[[signal_col]], as.character(universe$entrez_id))
  ranked <- sort(ranked, decreasing = TRUE)
  cat(sprintf("  %s: %s genes ranked, %.3f down to %.3f\n",
              signal_col, format(length(ranked), big.mark = ","),
              max(ranked), min(ranked)))
  ranked
}

#' Run GSEA over genes ranked by absolute signal and count neuronal terms among
#' the leading terms.
run_one_gsea <- function(ranked, condition, fdr_threshold, out_dir) {
  set.seed(1)
  gsea <- gseGO(
    geneList = ranked,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    keyType = "ENTREZID",
    minGSSize = 15,
    maxGSSize = 500,
    pvalueCutoff = fdr_threshold,
    pAdjustMethod = "BH",
    scoreType = "pos",
    seed = TRUE,
    eps = 0,
    verbose = FALSE
  )
  if (is.null(gsea)) {
    stop("gseGO() returned NULL for the ", condition, " K119ub signal ranking ",
         "of ", length(ranked), " genes. None of their Entrez IDs mapped to a ",
         "GO biological process gene set in org.Mm.eg.db.")
  }

  full_result <- gsea@result
  full_result$is_neuronal_term <- grepl(NEURONAL_PATTERN, full_result$Description,
                                        ignore.case = TRUE)
  write_section_table(full_result,
                  file.path(out_dir, sprintf("70_01f_gsea_%s.tsv",
                                             tolower(condition))))

  significant <- full_result[!is.na(full_result$p.adjust) &
                               full_result$p.adjust < fdr_threshold, , drop = FALSE]
  if (nrow(significant) == 0) {
    stop("GSEA on the ", condition, " K119ub signal ranking returned no term ",
         "at FDR < ", fdr_threshold,
         ". The full result is in 70_01f_gsea_", tolower(condition), ".tsv.")
  }

  significant <- significant[order(-significant$NES), , drop = FALSE]
  leading <- head(significant, LEADING_TERM_COUNT)

  cat(sprintf("  %s: %s significant terms, %s neuronal (%.1f%%)\n",
              condition, format(nrow(significant), big.mark = ","),
              format(sum(significant$is_neuronal_term), big.mark = ","),
              100 * mean(significant$is_neuronal_term)))
  cat(sprintf("  %s: %d of the %d leading terms are neuronal (%.1f%%)\n",
              condition, sum(leading$is_neuronal_term), nrow(leading),
              100 * mean(leading$is_neuronal_term)))
  cat("  Top 8 leading terms by NES:\n")
  for (i in seq_len(min(8, nrow(leading)))) {
    cat(sprintf("    %-52s NES = %+.2f  q = %.2e%s\n",
                substr(leading$Description[i], 1, 52), leading$NES[i],
                leading$p.adjust[i],
                ifelse(leading$is_neuronal_term[i], "  [neuronal]", "")))
  }

  list(gsea = gsea, significant = significant, leading = leading,
       condition = condition)
}

plot_gsea_dotplot <- function(gsea_run, out_dir) {
  n_show <- min(DOTPLOT_CATEGORIES, nrow(gsea_run$significant))
  figure <- dotplot(gsea_run$gsea, showCategory = n_show) +
    labs(
      title = sprintf("GSEA on absolute K119ub signal (%s)", gsea_run$condition),
      subtitle = sprintf(paste("Genes ranked by gene-body signal, not by",
                               "fold change | %s significant terms, %d of the",
                               "%d leading terms neuronal"),
                         format(nrow(gsea_run$significant), big.mark = ","),
                         sum(gsea_run$leading$is_neuronal_term),
                         nrow(gsea_run$leading))
    ) +
    theme_emseq(base_size = 10) +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"))

  save_multiformat_ggplot(
    figure,
    file.path(out_dir, sprintf("70_01f_gsea_%s_dotplot",
                               tolower(gsea_run$condition))),
    width = 13, height = 11
  )
  invisible(figure)
}

plot_gsea_running_score <- function(gsea_run, out_dir) {
  neuronal_terms <- gsea_run$significant[gsea_run$significant$is_neuronal_term, ,
                                         drop = FALSE]
  if (nrow(neuronal_terms) == 0) {
    stop("GSEA on the ", gsea_run$condition, " K119ub ranking found ",
         nrow(gsea_run$significant), " significant terms and none of them is ",
         "neuronal, so the running-score figure has no term to draw. ",
         "This is a negative result for the section hypothesis: see ",
         "70_01f_gsea_", tolower(gsea_run$condition), ".tsv.")
  }

  top <- neuronal_terms[which.max(neuronal_terms$NES), ]
  cat(sprintf("  Running score term: %s (NES = %+.2f, q = %.2e)\n",
              top$Description, top$NES, top$p.adjust))

  figure <- gseaplot2(gsea_run$gsea, geneSetID = top$ID,
                      title = sprintf("%s (%s K119ub signal)",
                                      top$Description, gsea_run$condition),
                      pvalue_table = TRUE)

  slug <- gsub("[^A-Za-z0-9]+", "_", substr(top$Description, 1, 40))
  slug <- gsub("^_|_$", "", slug)
  save_multiformat_ggplot(
    figure,
    file.path(out_dir, sprintf("70_01g_gsea_running_score_%s_%s",
                               tolower(gsea_run$condition), slug)),
    width = 11, height = 8
  )
  top
}

run_gsea_analysis <- function(universe, fdr_threshold, out_dir) {
  cat("--- Figures 70_01f and 70_01g: GSEA on absolute K119ub signal ---\n")

  rank_ctrl <- build_signal_ranking(universe, "gb_ctrl_signal")
  rank_mut <- build_signal_ranking(universe, "gb_mut_signal")

  gsea_ctrl <- run_one_gsea(rank_ctrl, "Control", fdr_threshold, out_dir)
  gsea_mut <- run_one_gsea(rank_mut, "Mutant", fdr_threshold, out_dir)

  plot_gsea_dotplot(gsea_ctrl, out_dir)
  plot_gsea_dotplot(gsea_mut, out_dir)
  plot_gsea_running_score(gsea_ctrl, out_dir)
  plot_gsea_running_score(gsea_mut, out_dir)

  summary_df <- do.call(rbind, lapply(list(gsea_ctrl, gsea_mut), function(run) {
    data.frame(
      condition = run$condition,
      n_terms_significant = nrow(run$significant),
      n_terms_neuronal = sum(run$significant$is_neuronal_term),
      frac_terms_neuronal = mean(run$significant$is_neuronal_term),
      n_leading = nrow(run$leading),
      n_leading_neuronal = sum(run$leading$is_neuronal_term),
      frac_leading_neuronal = mean(run$leading$is_neuronal_term),
      top_neuronal_term = run$significant$Description[
        run$significant$is_neuronal_term][1],
      top_neuronal_nes = run$significant$NES[run$significant$is_neuronal_term][1],
      fdr_threshold = fdr_threshold,
      stringsAsFactors = FALSE
    )
  }))
  write_section_table(summary_df,
                  file.path(out_dir, "70_01f_gsea_neuronal_summary.tsv"))
  cat("\n")

  list(ctrl = gsea_ctrl, mut = gsea_mut, summary = summary_df)
}

# =============================================================================
# FIGURE 70_01h: GENE LENGTH, NEURONAL AGAINST OTHER
# =============================================================================

plot_gene_length_comparison <- function(universe, out_dir) {
  cat("--- Figure 70_01h: gene length, neuronal against other ---\n")

  neuronal <- universe$is_neuronal
  length_test <- wilcox.test(universe$width[neuronal],
                             universe$width[!neuronal],
                             alternative = "greater")
  length_signal_cor <- cor.test(universe$width, universe$gb_ctrl_signal,
                                method = "spearman", exact = FALSE)

  cat(sprintf("  Median length: neuronal %s bp, other %s bp, %s\n",
              format(round(median(universe$width[neuronal])), big.mark = ","),
              format(round(median(universe$width[!neuronal])), big.mark = ","),
              format_p(length_test$p.value)))
  cat(sprintf("  Gene length against ctrl K119ub signal: rho = %.3f, %s\n",
              unname(length_signal_cor$estimate),
              format_p(length_signal_cor$p.value)))

  labels <- summarise_groups(universe, "group", "width")
  labels$label <- group_label(labels, digits = 0)
  labels$label_y <- max(universe$width) * 2.2

  figure <- ggplot(universe, aes(x = group, y = width, fill = group)) +
    geom_violin(scale = "width", alpha = 0.65, linewidth = 0.3,
                show.legend = FALSE) +
    geom_boxplot(width = 0.14, outlier.size = 0.2, outlier.alpha = 0.2,
                 show.legend = FALSE) +
    geom_text(data = labels, aes(x = group, y = label_y, label = label),
              inherit.aes = FALSE, size = 3.4, lineheight = 1.1) +
    scale_fill_manual(values = GROUP_COLORS) +
    scale_y_log10(labels = scales::comma,
                  expand = expansion(mult = c(0.05, 0.2))) +
    labs(
      title = "Gene length of neuronal and other genes",
      subtitle = sprintf(paste("One-sided Wilcoxon, neuronal longer: %s |",
                               "length against ctrl K119ub signal:",
                               "Spearman rho = %.3f"),
                         format_p(length_test$p.value),
                         unname(length_signal_cor$estimate)),
      x = "", y = "Gene length in bp (log10 scale)"
    ) +
    theme_emseq() +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"))

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "70_01h_gene_length_neuronal_vs_other"),
                          width = 9, height = 7)

  stats_df <- labels[, c("group", "n", "median", "mean")]
  colnames(stats_df) <- c("group", "n_genes", "median_width_bp", "mean_width_bp")
  stats_df$wilcoxon_p <- length_test$p.value
  stats_df$wilcoxon_alternative <- "neuronal longer"
  stats_df$length_vs_ctrl_signal_rho <- unname(length_signal_cor$estimate)
  stats_df$length_vs_ctrl_signal_p <- length_signal_cor$p.value
  write_section_table(stats_df,
                  file.path(out_dir, "70_01h_gene_length_comparison.tsv"))
  cat("\n")

  list(wilcoxon = length_test, correlation = length_signal_cor, stats = stats_df)
}

# =============================================================================
# FIGURE 70_01i: K119ub COMPARISON INSIDE GENE-LENGTH STRATA
# =============================================================================

#' Test neuronal against other genes inside one gene-length stratum, using a
#' top-quartile threshold computed within that stratum.
test_one_length_stratum <- function(universe, stratum, out_dir) {
  in_stratum <- universe$length_stratum == stratum
  df <- universe[in_stratum, , drop = FALSE]

  if (length(unique(df$is_neuronal)) < 2) {
    stop("Gene-length stratum ", stratum, " holds only ",
         ifelse(all(df$is_neuronal), "neuronal", "non-neuronal"),
         " genes, so the within-stratum comparison has one class.")
  }

  threshold <- quantile(df$gb_ctrl_signal, 0.75)
  df$stratum_top_quartile <- df$gb_ctrl_signal >= threshold

  gene_df <- data.frame(
    gene_name = df$gene_name,
    chr = df$chr,
    is_neuronal = df$is_neuronal,
    stratum_top_quartile = df$stratum_top_quartile,
    stringsAsFactors = FALSE
  )

  ft <- register_fisher_test(
    section = SECTION_ID,
    test_id = sprintf("length_stratum_%d", stratum),
    description = sprintf(paste("Within gene-length stratum %d of %d, derived",
                                "neuronal genes are enriched among the top",
                                "quartile of control H2AK119ub signal of that",
                                "stratum."), stratum, N_LENGTH_STRATA),
    gene_df = gene_df, row_var = "is_neuronal", col_var = "stratum_top_quartile",
    output_dir = out_dir
  )

  wilcox <- wilcox.test(df$gb_ctrl_signal[df$is_neuronal],
                        df$gb_ctrl_signal[!df$is_neuronal],
                        alternative = "greater")

  data.frame(
    stratum = stratum,
    width_min = min(df$width),
    width_max = max(df$width),
    n_genes = nrow(df),
    n_neuronal = sum(df$is_neuronal),
    frac_neuronal = mean(df$is_neuronal),
    top_quartile_threshold = unname(threshold),
    median_signal_neuronal = median(df$gb_ctrl_signal[df$is_neuronal]),
    median_signal_other = median(df$gb_ctrl_signal[!df$is_neuronal]),
    wilcoxon_p = wilcox$p.value,
    odds_ratio = unname(ft$estimate),
    or_ci_lower = ft$conf.int[1],
    or_ci_upper = ft$conf.int[2],
    fisher_p = ft$p.value,
    stringsAsFactors = FALSE
  )
}

plot_length_stratified <- function(universe, out_dir) {
  cat("--- Figure 70_01i: K119ub comparison inside gene-length strata ---\n")

  strata <- do.call(rbind, lapply(seq_len(N_LENGTH_STRATA), function(s) {
    test_one_length_stratum(universe, s, out_dir)
  }))
  strata$wilcoxon_q <- p.adjust(strata$wilcoxon_p, method = "BH")
  strata$fisher_q <- p.adjust(strata$fisher_p, method = "BH")
  strata$log2_odds_ratio <- log2(strata$odds_ratio)
  write_section_table(strata, file.path(out_dir, "70_01i_length_stratified.tsv"))

  cat("  Within-stratum results:\n")
  for (i in seq_len(nrow(strata))) {
    r <- strata[i, ]
    cat(sprintf("    L%d [%s-%s bp]: neuronal median %.3f against other %.3f, ",
                r$stratum, format(r$width_min, big.mark = ","),
                format(r$width_max, big.mark = ","),
                r$median_signal_neuronal, r$median_signal_other))
    cat(sprintf("OR = %.2f, %s\n", r$odds_ratio, format_p(r$fisher_p)))
  }

  stratum_names <- sprintf("L%d\n%s-%s kb", strata$stratum,
                           format(round(strata$width_min / 1000, 1), trim = TRUE),
                           format(round(strata$width_max / 1000, 1), trim = TRUE))
  names(stratum_names) <- as.character(strata$stratum)

  plot_df <- universe[, c("gene_name", "group", "gb_ctrl_signal",
                          "length_stratum")]
  plot_df$stratum_label <- factor(stratum_names[as.character(plot_df$length_stratum)],
                                  levels = stratum_names)

  labels <- summarise_groups_by(plot_df, "stratum_label", "group",
                                "gb_ctrl_signal")
  labels$label <- group_label(labels, digits = 3)
  labels$label_y <- max(plot_df$gb_ctrl_signal) * 1.9

  p_violin <- ggplot(plot_df, aes(x = group, y = gb_ctrl_signal, fill = group)) +
    geom_violin(scale = "width", alpha = 0.65, linewidth = 0.3,
                show.legend = FALSE) +
    geom_boxplot(width = 0.14, outlier.size = 0.15, outlier.alpha = 0.15,
                 show.legend = FALSE) +
    geom_text(data = labels, aes(x = group, y = label_y, label = label),
              inherit.aes = FALSE, size = 2.6, lineheight = 1.05) +
    facet_wrap(~ stratum_label, nrow = 1) +
    scale_fill_manual(values = GROUP_COLORS) +
    scale_y_log10(expand = expansion(mult = c(0.05, 0.2))) +
    labs(
      title = "Control K119ub signal inside gene-length strata",
      subtitle = "Equal-count strata of gene length, so length cannot drive a within-panel difference",
      x = "", y = "K119ub gene-body signal (log10 scale)"
    ) +
    theme_emseq(base_size = 11) +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"),
          axis.text.x = element_text(size = 8))

  strata$stratum_label <- factor(stratum_names, levels = stratum_names)
  strata$significance <- ifelse(strata$fisher_q < Q_THRESHOLD,
                                "Significant", "Not Significant")

  p_odds <- ggplot(strata, aes(x = stratum_label, y = log2_odds_ratio)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_errorbar(aes(ymin = log2(or_ci_lower), ymax = log2(or_ci_upper)),
                  width = 0.2, linewidth = 0.6, color = "grey45") +
    geom_point(aes(color = significance), size = 3.5) +
    geom_text(aes(label = sprintf("OR = %.2f\n%s", odds_ratio,
                                  format_p(fisher_p))),
              vjust = -0.9, size = 2.8, lineheight = 1.05, color = "grey20") +
    scale_color_manual(values = COLORS$significant, name = "") +
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.35))) +
    labs(
      title = "Neuronal enrichment in the within-stratum K119ub top quartile",
      x = "Gene-length stratum", y = expression(log[2](Odds~Ratio))
    ) +
    theme_emseq(base_size = 11) +
    theme(legend.position = "top", axis.text.x = element_text(size = 8))

  figure <- p_violin / p_odds +
    plot_layout(heights = c(1, 0.8)) +
    plot_annotation(
      title = "K119ub enrichment at neuronal genes holds inside gene-length strata",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = 16))
    )

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "70_01i_length_stratified_k119ub"),
                          width = 16, height = 12)
  cat("\n")
  invisible(strata)
}

# =============================================================================
# FIGURE 70_01j: CONSTITUTIVE SIGNAL AGAINST FOLD CHANGE
# =============================================================================

plot_signal_against_log2fc <- function(universe, out_dir) {
  cat("--- Figure 70_01j: constitutive signal against K119ub log2FC ---\n")

  correlation <- cor.test(universe$gb_ctrl_signal, universe$gb_log2fc,
                          method = "spearman", exact = FALSE)
  key_df <- universe[universe$gene_name %in% KEY_GENES, , drop = FALSE]
  cat(sprintf("  Ctrl signal against log2FC: rho = %.3f, %s\n",
              unname(correlation$estimate), format_p(correlation$p.value)))
  cat(sprintf("  Key genes in the universe: %d of %d (%s)\n",
              nrow(key_df), length(KEY_GENES),
              paste(key_df$gene_name, collapse = ", ")))

  figure <- ggplot(universe, aes(x = gb_ctrl_signal, y = gb_log2fc)) +
    geom_point(aes(color = group), size = 0.4, alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55",
               linewidth = 0.3) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "black",
                linewidth = 0.8) +
    scale_color_manual(values = GROUP_COLORS, name = "",
                       guide = guide_legend(override.aes = list(size = 3,
                                                                alpha = 1))) +
    scale_x_log10() +
    labs(
      title = "Constitutive K119ub level against the BAP1-driven change",
      subtitle = sprintf("Spearman rho = %.3f, %s, n = %s genes",
                         unname(correlation$estimate),
                         format_p(correlation$p.value),
                         format(nrow(universe), big.mark = ",")),
      x = "Control K119ub gene-body signal (log10 scale)",
      y = "K119ub log2FC (mut / ctrl)"
    ) +
    theme_emseq() +
    theme(legend.position = "top",
          plot.subtitle = element_text(size = 9, color = "grey40"))

  if (nrow(key_df) > 0) {
    figure <- figure +
      geom_point(data = key_df, size = 1.8, color = "black") +
      ggrepel::geom_text_repel(data = key_df, aes(label = gene_name),
                               size = 3, min.segment.length = 0,
                               max.overlaps = Inf, color = "black")
  }

  save_multiformat_ggplot(figure,
                          file.path(out_dir, "70_01j_signal_vs_log2fc"),
                          width = 10, height = 8)

  key_out <- key_df[, c("gene_name", "chr", "width", "gb_ctrl_signal",
                        "gb_mut_signal", "gb_log2fc", "is_neuronal",
                        "k119ub_decile")]
  key_out$spearman_rho_universe <- unname(correlation$estimate)
  key_out$spearman_p_universe <- correlation$p.value
  write_section_table(key_out, file.path(out_dir, "70_01j_key_genes.tsv"))
  cat("\n")
  invisible(correlation)
}

# =============================================================================
# HANDOFF AND PER-GENE TABLES
# =============================================================================

#' Write the union neuronal gene set with its K119ub level, in the exact column
#' order that sections 70_02, 70_03 and 70_04 read.
#'
#' The K119ub columns come from the quantifiable universe, which is the same
#' set every analysis in this section uses. A neuronal gene outside that set
#' carries NA in all three K119ub columns.
write_neuronal_handoff <- function(derived, external, universe, out_dir) {
  cat("--- Writing the neuronal gene set handoff ---\n")

  union_genes <- sort(union(derived, external))
  handoff <- data.frame(gene = union_genes, stringsAsFactors = FALSE)
  handoff$is_derived <- handoff$gene %in% derived
  handoff$is_external <- handoff$gene %in% external
  handoff$source <- ifelse(handoff$is_derived & handoff$is_external, "both",
                           ifelse(handoff$is_derived, "derived", "external"))

  signal_lookup <- universe[, c("gene_name", "gb_ctrl_signal", "gb_log2fc",
                                "k119ub_decile")]
  colnames(signal_lookup) <- c("gene", "k119ub_ctrl_signal", "k119ub_log2fc",
                               "k119ub_decile")
  handoff <- dplyr::left_join(handoff, signal_lookup, by = "gene")

  if (nrow(handoff) != length(union_genes)) {
    stop("The signal join changed the handoff row count: ", length(union_genes),
         " genes in, ", nrow(handoff), " out.")
  }
  missing <- setdiff(HANDOFF_COLUMNS, colnames(handoff))
  if (length(missing) > 0) {
    stop("The handoff table is missing contract columns: ",
         paste(missing, collapse = ", "))
  }

  handoff <- handoff[, HANDOFF_COLUMNS]
  n_with_signal <- sum(!is.na(handoff$k119ub_ctrl_signal))
  cat(sprintf("  Union set: %s genes (derived %s, external %s, both %s)\n",
              format(nrow(handoff), big.mark = ","),
              format(sum(handoff$source == "derived"), big.mark = ","),
              format(sum(handoff$source == "external"), big.mark = ","),
              format(sum(handoff$source == "both"), big.mark = ",")))
  cat(sprintf("  With a K119ub level: %s (%.1f%%); the rest carry NA\n",
              format(n_with_signal, big.mark = ","),
              100 * n_with_signal / nrow(handoff)))

  path <- file.path(out_dir, basename(HANDOFF_PATHS$neuronal_gene_set))
  write_section_table(handoff, path)

  if (normalizePath(out_dir, mustWork = FALSE) !=
      normalizePath(OUTPUT_PATHS$neuronal, mustWork = FALSE)) {
    cat(sprintf("  NOTE: sections 70_02, 70_03 and 70_04 read %s\n",
                HANDOFF_PATHS$neuronal_gene_set))
  }
  cat("\n")
  path
}

write_per_gene_table <- function(universe, out_dir) {
  cols <- c("gene_name", "entrez_id", "chr", "start", "end", "width",
            "gb_ctrl_signal", "gb_mut_signal", "gb_log2fc",
            "pr_ctrl_signal", "pr_mut_signal", "pr_log2fc",
            "is_neuronal", "is_neuronal_external",
            "k119ub_decile", "length_stratum",
            "ctrl_top_quartile", "ctrl_top_decile",
            "mut_top_quartile", "mut_top_decile")
  out <- universe[order(-universe$gb_ctrl_signal), cols]
  write_section_table(out,
                  file.path(out_dir, "70_01_per_gene_k119ub_neuronal.tsv"))
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
  cat("  SECTION 70_01: H2AK119ub AT NEURONAL GENES\n")
  cat("===========================================================================\n")
  cat("  Output directory: ", out_dir, "\n", sep = "")
  cat("  GO and GSEA FDR threshold: ", fdr_threshold, "\n", sep = "")
  cat("  Signal used: absolute gene-body K119ub, not fold change\n\n")

  quantifiable <- load_k119ub_signal()

  derived_genes <- derive_neuronal_gene_set(out_dir)
  cat("--- Loading the external neuronal gene set ---\n")
  external_genes <- load_gene_set(GENESET_PATHS$neuronal,
                                  "External neuronal gene set")
  cat("\n")

  overlap <- report_gene_set_overlap(derived_genes, external_genes,
                                     quantifiable, out_dir)

  universe <- annotate_universe(quantifiable, derived_genes, external_genes,
                                out_dir)

  distribution <- plot_signal_distribution(universe, out_dir)
  enrichment <- plot_high_signal_enrichment(universe, out_dir)
  dose_response <- plot_decile_dose_response(universe, out_dir)
  go_result <- run_go_enrichment(universe, fdr_threshold, out_dir)
  gsea_result <- run_gsea_analysis(universe, fdr_threshold, out_dir)
  length_result <- plot_gene_length_comparison(universe, out_dir)
  strata <- plot_length_stratified(universe, out_dir)
  plot_signal_against_log2fc(universe, out_dir)

  handoff_path <- write_neuronal_handoff(derived_genes, external_genes,
                                         universe, out_dir)
  write_per_gene_table(universe, out_dir)

  cat("\n---------------------------------------------------------------------------\n")
  cat("  SUMMARY\n")
  cat("---------------------------------------------------------------------------\n")
  cat(sprintf("  Universe: %s quantifiable genes, %s neuronal (%.1f%%)\n",
              format(nrow(universe), big.mark = ","),
              format(sum(universe$is_neuronal), big.mark = ","),
              100 * mean(universe$is_neuronal)))
  cat(sprintf("  Gene sets: derived %s, external %s, both %s, Jaccard %.3f\n",
              format(length(derived_genes), big.mark = ","),
              format(length(external_genes), big.mark = ","),
              format(overlap$n_both, big.mark = ","), overlap$jaccard))
  cat(sprintf("  K119ub signal, control: neuronal median %.3f against other %.3f, %s\n",
              distribution$stats$median_signal[
                distribution$stats$condition == "Control" &
                  distribution$stats$group == "Neuronal"],
              distribution$stats$median_signal[
                distribution$stats$condition == "Control" &
                  distribution$stats$group == "Other"],
              format_p(distribution$ctrl_p)))
  cat("  Fisher odds ratios:\n")
  for (i in seq_len(nrow(enrichment))) {
    cat(sprintf("    %-32s OR = %.2f [%.2f, %.2f], %s\n",
                enrichment$test[i], enrichment$odds_ratio[i],
                enrichment$ci_lower[i], enrichment$ci_upper[i],
                format_p(enrichment$p_value[i])))
  }
  cat(sprintf("  Decile dose-response: D1 %.1f%% to D10 %.1f%% neuronal, ",
              100 * dose_response$deciles$frac_neuronal[1],
              100 * dose_response$deciles$frac_neuronal[N_DECILES]))
  cat(sprintf("Spearman rho = %.3f, trend %s\n",
              unname(dose_response$spearman$estimate),
              format_p(dose_response$trend$p.value)))
  cat(sprintf("  GO BP: %s significant terms, %s neuronal\n",
              format(go_result$n_significant, big.mark = ","),
              format(go_result$n_neuronal, big.mark = ",")))
  for (i in seq_len(nrow(gsea_result$summary))) {
    r <- gsea_result$summary[i, ]
    cat(sprintf("  GSEA %s: %s significant terms, %d of the %d leading terms neuronal\n",
                r$condition, format(r$n_terms_significant, big.mark = ","),
                r$n_leading_neuronal, r$n_leading))
  }
  cat(sprintf("  Gene length: neuronal median %s bp against other %s bp, %s\n",
              format(round(length_result$stats$median_width_bp[
                length_result$stats$group == "Neuronal"]), big.mark = ","),
              format(round(length_result$stats$median_width_bp[
                length_result$stats$group == "Other"]), big.mark = ","),
              format_p(length_result$wilcoxon$p.value)))
  cat(sprintf("  Length-stratified odds ratios: %s\n",
              paste(sprintf("L%d %.2f", strata$stratum, strata$odds_ratio),
                    collapse = ", ")))
  cat(sprintf("  Handoff: %s\n", handoff_path))
  cat("\n=== Section 70_01 complete ===\n\n")
}

main()
