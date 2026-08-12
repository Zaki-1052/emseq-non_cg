# scripts/sections/60_mecp2/60_04_aging.R
#
# Section 60_04: MeCP2 developmental trajectory, young against adult.
#
# What this tests
#   MeCP2 occupancy rises across postnatal neuronal maturation. Two DiffBind
#   comparisons measure that rise separately in each genotype: adult against
#   young in control, and adult against young in mutant. This section asks
#   whether the mutant gains MeCP2 with age at loci the control does not, and
#   whether those loci carry the adult mCH change measured by this pipeline.
#
# Analyses
#   1. Counts of aging-up and aging-down peaks in each genotype, with a
#      peak-level Fisher test of the aging-up proportion between genotypes.
#   2. Gene-level overlap of the two aging-up gene sets, drawn as a Venn.
#   3. Peaks shared by both genotypes: Spearman correlation and paired Wilcoxon
#      test of the control aging fold against the mutant aging fold.
#   4. Mutant-specific aging-up genes: significant up in mutant, not significant
#      in control.
#   5. GO biological process enrichment on the mutant-specific genes, with all
#      peak-carrying genes as the universe.
#   6. Three registered gene-level Fisher tests of mutant-specific aging genes
#      against the adult mCH result, plus a Wilcoxon of mch_diff between the
#      mutant-specific genes and the rest.
#
# Reads
#   MECP2_PATHS$ctrl_aging   adult-vs-young MeCP2 DiffBind table, control
#   MECP2_PATHS$mut_aging    adult-vs-young MeCP2 DiffBind table, mutant
#   mch_results              gene-level adult mCH differential results (config)
#
# Writes
#   Five multi-format figures and thirteen TSV tables into OUTPUT_PATHS$mecp2
#   (override with --output-dir), plus three registered Fisher gene tables
#   under fisher_tables/ and three rows in HANDOFF_PATHS$fisher_registry.
#
# Adapted from Biomodal section 77 (MeCP2 aging trajectory, young vs adult).

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

SECTION_ID <- "60_04"

# Aging direction labels. "Aging up" means higher MeCP2 in adult than in young.
AGING_UP     <- "Aging up"
AGING_DOWN   <- "Aging down"
AGING_NONSIG <- "Not significant"

AGING_STATUS_ORDER <- c(AGING_UP, AGING_DOWN, AGING_NONSIG)

AGING_STATUS_COLORS <- c(
  "Aging up"        = unname(COLORS$mecp2["MeCP2 Up"]),
  "Aging down"      = unname(COLORS$mecp2["MeCP2 Down"]),
  "Not significant" = "grey70"
)

GENOTYPE_ORDER <- c("Control", "Mutant")

# Gene groups for the mCH comparison.
GROUP_MUT_SPECIFIC <- "Mutant-specific aging up"
GROUP_OTHER        <- "Other peak-carrying genes"

AGING_GROUP_ORDER <- c(GROUP_MUT_SPECIFIC, GROUP_OTHER)

AGING_GROUP_COLORS <- c(
  "Mutant-specific aging up"  = "#D95F02",
  "Other peak-carrying genes" = "grey70"
)

# Columns kept from each 53 MB aging DiffBind table.
AGING_KEEP_COLS <- c("Chr", "Start", "End", "Conc", "Conc_young", "Conc_adult",
                     "Fold", "p.value", "FDR")

# Columns kept from each ChIPseeker annotation table.
ANNOTATED_KEEP_COLS <- c("seqnames", "start", "end", "Fold", "FDR", "SYMBOL")

# Lowest correlation accepted between the negated Fold and Conc_adult minus
# Conc_young. A sign error in the file layout drives this strongly negative.
MIN_ORIENTATION_COR <- 0.5

# =============================================================================
# COMMAND LINE
# =============================================================================

parse_options <- function() {
  option_list <- list(
    make_option("--output-dir", type = "character", default = OUTPUT_PATHS$mecp2,
                dest = "output_dir",
                help = "Directory for figures and tables [default: %default]"),
    make_option("--fdr-threshold", type = "double", default = Q_THRESHOLD,
                dest = "fdr_threshold",
                help = "FDR cutoff for an aging peak [default: %default]"),
    make_option("--go-top-terms", type = "integer", default = 25,
                dest = "go_top_terms",
                help = "GO BP terms drawn in the dot plot [default: %default]")
  )
  parse_args(OptionParser(option_list = option_list))
}

# =============================================================================
# SMALL UTILITIES
# =============================================================================

#' Format a p-value for a figure subtitle.
fmt_p <- function(p) {
  if (length(p) != 1 || is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.3g", p)
}

# =============================================================================
# INPUT LOADING
# =============================================================================

#' Read one adult-versus-young MeCP2 DiffBind table.
#'
#' The file reports Fold as log2(young / adult), together with Conc_young and
#' Conc_adult. This function negates Fold so that a positive value means more
#' MeCP2 in adult than in young, then checks the negated Fold against
#' Conc_adult minus Conc_young and stops when the two disagree in sign.
#'
#' @param filepath Path to the DiffBind table.
#' @param label Display name used in messages.
#' @param fdr_threshold FDR cutoff for the aging status.
#' @return data.frame with AGING_KEEP_COLS plus aging_status and direction
load_aging_diffbind <- function(filepath, label, fdr_threshold) {
  db <- load_diffbind_flex(filepath, label, fdr_threshold = fdr_threshold)
  cat(sprintf("  %s: the counts above use the file's log2(young/adult) Fold.\n",
              label))

  missing <- setdiff(AGING_KEEP_COLS, colnames(db))
  if (length(missing) > 0) {
    stop(label, " aging table is missing columns: ",
         paste(missing, collapse = ", "), " (", filepath, ")")
  }

  db <- db[, AGING_KEEP_COLS, drop = FALSE]
  db$Fold <- -db$Fold

  conc_difference <- db$Conc_adult - db$Conc_young
  orientation_cor <- cor(db$Fold, conc_difference, use = "complete.obs")
  if (!is.finite(orientation_cor) || orientation_cor < MIN_ORIENTATION_COR) {
    stop(label, ": the negated Fold correlates with Conc_adult - Conc_young at ",
         sprintf("%.3f", orientation_cor), ", below the required ",
         MIN_ORIENTATION_COR, ". The file's fold orientation is not ",
         "log2(young/adult) as this section assumes. File: ", filepath)
  }

  db$aging_status <- AGING_NONSIG
  db$aging_status[db$FDR < fdr_threshold & db$Fold > 0] <- AGING_UP
  db$aging_status[db$FDR < fdr_threshold & db$Fold < 0] <- AGING_DOWN
  db$direction <- db$aging_status

  cat(sprintf("  %s aging fold (adult / young): %s peaks, %s up, %s down, %s not significant\n",
              label, format(nrow(db), big.mark = ","),
              format(sum(db$aging_status == AGING_UP), big.mark = ","),
              format(sum(db$aging_status == AGING_DOWN), big.mark = ","),
              format(sum(db$aging_status == AGING_NONSIG), big.mark = ",")))
  cat(sprintf("  %s orientation check: cor(Fold, Conc_adult - Conc_young) = %.3f\n",
              label, orientation_cor))
  db
}

#' Annotate an aging table to genes and keep the columns this section uses.
#'
#' @param db data.frame from load_aging_diffbind().
#' @param label Display name used in messages.
#' @return data.frame with ANNOTATED_KEEP_COLS
annotate_aging_peaks <- function(db, label) {
  annotated <- annotate_peaks_to_genes(db, label)

  missing <- setdiff(ANNOTATED_KEEP_COLS, colnames(annotated))
  if (length(missing) > 0) {
    stop(label, " annotation is missing columns: ",
         paste(missing, collapse = ", "))
  }

  annotated <- annotated[, ANNOTATED_KEEP_COLS, drop = FALSE]
  n_with_symbol <- sum(!is.na(annotated$SYMBOL) & nzchar(annotated$SYMBOL))
  cat(sprintf("  %s: %s of %s annotated peaks carry a gene symbol\n",
              label, format(n_with_symbol, big.mark = ","),
              format(nrow(annotated), big.mark = ",")))

  annotated[!is.na(annotated$SYMBOL) & nzchar(annotated$SYMBOL), , drop = FALSE]
}

# =============================================================================
# PEAK-LEVEL COUNTS
# =============================================================================

#' Count peaks by aging status in both genotypes.
#'
#' @param ctrl_db Control aging table.
#' @param mut_db Mutant aging table.
#' @return data.frame with genotype, aging_status, n_peaks, pct_of_peaks
count_aging_peaks <- function(ctrl_db, mut_db) {
  count_one <- function(db, genotype) {
    data.frame(
      genotype = genotype,
      aging_status = AGING_STATUS_ORDER,
      n_peaks = as.integer(table(factor(db$aging_status,
                                        levels = AGING_STATUS_ORDER))),
      n_peaks_total = nrow(db),
      stringsAsFactors = FALSE
    )
  }

  out <- rbind(count_one(ctrl_db, "Control"), count_one(mut_db, "Mutant"))
  out$pct_of_peaks <- 100 * out$n_peaks / out$n_peaks_total
  out
}

#' Peak-level Fisher test of the aging-up proportion between genotypes.
#'
#' Each peak is one observation. This is a peak-level test, so it calls
#' fisher.test() directly rather than register_fisher_test().
#'
#' @param ctrl_db Control aging table.
#' @param mut_db Mutant aging table.
#' @param status Aging status counted as the event.
#' @return data.frame with counts, odds ratio, and p-value
aging_proportion_fisher <- function(ctrl_db, mut_db, status) {
  n_mut_event <- sum(mut_db$aging_status == status)
  n_ctrl_event <- sum(ctrl_db$aging_status == status)

  tab <- matrix(
    c(n_mut_event, nrow(mut_db) - n_mut_event,
      n_ctrl_event, nrow(ctrl_db) - n_ctrl_event),
    nrow = 2,
    dimnames = list(c("event", "other"), c("mutant", "control"))
  )
  ft <- fisher.test(tab)

  data.frame(
    aging_status = status,
    n_ctrl_peaks = nrow(ctrl_db),
    n_ctrl_event = n_ctrl_event,
    pct_ctrl_event = 100 * n_ctrl_event / nrow(ctrl_db),
    n_mut_peaks = nrow(mut_db),
    n_mut_event = n_mut_event,
    pct_mut_event = 100 * n_mut_event / nrow(mut_db),
    odds_ratio = unname(ft$estimate),
    conf_low = ft$conf.int[1],
    conf_high = ft$conf.int[2],
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

#' Grouped bar chart of aging peak counts.
#'
#' The not-significant peaks outnumber the significant ones by two orders of
#' magnitude, so the count axis carries a square-root transform and every bar
#' also prints its exact count and percentage.
#'
#' @param count_df data.frame from count_aging_peaks().
#' @param subtitle Figure subtitle.
#' @return ggplot object
plot_aging_peak_counts <- function(count_df, subtitle) {
  plot_df <- count_df
  plot_df$genotype <- factor(plot_df$genotype, levels = GENOTYPE_ORDER)
  plot_df$aging_status <- factor(plot_df$aging_status, levels = AGING_STATUS_ORDER)
  plot_df$bar_label <- sprintf("%s\n(%.2f%%)",
                               format(plot_df$n_peaks, big.mark = ",", trim = TRUE),
                               plot_df$pct_of_peaks)

  ggplot(plot_df, aes(x = genotype, y = n_peaks, fill = aging_status)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.7,
             color = "black", linewidth = 0.3) +
    geom_text(aes(label = bar_label), position = position_dodge(width = 0.78),
              vjust = -0.3, size = 3.1, lineheight = 1.05) +
    scale_fill_manual(values = AGING_STATUS_COLORS, name = "MeCP2 aging status") +
    scale_y_sqrt(labels = scales::comma,
                 expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = "MeCP2 binding change from young to adult",
      subtitle = subtitle,
      x = "Genotype",
      y = "MeCP2 peaks (square-root axis)"
    ) +
    theme_emseq() +
    theme(legend.position = "top")
}

# =============================================================================
# GENE-LEVEL AGING TABLES
# =============================================================================

#' Count significant aging peaks per gene in each direction.
#'
#' aggregate_diffbind_by_gene() reports one fold and a total significant count
#' per gene. This adds the split of that count into up and down peaks, plus the
#' extreme folds, which the aging-up gene sets are built from.
#'
#' @param annotated data.frame from annotate_aging_peaks().
#' @param fdr_threshold FDR cutoff for a significant peak.
#' @param prefix Column-name prefix.
#' @return data.frame with one row per gene
count_gene_aging_peaks <- function(annotated, fdr_threshold, prefix) {
  out <- annotated %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(
      n_sig_up = sum(FDR < fdr_threshold & Fold > 0, na.rm = TRUE),
      n_sig_down = sum(FDR < fdr_threshold & Fold < 0, na.rm = TRUE),
      max_fold = max(Fold, na.rm = TRUE),
      min_fold = min(Fold, na.rm = TRUE),
      mean_fold = mean(Fold, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.frame()

  colnames(out) <- c("gene_name",
                     paste0(prefix, c("_n_sig_up", "_n_sig_down", "_max_fold",
                                      "_min_fold", "_mean_fold")))
  out
}

#' Build the one-row-per-gene aging table for one genotype.
#'
#' Collapses peaks with the median-of-significant rule, then joins the
#' directional peak counts.
#'
#' @param annotated data.frame from annotate_aging_peaks().
#' @param fdr_threshold FDR cutoff for a significant peak.
#' @param prefix Column-name prefix.
#' @return data.frame with one row per gene
build_genotype_gene_table <- function(annotated, fdr_threshold, prefix) {
  collapsed <- aggregate_diffbind_by_gene(annotated,
                                          method = "median_significant",
                                          fdr_threshold = fdr_threshold,
                                          prefix = prefix)
  directions <- count_gene_aging_peaks(annotated, fdr_threshold, prefix)
  out <- dplyr::left_join(collapsed, directions, by = "gene_name")

  cat(sprintf("  %s: %s genes carry a peak, %s have an aging-up peak, %s an aging-down peak\n",
              prefix, format(nrow(out), big.mark = ","),
              format(sum(out[[paste0(prefix, "_n_sig_up")]] > 0), big.mark = ","),
              format(sum(out[[paste0(prefix, "_n_sig_down")]] > 0), big.mark = ",")))
  out
}

#' Join the two genotype gene tables and label the aging gene classes.
#'
#' The analysis universe is the genes that carry at least one peak in both
#' genotypes, so that each gene has an aging test in each genotype.
#'
#' @param ctrl_gene Control gene table.
#' @param mut_gene Mutant gene table.
#' @return data.frame with one row per gene in the union of both tables
build_aging_gene_table <- function(ctrl_gene, mut_gene) {
  tbl <- dplyr::full_join(ctrl_gene, mut_gene, by = "gene_name")

  tbl$has_ctrl_peak <- !is.na(tbl$ctrl_aging_n_peaks)
  tbl$has_mut_peak <- !is.na(tbl$mut_aging_n_peaks)
  tbl$in_universe <- tbl$has_ctrl_peak & tbl$has_mut_peak

  tbl$ctrl_aging_up <- tbl$has_ctrl_peak & tbl$ctrl_aging_n_sig_up > 0
  tbl$mut_aging_up <- tbl$has_mut_peak & tbl$mut_aging_n_sig_up > 0
  tbl$ctrl_aging_down <- tbl$has_ctrl_peak & tbl$ctrl_aging_n_sig_down > 0
  tbl$mut_aging_down <- tbl$has_mut_peak & tbl$mut_aging_n_sig_down > 0
  tbl$ctrl_aging_sig <- tbl$has_ctrl_peak & tbl$ctrl_aging_n_sig > 0
  tbl$mut_aging_sig <- tbl$has_mut_peak & tbl$mut_aging_n_sig > 0

  tbl$mut_specific_aging_up <- tbl$in_universe & tbl$mut_aging_up &
    !tbl$ctrl_aging_sig
  tbl$shared_aging_up <- tbl$in_universe & tbl$mut_aging_up & tbl$ctrl_aging_up
  tbl$ctrl_specific_aging_up <- tbl$in_universe & tbl$ctrl_aging_up &
    !tbl$mut_aging_sig

  cat(sprintf("  Genes with a peak in both genotypes: %s of %s in the union\n",
              format(sum(tbl$in_universe), big.mark = ","),
              format(nrow(tbl), big.mark = ",")))
  cat(sprintf("  Control only: %s genes | mutant only: %s genes\n",
              format(sum(tbl$has_ctrl_peak & !tbl$has_mut_peak), big.mark = ","),
              format(sum(tbl$has_mut_peak & !tbl$has_ctrl_peak), big.mark = ",")))
  tbl
}

#' Counts behind the aging-up gene Venn.
#'
#' @param universe Rows of the aging gene table inside the analysis universe.
#' @return data.frame with one row per category
summarise_aging_up_overlap <- function(universe) {
  data.frame(
    category = c("Genes in universe",
                 "Control aging-up genes",
                 "Mutant aging-up genes",
                 "Aging-up in both",
                 "Control-only aging-up",
                 "Mutant-only aging-up",
                 "Mutant-specific aging-up (control not significant)"),
    n_genes = c(nrow(universe),
                sum(universe$ctrl_aging_up),
                sum(universe$mut_aging_up),
                sum(universe$ctrl_aging_up & universe$mut_aging_up),
                sum(universe$ctrl_aging_up & !universe$mut_aging_up),
                sum(universe$mut_aging_up & !universe$ctrl_aging_up),
                sum(universe$mut_specific_aging_up)),
    stringsAsFactors = FALSE
  )
}

#' Two-set Venn of the aging-up gene sets.
#'
#' @param ctrl_genes Control aging-up gene symbols.
#' @param mut_genes Mutant aging-up gene symbols.
#' @param subtitle Figure subtitle.
#' @return ggplot object
plot_aging_up_venn <- function(ctrl_genes, mut_genes, subtitle) {
  venn_sets <- list("Control aging-up" = ctrl_genes,
                    "Mutant aging-up" = mut_genes)

  ggVennDiagram(venn_sets, label = "count", label_alpha = 0,
                set_size = 4.2, label_size = 4.2) +
    scale_fill_gradient(low = "white", high = AGING_STATUS_COLORS[[AGING_UP]],
                        guide = "none") +
    scale_x_continuous(expand = expansion(mult = 0.12)) +
    labs(
      title = "Genes that gain MeCP2 with age",
      subtitle = subtitle
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
}

# =============================================================================
# SHARED PEAKS
# =============================================================================

#' Pair every control peak with the mutant peak it overlaps most.
#'
#' Overlapping pairs are ordered by shared base pairs, then reduced to one
#' mutant peak per control peak and one control peak per mutant peak, which
#' makes the pairing one to one.
#'
#' @param ctrl_db Control aging table.
#' @param mut_db Mutant aging table.
#' @return data.frame with one row per shared peak pair
pair_shared_peaks <- function(ctrl_db, mut_db) {
  ctrl_gr <- GRanges(seqnames = ctrl_db$Chr,
                     ranges = IRanges(start = ctrl_db$Start, end = ctrl_db$End))
  mut_gr <- GRanges(seqnames = mut_db$Chr,
                    ranges = IRanges(start = mut_db$Start, end = mut_db$End))

  hits <- findOverlaps(ctrl_gr, mut_gr)
  if (length(hits) == 0) {
    stop("No control aging peak overlaps any mutant aging peak. The two ",
         "DiffBind tables do not share a peak set.")
  }

  overlap_bp <- width(pintersect(ctrl_gr[queryHits(hits)],
                                 mut_gr[subjectHits(hits)]))
  pairs <- data.table::data.table(
    ctrl_idx = queryHits(hits),
    mut_idx = subjectHits(hits),
    overlap_bp = as.integer(overlap_bp)
  )
  data.table::setorder(pairs, -overlap_bp)
  pairs <- as.data.frame(pairs)
  pairs <- pairs[!duplicated(pairs$ctrl_idx), , drop = FALSE]
  pairs <- pairs[!duplicated(pairs$mut_idx), , drop = FALSE]

  out <- data.frame(
    chr = ctrl_db$Chr[pairs$ctrl_idx],
    ctrl_start = ctrl_db$Start[pairs$ctrl_idx],
    ctrl_end = ctrl_db$End[pairs$ctrl_idx],
    mut_start = mut_db$Start[pairs$mut_idx],
    mut_end = mut_db$End[pairs$mut_idx],
    overlap_bp = pairs$overlap_bp,
    ctrl_fold = ctrl_db$Fold[pairs$ctrl_idx],
    ctrl_fdr = ctrl_db$FDR[pairs$ctrl_idx],
    ctrl_status = ctrl_db$aging_status[pairs$ctrl_idx],
    mut_fold = mut_db$Fold[pairs$mut_idx],
    mut_fdr = mut_db$FDR[pairs$mut_idx],
    mut_status = mut_db$aging_status[pairs$mut_idx],
    stringsAsFactors = FALSE
  )
  out$fold_difference <- out$mut_fold - out$ctrl_fold
  out$up_in_both <- out$ctrl_status == AGING_UP & out$mut_status == AGING_UP

  cat(sprintf("  Shared peak pairs: %s (control %s peaks, mutant %s peaks)\n",
              format(nrow(out), big.mark = ","),
              format(nrow(ctrl_db), big.mark = ","),
              format(nrow(mut_db), big.mark = ",")))
  cat(sprintf("  Pairs aging-up in both genotypes: %s\n",
              format(sum(out$up_in_both), big.mark = ",")))
  out
}

#' Attach the control peak's gene symbol to each shared pair.
#'
#' Matches on the control peak coordinates, which ChIPseeker carries through
#' unchanged. Pairs whose peak was dropped by the annotation keep NA.
#'
#' @param pairs data.frame from pair_shared_peaks().
#' @param annotated Control annotation from annotate_aging_peaks().
#' @return pairs with an added gene_name column
attach_pair_symbols <- function(pairs, annotated) {
  annotated_key <- paste(as.character(annotated$seqnames),
                         as.integer(annotated$start),
                         as.integer(annotated$end), sep = ":")
  if (anyDuplicated(annotated_key) > 0) {
    stop("The control annotation holds duplicated peak coordinates, so gene ",
         "symbols cannot be matched to shared pairs unambiguously.")
  }

  pair_key <- paste(pairs$chr, as.integer(pairs$ctrl_start),
                    as.integer(pairs$ctrl_end), sep = ":")
  idx <- match(pair_key, annotated_key)
  pairs$gene_name <- annotated$SYMBOL[idx]

  cat(sprintf("  Shared pairs with a control gene symbol: %s of %s\n",
              format(sum(!is.na(pairs$gene_name)), big.mark = ","),
              format(nrow(pairs), big.mark = ",")))
  pairs
}

#' Correlation and paired Wilcoxon test on the shared peak folds.
#'
#' @param pairs data.frame from pair_shared_peaks().
#' @return one-row data.frame of the statistics
test_shared_peak_folds <- function(pairs) {
  spearman <- cor.test(pairs$ctrl_fold, pairs$mut_fold, method = "spearman",
                       exact = FALSE)
  pearson <- cor.test(pairs$ctrl_fold, pairs$mut_fold, method = "pearson")
  wilcox <- wilcox.test(pairs$mut_fold, pairs$ctrl_fold, paired = TRUE,
                        exact = FALSE)

  out <- data.frame(
    n_pairs = nrow(pairs),
    median_ctrl_fold = median(pairs$ctrl_fold),
    median_mut_fold = median(pairs$mut_fold),
    median_fold_difference = median(pairs$fold_difference),
    mean_fold_difference = mean(pairs$fold_difference),
    n_mut_above_identity = sum(pairs$fold_difference > 0),
    pct_mut_above_identity = 100 * sum(pairs$fold_difference > 0) / nrow(pairs),
    spearman_rho = unname(spearman$estimate),
    spearman_p = spearman$p.value,
    pearson_r = unname(pearson$estimate),
    pearson_p = pearson$p.value,
    wilcoxon_V = unname(wilcox$statistic),
    wilcoxon_p = wilcox$p.value,
    stringsAsFactors = FALSE
  )

  cat(sprintf("  Median aging fold: control %.4f, mutant %.4f, difference %.4f\n",
              out$median_ctrl_fold, out$median_mut_fold,
              out$median_fold_difference))
  cat(sprintf("  Mutant above the identity line: %s of %s pairs (%.1f%%)\n",
              format(out$n_mut_above_identity, big.mark = ","),
              format(out$n_pairs, big.mark = ","), out$pct_mut_above_identity))
  cat(sprintf("  Spearman rho = %.3f, %s\n", out$spearman_rho,
              fmt_p(out$spearman_p)))
  cat(sprintf("  Paired Wilcoxon (mutant vs control fold): %s\n",
              fmt_p(out$wilcoxon_p)))
  out
}

#' Scatter of the control aging fold against the mutant aging fold.
#'
#' Every shared pair enters the binned density layer. Pairs that are aging-up
#' in both genotypes are drawn again as points.
#'
#' @param pairs data.frame from attach_pair_symbols().
#' @param label_df Rows of pairs carrying a key-gene label.
#' @param stats One-row data.frame from test_shared_peak_folds().
#' @return ggplot object
plot_shared_peak_scatter <- function(pairs, label_df, stats) {
  axis_limit <- max(abs(c(pairs$ctrl_fold, pairs$mut_fold)))
  both_up <- pairs[pairs$up_in_both, , drop = FALSE]

  subtitle <- sprintf(
    paste("n = %s shared peaks | Spearman rho = %.3f, %s | paired Wilcoxon %s",
          "| median fold difference = %.3f | %.1f%% of pairs above the identity line"),
    format(stats$n_pairs, big.mark = ","), stats$spearman_rho,
    fmt_p(stats$spearman_p), fmt_p(stats$wilcoxon_p),
    stats$median_fold_difference, stats$pct_mut_above_identity)

  ggplot(pairs, aes(x = ctrl_fold, y = mut_fold)) +
    geom_bin2d(aes(fill = after_stat(log10(count))), bins = 140) +
    geom_point(data = both_up, aes(x = ctrl_fold, y = mut_fold),
               inherit.aes = FALSE, color = AGING_STATUS_COLORS[[AGING_UP]],
               alpha = 0.22, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black",
                linewidth = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, color = "#1B7837",
                linewidth = 0.8, se = TRUE, alpha = 0.25) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey45",
               linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey45",
               linewidth = 0.35) +
    geom_text_repel(data = label_df, aes(label = gene_name), size = 3,
                    fontface = "italic", color = "grey10", max.overlaps = 30,
                    segment.color = "grey55", segment.size = 0.3,
                    min.segment.length = 0) +
    scale_fill_gradientn(colors = c("#F7FBFF", "#9ECAE1", "#2171B5", "#08306B"),
                         name = "log10 shared\npeaks per bin") +
    coord_fixed(xlim = c(-axis_limit, axis_limit),
                ylim = c(-axis_limit, axis_limit)) +
    labs(
      title = "MeCP2 aging fold at shared peaks: control against mutant",
      subtitle = subtitle,
      x = "Control aging fold, log2(adult / young)",
      y = "Mutant aging fold, log2(adult / young)"
    ) +
    theme_emseq() +
    theme(legend.position = "right")
}

# =============================================================================
# GO ENRICHMENT
# =============================================================================

#' Map gene symbols to Entrez identifiers.
#'
#' @param symbols Character vector of gene symbols.
#' @param label Display name used in messages.
#' @return character vector of unique Entrez identifiers
map_symbols_to_entrez <- function(symbols, label) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  mapping <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)
  entrez <- unique(mapping$ENTREZID[!is.na(mapping$ENTREZID)])

  cat(sprintf("  %s: %s symbols mapped to %s Entrez identifiers\n",
              label, format(length(symbols), big.mark = ","),
              format(length(entrez), big.mark = ",")))
  if (length(entrez) == 0) {
    stop(label, ": no symbol mapped to an Entrez identifier, so GO enrichment ",
         "cannot run.")
  }
  entrez
}

#' GO biological process enrichment of the mutant-specific aging genes.
#'
#' @param gene_entrez Entrez identifiers of the mutant-specific genes.
#' @param universe_entrez Entrez identifiers of all peak-carrying genes.
#' @return enrichResult object
run_go_enrichment <- function(gene_entrez, universe_entrez) {
  ego <- enrichGO(
    gene = gene_entrez,
    universe = universe_entrez,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat("  WARNING: GO BP enrichment returned no term at q < 0.05 for the ",
        "mutant-specific aging genes against the peak-carrying universe.\n")
    return(NULL)
  }

  cat(sprintf("  Significant GO BP terms: %s\n",
              format(nrow(as.data.frame(ego)), big.mark = ",")))
  ego
}

#' Dot plot of the enriched GO biological process terms.
#'
#' @param ego enrichResult from run_go_enrichment().
#' @param n_terms Terms drawn.
#' @param n_genes Genes tested.
#' @param n_universe Genes in the universe.
#' @return ggplot object
plot_go_dotplot <- function(ego, n_terms, n_genes, n_universe) {
  enrichplot::dotplot(ego, showCategory = n_terms) +
    labs(
      title = "GO biological process: mutant-specific aging genes",
      subtitle = sprintf(
        "%s genes that gain MeCP2 with age only in the mutant | universe = %s genes with a peak in both genotypes | top %d terms by adjusted p",
        format(n_genes, big.mark = ","), format(n_universe, big.mark = ","),
        n_terms)
    ) +
    theme_emseq(base_size = 11) +
    theme(axis.text.y = element_text(size = 9))
}

# =============================================================================
# mCH INTEGRATION
# =============================================================================

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

#' Join the aging universe to the deduplicated adult mCH results.
#'
#' @param universe Rows of the aging gene table inside the analysis universe.
#' @param mch mch_results data.frame.
#' @return data.frame with one row per gene carrying both aging and mCH data
join_aging_to_mch <- function(universe, mch) {
  keep_idx <- deduplicate_mch_row_indices(mch)
  mch_gene <- mch[keep_idx, c("gene_name", "chr", "gene_length", "mch_ctrl",
                              "mch_mut", "mch_diff", "edger_logFC", "edger_fdr",
                              "mch_sig", "mch_hyper", "mch_hypo",
                              "mch_direction"), drop = FALSE]
  cat(sprintf("  Deduplicated %s mCH rows to %s gene names\n",
              format(nrow(mch), big.mark = ","),
              format(nrow(mch_gene), big.mark = ",")))

  joined <- dplyr::inner_join(universe, mch_gene, by = "gene_name")
  if (anyNA(joined$mch_diff)) {
    stop(sum(is.na(joined$mch_diff)), " genes carry a missing mch_diff after ",
         "the join with the aging universe, so the mCH comparison cannot run.")
  }

  joined$mch_diff_pct <- 100 * joined$mch_diff
  joined$aging_group <- factor(
    ifelse(joined$mut_specific_aging_up, GROUP_MUT_SPECIFIC, GROUP_OTHER),
    levels = AGING_GROUP_ORDER)

  cat(sprintf("  Genes with both a MeCP2 aging test and an mCH test: %s\n",
              format(nrow(joined), big.mark = ",")))
  cat(sprintf("  Of these, %s are mutant-specific aging-up and %s are mCH significant\n",
              format(sum(joined$mut_specific_aging_up), big.mark = ","),
              format(sum(joined$mch_sig), big.mark = ",")))
  joined
}

#' Register the three mCH overlap Fisher tests.
#'
#' @param joined data.frame from join_aging_to_mch().
#' @param out_dir Section output directory.
#' @return data.frame with one row per test
run_mch_fisher_tests <- function(joined, out_dir) {
  specs <- list(
    list(test_id = "mut_specific_aging_vs_mch_sig", col_var = "mch_sig",
         description = paste("Are genes that gain MeCP2 with age only in the",
                             "mutant more often differentially methylated for",
                             "mCH in the adult?")),
    list(test_id = "mut_specific_aging_vs_mch_hyper", col_var = "mch_hyper",
         description = paste("Are genes that gain MeCP2 with age only in the",
                             "mutant more often mCH hypermethylated in the",
                             "adult mutant?")),
    list(test_id = "mut_specific_aging_vs_mch_hypo", col_var = "mch_hypo",
         description = paste("Are genes that gain MeCP2 with age only in the",
                             "mutant more often mCH hypomethylated in the",
                             "adult mutant?"))
  )

  rows <- lapply(specs, function(spec) {
    ft <- register_fisher_test(
      section = SECTION_ID, test_id = spec$test_id,
      description = spec$description,
      gene_df = joined, row_var = "mut_specific_aging_up",
      col_var = spec$col_var, output_dir = out_dir)

    data.frame(
      test_id = spec$test_id,
      row_var = "mut_specific_aging_up",
      col_var = spec$col_var,
      n_genes = nrow(joined),
      n_row_true = sum(joined$mut_specific_aging_up),
      n_col_true = sum(joined[[spec$col_var]]),
      n_both_true = sum(joined$mut_specific_aging_up & joined[[spec$col_var]]),
      odds_ratio = unname(ft$estimate),
      conf_low = ft$conf.int[1],
      conf_high = ft$conf.int[2],
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Violin of mch_diff for mutant-specific aging genes against the rest.
#'
#' Each group carries its gene count and median on the figure.
#'
#' @param joined data.frame from join_aging_to_mch().
#' @param stats data.frame from summarise_groups() on the same columns, with a
#'   label column from group_label().
#' @param wilcox htest from wilcox.test().
#' @return ggplot object
plot_mch_by_aging_group <- function(joined, stats, wilcox) {
  y_min <- min(joined$mch_diff_pct)
  y_max <- max(joined$mch_diff_pct)
  y_span <- y_max - y_min

  stats$label_y <- y_max + 0.10 * y_span
  y_stat <- y_max + 0.26 * y_span

  ggplot(joined, aes(x = aging_group, y = mch_diff_pct, fill = aging_group)) +
    geom_violin(alpha = 0.6, scale = "width",
                draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 0.35) +
    geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white", alpha = 0.85,
                 linewidth = 0.35) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black",
               linewidth = 0.4) +
    geom_text(data = stats, aes(x = aging_group, y = label_y, label = label),
              inherit.aes = FALSE, size = 3.3, lineheight = 1.1) +
    annotate("text", x = 1.5, y = y_stat,
             label = sprintf("Wilcoxon rank sum: %s", fmt_p(wilcox$p.value)),
             size = 3.3, fontface = "italic") +
    scale_fill_manual(values = AGING_GROUP_COLORS) +
    scale_y_continuous(limits = c(y_min - 0.05 * y_span, y_stat + 0.12 * y_span)) +
    labs(
      title = "Adult mCH change at mutant-specific MeCP2 aging genes",
      subtitle = "Genes with a MeCP2 peak in both genotypes and an mCH test",
      x = "MeCP2 aging group",
      y = "mCH change (mutant - control, percentage points)"
    ) +
    theme_emseq() +
    theme(legend.position = "none")
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opt <- parse_options()
  out_dir <- opt$output_dir
  fdr_threshold <- opt$fdr_threshold
  go_top_terms <- opt$go_top_terms

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("================================================================================\n")
  cat("SECTION 60_04: MeCP2 DEVELOPMENTAL TRAJECTORY (YOUNG vs ADULT)\n")
  cat("================================================================================\n")
  cat("Output dir:    ", out_dir, "\n", sep = "")
  cat("FDR threshold: ", fdr_threshold, "\n", sep = "")
  cat("GO terms drawn:", go_top_terms, "\n")
  cat("\n")

  stopifnot(
    "Control MeCP2 aging DiffBind file not found" =
      file.exists(MECP2_PATHS$ctrl_aging),
    "Mutant MeCP2 aging DiffBind file not found" =
      file.exists(MECP2_PATHS$mut_aging)
  )

  # --- Load the two aging tables -------------------------------------------
  cat("--- Loading MeCP2 aging DiffBind tables ---\n")
  ctrl_db <- load_aging_diffbind(MECP2_PATHS$ctrl_aging, "Control aging",
                                 fdr_threshold)
  mut_db <- load_aging_diffbind(MECP2_PATHS$mut_aging, "Mutant aging",
                                fdr_threshold)

  # --- Analysis 1: peak counts by aging status ------------------------------
  cat("\n--- Analysis 1: aging peak counts by genotype ---\n")
  count_df <- count_aging_peaks(ctrl_db, mut_db)
  for (i in seq_len(nrow(count_df))) {
    cat(sprintf("  %-8s %-16s %s peaks (%.2f%% of the genotype's peaks)\n",
                count_df$genotype[i], count_df$aging_status[i],
                format(count_df$n_peaks[i], big.mark = ","),
                count_df$pct_of_peaks[i]))
  }

  peak_fisher <- rbind(
    aging_proportion_fisher(ctrl_db, mut_db, AGING_UP),
    aging_proportion_fisher(ctrl_db, mut_db, AGING_DOWN)
  )
  for (i in seq_len(nrow(peak_fisher))) {
    cat(sprintf("  Peak-level Fisher, %s in mutant vs control: OR = %.3f, %s\n",
                peak_fisher$aging_status[i], peak_fisher$odds_ratio[i],
                fmt_p(peak_fisher$p_value[i])))
  }

  up_row <- peak_fisher[peak_fisher$aging_status == AGING_UP, ]
  count_subtitle <- sprintf(
    "Aging up in mutant %.2f%% vs control %.2f%% of peaks | peak-level Fisher OR = %.2f, %s",
    up_row$pct_mut_event, up_row$pct_ctrl_event, up_row$odds_ratio,
    fmt_p(up_row$p_value))

  p_counts <- plot_aging_peak_counts(count_df, count_subtitle)
  save_multiformat_ggplot(p_counts,
                          file.path(out_dir, "60_04a_aging_peak_counts"),
                          width = 10, height = 8)

  # --- Annotate both peak sets to genes -------------------------------------
  cat("\n--- Annotating aging peaks to genes ---\n")
  ctrl_annotated <- annotate_aging_peaks(ctrl_db, "Control aging")
  mut_annotated <- annotate_aging_peaks(mut_db, "Mutant aging")

  cat("\n--- Collapsing aging peaks to genes ---\n")
  ctrl_gene <- build_genotype_gene_table(ctrl_annotated, fdr_threshold,
                                         "ctrl_aging")
  mut_gene <- build_genotype_gene_table(mut_annotated, fdr_threshold,
                                        "mut_aging")
  gene_table <- build_aging_gene_table(ctrl_gene, mut_gene)
  universe <- gene_table[gene_table$in_universe, , drop = FALSE]

  # --- Analysis 2: aging-up gene overlap ------------------------------------
  cat("\n--- Analysis 2: aging-up gene overlap between genotypes ---\n")
  overlap_df <- summarise_aging_up_overlap(universe)
  for (i in seq_len(nrow(overlap_df))) {
    cat(sprintf("  %-52s %s\n", overlap_df$category[i],
                format(overlap_df$n_genes[i], big.mark = ",")))
  }

  ctrl_up_genes <- universe$gene_name[universe$ctrl_aging_up]
  mut_up_genes <- universe$gene_name[universe$mut_aging_up]

  venn_subtitle <- sprintf(
    "Genes with a MeCP2 peak in both genotypes (n = %s) | control %s, mutant %s, shared %s",
    format(nrow(universe), big.mark = ","),
    format(length(ctrl_up_genes), big.mark = ","),
    format(length(mut_up_genes), big.mark = ","),
    format(length(intersect(ctrl_up_genes, mut_up_genes)), big.mark = ","))

  p_venn <- plot_aging_up_venn(ctrl_up_genes, mut_up_genes, venn_subtitle)
  save_multiformat_ggplot(p_venn,
                          file.path(out_dir, "60_04b_aging_up_gene_venn"),
                          width = 9, height = 7)

  gene_set_df <- data.frame(
    gene_name = universe$gene_name,
    ctrl_aging_up = universe$ctrl_aging_up,
    mut_aging_up = universe$mut_aging_up,
    ctrl_aging_sig = universe$ctrl_aging_sig,
    mut_aging_sig = universe$mut_aging_sig,
    mut_specific_aging_up = universe$mut_specific_aging_up,
    shared_aging_up = universe$shared_aging_up,
    ctrl_specific_aging_up = universe$ctrl_specific_aging_up,
    stringsAsFactors = FALSE
  )

  # --- Analysis 3: shared peaks ---------------------------------------------
  cat("\n--- Analysis 3: aging fold at peaks shared by both genotypes ---\n")
  pairs <- pair_shared_peaks(ctrl_db, mut_db)
  pairs <- attach_pair_symbols(pairs, ctrl_annotated)
  shared_stats <- test_shared_peak_folds(pairs)

  key_pairs <- pairs[!is.na(pairs$gene_name) &
                       pairs$gene_name %in% KEY_GENES, , drop = FALSE]
  key_pairs <- key_pairs[order(-abs(key_pairs$fold_difference)), , drop = FALSE]
  label_df <- key_pairs[!duplicated(key_pairs$gene_name), , drop = FALSE]
  cat(sprintf("  Labelling %d of %d key genes, at the shared peak with the largest fold difference\n",
              nrow(label_df), length(KEY_GENES)))

  p_scatter <- plot_shared_peak_scatter(pairs, label_df, shared_stats)
  save_multiformat_ggplot(p_scatter,
                          file.path(out_dir, "60_04c_shared_peak_aging_fold"),
                          width = 11, height = 10)

  # --- Analysis 4: mutant-specific aging genes ------------------------------
  cat("\n--- Analysis 4: mutant-specific aging-up genes ---\n")
  mut_specific <- universe[universe$mut_specific_aging_up, , drop = FALSE]
  cat(sprintf("  Mutant-specific aging-up genes: %s of %s in the universe (%.1f%%)\n",
              format(nrow(mut_specific), big.mark = ","),
              format(nrow(universe), big.mark = ","),
              100 * nrow(mut_specific) / nrow(universe)))
  cat(sprintf("  Median mutant aging fold at these genes: %.4f\n",
              median(mut_specific$mut_aging_fold)))

  # --- Analysis 5: GO enrichment --------------------------------------------
  cat("\n--- Analysis 5: GO biological process enrichment ---\n")
  gene_entrez <- map_symbols_to_entrez(mut_specific$gene_name,
                                       "Mutant-specific aging genes")
  universe_entrez <- map_symbols_to_entrez(universe$gene_name,
                                           "Peak-carrying universe")
  ego <- run_go_enrichment(gene_entrez, universe_entrez)

  if (is.null(ego)) {
    go_df <- data.frame(ID = character(0), Description = character(0),
                        GeneRatio = character(0), pvalue = numeric(0),
                        qvalue = numeric(0), stringsAsFactors = FALSE)
    cat("  Skipping GO dot plot (no significant terms).\n")
  } else {
    go_df <- as.data.frame(ego)
    n_terms_drawn <- min(go_top_terms, nrow(go_df))
    p_go <- plot_go_dotplot(ego, n_terms_drawn, nrow(mut_specific), nrow(universe))
    save_multiformat_ggplot(p_go,
                            file.path(out_dir, "60_04d_mut_specific_go_bp"),
                            width = 11, height = 11)

    for (i in seq_len(min(10, nrow(go_df)))) {
      cat(sprintf("  %-58s %s genes, q = %.3g\n",
                  substr(go_df$Description[i], 1, 58), go_df$GeneRatio[i],
                  go_df$qvalue[i]))
    }
  }

  # --- Analysis 6: overlap with the adult mCH result ------------------------
  cat("\n--- Analysis 6: overlap with the adult mCH result ---\n")
  joined <- join_aging_to_mch(universe, mch_results)

  mch_fisher <- run_mch_fisher_tests(joined, out_dir)

  mch_group_stats <- summarise_groups(joined, "aging_group", "mch_diff_pct")
  for (i in seq_len(nrow(mch_group_stats))) {
    cat(sprintf("  %-26s n = %s, median mCH change = %.4f pp, mean = %.4f pp\n",
                as.character(mch_group_stats$aging_group[i]),
                format(mch_group_stats$n[i], big.mark = ","),
                mch_group_stats$median[i], mch_group_stats$mean[i]))
  }

  mch_wilcox <- wilcox.test(
    joined$mch_diff_pct[joined$aging_group == GROUP_MUT_SPECIFIC],
    joined$mch_diff_pct[joined$aging_group == GROUP_OTHER],
    exact = FALSE)
  cat(sprintf("  Wilcoxon mCH change, mutant-specific vs other: W = %.0f, %s\n",
              unname(mch_wilcox$statistic), fmt_p(mch_wilcox$p.value)))

  mch_wilcox_tbl <- data.frame(
    group_1 = GROUP_MUT_SPECIFIC,
    group_2 = GROUP_OTHER,
    n_1 = sum(joined$aging_group == GROUP_MUT_SPECIFIC),
    n_2 = sum(joined$aging_group == GROUP_OTHER),
    median_1 = median(joined$mch_diff_pct[joined$aging_group == GROUP_MUT_SPECIFIC]),
    median_2 = median(joined$mch_diff_pct[joined$aging_group == GROUP_OTHER]),
    W = unname(mch_wilcox$statistic),
    p_value = mch_wilcox$p.value,
    stringsAsFactors = FALSE
  )
  mch_wilcox_tbl$median_difference <- mch_wilcox_tbl$median_1 -
    mch_wilcox_tbl$median_2

  # Figure text stays out of mch_group_stats, which is written as a table below.
  mch_group_labels <- mch_group_stats
  mch_group_labels$label <- group_label(mch_group_stats, digits = 4)

  p_mch <- plot_mch_by_aging_group(joined, mch_group_labels, mch_wilcox)
  save_multiformat_ggplot(p_mch,
                          file.path(out_dir, "60_04e_mch_diff_by_aging_group"),
                          width = 9, height = 8)

  # --- Tables ---------------------------------------------------------------
  cat("\n--- Writing tables ---\n")
  write_section_table(count_df,
                      file.path(out_dir, "60_04_aging_peak_counts.tsv"))
  write_section_table(peak_fisher,
                      file.path(out_dir, "60_04_aging_peak_proportion_fisher.tsv"))
  write_section_table(gene_table,
                      file.path(out_dir, "60_04_aging_gene_table.tsv"))
  write_section_table(overlap_df,
                      file.path(out_dir, "60_04_aging_up_gene_overlap.tsv"))
  write_section_table(gene_set_df,
                      file.path(out_dir, "60_04_aging_up_gene_sets.tsv"))
  write_section_table(pairs,
                      file.path(out_dir, "60_04_shared_peak_pairs.tsv"))
  write_section_table(shared_stats,
                      file.path(out_dir, "60_04_shared_peak_fold_stats.tsv"))
  write_section_table(mut_specific,
                      file.path(out_dir, "60_04_mut_specific_aging_genes.tsv"))
  write_section_table(go_df,
                      file.path(out_dir, "60_04_mut_specific_go_bp.tsv"))
  write_section_table(joined,
                      file.path(out_dir, "60_04_aging_mch_gene_level.tsv"))
  write_section_table(mch_fisher,
                      file.path(out_dir, "60_04_mch_overlap_fisher_summary.tsv"))
  write_section_table(mch_group_stats,
                      file.path(out_dir,
                                "60_04_mch_diff_by_aging_group_summary.tsv"))
  write_section_table(mch_wilcox_tbl,
                      file.path(out_dir,
                                "60_04_mch_diff_by_aging_group_wilcoxon.tsv"))

  # --- Summary --------------------------------------------------------------
  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 60_04 SUMMARY\n")
  cat("================================================================================\n")
  cat(sprintf("Control aging peaks:            %s (%s up, %s down)\n",
              format(nrow(ctrl_db), big.mark = ","),
              format(sum(ctrl_db$aging_status == AGING_UP), big.mark = ","),
              format(sum(ctrl_db$aging_status == AGING_DOWN), big.mark = ",")))
  cat(sprintf("Mutant aging peaks:             %s (%s up, %s down)\n",
              format(nrow(mut_db), big.mark = ","),
              format(sum(mut_db$aging_status == AGING_UP), big.mark = ","),
              format(sum(mut_db$aging_status == AGING_DOWN), big.mark = ",")))
  cat(sprintf("Peak-level Fisher, aging up:    OR = %.3f, %s\n",
              up_row$odds_ratio, fmt_p(up_row$p_value)))
  cat(sprintf("Genes with a peak in both:      %s\n",
              format(nrow(universe), big.mark = ",")))
  cat(sprintf("Aging-up genes:                 control %s, mutant %s, shared %s\n",
              format(length(ctrl_up_genes), big.mark = ","),
              format(length(mut_up_genes), big.mark = ","),
              format(length(intersect(ctrl_up_genes, mut_up_genes)), big.mark = ",")))
  cat(sprintf("Mutant-specific aging-up genes: %s\n",
              format(nrow(mut_specific), big.mark = ",")))
  cat(sprintf("Shared peak pairs:              %s, Spearman rho = %.3f, paired Wilcoxon %s\n",
              format(shared_stats$n_pairs, big.mark = ","),
              shared_stats$spearman_rho, fmt_p(shared_stats$wilcoxon_p)))
  cat(sprintf("Median aging fold difference:   %.4f (mutant minus control)\n",
              shared_stats$median_fold_difference))
  cat(sprintf("GO BP terms at q < 0.05:        %s\n",
              format(nrow(go_df), big.mark = ",")))
  for (i in seq_len(nrow(mch_fisher))) {
    cat(sprintf("Fisher %-32s OR = %.3f, %s\n",
                mch_fisher$test_id[i], mch_fisher$odds_ratio[i],
                fmt_p(mch_fisher$p_value[i])))
  }
  cat(sprintf("Wilcoxon mCH change by group:   %s\n", fmt_p(mch_wilcox$p.value)))
  cat("\nSection 60_04 complete.\n\n")
}

main()
