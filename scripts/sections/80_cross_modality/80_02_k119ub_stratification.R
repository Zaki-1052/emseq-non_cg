# scripts/sections/80_cross_modality/80_02_k119ub_stratification.R
#
# Section 80_02: K119ub intragenic vs intergenic stratification and quadrant
# reorientation.
#
# What this tests
#   Step 07's quadrant analysis dropped 28% of K119ub DiffBind peaks because
#   they fell outside gene bodies. This section asks: does the intergenic
#   fraction behave differently? And does the relationship between K119ub and
#   MeCP2 look different when we separate peaks by genomic location?
#
# Analyses
#   1. Axis-swapped quadrant scatter (H2AK119Ub x, MeCP2 y) to match the
#      causal model direction (H2AK119Ub upstream of MeCP2).
#   2. K119ub peak classification: intragenic vs intergenic using full peak
#      coordinates (not 400bp summit windows). Direction distribution for each.
#   3. MeCP2 peak classification: same binary split.
#   4. mCH at genes with intragenic-only, intergenic-only, or both K119ub peaks.
#   5. K119ub-up / mCH-down loci: genes where ubiquitination increases but
#      non-CG methylation decreases. Chromatin and MeCP2 characterisation.
#   6. Consensus peak overlap: intragenic vs intergenic, enhancer overlap.
#
# Reads
#   mch_results, k119ub_diffbind, mecp2_diffbind, k119ub_consensus,
#   mecp2_consensus, gene_bodies  (shared config)
#   gene_level_all_marks.tsv from 20_02  (HANDOFF_PATHS$gene_level_all_marks)
#   gene_chromatin_state.tsv from 10_01  (HANDOFF_PATHS$chromatin_state)
#   data/k119ub_gene_signal.tsv
#   data/h2aub_diffbind.txt  (raw, for full peak coordinates)
#   data/chromatin/activeenhancer.bed
#
# Writes
#   Figures and tables into OUTPUT_PATHS$cross_modality.
#   Fisher test shards into HANDOFF_PATHS$fisher_registry.

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)

SECTION_ID <- "80_02"

INTRAGENIC_PATTERN <- "^(Promoter|5' UTR|3' UTR|Exon|Intron)"

# =============================================================================
# OPTIONS
# =============================================================================

parse_options <- function() {
  option_list <- list(
    make_option("--output-dir", type = "character", default = OUTPUT_PATHS$cross_modality,
                dest = "output_dir",
                help = "Directory for figures and tables [default: %default]"),
    make_option("--fdr-threshold", type = "double", default = Q_THRESHOLD,
                dest = "fdr_threshold",
                help = "FDR cutoff for significance calls [default: %default]")
  )
  parse_args(OptionParser(option_list = option_list))
}

# =============================================================================
# DATA LOADING
# =============================================================================

fmt_int <- function(n) format(n, big.mark = ",", trim = TRUE)

#' Load the raw K119ub DiffBind file with full peak coordinates.
#'
#' The shared config pre-loads k119ub_diffbind using Summit_Chr (400bp windows).
#' This function reads the raw file and uses seqnames/start/end (full peaks).
load_k119ub_raw <- function(filepath) {
  stopifnot(file.exists(filepath))
  df <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "", fill = TRUE)

  required <- c("seqnames", "start", "end", "Fold", "FDR")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("Raw K119ub DiffBind is missing columns: ", paste(missing, collapse = ", "))
  }

  raw_chr <- as.character(df$seqnames)
  df$Chr <- ifelse(grepl("^chr", raw_chr), raw_chr, paste0("chr", raw_chr))
  df$peak_start <- df$start
  df$peak_end <- df$end
  df$peak_width <- df$end - df$start

  na_rows <- is.na(df$Chr) | is.na(df$peak_start) | is.na(df$peak_end)
  if (any(na_rows)) {
    cat(sprintf("  Dropping %d peaks with NA coordinates\n", sum(na_rows)))
    df <- df[!na_rows, , drop = FALSE]
  }

  df$direction <- "Unchanged"
  df$direction[df$FDR < Q_THRESHOLD & df$Fold > 0] <- "Gained"
  df$direction[df$FDR < Q_THRESHOLD & df$Fold < 0] <- "Lost"

  cat(sprintf("  K119ub raw: %s peaks (%s gained, %s lost) using full coordinates\n",
              fmt_int(nrow(df)), fmt_int(sum(df$direction == "Gained")),
              fmt_int(sum(df$direction == "Lost"))))
  df
}

#' Load the K119ub gene-body signal table.
load_k119ub_gene_signal <- function(filepath) {
  stopifnot(file.exists(filepath))
  df <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")
  required <- c("symbol", "gb_ctrl_signal", "gb_mut_signal", "gb_log2fc",
                 "gb_signal_class")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("K119ub gene signal missing columns: ", paste(missing, collapse = ", "))
  }

  df <- df[!is.na(df$symbol) & nzchar(df$symbol), , drop = FALSE]
  df$gene_span <- df$end - df$start
  ord <- order(df$symbol, -df$gene_span, df$entrez_id)
  keep <- sort(ord[!duplicated(df$symbol[ord])])
  df <- df[keep, , drop = FALSE]

  out <- data.frame(
    gene_name = df$symbol,
    k119ub_gb_log2fc = as.numeric(df$gb_log2fc),
    k119ub_gb_signal_class = as.character(df$gb_signal_class),
    stringsAsFactors = FALSE
  )

  cat(sprintf("  K119ub gene signal: %s genes (%s quantifiable)\n",
              fmt_int(nrow(out)),
              fmt_int(sum(out$k119ub_gb_signal_class == "quantifiable"))))
  out
}

#' Collapse MeCP2 DiffBind to gene-level using ChIPseeker + nearest-TSS rule.
aggregate_mecp2_to_genes <- function(diffbind, fdr_threshold) {
  annotated <- annotate_peaks_to_genes(diffbind, "MeCP2")
  usable <- !is.na(annotated$SYMBOL) & nzchar(annotated$SYMBOL) &
    !is.na(annotated$Fold) & !is.na(annotated$FDR)
  annotated <- annotated[usable, , drop = FALSE]
  out <- aggregate_diffbind_by_gene(annotated, method = "nearest_tss",
                                    fdr_threshold = fdr_threshold,
                                    prefix = "mecp2")
  cat(sprintf("  MeCP2 gene-level: %s genes\n", fmt_int(nrow(out))))
  out
}

# =============================================================================
# PEAK CLASSIFICATION
# =============================================================================

#' Annotate peaks and classify as intragenic or intergenic.
#'
#' Uses ChIPseeker with full peak coordinates (not 400bp summits).
#' Intragenic: promoter, UTR, exon, intron.
#' Intergenic: everything else (Distal Intergenic, Downstream).
classify_peak_location <- function(peaks_df, mark_name) {
  gr <- GRanges(seqnames = peaks_df$Chr,
                ranges = IRanges(start = peaks_df$peak_start,
                                 end = peaks_df$peak_end))
  mcols(gr)$Fold <- peaks_df$Fold
  mcols(gr)$FDR <- peaks_df$FDR
  mcols(gr)$direction <- peaks_df$direction

  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  cat(sprintf("  Annotating %s peaks with ChIPseeker...\n", mark_name))
  anno <- ChIPseeker::annotatePeak(gr, TxDb = txdb, tssRegion = c(-3000, 3000),
                                   annoDb = "org.Mm.eg.db", verbose = FALSE)
  anno_df <- as.data.frame(anno)

  anno_df$location <- ifelse(grepl(INTRAGENIC_PATTERN, anno_df$annotation),
                             "Intragenic", "Intergenic")
  anno_df$gene_symbol <- anno_df$SYMBOL

  n_intra <- sum(anno_df$location == "Intragenic")
  n_inter <- sum(anno_df$location == "Intergenic")
  cat(sprintf("  %s: %s intragenic (%.1f%%), %s intergenic (%.1f%%)\n",
              mark_name, fmt_int(n_intra), 100 * n_intra / nrow(anno_df),
              fmt_int(n_inter), 100 * n_inter / nrow(anno_df)))
  anno_df
}

# =============================================================================
# PLOTS
# =============================================================================

#' Quadrant scatter with H2AK119Ub on x-axis and MeCP2 on y-axis.
quadrant_scatter <- function(df, x_col, y_col, x_lab, y_lab, sig_col,
                             rho, n, key_genes) {
  df$sig_label <- ifelse(df[[sig_col]], "Significant", "Not Significant")
  df$sig_label <- factor(df$sig_label, levels = names(COLORS$significant))

  subtitle <- sprintf("Spearman rho = %.3f, n = %s", rho, fmt_int(n))

  key_df <- df[df$gene_name %in% key_genes, , drop = FALSE]

  q_counts <- table(
    x_pos = df[[x_col]] > 0,
    y_pos = df[[y_col]] > 0
  )

  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(aes(colour = sig_label), size = 0.5, alpha = 0.3) +
    scale_colour_manual(values = COLORS$significant, name = NULL) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    labs(x = x_lab, y = y_lab, subtitle = subtitle) +
    theme_emseq() +
    theme(legend.position = "bottom")

  # Quadrant count annotations
  xlim <- layer_scales(p)$x$range$range
  ylim <- layer_scales(p)$y$range$range
  if (is.null(xlim)) xlim <- range(df[[x_col]], na.rm = TRUE)
  if (is.null(ylim)) ylim <- range(df[[y_col]], na.rm = TRUE)

  q_labels <- data.frame(
    x = c(xlim[2] * 0.9, xlim[1] * 0.9, xlim[1] * 0.9, xlim[2] * 0.9),
    y = c(ylim[2] * 0.9, ylim[2] * 0.9, ylim[1] * 0.9, ylim[1] * 0.9),
    label = sprintf("Q%d\nn=%s", 1:4,
                    fmt_int(c(
                      sum(df[[x_col]] > 0 & df[[y_col]] > 0),
                      sum(df[[x_col]] < 0 & df[[y_col]] > 0),
                      sum(df[[x_col]] < 0 & df[[y_col]] < 0),
                      sum(df[[x_col]] > 0 & df[[y_col]] < 0)
                    ))),
    stringsAsFactors = FALSE
  )
  p <- p + geom_text(data = q_labels, aes(x = x, y = y, label = label),
                     inherit.aes = FALSE, size = 3, colour = "grey40", fontface = "bold")

  if (nrow(key_df) > 0) {
    p <- p +
      geom_point(data = key_df, colour = "black", size = 1.5, shape = 21,
                 fill = COLORS$direction[["Hypermethylated"]]) +
      geom_text_repel(data = key_df, aes(label = gene_name),
                      size = 2.8, max.overlaps = 20, seed = 42,
                      min.segment.length = 0)
  }
  p
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opts <- parse_options()
  out_dir <- opts$output_dir
  fdr <- opts$fdr_threshold
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("=== Section 80_02: K119ub stratification and quadrant reorientation ===\n\n")

  # =========================================================================
  # STEP 1: Axis-swapped quadrant scatter
  # =========================================================================

  cat("STEP 1: Axis-swapped quadrant scatter\n")

  k119ub_signal_path <- file.path(CODE_DIR, "data", "k119ub_gene_signal.tsv")
  k119ub_sig <- load_k119ub_gene_signal(k119ub_signal_path)

  mecp2_gene <- aggregate_mecp2_to_genes(mecp2_diffbind, fdr)

  quad_df <- merge(k119ub_sig, mecp2_gene, by = "gene_name", all = FALSE)
  quad_df <- quad_df[quad_df$k119ub_gb_signal_class == "quantifiable", , drop = FALSE]
  quad_df <- quad_df[!is.na(quad_df$mecp2_fold), , drop = FALSE]

  mch_sub <- data.frame(
    gene_name = mch_results$gene_name,
    mch_logfc = mch_results$edger_logFC,
    mch_sig = mch_results$mch_sig,
    mch_direction = mch_results$mch_direction,
    stringsAsFactors = FALSE
  )
  quad_df <- merge(quad_df, mch_sub, by = "gene_name", all.x = TRUE)
  quad_df$both_sig <- !is.na(quad_df$mecp2_has_sig) & quad_df$mecp2_has_sig &
    quad_df$k119ub_gb_log2fc != 0

  cat(sprintf("  Genes with both K119ub gene-body signal and MeCP2 fold: %s\n",
              fmt_int(nrow(quad_df))))

  rho_quad <- cor(quad_df$k119ub_gb_log2fc, quad_df$mecp2_fold,
                  method = "spearman", use = "complete.obs")

  p_quad <- quadrant_scatter(
    quad_df, "k119ub_gb_log2fc", "mecp2_fold",
    x_lab = expression("H2AK119Ub gene-body " * log[2] * "FC"),
    y_lab = expression("MeCP2 " * log[2] * "FC (nearest-TSS)"),
    sig_col = "mch_sig",
    rho = rho_quad, n = nrow(quad_df),
    key_genes = KEY_GENES
  ) + ggtitle("H2AK119Ub vs MeCP2 (causal axis orientation)")

  save_multiformat_ggplot(p_quad, file.path(out_dir, "80_02a_quadrant_reoriented"),
                          width = 9, height = 8)

  quad_summary <- data.frame(
    quadrant = paste0("Q", 1:4),
    description = c("K119ub_up + MeCP2_up", "K119ub_down + MeCP2_up",
                     "K119ub_down + MeCP2_down", "K119ub_up + MeCP2_down"),
    n_genes = c(
      sum(quad_df$k119ub_gb_log2fc > 0 & quad_df$mecp2_fold > 0),
      sum(quad_df$k119ub_gb_log2fc < 0 & quad_df$mecp2_fold > 0),
      sum(quad_df$k119ub_gb_log2fc < 0 & quad_df$mecp2_fold < 0),
      sum(quad_df$k119ub_gb_log2fc > 0 & quad_df$mecp2_fold < 0)
    ),
    spearman_rho = rho_quad,
    stringsAsFactors = FALSE
  )
  write_section_table(quad_summary, file.path(out_dir, "80_02_quadrant_summary.tsv"))

  # =========================================================================
  # STEP 2: K119ub peak location classification
  # =========================================================================

  cat("\nSTEP 2: K119ub peak location classification (full coordinates)\n")

  k119ub_raw_path <- file.path(CODE_DIR, "data", "h2aub_diffbind.txt")
  k119ub_raw <- load_k119ub_raw(k119ub_raw_path)
  k119ub_anno <- classify_peak_location(k119ub_raw, "K119ub")

  # Summary table
  loc_summary <- k119ub_anno %>%
    dplyr::group_by(location) %>%
    dplyr::summarise(
      n_peaks = dplyr::n(),
      median_width = median(width, na.rm = TRUE),
      n_gained = sum(direction == "Gained"),
      n_lost = sum(direction == "Lost"),
      n_unchanged = sum(direction == "Unchanged"),
      pct_gained = 100 * n_gained / dplyr::n(),
      pct_lost = 100 * n_lost / dplyr::n(),
      median_fold = median(Fold, na.rm = TRUE),
      .groups = "drop"
    )
  write_section_table(as.data.frame(loc_summary),
                      file.path(out_dir, "80_02_k119ub_location_summary.tsv"))

  cat(sprintf("  Summary written: %s rows\n", nrow(loc_summary)))

  # =========================================================================
  # STEP 3: K119ub direction x location
  # =========================================================================

  cat("\nSTEP 3: K119ub direction x location analysis\n")

  sig_peaks <- k119ub_anno[k119ub_anno$direction %in% c("Gained", "Lost"), , drop = FALSE]

  if (nrow(sig_peaks) > 0) {
    ct_dir_loc <- table(
      direction = sig_peaks$direction,
      location = sig_peaks$location
    )
    fisher_dir_loc <- fisher.test(ct_dir_loc)
    cat(sprintf("  Direction x location Fisher: OR = %.2f, p = %.2e\n",
                fisher_dir_loc$estimate, fisher_dir_loc$p.value))

    dir_loc_table <- as.data.frame.matrix(ct_dir_loc)
    dir_loc_table$direction <- rownames(dir_loc_table)
    write_section_table(dir_loc_table,
                        file.path(out_dir, "80_02_k119ub_direction_by_location.tsv"))
  }

  # Assign intergenic peaks to nearest gene, compare mCH
  cat("  Assigning intergenic peaks to nearest gene...\n")

  intra_genes <- unique(k119ub_anno$gene_symbol[k119ub_anno$location == "Intragenic" &
                                                   !is.na(k119ub_anno$gene_symbol)])
  inter_genes <- unique(k119ub_anno$gene_symbol[k119ub_anno$location == "Intergenic" &
                                                   !is.na(k119ub_anno$gene_symbol)])
  both_genes <- intersect(intra_genes, inter_genes)
  intra_only <- setdiff(intra_genes, inter_genes)
  inter_only <- setdiff(inter_genes, intra_genes)

  cat(sprintf("  Gene counts: %s intragenic-only, %s intergenic-only, %s both\n",
              fmt_int(length(intra_only)), fmt_int(length(inter_only)),
              fmt_int(length(both_genes))))

  mch_for_strat <- data.frame(
    gene_name = mch_results$gene_name,
    mch_logfc = mch_results$edger_logFC,
    mch_sig = mch_results$mch_sig,
    stringsAsFactors = FALSE
  )

  strat_list <- list(
    "Intragenic only" = mch_for_strat[mch_for_strat$gene_name %in% intra_only, ],
    "Intergenic only" = mch_for_strat[mch_for_strat$gene_name %in% inter_only, ],
    "Both" = mch_for_strat[mch_for_strat$gene_name %in% both_genes, ],
    "No K119ub peak" = mch_for_strat[!mch_for_strat$gene_name %in%
                                       c(intra_genes, inter_genes), ]
  )

  strat_stats <- do.call(rbind, lapply(names(strat_list), function(grp) {
    sub <- strat_list[[grp]]
    if (nrow(sub) == 0) return(NULL)
    data.frame(
      group = grp,
      n_genes = nrow(sub),
      n_sig = sum(sub$mch_sig, na.rm = TRUE),
      pct_sig = 100 * sum(sub$mch_sig, na.rm = TRUE) / nrow(sub),
      median_logfc = median(sub$mch_logfc, na.rm = TRUE),
      mean_logfc = mean(sub$mch_logfc, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  write_section_table(strat_stats,
                      file.path(out_dir, "80_02_mch_by_k119ub_location.tsv"))

  # Violin plot
  strat_plot_df <- do.call(rbind, lapply(names(strat_list), function(grp) {
    sub <- strat_list[[grp]]
    if (nrow(sub) == 0) return(NULL)
    data.frame(group = grp, mch_logfc = sub$mch_logfc, stringsAsFactors = FALSE)
  }))
  strat_plot_df$group <- factor(strat_plot_df$group,
                                levels = c("Intragenic only", "Intergenic only",
                                           "Both", "No K119ub peak"))

  p_strat <- ggplot(strat_plot_df, aes(x = group, y = mch_logfc, fill = group)) +
    geom_violin(draw_quantiles = 0.5, scale = "width", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(x = "K119ub peak location relative to gene", y = "mCH edgeR logFC",
         title = "mCH change by K119ub peak location") +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  ann_df <- summarise_groups(strat_plot_df, "group", "mch_logfc")
  ann_text <- group_label(ann_df, digits = 3)
  p_strat <- p_strat +
    geom_text(data = data.frame(group = ann_df$group,
                                y = max(strat_plot_df$mch_logfc, na.rm = TRUE) * 1.05,
                                label = ann_text),
              aes(x = group, y = y, label = label),
              size = 2.5, vjust = 0, inherit.aes = FALSE)

  save_multiformat_ggplot(p_strat, file.path(out_dir, "80_02b_mch_by_k119ub_location"),
                          width = 9, height = 7)

  # Pairwise Wilcoxon
  if (length(strat_list[["Intragenic only"]]$mch_logfc) > 0 &&
      length(strat_list[["Intergenic only"]]$mch_logfc) > 0) {
    wt <- wilcox.test(strat_list[["Intragenic only"]]$mch_logfc,
                      strat_list[["Intergenic only"]]$mch_logfc)
    cat(sprintf("  Intragenic vs Intergenic mCH logFC Wilcoxon: p = %.2e\n", wt$p.value))
  }

  # =========================================================================
  # STEP 4: MeCP2 peak location classification
  # =========================================================================

  cat("\nSTEP 4: MeCP2 peak location classification\n")

  mecp2_raw_path <- file.path(CODE_DIR, "data", "mecp2_diffbind.txt")
  mecp2_raw_df <- read.table(mecp2_raw_path, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, quote = "", fill = TRUE)

  ensure_chr_prefix <- function(x) {
    x <- as.character(x)
    ifelse(grepl("^chr", x), x, paste0("chr", x))
  }

  if ("seqnames" %in% colnames(mecp2_raw_df)) {
    mecp2_raw_df$Chr <- ensure_chr_prefix(mecp2_raw_df$seqnames)
    mecp2_raw_df$peak_start <- mecp2_raw_df$start
    mecp2_raw_df$peak_end <- mecp2_raw_df$end
  } else if ("Summit_Chr" %in% colnames(mecp2_raw_df)) {
    mecp2_raw_df$Chr <- ensure_chr_prefix(mecp2_raw_df$Summit_Chr)
    mecp2_raw_df$peak_start <- mecp2_raw_df$Summit_Start
    mecp2_raw_df$peak_end <- mecp2_raw_df$Summit_End
  } else {
    mecp2_raw_df$Chr <- ensure_chr_prefix(mecp2_raw_df$Chr)
    mecp2_raw_df$peak_start <- mecp2_raw_df$Start
    mecp2_raw_df$peak_end <- mecp2_raw_df$End
  }

  mecp2_raw_df$peak_width <- mecp2_raw_df$peak_end - mecp2_raw_df$peak_start
  mecp2_raw_df$direction <- "Unchanged"
  mecp2_raw_df$direction[mecp2_raw_df$FDR < fdr & mecp2_raw_df$Fold > 0] <- "Gained"
  mecp2_raw_df$direction[mecp2_raw_df$FDR < fdr & mecp2_raw_df$Fold < 0] <- "Lost"

  cat(sprintf("  MeCP2 raw: %s peaks\n", fmt_int(nrow(mecp2_raw_df))))

  mecp2_anno <- classify_peak_location(mecp2_raw_df, "MeCP2")

  mecp2_loc_summary <- mecp2_anno %>%
    dplyr::group_by(location) %>%
    dplyr::summarise(
      n_peaks = dplyr::n(),
      median_width = median(width, na.rm = TRUE),
      n_gained = sum(direction == "Gained"),
      n_lost = sum(direction == "Lost"),
      pct_gained = 100 * n_gained / dplyr::n(),
      .groups = "drop"
    )
  write_section_table(as.data.frame(mecp2_loc_summary),
                      file.path(out_dir, "80_02_mecp2_location_summary.tsv"))

  sig_mecp2 <- mecp2_anno[mecp2_anno$direction %in% c("Gained", "Lost"), , drop = FALSE]
  if (nrow(sig_mecp2) > 0) {
    ct_mecp2 <- table(direction = sig_mecp2$direction, location = sig_mecp2$location)
    fisher_mecp2 <- fisher.test(ct_mecp2)
    cat(sprintf("  MeCP2 direction x location Fisher: OR = %.2f, p = %.2e\n",
                fisher_mecp2$estimate, fisher_mecp2$p.value))
  }

  # =========================================================================
  # STEP 5: K119ub-up / mCH-down loci
  # =========================================================================

  cat("\nSTEP 5: K119ub-up / mCH-down loci characterisation\n")

  marks_path <- HANDOFF_PATHS$gene_level_all_marks
  chrom_path <- HANDOFF_PATHS$chromatin_state
  stopifnot(file.exists(marks_path))
  stopifnot(file.exists(chrom_path))

  marks <- read.table(marks_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                      quote = "")
  chrom <- read.table(chrom_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                      quote = "")

  div_df <- merge(
    data.frame(gene_name = mch_results$gene_name,
               mch_logfc = mch_results$edger_logFC,
               mch_sig = mch_results$mch_sig,
               mch_direction = mch_results$mch_direction,
               stringsAsFactors = FALSE),
    marks[, c("gene_name", "k119ub_fold", "k119ub_fdr", "k119ub_n_peaks")],
    by = "gene_name", all.x = TRUE
  )
  div_df <- merge(div_df,
                  chrom[, c("gene_name", "promoter_state", "body_state")],
                  by = "gene_name", all.x = TRUE)

  div_df$k119ub_gained <- !is.na(div_df$k119ub_fdr) & div_df$k119ub_fdr < fdr &
    div_df$k119ub_fold > 0
  div_df$mch_hypo <- div_df$mch_sig & div_df$mch_logfc < 0
  div_df$mch_hyper <- div_df$mch_sig & div_df$mch_logfc > 0

  # The divergent set: K119ub gained AND mCH down
  divergent <- div_df[div_df$k119ub_gained & div_df$mch_hypo, , drop = FALSE]
  concordant <- div_df[div_df$k119ub_gained & div_df$mch_hyper, , drop = FALSE]

  cat(sprintf("  K119ub gained + mCH hypo: %s genes\n", fmt_int(nrow(divergent))))
  cat(sprintf("  K119ub gained + mCH hyper: %s genes\n", fmt_int(nrow(concordant))))

  # Chromatin state distribution of divergent vs concordant
  if (nrow(divergent) > 0) {
    div_body <- table(divergent$body_state)
    con_body <- if (nrow(concordant) > 0) table(concordant$body_state) else integer(0)

    chrom_comparison <- data.frame(
      body_state = union(names(div_body), names(con_body)),
      stringsAsFactors = FALSE
    )
    chrom_comparison$n_divergent <- as.integer(div_body[chrom_comparison$body_state])
    chrom_comparison$n_concordant <- as.integer(con_body[chrom_comparison$body_state])
    chrom_comparison[is.na(chrom_comparison)] <- 0
    chrom_comparison$pct_divergent <- 100 * chrom_comparison$n_divergent / nrow(divergent)
    chrom_comparison$pct_concordant <- if (nrow(concordant) > 0) {
      100 * chrom_comparison$n_concordant / nrow(concordant)
    } else { NA }

    write_section_table(chrom_comparison,
                        file.path(out_dir, "80_02_divergent_chromatin_state.tsv"))
  }

  # MeCP2 direction at divergent genes
  div_mecp2 <- merge(divergent[, c("gene_name", "mch_logfc", "k119ub_fold")],
                     mecp2_gene, by = "gene_name", all.x = TRUE)
  div_mecp2$mecp2_direction <- "No peak"
  div_mecp2$mecp2_direction[!is.na(div_mecp2$mecp2_has_sig) &
                              div_mecp2$mecp2_has_sig &
                              div_mecp2$mecp2_fold > 0] <- "MeCP2 Up"
  div_mecp2$mecp2_direction[!is.na(div_mecp2$mecp2_has_sig) &
                              div_mecp2$mecp2_has_sig &
                              div_mecp2$mecp2_fold < 0] <- "MeCP2 Down"

  cat(sprintf("  MeCP2 at divergent genes: %s Up, %s Down, %s No peak\n",
              fmt_int(sum(div_mecp2$mecp2_direction == "MeCP2 Up")),
              fmt_int(sum(div_mecp2$mecp2_direction == "MeCP2 Down")),
              fmt_int(sum(div_mecp2$mecp2_direction == "No peak"))))

  write_section_table(div_mecp2, file.path(out_dir, "80_02_divergent_genes.tsv"))

  # Fisher: divergent vs concordant as gene set property
  has_k119ub <- div_df[div_df$k119ub_gained, , drop = FALSE]
  if (nrow(has_k119ub) > 0 && sum(has_k119ub$mch_sig) > 0) {
    gene_df_div <- data.frame(
      gene_name = has_k119ub$gene_name,
      chr = mch_results$chr[match(has_k119ub$gene_name, mch_results$gene_name)],
      mch_hyper = has_k119ub$mch_hyper,
      k119ub_gained = TRUE,
      stringsAsFactors = FALSE
    )
    gene_df_div <- gene_df_div[gene_df_div$mch_sig <- has_k119ub$mch_sig, ]
    gene_df_div <- gene_df_div[has_k119ub$mch_sig, , drop = FALSE]

    register_fisher_test(SECTION_ID, "k119ub_gained_mch_direction",
                         "Among K119ub-gained mCH-significant genes: hyper vs hypo",
                         data.frame(
                           gene_name = has_k119ub$gene_name[has_k119ub$mch_sig],
                           chr = mch_results$chr[match(
                             has_k119ub$gene_name[has_k119ub$mch_sig],
                             mch_results$gene_name)],
                           mch_hyper = has_k119ub$mch_hyper[has_k119ub$mch_sig],
                           k119ub_gained = TRUE,
                           stringsAsFactors = FALSE
                         ),
                         "mch_hyper", "k119ub_gained", out_dir)
  }

  # GO enrichment of divergent genes
  if (nrow(divergent) >= 10) {
    cat("  Running GO enrichment on divergent gene set...\n")
    all_genes_entrez <- tryCatch(
      bitr(mch_results$gene_name, fromType = "SYMBOL", toType = "ENTREZID",
           OrgDb = org.Mm.eg.db),
      error = function(e) NULL
    )
    div_entrez <- tryCatch(
      bitr(divergent$gene_name, fromType = "SYMBOL", toType = "ENTREZID",
           OrgDb = org.Mm.eg.db),
      error = function(e) NULL
    )

    if (!is.null(all_genes_entrez) && !is.null(div_entrez) && nrow(div_entrez) >= 5) {
      ego <- enrichGO(gene = div_entrez$ENTREZID,
                      universe = all_genes_entrez$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.1,
                      readable = TRUE)

      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        ego_df <- as.data.frame(ego)
        write_section_table(ego_df,
                            file.path(out_dir, "80_02_divergent_go_enrichment.tsv"))
        cat(sprintf("    %s GO terms at q < 0.1\n", nrow(ego_df)))

        if (nrow(ego_df) >= 3) {
          p_go <- dotplot(ego, showCategory = min(20, nrow(ego_df))) +
            ggtitle("GO BP: K119ub-gained + mCH-hypo genes") +
            theme_emseq()
          save_multiformat_ggplot(p_go,
                                 file.path(out_dir, "80_02e_divergent_go"),
                                 width = 10, height = 8)
        }
      } else {
        cat("    No GO terms enriched at q < 0.1\n")
      }
    }
  }

  # =========================================================================
  # STEP 6: Consensus peak analysis
  # =========================================================================

  cat("\nSTEP 6: Consensus peak analysis\n")

  gene_gr <- gene_bodies[seqnames(gene_bodies) %in% CANONICAL_CHRS]

  # K119ub consensus
  k119ub_cons_gr <- k119ub_consensus[seqnames(k119ub_consensus) %in% CANONICAL_CHRS]
  k119ub_cons_intra <- subsetByOverlaps(k119ub_cons_gr, gene_gr)
  k119ub_cons_inter <- k119ub_cons_gr[!overlapsAny(k119ub_cons_gr, gene_gr)]

  cat(sprintf("  K119ub consensus: %s total, %s intragenic (%.1f%%), %s intergenic (%.1f%%)\n",
              fmt_int(length(k119ub_cons_gr)),
              fmt_int(length(k119ub_cons_intra)),
              100 * length(k119ub_cons_intra) / length(k119ub_cons_gr),
              fmt_int(length(k119ub_cons_inter)),
              100 * length(k119ub_cons_inter) / length(k119ub_cons_gr)))

  # Enhancer overlap for intergenic consensus peaks
  enhancer_path <- file.path(CODE_DIR, "data", "chromatin", "activeenhancer.bed")
  if (file.exists(enhancer_path)) {
    enhancers <- load_chip_peaks(enhancer_path, "Active enhancer")
    n_enh <- sum(overlapsAny(k119ub_cons_inter, enhancers))
    cat(sprintf("  Intergenic K119ub consensus peaks overlapping active enhancers: %s (%.1f%%)\n",
                fmt_int(n_enh), 100 * n_enh / length(k119ub_cons_inter)))

    cons_summary <- data.frame(
      category = c("K119ub consensus total", "Intragenic", "Intergenic",
                   "Intergenic at enhancer"),
      n_peaks = c(length(k119ub_cons_gr), length(k119ub_cons_intra),
                  length(k119ub_cons_inter), n_enh),
      stringsAsFactors = FALSE
    )
  } else {
    cat("  Active enhancer BED not found, skipping enhancer overlap\n")
    cons_summary <- data.frame(
      category = c("K119ub consensus total", "Intragenic", "Intergenic"),
      n_peaks = c(length(k119ub_cons_gr), length(k119ub_cons_intra),
                  length(k119ub_cons_inter)),
      stringsAsFactors = FALSE
    )
  }
  write_section_table(cons_summary,
                      file.path(out_dir, "80_02_consensus_peak_summary.tsv"))

  # MeCP2 consensus
  mecp2_cons_gr <- mecp2_consensus[seqnames(mecp2_consensus) %in% CANONICAL_CHRS]
  mecp2_cons_intra <- subsetByOverlaps(mecp2_cons_gr, gene_gr)
  mecp2_cons_inter <- mecp2_cons_gr[!overlapsAny(mecp2_cons_gr, gene_gr)]

  cat(sprintf("  MeCP2 consensus: %s total, %s intragenic (%.1f%%), %s intergenic (%.1f%%)\n",
              fmt_int(length(mecp2_cons_gr)),
              fmt_int(length(mecp2_cons_intra)),
              100 * length(mecp2_cons_intra) / length(mecp2_cons_gr),
              fmt_int(length(mecp2_cons_inter)),
              100 * length(mecp2_cons_inter) / length(mecp2_cons_gr)))

  # Distance from intergenic K119ub consensus peaks to nearest gene
  if (length(k119ub_cons_inter) > 0) {
    nearest_idx <- nearest(k119ub_cons_inter, gene_gr)
    has_nearest <- !is.na(nearest_idx)
    if (sum(!has_nearest) > 0) {
      cat(sprintf("  %s intergenic peaks on chromosomes with no genes (dropped from distance calc)\n",
                  fmt_int(sum(!has_nearest))))
    }
    dists <- distance(k119ub_cons_inter[has_nearest], gene_gr[nearest_idx[has_nearest]])
    cat(sprintf("  Intergenic K119ub peak distance to nearest gene: median %s bp, max %s bp (n = %s)\n",
                fmt_int(median(dists, na.rm = TRUE)),
                fmt_int(max(dists, na.rm = TRUE)),
                fmt_int(length(dists))))

    dist_summary <- data.frame(
      metric = c("n_peaks", "n_with_nearest", "n_no_nearest",
                 "median_distance_bp", "mean_distance_bp",
                 "q25_distance_bp", "q75_distance_bp", "max_distance_bp"),
      value = c(length(nearest_idx), sum(has_nearest), sum(!has_nearest),
                median(dists, na.rm = TRUE), mean(dists, na.rm = TRUE),
                quantile(dists, 0.25, na.rm = TRUE), quantile(dists, 0.75, na.rm = TRUE),
                max(dists, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
    write_section_table(dist_summary,
                        file.path(out_dir, "80_02_intergenic_distance_to_gene.tsv"))
  }

  cat(sprintf("\n=== Section 80_02 complete ===\n"))
}

main()
