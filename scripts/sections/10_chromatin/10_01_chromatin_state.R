# scripts/sections/10_chromatin/10_01_chromatin_state.R
#
# Section 10_01 -- chromatin state of the promoter and of the gene body for
# every gene with an mCH measurement.
#
# The analysis asks whether genes that change non-CG methylation in the BAP1
# knockout sit in particular chromatin states. Each gene carries two states.
# The promoter state comes from histone mark overlaps over the TSS window
# (gene_promoters) and takes a level of PROMOTER_STATE_ORDER. The gene-body
# state comes from the same marks over the gene body (gene_bodies) and takes a
# level of BODY_STATE_ORDER. A gene can hold an active promoter over a
# Polycomb-marked body, so the two states are reported separately and their
# cross-tabulation is one of the figures.
#
# Reads (all through the shared config, nothing is re-read from disk here):
#   CHIP_PATHS      CTCF, H3K27ac, H3K27me3, H3K4me1, H3K4me3, bivalent BEDs
#   mch_results     gene-level mCH differential results
#   gene_bodies     GRanges of the same genes, in the same row order
#   gene_promoters  TSS +/- TSS_THRESHOLD windows, in the same row order
#
# Writes into OUT_DIR (default results/sections/10_chromatin/):
#   gene_chromatin_state.tsv        handoff table read by sections 10_03,
#                                   20_02, 40_02 and 40_03
#   gene_chromatin_state_full.tsv   the same genes with their mCH statistics
#   state_distribution_by_gene_set.tsv, state_significance_rate.tsv,
#   state_enrichment_fisher.tsv, mark_overlap_summary.tsv,
#   promoter_body_state_crosstab.tsv, mch_diff_by_state.tsv,
#   gene_length_by_state.tsv, key_genes_chromatin_state.tsv
#   twelve figures in PDF, SVG, PNG, and JPEG
#   fisher_tables/ plus one row per test in the shared Fisher registry
#
# Adapted from the Biomodal script section_10_chromatin_state.R. That script
# classified significant 5mC DMRs; this one classifies genes by mCH. The
# Biomodal panels 10d and 10e select mC-up / hmC-down gene sets, which need two
# methylation modalities, and are not ported.

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)

# =============================================================================
# SECTION CONSTANTS
# =============================================================================

SECTION_ID <- "10_01"

# Names on the left are the CHIP_PATHS keys and the prefixes that
# compute_chip_overlaps() turns into the *_overlap columns that
# classify_promoter_state() and classify_body_state() read. Names on the right
# are for display.
MARK_DISPLAY <- c(
  ctcf     = "CTCF",
  h3k27ac  = "H3K27ac",
  h3k27me3 = "H3K27me3",
  h3k4me1  = "H3K4me1",
  h3k4me3  = "H3K4me3",
  bivalent = "Bivalent"
)

MARK_KEYS <- names(MARK_DISPLAY)

# Overlap columns are stored twice per gene: prom_ over the promoter window,
# body_ over the gene body.
PROM_OVERLAP_COLS <- paste0("prom_", MARK_KEYS, "_overlap")
BODY_OVERLAP_COLS <- paste0("body_", MARK_KEYS, "_overlap")

# The handoff schema read by sections 10_03, 20_02, 40_02 and 40_03.
HANDOFF_COLS <- c("gene_name", "gene_id", "chr", "start", "end",
                  "promoter_state", "body_state",
                  PROM_OVERLAP_COLS, BODY_OVERLAP_COLS)

# The two regions each mark is tested over.
REGIONS <- list(
  list(id = "prom", label = "Promoter window"),
  list(id = "body", label = "Gene body")
)

REGION_LABELS <- vapply(REGIONS, function(r) r$label, character(1))

# The two state columns. Every summary table and every state figure is built
# once per entry.
STATE_KINDS <- list(
  list(id     = "promoter",
       col    = "promoter_state",
       levels = PROMOTER_STATE_ORDER,
       colors = PROMOTER_STATE_COLORS,
       title  = "Promoter state",
       phrase = "promoter state"),
  list(id     = "body",
       col    = "body_state",
       levels = BODY_STATE_ORDER,
       colors = BODY_STATE_COLORS,
       title  = "Gene-body state",
       phrase = "gene-body state")
)

STATE_KIND_TITLES <- vapply(STATE_KINDS, function(k) k$title, character(1))

# The promoter and body level names do not collide, so one level order serves
# the tables that hold both kinds.
STATE_LEVELS_ALL <- c(PROMOTER_STATE_ORDER, BODY_STATE_ORDER)

GENE_SET_ORDER <- c("All tested genes", "Significant mCH",
                    "Hypermethylated", "Hypomethylated")

# Grey for the background set, red for the significant set, taken from the
# shared significance palette.
GENE_SET_COMPARE_COLORS <- c(
  "All tested genes" = unname(COLORS$significant["Not Significant"]),
  "Significant mCH"  = unname(COLORS$significant["Significant"])
)

# =============================================================================
# SMALL HELPERS
# =============================================================================

#' Format an integer for a plot label or a log line.
fmt_int <- function(x) format(x, big.mark = ",", trim = TRUE)

#' Overlap column name for one mark in one region.
overlap_col <- function(region_id, mark_key) {
  paste0(region_id, "_", mark_key, "_overlap")
}

# =============================================================================
# OPTIONS AND INPUT CHECKS
# =============================================================================

parse_options <- function() {
  option_list <- list(
    make_option("--output-dir", dest = "output_dir", type = "character",
                default = OUTPUT_PATHS$chromatin,
                help = "Directory for figures and tables [default: %default]")
  )
  parse_args(OptionParser(
    option_list = option_list,
    description = "Promoter and gene-body chromatin state of genes with differential mCH."
  ))
}

#' Stop unless all six histone mark BED files exist.
check_chip_inputs <- function() {
  paths <- unlist(CHIP_PATHS[MARK_KEYS])
  missing <- names(paths)[!file.exists(paths)]
  if (length(missing) > 0) {
    stop("Histone mark BED files not found for: ",
         paste(missing, collapse = ", "), "\n",
         "  Expected in: ", dirname(CHIP_PATHS$ctcf))
  }
  invisible(TRUE)
}

#' Stop unless gene_bodies, gene_promoters and mch_results describe the same
#' genes in the same order, and unless both methylation directions are present.
check_data_assumptions <- function() {
  if (length(gene_bodies) != nrow(mch_results)) {
    stop("gene_bodies has ", length(gene_bodies), " ranges but mch_results has ",
         nrow(mch_results), " rows.")
  }
  if (length(gene_promoters) != nrow(mch_results)) {
    stop("gene_promoters has ", length(gene_promoters),
         " ranges but mch_results has ", nrow(mch_results), " rows.")
  }
  if (!identical(as.character(mcols(gene_bodies)$gene_name),
                 as.character(mch_results$gene_name))) {
    stop("gene_bodies and mch_results are not in the same gene order.")
  }
  if (!identical(as.character(mcols(gene_promoters)$gene_name),
                 as.character(mch_results$gene_name))) {
    stop("gene_promoters and mch_results are not in the same gene order.")
  }
  if (sum(mch_results$mch_hyper) == 0 || sum(mch_results$mch_hypo) == 0) {
    stop("The direction analyses need both hypermethylated and hypomethylated ",
         "genes. Found ", sum(mch_results$mch_hyper), " hyper and ",
         sum(mch_results$mch_hypo), " hypo.")
  }
  invisible(TRUE)
}

# =============================================================================
# CLASSIFICATION
# =============================================================================

#' Load the six histone mark peak sets as a named list of GRanges.
load_histone_marks <- function() {
  cat("Loading histone mark peak sets...\n")
  peaks <- lapply(MARK_KEYS, function(mark) {
    load_chip_peaks(CHIP_PATHS[[mark]], MARK_DISPLAY[[mark]])
  })
  names(peaks) <- MARK_KEYS
  peaks
}

#' Mark overlaps over one query GRanges, logged per mark.
region_overlaps <- function(query_gr, chip_peaks, region_label) {
  cat(sprintf("\nComputing histone mark overlaps over %s...\n",
              tolower(region_label)))
  overlaps <- compute_chip_overlaps(query_gr, chip_peaks)
  for (mark in MARK_KEYS) {
    col <- paste0(mark, "_overlap")
    cat(sprintf("  %-9s %s genes (%.1f%%)\n", MARK_DISPLAY[[mark]],
                fmt_int(sum(overlaps[[col]])), 100 * mean(overlaps[[col]])))
  }
  overlaps
}

#' Build the gene-level table: mCH statistics, promoter and gene-body mark
#' overlaps, and the two chromatin states, one row per tested gene.
build_gene_table <- function(chip_peaks) {
  prom_overlaps <- region_overlaps(gene_promoters, chip_peaks,
                                   "promoter windows")
  body_overlaps <- region_overlaps(gene_bodies, chip_peaks, "gene bodies")

  cat("\nClassifying promoter and gene-body chromatin states...\n")
  promoter_state <- classify_promoter_state(prom_overlaps)
  body_state <- classify_body_state(body_overlaps)

  colnames(prom_overlaps) <- paste0("prom_", colnames(prom_overlaps))
  colnames(body_overlaps) <- paste0("body_", colnames(body_overlaps))

  genes <- data.frame(
    gene_name      = mch_results$gene_name,
    gene_id        = mch_results$gene_id,
    chr            = mch_results$chr,
    start          = mch_results$start,
    end            = mch_results$end,
    gene_length    = mch_results$gene_length,
    mch_ctrl       = mch_results$mch_ctrl,
    mch_mut        = mch_results$mch_mut,
    mch_diff       = mch_results$mch_diff,
    edger_logFC    = mch_results$edger_logFC,
    edger_fdr      = mch_results$edger_fdr,
    mch_sig        = mch_results$mch_sig,
    mch_hyper      = mch_results$mch_hyper,
    mch_hypo       = mch_results$mch_hypo,
    mch_direction  = mch_results$mch_direction,
    promoter_state = promoter_state,
    body_state     = body_state,
    stringsAsFactors = FALSE
  )
  genes <- cbind(genes, prom_overlaps, body_overlaps)

  for (kind in STATE_KINDS) {
    cat(sprintf("\n  %s of all tested genes:\n", kind$title))
    for (s in kind$levels) {
      n <- sum(genes[[kind$col]] == s)
      cat(sprintf("    %-20s %6s (%5.1f%%)\n", s, fmt_int(n),
                  100 * n / nrow(genes)))
    }
  }
  genes
}

# =============================================================================
# SUMMARY TABLES
# =============================================================================

#' Logical row index for each of the four gene sets.
gene_set_index <- function(genes) {
  list(
    "All tested genes" = rep(TRUE, nrow(genes)),
    "Significant mCH"  = genes$mch_sig,
    "Hypermethylated"  = genes$mch_hyper,
    "Hypomethylated"   = genes$mch_hypo
  )
}

#' Count and percentage of each state within each gene set, for one state kind.
summarise_state_distribution <- function(genes, kind) {
  idx <- gene_set_index(genes)
  rows <- lapply(names(idx), function(set_name) {
    tab <- table(genes[[kind$col]][idx[[set_name]]])
    data.frame(
      state_kind = kind$title,
      gene_set   = set_name,
      state      = names(tab),
      count      = as.integer(tab),
      n_set      = sum(tab),
      percentage = 100 * as.integer(tab) / sum(tab),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$gene_set <- factor(out$gene_set, levels = GENE_SET_ORDER)
  out$state <- factor(out$state, levels = kind$levels)
  out[order(out$gene_set, out$state), ]
}

#' Gene counts and significance rate within each state, for one state kind.
summarise_significance_rate <- function(genes, kind) {
  rows <- lapply(kind$levels, function(s) {
    sub <- genes[genes[[kind$col]] == s, , drop = FALSE]
    data.frame(
      state_kind = kind$title,
      state      = s,
      n_genes    = nrow(sub),
      n_sig      = sum(sub$mch_sig),
      n_hyper    = sum(sub$mch_hyper),
      n_hypo     = sum(sub$mch_hypo),
      pct_sig    = 100 * sum(sub$mch_sig) / nrow(sub),
      pct_hyper  = 100 * sum(sub$mch_hyper) / nrow(sub),
      pct_hypo   = 100 * sum(sub$mch_hypo) / nrow(sub),
      median_mch_diff = median(sub$mch_diff, na.rm = TRUE),
      median_edger_logFC = median(sub$edger_logFC, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$state <- factor(out$state, levels = kind$levels)
  out
}

#' Percentage of genes in each gene set overlapping each mark, for one region.
summarise_mark_overlap <- function(genes, region) {
  cols <- vapply(MARK_KEYS, function(m) overlap_col(region$id, m), character(1))
  idx <- gene_set_index(genes)
  rows <- lapply(names(idx), function(set_name) {
    sub <- genes[idx[[set_name]], , drop = FALSE]
    data.frame(
      region      = region$label,
      gene_set    = set_name,
      n_genes     = nrow(sub),
      mark        = unname(MARK_DISPLAY),
      n_overlap   = vapply(cols, function(col) sum(sub[[col]]), integer(1)),
      pct_overlap = vapply(cols, function(col) 100 * mean(sub[[col]]), numeric(1)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })
  out <- do.call(rbind, rows)
  out$gene_set <- factor(out$gene_set, levels = GENE_SET_ORDER)
  out$mark <- factor(out$mark, levels = unname(MARK_DISPLAY))
  out
}

#' Gene counts and mCH outcome for every promoter-state by body-state pair.
#'
#' This is the table behind the cross-tabulation heatmap. Every pair of levels
#' appears, including pairs that hold no gene.
summarise_state_crosstab <- function(genes) {
  out <- genes %>%
    dplyr::group_by(promoter_state, body_state, .drop = FALSE) %>%
    dplyr::summarise(
      n_genes = dplyr::n(),
      n_sig   = sum(mch_sig),
      n_hyper = sum(mch_hyper),
      n_hypo  = sum(mch_hypo),
      median_mch_diff = median(mch_diff),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      pct_of_tested = 100 * n_genes / nrow(genes),
      pct_sig = ifelse(n_genes > 0, 100 * n_sig / n_genes, NA_real_)
    ) %>%
    as.data.frame()
  out$promoter_state <- factor(as.character(out$promoter_state),
                               levels = PROMOTER_STATE_ORDER)
  out$body_state <- factor(as.character(out$body_state),
                           levels = BODY_STATE_ORDER)
  out
}

#' Stack the per-kind tables of a summariser into one frame.
stack_over_kinds <- function(genes, summariser) {
  rows <- lapply(STATE_KINDS, function(kind) {
    df <- summariser(genes, kind)
    df$state <- as.character(df$state)
    df
  })
  out <- do.call(rbind, rows)
  out$state_kind <- factor(out$state_kind, levels = STATE_KIND_TITLES)
  out$state <- factor(out$state, levels = STATE_LEVELS_ALL)
  out
}

# =============================================================================
# FISHER TESTS
# =============================================================================

#' One gene-level Fisher test per promoter state and per gene-body state.
#'
#' Each test asks whether genes with significant differential mCH are enriched
#' in one state against every other tested gene. A state that holds no gene
#' carries no testable 2x2 table and is reported instead of tested.
run_state_fisher_tests <- function(genes, out_dir) {
  rows <- list()

  for (kind in STATE_KINDS) {
    state_counts <- table(genes[[kind$col]])
    tested_states <- names(state_counts)[state_counts > 0]
    empty_states <- names(state_counts)[state_counts == 0]
    if (length(empty_states) > 0) {
      cat(sprintf("  No gene falls in %s: %s. No Fisher test for these states.\n",
                  kind$phrase, paste(empty_states, collapse = ", ")))
    }

    for (state in tested_states) {
      df <- data.frame(gene_name = genes$gene_name, chr = genes$chr,
                       mch_sig = genes$mch_sig, stringsAsFactors = FALSE)
      df$in_state <- genes[[kind$col]] == state

      ft <- register_fisher_test(
        section = SECTION_ID,
        test_id = paste0(kind$id, "_", state),
        description = sprintf(
          paste("Genes with significant differential mCH versus all other",
                "tested genes, enrichment in %s %s."),
          kind$phrase, state),
        gene_df = df,
        row_var = "mch_sig",
        col_var = "in_state",
        output_dir = out_dir
      )

      in_set <- df$mch_sig
      rows[[length(rows) + 1]] <- data.frame(
        state_kind      = kind$title,
        state           = state,
        n_genes         = nrow(df),
        n_significant   = sum(in_set),
        n_in_state      = sum(df$in_state),
        n_sig_in_state  = sum(in_set & df$in_state),
        pct_of_sig_in_state        = 100 * sum(in_set & df$in_state) / sum(in_set),
        pct_of_background_in_state = 100 * sum(!in_set & df$in_state) / sum(!in_set),
        odds_ratio = unname(ft$estimate),
        ci_low     = ft$conf.int[1],
        ci_high    = ft$conf.int[2],
        p_value    = ft$p.value,
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, rows)
  out$p_adj_bh <- p.adjust(out$p_value, method = "BH")
  out$state_kind <- factor(out$state_kind, levels = STATE_KIND_TITLES)
  out$state <- factor(out$state, levels = STATE_LEVELS_ALL)
  out[order(out$state_kind, out$state), ]
}

# =============================================================================
# FIGURES
# =============================================================================

#' Two per-kind panels side by side under one title.
combine_kind_panels <- function(panels, title, subtitle = NULL,
                                widths = c(1, 1)) {
  (panels[[1]] | panels[[2]]) +
    plot_layout(widths = widths) +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 11)
      )
    )
}

#' Bar of the state distribution over all tested genes, for one state kind.
plot_distribution_panel <- function(distribution, kind) {
  d <- distribution[distribution$gene_set == "All tested genes", ]
  d$state <- factor(as.character(d$state), levels = kind$levels)

  ggplot(d, aes(x = state, y = percentage, fill = state)) +
    geom_col(color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%\n(n = %s)", percentage,
                                  fmt_int(count))),
              vjust = -0.2, size = 2.8, lineheight = 0.9) +
    scale_fill_manual(values = kind$colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, max(d$percentage) * 1.28),
                       expand = c(0, 0)) +
    labs(title = kind$title,
         subtitle = sprintf("All %s tested genes", fmt_int(unique(d$n_set))),
         x = kind$title, y = "Percentage of tested genes (%)") +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
}

#' Bar of the state distribution per mCH direction, beside a pie of the
#' distribution over all significant genes, for one state kind.
plot_direction_panel <- function(distribution, kind) {
  by_dir <- distribution[distribution$gene_set %in%
                           c("Hypermethylated", "Hypomethylated"), ]
  by_dir$gene_set <- droplevels(by_dir$gene_set)
  by_dir$state <- factor(as.character(by_dir$state), levels = kind$levels)
  sig <- distribution[distribution$gene_set == "Significant mCH", ]
  sig$state <- factor(as.character(sig$state), levels = kind$levels)

  n_hyper <- unique(by_dir$n_set[by_dir$gene_set == "Hypermethylated"])
  n_hypo  <- unique(by_dir$n_set[by_dir$gene_set == "Hypomethylated"])

  p_bar <- ggplot(by_dir, aes(x = state, y = percentage, fill = state)) +
    geom_col(color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%\n(n = %s)", percentage,
                                  fmt_int(count))),
              vjust = -0.2, size = 2.5, lineheight = 0.9) +
    facet_wrap(~gene_set, ncol = 2) +
    scale_fill_manual(values = kind$colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, max(by_dir$percentage) * 1.3),
                       expand = c(0, 0)) +
    labs(
      title = sprintf("%s by mCH direction", kind$title),
      subtitle = sprintf("Hypermethylated n = %s | Hypomethylated n = %s",
                         fmt_int(n_hyper), fmt_int(n_hypo)),
      x = kind$title, y = "Percentage of the gene set (%)"
    ) +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  p_pie <- ggplot(sig, aes(x = "", y = percentage, fill = state)) +
    geom_col(width = 1, color = "white", linewidth = 0.3) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = kind$colors, drop = FALSE, name = kind$title) +
    geom_text(aes(label = ifelse(percentage > 4,
                                 sprintf("%.1f%%", percentage), "")),
              position = position_stack(vjust = 0.5), size = 2.6) +
    labs(title = sprintf("All significant mCH genes (n = %s)",
                         fmt_int(unique(sig$n_set)))) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          legend.position = "right")

  (p_bar | p_pie) + plot_layout(widths = c(2, 1))
}

#' Grouped bar comparing significant genes against all tested genes.
plot_significant_vs_all_panel <- function(distribution, kind) {
  cmp <- distribution[distribution$gene_set %in%
                        c("All tested genes", "Significant mCH"), ]
  cmp$gene_set <- droplevels(cmp$gene_set)
  cmp$state <- factor(as.character(cmp$state), levels = kind$levels)
  set_sizes <- unique(cmp[, c("gene_set", "n_set")])

  ggplot(cmp, aes(x = state, y = percentage, fill = gene_set)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72,
             color = "black", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)),
              position = position_dodge(width = 0.8), vjust = -0.35,
              size = 2.7) +
    scale_fill_manual(values = GENE_SET_COMPARE_COLORS, name = "Gene set") +
    scale_y_continuous(limits = c(0, max(cmp$percentage) * 1.18),
                       expand = c(0, 0)) +
    labs(
      title = kind$title,
      subtitle = paste(sprintf("%s n = %s", as.character(set_sizes$gene_set),
                               fmt_int(set_sizes$n_set)), collapse = " | "),
      x = kind$title, y = "Percentage of the gene set (%)"
    ) +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
}

#' Stacked proportion bar of state for the four gene sets.
plot_stacked_proportion_panel <- function(distribution, kind) {
  d <- distribution
  d$state <- factor(as.character(d$state), levels = kind$levels)
  set_sizes <- unique(d[, c("gene_set", "n_set")])

  ggplot(d, aes(x = gene_set, y = percentage, fill = state)) +
    geom_col(position = "stack", color = "white", linewidth = 0.25) +
    geom_text(aes(label = ifelse(percentage > 3,
                                 sprintf("%.1f%%", percentage), "")),
              position = position_stack(vjust = 0.5), size = 3,
              fontface = "bold") +
    geom_text(data = set_sizes,
              aes(x = gene_set, y = 103,
                  label = sprintf("n = %s", fmt_int(n_set))),
              inherit.aes = FALSE, size = 3.2, fontface = "italic") +
    scale_fill_manual(values = kind$colors, drop = FALSE, name = kind$title) +
    scale_y_continuous(limits = c(0, 108), expand = c(0, 0)) +
    labs(title = kind$title, x = "Gene set",
         y = "Percentage of the gene set (%)") +
    theme_emseq() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"))
}

#' Percentage of genes that reach mCH significance within each state, with the
#' significant and total gene counts printed on the bars.
plot_significance_rate_panel <- function(rate, kind) {
  drawn <- rate[rate$n_genes > 0, ]
  drawn$state <- factor(as.character(drawn$state), levels = kind$levels)

  ggplot(drawn, aes(x = state, y = pct_sig, fill = state)) +
    geom_col(color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%\n%s / %s", pct_sig,
                                  fmt_int(n_sig), fmt_int(n_genes))),
              vjust = -0.15, size = 2.8, lineheight = 0.9) +
    scale_fill_manual(values = kind$colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, max(drawn$pct_sig) * 1.3),
                       expand = c(0, 0)) +
    labs(title = kind$title, x = kind$title,
         y = "Genes with significant mCH (%)") +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
}

#' Heatmap of promoter state against gene-body state.
#'
#' The tile fill is the share of all tested genes in the pair. The text gives
#' the gene count, that share, and the mCH significance rate inside the pair.
plot_state_crosstab <- function(crosstab, n_genes) {
  ggplot(crosstab, aes(x = body_state, y = promoter_state,
                       fill = pct_of_tested)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("n = %s\n%.1f%% of genes\n%s",
                                  fmt_int(n_genes),
                                  pct_of_tested,
                                  ifelse(is.na(pct_sig), "no genes",
                                         sprintf("%.1f%% sig", pct_sig)))),
              size = 2.7, lineheight = 0.95) +
    scale_fill_gradient(low = "white", high = "#3182bd",
                        name = "% of tested genes") +
    scale_y_discrete(limits = rev(PROMOTER_STATE_ORDER)) +
    labs(
      title = "Promoter state against gene-body state",
      subtitle = sprintf(paste("All %s tested genes. The two states are",
                               "classified over different intervals, so a gene",
                               "can hold any pair."),
                         fmt_int(n_genes)),
      x = "Gene-body state", y = "Promoter state"
    ) +
    theme_emseq() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"))
}

#' Heatmap of histone mark overlap percentage per gene set, for one region.
plot_mark_overlap_panel <- function(mark_overlap, region) {
  d <- mark_overlap[mark_overlap$region == region$label, ]

  ggplot(d, aes(x = mark, y = gene_set, fill = pct_overlap)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%\n(%s)", pct_overlap,
                                  fmt_int(n_overlap))),
              size = 2.8, lineheight = 0.9) +
    scale_fill_gradient2(low = "white", mid = "#fee090", high = "#d73027",
                         midpoint = 50, name = "% of genes",
                         limits = c(0, 100)) +
    scale_y_discrete(limits = rev(GENE_SET_ORDER)) +
    labs(title = region$label,
         subtitle = "Percentage of genes in each set whose interval overlaps the mark",
         x = "Histone mark", y = "") +
    theme_emseq() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(face = "bold"))
}

#' Odds ratio per state on a log2 axis, for one state kind.
plot_state_enrichment_panel <- function(enrichment, kind) {
  d <- enrichment[enrichment$state_kind == kind$title, ]
  drawn <- d[is.finite(d$odds_ratio) & d$odds_ratio > 0 &
               is.finite(d$ci_low) & d$ci_low > 0 &
               is.finite(d$ci_high), ]
  drawn$state <- factor(as.character(drawn$state), levels = kind$levels)
  drawn$log2_or <- log2(drawn$odds_ratio)
  drawn$log2_low <- log2(drawn$ci_low)
  drawn$log2_high <- log2(drawn$ci_high)
  drawn$star <- ifelse(drawn$p_adj_bh < Q_THRESHOLD, "*", "")
  n_dropped <- nrow(d) - nrow(drawn)

  ggplot(drawn, aes(x = state, y = log2_or, color = state)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_errorbar(aes(ymin = log2_low, ymax = log2_high), width = 0.22,
                  linewidth = 0.5) +
    geom_point(size = 2.6) +
    geom_text(aes(y = log2_high, label = star), vjust = -0.3, size = 5,
              show.legend = FALSE) +
    scale_color_manual(values = kind$colors, drop = FALSE) +
    coord_flip() +
    labs(
      title = kind$title,
      subtitle = sprintf(paste("* marks BH q < %.2f.",
                               "%d test(s) with a non-finite odds ratio or",
                               "interval are in the table only."),
                         Q_THRESHOLD, n_dropped),
      x = kind$title, y = "log2 odds ratio"
    ) +
    theme_emseq() +
    theme(legend.position = "none")
}

#' Violin and box of a numeric gene-level value across the states of one kind,
#' with n and the median printed above each group.
plot_value_by_state_panel <- function(genes, value_col, stats, kind, y_label) {
  values <- genes[[value_col]]
  span <- diff(range(values, na.rm = TRUE))
  label_y <- max(values, na.rm = TRUE) + 0.02 * span

  ggplot(genes, aes(x = .data[[kind$col]], y = .data[[value_col]],
                    fill = .data[[kind$col]])) +
    geom_violin(scale = "width", trim = FALSE, color = "grey30",
                linewidth = 0.3, alpha = 0.75) +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white",
                 color = "grey20", linewidth = 0.3) +
    geom_text(data = stats,
              aes(x = .data[[kind$col]], y = label_y, label = label),
              inherit.aes = FALSE, size = 2.7, vjust = 0, lineheight = 0.9) +
    scale_fill_manual(values = kind$colors, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    labs(title = kind$title, x = kind$title, y = y_label) +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
}

#' Horizontal bar of the mCH log fold change for the project key genes,
#' coloured by gene-body state and labelled with both states.
plot_key_genes <- function(key_genes) {
  key_genes <- key_genes[order(key_genes$edger_logFC), ]
  key_genes$plot_label <- factor(key_genes$plot_label,
                                 levels = key_genes$plot_label)

  ggplot(key_genes, aes(x = plot_label, y = edger_logFC, fill = body_state)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%s / %s  (FDR = %.2g)",
                                  as.character(promoter_state),
                                  as.character(body_state), edger_fdr),
                  hjust = ifelse(edger_logFC >= 0, -0.05, 1.05)),
              size = 2.7, color = "grey20") +
    geom_hline(yintercept = 0, color = "grey40") +
    scale_fill_manual(values = BODY_STATE_COLORS, drop = FALSE,
                      name = "Gene-body state") +
    scale_y_continuous(expand = expansion(mult = c(0.45, 0.45))) +
    coord_flip() +
    labs(
      title = "Chromatin state of the key genes",
      subtitle = sprintf(paste("%d of %d key genes are in the tested set.",
                               "Labels give promoter state / gene-body state."),
                         nrow(key_genes), length(KEY_GENES)),
      x = "", y = "mCH edgeR log fold change (mutant over control)"
    ) +
    theme_emseq() +
    theme(legend.position = "bottom")
}

#' Gene length against mCH log fold change, coloured by one state kind, with
#' the key genes labelled.
plot_length_vs_logfc_panel <- function(genes, key_genes, kind) {
  ggplot(genes, aes(x = log10(gene_length), y = edger_logFC,
                    color = .data[[kind$col]])) +
    geom_point(size = 0.5, alpha = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(data = key_genes, size = 2, alpha = 1, color = "black") +
    ggrepel::geom_text_repel(data = key_genes, aes(label = gene_name),
                             color = "black", size = 3, min.segment.length = 0,
                             box.padding = 0.4, max.overlaps = Inf,
                             show.legend = FALSE) +
    scale_color_manual(values = kind$colors, drop = FALSE, name = kind$title) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(
      title = kind$title,
      subtitle = sprintf("%s tested genes. Key genes are labelled in black.",
                         fmt_int(nrow(genes))),
      x = "log10 gene length (bp)",
      y = "mCH edgeR log fold change (mutant over control)"
    ) +
    theme_emseq() +
    theme(legend.position = "bottom")
}

# =============================================================================
# OUTPUT
# =============================================================================

#' Write the cross-section handoff table read by 10_03, 20_02, 40_02 and 40_03.
#'
#' The column set and its order are the handoff contract, so the frame is
#' checked against HANDOFF_COLS before it is written.
write_handoff_table <- function(genes) {
  missing <- setdiff(HANDOFF_COLS, colnames(genes))
  if (length(missing) > 0) {
    stop("Gene table is missing handoff columns: ",
         paste(missing, collapse = ", "))
  }

  handoff <- genes[, HANDOFF_COLS, drop = FALSE]
  if (!identical(colnames(handoff), HANDOFF_COLS)) {
    stop("Handoff table carries columns ",
         paste(colnames(handoff), collapse = ", "), "\n  Expected: ",
         paste(HANDOFF_COLS, collapse = ", "))
  }

  write_section_table(handoff, HANDOFF_PATHS$chromatin_state)
  cat(sprintf("  Handoff path: %s\n", HANDOFF_PATHS$chromatin_state))
  invisible(handoff)
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opt <- parse_options()
  out_dir <- opt$output_dir

  cat("================================================================================\n")
  cat("SECTION 10_01: PROMOTER AND GENE-BODY CHROMATIN STATE WITH DIFFERENTIAL mCH\n")
  cat("================================================================================\n")
  cat("Output dir:      ", out_dir, "\n", sep = "")
  cat("Promoter window: TSS +/- ", TSS_THRESHOLD, " bp\n\n", sep = "")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  check_chip_inputs()
  check_data_assumptions()

  chip_peaks <- load_histone_marks()
  genes <- build_gene_table(chip_peaks)
  genes$log10_gene_length <- log10(genes$gene_length)

  cat("\nBuilding summary tables...\n")
  distribution <- stack_over_kinds(genes, summarise_state_distribution)
  rate <- stack_over_kinds(genes, summarise_significance_rate)
  mark_overlap <- do.call(rbind, lapply(REGIONS, function(region) {
    summarise_mark_overlap(genes, region)
  }))
  mark_overlap$region <- factor(mark_overlap$region, levels = REGION_LABELS)
  crosstab <- summarise_state_crosstab(genes)

  diff_stats <- lapply(STATE_KINDS, function(kind) {
    summarise_groups(genes, kind$col, "mch_diff")
  })
  length_stats <- lapply(STATE_KINDS, function(kind) {
    summarise_groups(genes, kind$col, "log10_gene_length")
  })
  names(diff_stats) <- vapply(STATE_KINDS, function(k) k$id, character(1))
  names(length_stats) <- names(diff_stats)

  cat("\nRunning gene-level Fisher tests...\n")
  enrichment <- run_state_fisher_tests(genes, out_dir)

  key_genes <- genes[genes$gene_name %in% KEY_GENES, , drop = FALSE]
  if (nrow(key_genes) == 0) {
    stop("None of the KEY_GENES are in the tested gene set: ",
         paste(KEY_GENES, collapse = ", "))
  }
  repeated <- key_genes$gene_name[duplicated(key_genes$gene_name)]
  key_genes$plot_label <- ifelse(key_genes$gene_name %in% repeated,
                                 paste0(key_genes$gene_name, " (",
                                        key_genes$gene_id, ")"),
                                 key_genes$gene_name)

  cat("\nWriting tables...\n")
  write_handoff_table(genes)
  write_section_table(genes, file.path(out_dir, "gene_chromatin_state_full.tsv"))
  write_section_table(distribution,
                      file.path(out_dir, "state_distribution_by_gene_set.tsv"))
  write_section_table(rate, file.path(out_dir, "state_significance_rate.tsv"))
  write_section_table(enrichment,
                      file.path(out_dir, "state_enrichment_fisher.tsv"))
  write_section_table(mark_overlap,
                      file.path(out_dir, "mark_overlap_summary.tsv"))
  write_section_table(crosstab,
                      file.path(out_dir, "promoter_body_state_crosstab.tsv"))

  # Both tables take only the data columns of the group summaries. The two-line
  # figure text is built at the plot, so no column here holds a newline.
  diff_table <- do.call(rbind, lapply(STATE_KINDS, function(kind) {
    df <- diff_stats[[kind$id]]
    data.frame(state_kind = kind$title, state = as.character(df[[kind$col]]),
               n = df$n, median = df$median, mean = df$mean,
               stringsAsFactors = FALSE)
  }))
  length_table <- do.call(rbind, lapply(STATE_KINDS, function(kind) {
    df <- length_stats[[kind$id]]
    data.frame(state_kind = kind$title, state = as.character(df[[kind$col]]),
               n = df$n, median = df$median, mean = df$mean,
               stringsAsFactors = FALSE)
  }))
  write_section_table(diff_table, file.path(out_dir, "mch_diff_by_state.tsv"))
  write_section_table(length_table,
                      file.path(out_dir, "gene_length_by_state.tsv"))
  write_section_table(key_genes,
                      file.path(out_dir, "key_genes_chromatin_state.tsv"))

  cat("\nCreating figures...\n")

  save_multiformat_ggplot(
    combine_kind_panels(
      lapply(STATE_KINDS, function(kind) {
        plot_distribution_panel(distribution[distribution$state_kind == kind$title, ],
                                kind)
      }),
      title = "Chromatin state of all tested genes",
      subtitle = "Promoter state over the TSS window, gene-body state over the gene body",
      widths = c(4, 5)),
    file.path(out_dir, "10_01a_state_distribution_all_genes"),
    width = 14, height = 7)

  save_multiformat_ggplot(
    (plot_direction_panel(distribution[distribution$state_kind == STATE_KIND_TITLES[1], ],
                          STATE_KINDS[[1]]) /
       plot_direction_panel(distribution[distribution$state_kind == STATE_KIND_TITLES[2], ],
                            STATE_KINDS[[2]])) +
      plot_annotation(
        title = "Chromatin state by mCH direction",
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                                size = 15))),
    file.path(out_dir, "10_01b_state_distribution_by_direction"),
    width = 15, height = 14)

  save_multiformat_ggplot(
    combine_kind_panels(
      lapply(STATE_KINDS, function(kind) {
        plot_significant_vs_all_panel(
          distribution[distribution$state_kind == kind$title, ], kind)
      }),
      title = "Chromatin state: significant mCH genes against all tested genes",
      widths = c(4, 5)),
    file.path(out_dir, "10_01c_state_significant_vs_all_tested"),
    width = 15, height = 7)

  save_multiformat_ggplot(
    combine_kind_panels(
      lapply(STATE_KINDS, function(kind) {
        plot_stacked_proportion_panel(
          distribution[distribution$state_kind == kind$title, ], kind)
      }),
      title = "Chromatin state composition of each gene set"),
    file.path(out_dir, "10_01d_state_stacked_proportion"),
    width = 16, height = 9)

  save_multiformat_ggplot(
    combine_kind_panels(
      lapply(STATE_KINDS, function(kind) {
        plot_significance_rate_panel(rate[rate$state_kind == kind$title, ], kind)
      }),
      title = "Share of genes with significant differential mCH, by chromatin state",
      subtitle = sprintf("FDR < %.2f. Labels give significant genes over tested genes.",
                         Q_THRESHOLD),
      widths = c(4, 5)),
    file.path(out_dir, "10_01e_state_significance_rate"),
    width = 15, height = 7)

  save_multiformat_ggplot(
    plot_state_crosstab(crosstab, nrow(genes)),
    file.path(out_dir, "10_01f_promoter_body_state_crosstab"),
    width = 12, height = 8)

  save_multiformat_ggplot(
    combine_kind_panels(
      lapply(REGIONS, function(region) {
        plot_mark_overlap_panel(mark_overlap, region)
      }),
      title = "Histone mark overlap over promoter windows and gene bodies"),
    file.path(out_dir, "10_01g_histone_mark_overlap_heatmap"),
    width = 18, height = 6)

  save_multiformat_ggplot(
    combine_kind_panels(
      lapply(STATE_KINDS, function(kind) {
        plot_state_enrichment_panel(enrichment, kind)
      }),
      title = "Enrichment of mCH-significant genes in each chromatin state",
      subtitle = "Fisher odds ratio against all other tested genes, 95% CI",
      widths = c(4, 5)),
    file.path(out_dir, "10_01h_state_enrichment_odds_ratio"),
    width = 15, height = 7)

  # Figure text for the two violin panels, added after every table is written.
  diff_stats <- lapply(diff_stats, function(stats) {
    stats$label <- group_label(stats, digits = 5)
    stats
  })
  length_stats <- lapply(length_stats, function(stats) {
    stats$label <- group_label(stats, digits = 2)
    stats
  })

  save_multiformat_ggplot(
    combine_kind_panels(
      lapply(STATE_KINDS, function(kind) {
        plot_value_by_state_panel(genes, "mch_diff", diff_stats[[kind$id]],
                                  kind, "mCH difference (mutant minus control)")
      }),
      title = "mCH difference across chromatin states",
      subtitle = sprintf("All %s tested genes", fmt_int(nrow(genes))),
      widths = c(4, 5)),
    file.path(out_dir, "10_01i_mch_diff_by_state"),
    width = 15, height = 7)

  save_multiformat_ggplot(
    combine_kind_panels(
      lapply(STATE_KINDS, function(kind) {
        plot_value_by_state_panel(genes, "log10_gene_length",
                                  length_stats[[kind$id]], kind,
                                  "log10 gene length (bp)")
      }),
      title = "Gene length across chromatin states",
      subtitle = sprintf("All %s tested genes", fmt_int(nrow(genes))),
      widths = c(4, 5)),
    file.path(out_dir, "10_01j_gene_length_by_state"),
    width = 15, height = 7)

  save_multiformat_ggplot(
    plot_key_genes(key_genes),
    file.path(out_dir, "10_01k_key_genes_chromatin_state"),
    width = 12, height = 7)

  save_multiformat_ggplot(
    combine_kind_panels(
      lapply(STATE_KINDS, function(kind) {
        plot_length_vs_logfc_panel(genes, key_genes, kind)
      }),
      title = "Gene length against differential mCH, by chromatin state"),
    file.path(out_dir, "10_01l_gene_length_vs_mch_logfc"),
    width = 18, height = 9)

  cat("\nSection 10_01 complete.\n")
  cat(sprintf("  Genes classified: %s\n", fmt_int(nrow(genes))))
  cat(sprintf("  Fisher tests registered: %s\n", fmt_int(nrow(enrichment))))
  cat(sprintf("  Output: %s\n\n", out_dir))
}

main()
