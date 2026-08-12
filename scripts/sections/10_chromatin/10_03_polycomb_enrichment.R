# scripts/sections/10_chromatin/10_03_polycomb_enrichment.R
#
# Section 10_03: Polycomb target enrichment in differential gene-body mCH.
#
# What this tests
#   Whether Polycomb-marked genes are depleted from mCH hypermethylation in the
#   BAP1 knockout. An inaccessible-heterochromatin model predicts depletion:
#   DNMT3A cannot reach Polycomb-compacted chromatin, so the genes that gain
#   mCH should be normally active genes that acquire ectopic H2AK119ub.
#
#   Five Polycomb definitions are tested, so the conclusion does not rest on
#   one gene set:
#     1. Polycomb gene-body state
#     2. Polycomb or Mixed gene-body state
#     3. Repressed_Promoter or Bivalent_Promoter promoter state
#     4. H3K27me3 consensus peak overlapping the gene body
#     5. Top decile of control H2AK119ub gene-body signal
#
#   Definitions 1 to 4 read the two state columns and the gene-body overlap
#   columns of the section 10_01 handoff table. Each definition is tested
#   against mCH hypermethylation and against mCH hypomethylation. Every
#   gene-body state and every promoter state is also tested against both
#   directions, and the two families are reported as two panels.
#
# Reads
#   HANDOFF_PATHS$chromatin_state       gene-level promoter state, gene-body
#                                       state and mark overlaps, from 10_01
#   DIFFBIND_PATHS$k119ub_gene_signal   per-gene H2AK119ub signal
#   mch_results                         loaded by _shared_config.R
#
# Writes (OUT_DIR defaults to results/sections/10_chromatin/)
#   Figures 10_03a .. 10_03h, each in a multi-format subdirectory
#   Tables prefixed 10_03_
#   Fisher gene tables under fisher_tables/ and rows in the shared registry,
#   which section 40_04 validates by permutation
#
# Adapted from Biomodal section 30 (section_30_polycomb_target_enrichment.R).
# EM-seq measures mCH as one modality, so the second-modality Fisher tests and
# the second-modality magnitude violin (30e) of that section have no analogue
# here and are not carried over.

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)

SECTION_ID <- "10_03"

# Fill colours for the two-group comparisons in this section.
POLYCOMB_GROUP_COLORS <- c(
  "Non-Polycomb" = "#999999",
  "Polycomb"     = "#4daf4a"
)

POLYCOMB_GROUP_LEVELS <- c("Non-Polycomb", "Polycomb")

MCH_STATUS_LEVELS <- c("Hypermethylated", "Hypomethylated", "Not significant")

MCH_STATUS_COLORS <- c(
  "Hypermethylated" = unname(COLORS$direction["Hypermethylated"]),
  "Hypomethylated"  = unname(COLORS$direction["Hypomethylated"]),
  "Not significant" = "grey80"
)

# Section 10_01 records two chromatin states per gene. Every per-state test,
# table and figure in this section is built once per family and the two are
# reported as two panels.
state_families <- function() {
  list(
    list(id     = "body",
         col    = "body_state",
         levels = BODY_STATE_ORDER,
         colors = BODY_STATE_COLORS,
         title  = "Gene-body state"),
    list(id     = "promoter",
         col    = "promoter_state",
         levels = PROMOTER_STATE_ORDER,
         colors = PROMOTER_STATE_COLORS,
         title  = "Promoter state")
  )
}

STATE_FAMILY_TITLES <- vapply(state_families(), function(f) f$title,
                              character(1))

# The promoter and body level names do not collide, so one palette and one
# level order serve the figures that draw both families.
STATE_LEVELS_ALL <- c(BODY_STATE_ORDER, PROMOTER_STATE_ORDER)
STATE_COLORS_ALL <- c(BODY_STATE_COLORS, PROMOTER_STATE_COLORS)

# =============================================================================
# COMMAND LINE
# =============================================================================

parse_section_options <- function() {
  option_list <- list(
    make_option("--output-dir", dest = "output_dir", type = "character",
                default = OUTPUT_PATHS$chromatin,
                help = "Directory for figures and tables [default: %default]"),
    make_option("--fdr-threshold", dest = "fdr_threshold", type = "double",
                default = Q_THRESHOLD,
                help = "edgeR FDR cutoff calling a gene differentially methylated [default: %default]"),
    make_option("--k119ub-top-fraction", dest = "k119ub_top_fraction",
                type = "double", default = 0.10,
                help = "Top fraction of control H2AK119ub gene-body signal defining the K119ub Polycomb set [default: %default]")
  )
  opt <- parse_args(OptionParser(option_list = option_list))

  if (opt$fdr_threshold <= 0 || opt$fdr_threshold >= 1) {
    stop("--fdr-threshold must be between 0 and 1, got ", opt$fdr_threshold)
  }
  if (opt$k119ub_top_fraction <= 0 || opt$k119ub_top_fraction >= 1) {
    stop("--k119ub-top-fraction must be between 0 and 1, got ",
         opt$k119ub_top_fraction)
  }
  opt
}

# =============================================================================
# SMALL FORMATTERS
# =============================================================================

fmt_p <- function(p) {
  if (length(p) != 1 || is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

sig_stars <- function(q) {
  vapply(q, function(x) {
    if (is.na(x)) return("")
    if (x < 0.001) return("***")
    if (x < 0.01) return("**")
    if (x < 0.05) return("*")
    "ns"
  }, character(1))
}

# =============================================================================
# INPUT LOADING
# =============================================================================

#' Read the gene-level chromatin state table written by section 10_01.
#'
#' Section 10_01 writes one row per row of mch_results, so gene_id is the key
#' that joins the two tables one to one. Only the columns this section reads
#' are kept: the two state columns and the gene-body H3K27me3 overlap.
load_chromatin_state_table <- function(path) {
  if (!file.exists(path)) {
    stop("Gene chromatin state table not found: ", path,
         "\nRun section 10_01 first; it writes this file.")
  }

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")

  required <- c("gene_name", "gene_id", "promoter_state", "body_state",
                "body_h3k27me3_overlap")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("Chromatin state table from section 10_01 is missing columns: ",
         paste(missing, collapse = ", "), "\nFile: ", path)
  }

  unknown_promoter <- setdiff(unique(df$promoter_state), PROMOTER_STATE_ORDER)
  if (length(unknown_promoter) > 0) {
    stop("Chromatin state table holds promoter states outside ",
         "PROMOTER_STATE_ORDER: ", paste(unknown_promoter, collapse = ", "))
  }
  unknown_body <- setdiff(unique(df$body_state), BODY_STATE_ORDER)
  if (length(unknown_body) > 0) {
    stop("Chromatin state table holds body states outside BODY_STATE_ORDER: ",
         paste(unknown_body, collapse = ", "))
  }

  df <- df[, required, drop = FALSE]
  df$body_h3k27me3_overlap <- as.logical(df$body_h3k27me3_overlap)
  if (any(is.na(df$body_h3k27me3_overlap))) {
    stop("body_h3k27me3_overlap has values that are neither TRUE nor FALSE in ",
         path)
  }

  if (any(duplicated(df$gene_id))) {
    stop("gene_id is not unique in the section 10_01 handoff table: ",
         sum(duplicated(df$gene_id)), " repeated identifiers in ", path)
  }

  cat(sprintf("  Chromatin state table: %s genes\n",
              format(nrow(df), big.mark = ",")))
  df
}

#' Read the per-gene H2AK119ub signal table and keep one row per symbol.
load_k119ub_gene_signal <- function(path) {
  if (!file.exists(path)) {
    stop("H2AK119ub gene signal table not found: ", path)
  }

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")

  required <- c("symbol", "gb_ctrl_signal", "gb_mut_signal", "gb_log2fc",
                "gb_signal_class")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("H2AK119ub gene signal table is missing columns: ",
         paste(missing, collapse = ", "), "\nFile: ", path)
  }

  df <- df[, required, drop = FALSE]
  df <- df[order(df$symbol, -df$gb_ctrl_signal), , drop = FALSE]

  n_before <- nrow(df)
  df <- df[!duplicated(df$symbol), , drop = FALSE]

  cat(sprintf("  H2AK119ub gene signal: %s rows, %s unique symbols\n",
              format(n_before, big.mark = ","),
              format(nrow(df), big.mark = ",")))
  df
}

#' Reduce mch_results to one row per gene name.
#'
#' Some gene names carry more than one ENSMUSG identifier. The row with the
#' largest absolute edgeR log fold change is kept; ties break on the smaller
#' FDR, then on gene_id.
deduplicate_mch_results <- function(mch, fdr_threshold) {
  keep_cols <- c("gene_name", "gene_id", "chr", "start", "end", "gene_length",
                 "mch_ctrl", "mch_mut", "mch_diff", "edger_logFC", "edger_fdr")
  df <- mch[, keep_cols, drop = FALSE]

  if (any(is.na(df$edger_fdr))) {
    stop("edger_fdr is NA for ", sum(is.na(df$edger_fdr)),
         " genes in mch_differential_results.tsv.")
  }
  if (any(is.na(df$mch_diff))) {
    stop("mch_diff is NA for ", sum(is.na(df$mch_diff)),
         " genes in mch_differential_results.tsv.")
  }

  df <- df[order(-abs(df$edger_logFC), df$edger_fdr, df$gene_id), , drop = FALSE]

  n_before <- nrow(df)
  df <- df[!duplicated(df$gene_name), , drop = FALSE]

  df$mch_sig   <- df$edger_fdr < fdr_threshold
  df$mch_hyper <- df$mch_sig & df$mch_diff > 0
  df$mch_hypo  <- df$mch_sig & df$mch_diff < 0

  cat(sprintf("  mCH results: %s rows, %s unique genes at FDR < %.3f (%s hyper, %s hypo)\n",
              format(n_before, big.mark = ","),
              format(nrow(df), big.mark = ","),
              fdr_threshold,
              format(sum(df$mch_hyper), big.mark = ","),
              format(sum(df$mch_hypo), big.mark = ",")))
  df
}

# =============================================================================
# ANALYSIS UNIVERSE AND POLYCOMB DEFINITIONS
# =============================================================================

#' Join the mCH results to the chromatin state table on gene_id.
#'
#' gene_id is unique on both sides, so the join is one to one and no gene name
#' with two identifiers has to be collapsed. gene_name comes from the mCH side.
build_analysis_universe <- function(mch_unique, chromatin_states) {
  state_cols <- chromatin_states[, c("gene_id", "promoter_state", "body_state",
                                     "body_h3k27me3_overlap"), drop = FALSE]
  genes <- dplyr::inner_join(mch_unique, state_cols, by = "gene_id")

  if (nrow(genes) == 0) {
    stop("No gene_id is shared between mch_differential_results.tsv and the ",
         "chromatin state table from section 10_01.")
  }

  genes$promoter_state <- factor(genes$promoter_state,
                                 levels = PROMOTER_STATE_ORDER)
  genes$body_state <- factor(genes$body_state, levels = BODY_STATE_ORDER)

  cat(sprintf("  Analysis universe: %s genes with both an mCH result and a chromatin state\n",
              format(nrow(genes), big.mark = ",")))
  cat(sprintf("    %s mCH genes have no chromatin state and are dropped\n",
              format(nrow(mch_unique) - nrow(genes), big.mark = ",")))
  for (family in state_families()) {
    cat(sprintf("  %s counts:\n", family$title))
    for (state in family$levels) {
      n_state <- sum(genes[[family$col]] == state)
      cat(sprintf("    %-20s %6s (%.1f%%)\n", state,
                  format(n_state, big.mark = ","),
                  100 * n_state / nrow(genes)))
    }
  }
  genes
}

#' Add the chromatin-state Polycomb flags read from the section 10_01 handoff.
add_state_polycomb_flags <- function(genes) {
  genes$pc_body_polycomb <- genes$body_state == "Polycomb"
  genes$pc_body_polycomb_mixed <- genes$body_state %in% c("Polycomb", "Mixed")
  genes$pc_repressed_bivalent <- genes$promoter_state %in%
    c("Repressed_Promoter", "Bivalent_Promoter")
  genes$pc_k27me3_peak <- genes$body_h3k27me3_overlap

  cat(sprintf("  Polycomb gene-body state: %s genes (%.1f%%)\n",
              format(sum(genes$pc_body_polycomb), big.mark = ","),
              100 * mean(genes$pc_body_polycomb)))
  cat(sprintf("  Polycomb or Mixed gene-body state: %s genes (%.1f%%)\n",
              format(sum(genes$pc_body_polycomb_mixed), big.mark = ","),
              100 * mean(genes$pc_body_polycomb_mixed)))
  cat(sprintf("  Repressed or Bivalent promoter state: %s genes (%.1f%%)\n",
              format(sum(genes$pc_repressed_bivalent), big.mark = ","),
              100 * mean(genes$pc_repressed_bivalent)))
  cat(sprintf("  Gene body overlapping an H3K27me3 peak: %s genes (%.1f%%)\n",
              format(sum(genes$pc_k27me3_peak), big.mark = ","),
              100 * mean(genes$pc_k27me3_peak)))
  genes
}

#' Attach control H2AK119ub gene-body signal and flag the top fraction.
#'
#' The flag uses control signal, so it marks genes that are Polycomb-marked
#' before BAP1 loss. Genes with no entry in the signal table keep NA, and
#' register_fisher_test() drops them from the tests that use this flag.
add_k119ub_top_flag <- function(genes, k119ub_signal, top_fraction) {
  genes <- dplyr::left_join(genes, k119ub_signal,
                            by = c("gene_name" = "symbol"))

  cutoff <- quantile(genes$gb_ctrl_signal, probs = 1 - top_fraction,
                     na.rm = TRUE, names = FALSE)
  genes$pc_k119ub_top <- genes$gb_ctrl_signal >= cutoff
  genes$k119ub_cutoff <- cutoff

  n_with_signal <- sum(!is.na(genes$gb_ctrl_signal))
  cat(sprintf("  H2AK119ub signal matched for %s of %s genes\n",
              format(n_with_signal, big.mark = ","),
              format(nrow(genes), big.mark = ",")))
  cat(sprintf("  Top %.0f%% control signal cutoff = %.4f, %s genes flagged\n",
              100 * top_fraction, cutoff,
              format(sum(genes$pc_k119ub_top, na.rm = TRUE), big.mark = ",")))
  genes
}

state_column <- function(family_id, state) paste0(family_id, "_state_", state)

add_state_indicator_columns <- function(genes) {
  for (family in state_families()) {
    for (state in family$levels) {
      genes[[state_column(family$id, state)]] <- genes[[family$col]] == state
    }
  }
  genes
}

#' The Polycomb definitions tested in this section.
polycomb_definitions <- function(top_fraction) {
  list(
    list(id = "body_polycomb", col = "pc_body_polycomb",
         label = "Polycomb body state"),
    list(id = "body_polycomb_mixed", col = "pc_body_polycomb_mixed",
         label = "Polycomb or Mixed body"),
    list(id = "repressed_bivalent", col = "pc_repressed_bivalent",
         label = "Repressed/Bivalent promoter"),
    list(id = "k27me3_peak", col = "pc_k27me3_peak",
         label = "H3K27me3 body overlap"),
    list(id = "k119ub_top", col = "pc_k119ub_top",
         label = sprintf("K119ub top %.0f%%", 100 * top_fraction))
  )
}

mch_directions <- function() {
  list(
    list(id = "hyper", col = "mch_hyper", label = "Hypermethylated"),
    list(id = "hypo",  col = "mch_hypo",  label = "Hypomethylated")
  )
}

# =============================================================================
# FISHER TESTS
# =============================================================================

#' Run one gene-level 2x2 Fisher test through the shared registry.
#'
#' Returns a one-row data.frame holding the 2x2 counts, the conditional odds
#' ratio with its confidence interval, and the sample odds ratio.
run_registered_fisher <- function(genes, set_col, dmr_col, test_id, description,
                                  set_label, direction_label, out_dir) {
  needed <- c("gene_name", "chr", set_col, dmr_col)
  missing <- setdiff(needed, colnames(genes))
  if (length(missing) > 0) {
    stop("run_registered_fisher(): genes is missing columns: ",
         paste(missing, collapse = ", "))
  }

  df <- genes[, needed, drop = FALSE]
  df <- df[!is.na(df[[set_col]]) & !is.na(df[[dmr_col]]), , drop = FALSE]

  ft <- register_fisher_test(
    section = SECTION_ID, test_id = test_id, description = description,
    gene_df = df, row_var = set_col, col_var = dmr_col, output_dir = out_dir)

  in_set <- df[[set_col]]
  is_dmr <- df[[dmr_col]]
  a <- sum(in_set & is_dmr)
  b <- sum(!in_set & is_dmr)
  c_count <- sum(in_set & !is_dmr)
  d <- sum(!in_set & !is_dmr)

  data.frame(
    section = SECTION_ID,
    test_id = test_id,
    gene_set = set_label,
    direction = direction_label,
    set_col = set_col,
    dmr_col = dmr_col,
    n_genes = nrow(df),
    n_in_set = a + c_count,
    n_dmr = a + b,
    a_set_dmr = a,
    b_notset_dmr = b,
    c_set_notdmr = c_count,
    d_notset_notdmr = d,
    pct_dmr_in_set = ifelse((a + c_count) > 0, 100 * a / (a + c_count), NA_real_),
    pct_dmr_out_set = ifelse((b + d) > 0, 100 * b / (b + d), NA_real_),
    odds_ratio = unname(ft$estimate),
    ci_lower = ft$conf.int[1],
    ci_upper = ft$conf.int[2],
    sample_odds_ratio = (a * d) / (b * c_count),
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

#' Every Polycomb definition against both mCH directions.
run_definition_tests <- function(genes, definitions, out_dir) {
  directions <- mch_directions()
  rows <- list()

  for (def in definitions) {
    for (dir_spec in directions) {
      description <- sprintf(
        "Are genes in the %s Polycomb set enriched among mCH %s genes?",
        def$label, tolower(dir_spec$label))
      rows[[length(rows) + 1]] <- run_registered_fisher(
        genes = genes,
        set_col = def$col,
        dmr_col = dir_spec$col,
        test_id = sprintf("%s_%s", def$id, dir_spec$id),
        description = description,
        set_label = def$label,
        direction_label = dir_spec$label,
        out_dir = out_dir)
    }
  }
  dplyr::bind_rows(rows)
}

#' Every gene-body state and every promoter state against both mCH directions.
run_state_tests <- function(genes, out_dir) {
  directions <- mch_directions()
  rows <- list()

  for (family in state_families()) {
    for (state in family$levels) {
      for (dir_spec in directions) {
        description <- sprintf(
          "Are genes with the %s %s enriched among mCH %s genes?",
          state, tolower(family$title), tolower(dir_spec$label))
        row <- run_registered_fisher(
          genes = genes,
          set_col = state_column(family$id, state),
          dmr_col = dir_spec$col,
          test_id = sprintf("%s_state_%s_%s", family$id, state, dir_spec$id),
          description = description,
          set_label = state,
          direction_label = dir_spec$label,
          out_dir = out_dir)
        row$state_family <- family$title
        rows[[length(rows) + 1]] <- row
      }
    }
  }
  dplyr::bind_rows(rows)
}

# =============================================================================
# OVERLAP, RATE, AND MAGNITUDE TABLES
# =============================================================================

#' Sizes of the Polycomb definitions and their pairwise Jaccard overlap.
build_definition_tables <- function(genes, definitions) {
  sizes <- dplyr::bind_rows(lapply(definitions, function(def) {
    flag <- genes[[def$col]]
    tested <- !is.na(flag)
    data.frame(
      gene_set = def$label,
      column = def$col,
      n_tested = sum(tested),
      n_in_set = sum(flag, na.rm = TRUE),
      pct_in_set = 100 * sum(flag, na.rm = TRUE) / sum(tested),
      n_missing = sum(!tested),
      stringsAsFactors = FALSE
    )
  }))

  pairs <- list()
  for (i in seq_along(definitions)) {
    for (j in seq_along(definitions)) {
      if (j <= i) next
      flag_i <- genes[[definitions[[i]]$col]]
      flag_j <- genes[[definitions[[j]]$col]]
      both_known <- !is.na(flag_i) & !is.na(flag_j)
      fi <- flag_i[both_known]
      fj <- flag_j[both_known]
      union_size <- sum(fi | fj)
      pairs[[length(pairs) + 1]] <- data.frame(
        set_a = definitions[[i]]$label,
        set_b = definitions[[j]]$label,
        n_compared = sum(both_known),
        n_a = sum(fi),
        n_b = sum(fj),
        n_intersection = sum(fi & fj),
        n_union = union_size,
        jaccard = ifelse(union_size > 0, sum(fi & fj) / union_size, NA_real_),
        stringsAsFactors = FALSE
      )
    }
  }

  list(sizes = sizes, overlap = dplyr::bind_rows(pairs))
}

#' Percent hyper, hypo, and any-significant inside and outside each definition.
build_significance_rate_table <- function(genes, definitions) {
  dplyr::bind_rows(lapply(definitions, function(def) {
    flag <- genes[[def$col]]
    tested <- !is.na(flag)
    sub <- genes[tested, , drop = FALSE]
    group <- ifelse(flag[tested], "Polycomb", "Non-Polycomb")

    dplyr::bind_rows(lapply(POLYCOMB_GROUP_LEVELS, function(g) {
      idx <- group == g
      n <- sum(idx)
      data.frame(
        gene_set = def$label,
        group = g,
        n_genes = n,
        n_hyper = sum(sub$mch_hyper[idx]),
        n_hypo = sum(sub$mch_hypo[idx]),
        n_sig = sum(sub$mch_sig[idx]),
        pct_hyper = 100 * sum(sub$mch_hyper[idx]) / n,
        pct_hypo = 100 * sum(sub$mch_hypo[idx]) / n,
        pct_sig = 100 * sum(sub$mch_sig[idx]) / n,
        stringsAsFactors = FALSE
      )
    }))
  }))
}

#' Percent hyper, hypo, and any-significant for each state of both families.
#'
#' Keeps every level of both orders, including states with no genes, and
#' returns state as character so it joins to the Fisher results.
build_state_rate_table <- function(genes) {
  rows <- lapply(state_families(), function(family) {
    out <- genes %>%
      dplyr::group_by(.data[[family$col]], .drop = FALSE) %>%
      dplyr::summarise(
        n_genes = dplyr::n(),
        n_hyper = sum(mch_hyper),
        n_hypo = sum(mch_hypo),
        n_sig = sum(mch_sig),
        pct_hyper = 100 * sum(mch_hyper) / dplyr::n(),
        pct_hypo = 100 * sum(mch_hypo) / dplyr::n(),
        pct_sig = 100 * sum(mch_sig) / dplyr::n(),
        median_mch_diff = median(mch_diff),
        median_abs_mch_diff = median(abs(mch_diff)),
        .groups = "drop"
      ) %>%
      as.data.frame()
    names(out)[1] <- "state"
    out$state <- as.character(out$state)
    out$state_family <- family$title
    out[, c("state_family", setdiff(colnames(out), "state_family"))]
  })
  do.call(rbind, rows)
}

#' Long frame of |mch_diff| for one definition across three gene sets.
build_magnitude_frame <- function(genes, def_col) {
  flag <- genes[[def_col]]
  keep <- !is.na(flag)
  base <- data.frame(
    gene_name = genes$gene_name[keep],
    abs_mch_diff = abs(genes$mch_diff[keep]),
    polycomb_group = ifelse(flag[keep], "Polycomb", "Non-Polycomb"),
    mch_hyper = genes$mch_hyper[keep],
    mch_hypo = genes$mch_hypo[keep],
    stringsAsFactors = FALSE
  )

  set_levels <- c("All tested genes", "Hypermethylated (significant)",
                  "Hypomethylated (significant)")
  out <- rbind(
    transform(base, gene_set = set_levels[1]),
    transform(base[base$mch_hyper, , drop = FALSE], gene_set = set_levels[2]),
    transform(base[base$mch_hypo, , drop = FALSE], gene_set = set_levels[3])
  )
  out$gene_set <- factor(out$gene_set, levels = set_levels)
  out$polycomb_group <- factor(out$polycomb_group, levels = POLYCOMB_GROUP_LEVELS)
  out
}

#' Long frame of |mch_diff| for every definition, all tested genes.
build_definition_magnitude_frame <- function(genes, definitions) {
  out <- dplyr::bind_rows(lapply(definitions, function(def) {
    flag <- genes[[def$col]]
    keep <- !is.na(flag)
    data.frame(
      definition = def$label,
      abs_mch_diff = abs(genes$mch_diff[keep]),
      polycomb_group = ifelse(flag[keep], "Polycomb", "Non-Polycomb"),
      stringsAsFactors = FALSE
    )
  }))
  out$definition <- factor(out$definition,
                           levels = vapply(definitions, function(d) d$label,
                                           character(1)))
  out$polycomb_group <- factor(out$polycomb_group, levels = POLYCOMB_GROUP_LEVELS)
  out
}

#' Wilcoxon rank-sum comparisons of |mch_diff|, Polycomb versus non-Polycomb.
run_magnitude_tests <- function(genes, definitions) {
  rows <- list()

  for (def in definitions) {
    mag <- build_magnitude_frame(genes, def$col)
    for (set_name in levels(mag$gene_set)) {
      sub <- mag[mag$gene_set == set_name, , drop = FALSE]
      poly_vals <- sub$abs_mch_diff[sub$polycomb_group == "Polycomb"]
      other_vals <- sub$abs_mch_diff[sub$polycomb_group == "Non-Polycomb"]

      if (length(poly_vals) == 0 || length(other_vals) == 0) {
        stop("No genes in one group for definition '", def$label,
             "' and gene set '", set_name, "'. Polycomb n = ",
             length(poly_vals), ", non-Polycomb n = ", length(other_vals))
      }

      wt <- wilcox.test(poly_vals, other_vals)
      rows[[length(rows) + 1]] <- data.frame(
        gene_set = def$label,
        gene_subset = set_name,
        polycomb_n = length(poly_vals),
        polycomb_median = median(poly_vals),
        polycomb_mean = mean(poly_vals),
        polycomb_iqr = IQR(poly_vals),
        non_polycomb_n = length(other_vals),
        non_polycomb_median = median(other_vals),
        non_polycomb_mean = mean(other_vals),
        non_polycomb_iqr = IQR(other_vals),
        wilcox_W = unname(wt$statistic),
        wilcox_p = wt$p.value,
        stringsAsFactors = FALSE
      )
    }
  }

  out <- dplyr::bind_rows(rows)
  out$wilcox_q <- p.adjust(out$wilcox_p, method = "BH")
  out
}

# =============================================================================
# PLOT HELPERS
# =============================================================================

#' Per-facet n and median labels for a violin or boxplot.
#'
#' Returns one row per facet and group with the label text and the y position
#' the label is drawn at. The result is figure text, so it feeds a geom_text
#' layer and no written table.
facet_group_annotations <- function(df, facet_col, group_col, value_col,
                                    digits = 5) {
  facet_levels <- levels(df[[facet_col]])
  rows <- lapply(facet_levels, function(f) {
    sub <- df[df[[facet_col]] == f, , drop = FALSE]
    if (nrow(sub) == 0) {
      cat(sprintf("  Facet '%s' holds no genes; it gets no annotation.\n", f))
      return(NULL)
    }
    stats <- summarise_groups(sub, group_col, value_col)
    stats$label <- group_label(stats, digits = digits)
    stats[[facet_col]] <- factor(f, levels = facet_levels)
    stats$y_pos <- max(sub[[value_col]], na.rm = TRUE) * 1.06
    stats
  })
  dplyr::bind_rows(rows)
}

#' Prepare Fisher results for a forest plot on a log odds-ratio axis.
#'
#' Tests with an odds ratio of zero or infinity cannot sit on a log axis and
#' are left out of the plot; the count is printed and the table keeps them.
#' Open confidence bounds are drawn clipped to the axis range and flagged in
#' the ci_open column.
prepare_forest_frame <- function(fisher_df) {
  usable <- is.finite(fisher_df$odds_ratio) & fisher_df$odds_ratio > 0
  if (any(!usable)) {
    cat(sprintf("  %d test(s) have an odds ratio of 0 or infinity and are left off the forest plot.\n",
                sum(!usable)))
  }

  d <- fisher_df[usable, , drop = FALSE]
  if (nrow(d) == 0) {
    stop("Every Fisher test has an odds ratio of 0 or infinity; no forest plot.")
  }

  finite_lower <- d$ci_lower[is.finite(d$ci_lower) & d$ci_lower > 0]
  finite_upper <- d$ci_upper[is.finite(d$ci_upper) & d$ci_upper > 0]
  x_min <- min(c(finite_lower, d$odds_ratio)) * 0.8
  x_max <- max(c(finite_upper, d$odds_ratio)) * 1.25

  d$ci_open <- !is.finite(d$ci_lower) | d$ci_lower <= 0 | !is.finite(d$ci_upper)
  d$ci_lower_plot <- pmax(d$ci_lower, x_min)
  d$ci_upper_plot <- pmin(d$ci_upper, x_max)
  d$label_x <- x_max * 1.5
  d$axis_max <- x_max * 2.6
  d$sig_label <- sig_stars(d$q_value)
  d
}

forest_caption <- function(forest_df) {
  n_open <- sum(forest_df$ci_open)
  if (n_open == 0) {
    return("Points are conditional odds ratios; bars are 95% confidence intervals. Stars use BH-adjusted p-values.")
  }
  sprintf(paste0("Points are conditional odds ratios; bars are 95%% confidence intervals. ",
                 "%d interval(s) are open and drawn clipped to the axis; exact bounds are in the table. ",
                 "Stars use BH-adjusted p-values."), n_open)
}

# =============================================================================
# FIGURES
# =============================================================================

plot_status_stacked_bar <- function(genes, def, fisher_df) {
  flag <- genes[[def$col]]
  keep <- !is.na(flag)
  sub <- genes[keep, , drop = FALSE]

  sub$polycomb_group <- factor(ifelse(flag[keep], "Polycomb", "Non-Polycomb"),
                               levels = POLYCOMB_GROUP_LEVELS)
  sub$mch_status <- factor(
    dplyr::case_when(sub$mch_hyper ~ "Hypermethylated",
                     sub$mch_hypo ~ "Hypomethylated",
                     TRUE ~ "Not significant"),
    levels = MCH_STATUS_LEVELS)

  bar_data <- sub %>%
    dplyr::count(polycomb_group, mch_status, .drop = FALSE) %>%
    dplyr::group_by(polycomb_group) %>%
    dplyr::mutate(pct = 100 * n / sum(n), total = sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(x_label = sprintf("%s\n(n = %s)",
                                    as.character(polycomb_group),
                                    format(total, big.mark = ",", trim = TRUE)))

  hyper_row <- fisher_df[fisher_df$gene_set == def$label &
                           fisher_df$direction == "Hypermethylated", ]
  hypo_row <- fisher_df[fisher_df$gene_set == def$label &
                          fisher_df$direction == "Hypomethylated", ]
  if (nrow(hyper_row) != 1 || nrow(hypo_row) != 1) {
    stop("Expected one hyper and one hypo Fisher row for gene set '",
         def$label, "', found ", nrow(hyper_row), " and ", nrow(hypo_row), ".")
  }

  ggplot(bar_data, aes(x = x_label, y = pct, fill = mch_status)) +
    geom_col(position = "stack", width = 0.65, colour = "black", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%s\n(%.1f%%)",
                                  format(n, big.mark = ",", trim = TRUE), pct)),
              position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(values = MCH_STATUS_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = "mCH status of Polycomb and non-Polycomb genes",
      subtitle = sprintf("%s | hyper OR = %.2f (%s) | hypo OR = %.2f (%s)",
                         def$label,
                         hyper_row$odds_ratio, fmt_p(hyper_row$p_value),
                         hypo_row$odds_ratio, fmt_p(hypo_row$p_value)),
      x = NULL, y = "Percent of genes", fill = "mCH status"
    ) +
    theme_emseq() +
    theme(legend.position = "top")
}

plot_definition_forest <- function(forest_df) {
  d <- forest_df
  d$row_label <- sprintf("%s  |  %s", d$gene_set, d$direction)
  d$row_label <- factor(d$row_label,
                        levels = d$row_label[order(d$odds_ratio)])

  ggplot(d, aes(x = odds_ratio, y = row_label, colour = direction)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_errorbar(aes(xmin = ci_lower_plot, xmax = ci_upper_plot),
                  orientation = "y", width = 0.25, linewidth = 0.7) +
    geom_point(size = 3.2) +
    geom_text(aes(x = label_x,
                  label = sprintf("OR %.2f  %s  (%s/%s)",
                                  odds_ratio, sig_label,
                                  format(a_set_dmr, big.mark = ",", trim = TRUE),
                                  format(n_in_set, big.mark = ",", trim = TRUE))),
              hjust = 0, size = 3, colour = "black") +
    expand_limits(x = unique(d$axis_max)) +
    scale_x_log10(labels = scales::label_number(accuracy = 0.01)) +
    scale_colour_manual(values = COLORS$direction) +
    labs(
      title = "Polycomb enrichment in differential mCH",
      subtitle = "Odds ratio per Polycomb definition and mCH direction; counts are (direction genes in set / genes in set)",
      x = "Odds ratio (log scale)", y = NULL, colour = "mCH direction",
      caption = forest_caption(d)
    ) +
    theme_emseq() +
    theme(legend.position = "bottom",
          plot.caption = element_text(hjust = 0, size = 8))
}

plot_state_forest <- function(forest_df) {
  d <- forest_df
  d$state <- factor(d$gene_set, levels = rev(STATE_LEVELS_ALL))
  d$state_family <- factor(d$state_family, levels = STATE_FAMILY_TITLES)

  ggplot(d, aes(x = odds_ratio, y = state, colour = direction)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_errorbar(aes(xmin = ci_lower_plot, xmax = ci_upper_plot),
                  orientation = "y", width = 0.25, linewidth = 0.7) +
    geom_point(size = 3.2) +
    geom_text(aes(x = label_x,
                  label = sprintf("%.2f %s  (%s/%s)", odds_ratio, sig_label,
                                  format(a_set_dmr, big.mark = ",", trim = TRUE),
                                  format(n_in_set, big.mark = ",", trim = TRUE))),
              hjust = 0, size = 3, colour = "black") +
    facet_grid(state_family ~ direction, scales = "free_y", space = "free_y") +
    expand_limits(x = unique(d$axis_max)) +
    scale_x_log10(labels = scales::label_number(accuracy = 0.01)) +
    scale_colour_manual(values = COLORS$direction) +
    labs(
      title = "Differential mCH enrichment per chromatin state",
      subtitle = paste("Gene-body states in the top row, promoter states in the",
                       "bottom row. Each state against every other state;",
                       "counts are (direction genes in state / genes in state)"),
      x = "Odds ratio (log scale)", y = NULL, colour = "mCH direction",
      caption = forest_caption(d)
    ) +
    theme_emseq() +
    theme(legend.position = "none",
          plot.caption = element_text(hjust = 0, size = 8))
}

plot_state_rates <- function(state_rates, genome_hyper_pct, genome_hypo_pct) {
  long <- state_rates %>%
    dplyr::select(state_family, state, n_genes, pct_hyper, pct_hypo,
                  n_hyper, n_hypo) %>%
    dplyr::mutate(
      state = factor(state, levels = STATE_LEVELS_ALL),
      state_family = factor(state_family, levels = STATE_FAMILY_TITLES)) %>%
    tidyr::pivot_longer(
      cols = c("pct_hyper", "pct_hypo"),
      names_to = "direction", values_to = "pct") %>%
    dplyr::mutate(
      direction = ifelse(direction == "pct_hyper",
                         "Hypermethylated", "Hypomethylated"),
      n_direction = ifelse(direction == "Hypermethylated", n_hyper, n_hypo))

  reference <- data.frame(
    direction = c("Hypermethylated", "Hypomethylated"),
    genome_pct = c(genome_hyper_pct, genome_hypo_pct),
    stringsAsFactors = FALSE)

  ggplot(long, aes(x = state, y = pct, fill = state)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
    geom_hline(data = reference, aes(yintercept = genome_pct),
               linetype = "dashed", colour = "grey30", linewidth = 0.6) +
    geom_text(data = reference,
              aes(x = 0.6, y = genome_pct, label = sprintf("all genes: %.1f%%", genome_pct)),
              inherit.aes = FALSE, hjust = 0, vjust = -0.6, size = 3,
              fontface = "italic") +
    geom_text(aes(label = sprintf("%.1f%%\nn = %s / %s", pct,
                                  format(n_direction, big.mark = ",", trim = TRUE),
                                  format(n_genes, big.mark = ",", trim = TRUE))),
              vjust = -0.2, size = 2.8) +
    facet_grid(state_family ~ direction, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = STATE_COLORS_ALL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
    labs(
      title = "Differential mCH rate per chromatin state",
      subtitle = paste("Percent of genes in each state that are significant in",
                       "each direction. Gene-body states in the top row,",
                       "promoter states in the bottom row."),
      x = NULL, y = "Percent of genes in state", fill = "Chromatin state"
    ) +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "none")
}

plot_magnitude_violin <- function(mag_df, def_label, magnitude_tests) {
  annot <- facet_group_annotations(mag_df, "gene_set", "polycomb_group",
                                   "abs_mch_diff")

  subtitle_rows <- magnitude_tests[magnitude_tests$gene_set == def_label, ]
  subtitle <- paste(sprintf("%s: %s",
                            subtitle_rows$gene_subset,
                            vapply(subtitle_rows$wilcox_p, fmt_p, character(1))),
                    collapse = " | ")

  ggplot(mag_df, aes(x = polycomb_group, y = abs_mch_diff,
                     fill = polycomb_group)) +
    geom_violin(alpha = 0.75, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.14, fill = "white", outlier.size = 0.3,
                 outlier.alpha = 0.3) +
    geom_text(data = annot, aes(x = polycomb_group, y = y_pos, label = label),
              inherit.aes = FALSE, size = 2.9, vjust = 0) +
    facet_wrap(~ gene_set, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = POLYCOMB_GROUP_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.22))) +
    labs(
      title = sprintf("mCH change magnitude: %s", def_label),
      subtitle = sprintf("Wilcoxon rank-sum, Polycomb vs non-Polycomb — %s",
                         subtitle),
      x = NULL, y = "|mCH difference| (mutant - control)"
    ) +
    theme_emseq() +
    theme(legend.position = "none")
}

plot_magnitude_by_definition <- function(def_mag_df, magnitude_tests) {
  annot <- facet_group_annotations(def_mag_df, "definition", "polycomb_group",
                                   "abs_mch_diff")

  all_gene_tests <- magnitude_tests[magnitude_tests$gene_subset == "All tested genes", ]
  subtitle <- paste(sprintf("%s: %s", all_gene_tests$gene_set,
                            vapply(all_gene_tests$wilcox_p, fmt_p, character(1))),
                    collapse = " | ")

  ggplot(def_mag_df, aes(x = polycomb_group, y = abs_mch_diff,
                         fill = polycomb_group)) +
    geom_violin(alpha = 0.75, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.14, fill = "white", outlier.size = 0.3,
                 outlier.alpha = 0.3) +
    geom_text(data = annot, aes(x = polycomb_group, y = y_pos, label = label),
              inherit.aes = FALSE, size = 2.7, vjust = 0) +
    facet_wrap(~ definition, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = POLYCOMB_GROUP_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.22))) +
    labs(
      title = "mCH change magnitude across all Polycomb definitions",
      subtitle = sprintf("All tested genes. Wilcoxon rank-sum — %s", subtitle),
      x = NULL, y = "|mCH difference| (mutant - control)"
    ) +
    theme_emseq() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 1))
}

plot_significance_rate_by_definition <- function(rate_table) {
  definition_order <- unique(rate_table$gene_set)

  long <- rate_table %>%
    tidyr::pivot_longer(cols = c("pct_hyper", "pct_hypo", "pct_sig"),
                        names_to = "measure", values_to = "pct") %>%
    dplyr::mutate(
      n_measure = dplyr::case_when(
        measure == "pct_hyper" ~ n_hyper,
        measure == "pct_hypo" ~ n_hypo,
        TRUE ~ n_sig),
      measure = factor(dplyr::recode(measure,
                                     pct_hyper = "Hypermethylated",
                                     pct_hypo = "Hypomethylated",
                                     pct_sig = "Any significant"),
                       levels = c("Hypermethylated", "Hypomethylated",
                                  "Any significant")),
      gene_set = factor(gene_set, levels = definition_order),
      group = factor(group, levels = POLYCOMB_GROUP_LEVELS))

  ggplot(long, aes(x = gene_set, y = pct, fill = group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72,
             colour = "black", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.1f%%\nn = %s / %s", pct,
                                  format(n_measure, big.mark = ",", trim = TRUE),
                                  format(n_genes, big.mark = ",", trim = TRUE))),
              position = position_dodge(width = 0.8), vjust = -0.15, size = 2.5) +
    facet_wrap(~ measure, ncol = 3) +
    scale_fill_manual(values = POLYCOMB_GROUP_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    labs(
      title = "Differential mCH rate inside and outside each Polycomb definition",
      subtitle = "Percent of genes in each group that are significant, with the counts behind each bar",
      x = NULL, y = "Percent of genes", fill = NULL
    ) +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          legend.position = "top")
}

build_composite_panel <- function(p_stacked, p_forest, p_state_or, p_rate_def) {
  (p_stacked | p_forest) / (p_state_or | p_rate_def) +
    plot_annotation(
      title = "Polycomb target enrichment in BAP1-KO differential mCH",
      subtitle = "Inaccessible-heterochromatin prediction: Polycomb genes depleted from mCH hypermethylation",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic")
      )
    )
}

# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

print_fisher_summary <- function(fisher_df, heading) {
  cat(sprintf("\n%s\n", heading))
  cat(sprintf("  %-30s %-18s %8s %18s %12s %12s %6s\n",
              "Gene set", "Direction", "OR", "95% CI", "p", "q (BH)", "Sig"))
  for (i in seq_len(nrow(fisher_df))) {
    row <- fisher_df[i, ]
    cat(sprintf("  %-30s %-18s %8.3f %8.3f-%-9.3f %12.3g %12.3g %6s\n",
                row$gene_set, row$direction, row$odds_ratio,
                row$ci_lower, row$ci_upper, row$p_value, row$q_value,
                row$sig_label))
  }
}

#' Restate the odds ratio of each definition as enriched, depleted, or neither.
print_direction_of_effect <- function(fisher_df) {
  cat("\n--- Direction of effect at q < 0.05 ---\n")
  for (i in seq_len(nrow(fisher_df))) {
    row <- fisher_df[i, ]
    verdict <- if (is.na(row$q_value) || row$q_value >= 0.05) {
      "no significant difference"
    } else if (row$odds_ratio < 1) {
      "depleted"
    } else {
      "enriched"
    }
    cat(sprintf("  %-30s %-18s OR = %.3f, q = %.3g -> %s\n",
                row$gene_set, row$direction, row$odds_ratio, row$q_value,
                verdict))
  }
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opt <- parse_section_options()
  OUT_DIR <- opt$output_dir
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 10_03: POLYCOMB TARGET ENRICHMENT IN DIFFERENTIAL mCH\n")
  cat("================================================================================\n")
  cat(sprintf("Output directory: %s\n", OUT_DIR))
  cat(sprintf("FDR threshold:    %.3f\n", opt$fdr_threshold))
  cat(sprintf("K119ub top fraction: %.3f\n", opt$k119ub_top_fraction))
  cat("\n")

  # --- Step 1: inputs --------------------------------------------------------
  cat("STEP 1: Loading inputs\n")
  chromatin_states <- load_chromatin_state_table(HANDOFF_PATHS$chromatin_state)
  k119ub_signal <- load_k119ub_gene_signal(DIFFBIND_PATHS$k119ub_gene_signal)
  mch_unique <- deduplicate_mch_results(mch_results, opt$fdr_threshold)

  # --- Step 2: universe and Polycomb definitions -----------------------------
  cat("\nSTEP 2: Building the analysis universe and Polycomb definitions\n")
  genes <- build_analysis_universe(mch_unique, chromatin_states)
  genes <- add_state_polycomb_flags(genes)
  genes <- add_k119ub_top_flag(genes, k119ub_signal, opt$k119ub_top_fraction)
  genes <- add_state_indicator_columns(genes)

  definitions <- polycomb_definitions(opt$k119ub_top_fraction)

  def_tables <- build_definition_tables(genes, definitions)
  write_section_table(def_tables$sizes, file.path(OUT_DIR, "10_03_polycomb_definition_sizes.tsv"))
  write_section_table(def_tables$overlap, file.path(OUT_DIR, "10_03_polycomb_definition_overlap.tsv"))

  # --- Step 3: Fisher tests --------------------------------------------------
  cat("\nSTEP 3: Fisher exact tests through the shared registry\n")
  definition_fisher <- run_definition_tests(genes, definitions, OUT_DIR)
  state_fisher <- run_state_tests(genes, OUT_DIR)

  fisher_all <- dplyr::bind_rows(definition_fisher, state_fisher)
  fisher_all$q_value <- p.adjust(fisher_all$p_value, method = "BH")
  fisher_all$sig_label <- sig_stars(fisher_all$q_value)
  fisher_all$test_family <- c(rep("polycomb_definition", nrow(definition_fisher)),
                              rep("chromatin_state", nrow(state_fisher)))

  write_section_table(fisher_all, file.path(OUT_DIR, "10_03_polycomb_fisher_tests.tsv"))

  definition_fisher <- fisher_all[fisher_all$test_family == "polycomb_definition", ]
  state_fisher <- fisher_all[fisher_all$test_family == "chromatin_state", ]

  # --- Step 4: rate and magnitude tables -------------------------------------
  cat("\nSTEP 4: Rate and magnitude tables\n")
  rate_table <- build_significance_rate_table(genes, definitions)
  write_section_table(rate_table, file.path(OUT_DIR, "10_03_significance_rate_by_definition.tsv"))

  state_rates <- build_state_rate_table(genes)
  state_rates <- dplyr::left_join(
    state_rates,
    state_fisher %>%
      dplyr::select(state_family, state = gene_set, direction, odds_ratio,
                    ci_lower, ci_upper, p_value, q_value) %>%
      tidyr::pivot_wider(names_from = direction,
                         values_from = c(odds_ratio, ci_lower, ci_upper,
                                         p_value, q_value)),
    by = c("state_family", "state"))
  write_section_table(state_rates, file.path(OUT_DIR, "10_03_chromatin_state_enrichment.tsv"))

  magnitude_tests <- run_magnitude_tests(genes, definitions)
  write_section_table(magnitude_tests, file.path(OUT_DIR, "10_03_polycomb_magnitude_wilcoxon.tsv"))

  classification_cols <- c("gene_name", "gene_id", "chr", "start", "end",
                           "gene_length", "promoter_state", "body_state",
                           "mch_ctrl", "mch_mut", "mch_diff",
                           "edger_logFC", "edger_fdr",
                           "mch_sig", "mch_hyper", "mch_hypo",
                           "pc_body_polycomb", "pc_body_polycomb_mixed",
                           "pc_repressed_bivalent", "pc_k27me3_peak",
                           "pc_k119ub_top",
                           "gb_ctrl_signal", "gb_mut_signal", "gb_log2fc",
                           "gb_signal_class", "k119ub_cutoff")
  write_section_table(genes[, classification_cols, drop = FALSE],
                  file.path(OUT_DIR, "10_03_polycomb_gene_classification.tsv"))

  # --- Step 5: figures -------------------------------------------------------
  # The widest gene-body definition (Polycomb or Mixed) carries the stacked bar
  # and the magnitude violin.
  cat("\nSTEP 5: Figures\n")
  broad_def <- definitions[[which(vapply(
    definitions, function(d) d$col == "pc_body_polycomb_mixed", logical(1)))]]

  p_stacked <- plot_status_stacked_bar(genes, broad_def, definition_fisher)
  save_multiformat_ggplot(p_stacked,
                          file.path(OUT_DIR, "10_03a_polycomb_mch_status_stacked_bar"),
                          width = 8, height = 7)

  definition_forest <- prepare_forest_frame(definition_fisher)
  p_forest <- plot_definition_forest(definition_forest)
  save_multiformat_ggplot(p_forest,
                          file.path(OUT_DIR, "10_03b_polycomb_fisher_forest"),
                          width = 12, height = 8)

  state_forest <- prepare_forest_frame(state_fisher)
  p_state_or <- plot_state_forest(state_forest)
  save_multiformat_ggplot(p_state_or,
                          file.path(OUT_DIR, "10_03c_chromatin_state_odds_ratios"),
                          width = 13, height = 11)

  genome_hyper_pct <- 100 * sum(genes$mch_hyper) / nrow(genes)
  genome_hypo_pct <- 100 * sum(genes$mch_hypo) / nrow(genes)
  p_state_rates <- plot_state_rates(state_rates, genome_hyper_pct, genome_hypo_pct)
  save_multiformat_ggplot(p_state_rates,
                          file.path(OUT_DIR, "10_03d_chromatin_state_significance_rates"),
                          width = 14, height = 12)

  broad_magnitude <- build_magnitude_frame(genes, broad_def$col)
  p_violin <- plot_magnitude_violin(broad_magnitude, broad_def$label,
                                    magnitude_tests)
  save_multiformat_ggplot(p_violin,
                          file.path(OUT_DIR, "10_03e_mch_magnitude_violin"),
                          width = 13, height = 7)

  definition_magnitude <- build_definition_magnitude_frame(genes, definitions)
  p_violin_defs <- plot_magnitude_by_definition(definition_magnitude,
                                                magnitude_tests)
  save_multiformat_ggplot(p_violin_defs,
                          file.path(OUT_DIR, "10_03f_mch_magnitude_by_definition"),
                          width = 16, height = 7)

  p_rate_def <- plot_significance_rate_by_definition(rate_table)
  save_multiformat_ggplot(p_rate_def,
                          file.path(OUT_DIR, "10_03g_significance_rate_by_definition"),
                          width = 15, height = 7)

  p_composite <- build_composite_panel(p_stacked, p_forest, p_state_or, p_rate_def)
  save_multiformat_ggplot(p_composite,
                          file.path(OUT_DIR, "10_03h_composite_summary"),
                          width = 22, height = 17)

  # --- Step 6: console summary -----------------------------------------------
  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 10_03 SUMMARY\n")
  cat("================================================================================\n")
  cat(sprintf("Universe: %s genes with an mCH result and a chromatin state\n",
              format(nrow(genes), big.mark = ",")))
  cat(sprintf("mCH hypermethylated: %s (%.2f%%) | hypomethylated: %s (%.2f%%)\n",
              format(sum(genes$mch_hyper), big.mark = ","), genome_hyper_pct,
              format(sum(genes$mch_hypo), big.mark = ","), genome_hypo_pct))

  cat("\n--- Polycomb definition sizes ---\n")
  for (i in seq_len(nrow(def_tables$sizes))) {
    row <- def_tables$sizes[i, ]
    cat(sprintf("  %-30s %7s of %7s genes (%.1f%%), %s untested\n",
                row$gene_set,
                format(row$n_in_set, big.mark = ","),
                format(row$n_tested, big.mark = ","),
                row$pct_in_set,
                format(row$n_missing, big.mark = ",")))
  }

  print_fisher_summary(definition_fisher, "--- Polycomb definition tests ---")
  for (family in state_families()) {
    print_fisher_summary(
      state_fisher[state_fisher$state_family == family$title, ],
      sprintf("--- %s tests ---", family$title))
  }
  print_direction_of_effect(definition_fisher)

  cat("\n--- |mCH difference| medians, Polycomb vs non-Polycomb ---\n")
  for (i in seq_len(nrow(magnitude_tests))) {
    row <- magnitude_tests[i, ]
    cat(sprintf("  %-30s %-30s poly %.5f (n=%s) vs non-poly %.5f (n=%s), %s, q = %.3g\n",
                row$gene_set, row$gene_subset,
                row$polycomb_median, format(row$polycomb_n, big.mark = ","),
                row$non_polycomb_median, format(row$non_polycomb_n, big.mark = ","),
                fmt_p(row$wilcox_p), row$wilcox_q))
  }

  cat(sprintf("\nFigures and tables written to: %s\n", OUT_DIR))
  cat(sprintf("Fisher tests registered in: %s\n", HANDOFF_PATHS$fisher_registry))
  cat("Section 10_03 complete.\n")
}

main()
