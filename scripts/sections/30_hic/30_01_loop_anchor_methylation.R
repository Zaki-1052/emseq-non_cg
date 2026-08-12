# scripts/sections/30_hic/30_01_loop_anchor_methylation.R
#
# Section 30_01: gene-body mCH direction at differential Hi-C loop anchors.
#
# The question: does the direction of the mCH change track the direction of the
# Hi-C loop change after BAP1 loss? Genes are assigned to the anchors of gained
# and lost loops, then compared on hypermethylation rate, on the distribution of
# mch_diff, on convergence with H2AK119ub, and inside a logistic model of
# hypermethylation on anchor class plus K119ub gene-body fold change.
#
# Anchor-to-gene assignment runs two ways and both are reported:
#   direct   an anchor overlaps the gene body of a tested gene
#   great    an anchor overlaps a GREAT-style regulatory domain built from TxDb
#            (5 kb upstream, 1 kb downstream, extended to at most 100 kb toward
#            the neighbouring gene's basal domain)
#
# Reads:
#   HIC_PATHS$loops                    differential loop table
#   DIFFBIND_PATHS$k119ub_gene_signal  K119ub gene-body signal per gene symbol
#   k119ub_diffbind, mch_results, gene_bodies   pre-loaded by _shared_config.R
#   TxDb.Mmusculus.UCSC.mm10.knownGene gene models for the GREAT domains
#
# Writes: TSV tables and multi-format figures under OUTPUT_PATHS$hic
#   (results/sections/30_hic/), plus Fisher gene tables under fisher_tables/.
#
# Adapted from Biomodal section 27 (methylation x Hi-C loop anchor integration).
# Parts of section 27 need two separate methylation modalities, which EM-seq does
# not provide. This port keeps its loop-direction stratification (27b), its
# K119ub convergence test (27c), and its logistic model, all recast on mCH.

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)

# =============================================================================
# CONSTANTS
# =============================================================================

SECTION_ID <- "30_01"

# GREAT-style regulatory domain geometry, in base pairs.
GREAT_UPSTREAM      <- 5000
GREAT_DOWNSTREAM    <- 1000
GREAT_MAX_EXTENSION <- 100000

# Histone overlap columns carried per anchor in the loop table.
ANCHOR_HISTONE_MARKS <- c("H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3",
                          "Bivalent_Promoter")

# Loop direction coding used by the upstream Hi-C pipeline.
LOOP_DIRECTION_LABELS <- c(up_in_mutant = "Gained", down_in_mutant = "Lost")

# A gene can sit at anchors of both gained and lost loops, so the anchor class
# carries an explicit "Both anchors" level rather than forcing a priority rule.
ANCHOR_CLASS_ORDER <- c("Gained anchor", "Lost anchor", "Both anchors",
                        "No anchor")

ANCHOR_CLASS_COLORS <- c(
  "Gained anchor" = "#4575B4",
  "Lost anchor"   = "#D73027",
  "Both anchors"  = "#8C510A",
  "No anchor"     = "grey70"
)

# Association method key -> label used in tables, facets, and figure captions.
ASSOCIATION_METHODS <- c(great = "GREAT domain", direct = "Direct overlap")

MCH_DIGITS    <- 5
K119UB_DIGITS <- 3

# Separator used to build one grouping key from several columns. No group label
# in this section contains it, so the key splits back into its parts unchanged.
GROUP_KEY_SEP <- "@@"

# =============================================================================
# SMALL UTILITIES
# =============================================================================

#' Format one p-value for a figure subtitle.
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

#' Format a vector of p-values.
fmt_p_vec <- function(p) vapply(p, fmt_p, character(1))

#' Summarise a value column for every combination of several group columns.
#'
#' Builds one key from the group columns, calls the shared summarise_groups(),
#' then splits the key back into the original columns. The result holds data
#' only. The "n = N\nmed = M" text drawn on a figure comes from group_label()
#' called at the plot function.
annotate_group_stats <- function(df, group_cols, value_col) {
  key_parts <- lapply(group_cols, function(col) as.character(df[[col]]))
  df[[".group_key"]] <- do.call(paste, c(key_parts, list(sep = GROUP_KEY_SEP)))

  summ <- summarise_groups(df, ".group_key", value_col)
  parts <- do.call(rbind, strsplit(summ[[".group_key"]], GROUP_KEY_SEP,
                                   fixed = TRUE))
  colnames(parts) <- group_cols

  cbind(as.data.frame(parts, stringsAsFactors = FALSE),
        summ[, c("n", "median", "mean", "q25", "q75")])
}

#' Pairwise two-sided Wilcoxon rank-sum tests across the anchor classes present.
pairwise_wilcoxon <- function(df, group_col, value_col) {
  df <- df[!is.na(df[[value_col]]), , drop = FALSE]
  groups <- levels(droplevels(factor(df[[group_col]], levels = ANCHOR_CLASS_ORDER)))

  if (length(groups) < 2) {
    stop("pairwise_wilcoxon(): fewer than two groups present for ", group_col)
  }

  combos <- combn(groups, 2, simplify = FALSE)
  rows <- lapply(combos, function(pair) {
    x <- df[[value_col]][df[[group_col]] == pair[1]]
    y <- df[[value_col]][df[[group_col]] == pair[2]]
    wt <- wilcox.test(x, y)
    data.frame(
      group1 = pair[1], group2 = pair[2],
      n1 = length(x), n2 = length(y),
      median1 = median(x), median2 = median(y),
      median_difference = median(x) - median(y),
      W = unname(wt$statistic), p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out$p_adjusted <- p.adjust(out$p_value, method = "BH")
  out
}

# =============================================================================
# INPUT LOADING
# =============================================================================

#' Column names the loop table must carry.
#'
#' The names come from the header of extended_characterized_loops.tsv. Anchor
#' coordinates are chr1/start1/end1 and chr2/start2/end2, not anchor-prefixed.
#' The anchor annotation columns carry the "_extended" and "_ext" suffixes, and
#' the resolution column is resolution_kb.
required_loop_columns <- function() {
  base <- c("loop_id", "chr1", "start1", "end1", "chr2", "start2", "end2",
            "logFC", "FDR", "significant", "category", "direction",
            "resolution_kb", "loop_distance", "loop_type")
  per_anchor <- unlist(lapply(c(1L, 2L), function(i) {
    c(paste0("anchor", i, c("_type_extended", "_distance_to_tss_ext")),
      paste0("anchor", i, "_", ANCHOR_HISTONE_MARKS, "_overlap"))
  }))
  c(base, per_anchor)
}

#' Read the differential loop table and keep the significant loops.
read_significant_loops <- function(path) {
  if (!file.exists(path)) stop("Hi-C loop table not found: ", path)

  loops <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                      quote = "", comment.char = "")

  missing <- setdiff(required_loop_columns(), colnames(loops))
  if (length(missing) > 0) {
    stop("Loop table ", path, " is missing columns: ",
         paste(missing, collapse = ", "))
  }

  cat(sprintf("  Loop table: %d rows\n", nrow(loops)))

  loops$significant <- as.logical(loops$significant)
  loops <- loops[!is.na(loops$significant) & loops$significant, , drop = FALSE]
  if (nrow(loops) == 0) stop("No significant loops in ", path)

  unknown <- setdiff(unique(loops$direction), names(LOOP_DIRECTION_LABELS))
  if (length(unknown) > 0) {
    stop("Unexpected loop direction values in ", path, ": ",
         paste(unknown, collapse = ", "))
  }

  # loop_class is read off the direction column, so direction must agree with
  # the sign of logFC. A change of convention upstream would otherwise invert
  # every gained and lost label without any error.
  expected_gained <- unname(LOOP_DIRECTION_LABELS[loops$direction]) == "Gained"
  disagree <- expected_gained != (loops$logFC > 0)
  if (any(disagree)) {
    stop(sum(disagree), " of ", nrow(loops), " significant loops in ", path,
         " have a direction that disagrees with the sign of logFC.",
         "\nFirst offending loop_id: ", loops$loop_id[which(disagree)[1]],
         " (direction=", loops$direction[which(disagree)[1]],
         ", logFC=", loops$logFC[which(disagree)[1]], ")")
  }

  loops$loop_class <- unname(LOOP_DIRECTION_LABELS[loops$direction])
  cat(sprintf("  Significant loops: %d (%d gained, %d lost)\n",
              nrow(loops),
              sum(loops$loop_class == "Gained"),
              sum(loops$loop_class == "Lost")))
  loops
}

#' Read the K119ub gene-body signal table, one row per gene symbol.
read_k119ub_gene_signal <- function(path) {
  if (!file.exists(path)) stop("K119ub gene signal table not found: ", path)

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "", comment.char = "")

  required <- c("symbol", "gb_ctrl_signal", "gb_mut_signal", "gb_log2fc",
                "gb_signal_class")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("K119ub gene signal table ", path, " is missing columns: ",
         paste(missing, collapse = ", "))
  }

  # One row per symbol, preferring a row that carries a fold change.
  df <- df[order(df$symbol, is.na(df$gb_log2fc)), , drop = FALSE]
  df <- df[!duplicated(df$symbol), , drop = FALSE]

  out <- data.frame(
    gene_name           = df$symbol,
    k119ub_gb_ctrl      = df$gb_ctrl_signal,
    k119ub_gb_mut       = df$gb_mut_signal,
    k119ub_gb_log2fc    = df$gb_log2fc,
    k119ub_signal_class = df$gb_signal_class,
    stringsAsFactors = FALSE
  )

  cat(sprintf("  K119ub gene signal: %d genes, %d with a gene-body log2FC\n",
              nrow(out), sum(!is.na(out$k119ub_gb_log2fc))))
  out
}

# =============================================================================
# ANCHOR CONSTRUCTION
# =============================================================================

#' Reshape the loop table to one row per anchor (two rows per loop).
#'
#' Anchor i takes its coordinates from chr<i>/start<i>/end<i> and its annotation
#' from the anchor<i>_ columns. The start and end are copied unchanged here and
#' shifted to 1-based coordinates in anchor_table_to_granges().
build_anchor_table <- function(loops) {
  per_anchor <- lapply(c(1L, 2L), function(i) {
    prefix <- paste0("anchor", i, "_")

    df <- data.frame(
      anchor_id              = paste0(loops$loop_id, "_A", i),
      loop_id                = loops$loop_id,
      anchor_number          = i,
      chr                    = as.character(loops[[paste0("chr", i)]]),
      start                  = as.numeric(loops[[paste0("start", i)]]),
      end                    = as.numeric(loops[[paste0("end", i)]]),
      anchor_type            = as.character(loops[[paste0(prefix, "type_extended")]]),
      anchor_distance_to_tss = as.numeric(loops[[paste0(prefix, "distance_to_tss_ext")]]),
      loop_class             = loops$loop_class,
      loop_direction         = loops$direction,
      loop_logFC             = as.numeric(loops$logFC),
      loop_fdr               = as.numeric(loops$FDR),
      loop_category          = as.character(loops$category),
      loop_resolution_kb     = as.numeric(loops$resolution_kb),
      loop_distance          = as.numeric(loops$loop_distance),
      loop_type              = as.character(loops$loop_type),
      stringsAsFactors = FALSE
    )

    for (mark in ANCHOR_HISTONE_MARKS) {
      df[[paste0(mark, "_overlap")]] <-
        as.logical(loops[[paste0(prefix, mark, "_overlap")]])
    }
    df
  })

  anchors <- do.call(rbind, per_anchor)
  cat(sprintf("  Anchors: %d (%d gained, %d lost)\n",
              nrow(anchors),
              sum(anchors$loop_class == "Gained"),
              sum(anchors$loop_class == "Lost")))
  anchors
}

#' Convert the anchor table to GRanges, keeping the row order of the table.
#'
#' Anchor coordinates in the loop table are BED-style half-open intervals with
#' 0-based starts: end minus start equals the anchor width. The start is shifted
#' by one to match the 1-based ranges used everywhere else in the pipeline.
anchor_table_to_granges <- function(anchors) {
  GRanges(
    seqnames = anchors$chr,
    ranges = IRanges(start = anchors$start + 1L, end = anchors$end),
    anchor_id = anchors$anchor_id,
    loop_class = anchors$loop_class
  )
}

# =============================================================================
# GREAT-STYLE REGULATORY DOMAINS
# =============================================================================

#' Build GREAT-style regulatory domains from the mm10 TxDb gene models.
#'
#' Each gene gets a basal domain of GREAT_UPSTREAM bp upstream and
#' GREAT_DOWNSTREAM bp downstream of its TSS, strand-aware. The domain then
#' extends in both directions up to GREAT_MAX_EXTENSION bp, stopping at the
#' basal domain of the neighbouring gene on the same chromosome.
build_regulatory_domains <- function() {
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  genes_gr <- genes(txdb)
  genes_gr <- genes_gr[as.character(seqnames(genes_gr)) %in% CANONICAL_CHRS]

  tss_pos <- ifelse(as.character(strand(genes_gr)) == "-",
                    end(genes_gr), start(genes_gr))

  gene_info <- data.frame(
    entrez_id = names(genes_gr),
    chr = as.character(seqnames(genes_gr)),
    tss = tss_pos,
    strand = as.character(strand(genes_gr)),
    stringsAsFactors = FALSE
  )

  gene_info %>%
    dplyr::arrange(chr, tss) %>%
    dplyr::mutate(
      basal_start = ifelse(strand == "+", tss - GREAT_UPSTREAM, tss - GREAT_DOWNSTREAM),
      basal_end   = ifelse(strand == "+", tss + GREAT_DOWNSTREAM, tss + GREAT_UPSTREAM),
      max_start   = tss - GREAT_MAX_EXTENSION,
      max_end     = tss + GREAT_MAX_EXTENSION
    ) %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(
      prev_basal_end   = dplyr::lag(basal_end, default = -Inf),
      next_basal_start = dplyr::lead(basal_start, default = Inf),
      reg_start = pmax(max_start, prev_basal_end, 1),
      reg_end   = pmin(max_end, next_basal_start)
    ) %>%
    dplyr::ungroup() %>%
    # needs_basal is computed once, before either replacement. mutate() evaluates
    # its arguments in order, so testing reg_end < reg_start a second time would
    # read the already replaced reg_start and the two halves would disagree.
    dplyr::mutate(needs_basal = reg_end < reg_start) %>%
    dplyr::mutate(
      reg_start = ifelse(needs_basal, basal_start, reg_start),
      reg_end   = ifelse(needs_basal, basal_end, reg_end),
      reg_end   = pmax(reg_end, reg_start)
    ) %>%
    dplyr::select(entrez_id, chr, tss, strand, reg_start, reg_end) %>%
    as.data.frame()
}

#' Map Entrez gene IDs to gene symbols, one symbol per ID.
map_entrez_to_symbol <- function(entrez_ids) {
  keys <- unique(entrez_ids[!is.na(entrez_ids) & nzchar(entrez_ids)])
  mapping <- AnnotationDbi::select(org.Mm.eg.db, keys = keys,
                                   columns = "SYMBOL", keytype = "ENTREZID")
  mapping <- mapping[!is.na(mapping$SYMBOL), , drop = FALSE]
  mapping <- mapping[!duplicated(mapping$ENTREZID), , drop = FALSE]
  setNames(mapping$SYMBOL, mapping$ENTREZID)
}

#' GREAT regulatory domains as GRanges carrying a gene_name column.
build_great_domain_granges <- function() {
  domains <- build_regulatory_domains()
  cat(sprintf("  TxDb regulatory domains: %d genes (%d kb up, %d kb down, %d kb max extension)\n",
              nrow(domains), GREAT_UPSTREAM / 1000, GREAT_DOWNSTREAM / 1000,
              GREAT_MAX_EXTENSION / 1000))

  symbol_of <- map_entrez_to_symbol(domains$entrez_id)
  domains$gene_name <- unname(symbol_of[domains$entrez_id])
  domains <- domains[!is.na(domains$gene_name), , drop = FALSE]
  cat(sprintf("  Domains with a gene symbol: %d (%d unique symbols)\n",
              nrow(domains), length(unique(domains$gene_name))))

  GRanges(
    seqnames = domains$chr,
    ranges = IRanges(start = domains$reg_start, end = domains$reg_end),
    gene_name = domains$gene_name,
    entrez_id = domains$entrez_id
  )
}

# =============================================================================
# ANCHOR-TO-GENE ASSOCIATION
# =============================================================================

#' Overlap anchors with a target GRanges that carries a gene_name column.
#'
#' @return data.frame with anchor_idx (row index into the anchor table) and
#'   gene_name, one row per overlapping pair.
associate_anchors_to_genes <- function(anchor_gr, target_gr, method_label) {
  hits <- findOverlaps(anchor_gr, target_gr, ignore.strand = TRUE)
  if (length(hits) == 0) {
    stop("No anchor-to-gene overlaps for method '", method_label,
         "'. Check that anchor and gene coordinates use the same genome build.")
  }

  assoc <- data.frame(
    anchor_idx = queryHits(hits),
    gene_name = as.character(target_gr$gene_name[subjectHits(hits)]),
    stringsAsFactors = FALSE
  )
  assoc <- assoc[!is.na(assoc$gene_name) & nzchar(assoc$gene_name), , drop = FALSE]

  cat(sprintf("  %s: %d anchor-gene pairs, %d unique genes\n",
              method_label, nrow(assoc), length(unique(assoc$gene_name))))
  assoc
}

# =============================================================================
# GENE-LEVEL TABLE
# =============================================================================

#' One row per gene symbol from the mCH results.
#'
#' A few gene symbols carry more than one ENSMUSG entry. The row with the
#' largest absolute edgeR logFC is kept, so the gene-level Fisher tables hold
#' exactly one row per gene.
build_gene_universe <- function(mch) {
  deduped <- mch %>%
    dplyr::arrange(gene_name, dplyr::desc(abs(edger_logFC))) %>%
    dplyr::distinct(gene_name, .keep_all = TRUE) %>%
    dplyr::select(gene_name, gene_id, chr, start, end, gene_length,
                  mch_ctrl, mch_mut, mch_diff, edger_logFC, edger_fdr,
                  mch_sig, mch_hyper, mch_hypo, mch_direction) %>%
    as.data.frame()

  cat(sprintf("  mCH universe: %d genes after collapsing %d duplicate symbols\n",
              nrow(deduped), nrow(mch) - nrow(deduped)))
  deduped
}

#' Collapse anchor-gene pairs to one row per gene with anchor and K119ub flags.
summarise_anchor_membership <- function(assoc, anchors, k119ub_gained_at_anchor,
                                        k119ub_lost_at_anchor) {
  assoc$loop_class <- anchors$loop_class[assoc$anchor_idx]
  assoc$k119ub_gained <- k119ub_gained_at_anchor[assoc$anchor_idx]
  assoc$k119ub_lost <- k119ub_lost_at_anchor[assoc$anchor_idx]

  assoc %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(
      n_anchors               = dplyr::n(),
      n_gained_anchors        = sum(loop_class == "Gained"),
      n_lost_anchors          = sum(loop_class == "Lost"),
      at_gained_anchor        = any(loop_class == "Gained"),
      at_lost_anchor          = any(loop_class == "Lost"),
      k119ub_gained_at_anchor = any(k119ub_gained),
      k119ub_lost_at_anchor   = any(k119ub_lost),
      .groups = "drop"
    ) %>%
    as.data.frame()
}

#' Join anchor membership and K119ub signal onto the gene universe.
#'
#' Genes with no anchor overlap get FALSE flags and zero anchor counts, so the
#' returned table covers the whole universe of the association method.
build_gene_anchor_table <- function(universe, membership, k119ub_signal,
                                    method_key) {
  tbl <- universe %>%
    dplyr::left_join(membership, by = "gene_name") %>%
    dplyr::left_join(k119ub_signal, by = "gene_name")

  for (col in c("n_anchors", "n_gained_anchors", "n_lost_anchors")) {
    tbl[[col]][is.na(tbl[[col]])] <- 0
  }
  for (col in c("at_gained_anchor", "at_lost_anchor",
                "k119ub_gained_at_anchor", "k119ub_lost_at_anchor")) {
    tbl[[col]][is.na(tbl[[col]])] <- FALSE
  }

  tbl$anchor_class <- factor(
    ifelse(tbl$at_gained_anchor & tbl$at_lost_anchor, "Both anchors",
      ifelse(tbl$at_gained_anchor, "Gained anchor",
        ifelse(tbl$at_lost_anchor, "Lost anchor", "No anchor"))),
    levels = ANCHOR_CLASS_ORDER
  )

  tbl$at_any_anchor <- tbl$anchor_class != "No anchor"
  tbl$method_key <- method_key
  tbl$method <- unname(ASSOCIATION_METHODS[method_key])

  cat(sprintf("  %s gene table: %d genes (%d gained, %d lost, %d both, %d none)\n",
              unname(ASSOCIATION_METHODS[method_key]), nrow(tbl),
              sum(tbl$anchor_class == "Gained anchor"),
              sum(tbl$anchor_class == "Lost anchor"),
              sum(tbl$anchor_class == "Both anchors"),
              sum(tbl$anchor_class == "No anchor")))
  tbl
}

# =============================================================================
# ANALYSES
# =============================================================================

#' Hypermethylation and hypomethylation rates for every anchor class.
compute_class_rates <- function(gene_tbl) {
  gene_tbl %>%
    dplyr::group_by(method, anchor_class) %>%
    dplyr::summarise(
      n_genes         = dplyr::n(),
      n_sig           = sum(mch_sig),
      n_hyper         = sum(mch_hyper),
      n_hypo          = sum(mch_hypo),
      pct_sig         = 100 * mean(mch_sig),
      pct_hyper       = 100 * mean(mch_hyper),
      pct_hypo        = 100 * mean(mch_hypo),
      median_mch_diff = median(mch_diff, na.rm = TRUE),
      mean_mch_diff   = mean(mch_diff, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.frame()
}

#' Run and register the three gene-level Fisher tests for one method.
run_registered_fisher_tests <- function(gene_tbl, method_key, out_dir) {
  method_label <- unname(ASSOCIATION_METHODS[method_key])

  specs <- list(
    list(
      test_id = paste0("hyper_x_gained_anchor_", method_key),
      description = sprintf(
        "mCH hypermethylated genes against membership of gained Hi-C loop anchors (%s assignment)",
        method_label),
      col_var = "at_gained_anchor",
      gene_df = gene_tbl
    ),
    list(
      test_id = paste0("hyper_x_lost_anchor_", method_key),
      description = sprintf(
        "mCH hypermethylated genes against membership of lost Hi-C loop anchors (%s assignment)",
        method_label),
      col_var = "at_lost_anchor",
      gene_df = gene_tbl
    ),
    list(
      test_id = paste0("hyper_x_k119ub_gained_anchor_", method_key),
      description = sprintf(
        "Among anchor genes, mCH hypermethylation against a gained K119ub peak at the anchor (%s assignment)",
        method_label),
      col_var = "k119ub_gained_at_anchor",
      gene_df = gene_tbl[gene_tbl$at_any_anchor, , drop = FALSE]
    )
  )

  rows <- lapply(specs, function(spec) {
    ft <- register_fisher_test(
      section = SECTION_ID,
      test_id = spec$test_id,
      description = spec$description,
      gene_df = spec$gene_df,
      row_var = "mch_hyper",
      col_var = spec$col_var,
      output_dir = out_dir
    )
    data.frame(
      method = method_label,
      method_key = method_key,
      test_id = spec$test_id,
      description = spec$description,
      row_var = "mch_hyper",
      col_var = spec$col_var,
      n_genes = nrow(spec$gene_df),
      n_hyper = sum(spec$gene_df$mch_hyper),
      n_col_true = sum(spec$gene_df[[spec$col_var]]),
      n_both_true = sum(spec$gene_df$mch_hyper & spec$gene_df[[spec$col_var]]),
      odds_ratio = unname(ft$estimate),
      ci_lower = ft$conf.int[1],
      ci_upper = ft$conf.int[2],
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

#' Observed and expected counts for K119ub gained at anchor against mCH hyper.
compute_convergence_table <- function(gene_tbl) {
  anchor_genes <- gene_tbl[gene_tbl$at_any_anchor, , drop = FALSE]

  counts <- table(
    factor(ifelse(anchor_genes$k119ub_gained_at_anchor,
                  "K119ub gained", "K119ub not gained"),
           levels = c("K119ub gained", "K119ub not gained")),
    factor(ifelse(anchor_genes$mch_hyper,
                  "Hypermethylated", "Not hypermethylated"),
           levels = c("Hypermethylated", "Not hypermethylated"))
  )

  total <- sum(counts)
  row_totals <- rowSums(counts)
  col_totals <- colSums(counts)

  out <- expand.grid(
    k119ub_state = rownames(counts),
    mch_state = colnames(counts),
    stringsAsFactors = FALSE
  )
  row_idx <- match(out$k119ub_state, rownames(counts))
  col_idx <- match(out$mch_state, colnames(counts))

  out$observed <- as.vector(counts[cbind(row_idx, col_idx)])
  out$expected <- unname(row_totals[row_idx] * col_totals[col_idx] / total)
  out$observed_over_expected <- out$observed / out$expected
  out$method <- unique(gene_tbl$method)
  out$method_key <- unique(gene_tbl$method_key)
  out$n_anchor_genes <- total

  out[, c("method", "method_key", "k119ub_state", "mch_state", "observed",
          "expected", "observed_over_expected", "n_anchor_genes")]
}

#' Logistic model of mCH hypermethylation on anchor class and K119ub fold.
#'
#' Reference level for anchor_class is "No anchor". Genes without a K119ub
#' gene-body log2 fold change are not part of the fit; the count is returned.
fit_hyper_logistic <- function(gene_tbl) {
  model_df <- gene_tbl[!is.na(gene_tbl$k119ub_gb_log2fc), , drop = FALSE]
  model_df$anchor_class <- relevel(droplevels(model_df$anchor_class),
                                   ref = "No anchor")

  fit <- glm(mch_hyper ~ anchor_class + k119ub_gb_log2fc,
             family = binomial, data = model_df)

  method_label <- unique(gene_tbl$method)
  smry <- summary(fit)
  wald_ci <- confint.default(fit)

  coefficients <- data.frame(
    method = method_label,
    method_key = unique(gene_tbl$method_key),
    term = rownames(smry$coefficients),
    estimate = smry$coefficients[, 1],
    std_error = smry$coefficients[, 2],
    z_value = smry$coefficients[, 3],
    p_value = smry$coefficients[, 4],
    odds_ratio = exp(smry$coefficients[, 1]),
    or_lower = exp(wald_ci[, 1]),
    or_upper = exp(wald_ci[, 2]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  model_summary <- data.frame(
    method = method_label,
    method_key = unique(gene_tbl$method_key),
    formula = "mch_hyper ~ anchor_class + k119ub_gb_log2fc",
    n_genes_available = nrow(gene_tbl),
    n_genes_in_model = nrow(model_df),
    n_dropped_missing_k119ub = nrow(gene_tbl) - nrow(model_df),
    n_hyper_in_model = sum(model_df$mch_hyper),
    aic = AIC(fit),
    null_deviance = fit$null.deviance,
    residual_deviance = fit$deviance,
    df_null = fit$df.null,
    df_residual = fit$df.residual,
    stringsAsFactors = FALSE
  )

  cat(sprintf("  %s logistic model: n=%d genes, %d hypermethylated, AIC=%.1f\n",
              method_label, nrow(model_df), sum(model_df$mch_hyper), AIC(fit)))

  list(coefficients = coefficients, model_summary = model_summary)
}

#' Spearman correlation of mch_diff with the K119ub gene-body fold change.
compute_mch_k119ub_correlation <- function(gene_tbl) {
  rows <- lapply(levels(gene_tbl$anchor_class), function(cls) {
    sub <- gene_tbl[gene_tbl$anchor_class == cls &
                      !is.na(gene_tbl$k119ub_gb_log2fc) &
                      !is.na(gene_tbl$mch_diff), , drop = FALSE]
    if (nrow(sub) < 3) {
      cat(sprintf("  Anchor class '%s' holds %d usable genes; it gets no correlation.\n",
                  cls, nrow(sub)))
      return(NULL)
    }
    ct <- cor.test(sub$k119ub_gb_log2fc, sub$mch_diff,
                   method = "spearman", exact = FALSE)
    data.frame(
      method = unique(gene_tbl$method),
      method_key = unique(gene_tbl$method_key),
      anchor_class = cls,
      n_genes = nrow(sub),
      spearman_rho = unname(ct$estimate),
      p_value = ct$p.value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

# =============================================================================
# FIGURES
# =============================================================================

plot_hyper_rate_by_class <- function(rate_table, fisher_table, out_dir) {
  df <- rate_table
  df$anchor_class <- factor(df$anchor_class, levels = ANCHOR_CLASS_ORDER)
  df$bar_label <- sprintf("n = %s\n%.1f%%",
                          format(df$n_genes, big.mark = ","), df$pct_hyper)

  gained <- fisher_table[grepl("^hyper_x_gained_anchor_", fisher_table$test_id), ]
  lost <- fisher_table[grepl("^hyper_x_lost_anchor_", fisher_table$test_id), ]
  subtitle <- paste(
    sprintf("%s: gained-anchor OR = %.2f (%s), lost-anchor OR = %.2f (%s)",
            gained$method, gained$odds_ratio, fmt_p_vec(gained$p_value),
            lost$odds_ratio, fmt_p_vec(lost$p_value)),
    collapse = "\n")

  p <- ggplot(df, aes(x = anchor_class, y = pct_hyper, fill = anchor_class)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = bar_label), vjust = -0.25, size = 3.2,
              lineheight = 0.95) +
    facet_wrap(~ method, nrow = 1) +
    scale_fill_manual(values = ANCHOR_CLASS_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
    labs(
      title = "mCH hypermethylation rate by Hi-C loop anchor class",
      subtitle = subtitle,
      x = NULL,
      y = "Genes hypermethylated in mutant (%)"
    ) +
    theme_emseq() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1))

  save_multiformat_ggplot(p, file.path(out_dir, "30_01_hyper_rate_by_anchor_class"),
                          width = 11, height = 6)
}

plot_mch_diff_violin <- function(gene_tbl_all, wilcoxon_table, out_dir) {
  stats_df <- annotate_group_stats(gene_tbl_all, c("method", "anchor_class"),
                                   "mch_diff")
  stats_df$label <- group_label(stats_df, digits = MCH_DIGITS)
  stats_df$anchor_class <- factor(stats_df$anchor_class, levels = ANCHOR_CLASS_ORDER)
  stats_df$method <- factor(stats_df$method, levels = levels(gene_tbl_all$method))
  stats_df$label_y <- max(gene_tbl_all$mch_diff, na.rm = TRUE)

  key_pairs <- wilcoxon_table[
    wilcoxon_table$group1 == "Gained anchor" &
      wilcoxon_table$group2 %in% c("Lost anchor", "No anchor"), ]
  subtitle <- paste(sprintf("%s: %s vs %s %s",
                            key_pairs$method, key_pairs$group1, key_pairs$group2,
                            fmt_p_vec(key_pairs$p_value)),
                    collapse = " | ")

  p <- ggplot(gene_tbl_all, aes(x = anchor_class, y = mch_diff, fill = anchor_class)) +
    geom_violin(alpha = 0.55, scale = "width", linewidth = 0.3) +
    geom_boxplot(width = 0.14, outlier.size = 0.3, alpha = 0.85, linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_text(data = stats_df,
              aes(x = anchor_class, y = label_y, label = label),
              inherit.aes = FALSE, vjust = -0.1, size = 2.9, lineheight = 0.95) +
    facet_wrap(~ method, nrow = 1) +
    scale_fill_manual(values = ANCHOR_CLASS_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    labs(
      title = "Gene-body mCH difference by Hi-C loop anchor class",
      subtitle = subtitle,
      x = NULL,
      y = "mCH difference (mutant - control)"
    ) +
    theme_emseq() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1))

  save_multiformat_ggplot(p, file.path(out_dir, "30_01_mch_diff_by_anchor_class"),
                          width = 11, height = 6.5)
}

plot_logistic_forest <- function(coefficients, model_summaries, out_dir) {
  df <- coefficients[coefficients$term != "(Intercept)", , drop = FALSE]
  df$term_label <- sub("^anchor_class", "Anchor class: ", df$term)
  df$term_label <- sub("^k119ub_gb_log2fc$", "K119ub gene-body log2FC",
                       df$term_label)
  df$significance <- ifelse(df$p_value < 0.001, "***",
                     ifelse(df$p_value < 0.01, "**",
                     ifelse(df$p_value < 0.05, "*", "ns")))
  df$point_label <- sprintf("%.2f %s", df$odds_ratio, df$significance)

  subtitle <- paste(sprintf("%s: n = %s genes, AIC = %.1f",
                            model_summaries$method,
                            format(model_summaries$n_genes_in_model, big.mark = ","),
                            model_summaries$aic),
                    collapse = " | ")

  p <- ggplot(df, aes(x = odds_ratio, y = reorder(term_label, odds_ratio))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    geom_errorbar(aes(xmin = or_lower, xmax = or_upper), orientation = "y",
                  width = 0.18, linewidth = 0.7) +
    geom_point(size = 3, color = COLORS$direction[["Hypermethylated"]]) +
    geom_text(aes(label = point_label), vjust = -0.9, size = 3) +
    facet_wrap(~ method, nrow = 1) +
    scale_x_log10() +
    labs(
      title = "Logistic model of mCH hypermethylation",
      subtitle = paste0("mch_hyper ~ anchor_class + k119ub_gb_log2fc (reference: No anchor)\n",
                        subtitle),
      x = "Odds ratio (log scale, 95% Wald CI)",
      y = NULL
    ) +
    theme_emseq() +
    theme(axis.text.y = element_text(size = 9))

  save_multiformat_ggplot(p, file.path(out_dir, "30_01_logistic_forest"),
                          width = 12, height = 6.5)
}

plot_convergence_heatmap <- function(convergence_table, fisher_table, out_dir) {
  df <- convergence_table
  df$tile_label <- sprintf("O/E = %.2f\nn = %s", df$observed_over_expected,
                           format(df$observed, big.mark = ","))

  conv <- fisher_table[grepl("^hyper_x_k119ub_gained_anchor_", fisher_table$test_id), ]
  subtitle <- paste(sprintf("%s: OR = %.2f, %s (n = %s anchor genes)",
                            conv$method, conv$odds_ratio,
                            fmt_p_vec(conv$p_value),
                            format(conv$n_genes, big.mark = ",")),
                    collapse = " | ")

  p <- ggplot(df, aes(x = mch_state, y = k119ub_state,
                      fill = observed_over_expected)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = tile_label), size = 3.6, lineheight = 0.95) +
    facet_wrap(~ method, nrow = 1) +
    scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                         midpoint = 1, name = "Observed /\nexpected") +
    labs(
      title = "K119ub and mCH convergence at Hi-C loop anchors",
      subtitle = subtitle,
      x = "Gene-body mCH status",
      y = "K119ub peak status at the gene's anchors"
    ) +
    theme_emseq() +
    theme(panel.grid = element_blank())

  save_multiformat_ggplot(p, file.path(out_dir, "30_01_k119ub_convergence_heatmap"),
                          width = 12, height = 6)
}

plot_k119ub_signal_by_class <- function(signal_df, signal_stats, out_dir) {
  stats_df <- signal_stats
  stats_df$label <- group_label(stats_df, digits = K119UB_DIGITS)
  stats_df$anchor_class <- factor(stats_df$anchor_class, levels = ANCHOR_CLASS_ORDER)
  stats_df$method <- factor(stats_df$method, levels = levels(signal_df$method))
  stats_df$label_y <- max(signal_df$k119ub_gb_log2fc, na.rm = TRUE)

  p <- ggplot(signal_df, aes(x = anchor_class, y = k119ub_gb_log2fc,
                             fill = mch_group)) +
    geom_violin(alpha = 0.5, scale = "width", linewidth = 0.3,
                position = position_dodge(width = 0.85)) +
    geom_boxplot(width = 0.16, outlier.size = 0.3, alpha = 0.9, linewidth = 0.3,
                 position = position_dodge(width = 0.85)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_text(data = stats_df,
              aes(x = anchor_class, y = label_y, label = label,
                  group = mch_group),
              inherit.aes = FALSE, position = position_dodge(width = 0.85),
              vjust = -0.1, size = 2.6, lineheight = 0.95) +
    facet_wrap(~ method, nrow = 1) +
    scale_fill_manual(
      values = c("Hypermethylated" = COLORS$direction[["Hypermethylated"]],
                 "Not hypermethylated" = "grey70"),
      name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(
      title = "K119ub gene-body fold change by anchor class and mCH status",
      x = NULL,
      y = "K119ub gene-body log2FC (mutant / control)"
    ) +
    theme_emseq() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 20, hjust = 1))

  save_multiformat_ggplot(p, file.path(out_dir, "30_01_k119ub_signal_by_anchor_class"),
                          width = 12, height = 7)
}

plot_mch_vs_k119ub_scatter <- function(gene_tbl_all, correlation_table, out_dir) {
  df <- gene_tbl_all[!is.na(gene_tbl_all$k119ub_gb_log2fc) &
                       !is.na(gene_tbl_all$mch_diff), , drop = FALSE]
  label_df <- df[df$gene_name %in% KEY_GENES, , drop = FALSE]

  away <- correlation_table[correlation_table$anchor_class == "No anchor", ]
  at_gained <- correlation_table[correlation_table$anchor_class == "Gained anchor", ]
  subtitle <- paste(sprintf("%s: rho = %.3f at gained anchors, %.3f away from anchors",
                            at_gained$method, at_gained$spearman_rho,
                            away$spearman_rho[match(at_gained$method, away$method)]),
                    collapse = " | ")

  p <- ggplot(df, aes(x = k119ub_gb_log2fc, y = mch_diff, color = anchor_class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(size = 0.5, alpha = 0.35) +
    geom_text_repel(data = label_df, aes(label = gene_name), size = 2.8,
                    max.overlaps = 30, min.segment.length = 0,
                    show.legend = FALSE) +
    facet_wrap(~ method, nrow = 1) +
    scale_color_manual(values = ANCHOR_CLASS_COLORS, name = NULL) +
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(
      title = "Gene-body mCH change against K119ub gene-body change",
      subtitle = subtitle,
      x = "K119ub gene-body log2FC (mutant / control)",
      y = "mCH difference (mutant - control)"
    ) +
    theme_emseq() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p, file.path(out_dir, "30_01_mch_vs_k119ub_scatter"),
                          width = 12, height = 6.5)
}

# =============================================================================
# MAIN
# =============================================================================

parse_options <- function() {
  option_list <- list(
    make_option("--output-dir", type = "character", dest = "output_dir",
                default = OUTPUT_PATHS$hic,
                help = "Section output directory [default: %default]"),
    make_option("--fdr-threshold", type = "double", dest = "fdr_threshold",
                default = Q_THRESHOLD,
                help = "FDR cutoff for calling K119ub peaks gained or lost [default: %default]")
  )
  parse_args(OptionParser(option_list = option_list))
}

main <- function() {
  opt <- parse_options()
  out_dir <- opt$output_dir
  fdr_threshold <- opt$fdr_threshold
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("================================================================================\n")
  cat("SECTION 30_01: mCH DIRECTION AT DIFFERENTIAL Hi-C LOOP ANCHORS\n")
  cat("================================================================================\n")
  cat("Output directory: ", out_dir, "\n", sep = "")
  cat("K119ub peak FDR:  ", fdr_threshold, "\n\n", sep = "")

  # --- Inputs ---------------------------------------------------------------
  cat("--- Loading section inputs ---\n")
  loops <- read_significant_loops(HIC_PATHS$loops)
  k119ub_signal <- read_k119ub_gene_signal(DIFFBIND_PATHS$k119ub_gene_signal)

  loop_summary <- loops %>%
    dplyr::group_by(loop_class, category, resolution_kb) %>%
    dplyr::summarise(n_loops = dplyr::n(),
                     median_logFC = median(logFC),
                     median_loop_distance = median(loop_distance),
                     .groups = "drop") %>%
    as.data.frame()
  write_section_table(loop_summary, file.path(out_dir, "30_01_loop_summary.tsv"))

  # --- Anchors --------------------------------------------------------------
  cat("\n--- Building loop anchors ---\n")
  anchors <- build_anchor_table(loops)
  anchor_gr <- anchor_table_to_granges(anchors)
  write_section_table(anchors, file.path(out_dir, "30_01_anchor_table.tsv"))

  # --- K119ub peak status at each anchor ------------------------------------
  cat("\n--- Marking K119ub differential peaks at anchors ---\n")
  k119ub_gr <- diffbind_to_granges(k119ub_diffbind)
  k119ub_gained_peaks <- k119ub_gr[k119ub_gr$FDR < fdr_threshold & k119ub_gr$Fold > 0]
  k119ub_lost_peaks   <- k119ub_gr[k119ub_gr$FDR < fdr_threshold & k119ub_gr$Fold < 0]
  cat(sprintf("  K119ub peaks: %d gained, %d lost at FDR < %.2f\n",
              length(k119ub_gained_peaks), length(k119ub_lost_peaks),
              fdr_threshold))

  k119ub_gained_at_anchor <- countOverlaps(anchor_gr, k119ub_gained_peaks) > 0
  k119ub_lost_at_anchor   <- countOverlaps(anchor_gr, k119ub_lost_peaks) > 0
  cat(sprintf("  Anchors overlapping a gained K119ub peak: %d of %d\n",
              sum(k119ub_gained_at_anchor), length(anchor_gr)))
  cat(sprintf("  Anchors overlapping a lost K119ub peak:   %d of %d\n",
              sum(k119ub_lost_at_anchor), length(anchor_gr)))

  # --- Anchor-to-gene association, two ways ---------------------------------
  cat("\n--- Associating anchors with genes ---\n")
  mch_universe <- build_gene_universe(mch_results)

  great_domains_gr <- build_great_domain_granges()
  assoc <- list(
    great  = associate_anchors_to_genes(anchor_gr, great_domains_gr, "GREAT domain"),
    direct = associate_anchors_to_genes(anchor_gr, gene_bodies, "Direct overlap")
  )

  assoc_export <- do.call(rbind, lapply(names(assoc), function(key) {
    a <- assoc[[key]]
    data.frame(
      method = unname(ASSOCIATION_METHODS[key]),
      method_key = key,
      anchor_id = anchors$anchor_id[a$anchor_idx],
      loop_id = anchors$loop_id[a$anchor_idx],
      anchor_chr = anchors$chr[a$anchor_idx],
      anchor_start = anchors$start[a$anchor_idx],
      anchor_end = anchors$end[a$anchor_idx],
      anchor_type = anchors$anchor_type[a$anchor_idx],
      loop_class = anchors$loop_class[a$anchor_idx],
      loop_logFC = anchors$loop_logFC[a$anchor_idx],
      gene_name = a$gene_name,
      stringsAsFactors = FALSE
    )
  }))
  write_section_table(assoc_export, file.path(out_dir, "30_01_anchor_gene_associations.tsv"))

  # The GREAT universe holds only genes that have a TxDb regulatory domain,
  # because a gene without one can never be assigned to an anchor by that method.
  great_universe_genes <- unique(great_domains_gr$gene_name)
  universes <- list(
    great  = mch_universe[mch_universe$gene_name %in% great_universe_genes, , drop = FALSE],
    direct = mch_universe
  )
  cat(sprintf("  GREAT universe: %d of %d mCH genes have a TxDb regulatory domain\n",
              nrow(universes$great), nrow(mch_universe)))

  cat("\n--- Building gene-level anchor tables ---\n")
  gene_tables <- lapply(names(assoc), function(key) {
    membership <- summarise_anchor_membership(assoc[[key]], anchors,
                                              k119ub_gained_at_anchor,
                                              k119ub_lost_at_anchor)
    build_gene_anchor_table(universes[[key]], membership, k119ub_signal, key)
  })
  names(gene_tables) <- names(assoc)

  gene_tbl_all <- do.call(rbind, gene_tables)
  gene_tbl_all$method <- factor(gene_tbl_all$method,
                                levels = unname(ASSOCIATION_METHODS))
  write_section_table(gene_tbl_all, file.path(out_dir, "30_01_gene_anchor_class.tsv"))

  # --- Hypermethylation rate by anchor class --------------------------------
  cat("\n--- Hypermethylation rate by anchor class ---\n")
  rate_table <- compute_class_rates(gene_tbl_all)
  print(rate_table)
  write_section_table(rate_table, file.path(out_dir, "30_01_hyper_rate_by_anchor_class.tsv"))

  # --- Registered gene-level Fisher tests -----------------------------------
  cat("\n--- Gene-level Fisher tests ---\n")
  fisher_table <- do.call(rbind, lapply(names(gene_tables), function(key) {
    run_registered_fisher_tests(gene_tables[[key]], key, out_dir)
  }))
  write_section_table(fisher_table, file.path(out_dir, "30_01_fisher_anchor_class.tsv"))

  # --- Wilcoxon on mch_diff -------------------------------------------------
  cat("\n--- Wilcoxon tests on mch_diff by anchor class ---\n")
  wilcoxon_table <- do.call(rbind, lapply(names(gene_tables), function(key) {
    res <- pairwise_wilcoxon(gene_tables[[key]], "anchor_class", "mch_diff")
    cbind(method = unname(ASSOCIATION_METHODS[key]), method_key = key, res)
  }))
  print(wilcoxon_table)
  write_section_table(wilcoxon_table, file.path(out_dir, "30_01_wilcoxon_mch_diff.tsv"))

  # --- K119ub convergence ---------------------------------------------------
  cat("\n--- K119ub convergence at anchors ---\n")
  convergence_table <- do.call(rbind, lapply(gene_tables, compute_convergence_table))
  print(convergence_table)
  write_section_table(convergence_table, file.path(out_dir, "30_01_k119ub_convergence_counts.tsv"))

  signal_df <- gene_tbl_all[!is.na(gene_tbl_all$k119ub_gb_log2fc), , drop = FALSE]
  signal_df$mch_group <- ifelse(signal_df$mch_hyper, "Hypermethylated",
                                "Not hypermethylated")
  signal_stats <- annotate_group_stats(signal_df,
                                       c("method", "anchor_class", "mch_group"),
                                       "k119ub_gb_log2fc")
  write_section_table(signal_stats, file.path(out_dir, "30_01_k119ub_signal_by_anchor_class.tsv"))

  correlation_table <- do.call(rbind, lapply(gene_tables,
                                             compute_mch_k119ub_correlation))
  if (is.null(correlation_table) || nrow(correlation_table) == 0) {
    missing_classes <- setdiff(
      ANCHOR_CLASS_ORDER,
      if (!is.null(correlation_table)) unique(correlation_table$anchor_class) else character(0)
    )
    stop("compute_mch_k119ub_correlation(): zero rows returned. ",
         "Missing anchor classes: ", paste(missing_classes, collapse = ", "))
  }
  print(correlation_table)
  write_section_table(correlation_table, file.path(out_dir, "30_01_mch_k119ub_correlation.tsv"))

  # --- Logistic model -------------------------------------------------------
  cat("\n--- Logistic model of mCH hypermethylation ---\n")
  logistic_fits <- lapply(gene_tables, fit_hyper_logistic)
  logistic_coefficients <- do.call(rbind, lapply(logistic_fits, `[[`, "coefficients"))
  logistic_summaries <- do.call(rbind, lapply(logistic_fits, `[[`, "model_summary"))
  print(logistic_coefficients)
  write_section_table(logistic_coefficients, file.path(out_dir, "30_01_logistic_coefficients.tsv"))
  write_section_table(logistic_summaries, file.path(out_dir, "30_01_logistic_model_summary.tsv"))

  # --- Figures --------------------------------------------------------------
  cat("\n--- Figures ---\n")
  plot_hyper_rate_by_class(rate_table, fisher_table, out_dir)
  plot_mch_diff_violin(gene_tbl_all, wilcoxon_table, out_dir)
  plot_logistic_forest(logistic_coefficients, logistic_summaries, out_dir)
  plot_convergence_heatmap(convergence_table, fisher_table, out_dir)
  plot_k119ub_signal_by_class(signal_df, signal_stats, out_dir)
  plot_mch_vs_k119ub_scatter(gene_tbl_all, correlation_table, out_dir)

  cat("\n================================================================================\n")
  cat("SECTION 30_01 COMPLETE\n")
  cat("================================================================================\n\n")
}

main()
