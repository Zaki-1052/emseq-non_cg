# scripts/sections/30_hic/30_02_mecp2_loop_anchors.R
#
# Section 30_02: MeCP2 differential binding at Hi-C loop anchors.
#
# What this tests
#   MeCP2 reads methylated cytosines. In neurons the methylation it reads is
#   largely non-CG methylation (mCH) across gene bodies. This section asks
#   whether MeCP2 differential binding concentrates at the anchors of
#   differential Hi-C loops, which connects the methylation reader to 3D genome
#   organisation.
#
#   Five analyses:
#     1. MeCP2 peak density at loop anchors against a size-matched and
#        chromosome-matched background of non-anchor intervals.
#     2. Gene-level Fisher tests (registered for permutation validation in
#        section 40_04): MeCP2-gained genes against gained-anchor membership,
#        and against lost-anchor membership.
#     3. Wilcoxon tests of MeCP2 fold change by anchor class, at peak level and
#        at gene level.
#     4. The distance distribution from MeCP2 peaks to the nearest loop anchor.
#     5. Concentration-weighted MeCP2 fold per gene,
#        sum(Fold * Conc) / sum(Conc) over all peaks assigned to the gene.
#
# Reads
#   HIC_PATHS$loops   differential Hi-C loops with both anchor coordinates
#   mecp2_diffbind    MeCP2 differential binding table (pre-loaded by config)
#   mecp2_consensus   MeCP2 consensus peak set as GRanges (pre-loaded)
#   gene_bodies       gene-body GRanges of the mCH-tested gene universe
#   mch_results       gene-level mCH differential results
#
# Writes to results/sections/30_hic/ (OUTPUT_PATHS$hic, override with
#   --output-dir): five multi-format figures, fifteen TSV tables, and two
#   registered Fisher gene tables under fisher_tables/.
#
# Adapted from Biomodal section 31 (MeCP2 x Hi-C loop anchor integration).

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)

# =============================================================================
# SECTION CONSTANTS
# =============================================================================

SECTION_ID <- "30_02"

# Loop direction labels in HIC_PATHS$loops.
LOOP_GAINED <- "up_in_mutant"
LOOP_LOST   <- "down_in_mutant"

# Anchor classes. An anchor interval used by both gained and lost loops is
# "Mixed anchor". "Matched background" and "No anchor" are the comparison
# groups for peak-level and gene-level tests.
ANCHOR_CLASS_GAINED <- "Gained anchor"
ANCHOR_CLASS_LOST   <- "Lost anchor"
ANCHOR_CLASS_MIXED  <- "Mixed anchor"

PEAK_CLASS_ORDER <- c(ANCHOR_CLASS_GAINED, ANCHOR_CLASS_LOST, ANCHOR_CLASS_MIXED,
                      "Matched background", "Other genome")

REGION_CLASS_ORDER <- c(ANCHOR_CLASS_GAINED, ANCHOR_CLASS_LOST, ANCHOR_CLASS_MIXED,
                        "Matched background")

GENE_CLASS_ORDER <- c(ANCHOR_CLASS_GAINED, ANCHOR_CLASS_LOST, ANCHOR_CLASS_MIXED,
                      "No anchor")

ANCHOR_CLASS_COLORS <- c(
  "Gained anchor"      = "#4575B4",
  "Lost anchor"        = "#D73027",
  "Mixed anchor"       = "#8073AC",
  "Matched background" = "grey70",
  "No anchor"          = "grey85",
  "Other genome"       = "grey85"
)

MECP2_DIRECTION_ORDER <- c("Gained", "Lost", "Unchanged")

MECP2_DIRECTION_COLORS <- c(
  "Gained"    = unname(COLORS$mecp2["MeCP2 Up"]),
  "Lost"      = unname(COLORS$mecp2["MeCP2 Down"]),
  "Unchanged" = unname(COLORS$mecp2["Not Significant"])
)

# Columns the loop table must provide.
LOOP_REQUIRED_COLS <- c("loop_id", "chr1", "start1", "end1",
                        "chr2", "start2", "end2",
                        "logFC", "FDR", "significant", "direction")

# =============================================================================
# COMMAND LINE
# =============================================================================

parse_options <- function() {
  option_list <- list(
    make_option("--output-dir", type = "character", default = OUTPUT_PATHS$hic,
                help = "Directory for figures and tables [default: %default]"),
    make_option("--fdr-threshold", type = "double", default = Q_THRESHOLD,
                help = "FDR cutoff for MeCP2 peak significance [default: %default]"),
    make_option("--background-offset", type = "double", default = 1e6,
                help = paste("Distance in bp used to shift each anchor into the",
                             "matched background set [default: %default]"))
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
# LOOP ANCHORS
# =============================================================================

#' Read the differential Hi-C loop table and keep the significant loops.
#'
#' Every anchor class in this section is derived from the direction column, so
#' the table is reduced to loops that pass the significance flag and the
#' direction is checked against the sign of logFC.
#'
#' @param path TSV with one row per loop.
#' @return data.frame of significant loops with LOOP_REQUIRED_COLS present
read_loops <- function(path) {
  if (!file.exists(path)) {
    stop("Hi-C loop table not found: ", path,
         "\nRun scripts/utils/copy_reference_data.sh to populate data/hic/.")
  }

  loops <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                      quote = "", comment.char = "")

  missing <- setdiff(LOOP_REQUIRED_COLS, colnames(loops))
  if (length(missing) > 0) {
    stop("Loop table is missing columns: ", paste(missing, collapse = ", "),
         "\nFile: ", path)
  }

  n_all <- nrow(loops)
  loops$significant <- as.logical(loops$significant)
  loops <- loops[!is.na(loops$significant) & loops$significant, , drop = FALSE]
  if (nrow(loops) == 0) stop("No significant loops in ", path)
  cat(sprintf("  Loop table: %s rows, %s significant\n",
              format(n_all, big.mark = ","),
              format(nrow(loops), big.mark = ",")))

  unknown_direction <- setdiff(unique(loops$direction), c(LOOP_GAINED, LOOP_LOST))
  if (length(unknown_direction) > 0) {
    stop("Loop table has unexpected direction values: ",
         paste(unknown_direction, collapse = ", "),
         "\nExpected only '", LOOP_GAINED, "' and '", LOOP_LOST, "'.")
  }

  # Anchor classes come from direction alone, so direction must agree with the
  # sign of logFC. A change of convention upstream would otherwise relabel every
  # gained anchor as lost without any error.
  disagree <- (loops$direction == LOOP_GAINED) != (loops$logFC > 0)
  if (any(disagree)) {
    stop(sum(disagree), " of ", nrow(loops), " significant loops in ", path,
         " have a direction that disagrees with the sign of logFC.",
         "\nFirst offending loop_id: ", loops$loop_id[which(disagree)[1]],
         " (direction=", loops$direction[which(disagree)[1]],
         ", logFC=", loops$logFC[which(disagree)[1]], ")")
  }

  cat(sprintf("  Significant loops: %s (%s %s, %s %s)\n",
              format(nrow(loops), big.mark = ","),
              format(sum(loops$direction == LOOP_GAINED), big.mark = ","), LOOP_GAINED,
              format(sum(loops$direction == LOOP_LOST), big.mark = ","), LOOP_LOST))
  loops
}

#' Build one GRanges row per loop-anchor pair.
#'
#' Anchor coordinates in the loop table are BED-style half-open intervals with
#' 0-based starts, so the start is shifted by one to match the 1-based ranges
#' used everywhere else in the pipeline.
#'
#' @param loops data.frame from read_loops().
#' @return GRanges of length 2 * nrow(loops) carrying loop_id, loop_logfc,
#'   loop_fdr, loop_direction, anchor_slot
build_anchor_occurrences <- function(loops) {
  slot1 <- data.frame(
    loop_id = loops$loop_id,
    chr = loops$chr1, start = loops$start1, end = loops$end1,
    anchor_slot = "A1",
    stringsAsFactors = FALSE
  )
  slot2 <- data.frame(
    loop_id = loops$loop_id,
    chr = loops$chr2, start = loops$start2, end = loops$end2,
    anchor_slot = "A2",
    stringsAsFactors = FALSE
  )
  both <- rbind(slot1, slot2)
  both$loop_logfc <- rep(loops$logFC, 2)
  both$loop_fdr <- rep(loops$FDR, 2)
  both$loop_direction <- rep(loops$direction, 2)

  gr <- GRanges(
    seqnames = both$chr,
    ranges = IRanges(start = both$start + 1L, end = both$end),
    loop_id = both$loop_id,
    loop_logfc = both$loop_logfc,
    loop_fdr = both$loop_fdr,
    loop_direction = both$loop_direction,
    anchor_slot = both$anchor_slot
  )

  cat(sprintf("  Anchor occurrences: %s (2 per loop)\n",
              format(length(gr), big.mark = ",")))
  gr
}

#' Collapse anchor occurrences to unique intervals and classify them.
#'
#' An interval used only by up_in_mutant loops is a gained anchor, only by
#' down_in_mutant loops a lost anchor, and by both a mixed anchor.
#'
#' @param occurrences GRanges from build_anchor_occurrences().
#' @return GRanges of unique intervals with anchor_id, anchor_class, n_loops,
#'   n_gained_loops, n_lost_loops, mean_loop_logfc
build_anchor_intervals <- function(occurrences) {
  key <- paste0(as.character(seqnames(occurrences)), ":",
                start(occurrences), "-", end(occurrences))

  summary_df <- data.frame(
    anchor_key = key,
    chr = as.character(seqnames(occurrences)),
    start = start(occurrences),
    end = end(occurrences),
    loop_logfc = mcols(occurrences)$loop_logfc,
    is_gained = mcols(occurrences)$loop_direction == LOOP_GAINED,
    is_lost = mcols(occurrences)$loop_direction == LOOP_LOST,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(anchor_key, chr, start, end) %>%
    dplyr::summarise(
      n_loops = dplyr::n(),
      n_gained_loops = sum(is_gained),
      n_lost_loops = sum(is_lost),
      mean_loop_logfc = mean(loop_logfc),
      .groups = "drop"
    ) %>%
    as.data.frame()

  summary_df$anchor_class <- ifelse(
    summary_df$n_gained_loops > 0 & summary_df$n_lost_loops > 0, ANCHOR_CLASS_MIXED,
    ifelse(summary_df$n_gained_loops > 0, ANCHOR_CLASS_GAINED, ANCHOR_CLASS_LOST))

  gr <- GRanges(
    seqnames = summary_df$chr,
    ranges = IRanges(start = summary_df$start, end = summary_df$end),
    anchor_id = summary_df$anchor_key,
    anchor_class = summary_df$anchor_class,
    n_loops = summary_df$n_loops,
    n_gained_loops = summary_df$n_gained_loops,
    n_lost_loops = summary_df$n_lost_loops,
    mean_loop_logfc = summary_df$mean_loop_logfc
  )

  class_counts <- table(mcols(gr)$anchor_class)
  cat(sprintf("  Unique anchor intervals: %s\n",
              format(length(gr), big.mark = ",")))
  for (cls in names(class_counts)) {
    cat(sprintf("    %-16s %s\n", cls, format(class_counts[[cls]], big.mark = ",")))
  }
  gr
}

#' Build a size-matched and chromosome-matched background interval set.
#'
#' Each anchor is shifted along its own chromosome by offset_bp, keeping its
#' width. When the shifted interval would run past the chromosome end, the
#' shift is applied in the opposite direction. Shifted intervals that still
#' fall outside the chromosome, or that land on any anchor, are dropped.
#'
#' @param anchor_gr GRanges of unique anchor intervals.
#' @param offset_bp Shift distance in bp.
#' @param chrom_sizes Named integer vector of chromosome lengths.
#' @return GRanges of background intervals
build_matched_background <- function(anchor_gr, offset_bp, chrom_sizes) {
  chrs <- as.character(seqnames(anchor_gr))
  unknown <- setdiff(unique(chrs), names(chrom_sizes))
  if (length(unknown) > 0) {
    stop("Anchor chromosomes absent from the mm10 chromosome sizes: ",
         paste(unknown, collapse = ", "))
  }

  limits <- as.numeric(chrom_sizes[chrs])
  starts <- start(anchor_gr)
  ends <- end(anchor_gr)

  shifted_start <- starts + offset_bp
  shifted_end <- ends + offset_bp
  past_end <- shifted_end > limits
  shifted_start[past_end] <- starts[past_end] - offset_bp
  shifted_end[past_end] <- ends[past_end] - offset_bp

  in_bounds <- shifted_start >= 1 & shifted_end <= limits
  background <- GRanges(
    seqnames = chrs[in_bounds],
    ranges = IRanges(start = shifted_start[in_bounds], end = shifted_end[in_bounds])
  )
  background <- background[countOverlaps(background, anchor_gr) == 0]

  if (length(background) == 0) {
    stop("Matched background is empty at offset ", offset_bp,
         " bp. Choose a different --background-offset.")
  }

  cat(sprintf("  Matched background intervals: %s (from %s anchors, offset %s bp)\n",
              format(length(background), big.mark = ","),
              format(length(anchor_gr), big.mark = ","),
              format(offset_bp, big.mark = ",", scientific = FALSE)))
  background
}

# =============================================================================
# PEAK DENSITY
# =============================================================================

#' Count distinct peaks per Mb of covered sequence in a region set.
#'
#' Regions are merged with reduce() first, so overlapping regions inside one
#' class do not count their shared bases twice. Different classes are reduced
#' separately and can therefore share bases.
#'
#' @param regions GRanges for one class.
#' @param peaks GRanges of peaks.
#' @return data.frame with n_regions, covered_mb, n_peaks, peaks_per_mb
peak_density <- function(regions, peaks) {
  if (length(regions) == 0) {
    stop("peak_density() received an empty region set, so the density would ",
         "divide by zero.")
  }
  merged <- reduce(regions)
  covered_mb <- sum(as.numeric(width(merged))) / 1e6
  n_peaks <- length(subsetByOverlaps(peaks, merged))
  data.frame(
    n_regions = length(regions),
    n_merged_regions = length(merged),
    covered_mb = covered_mb,
    n_peaks = n_peaks,
    peaks_per_mb = n_peaks / covered_mb,
    stringsAsFactors = FALSE
  )
}

#' Peak density for every region class and every peak set.
#'
#' @param region_sets Named list of GRanges, one per region class.
#' @param peak_sets Named list of GRanges, one per peak set.
#' @return data.frame with one row per class x peak set
peak_density_table <- function(region_sets, peak_sets) {
  rows <- list()
  for (class_name in names(region_sets)) {
    for (peak_name in names(peak_sets)) {
      row <- peak_density(region_sets[[class_name]], peak_sets[[peak_name]])
      row$region_class <- class_name
      row$peak_set <- peak_name
      rows[[length(rows) + 1]] <- row
    }
  }
  out <- do.call(rbind, rows)
  out[, c("region_class", "peak_set", "n_regions", "n_merged_regions",
          "covered_mb", "n_peaks", "peaks_per_mb")]
}

#' Interval-level Fisher test: anchors against background for one peak set.
#'
#' Each interval is one observation and is scored for carrying at least one
#' peak. This is an interval-level test, so it calls fisher.test() directly
#' rather than register_fisher_test().
#'
#' @param anchor_gr GRanges of anchor intervals.
#' @param background_gr GRanges of background intervals.
#' @param peaks GRanges of the peak set under test.
#' @param peak_set_name Name recorded in the result row.
#' @return data.frame with counts, odds ratio, and p-value
anchor_vs_background_fisher <- function(anchor_gr, background_gr, peaks,
                                        peak_set_name) {
  anchor_hit <- countOverlaps(anchor_gr, peaks) > 0
  background_hit <- countOverlaps(background_gr, peaks) > 0

  tab <- matrix(
    c(sum(anchor_hit), sum(!anchor_hit),
      sum(background_hit), sum(!background_hit)),
    nrow = 2,
    dimnames = list(c("with_peak", "without_peak"), c("anchor", "background"))
  )
  ft <- fisher.test(tab)

  data.frame(
    peak_set = peak_set_name,
    n_anchor = length(anchor_gr),
    n_anchor_with_peak = sum(anchor_hit),
    pct_anchor_with_peak = 100 * sum(anchor_hit) / length(anchor_gr),
    n_background = length(background_gr),
    n_background_with_peak = sum(background_hit),
    pct_background_with_peak = 100 * sum(background_hit) / length(background_gr),
    odds_ratio = unname(ft$estimate),
    conf_low = ft$conf.int[1],
    conf_high = ft$conf.int[2],
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# MeCP2 PEAKS
# =============================================================================

#' Assign every MeCP2 peak to one region class.
#'
#' Priority: a peak on a mixed anchor, or on both a gained and a lost anchor,
#' is "Mixed anchor"; then gained-only, then lost-only; a peak on no anchor but
#' inside the matched background is "Matched background"; anything else is
#' "Other genome".
#'
#' @param peaks_gr GRanges of MeCP2 peaks.
#' @param anchor_gr GRanges of unique anchor intervals with anchor_class.
#' @param background_gr GRanges of matched background intervals.
#' @return character vector of PEAK_CLASS_ORDER labels
classify_peaks_by_anchor <- function(peaks_gr, anchor_gr, background_gr) {
  gained <- countOverlaps(peaks_gr,
                          anchor_gr[mcols(anchor_gr)$anchor_class == ANCHOR_CLASS_GAINED]) > 0
  lost <- countOverlaps(peaks_gr,
                        anchor_gr[mcols(anchor_gr)$anchor_class == ANCHOR_CLASS_LOST]) > 0
  mixed <- countOverlaps(peaks_gr,
                         anchor_gr[mcols(anchor_gr)$anchor_class == ANCHOR_CLASS_MIXED]) > 0
  in_background <- countOverlaps(peaks_gr, background_gr) > 0

  cls <- rep("Other genome", length(peaks_gr))
  cls[in_background & !gained & !lost & !mixed] <- "Matched background"
  cls[lost & !gained & !mixed] <- ANCHOR_CLASS_LOST
  cls[gained & !lost & !mixed] <- ANCHOR_CLASS_GAINED
  cls[mixed | (gained & lost)] <- ANCHOR_CLASS_MIXED
  cls
}

#' Distance from every peak to the nearest anchor interval.
#'
#' @param peaks_gr GRanges of peaks.
#' @param anchor_gr GRanges of anchor intervals.
#' @return numeric vector, NA where the peak's chromosome carries no anchor
distance_to_nearest_anchor <- function(peaks_gr, anchor_gr) {
  hits <- distanceToNearest(peaks_gr, anchor_gr, ignore.strand = TRUE)
  out <- rep(NA_real_, length(peaks_gr))
  out[queryHits(hits)] <- mcols(hits)$distance
  out
}

#' Attach the Conc column of the DiffBind table to ChIPseeker-annotated peaks.
#'
#' annotate_peaks_to_genes() carries only Fold, FDR, and direction through
#' ChIPseeker, so Conc is matched back on the peak coordinates. Coordinates are
#' cast to integer before they are pasted, so both sides format identically.
#'
#' @param annotated data.frame from annotate_peaks_to_genes().
#' @param diffbind data.frame from load_diffbind_flex(), holding Conc.
#' @return annotated with an added Conc column
attach_conc <- function(annotated, diffbind) {
  if (!"Conc" %in% colnames(diffbind)) {
    stop("MeCP2 DiffBind table has no Conc column, so the concentration-",
         "weighted fold cannot be computed. File: ", DIFFBIND_PATHS$mecp2)
  }

  peak_key <- paste(as.character(diffbind$Chr),
                    as.integer(diffbind$Start),
                    as.integer(diffbind$End), sep = ":")
  if (anyDuplicated(peak_key) > 0) {
    stop("MeCP2 DiffBind table has duplicated peak coordinates; ",
         "Conc cannot be matched unambiguously.")
  }

  annotated_key <- paste(as.character(annotated$seqnames),
                         as.integer(annotated$start),
                         as.integer(annotated$end), sep = ":")
  idx <- match(annotated_key, peak_key)
  if (anyNA(idx)) {
    stop(sum(is.na(idx)), " annotated MeCP2 peaks did not match a row of the ",
         "DiffBind table on coordinates.")
  }

  annotated$Conc <- diffbind$Conc[idx]
  annotated
}

#' Concentration-weighted MeCP2 fold change per gene.
#'
#' Weight each peak of a gene by its Conc, so abundant peaks dominate the gene
#' summary: sum(Fold * Conc) / sum(Conc).
#'
#' @param annotated data.frame with SYMBOL, Fold, and Conc.
#' @return data.frame with gene_name, mecp2_conc_weighted_fold, mecp2_total_conc,
#'   mecp2_n_peaks_weighted
conc_weighted_fold_by_gene <- function(annotated) {
  df <- annotated[!is.na(annotated$SYMBOL) & !is.na(annotated$Conc) &
                    !is.na(annotated$Fold), , drop = FALSE]
  if (any(df$Conc <= 0)) {
    stop("Non-positive Conc values in the MeCP2 table; the weighted fold ",
         "would not be a weighted mean.")
  }

  df %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(
      mecp2_conc_weighted_fold = sum(Fold * Conc) / sum(Conc),
      mecp2_total_conc = sum(Conc),
      mecp2_n_peaks_weighted = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::rename(gene_name = SYMBOL) %>%
    as.data.frame()
}

# =============================================================================
# GENE-LEVEL ANCHOR ASSIGNMENT
# =============================================================================

#' Gene-body overlap with each anchor class, plus the mean loop logFC.
#'
#' Genes whose name appears on more than one gene body are collapsed with any()
#' over the logical columns, so a gene counts as anchor-assigned when any of its
#' bodies overlaps an anchor.
#'
#' @param gene_gr gene_bodies GRanges.
#' @param anchor_gr GRanges of unique anchor intervals with anchor_class.
#' @param occurrences GRanges of anchor occurrences carrying loop_logfc.
#' @return data.frame with gene_name, chr, at_gained_anchor, at_lost_anchor,
#'   at_mixed_anchor, gene_anchor_class, n_anchors, mean_loop_logfc
assign_genes_to_anchors <- function(gene_gr, anchor_gr, occurrences) {
  gained_gr <- anchor_gr[mcols(anchor_gr)$anchor_class == ANCHOR_CLASS_GAINED]
  lost_gr   <- anchor_gr[mcols(anchor_gr)$anchor_class == ANCHOR_CLASS_LOST]
  mixed_gr  <- anchor_gr[mcols(anchor_gr)$anchor_class == ANCHOR_CLASS_MIXED]

  per_body <- data.frame(
    gene_name = mcols(gene_gr)$gene_name,
    chr = as.character(seqnames(gene_gr)),
    hits_gained = countOverlaps(gene_gr, gained_gr) > 0,
    hits_lost = countOverlaps(gene_gr, lost_gr) > 0,
    hits_mixed = countOverlaps(gene_gr, mixed_gr) > 0,
    n_anchors = countOverlaps(gene_gr, anchor_gr),
    stringsAsFactors = FALSE
  )

  loop_hits <- findOverlaps(gene_gr, occurrences)
  loop_logfc <- data.frame(
    gene_index = queryHits(loop_hits),
    loop_logfc = mcols(occurrences)$loop_logfc[subjectHits(loop_hits)],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(gene_index) %>%
    dplyr::summarise(mean_loop_logfc = mean(loop_logfc),
                     n_anchor_loops = dplyr::n(),
                     .groups = "drop") %>%
    as.data.frame()

  per_body$mean_loop_logfc <- NA_real_
  per_body$n_anchor_loops <- 0L
  per_body$mean_loop_logfc[loop_logfc$gene_index] <- loop_logfc$mean_loop_logfc
  per_body$n_anchor_loops[loop_logfc$gene_index] <- loop_logfc$n_anchor_loops

  collapsed <- per_body %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(
      chr = dplyr::first(chr),
      at_gained_anchor = any(hits_gained),
      at_lost_anchor = any(hits_lost),
      at_mixed_anchor = any(hits_mixed),
      n_anchors = sum(n_anchors),
      n_anchor_loops = sum(n_anchor_loops),
      mean_loop_logfc = if (all(is.na(mean_loop_logfc))) NA_real_ else
        mean(mean_loop_logfc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.frame()

  collapsed$gene_anchor_class <- "No anchor"
  collapsed$gene_anchor_class[collapsed$at_lost_anchor] <- ANCHOR_CLASS_LOST
  collapsed$gene_anchor_class[collapsed$at_gained_anchor] <- ANCHOR_CLASS_GAINED
  collapsed$gene_anchor_class[collapsed$at_mixed_anchor |
                                (collapsed$at_gained_anchor &
                                   collapsed$at_lost_anchor)] <- ANCHOR_CLASS_MIXED

  cat(sprintf("  Genes assigned to anchors: %s of %s\n",
              format(sum(collapsed$gene_anchor_class != "No anchor"), big.mark = ","),
              format(nrow(collapsed), big.mark = ",")))
  collapsed
}

# =============================================================================
# STATISTICS
# =============================================================================

#' Pairwise Wilcoxon rank sum tests between named groups.
#'
#' @param df data.frame holding the values.
#' @param class_col Grouping column name.
#' @param value_col Numeric column name.
#' @param comparisons List of length-2 character vectors naming the groups.
#' @return data.frame with one row per comparison
pairwise_wilcox <- function(df, class_col, value_col, comparisons) {
  rows <- lapply(comparisons, function(pair) {
    a <- df[[value_col]][df[[class_col]] == pair[1]]
    b <- df[[value_col]][df[[class_col]] == pair[2]]
    a <- a[!is.na(a)]
    b <- b[!is.na(b)]
    if (length(a) < 3 || length(b) < 3) {
      stop("Wilcoxon test '", pair[1], "' vs '", pair[2], "' on ", value_col,
           " has fewer than 3 values in a group (n = ", length(a), ", ",
           length(b), ").")
    }
    wt <- wilcox.test(a, b)
    data.frame(
      value = value_col,
      group1 = pair[1], group2 = pair[2],
      n1 = length(a), n2 = length(b),
      median1 = median(a), median2 = median(b),
      median_difference = median(a) - median(b),
      W = unname(wt$statistic),
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# =============================================================================
# FIGURES
# =============================================================================

#' Grouped bar chart of peak density by region class.
plot_peak_density <- function(density_df, subtitle) {
  plot_df <- density_df
  plot_df$region_class <- factor(plot_df$region_class, levels = REGION_CLASS_ORDER)
  plot_df <- plot_df[!is.na(plot_df$region_class), , drop = FALSE]

  ggplot(plot_df, aes(x = region_class, y = peaks_per_mb, fill = peak_set)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68,
             color = "black", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.1f\n(n=%s)", peaks_per_mb,
                                  format(n_peaks, big.mark = ",", trim = TRUE))),
              position = position_dodge(width = 0.75), vjust = -0.2, size = 2.8,
              lineheight = 1.0) +
    scale_fill_manual(values = c("All MeCP2 peaks" = "grey45",
                                 "MeCP2 gained" = MECP2_DIRECTION_COLORS[["Gained"]],
                                 "MeCP2 lost" = MECP2_DIRECTION_COLORS[["Lost"]]),
                      name = "Peak set") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(
      title = "MeCP2 peak density at Hi-C loop anchors",
      subtitle = subtitle,
      x = "Region class",
      y = "Distinct MeCP2 peaks per Mb of covered sequence"
    ) +
    theme_emseq() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 20, hjust = 1))
}

#' Violin plus box of a fold-change column by class, annotated with n and median.
#'
#' @param plot_df data.frame with the class column and the value column.
#' @param class_col Grouping column name.
#' @param value_col Numeric column name.
#' @param class_levels Factor level order.
#' @param stats data.frame from summarise_groups() for the same columns. It
#'   holds data only; the on-plot "n = N\nmed = M" text is built here with
#'   group_label().
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param y_label Y axis label.
#' @return ggplot object
plot_fold_violin <- function(plot_df, class_col, value_col, class_levels, stats,
                             title, subtitle, y_label) {
  plot_df <- plot_df[!is.na(plot_df[[value_col]]), , drop = FALSE]
  plot_df[[class_col]] <- factor(plot_df[[class_col]], levels = class_levels)
  plot_df <- plot_df[!is.na(plot_df[[class_col]]), , drop = FALSE]

  stats <- stats[stats[[class_col]] %in% class_levels, , drop = FALSE]
  stats[[class_col]] <- factor(stats[[class_col]], levels = class_levels)
  stats$label <- group_label(stats)

  y_limit <- as.numeric(quantile(abs(plot_df[[value_col]]), 0.99, na.rm = TRUE)) * 1.45
  stats$label_y <- y_limit * 0.78

  ggplot(plot_df, aes(x = .data[[class_col]], y = .data[[value_col]],
                      fill = .data[[class_col]])) +
    geom_violin(alpha = 0.55, scale = "width", draw_quantiles = c(0.25, 0.5, 0.75),
                linewidth = 0.3) +
    geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white", alpha = 0.85,
                 linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_text(data = stats,
              aes(x = .data[[class_col]], y = label_y, label = label),
              inherit.aes = FALSE, size = 3, lineheight = 1.1) +
    scale_fill_manual(values = ANCHOR_CLASS_COLORS, drop = FALSE) +
    coord_cartesian(ylim = c(-y_limit, y_limit)) +
    labs(title = title, subtitle = subtitle, x = "", y = y_label) +
    theme_emseq() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1))
}

#' Histogram of peak-to-anchor distance, faceted by MeCP2 direction.
plot_distance_histogram <- function(distance_df, stats_df, subtitle) {
  plot_df <- distance_df[!is.na(distance_df$distance_to_anchor), , drop = FALSE]
  plot_df$log10_distance <- log10(plot_df$distance_to_anchor + 1)
  plot_df$mecp2_direction <- factor(plot_df$mecp2_direction,
                                    levels = MECP2_DIRECTION_ORDER)

  label_df <- stats_df
  label_df$mecp2_direction <- factor(label_df$mecp2_direction,
                                     levels = MECP2_DIRECTION_ORDER)
  label_df$log10_median <- log10(label_df$median_distance + 1)
  label_df$label <- sprintf("n = %s\nmedian = %s bp\ninside anchor = %.1f%%",
                            format(label_df$n_peaks, big.mark = ",", trim = TRUE),
                            format(round(label_df$median_distance), big.mark = ",",
                                   trim = TRUE),
                            label_df$pct_inside_anchor)

  ggplot(plot_df, aes(x = log10_distance, fill = mecp2_direction)) +
    geom_histogram(bins = 60, color = "white", linewidth = 0.15) +
    geom_vline(data = label_df, aes(xintercept = log10_median),
               linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_text(data = label_df, aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = 1.05, vjust = 1.15, size = 3,
              lineheight = 1.1) +
    facet_wrap(~ mecp2_direction, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = MECP2_DIRECTION_COLORS, name = "MeCP2 peak") +
    labs(
      title = "Distance from MeCP2 peaks to the nearest Hi-C loop anchor",
      subtitle = subtitle,
      x = "log10(distance to nearest anchor + 1) in bp",
      y = "MeCP2 peaks"
    ) +
    theme_emseq() +
    theme(legend.position = "none")
}

#' Scatter of gene-level MeCP2 fold against mean loop logFC.
plot_mecp2_vs_loop <- function(scatter_df, label_df, quadrant_counts, subtitle) {
  x_max <- max(abs(scatter_df$mean_loop_logfc), na.rm = TRUE)
  y_max <- max(abs(scatter_df$mecp2_fold), na.rm = TRUE)

  quad_positions <- data.frame(
    quadrant = c("Q1", "Q2", "Q3", "Q4"),
    x = c(x_max * 0.72, -x_max * 0.72, -x_max * 0.72, x_max * 0.72),
    y = c(y_max * 0.92, y_max * 0.92, -y_max * 0.92, -y_max * 0.92),
    stringsAsFactors = FALSE
  )
  quad_positions$n <- as.integer(quadrant_counts[quad_positions$quadrant])
  quad_positions$label <- sprintf("%s\nn = %s", quad_positions$quadrant,
                                  format(quad_positions$n, big.mark = ",",
                                         trim = TRUE))

  ggplot(scatter_df, aes(x = mean_loop_logfc, y = mecp2_fold)) +
    geom_point(aes(color = gene_anchor_class), alpha = 0.55, size = 1.7) +
    geom_smooth(method = "lm", formula = y ~ x, color = "black",
                linewidth = 0.7, se = TRUE, alpha = 0.15) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
    geom_text_repel(data = label_df, aes(label = gene_name), size = 2.9,
                    fontface = "italic", color = "grey15", max.overlaps = 20,
                    segment.color = "grey60", segment.size = 0.3) +
    geom_text(data = quad_positions, aes(x = x, y = y, label = label,
                                         color = quadrant),
              inherit.aes = FALSE, size = 3.2, fontface = "bold",
              show.legend = FALSE) +
    scale_color_manual(values = c(ANCHOR_CLASS_COLORS, COLORS$quadrant),
                       breaks = GENE_CLASS_ORDER, name = "Gene anchor class") +
    labs(
      title = "MeCP2 binding change against loop strength change",
      subtitle = subtitle,
      x = "Mean loop log2 fold change over the anchors of the gene (mut / ctrl)",
      y = "MeCP2 log2 fold change, peak nearest the TSS (mut / ctrl)"
    ) +
    theme_emseq() +
    theme(legend.position = "top")
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opt <- parse_options()
  out_dir <- opt[["output-dir"]]
  fdr_threshold <- opt[["fdr-threshold"]]
  background_offset <- opt[["background-offset"]]

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("================================================================================\n")
  cat("SECTION 30_02: MeCP2 DIFFERENTIAL BINDING AT HI-C LOOP ANCHORS\n")
  cat("================================================================================\n")
  cat("Output dir:        ", out_dir, "\n", sep = "")
  cat("FDR threshold:     ", fdr_threshold, "\n", sep = "")
  cat("Background offset: ", format(background_offset, scientific = FALSE), " bp\n", sep = "")
  cat("\n")

  # --- Loop anchors ---------------------------------------------------------
  cat("Loading Hi-C loops...\n")
  loops <- read_loops(HIC_PATHS$loops)
  anchor_occurrences <- build_anchor_occurrences(loops)
  anchor_intervals <- build_anchor_intervals(anchor_occurrences)

  chrom_sizes <- seqlengths(TxDb.Mmusculus.UCSC.mm10.knownGene)
  background_intervals <- build_matched_background(anchor_intervals,
                                                   background_offset, chrom_sizes)

  # --- MeCP2 peaks ----------------------------------------------------------
  cat("\nPreparing MeCP2 peak sets...\n")
  mecp2_gr <- diffbind_to_granges(mecp2_diffbind)
  mcols(mecp2_gr)$direction <- ifelse(
    mecp2_diffbind$FDR < fdr_threshold & mecp2_diffbind$Fold > 0, "Gained",
    ifelse(mecp2_diffbind$FDR < fdr_threshold & mecp2_diffbind$Fold < 0, "Lost",
           "Unchanged"))

  mecp2_gained_gr <- mecp2_gr[mcols(mecp2_gr)$direction == "Gained"]
  mecp2_lost_gr   <- mecp2_gr[mcols(mecp2_gr)$direction == "Lost"]

  cat(sprintf("  MeCP2 peaks: %s total, %s gained, %s lost at FDR < %.3f\n",
              format(length(mecp2_gr), big.mark = ","),
              format(length(mecp2_gained_gr), big.mark = ","),
              format(length(mecp2_lost_gr), big.mark = ","),
              fdr_threshold))
  cat(sprintf("  MeCP2 consensus peaks: %s\n",
              format(length(mecp2_consensus), big.mark = ",")))

  # --- Analysis 1: peak density at anchors against matched background -------
  cat("\nAnalysis 1: MeCP2 peak density at anchors against matched background...\n")

  region_sets <- list(
    anchor_intervals[mcols(anchor_intervals)$anchor_class == ANCHOR_CLASS_GAINED],
    anchor_intervals[mcols(anchor_intervals)$anchor_class == ANCHOR_CLASS_LOST],
    anchor_intervals[mcols(anchor_intervals)$anchor_class == ANCHOR_CLASS_MIXED],
    background_intervals
  )
  names(region_sets) <- c(ANCHOR_CLASS_GAINED, ANCHOR_CLASS_LOST,
                          ANCHOR_CLASS_MIXED, "Matched background")

  peak_sets <- list(
    "All MeCP2 peaks" = mecp2_consensus,
    "MeCP2 gained" = mecp2_gained_gr,
    "MeCP2 lost" = mecp2_lost_gr
  )

  density_df <- peak_density_table(region_sets, peak_sets)
  for (i in seq_len(nrow(density_df))) {
    cat(sprintf("  %-18s %-16s %6.1f peaks/Mb (%s peaks over %.2f Mb)\n",
                density_df$region_class[i], density_df$peak_set[i],
                density_df$peaks_per_mb[i],
                format(density_df$n_peaks[i], big.mark = ","),
                density_df$covered_mb[i]))
  }

  fisher_intervals <- rbind(
    anchor_vs_background_fisher(anchor_intervals, background_intervals,
                                mecp2_consensus, "All MeCP2 peaks"),
    anchor_vs_background_fisher(anchor_intervals, background_intervals,
                                mecp2_gained_gr, "MeCP2 gained"),
    anchor_vs_background_fisher(anchor_intervals, background_intervals,
                                mecp2_lost_gr, "MeCP2 lost")
  )
  for (i in seq_len(nrow(fisher_intervals))) {
    cat(sprintf("  Fisher anchor vs background, %s: OR = %.2f, %s\n",
                fisher_intervals$peak_set[i], fisher_intervals$odds_ratio[i],
                fmt_p(fisher_intervals$p_value[i])))
  }

  density_subtitle <- sprintf(
    "Anchor vs matched background, MeCP2 gained peaks: OR = %.2f, %s | background = anchors shifted %s bp",
    fisher_intervals$odds_ratio[fisher_intervals$peak_set == "MeCP2 gained"],
    fmt_p(fisher_intervals$p_value[fisher_intervals$peak_set == "MeCP2 gained"]),
    format(background_offset, big.mark = ",", scientific = FALSE))

  p_density <- plot_peak_density(density_df, density_subtitle)
  save_multiformat_ggplot(p_density,
                          file.path(out_dir, "30_02a_mecp2_peak_density_at_anchors"),
                          width = 11, height = 8)

  # --- Analysis 2: MeCP2 fold by anchor class, peak level -------------------
  cat("\nAnalysis 2: MeCP2 peak fold change by anchor class...\n")

  peak_df <- data.frame(
    chr = as.character(seqnames(mecp2_gr)),
    start = start(mecp2_gr),
    end = end(mecp2_gr),
    mecp2_fold = mcols(mecp2_gr)$Fold,
    mecp2_fdr = mcols(mecp2_gr)$FDR,
    mecp2_direction = mcols(mecp2_gr)$direction,
    peak_class = classify_peaks_by_anchor(mecp2_gr, anchor_intervals,
                                          background_intervals),
    distance_to_anchor = distance_to_nearest_anchor(mecp2_gr, anchor_intervals),
    stringsAsFactors = FALSE
  )

  peak_class_counts <- table(peak_df$peak_class)
  for (cls in PEAK_CLASS_ORDER) {
    n_cls <- if (cls %in% names(peak_class_counts)) peak_class_counts[[cls]] else 0L
    cat(sprintf("  %-18s %s peaks\n", cls, format(n_cls, big.mark = ",")))
  }

  peak_stats <- summarise_groups(peak_df, "peak_class", "mecp2_fold")
  peak_comparisons <- list(
    c(ANCHOR_CLASS_GAINED, "Matched background"),
    c(ANCHOR_CLASS_LOST, "Matched background"),
    c(ANCHOR_CLASS_MIXED, "Matched background"),
    c(ANCHOR_CLASS_GAINED, ANCHOR_CLASS_LOST),
    c(ANCHOR_CLASS_GAINED, "Other genome"),
    c(ANCHOR_CLASS_LOST, "Other genome")
  )
  peak_wilcox <- pairwise_wilcox(peak_df, "peak_class", "mecp2_fold",
                                 peak_comparisons)
  for (i in seq_len(nrow(peak_wilcox))) {
    cat(sprintf("  Wilcoxon %s vs %s: median %.4f vs %.4f, %s\n",
                peak_wilcox$group1[i], peak_wilcox$group2[i],
                peak_wilcox$median1[i], peak_wilcox$median2[i],
                fmt_p(peak_wilcox$p_value[i])))
  }

  peak_violin_subtitle <- sprintf(
    "Wilcoxon %s vs matched background: %s | %s vs matched background: %s",
    ANCHOR_CLASS_GAINED,
    fmt_p(peak_wilcox$p_value[peak_wilcox$group1 == ANCHOR_CLASS_GAINED &
                                peak_wilcox$group2 == "Matched background"]),
    ANCHOR_CLASS_LOST,
    fmt_p(peak_wilcox$p_value[peak_wilcox$group1 == ANCHOR_CLASS_LOST &
                                peak_wilcox$group2 == "Matched background"]))

  p_peak_violin <- plot_fold_violin(
    peak_df, "peak_class", "mecp2_fold", PEAK_CLASS_ORDER, peak_stats,
    title = "MeCP2 peak fold change by loop anchor class",
    subtitle = peak_violin_subtitle,
    y_label = "MeCP2 log2 fold change (mut / ctrl)")
  save_multiformat_ggplot(p_peak_violin,
                          file.path(out_dir, "30_02b_mecp2_peak_fold_by_anchor_class"),
                          width = 10, height = 8)

  # --- Analysis 3: distance from MeCP2 peaks to the nearest anchor ----------
  cat("\nAnalysis 3: distance from MeCP2 peaks to the nearest anchor...\n")

  n_no_anchor_chr <- sum(is.na(peak_df$distance_to_anchor))
  cat(sprintf("  %s MeCP2 peaks lie on chromosomes with no loop anchor and have no distance\n",
              format(n_no_anchor_chr, big.mark = ",")))

  distance_stats <- peak_df %>%
    dplyr::filter(!is.na(distance_to_anchor)) %>%
    dplyr::group_by(mecp2_direction) %>%
    dplyr::summarise(
      n_peaks = dplyr::n(),
      median_distance = median(distance_to_anchor),
      mean_distance = mean(distance_to_anchor),
      q25_distance = unname(quantile(distance_to_anchor, 0.25)),
      q75_distance = unname(quantile(distance_to_anchor, 0.75)),
      n_inside_anchor = sum(distance_to_anchor == 0),
      pct_inside_anchor = 100 * sum(distance_to_anchor == 0) / dplyr::n(),
      .groups = "drop"
    ) %>%
    as.data.frame()

  for (i in seq_len(nrow(distance_stats))) {
    cat(sprintf("  %-10s n = %s, median = %s bp, inside anchor = %.1f%%\n",
                distance_stats$mecp2_direction[i],
                format(distance_stats$n_peaks[i], big.mark = ","),
                format(round(distance_stats$median_distance[i]), big.mark = ","),
                distance_stats$pct_inside_anchor[i]))
  }

  distance_wilcox <- pairwise_wilcox(
    peak_df[!is.na(peak_df$distance_to_anchor), , drop = FALSE],
    "mecp2_direction", "distance_to_anchor",
    list(c("Gained", "Unchanged"), c("Lost", "Unchanged"), c("Gained", "Lost")))
  for (i in seq_len(nrow(distance_wilcox))) {
    cat(sprintf("  Wilcoxon distance %s vs %s: %s\n",
                distance_wilcox$group1[i], distance_wilcox$group2[i],
                fmt_p(distance_wilcox$p_value[i])))
  }

  distance_subtitle <- sprintf(
    "Wilcoxon gained vs unchanged: %s | lost vs unchanged: %s | dashed line = group median",
    fmt_p(distance_wilcox$p_value[distance_wilcox$group1 == "Gained" &
                                    distance_wilcox$group2 == "Unchanged"]),
    fmt_p(distance_wilcox$p_value[distance_wilcox$group1 == "Lost" &
                                    distance_wilcox$group2 == "Unchanged"]))

  p_distance <- plot_distance_histogram(peak_df, distance_stats, distance_subtitle)
  save_multiformat_ggplot(p_distance,
                          file.path(out_dir, "30_02d_mecp2_distance_to_nearest_anchor"),
                          width = 9, height = 10)

  # --- Gene-level MeCP2 summaries ------------------------------------------
  cat("\nAggregating MeCP2 peaks to genes...\n")
  mecp2_annotated <- annotate_peaks_to_genes(mecp2_diffbind, "MeCP2")
  mecp2_annotated <- attach_conc(mecp2_annotated, mecp2_diffbind)

  mecp2_gene <- aggregate_diffbind_by_gene(mecp2_annotated, method = "nearest_tss",
                                           fdr_threshold = fdr_threshold,
                                           prefix = "mecp2")
  cat(sprintf("  %s genes carry a MeCP2 peak\n",
              format(nrow(mecp2_gene), big.mark = ",")))

  mecp2_conc_gene <- conc_weighted_fold_by_gene(mecp2_annotated)
  cat(sprintf("  %s genes have a concentration-weighted MeCP2 fold\n",
              format(nrow(mecp2_conc_gene), big.mark = ",")))

  # --- Gene-level anchor assignment ----------------------------------------
  cat("\nAssigning genes to loop anchors...\n")
  gene_anchor <- assign_genes_to_anchors(gene_bodies, anchor_intervals,
                                         anchor_occurrences)

  mch_gene <- mch_results %>%
    dplyr::group_by(gene_name) %>%
    dplyr::slice_max(order_by = abs(edger_logFC), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(gene_name, gene_length, mch_diff, edger_logFC, edger_fdr,
                  mch_sig, mch_direction) %>%
    as.data.frame()

  gene_table <- gene_anchor %>%
    dplyr::left_join(mch_gene, by = "gene_name") %>%
    dplyr::left_join(mecp2_gene, by = "gene_name") %>%
    dplyr::left_join(mecp2_conc_gene, by = "gene_name") %>%
    as.data.frame()

  gene_table$has_mecp2_peak <- !is.na(gene_table$mecp2_fold)
  gene_table$mecp2_gained <- gene_table$has_mecp2_peak &
    !is.na(gene_table$mecp2_fdr) &
    gene_table$mecp2_fdr < fdr_threshold & gene_table$mecp2_fold > 0
  gene_table$mecp2_lost <- gene_table$has_mecp2_peak &
    !is.na(gene_table$mecp2_fdr) &
    gene_table$mecp2_fdr < fdr_threshold & gene_table$mecp2_fold < 0

  cat(sprintf("  %s genes with a MeCP2 peak, %s MeCP2-gained, %s MeCP2-lost\n",
              format(sum(gene_table$has_mecp2_peak), big.mark = ","),
              format(sum(gene_table$mecp2_gained), big.mark = ","),
              format(sum(gene_table$mecp2_lost), big.mark = ",")))

  # --- Analysis 4: registered gene-level Fisher tests -----------------------
  cat("\nAnalysis 4: gene-level Fisher tests (registered for section 40_04)...\n")

  fisher_input <- gene_table[gene_table$has_mecp2_peak, , drop = FALSE]
  cat(sprintf("  Fisher gene universe: %s genes with both mCH testing and a MeCP2 peak\n",
              format(nrow(fisher_input), big.mark = ",")))

  ft_gained_anchor <- register_fisher_test(
    section = SECTION_ID, test_id = "mecp2_gained_at_gained_anchor",
    description = paste("Genes with a MeCP2 peak: is significant MeCP2 gain",
                        "associated with lying at a gained Hi-C loop anchor?"),
    gene_df = fisher_input, row_var = "mecp2_gained", col_var = "at_gained_anchor",
    output_dir = out_dir)

  ft_lost_anchor <- register_fisher_test(
    section = SECTION_ID, test_id = "mecp2_gained_at_lost_anchor",
    description = paste("Genes with a MeCP2 peak: is significant MeCP2 gain",
                        "associated with lying at a lost Hi-C loop anchor?"),
    gene_df = fisher_input, row_var = "mecp2_gained", col_var = "at_lost_anchor",
    output_dir = out_dir)

  gene_fisher_summary <- data.frame(
    test_id = c("mecp2_gained_at_gained_anchor", "mecp2_gained_at_lost_anchor"),
    anchor_class = c(ANCHOR_CLASS_GAINED, ANCHOR_CLASS_LOST),
    n_genes = nrow(fisher_input),
    n_mecp2_gained = sum(fisher_input$mecp2_gained),
    n_at_anchor = c(sum(fisher_input$at_gained_anchor),
                    sum(fisher_input$at_lost_anchor)),
    n_both = c(sum(fisher_input$mecp2_gained & fisher_input$at_gained_anchor),
               sum(fisher_input$mecp2_gained & fisher_input$at_lost_anchor)),
    odds_ratio = c(unname(ft_gained_anchor$estimate), unname(ft_lost_anchor$estimate)),
    conf_low = c(ft_gained_anchor$conf.int[1], ft_lost_anchor$conf.int[1]),
    conf_high = c(ft_gained_anchor$conf.int[2], ft_lost_anchor$conf.int[2]),
    p_value = c(ft_gained_anchor$p.value, ft_lost_anchor$p.value),
    stringsAsFactors = FALSE
  )

  # --- Analysis 5: gene-level MeCP2 fold by anchor class --------------------
  cat("\nAnalysis 5: gene-level MeCP2 fold change by anchor class...\n")

  gene_fold_df <- gene_table[gene_table$has_mecp2_peak, , drop = FALSE]

  gene_stats <- summarise_groups(gene_fold_df, "gene_anchor_class", "mecp2_fold")
  conc_stats <- summarise_groups(gene_fold_df, "gene_anchor_class",
                                 "mecp2_conc_weighted_fold")

  gene_comparisons <- list(
    c(ANCHOR_CLASS_GAINED, "No anchor"),
    c(ANCHOR_CLASS_LOST, "No anchor"),
    c(ANCHOR_CLASS_MIXED, "No anchor"),
    c(ANCHOR_CLASS_GAINED, ANCHOR_CLASS_LOST)
  )
  gene_wilcox <- pairwise_wilcox(gene_fold_df, "gene_anchor_class", "mecp2_fold",
                                 gene_comparisons)
  conc_wilcox <- pairwise_wilcox(gene_fold_df, "gene_anchor_class",
                                 "mecp2_conc_weighted_fold", gene_comparisons)

  for (i in seq_len(nrow(gene_wilcox))) {
    cat(sprintf("  Wilcoxon gene fold %s vs %s: median %.4f vs %.4f, %s\n",
                gene_wilcox$group1[i], gene_wilcox$group2[i],
                gene_wilcox$median1[i], gene_wilcox$median2[i],
                fmt_p(gene_wilcox$p_value[i])))
  }
  for (i in seq_len(nrow(conc_wilcox))) {
    cat(sprintf("  Wilcoxon Conc-weighted fold %s vs %s: median %.4f vs %.4f, %s\n",
                conc_wilcox$group1[i], conc_wilcox$group2[i],
                conc_wilcox$median1[i], conc_wilcox$median2[i],
                fmt_p(conc_wilcox$p_value[i])))
  }

  p_gene_nearest <- plot_fold_violin(
    gene_fold_df, "gene_anchor_class", "mecp2_fold", GENE_CLASS_ORDER, gene_stats,
    title = "Gene-level MeCP2 fold change by anchor class",
    subtitle = sprintf("Peak nearest the TSS | Wilcoxon %s vs no anchor: %s",
                       ANCHOR_CLASS_GAINED,
                       fmt_p(gene_wilcox$p_value[gene_wilcox$group1 == ANCHOR_CLASS_GAINED &
                                                   gene_wilcox$group2 == "No anchor"])),
    y_label = "MeCP2 log2 fold change (mut / ctrl)")

  p_gene_conc <- plot_fold_violin(
    gene_fold_df, "gene_anchor_class", "mecp2_conc_weighted_fold", GENE_CLASS_ORDER,
    conc_stats,
    title = "Concentration-weighted MeCP2 fold change by anchor class",
    subtitle = sprintf("sum(Fold x Conc) / sum(Conc) over all peaks of the gene | Wilcoxon %s vs no anchor: %s",
                       ANCHOR_CLASS_GAINED,
                       fmt_p(conc_wilcox$p_value[conc_wilcox$group1 == ANCHOR_CLASS_GAINED &
                                                   conc_wilcox$group2 == "No anchor"])),
    y_label = "Concentration-weighted MeCP2 log2 fold change")

  p_gene_violin <- p_gene_nearest + p_gene_conc + patchwork::plot_layout(ncol = 2)
  save_multiformat_ggplot(p_gene_violin,
                          file.path(out_dir, "30_02c_mecp2_gene_fold_by_anchor_class"),
                          width = 15, height = 8)

  # --- Analysis 6: MeCP2 fold against loop logFC ----------------------------
  cat("\nAnalysis 6: gene-level MeCP2 fold against loop logFC...\n")

  scatter_df <- gene_fold_df[!is.na(gene_fold_df$mean_loop_logfc) &
                               !is.na(gene_fold_df$mecp2_fold), , drop = FALSE]
  if (nrow(scatter_df) < 10) {
    stop("Only ", nrow(scatter_df), " anchor-assigned genes carry a MeCP2 fold; ",
         "the scatter cannot be built.")
  }

  spearman <- cor.test(scatter_df$mean_loop_logfc, scatter_df$mecp2_fold,
                       method = "spearman", exact = FALSE)
  cat(sprintf("  %s anchor-assigned genes, Spearman rho = %.3f, %s\n",
              format(nrow(scatter_df), big.mark = ","),
              unname(spearman$estimate), fmt_p(spearman$p.value)))

  scatter_df$quadrant <- assign_quadrant(scatter_df$mean_loop_logfc,
                                         scatter_df$mecp2_fold)
  quadrant_counts <- table(factor(scatter_df$quadrant,
                                  levels = c("Q1", "Q2", "Q3", "Q4")))
  for (q in names(quadrant_counts)) {
    cat(sprintf("    %s: %s genes\n", q,
                format(quadrant_counts[[q]], big.mark = ",")))
  }

  key_gene_rows <- scatter_df[scatter_df$gene_name %in% KEY_GENES, , drop = FALSE]
  top_fold_rows <- scatter_df[order(-abs(scatter_df$mecp2_fold)), , drop = FALSE]
  top_fold_rows <- head(top_fold_rows, 15)
  label_df <- unique(rbind(key_gene_rows, top_fold_rows))
  cat(sprintf("  Labelling %d genes (%d of %d KEY_GENES present)\n",
              nrow(label_df), nrow(key_gene_rows), length(KEY_GENES)))

  scatter_subtitle <- sprintf(
    "Spearman rho = %.3f, %s | %s genes at loop anchors with a MeCP2 peak",
    unname(spearman$estimate), fmt_p(spearman$p.value),
    format(nrow(scatter_df), big.mark = ","))

  p_scatter <- plot_mecp2_vs_loop(scatter_df, label_df, quadrant_counts,
                                  scatter_subtitle)
  save_multiformat_ggplot(p_scatter,
                          file.path(out_dir, "30_02e_mecp2_fold_vs_loop_logfc"),
                          width = 11, height = 9)

  scatter_stats <- data.frame(
    n_genes = nrow(scatter_df),
    spearman_rho = unname(spearman$estimate),
    spearman_p = spearman$p.value,
    n_Q1 = quadrant_counts[["Q1"]],
    n_Q2 = quadrant_counts[["Q2"]],
    n_Q3 = quadrant_counts[["Q3"]],
    n_Q4 = quadrant_counts[["Q4"]],
    stringsAsFactors = FALSE
  )

  # --- Per-loop anchor overlap table ---------------------------------------
  anchor1_gr <- GRanges(seqnames = loops$chr1,
                        ranges = IRanges(start = loops$start1 + 1L, end = loops$end1))
  anchor2_gr <- GRanges(seqnames = loops$chr2,
                        ranges = IRanges(start = loops$start2 + 1L, end = loops$end2))

  loop_overlap <- data.frame(
    loop_id = loops$loop_id,
    chr1 = loops$chr1, start1 = loops$start1, end1 = loops$end1,
    chr2 = loops$chr2, start2 = loops$start2, end2 = loops$end2,
    loop_logfc = loops$logFC, loop_fdr = loops$FDR, loop_direction = loops$direction,
    anchor1_mecp2_gained = countOverlaps(anchor1_gr, mecp2_gained_gr) > 0,
    anchor1_mecp2_lost = countOverlaps(anchor1_gr, mecp2_lost_gr) > 0,
    anchor2_mecp2_gained = countOverlaps(anchor2_gr, mecp2_gained_gr) > 0,
    anchor2_mecp2_lost = countOverlaps(anchor2_gr, mecp2_lost_gr) > 0,
    stringsAsFactors = FALSE
  )
  loop_overlap$any_mecp2_gained <- loop_overlap$anchor1_mecp2_gained |
    loop_overlap$anchor2_mecp2_gained
  loop_overlap$any_mecp2_lost <- loop_overlap$anchor1_mecp2_lost |
    loop_overlap$anchor2_mecp2_lost

  cat(sprintf("\n  Loops with a MeCP2-gained peak at either anchor: %s of %s (%.1f%%)\n",
              format(sum(loop_overlap$any_mecp2_gained), big.mark = ","),
              format(nrow(loop_overlap), big.mark = ","),
              100 * sum(loop_overlap$any_mecp2_gained) / nrow(loop_overlap)))
  cat(sprintf("  Loops with a MeCP2-lost peak at either anchor:   %s of %s (%.1f%%)\n",
              format(sum(loop_overlap$any_mecp2_lost), big.mark = ","),
              format(nrow(loop_overlap), big.mark = ","),
              100 * sum(loop_overlap$any_mecp2_lost) / nrow(loop_overlap)))

  # --- Tables ---------------------------------------------------------------
  cat("\nWriting tables...\n")

  anchor_export <- data.frame(
    anchor_id = mcols(anchor_intervals)$anchor_id,
    chr = as.character(seqnames(anchor_intervals)),
    start = start(anchor_intervals),
    end = end(anchor_intervals),
    width = width(anchor_intervals),
    anchor_class = mcols(anchor_intervals)$anchor_class,
    n_loops = mcols(anchor_intervals)$n_loops,
    n_gained_loops = mcols(anchor_intervals)$n_gained_loops,
    n_lost_loops = mcols(anchor_intervals)$n_lost_loops,
    mean_loop_logfc = mcols(anchor_intervals)$mean_loop_logfc,
    n_mecp2_peaks = countOverlaps(anchor_intervals, mecp2_consensus),
    n_mecp2_gained = countOverlaps(anchor_intervals, mecp2_gained_gr),
    n_mecp2_lost = countOverlaps(anchor_intervals, mecp2_lost_gr),
    stringsAsFactors = FALSE
  )

  background_export <- data.frame(
    chr = as.character(seqnames(background_intervals)),
    start = start(background_intervals),
    end = end(background_intervals),
    width = width(background_intervals),
    n_mecp2_peaks = countOverlaps(background_intervals, mecp2_consensus),
    n_mecp2_gained = countOverlaps(background_intervals, mecp2_gained_gr),
    n_mecp2_lost = countOverlaps(background_intervals, mecp2_lost_gr),
    stringsAsFactors = FALSE
  )

  write_section_table(anchor_export, file.path(out_dir, "30_02_anchor_intervals.tsv"))
  write_section_table(background_export, file.path(out_dir, "30_02_matched_background_intervals.tsv"))
  write_section_table(loop_overlap, file.path(out_dir, "30_02_loop_anchor_mecp2_overlap.tsv"))
  write_section_table(density_df, file.path(out_dir, "30_02_peak_density_by_class.tsv"))
  write_section_table(fisher_intervals, file.path(out_dir, "30_02_anchor_vs_background_fisher.tsv"))
  write_section_table(peak_df, file.path(out_dir, "30_02_mecp2_peaks_anchor_class.tsv"))
  write_section_table(peak_stats, file.path(out_dir, "30_02_peak_fold_by_class_summary.tsv"))
  write_section_table(peak_wilcox, file.path(out_dir, "30_02_peak_fold_by_class_wilcoxon.tsv"))
  write_section_table(distance_stats, file.path(out_dir, "30_02_distance_to_anchor_summary.tsv"))
  write_section_table(distance_wilcox, file.path(out_dir, "30_02_distance_to_anchor_wilcoxon.tsv"))
  write_section_table(gene_table, file.path(out_dir, "30_02_gene_anchor_table.tsv"))
  write_section_table(gene_stats, file.path(out_dir, "30_02_gene_fold_by_class_summary.tsv"))
  write_section_table(gene_wilcox, file.path(out_dir, "30_02_gene_fold_by_class_wilcoxon.tsv"))
  write_section_table(conc_stats, file.path(out_dir, "30_02_conc_weighted_fold_by_class_summary.tsv"))
  write_section_table(conc_wilcox, file.path(out_dir, "30_02_conc_weighted_fold_by_class_wilcoxon.tsv"))
  write_section_table(gene_fisher_summary, file.path(out_dir, "30_02_gene_fisher_summary.tsv"))
  write_section_table(scatter_stats, file.path(out_dir, "30_02_mecp2_vs_loop_logfc_stats.tsv"))

  # --- Summary --------------------------------------------------------------
  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 30_02 SUMMARY\n")
  cat("================================================================================\n")
  cat(sprintf("Loops:                    %s\n", format(nrow(loops), big.mark = ",")))
  cat(sprintf("Unique anchor intervals:  %s\n",
              format(length(anchor_intervals), big.mark = ",")))
  cat(sprintf("Matched background:       %s intervals\n",
              format(length(background_intervals), big.mark = ",")))
  cat(sprintf("MeCP2 peaks:              %s (%s gained, %s lost)\n",
              format(length(mecp2_gr), big.mark = ","),
              format(length(mecp2_gained_gr), big.mark = ","),
              format(length(mecp2_lost_gr), big.mark = ",")))
  cat(sprintf("Genes with a MeCP2 peak:  %s\n",
              format(sum(gene_table$has_mecp2_peak), big.mark = ",")))
  cat(sprintf("Fisher MeCP2 gain x gained anchor: OR = %.3f, %s\n",
              unname(ft_gained_anchor$estimate), fmt_p(ft_gained_anchor$p.value)))
  cat(sprintf("Fisher MeCP2 gain x lost anchor:   OR = %.3f, %s\n",
              unname(ft_lost_anchor$estimate), fmt_p(ft_lost_anchor$p.value)))
  cat(sprintf("Spearman MeCP2 fold vs loop logFC: rho = %.3f, %s\n",
              unname(spearman$estimate), fmt_p(spearman$p.value)))
  cat("\nSection 30_02 complete.\n\n")
}

main()
