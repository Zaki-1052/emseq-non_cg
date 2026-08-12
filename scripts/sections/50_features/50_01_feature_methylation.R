# scripts/sections/50_features/50_01_feature_methylation.R
#
# Section 50_01 -- non-CG methylation across sub-gene feature types.
#
# The analysis asks where inside a gene the BAP1 knockout changes mCH. Every
# protein-coding gene is split into six interval classes (5' UTR, exon, splice
# donor window, splice acceptor window, intron, 3' UTR). Read counts are pooled
# over the four samples of each genotype for every interval, so each interval
# carries one control rate, one mutant rate, and their difference. The script
# then compares those differences between feature types, along exon and intron
# rank position, across the scaled gene body, and between intervals that do and
# do not overlap an active enhancer.
#
# Reads:
#   DATA_PATHS$feature_dir            <sample>_feature_mch.tsv for the eight
#                                     samples in SAMPLE_META, written by step
#                                     02b (scripts/02b_mch_aggregate_features.R)
#                                     over the BEDs from
#                                     scripts/utils/generate_feature_beds.R
#   CHROMATIN_PATHS$active_enhancer   chromHMM active enhancer segmentation
#   mch_results                       gene-level mCH differential results, from
#                                     the shared config
#
# Writes into OUT_DIR (default results/sections/50_features/):
#   50_01_feature_intervals.tsv            one row per interval, both genotypes
#   50_01_gene_feature_delta.tsv           one row per gene and feature type
#   50_01_feature_type_summary.tsv         counts, bases, and pooled rates
#   50_01_coverage_threshold_sweep.tsv     intervals removed at each cutoff
#   50_01_delta_by_feature_stats.tsv       n, median, mean, quartiles, figure b
#   50_01_mch_level_by_feature_stats.tsv   n, median, mean, quartiles, figure h
#   50_01_kruskal_wallis.tsv               feature type effect on delta_mch
#   50_01_dunn_posthoc.tsv                 pairwise feature type comparisons
#   50_01_wilcoxon_vs_zero.tsv             delta_mch against zero per feature
#   50_01_exon_intron_rank_profile.tsv     mCH by rank position
#   50_01_metagene_profile.tsv             mCH across the scaled gene body
#   50_01_enhancer_overlap_summary.tsv     counts and rates per enhancer group
#   50_01_enhancer_split_stats.tsv         n, median, mean, quartiles, figure f
#   50_01_enhancer_wilcoxon.tsv            enhancer group comparison
#   50_01_key_gene_intervals.tsv           every interval of the key genes
#   50_01_gene_fisher_summary.tsv          the three registered gene-level tests
#   eight figures (50_01a to 50_01h) in PDF, SVG, PNG, and JPEG
#   fisher_tables/ plus one row per test in the shared Fisher registry
#
# Adapted from the Biomodal script section_52_cpg_resolution_gene_body.R. That
# script measured CG methylation with 5mC and 5hmC as two modalities. This one
# measures mCH as a single modality, so the Biomodal panel 52c and the second
# overlay series in panels 52e, 52f, and 52g carry no counterpart here. The
# Biomodal metagene came from a separate binned BED; this script builds the
# metagene from the exon and intron intervals themselves.

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)
library(dunn.test)

# =============================================================================
# SECTION CONSTANTS
# =============================================================================

SECTION_ID <- "50_01"

# Columns every per-sample table from step 02b must carry.
REQUIRED_FEATURE_COLS <- c("sample_id", "gene_name", "feature_type",
                           "feature_rank", "chr", "start", "end", "strand",
                           "feature_length", "n_ch_sites", "total_coverage",
                           "methylated_count", "mch_rate")

# Exons and introns tile the gene body without overlapping each other. UTR
# intervals sit inside exons and splice windows straddle exon boundaries, so
# the metagene and the gene body span use exons and introns only.
BODY_FEATURES <- c("Exon", "Intron")

# Coverage cutoffs reported in the threshold sweep. The cutoff actually applied
# comes from --min-coverage and is added to this list.
COVERAGE_SWEEP_THRESHOLDS <- c(0, 5, 10, 20, 30, 50, 100, 200)

FEATURE_COLORS <- c(
  "5UTR"                = "#66C2A5",
  "Exon"                = "#FC8D62",
  "SpliceSite_Donor"    = "#8DA0CB",
  "SpliceSite_Acceptor" = "#A6D854",
  "Intron"              = "#E78AC3",
  "3UTR"                = "#FFD92F"
)

ENHANCER_LEVELS <- c("Active enhancer", "No enhancer")

ENHANCER_COLORS <- c(
  "Active enhancer" = unname(BODY_STATE_COLORS["Enhancer_Marked"]),
  "No enhancer"     = "grey70"
)

GENOTYPE_LABELS <- c(ctrl = "Control", mut = "Mutant")

# =============================================================================
# OPTIONS AND INPUT CHECKS
# =============================================================================

parse_options <- function() {
  option_list <- list(
    make_option("--output-dir", dest = "output_dir", type = "character",
                default = OUTPUT_PATHS$features,
                help = "Directory for figures and tables [default: %default]"),
    make_option("--min-coverage", dest = "min_coverage", type = "double",
                default = 20,
                help = paste("Minimum pooled read coverage an interval needs in",
                             "each genotype before its rate is used",
                             "[default: %default]")),
    make_option("--max-rank", dest = "max_rank", type = "integer",
                default = 20,
                help = paste("Highest feature rank drawn in the exon and intron",
                             "rank profiles [default: %default]")),
    make_option("--metagene-bins", dest = "metagene_bins", type = "integer",
                default = 50,
                help = "Bins across the scaled gene body [default: %default]")
  )
  opt <- parse_args(OptionParser(
    option_list = option_list,
    description = "mCH across sub-gene feature types in BAP1-KO cerebellum."
  ))

  if (opt$min_coverage < 1) stop("--min-coverage must be at least 1")
  if (opt$max_rank < 2) stop("--max-rank must be at least 2")
  if (opt$metagene_bins < 5) stop("--metagene-bins must be at least 5")
  opt
}

#' Absolute path of one sample's step 02b feature table.
feature_table_path <- function(sample_id) {
  file.path(DATA_PATHS$feature_dir, paste0(sample_id, "_feature_mch.tsv"))
}

#' Stop unless every section input exists.
check_inputs <- function() {
  paths <- vapply(SAMPLE_META$sample_id, feature_table_path, character(1))
  missing <- SAMPLE_META$sample_id[!file.exists(paths)]
  if (length(missing) > 0) {
    stop("Per-sample feature tables not found for: ",
         paste(missing, collapse = ", "),
         "\n  Expected in: ", DATA_PATHS$feature_dir,
         "\n  Build the interval BEDs with scripts/utils/generate_feature_beds.R,",
         "\n  then run step 02b (scripts/02b_mch_aggregate_features.R) for every sample.")
  }

  if (!file.exists(CHROMATIN_PATHS$active_enhancer)) {
    stop("Active enhancer BED not found: ", CHROMATIN_PATHS$active_enhancer)
  }

  genotypes <- unique(SAMPLE_META$genotype)
  if (!setequal(genotypes, c("ctrl", "mut"))) {
    stop("SAMPLE_META must hold both ctrl and mut samples. Found: ",
         paste(genotypes, collapse = ", "))
  }
  invisible(TRUE)
}

#' Format an integer for a plot label or a log line.
fmt_int <- function(x) format(x, big.mark = ",", trim = TRUE)

# =============================================================================
# LOADING AND POOLING
# =============================================================================

#' Read one per-sample feature table from step 02b.
read_feature_table <- function(sample_id) {
  path <- feature_table_path(sample_id)
  dt <- data.table::fread(path, sep = "\t", showProgress = FALSE)

  missing <- setdiff(REQUIRED_FEATURE_COLS, colnames(dt))
  if (length(missing) > 0) {
    stop(basename(path), " is missing columns: ", paste(missing, collapse = ", "),
         "\n  Rerun step 02b (scripts/02b_mch_aggregate_features.R) for ", sample_id)
  }
  if (nrow(dt) == 0) stop(basename(path), " holds no intervals.")
  if (!all(dt$sample_id == sample_id)) {
    stop(basename(path), " carries sample_id values other than ", sample_id)
  }
  dt
}

#' Pool read counts over the samples of each genotype, one row per interval.
#'
#' Every sample table describes the same interval set in the same order,
#' because step 02b writes one row per BED interval. The interval key is
#' compared between samples and a mismatch stops the run.
#'
#' @return data.table with the interval metadata, per-genotype pooled counts,
#'   per-genotype rates, and delta_mch as mutant minus control
load_interval_counts <- function() {
  cat("Loading per-sample feature tables from:\n  ", DATA_PATHS$feature_dir, "\n",
      sep = "")

  intervals <- NULL
  key <- NULL
  pooled <- list(
    ctrl = list(coverage = NULL, methylated = NULL, sites = NULL),
    mut  = list(coverage = NULL, methylated = NULL, sites = NULL)
  )

  for (i in seq_len(nrow(SAMPLE_META))) {
    sample_id <- SAMPLE_META$sample_id[i]
    genotype  <- SAMPLE_META$genotype[i]
    dt <- read_feature_table(sample_id)

    this_key <- paste(dt$gene_name, dt$feature_type, dt$chr, dt$start, dt$end,
                      sep = "|")

    if (is.null(intervals)) {
      intervals <- dt[, .(gene_name, feature_type, feature_rank, chr, start, end,
                          strand, feature_length)]
      key <- this_key
      zeros <- numeric(nrow(dt))
      for (g in names(pooled)) {
        pooled[[g]]$coverage   <- zeros
        pooled[[g]]$methylated <- zeros
        pooled[[g]]$sites      <- zeros
      }
    } else if (!identical(this_key, key)) {
      stop("Interval rows in ", basename(feature_table_path(sample_id)),
           " do not match ", basename(feature_table_path(SAMPLE_META$sample_id[1])),
           ".\n  All samples must come from the same feature BEDs. Rerun step 02b",
           " with one feature directory.")
    }

    pooled[[genotype]]$coverage   <- pooled[[genotype]]$coverage + dt$total_coverage
    pooled[[genotype]]$methylated <- pooled[[genotype]]$methylated + dt$methylated_count
    pooled[[genotype]]$sites      <- pooled[[genotype]]$sites + dt$n_ch_sites

    cat(sprintf("  %-8s %s intervals  %s CH site observations  %s read coverage\n",
                sample_id, fmt_int(nrow(dt)),
                fmt_int(sum(dt$n_ch_sites)), fmt_int(sum(dt$total_coverage))))
  }

  intervals[, `:=`(
    ctrl_coverage   = pooled$ctrl$coverage,
    ctrl_methylated = pooled$ctrl$methylated,
    ctrl_ch_sites   = pooled$ctrl$sites,
    mut_coverage    = pooled$mut$coverage,
    mut_methylated  = pooled$mut$methylated,
    mut_ch_sites    = pooled$mut$sites
  )]

  intervals[, mch_ctrl := ifelse(ctrl_coverage > 0, ctrl_methylated / ctrl_coverage,
                                 NA_real_)]
  intervals[, mch_mut := ifelse(mut_coverage > 0, mut_methylated / mut_coverage,
                                NA_real_)]
  intervals[, delta_mch := mch_mut - mch_ctrl]

  unknown <- setdiff(unique(intervals$feature_type), FEATURE_TYPES)
  if (length(unknown) > 0) {
    stop("Feature types outside FEATURE_TYPES found: ",
         paste(unknown, collapse = ", "))
  }
  intervals[, feature_type := factor(feature_type, levels = FEATURE_TYPES)]

  cat(sprintf("\n  %s intervals across %s genes and %d feature types\n",
              fmt_int(nrow(intervals)), fmt_int(uniqueN(intervals$gene_name)),
              nlevels(intervals$feature_type)))
  intervals[]
}

# =============================================================================
# COVERAGE FILTER
# =============================================================================

#' Count the intervals each candidate coverage cutoff removes.
#'
#' An interval passes a cutoff when its pooled coverage reaches the cutoff in
#' both genotypes.
coverage_threshold_sweep <- function(intervals, thresholds) {
  rows <- lapply(thresholds, function(t) {
    pass <- intervals$ctrl_coverage >= t & intervals$mut_coverage >= t
    per_type <- data.frame(
      min_coverage = t,
      feature_type = levels(intervals$feature_type),
      n_intervals  = as.integer(table(intervals$feature_type)),
      n_pass       = as.integer(table(intervals$feature_type[pass])),
      stringsAsFactors = FALSE
    )
    overall <- data.frame(
      min_coverage = t,
      feature_type = "All",
      n_intervals  = nrow(intervals),
      n_pass       = sum(pass),
      stringsAsFactors = FALSE
    )
    rbind(per_type, overall)
  })
  out <- do.call(rbind, rows)
  out$n_removed <- out$n_intervals - out$n_pass
  out$pct_removed <- 100 * out$n_removed / out$n_intervals
  out
}

#' Keep the intervals whose pooled coverage reaches min_coverage in both
#' genotypes, and log what each feature type loses.
apply_coverage_filter <- function(intervals, min_coverage) {
  keep <- !is.na(intervals$mch_ctrl) & !is.na(intervals$mch_mut) &
    intervals$ctrl_coverage >= min_coverage & intervals$mut_coverage >= min_coverage

  cat(sprintf("\nApplying the coverage filter (>= %g pooled reads per genotype):\n",
              min_coverage))
  for (ft in levels(intervals$feature_type)) {
    in_type <- intervals$feature_type == ft
    n_total <- sum(in_type)
    n_kept <- sum(in_type & keep)
    cat(sprintf("  %-20s %10s of %10s kept (%5.1f%% removed)\n", ft,
                fmt_int(n_kept), fmt_int(n_total),
                100 * (n_total - n_kept) / n_total))
  }
  cat(sprintf("  %-20s %10s of %10s kept (%5.1f%% removed)\n", "All",
              fmt_int(sum(keep)), fmt_int(nrow(intervals)),
              100 * (nrow(intervals) - sum(keep)) / nrow(intervals)))

  kept <- intervals[keep]
  if (nrow(kept) == 0) {
    stop("No interval reaches ", min_coverage,
         " pooled reads in both genotypes. Lower --min-coverage.")
  }
  empty_types <- levels(kept$feature_type)[table(kept$feature_type) == 0]
  if (length(empty_types) > 0) {
    stop("No interval survives the coverage filter for: ",
         paste(empty_types, collapse = ", "),
         ". Lower --min-coverage.")
  }
  kept
}

# =============================================================================
# ENHANCER OVERLAP
# =============================================================================

#' Flag the intervals that overlap a chromHMM active enhancer segment.
add_enhancer_overlap <- function(intervals) {
  cat("\nLoading the active enhancer segmentation...\n")
  enhancers <- load_chip_peaks(CHROMATIN_PATHS$active_enhancer, "Active enhancer")

  gr <- GRanges(
    seqnames = intervals$chr,
    ranges = IRanges(start = intervals$start + 1L, end = intervals$end)
  )
  overlap <- countOverlaps(gr, enhancers) > 0

  if (sum(overlap) == 0) {
    stop("No feature interval overlaps ", CHROMATIN_PATHS$active_enhancer,
         ".\n  Check that the enhancer BED uses the same chromosome names as the",
         " feature BEDs (chr1 style).")
  }

  intervals[, enhancer_overlap := overlap]
  intervals[, enhancer_status := factor(
    ifelse(overlap, ENHANCER_LEVELS[1], ENHANCER_LEVELS[2]),
    levels = ENHANCER_LEVELS)]

  cat(sprintf("  %s of %s intervals overlap an active enhancer (%.1f%%)\n",
              fmt_int(sum(overlap)), fmt_int(nrow(intervals)),
              100 * mean(overlap)))
  intervals[]
}

# =============================================================================
# SUMMARY TABLES
# =============================================================================

#' Interval counts, bases, and pooled rates for every feature type.
feature_type_summary <- function(all_intervals, kept_intervals) {
  total <- all_intervals[, .(
    n_intervals_total = .N,
    total_bases_all = sum(as.numeric(feature_length))
  ), by = feature_type]

  kept <- kept_intervals[, .(
    n_intervals_pass = .N,
    n_genes = uniqueN(gene_name),
    total_bases_pass = sum(as.numeric(feature_length)),
    median_interval_length = median(feature_length),
    total_ch_sites = sum(ctrl_ch_sites + mut_ch_sites),
    pooled_mch_ctrl = sum(ctrl_methylated) / sum(ctrl_coverage),
    pooled_mch_mut = sum(mut_methylated) / sum(mut_coverage),
    mean_delta_mch = mean(delta_mch),
    median_delta_mch = median(delta_mch),
    n_enhancer_intervals = sum(enhancer_overlap)
  ), by = feature_type]

  out <- merge(total, kept, by = "feature_type", all = TRUE)
  out[, pooled_delta_mch := pooled_mch_mut - pooled_mch_ctrl]
  out[, pct_intervals_pass := 100 * n_intervals_pass / n_intervals_total]
  out[, pct_enhancer := 100 * n_enhancer_intervals / n_intervals_pass]
  setorder(out, feature_type)
  out[]
}

#' Pool the intervals of one feature type within one gene.
gene_feature_summary <- function(intervals) {
  gf <- intervals[, .(
    chr = chr[1],
    strand = strand[1],
    n_intervals = .N,
    total_bases = sum(as.numeric(feature_length)),
    n_enhancer_intervals = sum(enhancer_overlap),
    ctrl_coverage = sum(ctrl_coverage),
    ctrl_methylated = sum(ctrl_methylated),
    mut_coverage = sum(mut_coverage),
    mut_methylated = sum(mut_methylated)
  ), by = .(gene_name, feature_type)]

  gf[, mch_ctrl := ctrl_methylated / ctrl_coverage]
  gf[, mch_mut := mut_methylated / mut_coverage]
  gf[, delta_mch := mch_mut - mch_ctrl]
  gf[, has_enhancer := n_enhancer_intervals > 0]
  setorder(gf, gene_name, feature_type)
  gf[]
}

#' Gene-level mCH results with one row per gene symbol.
#'
#' mch_results carries a few gene symbols under more than one ENSMUSG id. The
#' row with the smallest FDR is kept, and the count of dropped rows is logged.
gene_level_mch <- function() {
  df <- mch_results[order(mch_results$edger_fdr,
                          -abs(mch_results$edger_logFC)), ]
  n_before <- nrow(df)
  df <- df[!duplicated(df$gene_name), ]
  cat(sprintf("\nGene-level mCH results: %s rows, %s unique gene symbols ",
              fmt_int(n_before), fmt_int(nrow(df))))
  cat(sprintf("(%s duplicate symbol rows dropped)\n", fmt_int(n_before - nrow(df))))

  data.frame(
    gene_name        = df$gene_name,
    gene_edger_logFC = df$edger_logFC,
    gene_edger_fdr   = df$edger_fdr,
    gene_mch_diff    = df$mch_diff,
    gene_mch_sig     = df$mch_sig,
    gene_mch_hyper   = df$mch_hyper,
    gene_mch_hypo    = df$mch_hypo,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# STATISTICS
# =============================================================================

#' Kruskal-Wallis test of delta_mch across the six feature types.
run_kruskal_wallis <- function(intervals) {
  kw <- kruskal.test(intervals$delta_mch, intervals$feature_type)
  chi_squared <- unname(kw$statistic)
  df <- as.integer(unname(kw$parameter))

  cat("\nKruskal-Wallis (delta_mch across feature types): ")
  cat(sprintf("chi-squared = %.2f, df = %d, p = %s\n",
              chi_squared, df, format.pval(kw$p.value, digits = 3)))

  list(
    test = kw,
    table = data.frame(
      test = "Kruskal-Wallis of delta_mch across feature types",
      n_intervals = nrow(intervals),
      n_groups = nlevels(intervals$feature_type),
      chi_squared = chi_squared,
      df = df,
      p_value = kw$p.value,
      stringsAsFactors = FALSE
    )
  )
}

#' Dunn's post-hoc pairwise comparison of delta_mch between feature types.
#'
#' The Z, P, and P.adjusted values are taken as dunn.test() returns them, with
#' the Benjamini-Hochberg adjustment requested by method = "bh".
run_dunn_posthoc <- function(intervals) {
  cat("\nDunn's post-hoc pairwise tests (BH-adjusted):\n")
  res <- dunn.test(intervals$delta_mch, as.character(intervals$feature_type),
                   method = "bh", kw = FALSE, table = FALSE, list = TRUE)

  out <- data.frame(
    comparison = res$comparisons,
    Z = res$Z,
    p_dunn = res$P,
    p_dunn_bh = res$P.adjusted,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$p_dunn_bh), ]
  n_sig <- sum(out$p_dunn_bh < Q_THRESHOLD)
  cat(sprintf("  %d of %d pairs reach BH q < %.2f\n", n_sig, nrow(out),
              Q_THRESHOLD))
  for (i in seq_len(nrow(out))) {
    cat(sprintf("    %-45s Z = %8.3f  q = %.3g\n",
                out$comparison[i], out$Z[i], out$p_dunn_bh[i]))
  }
  out
}

#' Wilcoxon signed rank test of delta_mch against zero inside each feature type.
run_wilcoxon_vs_zero <- function(intervals) {
  cat("\nWilcoxon signed rank tests of delta_mch against zero:\n")
  rows <- lapply(levels(intervals$feature_type), function(ft) {
    x <- intervals$delta_mch[intervals$feature_type == ft]
    x <- x[!is.na(x)]
    if (length(x) < 10) {
      stop("Feature type ", ft, " has only ", length(x),
           " intervals with a delta. Lower --min-coverage.")
    }
    wt <- wilcox.test(x, mu = 0, exact = FALSE, correct = TRUE)
    data.frame(
      feature_type = ft,
      n_intervals = length(x),
      median_delta_mch = median(x),
      mean_delta_mch = mean(x),
      pct_negative = 100 * mean(x < 0),
      statistic_V = unname(wt$statistic),
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$p_adj_bh <- p.adjust(out$p_value, method = "BH")
  out$feature_type <- factor(out$feature_type, levels = FEATURE_TYPES)
  out <- out[order(out$feature_type), ]

  for (i in seq_len(nrow(out))) {
    cat(sprintf("  %-20s n = %9s  median = %+.5f  q = %s\n",
                as.character(out$feature_type[i]), fmt_int(out$n_intervals[i]),
                out$median_delta_mch[i], format.pval(out$p_adj_bh[i], digits = 3)))
  }
  out
}

#' Wilcoxon rank sum test of delta_mch between enhancer groups, per feature type.
run_enhancer_wilcoxon <- function(intervals) {
  cat("\nWilcoxon rank sum tests, active enhancer against no enhancer:\n")
  rows <- lapply(levels(intervals$feature_type), function(ft) {
    sub <- intervals[intervals$feature_type == ft]
    marked <- sub$delta_mch[sub$enhancer_overlap]
    unmarked <- sub$delta_mch[!sub$enhancer_overlap]

    if (length(marked) == 0 || length(unmarked) == 0) {
      stop("Feature type ", ft, " has ", length(marked),
           " enhancer-overlapping and ", length(unmarked),
           " unmarked intervals, so the two groups cannot be compared.",
           "\n  Check that ", CHROMATIN_PATHS$active_enhancer,
           " covers the same chromosomes as the feature BEDs.")
    }

    wt <- wilcox.test(marked, unmarked, exact = FALSE, correct = TRUE)
    data.frame(
      feature_type = ft,
      n_enhancer = length(marked),
      n_no_enhancer = length(unmarked),
      median_delta_enhancer = median(marked),
      median_delta_no_enhancer = median(unmarked),
      median_difference = median(marked) - median(unmarked),
      statistic_W = unname(wt$statistic),
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$p_adj_bh <- p.adjust(out$p_value, method = "BH")
  out$feature_type <- factor(out$feature_type, levels = FEATURE_TYPES)
  out <- out[order(out$feature_type), ]

  for (i in seq_len(nrow(out))) {
    cat(sprintf("  %-20s enhancer median = %+.5f (n = %s)  unmarked median = %+.5f (n = %s)  q = %s\n",
                as.character(out$feature_type[i]),
                out$median_delta_enhancer[i], fmt_int(out$n_enhancer[i]),
                out$median_delta_no_enhancer[i], fmt_int(out$n_no_enhancer[i]),
                format.pval(out$p_adj_bh[i], digits = 3)))
  }
  out
}

#' Interval counts and pooled rates for each feature type and enhancer group.
enhancer_overlap_summary <- function(intervals) {
  out <- intervals[, .(
    n_intervals = .N,
    n_genes = uniqueN(gene_name),
    total_bases = sum(as.numeric(feature_length)),
    pooled_mch_ctrl = sum(ctrl_methylated) / sum(ctrl_coverage),
    pooled_mch_mut = sum(mut_methylated) / sum(mut_coverage),
    mean_delta_mch = mean(delta_mch),
    median_delta_mch = median(delta_mch)
  ), by = .(feature_type, enhancer_status)]
  out[, pooled_delta_mch := pooled_mch_mut - pooled_mch_ctrl]
  setorder(out, feature_type, enhancer_status)
  out[]
}

# =============================================================================
# RANK AND METAGENE PROFILES
# =============================================================================

#' mCH by exon and by intron rank position along the direction of transcription.
exon_intron_rank_profile <- function(intervals, max_rank) {
  sub <- intervals[feature_type %in% BODY_FEATURES & feature_rank <= max_rank]
  if (nrow(sub) == 0) {
    stop("No exon or intron interval has a rank at or below ", max_rank, ".")
  }

  out <- sub[, .(
    n_intervals = .N,
    n_genes = uniqueN(gene_name),
    median_interval_length = median(feature_length),
    pooled_mch_ctrl = sum(ctrl_methylated) / sum(ctrl_coverage),
    pooled_mch_mut = sum(mut_methylated) / sum(mut_coverage),
    mean_delta_mch = mean(delta_mch),
    median_delta_mch = median(delta_mch),
    se_delta_mch = sd(delta_mch) / sqrt(.N)
  ), by = .(feature_type, feature_rank)]
  out[, pooled_delta_mch := pooled_mch_mut - pooled_mch_ctrl]
  out[, feature_type := droplevels(feature_type)]
  setorder(out, feature_type, feature_rank)
  out[]
}

#' mCH across the scaled gene body, built from the exon and intron intervals.
#'
#' Each gene body runs from the first to the last base covered by its exon and
#' intron intervals. Every interval is placed on a 0 to 1 axis oriented 5' to
#' 3', then split between the bins it spans in proportion to the bases it shares
#' with each bin. Counts are pooled per gene and bin, and bins that reach the
#' coverage cutoff in both genotypes contribute one rate each.
metagene_profile <- function(intervals, n_bins, min_coverage) {
  body <- intervals[feature_type %in% BODY_FEATURES]
  if (nrow(body) == 0) {
    stop("No exon or intron interval survives the coverage filter, so the ",
         "metagene profile cannot be built.")
  }

  spans <- body[, .(
    body_start = min(start),
    body_end = max(end),
    n_strand = uniqueN(strand),
    gene_strand = strand[1]
  ), by = gene_name]

  bad_strand <- spans[n_strand > 1]
  if (nrow(bad_strand) > 0) {
    stop("Feature intervals disagree on strand for ", nrow(bad_strand),
         " genes, for example: ",
         paste(head(bad_strand$gene_name, 5), collapse = ", "))
  }
  spans <- spans[body_end > body_start]

  body <- merge(body, spans[, .(gene_name, body_start, body_end, gene_strand)],
                by = "gene_name")
  body[, body_length := body_end - body_start]
  body[, raw_start := (start - body_start) / body_length]
  body[, raw_end := (end - body_start) / body_length]
  body[, oriented_start := ifelse(gene_strand == "-", 1 - raw_end, raw_start)]
  body[, oriented_end := ifelse(gene_strand == "-", 1 - raw_start, raw_end)]
  body[, rel_start := pmax(0, pmin(1, oriented_start))]
  body[, rel_end := pmax(0, pmin(1, oriented_end))]
  body[, span := rel_end - rel_start]
  body <- body[span > 0]

  body[, bin_first := as.integer(pmin(n_bins, floor(rel_start * n_bins) + 1))]
  body[, bin_last := as.integer(pmax(1, pmin(n_bins, ceiling(rel_end * n_bins))))]
  body[, n_bins_spanned := bin_last - bin_first + 1L]

  # One row per interval and bin it reaches. sequence() counts 0-based inside
  # each interval's block of repeated rows.
  slim <- body[, .(gene_name, rel_start, rel_end, span, bin_first,
                   n_bins_spanned, ctrl_coverage, ctrl_methylated,
                   mut_coverage, mut_methylated)]
  expanded <- slim[rep(seq_len(nrow(slim)), slim$n_bins_spanned)]
  expanded[, bin := bin_first + sequence(slim$n_bins_spanned) - 1L]
  expanded[, bin_low := (bin - 1) / n_bins]
  expanded[, bin_high := bin / n_bins]
  expanded[, share := pmax(0, pmin(rel_end, bin_high) - pmax(rel_start, bin_low)) /
             span]

  gene_bin <- expanded[, .(
    ctrl_coverage = sum(ctrl_coverage * share),
    ctrl_methylated = sum(ctrl_methylated * share),
    mut_coverage = sum(mut_coverage * share),
    mut_methylated = sum(mut_methylated * share)
  ), by = .(gene_name, bin)]

  gene_bin <- gene_bin[ctrl_coverage >= min_coverage & mut_coverage >= min_coverage]
  if (nrow(gene_bin) == 0) {
    stop("No gene and bin combination reaches ", min_coverage,
         " pooled reads in both genotypes.")
  }
  gene_bin[, mch_ctrl := ctrl_methylated / ctrl_coverage]
  gene_bin[, mch_mut := mut_methylated / mut_coverage]
  gene_bin[, delta_mch := mch_mut - mch_ctrl]

  profile <- gene_bin[, .(
    n_genes = .N,
    mean_mch_ctrl = mean(mch_ctrl),
    se_mch_ctrl = sd(mch_ctrl) / sqrt(.N),
    mean_mch_mut = mean(mch_mut),
    se_mch_mut = sd(mch_mut) / sqrt(.N),
    mean_delta_mch = mean(delta_mch),
    median_delta_mch = median(delta_mch),
    se_delta_mch = sd(delta_mch) / sqrt(.N),
    pooled_mch_ctrl = sum(ctrl_methylated) / sum(ctrl_coverage),
    pooled_mch_mut = sum(mut_methylated) / sum(mut_coverage)
  ), by = bin]
  profile[, pooled_delta_mch := pooled_mch_mut - pooled_mch_ctrl]
  setorder(profile, bin)

  cat(sprintf("\nMetagene profile: %d bins, %s genes, mean %s genes per bin\n",
              n_bins, fmt_int(uniqueN(gene_bin$gene_name)),
              fmt_int(round(mean(profile$n_genes)))))
  list(profile = profile[], n_genes = uniqueN(gene_bin$gene_name))
}

# =============================================================================
# GENE-LEVEL FISHER TESTS
# =============================================================================

#' Build the gene table behind the registered Fisher tests.
#'
#' One row per gene that has both exon and intron intervals after the coverage
#' filter, with the gene-level mCH call joined on the gene symbol.
build_fisher_gene_table <- function(gene_feature, gene_mch) {
  exon <- gene_feature[feature_type == "Exon",
                       .(gene_name, chr,
                         exon_delta_mch = delta_mch,
                         exon_has_enhancer = has_enhancer)]
  intron <- gene_feature[feature_type == "Intron",
                         .(gene_name,
                           intron_delta_mch = delta_mch,
                           intron_has_enhancer = has_enhancer)]

  df <- merge(exon, intron, by = "gene_name")
  df <- merge(df, as.data.table(gene_mch), by = "gene_name", all.x = TRUE)

  df[, exon_hypo := exon_delta_mch < 0]
  df[, intron_hypo := intron_delta_mch < 0]
  df <- as.data.frame(df)

  cat(sprintf("\nFisher gene table: %s genes with exon and intron data, ",
              fmt_int(nrow(df))))
  cat(sprintf("%s of them in the gene-level mCH results\n",
              fmt_int(sum(!is.na(df$gene_mch_sig)))))
  df
}

#' Run and register the three gene-level 2x2 tests.
run_gene_fisher_tests <- function(fisher_genes, out_dir) {
  specs <- list(
    list(test_id = "exon_intron_hypo_concordance",
         row_var = "exon_hypo", col_var = "intron_hypo",
         description = paste("Genes losing mCH across their exons against genes",
                             "losing mCH across their introns.")),
    list(test_id = "intron_enhancer_vs_gene_sig",
         row_var = "intron_has_enhancer", col_var = "gene_mch_sig",
         description = paste("Genes with an intron over an active enhancer",
                             "against genes with significant gene-body mCH change.")),
    list(test_id = "intron_hypo_vs_gene_hypo",
         row_var = "intron_hypo", col_var = "gene_mch_hypo",
         description = paste("Genes losing mCH across their introns against genes",
                             "called significantly hypomethylated at the gene level."))
  )

  rows <- lapply(specs, function(spec) {
    ft <- register_fisher_test(
      section = SECTION_ID,
      test_id = spec$test_id,
      description = spec$description,
      gene_df = fisher_genes,
      row_var = spec$row_var,
      col_var = spec$col_var,
      output_dir = out_dir
    )

    keep <- !is.na(fisher_genes[[spec$row_var]]) &
      !is.na(fisher_genes[[spec$col_var]])
    row_true <- fisher_genes[[spec$row_var]][keep]
    col_true <- fisher_genes[[spec$col_var]][keep]

    data.frame(
      test_id = spec$test_id,
      description = spec$description,
      row_var = spec$row_var,
      col_var = spec$col_var,
      n_genes = sum(keep),
      n_row_true = sum(row_true),
      n_col_true = sum(col_true),
      n_both_true = sum(row_true & col_true),
      odds_ratio = unname(ft$estimate),
      ci_low = ft$conf.int[1],
      ci_high = ft$conf.int[2],
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# =============================================================================
# FIGURES
# =============================================================================

#' Axis limits that hide the extreme tails of a distribution.
clip_limits <- function(values, lower = 0.005, upper = 0.995) {
  unname(quantile(values, c(lower, upper), na.rm = TRUE))
}

#' Interval counts, bases, and coverage pass rate per feature type.
plot_feature_overview <- function(summary_df, min_coverage) {
  df <- as.data.frame(summary_df)

  p_count <- ggplot(df, aes(x = feature_type, y = n_intervals_pass,
                            fill = feature_type)) +
    geom_col(color = "black", linewidth = 0.3) +
    geom_text(aes(label = fmt_int(n_intervals_pass)), vjust = -0.35, size = 2.9) +
    scale_fill_manual(values = FEATURE_COLORS, drop = FALSE) +
    scale_y_continuous(limits = c(0, max(df$n_intervals_pass) * 1.18),
                       expand = c(0, 0), labels = scales::comma) +
    labs(title = "Intervals used", x = "Feature type", y = "Intervals") +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "none")

  p_bases <- ggplot(df, aes(x = feature_type, y = total_bases_pass / 1e6,
                            fill = feature_type)) +
    geom_col(color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f", total_bases_pass / 1e6)),
              vjust = -0.35, size = 2.9) +
    scale_fill_manual(values = FEATURE_COLORS, drop = FALSE) +
    scale_y_continuous(limits = c(0, max(df$total_bases_pass / 1e6) * 1.18),
                       expand = c(0, 0)) +
    labs(title = "Bases covered", x = "Feature type", y = "Total bases (Mb)") +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "none")

  p_pass <- ggplot(df, aes(x = feature_type, y = pct_intervals_pass,
                           fill = feature_type)) +
    geom_col(color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%\n%s / %s", pct_intervals_pass,
                                  fmt_int(n_intervals_pass),
                                  fmt_int(n_intervals_total))),
              vjust = -0.15, size = 2.6, lineheight = 0.9) +
    scale_fill_manual(values = FEATURE_COLORS, drop = FALSE) +
    scale_y_continuous(limits = c(0, 118), expand = c(0, 0)) +
    labs(title = "Intervals passing the coverage filter",
         x = "Feature type", y = "Intervals kept (%)") +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "none")

  p_count + p_bases + p_pass +
    plot_layout(ncol = 3) +
    plot_annotation(
      title = "Sub-gene feature intervals used for the mCH comparison",
      subtitle = sprintf("Coverage filter: at least %g pooled reads per genotype per interval",
                         min_coverage),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = 15),
                    plot.subtitle = element_text(hjust = 0.5, size = 11))
    )
}

#' Absolute mCH per feature type, control beside mutant.
plot_mch_level_by_feature <- function(level_df, level_stats) {
  limits <- clip_limits(level_df$mch_pct)
  label_y <- limits[2] + 0.05 * diff(limits)
  dodge <- position_dodge(width = 0.8)

  ggplot(level_df, aes(x = feature_type, y = mch_pct, fill = genotype)) +
    geom_boxplot(outlier.shape = NA, width = 0.68, position = dodge,
                 color = "grey20", linewidth = 0.3) +
    geom_text(data = level_stats,
              aes(x = feature_type, y = label_y, label = label, group = genotype),
              position = dodge, inherit.aes = FALSE, size = 2.4, vjust = 0,
              lineheight = 0.9) +
    scale_fill_manual(values = COLORS$genotype, labels = GENOTYPE_LABELS,
                      name = "Genotype") +
    coord_cartesian(ylim = c(limits[1], limits[2] + 0.32 * diff(limits))) +
    labs(
      title = "mCH level by sub-gene feature type",
      subtitle = paste("Read counts pooled over the four samples of each genotype.",
                       "Axis clipped to the 0.5 and 99.5 percentiles;",
                       "outliers are not drawn."),
      x = "Feature type", y = "mCH (%)"
    ) +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "top")
}

#' delta_mch by feature type, carrying the Kruskal-Wallis result.
plot_delta_by_feature <- function(plot_df, delta_stats, kw_table) {
  limits <- clip_limits(plot_df$delta_mch_pct)
  label_y <- limits[2] + 0.05 * diff(limits)

  ggplot(plot_df, aes(x = feature_type, y = delta_mch_pct, fill = feature_type)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_boxplot(outlier.shape = NA, width = 0.7, color = "grey20",
                 linewidth = 0.3) +
    geom_text(data = delta_stats,
              aes(x = feature_type, y = label_y, label = label),
              inherit.aes = FALSE, size = 2.7, vjust = 0, lineheight = 0.9) +
    scale_fill_manual(values = FEATURE_COLORS, drop = FALSE) +
    coord_cartesian(ylim = c(limits[1], limits[2] + 0.30 * diff(limits))) +
    labs(
      title = "mCH change by sub-gene feature type",
      subtitle = sprintf(paste("Mutant minus control, in percentage points.",
                               "Kruskal-Wallis chi-squared = %.1f, df = %d, p = %s.",
                               "\nAxis clipped to the 0.5 and 99.5 percentiles;",
                               "outliers are not drawn."),
                         kw_table$chi_squared, kw_table$df,
                         format.pval(kw_table$p_value, digits = 3)),
      x = "Feature type", y = "delta mCH (percentage points)"
    ) +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "none")
}

#' Square matrix of BH-adjusted Dunn q-values, one cell per feature pair.
build_dunn_matrix <- function(posthoc, feature_levels) {
  n <- length(feature_levels)
  mat <- matrix(NA_real_, nrow = n, ncol = n,
                dimnames = list(feature_levels, feature_levels))
  for (i in seq_len(nrow(posthoc))) {
    parts <- trimws(strsplit(posthoc$comparison[i], "-", fixed = TRUE)[[1]])
    if (length(parts) == 2 && all(parts %in% feature_levels)) {
      mat[parts[1], parts[2]] <- posthoc$p_dunn_bh[i]
      mat[parts[2], parts[1]] <- posthoc$p_dunn_bh[i]
    }
  }
  if (all(is.na(mat[upper.tri(mat)]))) {
    stop("No Dunn comparison label matched the feature type names. ",
         "Labels seen: ", paste(head(posthoc$comparison, 3), collapse = "; "))
  }
  mat
}

#' Lower-triangle heatmap of the Dunn post-hoc q-values.
plot_dunn_heatmap <- function(mat) {
  df <- expand.grid(feature_row = rownames(mat), feature_col = colnames(mat),
                    stringsAsFactors = FALSE)
  df$q_value <- mapply(function(r, c) mat[r, c], df$feature_row, df$feature_col)
  keep <- match(df$feature_row, rownames(mat)) > match(df$feature_col, colnames(mat))
  df <- df[keep & !is.na(df$q_value), ]

  df$neg_log10_q <- -log10(pmax(df$q_value, 1e-300))
  df$stars <- ifelse(df$q_value < 0.001, "***",
                     ifelse(df$q_value < 0.01, "**",
                            ifelse(df$q_value < 0.05, "*", "ns")))
  df$display <- sprintf("%.2g\n%s", df$q_value, df$stars)
  df$feature_row <- factor(df$feature_row, levels = rev(rownames(mat)))
  df$feature_col <- factor(df$feature_col, levels = colnames(mat))

  ggplot(df, aes(x = feature_col, y = feature_row, fill = neg_log10_q)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = display), size = 2.8, lineheight = 0.85) +
    scale_fill_gradient2(low = "grey90", mid = "#FDB863", high = "#B2182B",
                         midpoint = -log10(Q_THRESHOLD), limits = c(0, NA),
                         name = expression(-log[10] * "(q)")) +
    labs(
      title = "Dunn's post-hoc comparison of delta mCH between feature types",
      subtitle = "Benjamini-Hochberg q-values. * q < 0.05, ** q < 0.01, *** q < 0.001, ns not significant",
      x = NULL, y = NULL
    ) +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          panel.grid = element_blank())
}

#' Per-interval delta mCH along each key gene, drawn as lollipops.
plot_key_gene_loci <- function(key_intervals, n_key_genes) {
  ggplot(key_intervals, aes(x = midpoint_mb, y = delta_mch * 100,
                            color = feature_type)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.3) +
    geom_segment(aes(xend = midpoint_mb, yend = 0), linewidth = 0.4) +
    geom_point(size = 1.4) +
    facet_wrap(~ facet_label, scales = "free_x", ncol = 2) +
    scale_color_manual(values = FEATURE_COLORS, drop = FALSE, name = "Feature type") +
    labs(
      title = "mCH change along the key genes",
      subtitle = sprintf(paste("%d of %d key genes have intervals passing the",
                               "coverage filter. One point per interval,",
                               "mutant minus control."),
                         n_key_genes, length(KEY_GENES)),
      x = "Genomic position (Mb)", y = "delta mCH (percentage points)"
    ) +
    theme_emseq() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 8, face = "bold"))
}

#' mCH and delta mCH across the scaled gene body.
plot_metagene <- function(profile, n_bins, n_genes) {
  rate_long <- rbind(
    data.frame(bin = profile$bin, genotype = "ctrl",
               mean_rate = profile$mean_mch_ctrl * 100,
               se_rate = profile$se_mch_ctrl * 100),
    data.frame(bin = profile$bin, genotype = "mut",
               mean_rate = profile$mean_mch_mut * 100,
               se_rate = profile$se_mch_mut * 100)
  )

  break_positions <- c("5'" = 1, "25%" = round(n_bins / 4),
                       "50%" = round(n_bins / 2), "75%" = round(3 * n_bins / 4),
                       "3'" = n_bins)
  distinct <- !duplicated(break_positions)
  breaks <- unname(break_positions[distinct])
  labels <- names(break_positions)[distinct]

  p_rate <- ggplot(rate_long, aes(x = bin, y = mean_rate, color = genotype,
                                  fill = genotype)) +
    geom_ribbon(aes(ymin = mean_rate - se_rate, ymax = mean_rate + se_rate),
                alpha = 0.25, color = NA) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = COLORS$genotype, labels = GENOTYPE_LABELS,
                       name = "Genotype") +
    scale_fill_manual(values = COLORS$genotype, labels = GENOTYPE_LABELS,
                      name = "Genotype") +
    scale_x_continuous(breaks = breaks, labels = labels) +
    labs(title = "mCH across the scaled gene body",
         x = NULL, y = "mCH (%), mean and standard error") +
    theme_emseq() +
    theme(legend.position = "top")

  p_delta <- ggplot(profile, aes(x = bin, y = mean_delta_mch * 100)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_ribbon(aes(ymin = (mean_delta_mch - se_delta_mch) * 100,
                    ymax = (mean_delta_mch + se_delta_mch) * 100),
                alpha = 0.25, fill = "grey40") +
    geom_line(linewidth = 0.9, color = "grey15") +
    scale_x_continuous(breaks = breaks, labels = labels) +
    labs(title = "mCH change across the scaled gene body",
         x = "Relative position in the gene body (5' to 3')",
         y = "delta mCH (percentage points)") +
    theme_emseq()

  p_rate / p_delta +
    plot_annotation(
      title = "Metagene mCH profile",
      subtitle = sprintf(paste("%s genes, %d bins. Built from exon and intron",
                               "intervals, each split between the bins it spans",
                               "in proportion to shared bases."),
                         fmt_int(n_genes), n_bins),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = 15),
                    plot.subtitle = element_text(hjust = 0.5, size = 11))
    )
}

#' mCH and delta mCH by exon and intron rank position.
plot_rank_profile <- function(rank_profile, max_rank) {
  rate_long <- rbind(
    data.frame(feature_type = rank_profile$feature_type,
               feature_rank = rank_profile$feature_rank,
               genotype = "ctrl", rate = rank_profile$pooled_mch_ctrl * 100),
    data.frame(feature_type = rank_profile$feature_type,
               feature_rank = rank_profile$feature_rank,
               genotype = "mut", rate = rank_profile$pooled_mch_mut * 100)
  )

  p_rate <- ggplot(rate_long, aes(x = feature_rank, y = rate, color = genotype)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.4) +
    facet_wrap(~ feature_type, ncol = 2) +
    scale_color_manual(values = COLORS$genotype, labels = GENOTYPE_LABELS,
                       name = "Genotype") +
    scale_x_continuous(breaks = pretty(seq_len(max_rank))) +
    labs(title = "mCH by rank position",
         x = NULL, y = "Pooled mCH (%)") +
    theme_emseq() +
    theme(legend.position = "top")

  p_delta <- ggplot(rank_profile, aes(x = feature_rank, y = mean_delta_mch * 100)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_ribbon(aes(ymin = (mean_delta_mch - se_delta_mch) * 100,
                    ymax = (mean_delta_mch + se_delta_mch) * 100),
                alpha = 0.25, fill = "grey40") +
    geom_line(linewidth = 0.8, color = "grey15") +
    geom_point(size = 1.4, color = "grey15") +
    facet_wrap(~ feature_type, ncol = 2) +
    scale_x_continuous(breaks = pretty(seq_len(max_rank))) +
    labs(title = "mCH change by rank position",
         x = "Rank along the direction of transcription (1 is most 5')",
         y = "delta mCH (percentage points)") +
    theme_emseq()

  p_rate / p_delta +
    plot_annotation(
      title = "Exon and intron rank position profiles",
      subtitle = sprintf(paste("Ranks 1 to %d. Rates pool read counts over the",
                               "four samples of each genotype; the change line",
                               "shows the mean and standard error over intervals."),
                         max_rank),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = 15),
                    plot.subtitle = element_text(hjust = 0.5, size = 11))
    )
}

#' delta_mch split by active enhancer overlap, within each feature type.
plot_enhancer_split <- function(plot_df, enhancer_stats, wilcoxon) {
  limits <- clip_limits(plot_df$delta_mch_pct)
  label_y <- limits[2] + 0.05 * diff(limits)
  dodge <- position_dodge(width = 0.8)

  star_df <- wilcoxon
  star_df$stars <- ifelse(star_df$p_adj_bh < 0.001, "***",
                          ifelse(star_df$p_adj_bh < 0.01, "**",
                                 ifelse(star_df$p_adj_bh < 0.05, "*", "ns")))
  star_df$y <- limits[2] + 0.28 * diff(limits)

  ggplot(plot_df, aes(x = feature_type, y = delta_mch_pct,
                      fill = enhancer_status)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_boxplot(outlier.shape = NA, width = 0.68, position = dodge,
                 color = "grey20", linewidth = 0.3) +
    geom_text(data = enhancer_stats,
              aes(x = feature_type, y = label_y, label = label,
                  group = enhancer_status),
              position = dodge, inherit.aes = FALSE, size = 2.3, vjust = 0,
              lineheight = 0.9) +
    geom_text(data = star_df, aes(x = feature_type, y = y, label = stars),
              inherit.aes = FALSE, size = 3.4, vjust = 0) +
    scale_fill_manual(values = ENHANCER_COLORS, name = "Interval overlaps") +
    coord_cartesian(ylim = c(limits[1], limits[2] + 0.42 * diff(limits))) +
    labs(
      title = "mCH change by active enhancer overlap",
      subtitle = paste("Wilcoxon rank sum test per feature type, BH-adjusted:",
                       "* q < 0.05, ** q < 0.01, *** q < 0.001, ns not significant.",
                       "\nAxis clipped to the 0.5 and 99.5 percentiles;",
                       "outliers are not drawn."),
      x = "Feature type", y = "delta mCH (percentage points)"
    ) +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "top")
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opt <- parse_options()
  out_dir <- opt$output_dir

  cat("================================================================================\n")
  cat("SECTION 50_01: mCH ACROSS SUB-GENE FEATURE TYPES\n")
  cat("================================================================================\n")
  cat("Output dir:    ", out_dir, "\n", sep = "")
  cat("Feature dir:   ", DATA_PATHS$feature_dir, "\n", sep = "")
  cat("Min coverage:  ", opt$min_coverage, " pooled reads per genotype\n", sep = "")
  cat("Max rank:      ", opt$max_rank, "\n", sep = "")
  cat("Metagene bins: ", opt$metagene_bins, "\n\n", sep = "")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  check_inputs()

  all_intervals <- load_interval_counts()

  cat("\nCoverage threshold sweep:\n")
  sweep_thresholds <- sort(unique(c(COVERAGE_SWEEP_THRESHOLDS, opt$min_coverage)))
  sweep <- coverage_threshold_sweep(all_intervals, sweep_thresholds)
  sweep_all <- sweep[sweep$feature_type == "All", ]
  for (i in seq_len(nrow(sweep_all))) {
    cat(sprintf("  >= %6g reads: %10s kept, %10s removed (%5.1f%%)\n",
                sweep_all$min_coverage[i], fmt_int(sweep_all$n_pass[i]),
                fmt_int(sweep_all$n_removed[i]), sweep_all$pct_removed[i]))
  }

  intervals <- apply_coverage_filter(all_intervals, opt$min_coverage)
  intervals <- add_enhancer_overlap(intervals)

  cat("\nBuilding summary tables...\n")
  type_summary <- feature_type_summary(all_intervals, intervals)
  gene_feature <- gene_feature_summary(intervals)
  gene_mch <- gene_level_mch()

  plot_df <- data.frame(
    feature_type = intervals$feature_type,
    delta_mch_pct = intervals$delta_mch * 100,
    enhancer_status = intervals$enhancer_status,
    stringsAsFactors = FALSE
  )
  delta_stats <- summarise_groups(plot_df, "feature_type", "delta_mch_pct")

  level_df <- rbind(
    data.frame(feature_type = intervals$feature_type, genotype = "ctrl",
               mch_pct = intervals$mch_ctrl * 100, stringsAsFactors = FALSE),
    data.frame(feature_type = intervals$feature_type, genotype = "mut",
               mch_pct = intervals$mch_mut * 100, stringsAsFactors = FALSE)
  )
  level_stats <- do.call(rbind, lapply(c("ctrl", "mut"), function(g) {
    stats <- summarise_groups(level_df[level_df$genotype == g, ], "feature_type",
                              "mch_pct")
    stats$genotype <- g
    stats
  }))

  enhancer_stats <- do.call(rbind, lapply(ENHANCER_LEVELS, function(status) {
    stats <- summarise_groups(plot_df[plot_df$enhancer_status == status, ],
                              "feature_type", "delta_mch_pct")
    stats$enhancer_status <- factor(status, levels = ENHANCER_LEVELS)
    stats
  }))

  kw <- run_kruskal_wallis(intervals)
  dunn_posthoc <- run_dunn_posthoc(intervals)
  wilcox_zero <- run_wilcoxon_vs_zero(intervals)
  enhancer_summary <- enhancer_overlap_summary(intervals)
  enhancer_wilcox <- run_enhancer_wilcoxon(intervals)

  rank_profile <- exon_intron_rank_profile(intervals, opt$max_rank)
  metagene <- metagene_profile(intervals, opt$metagene_bins, opt$min_coverage)

  cat("\nRunning gene-level Fisher tests...\n")
  fisher_genes <- build_fisher_gene_table(gene_feature, gene_mch)
  fisher_summary <- run_gene_fisher_tests(fisher_genes, out_dir)

  key_intervals <- intervals[gene_name %in% KEY_GENES]
  if (nrow(key_intervals) == 0) {
    stop("No KEY_GENES interval survives the coverage filter: ",
         paste(KEY_GENES, collapse = ", "))
  }
  key_intervals <- merge(key_intervals, as.data.table(gene_mch), by = "gene_name",
                         all.x = TRUE)
  key_intervals[, midpoint_mb := (start + end) / 2e6]
  key_intervals[, facet_label := sprintf("%s (%s, gene logFC = %.2f, FDR = %.2g)",
                                         gene_name, chr, gene_edger_logFC,
                                         gene_edger_fdr)]
  setorder(key_intervals, gene_name, start)
  n_key_genes <- uniqueN(key_intervals$gene_name)
  cat(sprintf("\nKey genes: %d of %d have intervals passing the filter (%s intervals)\n",
              n_key_genes, length(KEY_GENES), fmt_int(nrow(key_intervals))))

  cat("\nWriting tables...\n")
  write_section_table(intervals, file.path(out_dir, "50_01_feature_intervals.tsv"))
  write_section_table(gene_feature, file.path(out_dir, "50_01_gene_feature_delta.tsv"))
  write_section_table(type_summary, file.path(out_dir, "50_01_feature_type_summary.tsv"))
  write_section_table(sweep, file.path(out_dir, "50_01_coverage_threshold_sweep.tsv"))
  write_section_table(delta_stats, file.path(out_dir, "50_01_delta_by_feature_stats.tsv"))
  write_section_table(level_stats, file.path(out_dir, "50_01_mch_level_by_feature_stats.tsv"))
  write_section_table(kw$table, file.path(out_dir, "50_01_kruskal_wallis.tsv"))
  write_section_table(dunn_posthoc, file.path(out_dir, "50_01_dunn_posthoc.tsv"))
  write_section_table(wilcox_zero, file.path(out_dir, "50_01_wilcoxon_vs_zero.tsv"))
  write_section_table(rank_profile, file.path(out_dir, "50_01_exon_intron_rank_profile.tsv"))
  write_section_table(metagene$profile, file.path(out_dir, "50_01_metagene_profile.tsv"))
  write_section_table(enhancer_summary, file.path(out_dir, "50_01_enhancer_overlap_summary.tsv"))
  write_section_table(enhancer_stats, file.path(out_dir, "50_01_enhancer_split_stats.tsv"))
  write_section_table(enhancer_wilcox, file.path(out_dir, "50_01_enhancer_wilcoxon.tsv"))
  write_section_table(key_intervals, file.path(out_dir, "50_01_key_gene_intervals.tsv"))
  write_section_table(fisher_summary, file.path(out_dir, "50_01_gene_fisher_summary.tsv"))

  cat("\nCreating figures...\n")

  # Figure text lives beside the plot call, never in the tables written above.
  delta_labels <- delta_stats
  delta_labels$label <- group_label(delta_stats)
  level_labels <- level_stats
  level_labels$label <- group_label(level_stats)
  enhancer_labels <- enhancer_stats
  enhancer_labels$label <- group_label(enhancer_stats)

  save_multiformat_ggplot(
    plot_feature_overview(type_summary, opt$min_coverage),
    file.path(out_dir, "50_01a_feature_distribution"),
    width = 16, height = 6)

  save_multiformat_ggplot(
    plot_delta_by_feature(plot_df, delta_labels, kw$table),
    file.path(out_dir, "50_01b_delta_mch_by_feature"),
    width = 11, height = 7)

  save_multiformat_ggplot(
    plot_dunn_heatmap(build_dunn_matrix(dunn_posthoc, FEATURE_TYPES)),
    file.path(out_dir, "50_01c_dunn_posthoc_heatmap"),
    width = 10, height = 8)

  save_multiformat_ggplot(
    plot_key_gene_loci(key_intervals, n_key_genes),
    file.path(out_dir, "50_01d_key_gene_feature_lollipop"),
    width = 14, height = 3.4 * ceiling(n_key_genes / 2))

  save_multiformat_ggplot(
    plot_metagene(metagene$profile, opt$metagene_bins, metagene$n_genes),
    file.path(out_dir, "50_01e_metagene_profile"),
    width = 11, height = 9)

  save_multiformat_ggplot(
    plot_enhancer_split(plot_df, enhancer_labels, enhancer_wilcox),
    file.path(out_dir, "50_01f_enhancer_split_delta"),
    width = 12, height = 7)

  save_multiformat_ggplot(
    plot_rank_profile(rank_profile, opt$max_rank),
    file.path(out_dir, "50_01g_exon_intron_rank_profile"),
    width = 12, height = 9)

  save_multiformat_ggplot(
    plot_mch_level_by_feature(level_df, level_labels),
    file.path(out_dir, "50_01h_mch_level_by_feature"),
    width = 12, height = 7)

  cat("\nSection 50_01 complete.\n")
  cat(sprintf("  Intervals loaded:  %s\n", fmt_int(nrow(all_intervals))))
  cat(sprintf("  Intervals used:    %s\n", fmt_int(nrow(intervals))))
  cat(sprintf("  Genes used:        %s\n", fmt_int(uniqueN(intervals$gene_name))))
  cat(sprintf("  Fisher tests registered: %d\n", nrow(fisher_summary)))
  cat(sprintf("  Output: %s\n\n", out_dir))
}

main()
