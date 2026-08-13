# scripts/sections/40_permutation/40_03_domains.R
#
# Section 40_03: Genomic permutation validation of mCH direction against large
# chromatin domains.
#
# What this tests
#   Gene-level Fisher tests assume that genes are independent draws. Megabase
#   A/B compartments and Polycomb domains break that assumption: neighbouring
#   genes share a domain, so an enrichment can arise from position alone. This
#   section repeats the domain questions with regioneReloaded permutation
#   tests. Each mCH direction set is randomised across the genome with
#   randomizeRegions, which keeps interval count, interval size, and chromosome
#   assignment, and the observed overlap is scored against that null.
#
#   Two sub-analyses, each with the same two-row Alist (mCH hypermethylated
#   gene bodies, mCH hypomethylated gene bodies):
#     Compartment  A compartment, B compartment, B to A shift, A to B shift
#     Polycomb     merged H3K27me3 domains, Polycomb-state intervals
#
#   Every permutation cell also gets a plain gene-level Fisher test on the same
#   intervals, so the odds ratio and the permutation z-score sit in one table.
#
# Reads
#   mch_results                     gene-level mCH results, from _shared_config.R
#   HIC_PATHS$compartments          HOMER getDiffExpression output. 25 kb bins,
#                                   six per-sample PC1 columns (ctrl_M1..mut_M3),
#                                   a "ctrl vs. mut Difference" column, and a
#                                   "ctrl vs. mut adj. p-value" column.
#   CHIP_PATHS$h3k27me3             H3K27me3 consensus peaks, merged into domains
#   HANDOFF_PATHS$chromatin_state   gene-level chromatin state written by 10_01.
#                                   The Polycomb-state intervals are the gene
#                                   bodies whose body_state is "Polycomb".
#
# Writes (into OUT_DIR, default OUTPUT_PATHS$permutation)
#   Caches   40_03_perm_compartment.rds, 40_03_perm_polycomb.rds,
#            40_03_localz_compartment.rds, 40_03_localz_polycomb.rds
#            Set FORCE_RERUN to any non-empty value to recompute them.
#   Tables   40_03_association_compartment.tsv, 40_03_association_polycomb.tsv,
#            40_03_fisher_tests.tsv, 40_03_fisher_vs_permutation.tsv,
#            40_03_region_set_summary.tsv, 40_03_compartment_bin_summary.tsv,
#            40_03_polycomb_domain_summary.tsv,
#            40_03_key_genes_domain_membership.tsv,
#            40_03_local_zscore_compartment.tsv,
#            40_03_local_zscore_polycomb.tsv,
#            40_03_permutation_parameters.tsv
#   Figures  40_03a_crosswise_mch_x_compartment, 40_03b_forest_mch_x_compartment,
#            40_03c_crosswise_mch_x_polycomb, 40_03d_forest_mch_x_polycomb,
#            40_03e_fisher_vs_permutation, 40_03f_observed_vs_expected_overlaps,
#            40_03g_local_zscore_compartment, 40_03h_local_zscore_polycomb
#   Fisher gene tables under fisher_tables/ and rows in the shared registry,
#   which section 40_04 validates by label shuffle.
#
# Adapted from Biomodal section 36 (section_36_permutation_domains.R). That
# script carried four Alist rows (mC and hmC, both directions); EM-seq measures
# one methylation modality, so the Alist here has two rows.
#
# Run:
#   Rscript scripts/sections/40_permutation/40_03_domains.R 5000 2>&1 | tee logs/40_03.txt

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)
library(regioneR)
library(regioneReloaded)
library(BSgenome.Mmusculus.UCSC.mm10)

# chooseHclustMet picks a clustering method by cophenetic correlation, which is
# undefined for the two-row matrix this section builds.
patch_chooseHclustMet()

# =============================================================================
# SECTION CONSTANTS
# =============================================================================

SECTION_ID <- "40_03"

# Permutation runs on the assembled autosomes plus chrX. chrY is left out: the
# cohort is mixed sex and the mCH results carry both sexes on one axis.
PERM_CHRS <- paste0("chr", c(1:19, "X"))

PERM_SEED <- 42L
PERM_RANFUN <- "randomizeRegions"
PERM_EVFUN <- "numOverlaps"
PERM_PER_CHROMOSOME <- TRUE
DEFAULT_NTIMES <- 5000L
DEFAULT_CORES <- 8L

# Order of the six per-sample PC1 columns in the HOMER compartment header.
PC1_SAMPLE_ORDER <- c("ctrl_M1", "ctrl_M2", "ctrl_M3",
                      "mut_M1", "mut_M2", "mut_M3")

COMPARTMENT_LEVELS <- c("A", "B")
SHIFT_LEVELS <- c("B to A", "A to B", "Stable")

# The two Alist rows, shared by both sub-analyses.
HYPER_SET <- "mCH Hyper genes"
HYPO_SET <- "mCH Hypo genes"
DIRECTION_SETS <- c(HYPER_SET, HYPO_SET)

DIRECTION_SET_COLORS <- c(
  "mCH Hyper genes" = unname(COLORS$direction["Hypermethylated"]),
  "mCH Hypo genes"  = unname(COLORS$direction["Hypomethylated"])
)

SUB_ANALYSIS_LEVELS <- c("Compartment", "Polycomb")

CONCORDANCE_LEVELS <- c("Confirmed", "Weakened", "Strengthened",
                        "Concordant NS")

CONCORDANCE_COLORS <- c(
  "Confirmed"     = "#2CA02C",
  "Weakened"      = "#D62728",
  "Strengthened"  = "#FF7F0E",
  "Concordant NS" = "grey60"
)

# Columns the regioneReloaded multiOverlaps data.frames must carry.
PERM_RESULT_COLUMNS <- c("name", "n_regionA", "n_regionB", "z_score",
                         "p_value", "n_overlaps", "mean_perm_test",
                         "sd_perm_test", "norm_zscore", "adj.p_value")

# Columns section 10_01 writes into gene_chromatin_state.tsv.
CHROMATIN_STATE_COLUMNS <- c(
  "gene_name", "gene_id", "chr", "start", "end",
  "promoter_state", "body_state",
  "prom_ctcf_overlap", "prom_h3k27ac_overlap", "prom_h3k27me3_overlap",
  "prom_h3k4me1_overlap", "prom_h3k4me3_overlap", "prom_bivalent_overlap",
  "body_ctcf_overlap", "body_h3k27ac_overlap", "body_h3k27me3_overlap",
  "body_h3k4me1_overlap", "body_h3k4me3_overlap", "body_bivalent_overlap")

# Two-sided normal cutoff drawn on the forest and comparison plots.
Z_CUTOFF <- 1.96

# =============================================================================
# OPTIONS
# =============================================================================

parse_options <- function() {
  option_list <- list(
    make_option("--output-dir", dest = "output_dir", type = "character",
                default = OUTPUT_PATHS$permutation,
                help = "Directory for caches, figures, and tables [default: %default]"),
    make_option("--fdr-threshold", dest = "fdr_threshold", type = "double",
                default = Q_THRESHOLD,
                help = "edgeR FDR cutoff calling a gene differentially methylated [default: %default]"),
    make_option("--shift-fdr", dest = "shift_fdr", type = "double",
                default = 0.05,
                help = "Adjusted p cutoff for a compartment shift [default: %default]"),
    make_option("--shift-diff", dest = "shift_diff", type = "double",
                default = 0.30,
                help = "Minimum |PC1 difference| for a compartment shift [default: %default]"),
    make_option("--k27me3-merge-gap", dest = "k27me3_merge_gap", type = "integer",
                default = 5000L,
                help = "Gap in bp below which neighbouring H3K27me3 peaks join one domain [default: %default]"),
    make_option("--lz-window", dest = "lz_window", type = "integer",
                default = 50000L,
                help = "Half-width in bp of the local z-score window [default: %default]"),
    make_option("--lz-step", dest = "lz_step", type = "integer",
                default = 1000L,
                help = "Step in bp between local z-score offsets [default: %default]")
  )

  parsed <- parse_args(
    OptionParser(
      option_list = option_list,
      usage = "%prog [options] [ntimes]",
      description = sprintf(paste(
        "Permutation validation of mCH direction against A/B compartments and",
        "Polycomb domains. The optional positional argument is the number of",
        "permutations (default %d)."), DEFAULT_NTIMES)
    ),
    positional_arguments = TRUE
  )

  opt <- parsed$options
  opt$ntimes <- resolve_ntimes(parsed$args)
  opt$cores <- resolve_cores()
  opt$force_rerun <- nzchar(Sys.getenv("FORCE_RERUN", unset = ""))

  if (opt$fdr_threshold <= 0 || opt$fdr_threshold >= 1) {
    stop("--fdr-threshold must be between 0 and 1, got ", opt$fdr_threshold)
  }
  if (opt$shift_fdr <= 0 || opt$shift_fdr >= 1) {
    stop("--shift-fdr must be between 0 and 1, got ", opt$shift_fdr)
  }
  if (opt$shift_diff <= 0) {
    stop("--shift-diff must be greater than 0, got ", opt$shift_diff)
  }
  if (opt$k27me3_merge_gap < 1) {
    stop("--k27me3-merge-gap must be at least 1 bp, got ", opt$k27me3_merge_gap)
  }
  if (opt$lz_window < 1 || opt$lz_step < 1) {
    stop("--lz-window and --lz-step must both be at least 1 bp, got ",
         opt$lz_window, " and ", opt$lz_step)
  }
  if (opt$lz_step > opt$lz_window) {
    stop("--lz-step (", opt$lz_step, ") is larger than --lz-window (",
         opt$lz_window, "), which leaves fewer than two offsets.")
  }
  opt
}

#' Read the permutation count from the positional command-line argument.
resolve_ntimes <- function(positional) {
  if (length(positional) == 0) return(DEFAULT_NTIMES)
  if (length(positional) > 1) {
    stop("Only one positional argument is allowed (ntimes). Got: ",
         paste(positional, collapse = " "))
  }
  n <- suppressWarnings(as.integer(positional[1]))
  if (is.na(n) || n < 1) {
    stop("ntimes must be a positive whole number, got '", positional[1], "'")
  }
  n
}

#' Read the core count from SLURM_CPUS_PER_TASK.
resolve_cores <- function() {
  slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
  if (!nzchar(slurm_cpus)) return(DEFAULT_CORES)
  n <- suppressWarnings(as.integer(slurm_cpus))
  if (is.na(n) || n < 1) {
    stop("SLURM_CPUS_PER_TASK is set to '", slurm_cpus,
         "', which is not a positive whole number.")
  }
  n
}

# =============================================================================
# SMALL HELPERS
# =============================================================================

fmt_int <- function(x) format(x, big.mark = ",", trim = TRUE)

fmt_p <- function(p) {
  if (length(p) != 1 || is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

sig_stars <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return("")
    if (x < 0.001) return("***")
    if (x < 0.01) return("**")
    if (x < 0.05) return("*")
    "ns"
  }, character(1))
}

#' Turn a region set display name into a column-safe identifier.
region_set_slug <- function(x) {
  slug <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("^_+|_+$", "", slug)
}

#' Interval count, total bp, and median width of a region set.
region_set_stats <- function(gr, set_name, sub_analysis, role) {
  data.frame(
    sub_analysis = sub_analysis,
    role = role,
    region_set = set_name,
    n_intervals = length(gr),
    total_bp = sum(as.numeric(width(gr))),
    median_width_bp = median(width(gr)),
    mean_width_bp = mean(width(gr)),
    max_width_bp = max(width(gr)),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# GENOME AND REGION SET PREPARATION
# =============================================================================

#' Build the mm10 permutation genome, restricted to PERM_CHRS.
build_permutation_genome <- function() {
  genome <- getGenomeAndMask("mm10")$genome
  genome <- genome[as.character(seqnames(genome)) %in% PERM_CHRS]
  seqlevels(genome) <- PERM_CHRS

  if (length(genome) != length(PERM_CHRS)) {
    stop("mm10 genome has ", length(genome), " ranges after restricting to ",
         length(PERM_CHRS), " chromosomes. Expected one range per chromosome.")
  }

  cat(sprintf("  mm10 permutation genome: %d chromosomes, %s bp\n",
              length(genome), fmt_int(sum(as.numeric(width(genome))))))
  genome
}

#' Restrict a region set to PERM_CHRS and align its seqlevels to the genome.
#'
#' Metadata columns are dropped: crosswisePermTest only uses the coordinates,
#' and the region sets are cached as RDS.
restrict_to_permutation_chrs <- function(gr, set_name) {
  n_before <- length(gr)
  gr <- gr[as.character(seqnames(gr)) %in% PERM_CHRS]
  seqlevels(gr) <- PERM_CHRS
  mcols(gr) <- NULL

  if (length(gr) == 0) {
    stop("Region set '", set_name, "' has no interval on ",
         paste(PERM_CHRS, collapse = ", "))
  }

  cat(sprintf("    %-26s %8s intervals (%s off-target dropped), median width %s bp\n",
              set_name, fmt_int(length(gr)), fmt_int(n_before - length(gr)),
              fmt_int(median(width(gr)))))
  gr
}

# =============================================================================
# GENE UNIVERSE
# =============================================================================

#' Keep one row per gene symbol, the one with the largest |edger_logFC|.
#'
#' Some symbols carry more than one ENSMUSG identifier. One interval per symbol
#' keeps the permutation Alist and the gene-level Fisher tables in agreement.
dedup_mch_by_gene <- function(mch) {
  n_dup <- sum(duplicated(mch$gene_name))
  out <- mch %>%
    dplyr::group_by(gene_name) %>%
    dplyr::slice_max(abs(edger_logFC), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    as.data.frame()
  cat(sprintf("  mCH genes: %s tested, %d duplicate symbols collapsed, %s kept\n",
              fmt_int(nrow(mch)), n_dup, fmt_int(nrow(out))))
  out
}

#' Build the gene universe: one row per symbol, on PERM_CHRS, with the
#' direction flags recomputed at the section FDR threshold.
build_gene_universe <- function(mch, fdr_threshold) {
  if (any(is.na(mch$edger_fdr))) {
    stop("edger_fdr is NA for ", sum(is.na(mch$edger_fdr)), " genes.")
  }
  if (any(is.na(mch$mch_diff))) {
    stop("mch_diff is NA for ", sum(is.na(mch$mch_diff)), " genes.")
  }

  genes <- dedup_mch_by_gene(mch)

  n_before <- nrow(genes)
  genes <- genes[genes$chr %in% PERM_CHRS, , drop = FALSE]
  cat(sprintf("  Genes on the %d permutation chromosomes: %s of %s\n",
              length(PERM_CHRS), fmt_int(nrow(genes)), fmt_int(n_before)))

  genes$mch_sig <- genes$edger_fdr < fdr_threshold
  genes$mch_hyper <- genes$mch_sig & genes$edger_logFC > 0
  genes$mch_hypo <- genes$mch_sig & genes$edger_logFC < 0
  genes$mch_direction <- ifelse(genes$edger_logFC > 0,
                                "Hypermethylated", "Hypomethylated")

  if (sum(genes$mch_hyper) == 0 || sum(genes$mch_hypo) == 0) {
    stop("Both mCH directions are needed. Found ", sum(genes$mch_hyper),
         " hypermethylated and ", sum(genes$mch_hypo), " hypomethylated genes ",
         "at FDR < ", fdr_threshold)
  }

  cat(sprintf("  Significant at FDR < %.3f: %s genes (%s hyper, %s hypo)\n",
              fdr_threshold, fmt_int(sum(genes$mch_sig)),
              fmt_int(sum(genes$mch_hyper)), fmt_int(sum(genes$mch_hypo))))
  genes
}

#' Gene-body GRanges for the universe.
#'
#' mch_results coordinates are 0-based half-open, matching the construction of
#' gene_bodies in _shared_config.R.
gene_universe_granges <- function(genes) {
  GRanges(
    seqnames = genes$chr,
    ranges = IRanges(start = genes$start + 1L, end = genes$end)
  )
}

#' The two-row Alist: hypermethylated and hypomethylated gene bodies.
build_direction_region_sets <- function(genes, gene_gr) {
  cat("  Alist (mCH direction sets):\n")
  sets <- list()
  sets[[HYPER_SET]] <- restrict_to_permutation_chrs(gene_gr[genes$mch_hyper],
                                                    HYPER_SET)
  sets[[HYPO_SET]] <- restrict_to_permutation_chrs(gene_gr[genes$mch_hypo],
                                                   HYPO_SET)
  sets
}

# =============================================================================
# SUB-ANALYSIS A: A/B COMPARTMENTS
# =============================================================================

#' Read the HOMER compartment table and classify every bin.
#'
#' Column names hold full command lines and file paths, so the reader keeps them
#' verbatim (check.names = FALSE) and finds the columns by pattern. Compartment
#' A is mean control PC1 > 0. The shift category comes from the mutant-minus-
#' control PC1 difference and its adjusted p-value. This matches the parsing in
#' section 10_02.
load_compartment_bins <- function(path, shift_fdr, shift_diff) {
  if (!file.exists(path)) stop("Compartment file not found: ", path)

  raw <- read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE, comment.char = "", quote = "")
  cat(sprintf("  Loaded %s bins from %s\n", fmt_int(nrow(raw)), basename(path)))

  pc1_cols <- grep("bedGraph avg over given bp", names(raw), value = TRUE)
  if (length(pc1_cols) != 6) {
    stop("Expected 6 per-sample PC1 columns matching ",
         "'bedGraph avg over given bp', found ", length(pc1_cols))
  }
  for (i in seq_along(PC1_SAMPLE_ORDER)) {
    if (!grepl(PC1_SAMPLE_ORDER[i], pc1_cols[i], fixed = TRUE)) {
      stop("PC1 column ", i, " does not name sample ", PC1_SAMPLE_ORDER[i],
           ". Column header: ", pc1_cols[i])
    }
  }

  difference_col <- grep("ctrl vs\\. mut Difference", names(raw), value = TRUE)
  if (length(difference_col) != 1) {
    stop("Expected exactly 1 'ctrl vs. mut Difference' column, found ",
         length(difference_col))
  }
  adj_pvalue_col <- grep("ctrl vs\\. mut adj\\. p-value", names(raw), value = TRUE)
  if (length(adj_pvalue_col) != 1) {
    stop("Expected exactly 1 'ctrl vs. mut adj. p-value' column, found ",
         length(adj_pvalue_col))
  }

  missing_coords <- setdiff(c("Chr", "Start", "End"), names(raw))
  if (length(missing_coords) > 0) {
    stop("Compartment file is missing coordinate columns: ",
         paste(missing_coords, collapse = ", "))
  }

  cat(sprintf("  Control PC1 columns: %s\n",
              paste(PC1_SAMPLE_ORDER[1:3], collapse = ", ")))
  cat(sprintf("  Difference column:   %s\n", difference_col))
  cat(sprintf("  Adjusted p column:   %s\n", adj_pvalue_col))

  ctrl_pc1 <- vapply(pc1_cols[1:3], function(cl) as.numeric(raw[[cl]]),
                     numeric(nrow(raw)))
  mut_pc1 <- vapply(pc1_cols[4:6], function(cl) as.numeric(raw[[cl]]),
                    numeric(nrow(raw)))

  bins <- data.frame(
    bin_chr = as.character(raw[["Chr"]]),
    bin_start = as.integer(raw[["Start"]]),
    bin_end = as.integer(raw[["End"]]),
    mean_ctrl_pc1 = rowMeans(ctrl_pc1),
    mean_mut_pc1 = rowMeans(mut_pc1),
    pc1_difference = as.numeric(raw[[difference_col]]),
    pc1_adj_pvalue = as.numeric(raw[[adj_pvalue_col]]),
    stringsAsFactors = FALSE
  )

  n_na <- sum(!complete.cases(bins))
  if (n_na > 0) {
    stop(n_na, " compartment bins hold NA in a coordinate, PC1, difference, ",
         "or adjusted p-value field.")
  }

  bins$compartment <- factor(ifelse(bins$mean_ctrl_pc1 > 0, "A", "B"),
                             levels = COMPARTMENT_LEVELS)

  shift <- rep("Stable", nrow(bins))
  passes <- bins$pc1_adj_pvalue < shift_fdr
  shift[passes & bins$pc1_difference > shift_diff] <- "B to A"
  shift[passes & bins$pc1_difference < -shift_diff] <- "A to B"
  bins$shift <- factor(shift, levels = SHIFT_LEVELS)

  n_before <- nrow(bins)
  bins <- bins[bins$bin_chr %in% PERM_CHRS, , drop = FALSE]
  cat(sprintf("  Bins on the permutation chromosomes: %s of %s\n",
              fmt_int(nrow(bins)), fmt_int(n_before)))

  cat(sprintf("  Compartment A: %s bins, B: %s bins\n",
              fmt_int(sum(bins$compartment == "A")),
              fmt_int(sum(bins$compartment == "B"))))
  cat(sprintf("  Shift (adj p < %.2f, |difference| > %.2f): B to A %s, A to B %s, Stable %s\n",
              shift_fdr, shift_diff,
              fmt_int(sum(bins$shift == "B to A")),
              fmt_int(sum(bins$shift == "A to B")),
              fmt_int(sum(bins$shift == "Stable"))))
  bins
}

#' Bin counts per compartment and shift category.
compartment_bin_summary <- function(bins) {
  counts <- as.data.frame(table(compartment = bins$compartment,
                                shift = bins$shift),
                          stringsAsFactors = FALSE)
  names(counts)[3] <- "n_bins"
  counts$pct_of_bins <- 100 * counts$n_bins / nrow(bins)
  counts
}

#' Build the domain region set for one bin class.
#'
#' HOMER bins are half-open, so adding 1 to Start keeps neighbouring bins
#' disjoint. reduce() then merges bins of the same class that touch or overlap,
#' so a gene crossing a run of same-class bins counts as one overlap rather than
#' one overlap per bin.
bins_to_domains <- function(bins, keep, set_name) {
  selected <- bins[keep, , drop = FALSE]
  if (nrow(selected) == 0) {
    stop("Compartment class '", set_name, "' holds no bin. ",
         "Check --shift-fdr and --shift-diff.")
  }

  gr <- GRanges(
    seqnames = selected$bin_chr,
    ranges = IRanges(start = selected$bin_start + 1L, end = selected$bin_end)
  )
  domains <- GenomicRanges::reduce(gr)

  cat(sprintf("    %-26s %8s bins -> %s domains\n",
              set_name, fmt_int(nrow(selected)), fmt_int(length(domains))))
  domains
}

#' The four compartment region sets, each merged into contiguous domains.
build_compartment_region_sets <- function(bins) {
  cat("  Blist (compartment domains):\n")
  raw_sets <- list(
    "A Compartment" = bins_to_domains(bins, bins$compartment == "A",
                                      "A Compartment"),
    "B Compartment" = bins_to_domains(bins, bins$compartment == "B",
                                      "B Compartment"),
    "B to A Shift"  = bins_to_domains(bins, bins$shift == "B to A",
                                      "B to A Shift"),
    "A to B Shift"  = bins_to_domains(bins, bins$shift == "A to B",
                                      "A to B Shift")
  )

  sets <- lapply(names(raw_sets), function(nm) {
    restrict_to_permutation_chrs(raw_sets[[nm]], nm)
  })
  names(sets) <- names(raw_sets)

  for (nm in names(sets)) {
    if (length(sets[[nm]]) < 50) {
      cat(sprintf(paste("    Note: '%s' has %d domains. A permutation z-score",
                        "on that few intervals carries wide sampling error.\n"),
                  nm, length(sets[[nm]])))
    }
  }
  sets
}

# =============================================================================
# SUB-ANALYSIS B: POLYCOMB DOMAINS
# =============================================================================

#' Merge H3K27me3 consensus peaks into domains.
#'
#' Peaks separated by less than merge_gap bp join one domain, so a run of
#' neighbouring H3K27me3 peaks is scored as the single broad region it marks.
build_k27me3_domains <- function(path, merge_gap) {
  peaks <- load_chip_peaks(path, "H3K27me3 consensus")
  domains <- GenomicRanges::reduce(peaks, min.gapwidth = merge_gap)

  cat(sprintf("    %s peaks merged at a %s bp gap -> %s domains, median width %s bp\n",
              fmt_int(length(peaks)), fmt_int(merge_gap),
              fmt_int(length(domains)), fmt_int(median(width(domains)))))

  list(peaks = peaks, domains = domains)
}

#' Read the Polycomb-body gene intervals from the section 10_01 handoff table.
#'
#' body_state is the chromatin state of the gene body. It is the column that
#' answers the Polycomb question; promoter_state describes the TSS window
#' instead.
#'
#' The handoff coordinates come from mch_results, so they are 0-based half-open
#' and get the same +1 as gene_bodies. Overlapping gene bodies are reduced to
#' disjoint intervals.
load_polycomb_state_intervals <- function(path) {
  if (!file.exists(path)) {
    stop("Gene chromatin state table not found: ", path,
         "\nRun section 10_01 first; it writes this file.")
  }

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")

  missing <- setdiff(CHROMATIN_STATE_COLUMNS, colnames(df))
  if (length(missing) > 0) {
    stop("Chromatin state table from section 10_01 is missing columns: ",
         paste(missing, collapse = ", "), "\nFile: ", path)
  }

  unknown_states <- setdiff(unique(df$body_state), BODY_STATE_ORDER)
  if (length(unknown_states) > 0) {
    stop("Chromatin state table holds body states outside BODY_STATE_ORDER: ",
         paste(unknown_states, collapse = ", "))
  }

  polycomb <- df[df$body_state == "Polycomb", , drop = FALSE]
  if (nrow(polycomb) == 0) {
    stop("No gene carries the Polycomb body state in ", path)
  }

  gr <- GRanges(
    seqnames = polycomb$chr,
    ranges = IRanges(start = polycomb$start + 1L, end = polycomb$end)
  )
  intervals <- GenomicRanges::reduce(gr)

  cat(sprintf("    %s Polycomb-body genes -> %s disjoint intervals, median width %s bp\n",
              fmt_int(nrow(polycomb)), fmt_int(length(intervals)),
              fmt_int(median(width(intervals)))))

  list(genes = polycomb, intervals = intervals)
}

#' The two Polycomb region sets.
build_polycomb_region_sets <- function(k27me3_path, merge_gap, state_path) {
  cat("  Blist (Polycomb domains):\n")
  k27me3 <- build_k27me3_domains(k27me3_path, merge_gap)
  state <- load_polycomb_state_intervals(state_path)

  sets <- list(
    "H3K27me3 Domains" = restrict_to_permutation_chrs(k27me3$domains,
                                                      "H3K27me3 Domains"),
    "Polycomb State Intervals" = restrict_to_permutation_chrs(
      state$intervals, "Polycomb State Intervals")
  )

  summary_table <- data.frame(
    region_set = c("H3K27me3 Domains", "Polycomb State Intervals"),
    source = c(basename(k27me3_path), basename(state_path)),
    n_input_intervals = c(length(k27me3$peaks), nrow(state$genes)),
    merge_gap_bp = c(merge_gap, 1L),
    n_merged_intervals = c(length(k27me3$domains), length(state$intervals)),
    n_after_chr_filter = c(length(sets[["H3K27me3 Domains"]]),
                           length(sets[["Polycomb State Intervals"]])),
    stringsAsFactors = FALSE
  )

  list(sets = sets, summary = summary_table)
}

# =============================================================================
# PERMUTATION TESTS
# =============================================================================

#' Load a cached permutation object, or build it and cache it with its count.
#'
#' The permutation count is stored beside the object. Without it, a cache built
#' at one --ntimes would be reported under the --ntimes of a later run, and
#' every figure subtitle and parameter table would name a count that never ran.
#'
#' @param cache_path RDS file holding list(object, ntimes).
#' @param label Display name used in the messages.
#' @param opt Parsed options; supplies ntimes and force_rerun.
#' @param compute Function of no arguments returning the object to cache.
#' @return list with the object and the count it was built at
load_or_run_cached <- function(cache_path, label, opt, compute) {
  if (file.exists(cache_path) && !opt$force_rerun) {
    cat(sprintf("  Cached %s found; loading %s\n", label, basename(cache_path)))
    cached <- readRDS(cache_path)

    missing <- setdiff(c("object", "ntimes"), names(cached))
    if (length(missing) > 0) {
      stop("Cache ", cache_path, " is missing entries: ",
           paste(missing, collapse = ", "),
           ". Set FORCE_RERUN=1 to rebuild it.")
    }

    cat(sprintf("    Cache holds %s permutations.\n", fmt_int(cached$ntimes)))
    if (cached$ntimes != opt$ntimes) {
      cat(sprintf(paste("    This run asked for %s. The cached count is what",
                        "gets reported. Set FORCE_RERUN=1 to rerun.\n"),
                  fmt_int(opt$ntimes)))
    }
    return(cached)
  }

  cached <- list(object = compute(), ntimes = opt$ntimes)
  saveRDS(cached, cache_path)
  cat(sprintf("  Saved cache: %s\n", cache_path))
  cached
}

#' Set opt$ntimes to the count the cached objects actually hold.
#'
#' Every cached object records its own count. A run that mixes counts cannot be
#' summarised by one number, so it stops instead.
reconcile_cached_ntimes <- function(opt, cached_objects) {
  counts <- unique(vapply(cached_objects,
                          function(x) as.numeric(x$ntimes), numeric(1)))

  if (length(counts) > 1) {
    stop("The cached permutation objects were built at different permutation ",
         "counts: ", paste(fmt_int(sort(counts)), collapse = ", "),
         ". Set FORCE_RERUN=1 to rebuild them all at ", fmt_int(opt$ntimes), ".")
  }

  if (counts != opt$ntimes) {
    cat(sprintf(paste("  Reporting %s permutations, the count every cache holds,",
                      "not the %s this run asked for.\n"),
                fmt_int(counts), fmt_int(opt$ntimes)))
    opt$ntimes <- counts
  }
  opt
}

#' Run one crosswise permutation test, or load it from the RDS cache.
#'
#' Returns the cache record, so the caller can read the count the object was
#' built at. The genoMatriXeR object is in the "object" entry.
run_crosswise <- function(Alist, Blist, genome, opt, cache_path, label) {
  load_or_run_cached(
    cache_path, sprintf("%s permutation", label), opt,
    function() {
      cat(sprintf("  Running crosswisePermTest for %s: %d x %d = %d cells\n",
                  label, length(Alist), length(Blist),
                  length(Alist) * length(Blist)))
      cat(sprintf("    ntimes = %s, cores = %d, ranFUN = %s, evFUN = %s, per.chromosome = %s, seed = %d\n",
                  fmt_int(opt$ntimes), opt$cores, PERM_RANFUN, PERM_EVFUN,
                  PERM_PER_CHROMOSOME, PERM_SEED))

      options(mc.cores = opt$cores)
      set.seed(PERM_SEED)
      cw <- crosswisePermTest(
        Alist          = Alist,
        Blist          = Blist,
        genome         = genome,
        ranFUN         = PERM_RANFUN,
        evFUN          = PERM_EVFUN,
        ntimes         = opt$ntimes,
        force.parallel = TRUE,
        per.chromosome = PERM_PER_CHROMOSOME
      )

      # symm_matrix = FALSE: the design is 2 x N, not square. The patched
      # chooseHclustMet handles the two-row clustering.
      makeCrosswiseMatrix(cw, pvcut = 1, symm_matrix = FALSE,
                          hc.method = "average")
    })
}

#' Flatten a genoMatriXeR object into one row per Alist x Blist cell.
extract_association_table <- function(cw, sub_analysis) {
  evaluation <- getMultiEvaluation(cw)

  if (length(evaluation) == 0) {
    stop("getMultiEvaluation() returned no result for sub-analysis ",
         sub_analysis)
  }

  rows <- lapply(names(evaluation), function(a_name) {
    df <- as.data.frame(evaluation[[a_name]])
    missing <- setdiff(PERM_RESULT_COLUMNS, colnames(df))
    if (length(missing) > 0) {
      stop("regioneReloaded result for '", a_name, "' is missing columns: ",
           paste(missing, collapse = ", "),
           ". Columns present: ", paste(colnames(df), collapse = ", "))
    }
    data.frame(
      sub_analysis = sub_analysis,
      region_set_a = a_name,
      region_set_b = as.character(df$name),
      n_region_a = df$n_regionA,
      n_region_b = df$n_regionB,
      observed_hits = df$n_overlaps,
      expected_hits = df$mean_perm_test,
      sd_perm_hits = df$sd_perm_test,
      perm_z_score = df$z_score,
      perm_norm_zscore = df$norm_zscore,
      perm_p_value = df$p_value,
      perm_adj_p_value = df[["adj.p_value"]],
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out$test_id <- paste(out$region_set_a, "x", out$region_set_b)
  out$perm_sig_label <- sig_stars(out$perm_adj_p_value)
  rownames(out) <- NULL
  out
}

#' Run one local z-score profile, or load it from the RDS cache.
#'
#' Returns the cache record. The multiLocalZScore object is in the "object"
#' entry and the count it was built at is in "ntimes".
run_local_zscore <- function(a_gr, Blist, genome, opt, cache_path, label) {
  load_or_run_cached(
    cache_path, sprintf("%s local z-score", label), opt,
    function() {
      cat(sprintf("  Running multiLocalZscore for %s: window = %s bp, step = %s bp, ntimes = %s\n",
                  label, fmt_int(opt$lz_window), fmt_int(opt$lz_step),
                  fmt_int(opt$ntimes)))

      options(mc.cores = opt$cores)
      set.seed(PERM_SEED)
      mlz <- multiLocalZscore(
        A        = a_gr,
        Blist    = Blist,
        ranFUN   = PERM_RANFUN,
        evFUN    = PERM_EVFUN,
        genome   = genome,
        window   = opt$lz_window,
        step     = opt$lz_step,
        force.parallel = TRUE,
        ntimes   = opt$ntimes
      )
      makeLZMatrix(mlz)
    })
}

#' Pull the local z-score summary table out of a multiLocalZScore object.
local_zscore_table <- function(mlz, sub_analysis) {
  evaluation <- getMultiEvaluation(mlz)
  if (!"resumeTable" %in% names(evaluation)) {
    stop("multiLocalZscore evaluation has no resumeTable. Names present: ",
         paste(names(evaluation), collapse = ", "))
  }
  tab <- as.data.frame(evaluation$resumeTable)
  tab$sub_analysis <- sub_analysis
  tab
}

# =============================================================================
# GENE-LEVEL FISHER TESTS
# =============================================================================

#' Overlap flags for every gene against every region set.
#'
#' Columns are named by region_set_slug(), one logical column per region set.
compute_membership <- function(gene_gr, region_sets) {
  flags <- lapply(region_sets, function(rs) countOverlaps(gene_gr, rs) > 0)
  names(flags) <- region_set_slug(names(region_sets))
  as.data.frame(flags)
}

#' One registered gene-level Fisher test per direction and region set.
#'
#' Each test asks whether the mCH direction set is enriched among the genes that
#' overlap the region set, against every other gene in the same universe the
#' permutation randomises over.
run_domain_fisher_tests <- function(genes, membership, region_sets,
                                    sub_analysis, out_dir) {
  directions <- list(
    list(set = HYPER_SET, col = "mch_hyper", label = "hypermethylated"),
    list(set = HYPO_SET, col = "mch_hypo", label = "hypomethylated")
  )

  rows <- list()
  for (target_name in names(region_sets)) {
    target_col <- region_set_slug(target_name)
    in_target <- membership[[target_col]]

    for (spec in directions) {
      gene_df <- data.frame(gene_name = genes$gene_name, chr = genes$chr,
                            stringsAsFactors = FALSE)
      gene_df[[spec$col]] <- genes[[spec$col]]
      gene_df$in_target <- in_target

      test_id <- sprintf("%s_x_%s", region_set_slug(spec$set), target_col)
      ft <- register_fisher_test(
        section = SECTION_ID,
        test_id = test_id,
        description = sprintf(
          "Are mCH %s genes enriched among genes overlapping %s?",
          spec$label, target_name),
        gene_df = gene_df,
        row_var = spec$col,
        col_var = "in_target",
        output_dir = out_dir
      )

      in_direction <- gene_df[[spec$col]]
      n_both <- sum(in_direction & in_target)
      n_direction_only <- sum(in_direction & !in_target)
      n_target_only <- sum(!in_direction & in_target)
      n_neither <- sum(!in_direction & !in_target)

      rows[[length(rows) + 1]] <- data.frame(
        sub_analysis = sub_analysis,
        region_set_a = spec$set,
        region_set_b = target_name,
        fisher_test_id = test_id,
        n_genes = nrow(gene_df),
        n_direction = sum(in_direction),
        n_in_target = sum(in_target),
        n_direction_in_target = n_both,
        n_direction_outside_target = n_direction_only,
        n_other_in_target = n_target_only,
        n_other_outside_target = n_neither,
        pct_direction_in_target = 100 * n_both / sum(in_direction),
        pct_other_in_target = 100 * n_target_only / (n_target_only + n_neither),
        fisher_or = unname(ft$estimate),
        fisher_ci_lower = ft$conf.int[1],
        fisher_ci_upper = ft$conf.int[2],
        fisher_p = ft$p.value,
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Join the Fisher and permutation results cell by cell.
merge_fisher_permutation <- function(fisher_df, perm_df) {
  merged <- dplyr::inner_join(
    fisher_df, perm_df,
    by = c("sub_analysis", "region_set_a", "region_set_b"))

  if (nrow(merged) != nrow(fisher_df)) {
    stop("Fisher rows (", nrow(fisher_df), ") and merged rows (", nrow(merged),
         ") disagree. A permutation cell has no matching Fisher test, or the ",
         "region set names differ between the two.")
  }

  merged$concordance <- factor(
    dplyr::case_when(
      merged$fisher_p < Q_THRESHOLD & merged$perm_adj_p_value < Q_THRESHOLD ~ "Confirmed",
      merged$fisher_p < Q_THRESHOLD & merged$perm_adj_p_value >= Q_THRESHOLD ~ "Weakened",
      merged$fisher_p >= Q_THRESHOLD & merged$perm_adj_p_value < Q_THRESHOLD ~ "Strengthened",
      TRUE ~ "Concordant NS"
    ),
    levels = CONCORDANCE_LEVELS)

  merged$fisher_sig_label <- sig_stars(merged$fisher_p)
  merged$sub_analysis <- factor(merged$sub_analysis,
                                levels = SUB_ANALYSIS_LEVELS)
  merged$region_set_a <- factor(merged$region_set_a, levels = DIRECTION_SETS)
  merged
}

# =============================================================================
# FIGURES
# =============================================================================

#' Crosswise association heatmap for one sub-analysis.
plot_crosswise_heatmap <- function(cw, title, subtitle) {
  plotCrosswiseMatrix(cw, matrix_type = "association") +
    labs(title = title, subtitle = subtitle) +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

#' Place the annotation column of a forest plot to the right of the points.
forest_label_positions <- function(values) {
  span <- range(c(0, values), na.rm = TRUE)
  pad <- max(diff(span), 1) * 0.10
  list(label_x = span[2] + pad, axis_max = span[2] + pad * 7)
}

#' Permutation z-score per cell, annotated with the plain Fisher odds ratio.
plot_forest <- function(merged, sub_analysis, ntimes) {
  d <- merged[merged$sub_analysis == sub_analysis, , drop = FALSE]
  if (nrow(d) == 0) stop("No merged result for sub-analysis ", sub_analysis)

  positions <- forest_label_positions(d$perm_norm_zscore)
  d$label_x <- positions$label_x
  d$row_label <- factor(d$test_id,
                        levels = d$test_id[order(d$perm_norm_zscore)])
  d$point_label <- sprintf(
    "z = %.2f %s | Fisher OR = %.2f %s | %s/%s genes",
    d$perm_norm_zscore, d$perm_sig_label, d$fisher_or, d$fisher_sig_label,
    vapply(d$n_direction_in_target, fmt_int, character(1)),
    vapply(d$n_direction, fmt_int, character(1)))

  ggplot(d, aes(x = perm_norm_zscore, y = row_label, colour = region_set_a)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = c(-Z_CUTOFF, Z_CUTOFF), linetype = "dotted",
               colour = "grey70") +
    geom_point(size = 3.4) +
    geom_text(aes(x = label_x, label = point_label),
              hjust = 0, size = 3, colour = "grey20") +
    expand_limits(x = positions$axis_max) +
    scale_colour_manual(values = DIRECTION_SET_COLORS, drop = FALSE) +
    labs(
      title = sprintf("Permutation association: mCH direction against %s domains",
                      sub_analysis),
      subtitle = sprintf(paste("%s %s permutations, per.chromosome = %s.",
                               "Dotted lines mark z = +/- %.2f."),
                         fmt_int(ntimes), PERM_RANFUN, PERM_PER_CHROMOSOME,
                         Z_CUTOFF),
      x = "Permutation normalised z-score", y = NULL,
      colour = "mCH direction set",
      caption = paste("Labels give the permutation z-score, the plain gene-level",
                      "Fisher odds ratio, and the count of direction genes",
                      "overlapping the region set.")
    ) +
    theme_emseq() +
    theme(legend.position = "top",
          plot.caption = element_text(hjust = 0, size = 8),
          axis.text.y = element_text(size = 9))
}

#' Fisher odds ratio against permutation z-score, coloured by concordance.
plot_fisher_vs_permutation <- function(merged, ntimes) {
  d <- merged
  positions <- forest_label_positions(d$perm_norm_zscore)
  d$label_x <- positions$label_x
  d$row_label <- factor(d$test_id,
                        levels = d$test_id[order(d$perm_norm_zscore)])
  d$point_label <- sprintf(
    "OR = %.2f (%s) | perm %s",
    d$fisher_or,
    vapply(d$fisher_p, fmt_p, character(1)),
    vapply(d$perm_adj_p_value, fmt_p, character(1)))

  concordance_counts <- table(d$concordance)
  subtitle <- paste(sprintf("%s %d", names(concordance_counts),
                            as.integer(concordance_counts)),
                    collapse = " | ")

  ggplot(d, aes(x = perm_norm_zscore, y = row_label, colour = concordance)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = c(-Z_CUTOFF, Z_CUTOFF), linetype = "dotted",
               colour = "grey70") +
    geom_point(size = 3.4) +
    geom_text(aes(x = label_x, label = point_label),
              hjust = 0, size = 2.9, colour = "grey20") +
    expand_limits(x = positions$axis_max) +
    facet_grid(sub_analysis ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(values = CONCORDANCE_COLORS, drop = FALSE) +
    labs(
      title = "Gene-level Fisher test against genomic permutation",
      subtitle = sprintf("%s permutations. Concordance: %s",
                         fmt_int(ntimes), subtitle),
      x = "Permutation normalised z-score", y = NULL,
      colour = "Concordance",
      caption = sprintf(paste(
        "Confirmed: both significant. Weakened: Fisher only.",
        "Strengthened: permutation only. Concordant NS: neither.",
        "Fisher uses its raw p, permutation its adjusted p, both at %.2f."),
        Q_THRESHOLD)
    ) +
    theme_emseq() +
    theme(legend.position = "top",
          strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
          plot.caption = element_text(hjust = 0, size = 8),
          axis.text.y = element_text(size = 9))
}

#' Observed overlap count against the permuted mean, one panel per region set.
plot_observed_expected <- function(perm_all, ntimes) {
  d <- perm_all
  d$region_set_a <- factor(d$region_set_a, levels = rev(DIRECTION_SETS))
  d$count_label <- sprintf(
    "obs %s\nexp %s +/- %s",
    vapply(d$observed_hits, fmt_int, character(1)),
    vapply(round(d$expected_hits), fmt_int, character(1)),
    vapply(round(d$sd_perm_hits), fmt_int, character(1)))

  ggplot(d, aes(y = region_set_a)) +
    geom_errorbarh(aes(xmin = expected_hits - sd_perm_hits,
                       xmax = expected_hits + sd_perm_hits),
                   height = 0.18, colour = "grey45", linewidth = 0.6) +
    geom_point(aes(x = expected_hits, shape = "Permuted mean"),
               size = 2.6, colour = "grey30") +
    geom_point(aes(x = observed_hits, colour = region_set_a,
                   shape = "Observed"), size = 3.2) +
    geom_text(aes(x = observed_hits, label = count_label),
              vjust = -0.7, size = 2.6, lineheight = 0.95, colour = "grey20") +
    facet_wrap(~ region_set_b, scales = "free_x", ncol = 2) +
    scale_shape_manual(values = c("Observed" = 16, "Permuted mean" = 17),
                       name = NULL) +
    scale_colour_manual(values = DIRECTION_SET_COLORS, guide = "none") +
    scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
    labs(
      title = "Observed overlap against the permutation null",
      subtitle = sprintf(paste("%s %s permutations.",
                               "Bars give one standard deviation of the",
                               "permuted overlap count."),
                         fmt_int(ntimes), PERM_RANFUN),
      x = "Overlap count (numOverlaps)", y = NULL
    ) +
    theme_emseq() +
    theme(legend.position = "top")
}

#' Local z-score profile of one direction set across the region sets.
plot_local_zscore <- function(mlz, region_set_names, title, subtitle) {
  plotSingleLZ(mLZ = mlz, RS = region_set_names, smoothing = TRUE) +
    labs(title = title, subtitle = subtitle) +
    theme_emseq()
}

# =============================================================================
# SUPPORTING TABLES
# =============================================================================

#' Region set sizes for both the Alist and the two Blists.
build_region_set_summary <- function(direction_sets, compartment_sets,
                                     polycomb_sets) {
  rows <- c(
    lapply(names(direction_sets), function(nm) {
      region_set_stats(direction_sets[[nm]], nm, "Both", "Alist")
    }),
    lapply(names(compartment_sets), function(nm) {
      region_set_stats(compartment_sets[[nm]], nm, "Compartment", "Blist")
    }),
    lapply(names(polycomb_sets), function(nm) {
      region_set_stats(polycomb_sets[[nm]], nm, "Polycomb", "Blist")
    })
  )
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Domain membership of the project key genes.
build_key_gene_membership <- function(genes, membership) {
  keep <- genes$gene_name %in% KEY_GENES
  if (!any(keep)) {
    stop("None of the KEY_GENES are in the tested gene universe: ",
         paste(KEY_GENES, collapse = ", "))
  }

  out <- cbind(
    genes[keep, c("gene_name", "gene_id", "chr", "start", "end", "gene_length",
                  "mch_diff", "edger_logFC", "edger_fdr", "mch_sig",
                  "mch_direction"), drop = FALSE],
    membership[keep, , drop = FALSE]
  )
  rownames(out) <- NULL

  cat(sprintf("  Key genes in the universe: %d of %d\n",
              nrow(out), length(KEY_GENES)))
  out
}

#' One row recording every parameter the permutation ran with.
build_parameter_table <- function(opt) {
  data.frame(
    section = SECTION_ID,
    ntimes = opt$ntimes,
    seed = PERM_SEED,
    ran_fun = PERM_RANFUN,
    ev_fun = PERM_EVFUN,
    per_chromosome = PERM_PER_CHROMOSOME,
    cores = opt$cores,
    chromosomes = paste(PERM_CHRS, collapse = ","),
    fdr_threshold = opt$fdr_threshold,
    shift_fdr = opt$shift_fdr,
    shift_diff = opt$shift_diff,
    k27me3_merge_gap_bp = opt$k27me3_merge_gap,
    lz_window_bp = opt$lz_window,
    lz_step_bp = opt$lz_step,
    force_rerun = opt$force_rerun,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

print_cell_summary <- function(merged) {
  cat(sprintf("  %-18s %-26s %10s %10s %12s %12s %-14s\n",
              "mCH set", "Region set", "perm z", "Fisher OR", "perm adj p",
              "Fisher p", "Concordance"))
  ordered <- merged[order(merged$sub_analysis, -merged$perm_norm_zscore), ]
  for (i in seq_len(nrow(ordered))) {
    row <- ordered[i, ]
    cat(sprintf("  %-18s %-26s %10.2f %10.3f %12.3g %12.3g %-14s\n",
                as.character(row$region_set_a), row$region_set_b,
                row$perm_norm_zscore, row$fisher_or,
                row$perm_adj_p_value, row$fisher_p,
                as.character(row$concordance)))
  }
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opt <- parse_options()
  out_dir <- opt$output_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 40_03: DOMAIN-LEVEL PERMUTATION (COMPARTMENTS AND POLYCOMB)\n")
  cat("================================================================================\n")
  cat(sprintf("Output dir:      %s\n", out_dir))
  cat(sprintf("Permutations:    %s (seed %d, %d cores)\n",
              fmt_int(opt$ntimes), PERM_SEED, opt$cores))
  cat(sprintf("Randomisation:   %s, per.chromosome = %s, evFUN = %s\n",
              PERM_RANFUN, PERM_PER_CHROMOSOME, PERM_EVFUN))
  cat(sprintf("Chromosomes:     %s\n", paste(PERM_CHRS, collapse = " ")))
  cat(sprintf("FORCE_RERUN:     %s\n", opt$force_rerun))
  cat("\n")

  # --- Step 1: genome and gene universe --------------------------------------
  cat("STEP 1: Genome and gene universe\n")
  genome <- build_permutation_genome()
  genes <- build_gene_universe(mch_results, opt$fdr_threshold)
  gene_gr <- gene_universe_granges(genes)
  direction_sets <- build_direction_region_sets(genes, gene_gr)
  cat("\n")

  # --- Step 2: compartment region sets ---------------------------------------
  cat("STEP 2: A/B compartment domains\n")
  bins <- load_compartment_bins(HIC_PATHS$compartments,
                                opt$shift_fdr, opt$shift_diff)
  compartment_sets <- build_compartment_region_sets(bins)
  bin_summary <- compartment_bin_summary(bins)
  cat("\n")

  # --- Step 3: Polycomb region sets ------------------------------------------
  cat("STEP 3: Polycomb domains\n")
  polycomb <- build_polycomb_region_sets(CHIP_PATHS$h3k27me3,
                                         opt$k27me3_merge_gap,
                                         HANDOFF_PATHS$chromatin_state)
  polycomb_sets <- polycomb$sets
  cat("\n")

  # --- Step 4: permutation tests ---------------------------------------------
  cat("STEP 4: Crosswise permutation tests\n")
  cached_compartment <- run_crosswise(
    direction_sets, compartment_sets, genome, opt,
    file.path(out_dir, "40_03_perm_compartment.rds"), "Compartment")
  cached_polycomb <- run_crosswise(
    direction_sets, polycomb_sets, genome, opt,
    file.path(out_dir, "40_03_perm_polycomb.rds"), "Polycomb")

  cw_compartment <- cached_compartment$object
  cw_polycomb <- cached_polycomb$object

  perm_compartment <- extract_association_table(cw_compartment, "Compartment")
  perm_polycomb <- extract_association_table(cw_polycomb, "Polycomb")
  perm_all <- rbind(perm_compartment, perm_polycomb)
  cat(sprintf("  Extracted %d permutation cells (%d compartment, %d Polycomb)\n",
              nrow(perm_all), nrow(perm_compartment), nrow(perm_polycomb)))
  cat("\n")

  # --- Step 5: gene-level Fisher tests ---------------------------------------
  cat("STEP 5: Gene-level Fisher tests on the same intervals\n")
  membership_compartment <- compute_membership(gene_gr, compartment_sets)
  membership_polycomb <- compute_membership(gene_gr, polycomb_sets)
  membership <- cbind(membership_compartment, membership_polycomb)

  fisher_compartment <- run_domain_fisher_tests(
    genes, membership_compartment, compartment_sets, "Compartment", out_dir)
  fisher_polycomb <- run_domain_fisher_tests(
    genes, membership_polycomb, polycomb_sets, "Polycomb", out_dir)
  fisher_all <- rbind(fisher_compartment, fisher_polycomb)

  merged <- merge_fisher_permutation(fisher_all, perm_all)
  cat("\n")

  # --- Step 6: local z-score profiles ----------------------------------------
  cat("STEP 6: Local z-score profiles for the mCH hypermethylated set\n")
  cached_lz_compartment <- run_local_zscore(
    direction_sets[[HYPER_SET]], compartment_sets, genome, opt,
    file.path(out_dir, "40_03_localz_compartment.rds"), "Compartment")
  cached_lz_polycomb <- run_local_zscore(
    direction_sets[[HYPER_SET]], polycomb_sets, genome, opt,
    file.path(out_dir, "40_03_localz_polycomb.rds"), "Polycomb")

  mlz_compartment <- cached_lz_compartment$object
  mlz_polycomb <- cached_lz_polycomb$object

  # Every table, figure subtitle, and parameter row below names opt$ntimes, so
  # it has to be the count the four cached objects were actually built at.
  opt <- reconcile_cached_ntimes(opt, list(cached_compartment, cached_polycomb,
                                           cached_lz_compartment,
                                           cached_lz_polycomb))
  cat("\n")

  # --- Step 7: tables --------------------------------------------------------
  cat("STEP 7: Tables\n")
  write_section_table(build_region_set_summary(direction_sets, compartment_sets,
                                           polycomb_sets),
                  file.path(out_dir, "40_03_region_set_summary.tsv"))
  write_section_table(bin_summary, file.path(out_dir, "40_03_compartment_bin_summary.tsv"))
  write_section_table(polycomb$summary, file.path(out_dir, "40_03_polycomb_domain_summary.tsv"))
  write_section_table(perm_compartment, file.path(out_dir, "40_03_association_compartment.tsv"))
  write_section_table(perm_polycomb, file.path(out_dir, "40_03_association_polycomb.tsv"))
  write_section_table(fisher_all, file.path(out_dir, "40_03_fisher_tests.tsv"))
  write_section_table(merged, file.path(out_dir, "40_03_fisher_vs_permutation.tsv"))
  write_section_table(local_zscore_table(mlz_compartment, "Compartment"),
                  file.path(out_dir, "40_03_local_zscore_compartment.tsv"))
  write_section_table(local_zscore_table(mlz_polycomb, "Polycomb"),
                  file.path(out_dir, "40_03_local_zscore_polycomb.tsv"))
  write_section_table(build_key_gene_membership(genes, membership),
                  file.path(out_dir, "40_03_key_genes_domain_membership.tsv"))
  write_section_table(build_parameter_table(opt),
                  file.path(out_dir, "40_03_permutation_parameters.tsv"))
  cat("\n")

  # --- Step 8: figures -------------------------------------------------------
  cat("STEP 8: Figures\n")
  perm_subtitle <- sprintf(
    "%s, %s permutations, per.chromosome = %s",
    PERM_RANFUN, fmt_int(opt$ntimes), PERM_PER_CHROMOSOME)

  save_multiformat_ggplot(
    plot_crosswise_heatmap(
      cw_compartment,
      "Permutation association: mCH direction against A/B compartments",
      perm_subtitle),
    file.path(out_dir, "40_03a_crosswise_mch_x_compartment"),
    width = 10, height = 6)

  save_multiformat_ggplot(
    plot_forest(merged, "Compartment", opt$ntimes),
    file.path(out_dir, "40_03b_forest_mch_x_compartment"),
    width = 15, height = 7)

  save_multiformat_ggplot(
    plot_crosswise_heatmap(
      cw_polycomb,
      "Permutation association: mCH direction against Polycomb domains",
      perm_subtitle),
    file.path(out_dir, "40_03c_crosswise_mch_x_polycomb"),
    width = 9, height = 6)

  save_multiformat_ggplot(
    plot_forest(merged, "Polycomb", opt$ntimes),
    file.path(out_dir, "40_03d_forest_mch_x_polycomb"),
    width = 15, height = 6)

  save_multiformat_ggplot(
    plot_fisher_vs_permutation(merged, opt$ntimes),
    file.path(out_dir, "40_03e_fisher_vs_permutation"),
    width = 15, height = 9)

  save_multiformat_ggplot(
    plot_observed_expected(perm_all, opt$ntimes),
    file.path(out_dir, "40_03f_observed_vs_expected_overlaps"),
    width = 13, height = 10)

  save_multiformat_ggplot(
    plot_local_zscore(
      mlz_compartment, names(compartment_sets),
      "Local z-score: mCH hypermethylated genes at compartment domains",
      sprintf(paste("Z-score at each offset of the gene bodies, +/- %s bp in",
                    "%s bp steps, %s permutations per offset."),
              fmt_int(opt$lz_window), fmt_int(opt$lz_step),
              fmt_int(opt$ntimes))),
    file.path(out_dir, "40_03g_local_zscore_compartment"),
    width = 11, height = 6)

  save_multiformat_ggplot(
    plot_local_zscore(
      mlz_polycomb, names(polycomb_sets),
      "Local z-score: mCH hypermethylated genes at Polycomb domains",
      sprintf(paste("Z-score at each offset of the gene bodies, +/- %s bp in",
                    "%s bp steps, %s permutations per offset."),
              fmt_int(opt$lz_window), fmt_int(opt$lz_step),
              fmt_int(opt$ntimes))),
    file.path(out_dir, "40_03h_local_zscore_polycomb"),
    width = 11, height = 6)
  cat("\n")

  # --- Step 9: summary -------------------------------------------------------
  cat("================================================================================\n")
  cat("SECTION 40_03 SUMMARY\n")
  cat("================================================================================\n")
  cat(sprintf("Gene universe: %s genes (%s hyper, %s hypo) at FDR < %.3f\n",
              fmt_int(nrow(genes)), fmt_int(sum(genes$mch_hyper)),
              fmt_int(sum(genes$mch_hypo)), opt$fdr_threshold))
  cat(sprintf("Compartment sub-analysis: %d x %d = %d cells\n",
              length(direction_sets), length(compartment_sets),
              length(direction_sets) * length(compartment_sets)))
  cat(sprintf("Polycomb sub-analysis:    %d x %d = %d cells\n",
              length(direction_sets), length(polycomb_sets),
              length(direction_sets) * length(polycomb_sets)))
  cat(sprintf("Fisher tests registered:  %d\n", nrow(fisher_all)))
  cat("\n")
  print_cell_summary(merged)

  cat("\nConcordance between Fisher and permutation:\n")
  for (level in CONCORDANCE_LEVELS) {
    n <- sum(merged$concordance == level)
    cat(sprintf("  %-16s %2d of %d (%.0f%%)\n", level, n, nrow(merged),
                100 * n / nrow(merged)))
  }

  cat(sprintf("\nFigures and tables written to: %s\n", out_dir))
  cat(sprintf("Fisher tests registered in:    %s\n",
              HANDOFF_PATHS$fisher_registry))
  cat("Section 40_03 complete.\n\n")
}

main()
