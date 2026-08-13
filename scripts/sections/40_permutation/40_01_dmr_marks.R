# scripts/sections/40_permutation/40_01_dmr_marks.R
#
# Section 40_01: genomic permutation test of mCH gene bodies against chromatin
# mark peak sets.
#
# What this tests
#   Gene-level Fisher tests ask whether mCH-significant genes overlap a peak set
#   more often than the other tested genes. That question ignores where the
#   regions sit in the genome. Gene bodies cluster, peaks cluster, and both
#   follow gene density, so a Fisher odds ratio can be large without any spatial
#   association.
#
#   This section answers the same question with a genomic null. Each mCH
#   direction set is randomly relocated within its own chromosome many times.
#   The observed overlap count with each peak set is compared against the
#   distribution of overlap counts from the relocated sets. A z-score far from
#   zero means the co-occurrence is more than the genomic arrangement alone
#   produces.
#
#   Two permutation runs are executed on the same region sets:
#     crosswisePermTest  regioneReloaded, all Alist x Blist cells at once.
#                        It keeps the z-score, the p-value, the observed count,
#                        and the mean and standard deviation of the null. It
#                        discards the individual permuted counts.
#     permTest           regioneR, one call per Alist set over the same
#                        evaluation function list. It keeps every permuted
#                        count, which the null-distribution figure draws.
#   Both runs report an observed overlap count. The script stops if the two
#   counts disagree for any cell.
#
# Reads (all through _shared_config.R)
#   gene_bodies                      tested gene universe as GRanges
#   mch_results                      gene-level mCH differential results
#   DIFFBIND_PATHS$atac_up/atac_down ATAC differential peak BED files
#   k119ub_diffbind, k27ac_diffbind  DiffBind tables, split by direction
#   MECP2_PATHS$up, MECP2_PATHS$down MeCP2 differential peak BED files
#   BSgenome.Mmusculus.UCSC.mm10     chromosome lengths for the randomisation
#
# Writes (OUT_DIR defaults to results/sections/40_permutation/)
#   Figures 40_01a .. 40_01e, each in a multi-format subdirectory
#   Tables prefixed 40_01_
#   The cached permutation objects as 40_01_permutation_objects.rds
#   Fisher gene tables under fisher_tables/ and rows in the shared registry,
#   which section 40_04 validates by label permutation
#
# Caching
#   The RDS cache is loaded instead of recomputed when it exists and the
#   FORCE_RERUN environment variable is empty. Figures and tables are rebuilt on
#   every run from whichever object is in memory.
#
# Adapted from Biomodal section 34
# (section_34_permutation_dmr_chromatin_marks.R). That section permuted mC DMR
# intervals; this port permutes mCH-significant gene bodies. Its condition
# specific peak sets (K119ub Ctrl/Mut, H3K27ac Ctrl/Mut) are replaced by
# differential sets, and MeCP2 up and down peaks are added.

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)
library(regioneR)
library(regioneReloaded)
library(BSgenome.Mmusculus.UCSC.mm10)

patch_chooseHclustMet()

# =============================================================================
# CONSTANTS
# =============================================================================

SECTION_ID <- "40_01"

OUT_DIR <- OUTPUT_PATHS$permutation

# Randomisation settings. Every region keeps its own chromosome and its own
# width; only its position moves.
PERM_SEED     <- 42L
PERM_PER_CHR  <- TRUE
PERM_RANFUN   <- "randomizeRegions"
PERM_EVFUN    <- "numOverlaps"
PERM_DEFAULT_NTIMES <- 5000L
PERM_DEFAULT_CORES  <- 8L

# The autosomes plus chrX. chrY is excluded because half the cohort is female.
STANDARD_CHRS <- paste0("chr", c(1:19, "X"))

CACHE_FILENAME <- "40_01_permutation_objects.rds"

# --- Region set naming, ordering, and colour --------------------------------

ALIST_ORDER <- c("mCH Hyper gene bodies", "mCH Hypo gene bodies")

ALIST_SHORT <- c(
  "mCH Hyper gene bodies" = "Hyper",
  "mCH Hypo gene bodies"  = "Hypo"
)

ALIST_COLORS <- c(
  "mCH Hyper gene bodies" = unname(COLORS$direction[["Hypermethylated"]]),
  "mCH Hypo gene bodies"  = unname(COLORS$direction[["Hypomethylated"]])
)

BLIST_ORDER <- c("ATAC Up", "ATAC Down",
                 "K119ub Gained", "K119ub Lost",
                 "K27ac Gained", "K27ac Lost",
                 "MeCP2 Up", "MeCP2 Down")

BLIST_GROUPS <- c(
  "ATAC Up"       = "ATAC",
  "ATAC Down"     = "ATAC",
  "K119ub Gained" = "H2AK119ub",
  "K119ub Lost"   = "H2AK119ub",
  "K27ac Gained"  = "H3K27ac",
  "K27ac Lost"    = "H3K27ac",
  "MeCP2 Up"      = "MeCP2",
  "MeCP2 Down"    = "MeCP2"
)

BLIST_GROUP_ORDER <- c("ATAC", "H2AK119ub", "H3K27ac", "MeCP2")

BLIST_COLORS <- c(
  "ATAC Up"       = unname(COLORS$atac[["ATAC Up"]]),
  "ATAC Down"     = unname(COLORS$atac[["ATAC Down"]]),
  "K119ub Gained" = unname(COLORS$k119ub[["K119ub Gained"]]),
  "K119ub Lost"   = unname(COLORS$k119ub[["K119ub Lost"]]),
  "K27ac Gained"  = unname(COLORS$h3k27ac[["H3K27ac Gained"]]),
  "K27ac Lost"    = unname(COLORS$h3k27ac[["H3K27ac Lost"]]),
  "MeCP2 Up"      = unname(COLORS$mecp2[["MeCP2 Up"]]),
  "MeCP2 Down"    = unname(COLORS$mecp2[["MeCP2 Down"]])
)

CONCORDANCE_ORDER <- c("Confirmed", "Weakened", "Strengthened", "Concordant NS")

CONCORDANCE_COLORS <- c(
  "Confirmed"     = "#2CA02C",
  "Weakened"      = "#D62728",
  "Strengthened"  = "#FF7F0E",
  "Concordant NS" = "grey60"
)

# Smallest p-value the z-score conversion accepts, so qnorm stays finite.
MIN_P_FOR_Z <- 1e-300

# =============================================================================
# COMMAND LINE
# =============================================================================

parse_section_options <- function() {
  option_list <- list(
    make_option("--output-dir", dest = "output_dir", type = "character",
                default = OUT_DIR,
                help = "Directory for figures, tables, and the RDS cache [default: %default]"),
    make_option("--ntimes", dest = "ntimes", type = "integer",
                default = PERM_DEFAULT_NTIMES,
                help = "Number of randomisations per region set [default: %default]")
  )

  parser <- OptionParser(
    option_list = option_list,
    usage = "Rscript 40_01_dmr_marks.R [--output-dir DIR] [--ntimes N] [N]"
  )
  parsed <- parse_args(parser, positional_arguments = TRUE)
  opt <- parsed$options

  # A bare positional number sets ntimes, matching the Biomodal call style.
  if (length(parsed$args) > 1) {
    stop("At most one positional argument is accepted (ntimes). Got: ",
         paste(parsed$args, collapse = ", "))
  }
  if (length(parsed$args) == 1) {
    positional_ntimes <- suppressWarnings(as.integer(parsed$args[1]))
    if (is.na(positional_ntimes)) {
      stop("The positional argument must be an integer number of ",
           "randomisations. Got: ", parsed$args[1])
    }
    opt$ntimes <- positional_ntimes
  }

  if (opt$ntimes < 100) {
    stop("--ntimes must be at least 100, got ", opt$ntimes)
  }
  opt
}

#' Number of cores for the permutation, from the SLURM allocation.
resolve_cores <- function() {
  slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
  if (!nzchar(slurm_cpus)) return(PERM_DEFAULT_CORES)

  cores <- suppressWarnings(as.integer(slurm_cpus))
  if (is.na(cores) || cores < 1) {
    stop("SLURM_CPUS_PER_TASK is set to '", slurm_cpus,
         "' which is not a positive integer.")
  }
  cores
}

force_rerun_requested <- function() {
  nzchar(Sys.getenv("FORCE_RERUN", unset = ""))
}

# =============================================================================
# SMALL UTILITIES
# =============================================================================

fmt_p <- function(p) {
  if (length(p) != 1 || is.na(p)) return("p = NA")
  if (p == 0) return("p < 1/ntimes")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

#' Round a count and group its thousands, without padding to a common width.
fmt_count <- function(x) {
  trimws(format(round(x), big.mark = ",", scientific = FALSE))
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

#' Return the one column name from candidates that df carries.
#'
#' Stops, and prints every column of df, when none or more than one candidate is
#' present.
resolve_column <- function(df, candidates, purpose) {
  present <- candidates[candidates %in% colnames(df)]
  if (length(present) != 1) {
    stop("Cannot resolve the ", purpose, " column of the permutation results. ",
         "Looked for: ", paste(candidates, collapse = ", "),
         ". Columns present: ", paste(colnames(df), collapse = ", "))
  }
  present
}

#' Convert a two-sided p-value and an odds ratio to a signed z-score.
#'
#' qnorm is called on the lower tail so that very small p-values stay finite.
fisher_p_to_signed_z <- function(p_value, odds_ratio) {
  magnitude <- qnorm(pmax(p_value, MIN_P_FOR_Z) / 2, lower.tail = FALSE)
  magnitude * sign(log2(odds_ratio))
}

#' Monte Carlo standard error of a permutation z-score.
#'
#' The z-score divides a fixed observed count by an estimated null mean and
#' standard deviation, both from ntimes draws. The delta method gives
#' se(z) = sqrt(1 / ntimes + z^2 / (2 * ntimes)).
zscore_standard_error <- function(z_score, ntimes) {
  sqrt(1 / ntimes + z_score^2 / (2 * ntimes))
}

overlap_column <- function(set_name) {
  paste0("ov_", gsub("[^A-Za-z0-9]+", "_", tolower(set_name)))
}

test_key <- function(alist_set, blist_set) {
  paste(alist_set, "x", blist_set)
}

# =============================================================================
# GENOME
# =============================================================================

#' mm10 chromosome ranges for chr1-chr19 and chrX.
build_permutation_genome <- function() {
  cat("Constructing the mm10 genome object...\n")
  genome_full <- getGenomeAndMask("mm10")$genome
  genome <- genome_full[as.character(seqnames(genome_full)) %in% STANDARD_CHRS]

  if (length(genome) != length(STANDARD_CHRS)) {
    stop("Expected ", length(STANDARD_CHRS), " chromosomes in the mm10 genome ",
         "object, got ", length(genome), ": ",
         paste(as.character(seqnames(genome)), collapse = ", "))
  }

  cat(sprintf("  Genome: %d chromosomes, %.2f Gb total\n",
              length(genome), sum(as.numeric(width(genome))) / 1e9))
  genome
}

#' Keep the standard chromosomes and give the GRanges the genome's seqlevels.
restrict_to_standard <- function(gr, genome, set_name) {
  n_before <- length(gr)
  gr <- gr[as.character(seqnames(gr)) %in% STANDARD_CHRS]
  if (length(gr) == 0) {
    stop("Region set '", set_name, "' has no regions on ",
         paste(STANDARD_CHRS, collapse = ", "))
  }
  seqlevels(gr) <- seqlevels(genome)

  cat(sprintf("  %-24s %7d regions (%d dropped off chr1-19/X), median width %d bp\n",
              set_name, length(gr), n_before - length(gr),
              as.integer(median(width(gr)))))
  gr
}

# =============================================================================
# REGION SETS
# =============================================================================

#' The two mCH direction sets, taken from the tested gene universe.
build_alist <- function(genome) {
  cat("\nConstructing Alist (mCH direction gene body sets)...\n")

  is_hyper <- gene_bodies$mch_sig & gene_bodies$edger_logFC > 0
  is_hypo  <- gene_bodies$mch_sig & gene_bodies$edger_logFC < 0

  if (sum(is_hyper) != sum(mch_results$mch_hyper) ||
      sum(is_hypo) != sum(mch_results$mch_hypo)) {
    stop("Gene body direction counts disagree with mch_results. ",
         "gene_bodies: ", sum(is_hyper), " hyper / ", sum(is_hypo), " hypo. ",
         "mch_results: ", sum(mch_results$mch_hyper), " hyper / ",
         sum(mch_results$mch_hypo), " hypo.")
  }

  alist <- list(
    "mCH Hyper gene bodies" = restrict_to_standard(gene_bodies[is_hyper], genome,
                                                   "mCH Hyper gene bodies"),
    "mCH Hypo gene bodies"  = restrict_to_standard(gene_bodies[is_hypo], genome,
                                                   "mCH Hypo gene bodies")
  )
  alist[ALIST_ORDER]
}

#' Split a DiffBind table into its gained and lost peak sets.
diffbind_direction_granges <- function(db, mark_name) {
  gr <- diffbind_to_granges(db)
  gained <- gr[gr$direction == "Gained"]
  lost   <- gr[gr$direction == "Lost"]

  if (length(gained) == 0 || length(lost) == 0) {
    stop(mark_name, " DiffBind table has ", length(gained), " gained and ",
         length(lost), " lost peaks. Both sets must be non-empty.")
  }
  list(gained = gained, lost = lost)
}

#' The eight chromatin mark peak sets tested against the mCH gene bodies.
build_blist <- function(genome) {
  cat("\nConstructing Blist (chromatin mark peak sets)...\n")

  atac_up   <- load_chip_peaks(DIFFBIND_PATHS$atac_up, "ATAC Up")
  atac_down <- load_chip_peaks(DIFFBIND_PATHS$atac_down, "ATAC Down")
  mecp2_up   <- load_chip_peaks(MECP2_PATHS$up, "MeCP2 Up")
  mecp2_down <- load_chip_peaks(MECP2_PATHS$down, "MeCP2 Down")

  k119ub <- diffbind_direction_granges(k119ub_diffbind, "H2AK119ub")
  k27ac  <- diffbind_direction_granges(k27ac_diffbind, "H3K27ac")

  raw <- list(
    "ATAC Up"       = atac_up,
    "ATAC Down"     = atac_down,
    "K119ub Gained" = k119ub$gained,
    "K119ub Lost"   = k119ub$lost,
    "K27ac Gained"  = k27ac$gained,
    "K27ac Lost"    = k27ac$lost,
    "MeCP2 Up"      = mecp2_up,
    "MeCP2 Down"    = mecp2_down
  )

  cat("\n  Restricting Blist to chr1-19 and chrX...\n")
  blist <- lapply(names(raw), function(nm) restrict_to_standard(raw[[nm]], genome, nm))
  names(blist) <- names(raw)
  blist[BLIST_ORDER]
}

#' Region counts and widths of every set entering the permutation.
build_region_set_table <- function(alist, blist) {
  sets <- c(alist, blist)
  roles <- c(rep("Alist", length(alist)), rep("Blist", length(blist)))

  do.call(rbind, lapply(seq_along(sets), function(i) {
    gr <- sets[[i]]
    data.frame(
      role = roles[i],
      region_set = names(sets)[i],
      n_regions = length(gr),
      total_bp = sum(as.numeric(width(gr))),
      median_width_bp = as.numeric(median(width(gr))),
      mean_width_bp = mean(as.numeric(width(gr))),
      n_chromosomes = length(unique(as.character(seqnames(gr)))),
      stringsAsFactors = FALSE
    )
  }))
}

# =============================================================================
# PERMUTATION
# =============================================================================

#' Run crosswisePermTest over every Alist x Blist pair.
run_crosswise <- function(alist, blist, genome, ntimes, cores) {
  cat(sprintf("\nRunning crosswisePermTest: %d Alist x %d Blist sets, ntimes=%d, cores=%d\n",
              length(alist), length(blist), ntimes, cores))

  options(mc.cores = cores)
  set.seed(PERM_SEED)
  cw <- crosswisePermTest(
    Alist          = alist,
    Blist          = blist,
    genome         = genome,
    ranFUN         = PERM_RANFUN,
    evFUN          = PERM_EVFUN,
    ntimes         = ntimes,
    force.parallel = TRUE,
    per.chromosome = PERM_PER_CHR
  )

  cat("  Building the crosswise association matrix...\n")
  makeCrosswiseMatrix(cw, pvcut = 1, symm_matrix = FALSE, hc.method = "average")
}

#' Run one regioneR permTest per Alist set, keeping every permuted count.
run_null_distributions <- function(alist, blist, genome, ntimes, cores) {
  cat(sprintf("\nRunning regioneR permTest for the null distributions (ntimes=%d, cores=%d)\n",
              ntimes, cores))

  func_list <- createFunctionsList(FUN = numOverlaps, param.name = "B",
                                   values = blist)

  options(mc.cores = cores)
  set.seed(PERM_SEED)
  out <- lapply(names(alist), function(a_name) {
    cat(sprintf("  %s ...\n", a_name))
    permTest(
      A                  = alist[[a_name]],
      evaluate.function  = func_list,
      randomize.function = randomizeRegions,
      genome             = genome,
      ntimes             = ntimes,
      force.parallel     = TRUE,
      per.chromosome     = PERM_PER_CHR
    )
  })
  names(out) <- names(alist)
  out
}

#' Load the cached permutation objects, or compute and cache them.
load_or_run_permutations <- function(alist, blist, genome, ntimes, cores, cache_path) {
  if (file.exists(cache_path) && !force_rerun_requested()) {
    cat("\nLoading cached permutation objects from: ", cache_path, "\n", sep = "")
    cached <- readRDS(cache_path)

    required <- c("crosswise", "nulls", "ntimes")
    missing <- setdiff(required, names(cached))
    if (length(missing) > 0) {
      stop("Cache ", cache_path, " is missing entries: ",
           paste(missing, collapse = ", "),
           ". Set FORCE_RERUN=1 to rebuild it.")
    }
    cat(sprintf("  Cache holds %d randomisations per region set.\n", cached$ntimes))
    return(cached)
  }

  if (force_rerun_requested()) {
    cat("\nFORCE_RERUN is set; recomputing the permutations.\n")
  }

  results <- list(
    crosswise = run_crosswise(alist, blist, genome, ntimes, cores),
    nulls     = run_null_distributions(alist, blist, genome, ntimes, cores),
    ntimes    = ntimes,
    seed      = PERM_SEED,
    ran_fun   = PERM_RANFUN,
    ev_fun    = PERM_EVFUN,
    per_chromosome = PERM_PER_CHR
  )

  saveRDS(results, cache_path)
  cat("  Saved cache: ", cache_path, "\n", sep = "")
  results
}

# =============================================================================
# RESULT EXTRACTION
# =============================================================================

#' One row per Alist x Blist cell from the crosswisePermTest object.
extract_crosswise_table <- function(cw, ntimes) {
  multi_overlaps <- cw@multiOverlaps

  missing_sets <- setdiff(ALIST_ORDER, names(multi_overlaps))
  if (length(missing_sets) > 0) {
    stop("crosswisePermTest results are missing Alist sets: ",
         paste(missing_sets, collapse = ", "))
  }

  rows <- lapply(ALIST_ORDER, function(a_name) {
    df <- multi_overlaps[[a_name]]

    hits_col <- resolve_column(df, c("n_hits", "n_overlaps"), "observed overlap")
    adj_col  <- resolve_column(df, c("adj_p_value", "adj.p_value"),
                               "adjusted p-value")
    for (col in c("name", "n_regionA", "n_regionB", "z_score", "p_value",
                  "mean_perm_test", "sd_perm_test", "norm_zscore")) {
      if (!col %in% colnames(df)) {
        stop("crosswisePermTest results for '", a_name, "' have no '", col,
             "' column. Columns present: ",
             paste(colnames(df), collapse = ", "))
      }
    }

    data.frame(
      alist_set          = a_name,
      blist_set          = as.character(df$name),
      n_region_a         = as.numeric(df$n_regionA),
      n_region_b         = as.numeric(df$n_regionB),
      observed_overlaps  = as.numeric(df[[hits_col]]),
      null_mean          = as.numeric(df$mean_perm_test),
      null_sd            = as.numeric(df$sd_perm_test),
      z_score            = as.numeric(df$z_score),
      norm_zscore        = as.numeric(df$norm_zscore),
      p_value            = as.numeric(df$p_value),
      adj_p_value        = as.numeric(df[[adj_col]]),
      stringsAsFactors   = FALSE
    )
  })

  out <- do.call(rbind, rows)

  unknown <- setdiff(out$blist_set, BLIST_ORDER)
  if (length(unknown) > 0) {
    stop("crosswisePermTest returned unexpected Blist names: ",
         paste(unknown, collapse = ", "))
  }

  out$mark_group <- unname(BLIST_GROUPS[out$blist_set])
  out$alist_label <- unname(ALIST_SHORT[out$alist_set])
  out$test_id <- test_key(out$alist_set, out$blist_set)
  out$observed_over_expected <- out$observed_overlaps / out$null_mean
  out$se_z <- zscore_standard_error(out$z_score, ntimes)
  out$z_ci_lower <- out$z_score - 1.96 * out$se_z
  out$z_ci_upper <- out$z_score + 1.96 * out$se_z
  out$ntimes <- ntimes
  out$sig_label <- sig_stars(out$p_value)

  out <- out[order(match(out$alist_set, ALIST_ORDER),
                   match(out$blist_set, BLIST_ORDER)), , drop = FALSE]
  rownames(out) <- NULL

  cat(sprintf("  Extracted %d pairwise permutation results\n", nrow(out)))
  cat(sprintf("  z-score range: [%.2f, %.2f]\n",
              min(out$z_score), max(out$z_score)))
  out
}

#' Every permuted overlap count in long form, one row per randomisation.
extract_null_draws <- function(nulls) {
  rows <- lapply(ALIST_ORDER, function(a_name) {
    pt <- nulls[[a_name]]
    if (is.null(pt)) {
      stop("Null distribution results are missing Alist set: ", a_name)
    }

    do.call(rbind, lapply(BLIST_ORDER, function(b_name) {
      res <- pt[[b_name]]
      if (is.null(res)) {
        stop("Null distribution results for '", a_name,
             "' are missing Blist set: ", b_name,
             ". Sets present: ", paste(names(pt), collapse = ", "))
      }
      data.frame(
        alist_set = a_name,
        blist_set = b_name,
        permutation = seq_along(res$permuted),
        permuted_overlaps = as.numeric(res$permuted),
        stringsAsFactors = FALSE
      )
    }))
  })

  out <- do.call(rbind, rows)
  out$test_id <- test_key(out$alist_set, out$blist_set)
  out
}

#' One row per cell summarising the retained null draws.
summarise_null_draws <- function(nulls, draws) {
  observed <- do.call(rbind, lapply(ALIST_ORDER, function(a_name) {
    do.call(rbind, lapply(BLIST_ORDER, function(b_name) {
      res <- nulls[[a_name]][[b_name]]
      data.frame(
        alist_set = a_name,
        blist_set = b_name,
        observed_overlaps = as.numeric(res$observed),
        null_run_z = as.numeric(res$zscore),
        null_run_p = as.numeric(res$pval),
        null_run_alternative = as.character(res$alternative),
        stringsAsFactors = FALSE
      )
    }))
  }))
  observed$test_id <- test_key(observed$alist_set, observed$blist_set)

  summaries <- draws %>%
    dplyr::group_by(test_id) %>%
    dplyr::summarise(
      n_permutations = dplyr::n(),
      null_mean_draws = mean(permuted_overlaps),
      null_sd_draws = sd(permuted_overlaps),
      null_median = median(permuted_overlaps),
      null_min = min(permuted_overlaps),
      null_max = max(permuted_overlaps),
      null_q025 = quantile(permuted_overlaps, 0.025, names = FALSE),
      null_q975 = quantile(permuted_overlaps, 0.975, names = FALSE),
      .groups = "drop"
    ) %>%
    as.data.frame()

  out <- dplyr::left_join(observed, summaries, by = "test_id")

  counts <- draws %>%
    dplyr::left_join(observed[, c("test_id", "observed_overlaps")], by = "test_id") %>%
    dplyr::group_by(test_id) %>%
    dplyr::summarise(
      n_null_ge_observed = sum(permuted_overlaps >= observed_overlaps),
      n_null_le_observed = sum(permuted_overlaps <= observed_overlaps),
      .groups = "drop"
    ) %>%
    as.data.frame()

  out <- dplyr::left_join(out, counts, by = "test_id")
  out$p_empirical_greater <- (1 + out$n_null_ge_observed) / (1 + out$n_permutations)
  out$p_empirical_less    <- (1 + out$n_null_le_observed) / (1 + out$n_permutations)
  out$mark_group <- unname(BLIST_GROUPS[out$blist_set])
  out$alist_label <- unname(ALIST_SHORT[out$alist_set])
  out
}

#' Stop when the two permutation runs disagree on an observed overlap count.
#'
#' The observed count does not depend on the randomisation, so the two runs must
#' report the same number for every cell.
check_observed_counts_agree <- function(crosswise, null_summary) {
  merged <- dplyr::inner_join(
    crosswise[, c("test_id", "observed_overlaps")],
    null_summary[, c("test_id", "observed_overlaps")],
    by = "test_id", suffix = c("_crosswise", "_nullrun"))

  if (nrow(merged) != nrow(crosswise)) {
    stop("Only ", nrow(merged), " of ", nrow(crosswise),
         " permutation cells matched between the two runs.")
  }

  mismatched <- merged[merged$observed_overlaps_crosswise !=
                         merged$observed_overlaps_nullrun, , drop = FALSE]
  if (nrow(mismatched) > 0) {
    print(mismatched)
    stop("crosswisePermTest and permTest report different observed overlap ",
         "counts for ", nrow(mismatched), " cell(s). The region sets differ ",
         "between the two runs.")
  }

  cat(sprintf("  Observed overlap counts agree across both runs for all %d cells.\n",
              nrow(merged)))
  invisible(TRUE)
}

#' Join the two runs into one association table.
build_association_table <- function(crosswise, null_summary) {
  keep <- c("test_id", "null_run_z", "null_run_p", "null_run_alternative",
            "n_permutations", "null_mean_draws", "null_sd_draws", "null_median",
            "null_min", "null_max", "null_q025", "null_q975",
            "n_null_ge_observed", "n_null_le_observed",
            "p_empirical_greater", "p_empirical_less")

  out <- dplyr::left_join(crosswise, null_summary[, keep], by = "test_id")
  out$z_difference <- out$z_score - out$null_run_z
  out$ran_fun <- PERM_RANFUN
  out$ev_fun <- PERM_EVFUN
  out$per_chromosome <- PERM_PER_CHR
  out$seed <- PERM_SEED
  out
}

#' The fixed settings of this run, as a one-column-per-setting table.
build_parameter_table <- function(genome, alist, blist, ntimes, cores) {
  data.frame(
    section = SECTION_ID,
    ntimes = ntimes,
    mc_cores = cores,
    seed = PERM_SEED,
    ran_fun = PERM_RANFUN,
    ev_fun = PERM_EVFUN,
    per_chromosome = PERM_PER_CHR,
    mask = "none",
    n_alist_sets = length(alist),
    n_blist_sets = length(blist),
    n_pairwise_tests = length(alist) * length(blist),
    n_chromosomes = length(genome),
    genome_bp = sum(as.numeric(width(genome))),
    chromosomes = paste(STANDARD_CHRS, collapse = ","),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# GENE-LEVEL FISHER TESTS
# =============================================================================

#' One row per gene symbol, with an overlap flag for every Blist set.
#'
#' The permutation universe is the whole genome; the Fisher universe is the set
#' of genes that carry an mCH measurement on chr1-19 or chrX. The two tests
#' therefore ask related but different questions, which is the point of the
#' comparison later in this section.
build_gene_overlap_table <- function(mch, blist) {
  genes <- mch %>%
    dplyr::arrange(gene_name, dplyr::desc(abs(edger_logFC)), edger_fdr) %>%
    dplyr::distinct(gene_name, .keep_all = TRUE) %>%
    as.data.frame()

  cat(sprintf("  mCH universe: %s genes after collapsing %s duplicate symbols\n",
              format(nrow(genes), big.mark = ","),
              format(nrow(mch) - nrow(genes), big.mark = ",")))

  genes <- genes[genes$chr %in% STANDARD_CHRS, , drop = FALSE]
  cat(sprintf("  On chr1-19 and chrX: %s genes (%s hyper, %s hypo)\n",
              format(nrow(genes), big.mark = ","),
              format(sum(genes$mch_hyper), big.mark = ","),
              format(sum(genes$mch_hypo), big.mark = ",")))

  gr <- GRanges(seqnames = genes$chr,
                ranges = IRanges(start = genes$start + 1L, end = genes$end))

  for (set_name in names(blist)) {
    flag <- countOverlaps(gr, blist[[set_name]]) > 0
    genes[[overlap_column(set_name)]] <- flag
    cat(sprintf("    %-16s overlaps %s of %s gene bodies\n", set_name,
                format(sum(flag), big.mark = ","),
                format(nrow(genes), big.mark = ",")))
  }
  genes
}

#' Run one registered gene-level Fisher test and return its counts as a row.
run_one_fisher <- function(gene_df, row_var, col_var, test_id, description,
                           alist_set, blist_set, universe_label, out_dir) {
  ft <- register_fisher_test(
    section = SECTION_ID, test_id = test_id, description = description,
    gene_df = gene_df, row_var = row_var, col_var = col_var,
    output_dir = out_dir)

  row_flag <- gene_df[[row_var]]
  col_flag <- gene_df[[col_var]]

  data.frame(
    section = SECTION_ID,
    test_id = test_id,
    description = description,
    alist_set = alist_set,
    blist_set = blist_set,
    universe = universe_label,
    row_var = row_var,
    col_var = col_var,
    n_genes = nrow(gene_df),
    n_row_true = sum(row_flag),
    n_col_true = sum(col_flag),
    n_both_true = sum(row_flag & col_flag),
    pct_overlap_in_row = 100 * sum(row_flag & col_flag) / sum(row_flag),
    pct_overlap_out_row = 100 * sum(!row_flag & col_flag) / sum(!row_flag),
    odds_ratio = unname(ft$estimate),
    ci_lower = ft$conf.int[1],
    ci_upper = ft$conf.int[2],
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

#' Every gene-level Fisher test of this section.
#'
#' Three families:
#'   hyper_x_<set>          all tested genes, hypermethylated against overlap
#'   hypo_x_<set>           all tested genes, hypomethylated against overlap
#'   hyper_vs_hypo_x_<set>  significant genes only, direction against overlap
run_gene_fisher_tests <- function(genes, blist, out_dir) {
  significant <- genes[genes$mch_sig, , drop = FALSE]
  if (nrow(significant) == 0) {
    stop("No mCH-significant genes remain for the direction contrast tests.")
  }

  rows <- list()
  for (set_name in names(blist)) {
    col_var <- overlap_column(set_name)
    key <- sub("^ov_", "", col_var)

    rows[[length(rows) + 1]] <- run_one_fisher(
      gene_df = genes, row_var = "mch_hyper", col_var = col_var,
      test_id = paste0("hyper_x_", key),
      description = sprintf(
        "Do mCH hypermethylated gene bodies overlap %s peaks more than other tested genes?",
        set_name),
      alist_set = "mCH Hyper gene bodies", blist_set = set_name,
      universe_label = "All tested genes", out_dir = out_dir)

    rows[[length(rows) + 1]] <- run_one_fisher(
      gene_df = genes, row_var = "mch_hypo", col_var = col_var,
      test_id = paste0("hypo_x_", key),
      description = sprintf(
        "Do mCH hypomethylated gene bodies overlap %s peaks more than other tested genes?",
        set_name),
      alist_set = "mCH Hypo gene bodies", blist_set = set_name,
      universe_label = "All tested genes", out_dir = out_dir)

    rows[[length(rows) + 1]] <- run_one_fisher(
      gene_df = significant, row_var = "mch_hyper", col_var = col_var,
      test_id = paste0("hyper_vs_hypo_x_", key),
      description = sprintf(
        "Among mCH-significant genes, do hypermethylated bodies overlap %s peaks more than hypomethylated bodies?",
        set_name),
      alist_set = "Significant genes", blist_set = set_name,
      universe_label = "mCH-significant genes", out_dir = out_dir)
  }

  out <- do.call(rbind, rows)
  out$q_value <- p.adjust(out$p_value, method = "BH")
  out$sig_label <- sig_stars(out$q_value)
  out
}

# =============================================================================
# FISHER VERSUS PERMUTATION
# =============================================================================

classify_concordance <- function(fisher_p, perm_p, threshold = Q_THRESHOLD) {
  out <- dplyr::case_when(
    fisher_p < threshold & perm_p < threshold ~ "Confirmed",
    fisher_p < threshold & perm_p >= threshold ~ "Weakened",
    fisher_p >= threshold & perm_p < threshold ~ "Strengthened",
    TRUE ~ "Concordant NS"
  )
  factor(out, levels = CONCORDANCE_ORDER)
}

#' Match every permutation cell to the gene-level Fisher test of the same pair.
#'
#' The Fisher columns are renamed before the join, so the merged table carries
#' one unprefixed set of permutation columns and one fisher_ prefixed set.
build_comparison_table <- function(association, fisher_table) {
  cells <- fisher_table[fisher_table$universe == "All tested genes", , drop = FALSE]

  fisher_slim <- data.frame(
    test_id                  = test_key(cells$alist_set, cells$blist_set),
    fisher_test_id           = cells$test_id,
    fisher_odds_ratio        = cells$odds_ratio,
    fisher_ci_lower          = cells$ci_lower,
    fisher_ci_upper          = cells$ci_upper,
    fisher_p                 = cells$p_value,
    fisher_q                 = cells$q_value,
    fisher_n_genes           = cells$n_genes,
    fisher_n_direction_genes = cells$n_row_true,
    fisher_n_overlap_genes   = cells$n_col_true,
    fisher_n_both            = cells$n_both_true,
    stringsAsFactors = FALSE
  )

  merged <- dplyr::left_join(association, fisher_slim, by = "test_id")

  if (any(is.na(merged$fisher_odds_ratio))) {
    unmatched <- merged$test_id[is.na(merged$fisher_odds_ratio)]
    stop("No matching gene-level Fisher test for permutation cell(s): ",
         paste(unmatched, collapse = ", "))
  }

  merged$permutation_p <- merged$p_value
  merged$fisher_z_equiv <- fisher_p_to_signed_z(merged$fisher_p,
                                                merged$fisher_odds_ratio)
  merged$concordance <- classify_concordance(merged$fisher_p, merged$p_value)
  merged
}

# =============================================================================
# FIGURES
# =============================================================================

plot_crosswise_heatmap <- function(association, ntimes, out_dir) {
  df <- association
  df$alist_set <- factor(df$alist_set, levels = rev(ALIST_ORDER))
  df$blist_set <- factor(df$blist_set, levels = BLIST_ORDER)
  df$tile_label <- sprintf("z = %.1f %s\nobs %s\nexp %s",
                           df$z_score, df$sig_label,
                           fmt_count(df$observed_overlaps),
                           fmt_count(df$null_mean))

  limit <- max(abs(df$z_score))

  p <- ggplot(df, aes(x = blist_set, y = alist_set, fill = z_score)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = tile_label), size = 2.9, lineheight = 0.95) +
    scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                         midpoint = 0, limits = c(-limit, limit),
                         name = "Permutation\nz-score") +
    labs(
      title = "mCH gene bodies against chromatin mark peaks",
      subtitle = sprintf(
        "crosswisePermTest, %s randomisations, %s per chromosome, evaluation %s. Stars use the permutation p-value.",
        format(ntimes, big.mark = ","), PERM_RANFUN, PERM_EVFUN),
      x = NULL, y = NULL
    ) +
    theme_emseq() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))

  save_multiformat_ggplot(p, file.path(out_dir, "40_01a_crosswise_zscore_heatmap"),
                          width = 13, height = 5.5)
}

plot_zscore_forest <- function(association, ntimes, out_dir) {
  df <- association
  df$blist_set <- factor(df$blist_set, levels = rev(BLIST_ORDER))
  df$alist_set <- factor(df$alist_set, levels = ALIST_ORDER)
  df$point_label <- sprintf("z = %.1f %s  (%s)", df$z_score, df$sig_label,
                            vapply(df$p_value, fmt_p, character(1)))

  # Span of the plotted z axis, floored at 1 so the labels stay off the points
  # even when every interval sits at nearly the same z.
  span <- max(diff(range(c(df$z_ci_lower, df$z_ci_upper))), 1)
  df$label_x <- max(df$z_ci_upper) + 0.08 * span
  axis_max <- df$label_x[1] + 0.55 * span

  p <- ggplot(df, aes(x = z_score, y = blist_set, colour = blist_set)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = c(-1.96, 1.96), linetype = "dotted", colour = "grey70") +
    geom_errorbarh(aes(xmin = z_ci_lower, xmax = z_ci_upper), height = 0.22,
                   linewidth = 0.7) +
    geom_point(size = 3.2) +
    geom_text(aes(x = label_x, label = point_label), hjust = 0, size = 2.9,
              colour = "black") +
    expand_limits(x = axis_max) +
    facet_wrap(~ alist_set, ncol = 1) +
    scale_colour_manual(values = BLIST_COLORS, guide = "none") +
    labs(
      title = "Permutation z-scores for mCH gene bodies against chromatin marks",
      subtitle = sprintf(
        "Bars are 95%% Monte Carlo intervals from %s randomisations. Dotted lines mark z = +/- 1.96.",
        format(ntimes, big.mark = ",")),
      x = "Permutation z-score", y = NULL
    ) +
    theme_emseq()

  save_multiformat_ggplot(p, file.path(out_dir, "40_01b_zscore_forest"),
                          width = 12, height = 9)
}

plot_null_histograms <- function(draws, null_summary, out_dir) {
  df <- draws
  df$facet_label <- sprintf("%s x %s",
                            unname(ALIST_SHORT[df$alist_set]), df$blist_set)

  annot <- null_summary
  annot$facet_label <- sprintf("%s x %s",
                               unname(ALIST_SHORT[annot$alist_set]), annot$blist_set)
  annot$text_label <- sprintf(
    "n perm = %s\nobs = %s\nnull med = %s\nnull sd = %.1f\nz = %.1f",
    fmt_count(annot$n_permutations),
    fmt_count(annot$observed_overlaps),
    fmt_count(annot$null_median),
    annot$null_sd_draws,
    annot$null_run_z)

  facet_levels <- unlist(lapply(ALIST_ORDER, function(a) {
    sprintf("%s x %s", unname(ALIST_SHORT[a]), BLIST_ORDER)
  }))
  df$facet_label <- factor(df$facet_label, levels = facet_levels)
  annot$facet_label <- factor(annot$facet_label, levels = facet_levels)

  df$alist_set <- factor(df$alist_set, levels = ALIST_ORDER)

  p <- ggplot(df, aes(x = permuted_overlaps, fill = alist_set)) +
    geom_histogram(bins = 40, colour = "grey35", linewidth = 0.15,
                   alpha = 0.75) +
    scale_fill_manual(values = ALIST_COLORS, name = NULL) +
    geom_vline(data = annot, aes(xintercept = null_median),
               linetype = "dashed", colour = "grey30", linewidth = 0.5) +
    geom_vline(data = annot, aes(xintercept = observed_overlaps),
               colour = "black", linewidth = 0.9) +
    geom_text(data = annot, aes(x = -Inf, y = Inf, label = text_label),
              inherit.aes = FALSE, hjust = -0.06, vjust = 1.12, size = 2.5,
              lineheight = 0.95) +
    facet_wrap(~ facet_label, ncol = 4, scales = "free") +
    scale_x_continuous(labels = scales::label_number(big.mark = ",")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(
      title = "Null overlap distributions with the observed count marked",
      subtitle = paste0("Bars are permuted overlap counts from randomised gene bodies. ",
                        "Black line is the observed count; dashed line is the null median."),
      x = "Overlapping regions per randomisation", y = "Randomisations"
    ) +
    theme_emseq() +
    theme(strip.text = element_text(size = 8.5),
          legend.position = "top")

  save_multiformat_ggplot(p, file.path(out_dir, "40_01c_null_distributions"),
                          width = 17, height = 9)
}

plot_fisher_vs_permutation <- function(comparison, ntimes, out_dir) {
  long <- comparison %>%
    dplyr::select(test_id, blist_set, alist_set, mark_group, concordance,
                  fisher_z_equiv, z_score) %>%
    tidyr::pivot_longer(cols = c("fisher_z_equiv", "z_score"),
                        names_to = "method", values_to = "z") %>%
    dplyr::mutate(
      method = ifelse(method == "fisher_z_equiv", "Fisher", "Permutation"),
      mark_group = factor(mark_group, levels = BLIST_GROUP_ORDER))

  p <- ggplot(long, aes(x = z, y = reorder(test_id, z), shape = method,
                        colour = concordance)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = c(-1.96, 1.96), linetype = "dotted", colour = "grey70") +
    geom_point(size = 3, alpha = 0.9) +
    scale_shape_manual(values = c("Fisher" = 16, "Permutation" = 17),
                       name = "Method") +
    scale_colour_manual(values = CONCORDANCE_COLORS, name = "Concordance",
                        drop = FALSE) +
    facet_grid(mark_group ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = "Gene-level Fisher tests against the genomic permutation",
      subtitle = sprintf(
        "Fisher z is qnorm of the two-sided p-value signed by the odds ratio; permutation z is from %s randomisations.",
        format(ntimes, big.mark = ",")),
      x = "z-score", y = NULL
    ) +
    theme_emseq() +
    theme(strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
          axis.text.y = element_text(size = 8))

  save_multiformat_ggplot(p, file.path(out_dir, "40_01d_fisher_vs_permutation"),
                          width = 14, height = 10)
}

plot_observed_vs_expected <- function(association, out_dir) {
  long <- association %>%
    dplyr::select(alist_set, blist_set, observed_overlaps, null_mean, null_sd,
                  null_q025, null_q975, observed_over_expected) %>%
    tidyr::pivot_longer(cols = c("observed_overlaps", "null_mean"),
                        names_to = "source", values_to = "count") %>%
    dplyr::mutate(
      source = ifelse(source == "observed_overlaps", "Observed", "Null mean"),
      source = factor(source, levels = c("Observed", "Null mean")),
      blist_set = factor(blist_set, levels = BLIST_ORDER),
      alist_set = factor(alist_set, levels = ALIST_ORDER),
      error_lower = ifelse(source == "Null mean", null_q025, NA_real_),
      error_upper = ifelse(source == "Null mean", null_q975, NA_real_))

  ratio_labels <- association
  ratio_labels$blist_set <- factor(ratio_labels$blist_set, levels = BLIST_ORDER)
  ratio_labels$alist_set <- factor(ratio_labels$alist_set, levels = ALIST_ORDER)
  ratio_labels$label <- sprintf("O/E = %.2f", ratio_labels$observed_over_expected)
  ratio_labels$y_pos <- pmax(ratio_labels$observed_overlaps, ratio_labels$null_q975)

  p <- ggplot(long, aes(x = blist_set, y = count, fill = source)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72,
             colour = "black", linewidth = 0.2) +
    geom_errorbar(aes(ymin = error_lower, ymax = error_upper),
                  position = position_dodge(width = 0.8), width = 0.2,
                  linewidth = 0.4, na.rm = TRUE) +
    geom_text(aes(label = fmt_count(count)),
              position = position_dodge(width = 0.8), vjust = -0.35, size = 2.5) +
    geom_text(data = ratio_labels,
              aes(x = blist_set, y = y_pos, label = label),
              inherit.aes = FALSE, vjust = -2.2, size = 2.9, fontface = "bold") +
    facet_wrap(~ alist_set, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c("Observed" = "#D7191C", "Null mean" = "grey70"),
                      name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(
      title = "Observed overlaps against the permutation null",
      subtitle = "Error bars span the 2.5th to 97.5th percentile of the null draws",
      x = NULL, y = "Gene bodies overlapping the peak set"
    ) +
    theme_emseq() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 25, hjust = 1))

  save_multiformat_ggplot(p, file.path(out_dir, "40_01e_observed_vs_expected"),
                          width = 12, height = 9)
}

# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

print_association_summary <- function(association) {
  cat("\n--- Permutation results ---\n")
  cat(sprintf("  %-24s %-16s %10s %10s %8s %14s %8s\n",
              "mCH set", "Peak set", "observed", "null mean", "z", "p", "O/E"))
  for (i in seq_len(nrow(association))) {
    row <- association[i, ]
    cat(sprintf("  %-24s %-16s %10s %10.1f %8.2f %14.3g %8.2f\n",
                row$alist_set, row$blist_set,
                format(row$observed_overlaps, big.mark = ","),
                row$null_mean, row$z_score, row$p_value,
                row$observed_over_expected))
  }
}

print_concordance_summary <- function(comparison) {
  cat("\n--- Fisher versus permutation concordance ---\n")
  for (level in CONCORDANCE_ORDER) {
    n <- sum(comparison$concordance == level)
    cat(sprintf("  %-16s %2d / %2d (%.0f%%)\n", paste0(level, ":"), n,
                nrow(comparison), 100 * n / nrow(comparison)))
  }
}

print_replicate_summary <- function(association) {
  cat("\n--- Agreement between the two permutation runs ---\n")
  cat(sprintf("  Largest z-score difference: %.3f\n",
              max(abs(association$z_difference))))
  cat(sprintf("  Median absolute difference: %.3f\n",
              median(abs(association$z_difference))))
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opt <- parse_section_options()
  out_dir <- opt$output_dir
  ntimes <- opt$ntimes
  cores <- resolve_cores()
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 40_01: mCH GENE BODY x CHROMATIN MARK GENOMIC PERMUTATION\n")
  cat("================================================================================\n")
  cat(sprintf("Output directory:   %s\n", out_dir))
  cat(sprintf("Randomisations:     %s\n", format(ntimes, big.mark = ",")))
  cat(sprintf("Cores:              %d\n", cores))
  cat(sprintf("Randomisation:      %s (per chromosome = %s)\n",
              PERM_RANFUN, PERM_PER_CHR))
  cat(sprintf("Evaluation:         %s\n", PERM_EVFUN))
  cat(sprintf("Seed:               %d\n", PERM_SEED))
  cat(sprintf("regioneR:           %s\n", packageVersion("regioneR")))
  cat(sprintf("regioneReloaded:    %s\n", packageVersion("regioneReloaded")))
  cat("\n")

  # --- Step 1: genome and region sets ---------------------------------------
  cat("STEP 1: Genome and region sets\n")
  genome <- build_permutation_genome()
  alist <- build_alist(genome)
  blist <- build_blist(genome)

  region_sets <- build_region_set_table(alist, blist)
  write_section_table(region_sets, file.path(out_dir, "40_01_region_set_sizes.tsv"))

  # --- Step 2: permutation ---------------------------------------------------
  cat("\nSTEP 2: Permutation\n")
  cache_path <- file.path(out_dir, CACHE_FILENAME)
  perm <- load_or_run_permutations(alist, blist, genome, ntimes, cores, cache_path)

  if (perm$ntimes != ntimes) {
    cat(sprintf("  The cache holds %s randomisations; --ntimes asked for %s. ",
                format(perm$ntimes, big.mark = ","),
                format(ntimes, big.mark = ",")))
    cat("Reporting the cached run.\n")
    ntimes <- perm$ntimes
  }

  parameters <- build_parameter_table(genome, alist, blist, ntimes, cores)
  write_section_table(parameters, file.path(out_dir, "40_01_permutation_parameters.tsv"))

  # --- Step 3: result extraction ---------------------------------------------
  cat("\nSTEP 3: Extracting permutation results\n")
  crosswise <- extract_crosswise_table(perm$crosswise, ntimes)
  null_draws <- extract_null_draws(perm$nulls)
  null_summary <- summarise_null_draws(perm$nulls, null_draws)
  check_observed_counts_agree(crosswise, null_summary)

  association <- build_association_table(crosswise, null_summary)
  write_section_table(association, file.path(out_dir, "40_01_permutation_association.tsv"))
  write_section_table(null_summary, file.path(out_dir, "40_01_null_distribution_summary.tsv"))
  write_section_table(null_draws, file.path(out_dir, "40_01_null_draws.tsv"))
  print_association_summary(association)
  print_replicate_summary(association)

  # --- Step 4: gene-level Fisher tests ---------------------------------------
  cat("\nSTEP 4: Gene-level Fisher tests\n")
  genes <- build_gene_overlap_table(mch_results, blist)
  fisher_table <- run_gene_fisher_tests(genes, blist, out_dir)
  write_section_table(fisher_table, file.path(out_dir, "40_01_fisher_gene_level.tsv"))

  gene_columns <- c("gene_name", "gene_id", "chr", "start", "end", "gene_length",
                    "mch_ctrl", "mch_mut", "mch_diff", "edger_logFC", "edger_fdr",
                    "mch_sig", "mch_hyper", "mch_hypo",
                    vapply(names(blist), overlap_column, character(1)))
  write_section_table(genes[, gene_columns, drop = FALSE],
                  file.path(out_dir, "40_01_gene_peak_overlaps.tsv"))

  # --- Step 5: Fisher versus permutation -------------------------------------
  cat("\nSTEP 5: Fisher versus permutation\n")
  comparison <- build_comparison_table(association, fisher_table)
  write_section_table(comparison, file.path(out_dir, "40_01_fisher_vs_permutation.tsv"))
  print_concordance_summary(comparison)

  # --- Step 6: figures -------------------------------------------------------
  cat("\nSTEP 6: Figures\n")
  plot_crosswise_heatmap(association, ntimes, out_dir)
  plot_zscore_forest(association, ntimes, out_dir)
  plot_null_histograms(null_draws, null_summary, out_dir)
  plot_fisher_vs_permutation(comparison, ntimes, out_dir)
  plot_observed_vs_expected(association, out_dir)

  # --- Step 7: summary -------------------------------------------------------
  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 40_01 SUMMARY\n")
  cat("================================================================================\n")
  cat(sprintf("Pairwise permutation tests: %d\n", nrow(association)))
  cat(sprintf("Significant at p < %.2f:     %d\n", Q_THRESHOLD,
              sum(association$p_value < Q_THRESHOLD)))

  strongest <- association[which.max(abs(association$z_score)), ]
  cat(sprintf("Strongest association:      %s (z = %.2f, %s, O/E = %.2f)\n",
              strongest$test_id, strongest$z_score, fmt_p(strongest$p_value),
              strongest$observed_over_expected))

  cat(sprintf("Gene-level Fisher tests registered: %d\n", nrow(fisher_table)))
  cat(sprintf("Registry: %s\n", HANDOFF_PATHS$fisher_registry))
  cat(sprintf("Cache:    %s\n", cache_path))
  cat(sprintf("Figures and tables written to: %s\n", out_dir))
  cat("Section 40_01 complete.\n")
}

main()
