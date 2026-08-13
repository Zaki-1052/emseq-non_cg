# scripts/sections/40_permutation/40_02_atac_loops.R
#
# Section 40_02: genomic permutation validation of chromatin accessibility and
# Hi-C loop anchor co-occurrence with chromatin features.
#
# The question: do differential ATAC peaks and differential Hi-C loop anchors
# sit on chromatin features more often than intervals of the same size, count,
# and chromosome placed at random? A Fisher's exact test on interval counts
# treats every interval as an independent draw. Genomic intervals are not
# independent, so this section repeats each contrast as a permutation test that
# randomises interval positions while it preserves size, count, and chromosome.
#
# Four sub-analyses, each a crosswise (many A by many B) permutation matrix:
#   40_02a  ATAC up and ATAC down peaks against six histone and CTCF marks
#   40_02b  ATAC up and ATAC down peaks against gene bodies grouped by the five
#           body chromatin states
#   40_02m  ATAC up and ATAC down peaks against gene bodies grouped by the four
#           promoter chromatin states
#   40_02c  Gained and lost Hi-C loop anchors against the six marks, ATAC peaks,
#           MeCP2 peaks, and mCH hypermethylated and hypomethylated gene bodies
#
# Section 10_01 records two chromatin states per gene: promoter_state over the
# TSS window and body_state over the gene body. Blocks b and m read one state
# column each. Both group the same gene-body intervals; only the label differs.
#
# Every permutation association is paired with the matching interval-level
# Fisher's exact test, so the two methods can be read side by side.
#
# Reads:
#   CHIP_PATHS$ctcf .. $bivalent        six mark BED files
#   DIFFBIND_PATHS$atac_up / $atac_down differential ATAC peak BED files
#   MECP2_PATHS$up / $down              differential MeCP2 peak BED files
#   HIC_PATHS$loops                     differential Hi-C loop table
#   HANDOFF_PATHS$chromatin_state       gene chromatin state table from 10_01
#   gene_bodies                         pre-loaded by _shared_config.R
#   BSgenome.Mmusculus.UCSC.mm10        mm10 chromosome lengths through regioneR
#
# Writes: TSV tables, RDS permutation caches, and multi-format figures under
#   OUTPUT_PATHS$permutation (results/sections/40_permutation/).
#
# Adapted from the Biomodal script section_35_permutation_atac_loops.R. That
# script tested mC hyper and hypo DMRs in sub-analysis C; this single-modality
# port tests mCH hypermethylated and hypomethylated gene bodies instead, and
# adds the six marks to the sub-analysis C target list.

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)
library(regioneR)
library(regioneReloaded)
library(BSgenome.Mmusculus.UCSC.mm10)

# =============================================================================
# SECTION CONSTANTS
# =============================================================================

SECTION_ID <- "40_02"

# Permutation runs on the chromosomes that the mm10 assembly and every input
# peak set share. Unplaced contigs, chrM, and chrY have too few peaks for a
# per-chromosome randomisation.
ANALYSIS_CHRS <- paste0("chr", c(1:19, "X"))

DEFAULT_NTIMES <- 5000L
DEFAULT_CORES  <- 8L
DEFAULT_SEED   <- 42L

PERM_RANFUN         <- "randomizeRegions"
PERM_EVFUN          <- "numOverlaps"
PERM_PER_CHROMOSOME <- TRUE

# Local z-score geometry for the strongest loop anchor association.
LZ_WINDOW <- 50000L
LZ_STEP   <- 1000L

# Names on the left are CHIP_PATHS keys; names on the right label the figures.
MARK_DISPLAY <- c(
  ctcf     = "CTCF",
  h3k27ac  = "H3K27ac",
  h3k27me3 = "H3K27me3",
  h3k4me1  = "H3K4me1",
  h3k4me3  = "H3K4me3",
  bivalent = "Bivalent"
)

# Loop direction coding used by the upstream Hi-C pipeline.
LOOP_DIRECTION_LABELS <- c(up_in_mutant = "Gained loop anchors",
                           down_in_mutant = "Lost loop anchors")

SUB_ANALYSIS_LABELS <- c(
  a = "40_02a  ATAC x marks",
  b = "40_02b  ATAC x body states",
  m = "40_02m  ATAC x promoter states",
  c = "40_02c  Loop anchors x features"
)

# Columns section 10_01 writes into gene_chromatin_state.tsv.
CHROMATIN_STATE_COLUMNS <- c(
  "gene_name", "gene_id", "chr", "start", "end",
  "promoter_state", "body_state",
  "prom_ctcf_overlap", "prom_h3k27ac_overlap", "prom_h3k27me3_overlap",
  "prom_h3k4me1_overlap", "prom_h3k4me3_overlap", "prom_bivalent_overlap",
  "body_ctcf_overlap", "body_h3k27ac_overlap", "body_h3k27me3_overlap",
  "body_h3k4me1_overlap", "body_h3k4me3_overlap", "body_bivalent_overlap")

CONCORDANCE_LEVELS <- c("Confirmed", "Weakened", "Strengthened", "Concordant NS")

CONCORDANCE_COLORS <- c(
  "Confirmed"     = "#1B7837",
  "Weakened"      = "#D7191C",
  "Strengthened"  = "#E08214",
  "Concordant NS" = "#808080"
)

OBSERVED_EXPECTED_COLORS <- c(
  "Observed"                 = "#B2182B",
  "Expected (permuted mean)" = "#4D4D4D"
)

# Two-sided normal cutoff drawn on the forest plots.
Z_CUTOFF <- 1.96

# =============================================================================
# SMALL UTILITIES
# =============================================================================

#' Format an integer with thousands separators.
fmt_int <- function(x) format(x, big.mark = ",", trim = TRUE)

#' Format a p-value for a figure subtitle.
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

#' Significance stars for a vector of p-values.
significance_stars <- function(p) {
  ifelse(is.na(p), "NA",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", "ns"))))
}

#' Stop unless a data.frame carries every named column.
require_columns <- function(df, required, what) {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(what, " is missing columns: ", paste(missing, collapse = ", "),
         "\n  Found: ", paste(colnames(df), collapse = ", "))
  }
  invisible(TRUE)
}

# =============================================================================
# OPTIONS AND INPUT CHECKS
# =============================================================================

#' Default core count, taken from the SLURM allocation when the job sets one.
default_cores <- function() {
  slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
  if (!nzchar(slurm_cpus)) return(DEFAULT_CORES)
  cores <- suppressWarnings(as.integer(slurm_cpus))
  if (is.na(cores) || cores < 1L) {
    stop("SLURM_CPUS_PER_TASK is set to '", slurm_cpus,
         "', which is not a positive integer.")
  }
  cores
}

#' Parse command line options.
#'
#' The permutation count can arrive either as --ntimes or as a single bare
#' positional argument, which is how the SLURM wrapper calls the permutation
#' sections (Rscript 40_02_atac_loops.R 5000).
parse_options <- function() {
  option_list <- list(
    make_option("--output-dir", dest = "output_dir", type = "character",
                default = OUTPUT_PATHS$permutation,
                help = "Directory for figures, tables, and caches [default: %default]"),
    make_option("--ntimes", dest = "ntimes", type = "integer",
                default = DEFAULT_NTIMES,
                help = "Permutations per pairwise test [default: %default]"),
    make_option("--cores", dest = "cores", type = "integer",
                default = default_cores(),
                help = "Cores for the permutation loop [default: %default]"),
    make_option("--seed", dest = "seed", type = "integer", default = DEFAULT_SEED,
                help = "Random seed set before each permutation run [default: %default]")
  )

  parsed <- parse_args(
    OptionParser(
      option_list = option_list,
      usage = "usage: %prog [options] [ntimes]",
      description = paste("Permutation validation of ATAC peak and Hi-C loop",
                          "anchor co-occurrence with chromatin features.")
    ),
    positional_arguments = c(0, 1)
  )

  opt <- parsed$options
  if (length(parsed$args) == 1) {
    ntimes <- suppressWarnings(as.integer(parsed$args[1]))
    if (is.na(ntimes) || ntimes < 1L) {
      stop("The positional ntimes argument must be a positive integer, got: ",
           parsed$args[1])
    }
    opt$ntimes <- ntimes
  }

  if (opt$ntimes < 1L) stop("--ntimes must be at least 1, got ", opt$ntimes)
  if (opt$cores < 1L) stop("--cores must be at least 1, got ", opt$cores)
  opt
}

#' Stop unless every section input file exists.
check_inputs <- function() {
  mark_paths <- unlist(CHIP_PATHS[names(MARK_DISPLAY)])
  names(mark_paths) <- paste0("mark_", names(MARK_DISPLAY))

  required <- c(
    mark_paths,
    c(atac_up         = DIFFBIND_PATHS$atac_up,
      atac_down       = DIFFBIND_PATHS$atac_down,
      mecp2_up        = MECP2_PATHS$up,
      mecp2_down      = MECP2_PATHS$down,
      loops           = HIC_PATHS$loops,
      chromatin_state = HANDOFF_PATHS$chromatin_state)
  )

  missing <- names(required)[!file.exists(required)]
  if (length(missing) == 0) return(invisible(TRUE))

  detail <- paste(sprintf("    %s: %s", missing, required[missing]),
                  collapse = "\n")
  extra <- ""
  if ("chromatin_state" %in% missing) {
    extra <- paste0("\n  The chromatin state table is the handoff from section",
                    " 10_01. Run 10_01 first.")
  }
  stop("Section ", SECTION_ID, " inputs not found:\n", detail, extra)
}

# =============================================================================
# GENOME AND REGION SET PREPARATION
# =============================================================================

#' Build the mm10 genome GRanges restricted to the analysis chromosomes.
build_analysis_genome <- function() {
  genome_full <- getGenomeAndMask("mm10")$genome
  genome <- genome_full[as.character(seqnames(genome_full)) %in% ANALYSIS_CHRS]

  found <- unique(as.character(seqnames(genome)))
  missing <- setdiff(ANALYSIS_CHRS, found)
  if (length(missing) > 0) {
    stop("mm10 genome object has no entry for: ", paste(missing, collapse = ", "))
  }

  cat(sprintf("  mm10 genome: %d chromosomes, %s bp\n",
              length(genome), fmt_int(sum(width(genome)))))
  genome
}

#' Drop regions outside the analysis chromosomes and report the loss.
#'
#' randomizeRegions with per.chromosome = TRUE needs every region to sit on a
#' chromosome that the genome object defines.
restrict_to_analysis_chrs <- function(gr, set_name) {
  n_before <- length(gr)
  gr <- gr[as.character(seqnames(gr)) %in% ANALYSIS_CHRS]
  if (length(gr) == 0) {
    stop("Region set '", set_name, "' has no regions on ",
         paste(ANALYSIS_CHRS, collapse = ", "))
  }
  seqlevels(gr) <- seqlevelsInUse(gr)

  cat(sprintf("  %-28s %8s regions kept, %6s dropped off chr1-19 and chrX\n",
              set_name, fmt_int(length(gr)), fmt_int(n_before - length(gr))))
  gr
}

#' Load one BED peak set and restrict it to the analysis chromosomes.
load_peak_set <- function(path, set_name) {
  gr <- load_chip_peaks(path, set_name)
  restrict_to_analysis_chrs(gr, set_name)
}

#' The six histone and CTCF mark peak sets, named for the figures.
load_mark_sets <- function() {
  cat("\nLoading mark peak sets...\n")
  sets <- lapply(names(MARK_DISPLAY), function(key) {
    load_peak_set(CHIP_PATHS[[key]], MARK_DISPLAY[[key]])
  })
  names(sets) <- unname(MARK_DISPLAY)
  sets
}

#' Differential ATAC peak sets.
load_atac_sets <- function() {
  cat("\nLoading differential ATAC peak sets...\n")
  list(
    "ATAC Up"   = load_peak_set(DIFFBIND_PATHS$atac_up, "ATAC Up"),
    "ATAC Down" = load_peak_set(DIFFBIND_PATHS$atac_down, "ATAC Down")
  )
}

#' Differential MeCP2 peak sets.
load_mecp2_sets <- function() {
  cat("\nLoading differential MeCP2 peak sets...\n")
  list(
    "MeCP2 Up"   = load_peak_set(MECP2_PATHS$up, "MeCP2 Up"),
    "MeCP2 Down" = load_peak_set(MECP2_PATHS$down, "MeCP2 Down")
  )
}

#' Gene bodies that gained and lost mCH in the mutant.
build_mch_sets <- function() {
  cat("\nBuilding mCH gene body sets...\n")
  # which() keeps a gene with a missing mch_sig or mch_diff out of both sets
  # instead of turning it into an NA-valued range.
  hyper <- gene_bodies[which(gene_bodies$mch_sig & gene_bodies$edger_logFC > 0)]
  hypo  <- gene_bodies[which(gene_bodies$mch_sig & gene_bodies$edger_logFC < 0)]

  list(
    "mCH Hyper genes" = restrict_to_analysis_chrs(hyper, "mCH Hyper genes"),
    "mCH Hypo genes"  = restrict_to_analysis_chrs(hypo, "mCH Hypo genes")
  )
}

#' Read the differential loop table and build gained and lost anchor sets.
#'
#' Each significant loop contributes both of its anchors. Overlapping anchors
#' of the same direction are merged with reduce(), so an anchor region counts
#' once however many loops share it.
build_loop_anchor_sets <- function(path) {
  cat("\nLoading Hi-C loop anchors...\n")
  loops <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                      quote = "", comment.char = "")

  col_map <- c(chr1 = "anchor1_chr", start1 = "anchor1_start", end1 = "anchor1_end",
               chr2 = "anchor2_chr", start2 = "anchor2_start", end2 = "anchor2_end")
  hit <- names(col_map) %in% colnames(loops)
  if (any(hit)) colnames(loops)[match(names(col_map)[hit], colnames(loops))] <- col_map[hit]

  require_columns(
    loops,
    c("loop_id", "direction", "significant",
      "anchor1_chr", "anchor1_start", "anchor1_end",
      "anchor2_chr", "anchor2_start", "anchor2_end"),
    paste0("Loop table ", path)
  )
  cat(sprintf("  Loop table: %s rows\n", fmt_int(nrow(loops))))

  loops$significant <- as.logical(loops$significant)
  loops <- loops[!is.na(loops$significant) & loops$significant, , drop = FALSE]
  if (nrow(loops) == 0) stop("No significant loops in ", path)

  unknown <- setdiff(unique(loops$direction), names(LOOP_DIRECTION_LABELS))
  if (length(unknown) > 0) {
    stop("Unexpected loop direction values in ", path, ": ",
         paste(unknown, collapse = ", "))
  }

  # Anchor coordinates are used exactly as the loop table writes them, which is
  # the same convention section 30_01 uses when it builds its anchor GRanges.
  anchors_of <- function(rows) {
    a1 <- GRanges(seqnames = rows$anchor1_chr,
                  ranges = IRanges(start = rows$anchor1_start,
                                   end = rows$anchor1_end))
    a2 <- GRanges(seqnames = rows$anchor2_chr,
                  ranges = IRanges(start = rows$anchor2_start,
                                   end = rows$anchor2_end))
    reduce(c(a1, a2))
  }

  sets <- lapply(names(LOOP_DIRECTION_LABELS), function(direction) {
    rows <- loops[loops$direction == direction, , drop = FALSE]
    label <- LOOP_DIRECTION_LABELS[[direction]]
    if (nrow(rows) == 0) {
      stop("No significant loops with direction '", direction, "' in ", path)
    }
    gr <- anchors_of(rows)
    cat(sprintf("  %-28s %8s merged anchors from %s loops\n",
                label, fmt_int(length(gr)), fmt_int(nrow(rows))))
    restrict_to_analysis_chrs(gr, label)
  })
  names(sets) <- unname(LOOP_DIRECTION_LABELS)
  sets
}

#' Read the section 10_01 gene chromatin state table.
load_chromatin_state_table <- function(path) {
  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")
  require_columns(df, CHROMATIN_STATE_COLUMNS,
                  paste0("Chromatin state table ", path))
  cat(sprintf("\nChromatin state table: %s genes\n", fmt_int(nrow(df))))
  df
}

#' Split the gene bodies of one state column into one region set per state.
#'
#' Adds one to the start column so these intervals match the gene_bodies
#' object, which _shared_config.R builds with mch_results$start + 1L.
#'
#' @param df table from load_chromatin_state_table().
#' @param state_column "body_state" or "promoter_state".
#' @param state_order the state levels, in the order the figures use.
#' @param label display name used in the messages and the region set summary.
#' @return named list of GRanges, one per state
build_state_region_sets <- function(df, state_column, state_order, label) {
  cat(sprintf("\nBuilding %s region sets from %s...\n", label, state_column))

  unknown <- setdiff(unique(df[[state_column]]), state_order)
  if (length(unknown) > 0) {
    stop("Unexpected ", state_column, " values: ",
         paste(unknown, collapse = ", "),
         "\n  Expected only: ", paste(state_order, collapse = ", "))
  }

  gr <- GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start + 1L, end = df$end),
    gene_name = df$gene_name,
    state = df[[state_column]]
  )
  gr <- restrict_to_analysis_chrs(gr, paste(label, "genes"))

  sets <- lapply(state_order, function(state) {
    state_gr <- gr[gr$state == state]
    cat(sprintf("  %-28s %8s gene bodies\n", state, fmt_int(length(state_gr))))
    state_gr
  })
  names(sets) <- state_order

  empty <- names(sets)[vapply(sets, length, integer(1)) == 0]
  if (length(empty) > 0) {
    stop(label, " with no gene body on the analysis chromosomes: ",
         paste(empty, collapse = ", "),
         "\n  A permutation test needs a non-empty region set for every state.")
  }
  sets
}

#' One summary row per region set: count, span, and width distribution.
summarise_region_sets <- function(sets, role, sub_analysis) {
  rows <- lapply(names(sets), function(set_name) {
    gr <- sets[[set_name]]
    data.frame(
      sub_analysis  = sub_analysis,
      role          = role,
      region_set    = set_name,
      n_regions     = length(gr),
      n_chromosomes = length(seqlevelsInUse(gr)),
      total_bp      = sum(as.numeric(width(gr))),
      median_width  = median(width(gr)),
      mean_width    = mean(width(gr)),
      min_width     = min(width(gr)),
      max_width     = max(width(gr)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# =============================================================================
# PERMUTATION TESTS
# =============================================================================

#' Run one crosswise permutation matrix, reusing an RDS cache when present.
#'
#' FORCE_RERUN in the environment makes the run ignore the cache.
run_crosswise <- function(Alist, Blist, genome, ntimes, cores, seed,
                          cache_path, label, force_rerun) {
  if (file.exists(cache_path) && !force_rerun) {
    cat(sprintf("[%s] Loading cached permutation result: %s\n", label, cache_path))
    return(readRDS(cache_path))
  }

  cat(sprintf("[%s] %d x %d = %d pairwise tests, %s permutations each\n",
              label, length(Alist), length(Blist),
              length(Alist) * length(Blist), fmt_int(ntimes)))

  options(mc.cores = cores)
  set.seed(seed)
  cw <- crosswisePermTest(
    Alist          = Alist,
    Blist          = Blist,
    genome         = genome,
    ranFUN         = PERM_RANFUN,
    evFUN          = PERM_EVFUN,
    ntimes         = ntimes,
    force.parallel = TRUE,
    per.chromosome = PERM_PER_CHROMOSOME
  )

  # symm_matrix = FALSE because the design is not square. The patched
  # chooseHclustMet handles the two-row clustering.
  cw <- makeCrosswiseMatrix(cw, pvcut = 1, symm_matrix = FALSE,
                            hc.method = "average")

  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cw, cache_path)
  cat(sprintf("[%s] Saved cache: %s\n", label, cache_path))
  cw
}

#' Flatten a crosswise result into one row per A-by-B pair.
#'
#' regioneReloaded names the adjusted p-value column "adj.p_value". The column
#' check below fails loudly if a future version renames it.
extract_permutation_table <- function(cw, sub_analysis) {
  evaluation <- getMultiEvaluation(cw)
  if (length(evaluation) == 0) {
    stop("Crosswise result for ", sub_analysis, " holds no evaluation tables. ",
         "One or more permutation tests failed inside crosswisePermTest().")
  }

  required <- c("name", "n_regionA", "n_regionB", "n_overlaps", "mean_perm_test",
                "sd_perm_test", "z_score", "norm_zscore", "p_value", "adj.p_value")

  rows <- lapply(names(evaluation), function(rs1) {
    df <- evaluation[[rs1]]
    require_columns(df, required,
                    paste0("regioneReloaded evaluation table for ", rs1))
    data.frame(
      sub_analysis      = sub_analysis,
      RS1               = rs1,
      RS2               = as.character(df$name),
      n_region_a        = df$n_regionA,
      n_region_b        = df$n_regionB,
      observed_overlaps = df$n_overlaps,
      expected_overlaps = df$mean_perm_test,
      sd_permuted       = df$sd_perm_test,
      z_score           = df$z_score,
      norm_zscore       = df$norm_zscore,
      perm_p_value      = df$p_value,
      perm_adj_p_value  = df[["adj.p_value"]],
      stringsAsFactors  = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out$observed_over_expected <- out$observed_overlaps / out$expected_overlaps
  rownames(out) <- NULL
  out
}

# =============================================================================
# INTERVAL-LEVEL FISHER TESTS
# =============================================================================

#' Fisher's exact test contrasting two query region sets on one target set.
#'
#' Each interval of a query set is one observation, scored for carrying at
#' least one target region. This is an interval-level test, so it calls
#' fisher.test() directly rather than register_fisher_test().
contrast_fisher <- function(query_gr, contrast_gr, target_gr,
                            query_name, contrast_name, target_name,
                            sub_analysis) {
  query_hit    <- sum(countOverlaps(query_gr, target_gr) > 0)
  contrast_hit <- sum(countOverlaps(contrast_gr, target_gr) > 0)

  tab <- matrix(
    c(query_hit, length(query_gr) - query_hit,
      contrast_hit, length(contrast_gr) - contrast_hit),
    nrow = 2,
    dimnames = list(c("with_target", "without_target"),
                    c(query_name, contrast_name))
  )
  ft <- fisher.test(tab)

  data.frame(
    sub_analysis          = sub_analysis,
    RS1                   = query_name,
    RS2                   = target_name,
    contrast_against      = contrast_name,
    n_query               = length(query_gr),
    n_query_with_target   = query_hit,
    pct_query_with_target = 100 * query_hit / length(query_gr),
    n_contrast            = length(contrast_gr),
    n_contrast_with_target = contrast_hit,
    pct_contrast_with_target = 100 * contrast_hit / length(contrast_gr),
    fisher_or             = unname(ft$estimate),
    fisher_ci_low         = ft$conf.int[1],
    fisher_ci_high        = ft$conf.int[2],
    fisher_p              = ft$p.value,
    stringsAsFactors      = FALSE
  )
}

#' Every interval-level Fisher test behind one sub-analysis.
#'
#' Alist holds exactly two direction sets. Each direction is tested against the
#' other on every target, so both rows of the permutation matrix carry a Fisher
#' odds ratio.
run_direction_fisher_tests <- function(Alist, Blist, sub_analysis) {
  if (length(Alist) != 2) {
    stop("run_direction_fisher_tests() expects two query sets, got ",
         length(Alist))
  }

  rows <- list()
  for (target_name in names(Blist)) {
    for (i in seq_along(Alist)) {
      j <- if (i == 1L) 2L else 1L
      rows[[length(rows) + 1L]] <- contrast_fisher(
        query_gr = Alist[[i]], contrast_gr = Alist[[j]],
        target_gr = Blist[[target_name]],
        query_name = names(Alist)[i], contrast_name = names(Alist)[j],
        target_name = target_name, sub_analysis = sub_analysis
      )
    }
  }
  do.call(rbind, rows)
}

#' Join the permutation table to the Fisher table and classify agreement.
merge_permutation_fisher <- function(perm_tbl, fisher_tbl, q_threshold) {
  merged <- dplyr::inner_join(perm_tbl, fisher_tbl,
                              by = c("sub_analysis", "RS1", "RS2"))
  if (nrow(merged) != nrow(perm_tbl)) {
    stop("Permutation and Fisher tables do not pair one to one: ",
         nrow(perm_tbl), " permutation rows, ", nrow(merged), " joined rows.")
  }

  merged$concordance <- factor(
    dplyr::case_when(
      merged$fisher_p < q_threshold & merged$perm_adj_p_value < q_threshold ~ "Confirmed",
      merged$fisher_p < q_threshold & merged$perm_adj_p_value >= q_threshold ~ "Weakened",
      merged$fisher_p >= q_threshold & merged$perm_adj_p_value < q_threshold ~ "Strengthened",
      TRUE ~ "Concordant NS"
    ),
    levels = CONCORDANCE_LEVELS
  )

  merged$test_id <- paste(merged$RS1, "x", merged$RS2)
  merged$perm_stars <- significance_stars(merged$perm_adj_p_value)
  merged$fisher_stars <- significance_stars(merged$fisher_p)
  merged[order(merged$RS1, -abs(merged$norm_zscore)), , drop = FALSE]
}

#' Count each concordance class within a sub-analysis.
summarise_concordance <- function(assoc) {
  counts <- table(factor(assoc$concordance, levels = CONCORDANCE_LEVELS))
  data.frame(
    sub_analysis = unique(assoc$sub_analysis),
    concordance  = names(counts),
    n_tests      = as.integer(counts),
    pct_tests    = 100 * as.integer(counts) / sum(counts),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# LOCAL Z-SCORE
# =============================================================================

#' The shift offsets that regioneR::localZScore evaluates.
#'
#' localZScore builds them as the step multiples up to floor(window / step) in
#' both directions, with zero in the middle.
build_shift_positions <- function(window, step) {
  n_steps <- floor(window / step)
  c(rev(-seq_len(n_steps) * step), 0L, seq_len(n_steps) * step)
}

#' Run the local z-score profile for one A set against a B list, with caching.
run_local_zscore <- function(A, Blist, genome, ntimes, cores, seed,
                             cache_path, force_rerun) {
  if (file.exists(cache_path) && !force_rerun) {
    cat(sprintf("[40_02 LZ] Loading cached local z-score: %s\n", cache_path))
    return(readRDS(cache_path))
  }

  cat(sprintf("[40_02 LZ] %d region sets, +/- %s bp window, %s bp step, %s permutations\n",
              length(Blist), fmt_int(LZ_WINDOW), fmt_int(LZ_STEP), fmt_int(ntimes)))

  set.seed(seed)
  mlz <- multiLocalZscore(
    A              = A,
    Blist          = Blist,
    genome         = genome,
    ranFUN         = PERM_RANFUN,
    evFUN          = PERM_EVFUN,
    ntimes         = ntimes,
    window         = LZ_WINDOW,
    step           = LZ_STEP,
    force.parallel = TRUE,
    per.chromosome = PERM_PER_CHROMOSOME
  )
  mlz <- makeLZMatrix(mlz)

  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(mlz, cache_path)
  cat(sprintf("[40_02 LZ] Saved cache: %s\n", cache_path))
  mlz
}

#' Long table of the shifted z-score profile, one row per region set and shift.
extract_local_zscore_table <- function(mlz, anchor_set_name) {
  evaluation <- getMultiEvaluation(mlz)
  shifted <- evaluation$shifts
  if (length(shifted) == 0) {
    stop("Local z-score object holds no shifted z-score profiles.")
  }

  shifts <- build_shift_positions(LZ_WINDOW, LZ_STEP)
  profile_lengths <- vapply(shifted, length, integer(1))
  if (any(profile_lengths != length(shifts))) {
    stop("Shifted z-score profiles have ", paste(unique(profile_lengths), collapse = "/"),
         " values but the reconstructed shift positions have ", length(shifts), ".")
  }

  rows <- lapply(names(shifted), function(set_name) {
    data.frame(
      anchor_set     = anchor_set_name,
      region_set     = set_name,
      shift_bp       = shifts,
      shifted_zscore = as.numeric(shifted[[set_name]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# =============================================================================
# FIGURES
# =============================================================================

#' Crosswise association matrix heatmap for one sub-analysis.
plot_crosswise_heatmap <- function(cw, title, subtitle, caption,
                                   out_dir, figure_name, width, height) {
  p <- plotCrosswiseMatrix(cw, matrix_type = "association") +
    labs(title = title, subtitle = subtitle, caption = caption,
         x = NULL, y = NULL, fill = "Normalised\nz-score") +
    theme_emseq() +
    theme(
      axis.text.x = element_text(angle = 40, hjust = 1),
      panel.grid = element_blank()
    )

  save_multiformat_ggplot(p, file.path(out_dir, figure_name),
                          width = width, height = height)
}

#' Forest plot of normalised permutation z-scores for one sub-analysis.
plot_association_forest <- function(assoc, title, subtitle, caption,
                                    out_dir, figure_name, width, height) {
  df <- assoc
  df$point_label <- sprintf("OR = %.2f  %s", df$fisher_or, df$perm_stars)

  p <- ggplot(df, aes(x = norm_zscore, y = reorder(test_id, norm_zscore),
                      color = concordance)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-Z_CUTOFF, Z_CUTOFF), linetype = "dotted",
               color = "grey65") +
    geom_point(size = 3) +
    geom_text(aes(label = point_label), hjust = -0.15, size = 2.6,
              show.legend = FALSE) +
    scale_color_manual(values = CONCORDANCE_COLORS, drop = FALSE,
                       name = "Fisher vs\npermutation") +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.32))) +
    labs(
      title = title,
      subtitle = subtitle,
      caption = caption,
      x = "Normalised permutation z-score",
      y = NULL
    ) +
    theme_emseq() +
    theme(axis.text.y = element_text(size = 8))

  save_multiformat_ggplot(p, file.path(out_dir, figure_name),
                          width = width, height = height)
}

#' Observed against permuted-mean overlap counts for one sub-analysis.
plot_observed_expected <- function(assoc, title, subtitle, caption,
                                   out_dir, figure_name, width, height) {
  long <- rbind(
    data.frame(RS1 = assoc$RS1, RS2 = assoc$RS2, kind = "Observed",
               overlaps = assoc$observed_overlaps, stringsAsFactors = FALSE),
    data.frame(RS1 = assoc$RS1, RS2 = assoc$RS2,
               kind = "Expected (permuted mean)",
               overlaps = assoc$expected_overlaps, stringsAsFactors = FALSE)
  )
  long$kind <- factor(long$kind, levels = names(OBSERVED_EXPECTED_COLORS))

  ratio_labels <- assoc
  ratio_labels$label_y <- pmax(ratio_labels$observed_overlaps,
                               ratio_labels$expected_overlaps)
  ratio_labels$ratio_label <- sprintf("O/E = %.2f", ratio_labels$observed_over_expected)

  p <- ggplot(long, aes(x = RS2, y = overlaps, fill = kind)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    geom_text(aes(label = fmt_int(round(overlaps))),
              position = position_dodge(width = 0.8),
              vjust = -0.35, size = 2.3) +
    geom_text(data = ratio_labels,
              aes(x = RS2, y = label_y, label = ratio_label),
              inherit.aes = FALSE, vjust = -1.9, size = 2.6, fontface = "bold") +
    facet_wrap(~ RS1, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = OBSERVED_EXPECTED_COLORS, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.24))) +
    labs(
      title = title,
      subtitle = subtitle,
      caption = caption,
      x = NULL,
      y = "Overlapping intervals"
    ) +
    theme_emseq() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 30, hjust = 1)
    )

  save_multiformat_ggplot(p, file.path(out_dir, figure_name),
                          width = width, height = height)
}

#' One figure comparing Fisher and permutation across all three sub-analyses.
plot_fisher_vs_permutation <- function(assoc_all, ntimes, out_dir) {
  df <- assoc_all
  df$point_label <- sprintf("OR = %.2f", df$fisher_or)

  counts <- table(factor(df$concordance, levels = CONCORDANCE_LEVELS))
  subtitle <- paste(sprintf("%s: %d", names(counts), as.integer(counts)),
                    collapse = " | ")

  p <- ggplot(df, aes(x = norm_zscore, y = reorder(test_id, norm_zscore),
                      color = concordance)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-Z_CUTOFF, Z_CUTOFF), linetype = "dotted",
               color = "grey65") +
    geom_point(size = 2.6) +
    geom_text(aes(label = point_label), hjust = -0.18, size = 2.2,
              show.legend = FALSE) +
    scale_color_manual(values = CONCORDANCE_COLORS, drop = FALSE,
                       name = "Fisher vs\npermutation") +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.34))) +
    facet_wrap(~ sub_analysis, ncol = 1, scales = "free_y") +
    labs(
      title = "Interval Fisher tests against genomic permutation",
      subtitle = paste0("Point = normalised permutation z-score, label = Fisher odds ratio\n",
                        subtitle),
      caption = sprintf("%s: %s permutations, per.chromosome = TRUE",
                        PERM_RANFUN, fmt_int(ntimes)),
      x = "Normalised permutation z-score",
      y = NULL
    ) +
    theme_emseq() +
    theme(axis.text.y = element_text(size = 7))

  save_multiformat_ggplot(p, file.path(out_dir, "40_02j_fisher_vs_permutation"),
                          width = 13, height = max(11, nrow(df) * 0.30))
}

#' Local z-score profile of the strongest loop anchor association.
plot_local_zscore <- function(mlz, anchor_set_name, target_set_name, ntimes,
                              out_dir) {
  p <- plotSingleLZ(mlz, RS = target_set_name, smoothing = TRUE) +
    labs(
      title = sprintf("Local z-score profile: %s at %s",
                      anchor_set_name, target_set_name),
      subtitle = sprintf("+/- %s bp window, %s bp step, %s permutations",
                         fmt_int(LZ_WINDOW), fmt_int(LZ_STEP), fmt_int(ntimes)),
      x = "Shift of the anchor set (bp)",
      y = "Local z-score"
    ) +
    theme_emseq()

  save_multiformat_ggplot(p, file.path(out_dir, "40_02k_local_zscore_loop_anchors"),
                          width = 10, height = 6.5)
}

#' Every local z-score profile in one faceted figure.
plot_local_zscore_all <- function(lz_table, anchor_set_name, ntimes, out_dir) {
  p <- ggplot(lz_table, aes(x = shift_bp, y = shifted_zscore)) +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey55") +
    geom_line(linewidth = 0.6, color = "#2166AC") +
    facet_wrap(~ region_set, scales = "free_y") +
    scale_x_continuous(labels = function(x) paste0(x / 1000, " kb")) +
    labs(
      title = sprintf("Local z-score profiles: %s against every feature set",
                      anchor_set_name),
      subtitle = sprintf("+/- %s bp window, %s bp step, %s permutations",
                         fmt_int(LZ_WINDOW), fmt_int(LZ_STEP), fmt_int(ntimes)),
      caption = "Zero shift is the observed anchor position.",
      x = "Shift of the anchor set",
      y = "Local z-score"
    ) +
    theme_emseq() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

  save_multiformat_ggplot(p, file.path(out_dir, "40_02l_local_zscore_all_features"),
                          width = 13, height = 9)
}

# =============================================================================
# SUB-ANALYSIS DRIVER
# =============================================================================

#' Run one sub-analysis end to end: permutation, Fisher, tables, figures.
#'
#' @param key short identifier of the sub-analysis, also part of the cache name.
#' @param figure_letters list with heatmap, forest, and observed entries. Each is
#'   the letter that follows "40_02" in that figure's file name.
run_sub_analysis <- function(key, Alist, Blist, genome, opt, out_dir, cache_dir,
                             force_rerun, figure_geometry, figure_letters,
                             titles) {
  label <- unname(SUB_ANALYSIS_LABELS[key])
  cat("\n--------------------------------------------------------------------------------\n")
  cat(label, "\n")
  cat("--------------------------------------------------------------------------------\n")

  # The permutation count is part of the cache file name, so a run with a new
  # --ntimes never reads a cache built at a different count.
  cache_path <- file.path(cache_dir,
                          sprintf("40_02%s_crosswise_n%d.rds", key, opt$ntimes))
  cw <- run_crosswise(Alist, Blist, genome, opt$ntimes, opt$cores, opt$seed,
                      cache_path, label, force_rerun)

  perm_tbl <- extract_permutation_table(cw, label)

  # A cache built from a different region set list produces a result that no
  # longer matches Blist, so the names are checked before the row count.
  absent <- setdiff(names(Blist), unique(perm_tbl$RS2))
  unexpected <- setdiff(unique(perm_tbl$RS2), names(Blist))
  if (length(absent) > 0 || length(unexpected) > 0) {
    stop(label, ": the permutation result does not cover the requested region ",
         "sets.\n  Requested but absent: ",
         if (length(absent) > 0) paste(absent, collapse = ", ") else "none",
         "\n  Present but not requested: ",
         if (length(unexpected) > 0) paste(unexpected, collapse = ", ") else "none",
         "\n  The cache at ", cache_path, " was built from a different region ",
         "set list. Delete it or set FORCE_RERUN=1.")
  }

  expected_pairs <- length(Alist) * length(Blist)
  if (nrow(perm_tbl) != expected_pairs) {
    stop(label, ": expected ", expected_pairs, " permutation rows but got ",
         nrow(perm_tbl),
         ". crosswisePermTest() dropped a region set after an internal error.")
  }
  cat(sprintf("  Permutation rows: %d\n", nrow(perm_tbl)))

  fisher_tbl <- run_direction_fisher_tests(Alist, Blist, label)
  cat(sprintf("  Interval Fisher tests: %d\n", nrow(fisher_tbl)))

  assoc <- merge_permutation_fisher(perm_tbl, fisher_tbl, Q_THRESHOLD)
  print(assoc[, c("RS1", "RS2", "observed_overlaps", "expected_overlaps",
                  "norm_zscore", "perm_adj_p_value", "fisher_or", "fisher_p",
                  "concordance")])

  write_section_table(assoc,
                  file.path(out_dir, sprintf("40_02%s_association_%s.tsv", key, titles$slug)))

  caption <- sprintf("%s, %s permutations, per.chromosome = TRUE",
                     PERM_RANFUN, fmt_int(opt$ntimes))

  plot_crosswise_heatmap(
    cw, titles$heatmap_title, titles$subtitle, caption, out_dir,
    sprintf("40_02%s_crosswise_%s", figure_letters$heatmap, titles$slug),
    figure_geometry$heatmap_width, figure_geometry$heatmap_height
  )

  plot_association_forest(
    assoc, titles$forest_title, titles$subtitle, caption, out_dir,
    sprintf("40_02%s_forest_%s", figure_letters$forest, titles$slug),
    figure_geometry$forest_width, figure_geometry$forest_height
  )

  plot_observed_expected(
    assoc, titles$observed_title, titles$subtitle, caption, out_dir,
    sprintf("40_02%s_observed_expected_%s", figure_letters$observed,
            titles$slug),
    figure_geometry$bar_width, figure_geometry$bar_height
  )

  list(association = assoc, crosswise = cw)
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opt <- parse_options()
  out_dir <- opt$output_dir
  cache_dir <- file.path(out_dir, "cache")
  force_rerun <- nzchar(Sys.getenv("FORCE_RERUN", unset = ""))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  cat("================================================================================\n")
  cat("SECTION 40_02: PERMUTATION VALIDATION OF ATAC PEAKS AND Hi-C LOOP ANCHORS\n")
  cat("================================================================================\n")
  cat("Output directory: ", out_dir, "\n", sep = "")
  cat("Cache directory:  ", cache_dir, "\n", sep = "")
  cat("Permutations:     ", fmt_int(opt$ntimes), "\n", sep = "")
  cat("Cores:            ", opt$cores, "\n", sep = "")
  cat("Seed:             ", opt$seed, "\n", sep = "")
  cat("Randomisation:    ", PERM_RANFUN, " (per.chromosome = ",
      PERM_PER_CHROMOSOME, ")\n", sep = "")
  cat("Evaluation:       ", PERM_EVFUN, "\n", sep = "")
  cat("FORCE_RERUN:      ", if (force_rerun) "set, caches ignored" else "not set",
      "\n\n", sep = "")

  check_inputs()
  patch_chooseHclustMet()

  # --- Genome and region sets ------------------------------------------------
  cat("\n--- Building the genome object ---\n")
  genome <- build_analysis_genome()

  cat("\n--- Loading region sets ---\n")
  marks <- load_mark_sets()
  atac <- load_atac_sets()
  mecp2 <- load_mecp2_sets()
  mch_sets <- build_mch_sets()
  loop_anchors <- build_loop_anchor_sets(HIC_PATHS$loops)

  state_table <- load_chromatin_state_table(HANDOFF_PATHS$chromatin_state)
  body_states <- build_state_region_sets(state_table, "body_state",
                                         BODY_STATE_ORDER, "body state")
  promoter_states <- build_state_region_sets(state_table, "promoter_state",
                                             PROMOTER_STATE_ORDER,
                                             "promoter state")

  loop_features <- c(marks, atac, mecp2, mch_sets)

  region_summary <- rbind(
    summarise_region_sets(atac, "query", unname(SUB_ANALYSIS_LABELS["a"])),
    summarise_region_sets(marks, "target", unname(SUB_ANALYSIS_LABELS["a"])),
    summarise_region_sets(atac, "query", unname(SUB_ANALYSIS_LABELS["b"])),
    summarise_region_sets(body_states, "target", unname(SUB_ANALYSIS_LABELS["b"])),
    summarise_region_sets(atac, "query", unname(SUB_ANALYSIS_LABELS["m"])),
    summarise_region_sets(promoter_states, "target", unname(SUB_ANALYSIS_LABELS["m"])),
    summarise_region_sets(loop_anchors, "query", unname(SUB_ANALYSIS_LABELS["c"])),
    summarise_region_sets(loop_features, "target", unname(SUB_ANALYSIS_LABELS["c"]))
  )
  write_section_table(region_summary, file.path(out_dir, "40_02_region_set_summary.tsv"))

  parameters <- data.frame(
    parameter = c("section", "ntimes", "cores", "seed", "randomisation_function",
                  "evaluation_function", "per_chromosome", "analysis_chromosomes",
                  "fdr_threshold", "local_zscore_window_bp", "local_zscore_step_bp",
                  "force_rerun"),
    value = c(SECTION_ID, opt$ntimes, opt$cores, opt$seed, PERM_RANFUN,
              PERM_EVFUN, PERM_PER_CHROMOSOME,
              paste(ANALYSIS_CHRS, collapse = ","), Q_THRESHOLD,
              LZ_WINDOW, LZ_STEP, force_rerun),
    stringsAsFactors = FALSE
  )
  write_section_table(parameters, file.path(out_dir, "40_02_permutation_parameters.tsv"))

  # --- Sub-analysis A: ATAC against the six marks ----------------------------
  result_a <- run_sub_analysis(
    key = "a", Alist = atac, Blist = marks, genome = genome, opt = opt,
    out_dir = out_dir, cache_dir = cache_dir, force_rerun = force_rerun,
    figure_geometry = list(heatmap_width = 8, heatmap_height = 6.5,
                           forest_width = 11, forest_height = 7,
                           bar_width = 11, bar_height = 8),
    figure_letters = list(heatmap = "a", forest = "d", observed = "g"),
    titles = list(
      slug = "atac_x_marks",
      subtitle = "Differential ATAC peaks against consensus mark peaks",
      heatmap_title = "Permutation association: ATAC peaks and chromatin marks",
      forest_title = "ATAC peaks against chromatin marks",
      observed_title = "Observed and expected ATAC peak overlaps with chromatin marks"
    )
  )

  # --- Sub-analysis B: ATAC against the five gene-body states ----------------
  result_b <- run_sub_analysis(
    key = "b", Alist = atac, Blist = body_states, genome = genome, opt = opt,
    out_dir = out_dir, cache_dir = cache_dir, force_rerun = force_rerun,
    figure_geometry = list(heatmap_width = 8, heatmap_height = 7,
                           forest_width = 11, forest_height = 7.5,
                           bar_width = 12, bar_height = 8),
    figure_letters = list(heatmap = "b", forest = "e", observed = "h"),
    titles = list(
      slug = "atac_x_body_state",
      subtitle = "Differential ATAC peaks against gene bodies grouped by body chromatin state",
      heatmap_title = "Permutation association: ATAC peaks and gene-body chromatin states",
      forest_title = "ATAC peaks against gene-body chromatin states",
      observed_title = "Observed and expected ATAC peak overlaps with gene-body chromatin states"
    )
  )

  # --- Sub-analysis M: ATAC against the four promoter states -----------------
  # The intervals are the same gene bodies as block B. Only the label changes:
  # here a gene is grouped by the chromatin state of its promoter window.
  result_m <- run_sub_analysis(
    key = "m", Alist = atac, Blist = promoter_states, genome = genome, opt = opt,
    out_dir = out_dir, cache_dir = cache_dir, force_rerun = force_rerun,
    figure_geometry = list(heatmap_width = 8, heatmap_height = 6.5,
                           forest_width = 11, forest_height = 7,
                           bar_width = 11, bar_height = 8),
    figure_letters = list(heatmap = "m", forest = "n", observed = "o"),
    titles = list(
      slug = "atac_x_promoter_state",
      subtitle = "Differential ATAC peaks against gene bodies grouped by promoter chromatin state",
      heatmap_title = "Permutation association: ATAC peaks and promoter chromatin states",
      forest_title = "ATAC peaks against promoter chromatin states",
      observed_title = "Observed and expected ATAC peak overlaps with promoter chromatin states"
    )
  )

  # --- Sub-analysis C: loop anchors against every feature set ----------------
  result_c <- run_sub_analysis(
    key = "c", Alist = loop_anchors, Blist = loop_features, genome = genome,
    opt = opt, out_dir = out_dir, cache_dir = cache_dir, force_rerun = force_rerun,
    figure_geometry = list(heatmap_width = 8, heatmap_height = 10,
                           forest_width = 12, forest_height = 10,
                           bar_width = 14, bar_height = 9),
    figure_letters = list(heatmap = "c", forest = "f", observed = "i"),
    titles = list(
      slug = "loop_anchors_x_features",
      subtitle = "Differential Hi-C loop anchors against marks, ATAC, MeCP2, and mCH gene bodies",
      heatmap_title = "Permutation association: loop anchors and chromatin features",
      forest_title = "Hi-C loop anchors against chromatin features",
      observed_title = "Observed and expected loop anchor overlaps with chromatin features"
    )
  )

  # --- Combined comparison ---------------------------------------------------
  cat("\n--------------------------------------------------------------------------------\n")
  cat("Fisher against permutation, all sub-analyses\n")
  cat("--------------------------------------------------------------------------------\n")

  assoc_all <- rbind(result_a$association, result_b$association,
                     result_m$association, result_c$association)
  write_section_table(assoc_all, file.path(out_dir, "40_02_fisher_vs_permutation.tsv"))

  concordance_summary <- rbind(
    summarise_concordance(result_a$association),
    summarise_concordance(result_b$association),
    summarise_concordance(result_m$association),
    summarise_concordance(result_c$association),
    data.frame(
      sub_analysis = "All sub-analyses",
      concordance = CONCORDANCE_LEVELS,
      n_tests = as.integer(table(factor(assoc_all$concordance,
                                        levels = CONCORDANCE_LEVELS))),
      pct_tests = 100 * as.integer(table(factor(assoc_all$concordance,
                                                levels = CONCORDANCE_LEVELS))) /
        nrow(assoc_all),
      stringsAsFactors = FALSE
    )
  )
  print(concordance_summary)
  write_section_table(concordance_summary, file.path(out_dir, "40_02_concordance_summary.tsv"))

  plot_fisher_vs_permutation(assoc_all, opt$ntimes, out_dir)

  # --- Local z-score for the strongest loop anchor association ---------------
  cat("\n--------------------------------------------------------------------------------\n")
  cat("Local z-score profile for the strongest loop anchor association\n")
  cat("--------------------------------------------------------------------------------\n")

  assoc_c <- result_c$association
  strongest <- assoc_c[which.max(abs(assoc_c$norm_zscore)), , drop = FALSE]
  cat(sprintf("  Strongest association: %s x %s (normalised z = %.3f, adj p = %.4g)\n",
              strongest$RS1, strongest$RS2, strongest$norm_zscore,
              strongest$perm_adj_p_value))

  mlz <- run_local_zscore(
    A = loop_anchors[[strongest$RS1]],
    Blist = loop_features,
    genome = genome,
    ntimes = opt$ntimes,
    cores = opt$cores,
    seed = opt$seed,
    cache_path = file.path(cache_dir,
                           sprintf("40_02c_local_zscore_n%d.rds", opt$ntimes)),
    force_rerun = force_rerun
  )

  lz_table <- extract_local_zscore_table(mlz, strongest$RS1)
  write_section_table(lz_table, file.path(out_dir, "40_02_local_zscore_shifts.tsv"))

  lz_summary <- getMultiEvaluation(mlz)$resumeTable
  require_columns(lz_summary, c("name", "z_score", "norm_zscore", "p_value",
                                "adj.p_value"),
                  "Local z-score summary table")
  lz_summary$anchor_set <- strongest$RS1
  write_section_table(lz_summary, file.path(out_dir, "40_02_local_zscore_summary.tsv"))

  plot_local_zscore(mlz, strongest$RS1, strongest$RS2, opt$ntimes, out_dir)
  plot_local_zscore_all(lz_table, strongest$RS1, opt$ntimes, out_dir)

  # --- Summary ---------------------------------------------------------------
  cat("\n================================================================================\n")
  cat("SECTION 40_02 SUMMARY\n")
  cat("================================================================================\n")
  cat(sprintf("Sub-analysis a: %2d x %2d = %3d tests\n", length(atac), length(marks),
              length(atac) * length(marks)))
  cat(sprintf("Sub-analysis b: %2d x %2d = %3d tests\n", length(atac),
              length(body_states), length(atac) * length(body_states)))
  cat(sprintf("Sub-analysis m: %2d x %2d = %3d tests\n", length(atac),
              length(promoter_states), length(atac) * length(promoter_states)))
  cat(sprintf("Sub-analysis c: %2d x %2d = %3d tests\n", length(loop_anchors),
              length(loop_features), length(loop_anchors) * length(loop_features)))
  cat(sprintf("Total pairwise permutation tests: %d\n", nrow(assoc_all)))
  cat(sprintf("Permutations per test: %s\n", fmt_int(opt$ntimes)))
  cat("\nAgreement between the interval Fisher tests and the permutation tests:\n")
  for (level in CONCORDANCE_LEVELS) {
    n <- sum(assoc_all$concordance == level)
    cat(sprintf("  %-16s %3d (%5.1f%%)\n", level, n, 100 * n / nrow(assoc_all)))
  }
  cat("\nSECTION 40_02 COMPLETE\n")
  cat("================================================================================\n\n")
}

main()
