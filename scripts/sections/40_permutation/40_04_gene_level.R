# scripts/sections/40_permutation/40_04_gene_level.R
#
# Section 40_04: label-shuffle validation of every gene-level Fisher test in the
# pipeline.
#
# What this tests
#   Every gene-level 2x2 Fisher test in this pipeline goes through
#   register_fisher_test(), which writes the gene table behind the test and one
#   row into a shared registry. Fisher's exact test assumes the genes are
#   exchangeable. They are not. mCH direction runs in long blocks along a
#   chromosome, and chromatin mark peaks follow gene density, so two gene sets
#   can share chromosomes and produce a large odds ratio without any gene-level
#   association.
#
#   This section reruns each registered test against a chromosome-stratified
#   null. The column labels are shuffled inside each chromosome, which holds the
#   per-chromosome count of column-TRUE genes fixed and holds the row labels
#   fixed. The odds ratio is recomputed on every shuffle to build an empirical
#   null. The observed odds ratio is then read against that null instead of
#   against the hypergeometric null Fisher assumes.
#
#   The case this section exists to catch is a test whose analytic Fisher p is
#   significant and whose empirical p is not. That test is reported as
#   "Weakened".
#
# Test statistic
#   log2 odds ratio with 0.5 added to each of the four cells (the
#   Haldane-Anscombe correction). The correction is applied identically to the
#   observed table and to every shuffled table, so the two are on one scale, and
#   it keeps the statistic finite when a cell is empty. No permutation is
#   dropped for being infinite.
#
#   Because the row labels stay fixed and the shuffle holds the column margin
#   fixed inside each chromosome, both margins of the shuffled table equal the
#   observed margins. The count of row-TRUE and column-TRUE genes therefore
#   fixes the whole table, and the other three cells follow from it.
#
# Registered gene-level Fisher tests audited here (read from the registry at run
# time; this list records which sections produce rows and how their test ids are
# built):
#   10_01  promoter_<promoter state>, body_<gene-body state>
#   10_02  hyper_in_compartment_a, hypo_in_compartment_b,
#          hyper_in_shift_b_to_a, hypo_in_shift_a_to_b
#   10_03  <polycomb definition id>_hyper / _hypo,
#          promoter_state_<promoter state>_hyper / _hypo,
#          body_state_<gene-body state>_hyper / _hypo
#   10_04  a_compartment_sig, b2_sig, label_changed_sig, k27me3_sig, k27ac_sig,
#          a_compartment_hyper
#   20_01  hyper_vs_mecp2_up, hypo_vs_mecp2_down
#   20_02  atac_direction, k27ac_direction, k27me3_direction, k119ub_direction,
#          convergence_2plus
#   20_03  k119ub_up_vs_mecp2_up, mecp2_up_vs_k27ac_up_euchromatin,
#          mecp2_up_vs_k27me3_up_heterochromatin, k119ub_up_vs_mch_hyper,
#          mecp2_up_vs_mch_hyper, k119ub_up_vs_mch_hyper_sig,
#          mecp2_up_vs_mch_hyper_sig
#   20_04  mch_gain_vs_mecp2_gain_<mark category slug>
#   30_01  hyper_x_gained_anchor_<great|direct>,
#          hyper_x_lost_anchor_<great|direct>,
#          hyper_x_k119ub_gained_anchor_<great|direct>
#   30_02  mecp2_gained_at_gained_anchor, mecp2_gained_at_lost_anchor
#   40_01  hyper_x_<peak set>, hypo_x_<peak set>, hyper_vs_hypo_x_<peak set>
#   40_03  mch_hyper_genes_x_<region set>, mch_hypo_genes_x_<region set>
#   50_01  exon_intron_hypo_concordance, intron_enhancer_vs_gene_sig,
#          intron_hypo_vs_gene_hypo
#   60_01  mch_sig_vs_mecp2_sig
#   60_02  mecp2_no_mch_vs_k119ub_signal_gain,
#          mecp2_no_mch_vs_k119ub_peak_gain
#   60_03  collapse_rule_gain_agreement, many_peaks_vs_gain
#   60_04  mut_specific_aging_vs_mch_sig / _mch_hyper / _mch_hypo
#   70_01  neuronal_ctrl_top_quartile, neuronal_ctrl_top_decile,
#          neuronal_mut_top_quartile, neuronal_mut_top_decile,
#          external_ctrl_top_quartile, decile_01 .. decile_10,
#          length_stratum_1 .. length_stratum_5
#   70_02  <mark>_predicted, all_four_predicted, three_plus_predicted,
#          two_plus_predicted_sig
#   70_03  neuronal_vs_mecp2_up, synapse_vs_mecp2_up
#   70_04  <set slug>_x_mecp2_up, <set slug>_both_x_<mark>_remodeled,
#          <set slug>_both_x_k119ub_top_decile,
#          <set slug>_both_x_k119ub_bottom_decile
#
# Reads
#   HANDOFF_PATHS$fisher_registry    the shared registry, written incrementally
#                                    by the sections listed above
#   registry$gene_table_path         one gene table per registered test, with
#                                    gene_name, chr, and the two logical columns
#
# Writes (OUT_DIR defaults to results/sections/40_permutation/)
#   Figures 40_04a .. 40_04e, each in a multi-format subdirectory
#   Tables prefixed 40_04_
#   The cached permutation results as 40_04_gene_level_permutation.rds
#
# Caching
#   The RDS cache is loaded instead of recomputed when it exists and the
#   FORCE_RERUN environment variable is empty. The cache records which tests it
#   covers; the script stops when the registry no longer matches it.
#
# This section runs last, after every section that registers a Fisher test.
#
# This section calls fisher.test() directly rather than register_fisher_test().
# It does not add a new gene-level test; it reproduces the analytic p-value of a
# test another section already registered, so that the analytic and the
# empirical verdict can be compared side by side.
#
# Adapted from Biomodal section 37 (section_37_permutation_gene_level.R). That
# section validated a hardcoded list of 15 tests. Here the list comes from the
# registry, so it always matches the tests that actually ran.

# --- Locate and source shared config ---
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
source(file.path(SECTIONS_DIR, "_shared_config.R"))

library(optparse)
library(gridExtra)

# regioneR is loaded by every 40_ permutation section. The shuffle in this
# section acts on gene labels rather than on genomic intervals, so no genome
# object and no interval randomisation are built here.
library(regioneR)

# =============================================================================
# CONSTANTS
# =============================================================================

SECTION_ID <- "40_04"

OUT_DIR <- OUTPUT_PATHS$permutation

PERM_SEED <- 42L
PERM_DEFAULT_NTIMES <- 10000L
PERM_MIN_NTIMES <- 100L

CACHE_FILENAME <- "40_04_gene_level_permutation.rds"

# Columns register_fisher_test() writes into the registry.
REGISTRY_COLUMNS <- c("section", "test_id", "description", "row_var", "col_var",
                      "gene_table_path", "n_genes", "n_row_true", "n_col_true",
                      "n_both_true", "odds_ratio", "p_value")

# Sections that call register_fisher_test(). The registry decides which tests
# run; this list only reports which of those sections wrote no registry row.
EXPECTED_SECTIONS <- c("10_01", "10_02", "10_03", "10_04",
                       "20_01", "20_02", "20_03", "20_04",
                       "30_01", "30_02",
                       "40_01", "40_03",
                       "50_01",
                       "60_01", "60_02", "60_03", "60_04",
                       "70_01", "70_02", "70_03", "70_04")

CONCORDANCE_ORDER <- c("Confirmed", "Weakened", "Strengthened", "Sign reversed",
                       "Concordant NS", "Degenerate")

CONCORDANCE_COLORS <- c(
  "Confirmed"     = "#2CA02C",
  "Weakened"      = "#D62728",
  "Strengthened"  = "#FF7F0E",
  "Sign reversed" = "#9467BD",
  "Concordant NS" = "grey60",
  "Degenerate"    = "grey30"
)

SERIES_ORDER <- c("Observed", "Permutation-corrected", "Null mean")

SERIES_COLORS <- c(
  "Observed"              = "#000000",
  "Permutation-corrected" = "#D7191C",
  "Null mean"             = "grey55"
)

SERIES_SHAPES <- c(
  "Observed"              = 16,
  "Permutation-corrected" = 17,
  "Null mean"             = 18
)

# Facets and table rows per figure page. Long registries are drawn across
# several pages so that no page grows past a readable size.
HISTOGRAM_FACETS_PER_PAGE <- 20L
HISTOGRAM_FACET_COLUMNS   <- 4L
VIOLINS_PER_PAGE          <- 20L
TABLE_ROWS_PER_PAGE       <- 30L

# Smallest p-value the log10 axis of figure 40_04d accepts.
MIN_P_FOR_LOG <- 1e-300

# Relative tolerance when the recomputed odds ratio is checked against the
# odds ratio the registry recorded.
OR_TOLERANCE <- 1e-6

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
                help = "Number of label shuffles per test [default: %default]")
  )

  parser <- OptionParser(
    option_list = option_list,
    usage = "Rscript 40_04_gene_level.R [--output-dir DIR] [--ntimes N] [N]"
  )
  parsed <- parse_args(parser, positional_arguments = TRUE)
  opt <- parsed$options

  # A bare positional number sets ntimes, matching the other 40_ sections.
  if (length(parsed$args) > 1) {
    stop("At most one positional argument is accepted (ntimes). Got: ",
         paste(parsed$args, collapse = ", "))
  }
  if (length(parsed$args) == 1) {
    positional_ntimes <- suppressWarnings(as.integer(parsed$args[1]))
    if (is.na(positional_ntimes)) {
      stop("The positional argument must be an integer number of shuffles. Got: ",
           parsed$args[1])
    }
    opt$ntimes <- positional_ntimes
  }

  if (opt$ntimes < PERM_MIN_NTIMES) {
    stop("--ntimes must be at least ", PERM_MIN_NTIMES, ", got ", opt$ntimes)
  }
  opt
}

force_rerun_requested <- function() {
  nzchar(Sys.getenv("FORCE_RERUN", unset = ""))
}

# =============================================================================
# SMALL UTILITIES
# =============================================================================

#' Write a section table into out_dir under filename.
#'
#' Joins the path and hands the frame to write_section_table(), which rejects
#' any column holding figure text.
write_tsv_table <- function(df, out_dir, filename) {
  write_section_table(df, file.path(out_dir, filename))
}

fmt_p <- function(p) {
  if (length(p) != 1 || is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

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

test_key <- function(section, test_id) {
  paste(section, test_id, sep = ":")
}

#' Split a vector into pages of at most per_page entries.
paginate <- function(x, per_page) {
  if (length(x) == 0) return(list())
  split(x, ceiling(seq_along(x) / per_page))
}

#' Shorten a description so it fits one row of the summary table figure.
truncate_text <- function(x, width = 58) {
  ifelse(nchar(x) > width, paste0(substr(x, 1, width - 3), "..."), x)
}

#' Benjamini-Hochberg adjustment over the p-values that exist.
#'
#' p.adjust() counts missing values in its default n, which would shrink the
#' adjusted values of the tests that do have a p-value.
adjust_bh <- function(p) {
  p.adjust(p, method = "BH", n = sum(!is.na(p)))
}

# =============================================================================
# REGISTRY
# =============================================================================

#' Read the shared Fisher registry and check every gene table behind it.
load_registry <- function(registry_path) {
  if (!file.exists(registry_path)) {
    stop("The Fisher test registry does not exist: ", registry_path,
         "\nSection 40_04 validates tests other sections registered, so those ",
         "sections must run first: ",
         paste(EXPECTED_SECTIONS, collapse = ", "))
  }

  registry <- read.table(registry_path, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, quote = "", comment.char = "")

  missing_cols <- setdiff(REGISTRY_COLUMNS, colnames(registry))
  if (length(missing_cols) > 0) {
    stop("The Fisher registry ", registry_path, " is missing columns: ",
         paste(missing_cols, collapse = ", "))
  }

  if (nrow(registry) == 0) {
    stop("The Fisher test registry is empty: ", registry_path,
         "\nSection 40_04 validates tests other sections registered, so those ",
         "sections must run first: ",
         paste(EXPECTED_SECTIONS, collapse = ", "))
  }

  registry$key <- test_key(registry$section, registry$test_id)

  duplicated_keys <- registry$key[duplicated(registry$key)]
  if (length(duplicated_keys) > 0) {
    stop("The Fisher registry holds repeated section:test_id keys: ",
         paste(unique(duplicated_keys), collapse = ", "))
  }

  absent <- registry$key[!file.exists(registry$gene_table_path)]
  if (length(absent) > 0) {
    stop("Gene tables named by the registry do not exist for: ",
         paste(absent, collapse = ", "),
         "\nRerun the sections that own these tests.")
  }

  registry <- registry[order(registry$section, registry$test_id), , drop = FALSE]
  rownames(registry) <- NULL

  cat(sprintf("  Registry: %d tests from %d sections\n",
              nrow(registry), length(unique(registry$section))))
  registry
}

#' Compare the sections present in the registry against the sections that
#' contain a register_fisher_test() call.
build_coverage_table <- function(registry) {
  present <- unique(registry$section)
  sections <- union(EXPECTED_SECTIONS, present)

  sections <- sort(sections)
  counts <- vapply(sections, function(s) sum(registry$section == s), integer(1))

  data.frame(
    section = sections,
    calls_register_fisher_test = sections %in% EXPECTED_SECTIONS,
    present_in_registry = sections %in% present,
    n_tests = unname(counts),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# GENE TABLES
# =============================================================================

#' Read the gene table behind one registered test.
#'
#' Returns gene_name, chr, and the two logical columns under the fixed names
#' row_flag and col_flag.
load_gene_table <- function(row) {
  df <- read.table(row$gene_table_path, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "", comment.char = "")

  required <- c("gene_name", "chr", row$row_var, row$col_var)
  missing_cols <- setdiff(required, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Gene table ", row$gene_table_path, " for ", row$key,
         " is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.logical(df[[row$row_var]]) || !is.logical(df[[row$col_var]])) {
    stop("Gene table ", row$gene_table_path, " for ", row$key,
         " must hold logical values in ", row$row_var, " and ", row$col_var,
         ". Read types: ", class(df[[row$row_var]])[1], " and ",
         class(df[[row$col_var]])[1])
  }

  if (anyNA(df[[row$row_var]]) || anyNA(df[[row$col_var]])) {
    stop("Gene table ", row$gene_table_path, " for ", row$key,
         " holds missing values in ", row$row_var, " or ", row$col_var,
         ". register_fisher_test() drops those rows before writing.")
  }

  if (anyNA(df$chr) || any(!nzchar(df$chr))) {
    stop("Gene table ", row$gene_table_path, " for ", row$key,
         " holds an empty chromosome name.")
  }

  if (nrow(df) != row$n_genes) {
    stop("Gene table ", row$gene_table_path, " for ", row$key, " holds ",
         nrow(df), " genes; the registry recorded ", row$n_genes,
         ". The section and the registry are out of step.")
  }

  data.frame(
    gene_name = df$gene_name,
    chr = as.character(df$chr),
    row_flag = df[[row$row_var]],
    col_flag = df[[row$col_var]],
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# CONTINGENCY TABLE AND TEST STATISTIC
# =============================================================================

#' The four cells of the 2x2 table, in the layout register_fisher_test() uses.
contingency_cells <- function(row_flag, col_flag) {
  a <- sum(row_flag & col_flag)
  list(
    a = a,
    b = sum(row_flag) - a,
    c = sum(col_flag) - a,
    d = length(row_flag) - sum(row_flag) - sum(col_flag) + a
  )
}

#' Rebuild the 2x2 matrix register_fisher_test() passed to fisher.test().
#'
#' Rows are row_var TRUE then FALSE; columns are col_var TRUE then FALSE.
contingency_matrix <- function(cells) {
  matrix(c(cells$a, cells$c, cells$b, cells$d), nrow = 2,
         dimnames = list(c("TRUE", "FALSE"), c("TRUE", "FALSE")))
}

#' log2 odds ratio with 0.5 added to every cell.
haldane_log2_or <- function(a, b, c_cell, d) {
  log2(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c_cell + 0.5)))
}

#' The other three cells implied by a, given both fixed margins.
cells_from_a <- function(a, n_row_true, n_col_true, n_genes) {
  list(
    a = a,
    b = n_row_true - a,
    c = n_col_true - a,
    d = n_genes - n_row_true - n_col_true + a
  )
}

#' Whether two odds ratios agree, treating NA and Inf explicitly.
odds_ratios_agree <- function(x, y) {
  if (is.na(x) && is.na(y)) return(TRUE)
  if (is.na(x) || is.na(y)) return(FALSE)
  if (is.infinite(x) || is.infinite(y)) return(identical(x, y))
  abs(x - y) <= OR_TOLERANCE * max(1, abs(x), abs(y))
}

# =============================================================================
# PERMUTATION
# =============================================================================

#' Shuffle the column labels inside each chromosome, ntimes over.
#'
#' The row labels never move and the number of column-TRUE genes on each
#' chromosome never changes, so each shuffle yields a 2x2 table with the same
#' margins as the observed table. Only the count of genes that are TRUE in both
#' columns varies, and that count is what this returns.
#'
#' @param row_flag logical vector, one entry per gene
#' @param col_flag logical vector, one entry per gene
#' @param chr character vector of chromosome names, one entry per gene
#' @param ntimes number of shuffles
#' @return integer vector of length ntimes holding the both-TRUE count
shuffle_within_chromosomes <- function(row_flag, col_flag, chr, ntimes) {
  chr_index <- split(seq_along(chr), chr)
  a_draws <- integer(ntimes)

  for (i in seq_len(ntimes)) {
    shuffled <- col_flag
    for (idx in chr_index) {
      shuffled[idx] <- col_flag[idx][sample.int(length(idx))]
    }
    a_draws[i] <- sum(row_flag & shuffled)
  }
  a_draws
}

#' Per-chromosome composition of one test.
#'
#' expected_a is the count of both-TRUE genes the within-chromosome shuffle
#' produces on average for that chromosome.
chromosome_composition <- function(genes, key, section, test_id) {
  parts <- split(genes, genes$chr)

  rows <- lapply(names(parts), function(chr_name) {
    part <- parts[[chr_name]]
    n_c <- nrow(part)
    r_c <- sum(part$row_flag)
    k_c <- sum(part$col_flag)
    data.frame(
      key = key,
      section = section,
      test_id = test_id,
      chr = chr_name,
      n_genes = n_c,
      n_row_true = r_c,
      n_col_true = k_c,
      n_both_true = sum(part$row_flag & part$col_flag),
      expected_a = r_c * k_c / n_c,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Reason a test cannot be permuted, or "Tested" when it can.
degeneracy_reason <- function(n_genes, n_row_true, n_col_true) {
  if (n_row_true == 0 || n_row_true == n_genes) {
    return("Degenerate: row label is constant")
  }
  if (n_col_true == 0 || n_col_true == n_genes) {
    return("Degenerate: column label is constant")
  }
  "Tested"
}

#' Run the observed test and the label-shuffle null for one registered test.
run_one_test <- function(row, ntimes) {
  genes <- load_gene_table(row)

  n_genes <- nrow(genes)
  n_row_true <- sum(genes$row_flag)
  n_col_true <- sum(genes$col_flag)

  cells <- contingency_cells(genes$row_flag, genes$col_flag)

  # Reproduces the analytic result the owning section registered. The odds ratio
  # is compared against the registry in check_against_registry().
  fisher_result <- fisher.test(contingency_matrix(cells))
  fisher_or <- unname(fisher_result$estimate)

  observed_log2_or <- haldane_log2_or(cells$a, cells$b, cells$c, cells$d)

  status <- degeneracy_reason(n_genes, n_row_true, n_col_true)

  composition <- chromosome_composition(genes, row$key, row$section, row$test_id)
  expected_a_stratified <- sum(composition$expected_a)
  expected_a_unstratified <- n_row_true * n_col_true / n_genes

  result <- list(
    key = row$key,
    section = row$section,
    test_id = row$test_id,
    description = row$description,
    row_var = row$row_var,
    col_var = row$col_var,
    n_genes = n_genes,
    n_row_true = n_row_true,
    n_col_true = n_col_true,
    n_chromosomes = length(unique(genes$chr)),
    cells = cells,
    fisher_or = fisher_or,
    fisher_p = fisher_result$p.value,
    registry_or = row$odds_ratio,
    registry_p = row$p_value,
    observed_log2_or = observed_log2_or,
    expected_a_stratified = expected_a_stratified,
    expected_a_unstratified = expected_a_unstratified,
    composition = composition,
    status = status,
    a_draws = integer(0),
    null_log2_or = numeric(0)
  )

  if (status != "Tested") {
    cat(sprintf("  %-44s %s -- no null distribution built\n", row$key, status))
    return(result)
  }

  set.seed(PERM_SEED)
  started <- Sys.time()
  a_draws <- shuffle_within_chromosomes(genes$row_flag, genes$col_flag,
                                        genes$chr, ntimes)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  null_cells <- cells_from_a(a_draws, n_row_true, n_col_true, n_genes)
  null_log2_or <- haldane_log2_or(null_cells$a, null_cells$b,
                                  null_cells$c, null_cells$d)

  result$a_draws <- a_draws
  result$null_log2_or <- null_log2_or
  result$elapsed_seconds <- elapsed

  if (sd(null_log2_or) == 0) {
    result$status <- "Degenerate: null has zero variance"
    cat(sprintf("  %-44s %s\n", row$key, result$status))
    return(result)
  }

  cat(sprintf("  %-44s n=%7s  obs log2OR=%7.3f  null mean=%7.3f  z=%7.2f  (%.1fs)\n",
              row$key, format(n_genes, big.mark = ","), observed_log2_or,
              mean(null_log2_or),
              (observed_log2_or - mean(null_log2_or)) / sd(null_log2_or),
              elapsed))
  result
}

#' Run every registered test.
run_all_tests <- function(registry, ntimes) {
  cat(sprintf("\nShuffling %s times per test, seed %d, one seed reset per test.\n",
              format(ntimes, big.mark = ","), PERM_SEED))

  results <- vector("list", nrow(registry))
  for (i in seq_len(nrow(registry))) {
    results[[i]] <- run_one_test(registry[i, ], ntimes)
  }
  names(results) <- registry$key
  results
}

#' Load the cached permutation results, or compute and cache them.
load_or_run_permutations <- function(registry, ntimes, cache_path) {
  if (file.exists(cache_path) && !force_rerun_requested()) {
    cat("\nLoading cached permutation results from: ", cache_path, "\n", sep = "")
    cached <- readRDS(cache_path)

    required <- c("results", "ntimes", "seed", "keys")
    missing_entries <- setdiff(required, names(cached))
    if (length(missing_entries) > 0) {
      stop("Cache ", cache_path, " is missing entries: ",
           paste(missing_entries, collapse = ", "),
           ". Set FORCE_RERUN=1 to rebuild it.")
    }

    added <- setdiff(registry$key, cached$keys)
    removed <- setdiff(cached$keys, registry$key)
    if (length(added) > 0 || length(removed) > 0) {
      stop("The cache at ", cache_path, " covers a different set of tests than ",
           "the registry.\n  In the registry but not the cache: ",
           if (length(added) > 0) paste(added, collapse = ", ") else "none",
           "\n  In the cache but not the registry: ",
           if (length(removed) > 0) paste(removed, collapse = ", ") else "none",
           "\nSet FORCE_RERUN=1 to rebuild the cache.")
    }

    cat(sprintf("  Cache holds %d tests at %s shuffles each.\n",
                length(cached$results), format(cached$ntimes, big.mark = ",")))
    return(cached)
  }

  if (force_rerun_requested()) {
    cat("\nFORCE_RERUN is set; recomputing every permutation.\n")
  }

  results <- run_all_tests(registry, ntimes)

  cached <- list(
    results = results,
    keys = registry$key,
    ntimes = ntimes,
    seed = PERM_SEED,
    generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  saveRDS(cached, cache_path)
  cat("  Saved cache: ", cache_path, "\n", sep = "")
  cached
}

# =============================================================================
# RESULT TABLES
# =============================================================================

#' Stop when a recomputed odds ratio disagrees with the registry.
check_against_registry <- function(results, registry) {
  cat("\nChecking recomputed odds ratios against the registry...\n")

  mismatched <- character(0)
  for (key in registry$key) {
    r <- results[[key]]
    if (!odds_ratios_agree(r$fisher_or, r$registry_or)) {
      mismatched <- c(mismatched, sprintf("%s (recomputed %.6g, registry %.6g)",
                                          key, r$fisher_or, r$registry_or))
    }
  }

  if (length(mismatched) > 0) {
    stop("The odds ratio recomputed from the gene table differs from the ",
         "registry for:\n  ", paste(mismatched, collapse = "\n  "),
         "\nThe gene tables and the registry come from different runs. Rerun ",
         "the owning sections.")
  }

  cat(sprintf("  All %d odds ratios reproduce the registry.\n", nrow(registry)))
  invisible(TRUE)
}

#' One summary row per test.
build_summary_table <- function(results, ntimes) {
  rows <- lapply(results, function(r) {
    tested <- r$status == "Tested"

    null_mean <- if (tested) mean(r$null_log2_or) else NA_real_
    null_sd <- if (tested) sd(r$null_log2_or) else NA_real_
    null_median <- if (tested) median(r$null_log2_or) else NA_real_
    null_q025 <- if (tested) unname(quantile(r$null_log2_or, 0.025)) else NA_real_
    null_q975 <- if (tested) unname(quantile(r$null_log2_or, 0.975)) else NA_real_

    z_score <- if (tested) (r$observed_log2_or - null_mean) / null_sd else NA_real_
    corrected_log2_or <- if (tested) r$observed_log2_or - null_mean else NA_real_

    deviation_obs <- if (tested) abs(r$observed_log2_or - null_mean) else NA_real_
    deviation_null <- if (tested) abs(r$null_log2_or - null_mean) else numeric(0)

    n_draws <- length(r$null_log2_or)
    empirical_p <- if (tested) {
      (sum(deviation_null >= deviation_obs) + 1) / (n_draws + 1)
    } else NA_real_
    empirical_p_right <- if (tested) {
      (sum(r$null_log2_or >= r$observed_log2_or) + 1) / (n_draws + 1)
    } else NA_real_
    empirical_p_left <- if (tested) {
      (sum(r$null_log2_or <= r$observed_log2_or) + 1) / (n_draws + 1)
    } else NA_real_
    empirical_p_doubled <- if (tested) {
      min(1, 2 * min(empirical_p_left, empirical_p_right))
    } else NA_real_

    data.frame(
      key = r$key,
      section = r$section,
      test_id = r$test_id,
      description = r$description,
      row_var = r$row_var,
      col_var = r$col_var,
      status = r$status,
      n_genes = r$n_genes,
      n_chromosomes = r$n_chromosomes,
      n_row_true = r$n_row_true,
      n_col_true = r$n_col_true,
      cell_a = r$cells$a,
      cell_b = r$cells$b,
      cell_c = r$cells$c,
      cell_d = r$cells$d,
      expected_a_unstratified = r$expected_a_unstratified,
      expected_a_stratified = r$expected_a_stratified,
      fisher_or = r$fisher_or,
      fisher_p = r$fisher_p,
      observed_or = 2^r$observed_log2_or,
      observed_log2_or = r$observed_log2_or,
      n_shuffles = n_draws,
      null_mean_log2_or = null_mean,
      null_sd_log2_or = null_sd,
      null_median_log2_or = null_median,
      null_q025_log2_or = null_q025,
      null_q975_log2_or = null_q975,
      null_mean_or = 2^null_mean,
      z_score = z_score,
      empirical_p = empirical_p,
      empirical_p_left = empirical_p_left,
      empirical_p_right = empirical_p_right,
      empirical_p_doubled = empirical_p_doubled,
      corrected_log2_or = corrected_log2_or,
      corrected_or = 2^corrected_log2_or,
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, rows)
  rownames(summary_df) <- NULL

  summary_df$empirical_p_resolution <- 1 / (summary_df$n_shuffles + 1)
  summary_df$fisher_p_bh <- adjust_bh(summary_df$fisher_p)
  summary_df$empirical_p_bh <- adjust_bh(summary_df$empirical_p)

  summary_df$concordance <- classify_concordance(summary_df)
  summary_df$concordance <- factor(summary_df$concordance,
                                   levels = CONCORDANCE_ORDER)
  summary_df$empirical_stars <- sig_stars(summary_df$empirical_p)
  summary_df$fisher_stars <- sig_stars(summary_df$fisher_p)
  summary_df$n_shuffles_requested <- ntimes

  summary_df
}

#' Compare the analytic and the empirical verdict for each test.
#'
#' "Weakened" is the case this section exists to catch: Fisher calls the test
#' significant, the chromosome-stratified null does not.
classify_concordance <- function(df) {
  fisher_sig <- df$fisher_p < Q_THRESHOLD
  empirical_sig <- df$empirical_p < Q_THRESHOLD
  sign_flip <- sign(df$observed_log2_or) != sign(df$corrected_log2_or)

  dplyr::case_when(
    df$status != "Tested" ~ "Degenerate",
    fisher_sig & empirical_sig & sign_flip ~ "Sign reversed",
    fisher_sig & empirical_sig ~ "Confirmed",
    fisher_sig & !empirical_sig ~ "Weakened",
    !fisher_sig & empirical_sig ~ "Strengthened",
    TRUE ~ "Concordant NS"
  )
}

#' Every shuffled statistic, one row per draw.
build_null_draws_table <- function(results) {
  rows <- lapply(results, function(r) {
    if (length(r$null_log2_or) == 0) return(NULL)
    data.frame(
      key = r$key,
      section = r$section,
      test_id = r$test_id,
      draw = seq_along(r$null_log2_or),
      both_true = r$a_draws,
      null_log2_or = r$null_log2_or,
      stringsAsFactors = FALSE
    )
  })

  kept <- rows[!vapply(rows, is.null, logical(1))]
  if (length(kept) == 0) {
    stop("No registered test produced a null distribution. Every test in the ",
         "registry has a constant row or column label, so none can be ",
         "validated by shuffling.")
  }
  do.call(rbind, kept)
}

#' Chromosome composition of every test, stacked.
build_composition_table <- function(results) {
  do.call(rbind, lapply(results, function(r) r$composition))
}

#' Count of tests in each concordance class.
build_concordance_counts <- function(summary_df) {
  counts <- table(summary_df$concordance)
  data.frame(
    concordance = names(counts),
    n_tests = as.integer(counts),
    fraction = as.numeric(counts) / nrow(summary_df),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# FIGURES
# =============================================================================

#' Long-format table for the forest plot: three series per test.
build_forest_data <- function(summary_df) {
  tested <- summary_df[summary_df$status == "Tested", , drop = FALSE]
  tested <- tested[order(tested$corrected_log2_or), , drop = FALSE]
  key_levels <- tested$key

  observed <- data.frame(
    key = tested$key,
    series = "Observed",
    log2_or = tested$observed_log2_or,
    low = NA_real_,
    high = NA_real_,
    stringsAsFactors = FALSE
  )

  corrected <- data.frame(
    key = tested$key,
    series = "Permutation-corrected",
    log2_or = tested$corrected_log2_or,
    low = NA_real_,
    high = NA_real_,
    stringsAsFactors = FALSE
  )

  null_mean <- data.frame(
    key = tested$key,
    series = "Null mean",
    log2_or = tested$null_mean_log2_or,
    low = tested$null_q025_log2_or,
    high = tested$null_q975_log2_or,
    stringsAsFactors = FALSE
  )

  out <- rbind(observed, corrected, null_mean)
  out$key <- factor(out$key, levels = key_levels)
  out$series <- factor(out$series, levels = SERIES_ORDER)
  out
}

#' Verdict text drawn to the right of each row of the forest plot.
build_forest_labels <- function(summary_df, key_levels, label_x) {
  tested <- summary_df[match(key_levels, summary_df$key), , drop = FALSE]
  data.frame(
    key = factor(tested$key, levels = key_levels),
    log2_or = label_x,
    label = sprintf("%s  (emp p = %.4f, z = %.1f)",
                    as.character(tested$concordance),
                    tested$empirical_p, tested$z_score),
    stringsAsFactors = FALSE
  )
}

plot_forest <- function(summary_df, ntimes, out_dir) {
  cat("--- Figure 40_04a: observed against permutation-corrected odds ratios ---\n")

  forest <- build_forest_data(summary_df)
  if (nrow(forest) == 0) {
    stop("No registered test could be permuted, so the forest plot has nothing ",
         "to draw. Every test in the registry is degenerate.")
  }

  key_levels <- levels(forest$key)
  n_tests <- length(key_levels)

  value_range <- range(c(forest$log2_or, forest$low, forest$high), na.rm = TRUE)
  span <- max(diff(value_range), 1)
  label_x <- value_range[2] + 0.06 * span
  labels <- build_forest_labels(summary_df, key_levels, label_x)

  p <- ggplot(forest, aes(x = log2_or, y = key)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_errorbarh(aes(xmin = low, xmax = high), height = 0.25,
                   colour = SERIES_COLORS[["Null mean"]], linewidth = 0.6,
                   na.rm = TRUE) +
    geom_point(aes(colour = series, shape = series), size = 2.6) +
    geom_text(data = labels, aes(x = log2_or, y = key, label = label),
              inherit.aes = FALSE, hjust = 0, size = 2.4, colour = "grey20") +
    expand_limits(x = label_x + 0.70 * span) +
    scale_colour_manual(values = SERIES_COLORS, name = NULL, drop = FALSE) +
    scale_shape_manual(values = SERIES_SHAPES, name = NULL, drop = FALSE) +
    labs(
      title = "Gene-level odds ratios before and after the chromosome-stratified null",
      subtitle = sprintf(paste("Corrected log2 odds ratio is observed minus the",
                               "null mean. Bars span the 2.5th to 97.5th",
                               "percentile of %s shuffles."),
                         format(ntimes, big.mark = ",")),
      x = "log2 odds ratio (0.5 added to each cell)",
      y = NULL
    ) +
    theme_emseq() +
    theme(legend.position = "top",
          axis.text.y = element_text(size = 7))

  save_multiformat_ggplot(p, file.path(out_dir, "40_04a_odds_ratio_forest"),
                          width = 15, height = max(7, n_tests * 0.30 + 3))
}

plot_null_histograms <- function(summary_df, draws, ntimes, out_dir) {
  cat("--- Figure 40_04b: null distributions with the observed value marked ---\n")

  tested <- summary_df[summary_df$status == "Tested", , drop = FALSE]
  tested <- tested[order(tested$section, tested$test_id), , drop = FALSE]
  pages <- paginate(tested$key, HISTOGRAM_FACETS_PER_PAGE)

  for (page_index in seq_along(pages)) {
    page_keys <- pages[[page_index]]
    page_draws <- draws[draws$key %in% page_keys, , drop = FALSE]
    page_stats <- tested[tested$key %in% page_keys, , drop = FALSE]

    page_draws$key <- factor(page_draws$key, levels = page_keys)
    page_stats$key <- factor(page_stats$key, levels = page_keys)

    page_stats$text_label <- sprintf(
      "n genes = %s\nobs = %.2f\nnull mean = %.2f\nz = %.2f\nemp %s",
      fmt_count(page_stats$n_genes),
      page_stats$observed_log2_or,
      page_stats$null_mean_log2_or,
      page_stats$z_score,
      vapply(page_stats$empirical_p, fmt_p, character(1)))

    p <- ggplot(page_draws, aes(x = null_log2_or)) +
      geom_histogram(bins = 50, fill = "grey75", colour = "grey40",
                     linewidth = 0.15) +
      geom_vline(data = page_stats, aes(xintercept = null_mean_log2_or),
                 linetype = "dashed", colour = "grey25", linewidth = 0.5) +
      geom_vline(data = page_stats, aes(xintercept = observed_log2_or),
                 colour = "#D7191C", linewidth = 0.9) +
      geom_text(data = page_stats, aes(x = -Inf, y = Inf, label = text_label),
                inherit.aes = FALSE, hjust = -0.06, vjust = 1.1, size = 2.3,
                lineheight = 0.95) +
      facet_wrap(~ key, ncol = HISTOGRAM_FACET_COLUMNS, scales = "free") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.40))) +
      labs(
        title = sprintf("Label-shuffle null distributions (page %d of %d)",
                        page_index, length(pages)),
        subtitle = sprintf(paste("Grey bars are %s within-chromosome shuffles.",
                                 "Red line is the observed log2 odds ratio;",
                                 "dashed line is the null mean."),
                           format(ntimes, big.mark = ",")),
        x = "log2 odds ratio per shuffle",
        y = "Shuffles"
      ) +
      theme_emseq() +
      theme(strip.text = element_text(size = 7))

    n_rows <- ceiling(length(page_keys) / HISTOGRAM_FACET_COLUMNS)
    save_multiformat_ggplot(
      p,
      file.path(out_dir, sprintf("40_04b_null_distributions_page%02d", page_index)),
      width = 17, height = max(6, n_rows * 2.7 + 2))
  }
}

plot_null_violins <- function(summary_df, draws, ntimes, out_dir) {
  cat("--- Figure 40_04c: null spread per test with the observed value marked ---\n")

  tested <- summary_df[summary_df$status == "Tested", , drop = FALSE]
  tested <- tested[order(tested$section, tested$test_id), , drop = FALSE]
  pages <- paginate(tested$key, VIOLINS_PER_PAGE)

  for (page_index in seq_along(pages)) {
    page_keys <- pages[[page_index]]
    page_draws <- draws[draws$key %in% page_keys, , drop = FALSE]
    page_stats <- tested[tested$key %in% page_keys, , drop = FALSE]

    page_draws$key <- factor(page_draws$key, levels = page_keys)
    page_stats$key <- factor(page_stats$key, levels = page_keys)

    # The annotation frame is plot data only. It carries the n and median that
    # every violin figure in this pipeline shows, and it is never written.
    annotation <- summarise_groups(page_draws, "key", "null_log2_or")
    annotation$key <- factor(annotation$key, levels = page_keys)
    annotation$label <- group_label(annotation)

    y_span <- range(c(page_draws$null_log2_or, page_stats$observed_log2_or))
    annotation$label_y <- y_span[2] + 0.16 * diff(y_span)

    p <- ggplot(page_draws, aes(x = key, y = null_log2_or)) +
      geom_violin(fill = "grey80", colour = "grey40", linewidth = 0.3,
                  scale = "width") +
      geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white",
                   linewidth = 0.3) +
      geom_point(data = page_stats, aes(x = key, y = observed_log2_or),
                 inherit.aes = FALSE, colour = "#D7191C", size = 2.6) +
      geom_text(data = annotation, aes(x = key, y = label_y, label = label),
                inherit.aes = FALSE, size = 2.4, lineheight = 0.95) +
      expand_limits(y = y_span[2] + 0.30 * diff(y_span)) +
      labs(
        title = sprintf("Null spread against the observed statistic (page %d of %d)",
                        page_index, length(pages)),
        subtitle = sprintf(paste("Violins are %s within-chromosome shuffles per",
                                 "test. Red point is the observed log2 odds",
                                 "ratio. n and median describe the shuffles."),
                           format(ntimes, big.mark = ",")),
        x = NULL,
        y = "log2 odds ratio"
      ) +
      theme_emseq() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

    save_multiformat_ggplot(
      p,
      file.path(out_dir, sprintf("40_04c_null_spread_page%02d", page_index)),
      width = max(10, length(page_keys) * 0.60 + 3), height = 9)
  }
}

plot_analytic_vs_empirical <- function(summary_df, ntimes, out_dir) {
  cat("--- Figure 40_04d: analytic against empirical p-value ---\n")

  df <- summary_df[summary_df$status == "Tested", , drop = FALSE]
  df$fisher_p_plot <- pmax(df$fisher_p, MIN_P_FOR_LOG)
  df$empirical_p_plot <- df$empirical_p

  disagreeing <- df[df$concordance %in%
                      c("Weakened", "Strengthened", "Sign reversed"), , drop = FALSE]

  resolution <- 1 / (ntimes + 1)

  p <- ggplot(df, aes(x = fisher_p_plot, y = empirical_p_plot,
                      colour = concordance)) +
    geom_hline(yintercept = Q_THRESHOLD, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = Q_THRESHOLD, linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = resolution, linetype = "dotted", colour = "grey60") +
    geom_abline(slope = 1, intercept = 0, colour = "grey75", linewidth = 0.4) +
    geom_point(size = 2.8, alpha = 0.9) +
    ggrepel::geom_text_repel(data = disagreeing, aes(label = key),
                             size = 2.5, max.overlaps = 40,
                             min.segment.length = 0, show.legend = FALSE) +
    scale_x_log10(labels = scales::label_scientific()) +
    scale_y_log10(labels = scales::label_scientific()) +
    scale_colour_manual(values = CONCORDANCE_COLORS, name = "Verdict",
                        drop = FALSE) +
    annotate("text", x = Q_THRESHOLD, y = resolution, hjust = -0.05, vjust = -0.6,
             size = 2.8, colour = "grey35",
             label = sprintf("empirical p floor = 1/(%s + 1)",
                             format(ntimes, big.mark = ","))) +
    labs(
      title = "Analytic Fisher p against the label-shuffle empirical p",
      subtitle = sprintf(paste("Dashed lines mark p = %.2f. Points below and",
                               "right of the crossing are significant by Fisher",
                               "only. Fisher p is floored at %.0e for the axis."),
                         Q_THRESHOLD, MIN_P_FOR_LOG),
      x = "Fisher exact p-value",
      y = "Empirical p-value from the chromosome-stratified shuffle"
    ) +
    theme_emseq()

  save_multiformat_ggplot(p, file.path(out_dir, "40_04d_analytic_vs_empirical_p"),
                          width = 12, height = 9)
}

#' Rows of the summary table figure.
build_table_figure_data <- function(summary_df) {
  data.frame(
    Section = summary_df$section,
    Test = summary_df$test_id,
    Description = truncate_text(summary_df$description),
    N = fmt_count(summary_df$n_genes),
    `Fisher OR` = ifelse(is.na(summary_df$fisher_or), "-",
                         sprintf("%.2f", summary_df$fisher_or)),
    `Fisher p` = ifelse(is.na(summary_df$fisher_p), "-",
                        sprintf("%.2e", summary_df$fisher_p)),
    `Obs log2OR` = sprintf("%.2f", summary_df$observed_log2_or),
    `Null mean` = ifelse(is.na(summary_df$null_mean_log2_or), "-",
                         sprintf("%.2f", summary_df$null_mean_log2_or)),
    `Corr log2OR` = ifelse(is.na(summary_df$corrected_log2_or), "-",
                           sprintf("%.2f", summary_df$corrected_log2_or)),
    Z = ifelse(is.na(summary_df$z_score), "-",
               sprintf("%.2f", summary_df$z_score)),
    `Emp p` = ifelse(is.na(summary_df$empirical_p), "-",
                     sprintf("%.4f", summary_df$empirical_p)),
    Verdict = as.character(summary_df$concordance),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

plot_summary_table <- function(summary_df, ntimes, out_dir) {
  cat("--- Figure 40_04e: summary table ---\n")

  ordered <- summary_df[order(summary_df$section, summary_df$test_id), ,
                        drop = FALSE]
  table_df <- build_table_figure_data(ordered)
  pages <- paginate(seq_len(nrow(table_df)), TABLE_ROWS_PER_PAGE)

  table_theme <- gridExtra::ttheme_default(
    core = list(fg_params = list(fontsize = 7, hjust = 0, x = 0.02)),
    colhead = list(fg_params = list(fontsize = 8, fontface = "bold"))
  )

  for (page_index in seq_along(pages)) {
    page_rows <- pages[[page_index]]
    page_df <- table_df[page_rows, , drop = FALSE]

    table_grob <- gridExtra::tableGrob(page_df, rows = NULL, theme = table_theme)

    p <- ggplot() +
      annotation_custom(table_grob) +
      labs(
        title = sprintf("Fisher against label shuffle (page %d of %d)",
                        page_index, length(pages)),
        subtitle = sprintf(paste("%s within-chromosome shuffles per test, seed",
                                 "%d. Corr log2OR is observed minus the null",
                                 "mean."),
                           format(ntimes, big.mark = ","), PERM_SEED)
      ) +
      theme_void() +
      theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
            plot.subtitle = element_text(size = 10, hjust = 0.5),
            plot.margin = margin(12, 12, 12, 12))

    save_multiformat_ggplot(
      p,
      file.path(out_dir, sprintf("40_04e_summary_table_page%02d", page_index)),
      width = 20, height = max(6, nrow(page_df) * 0.28 + 3))
  }
}

# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

print_coverage_summary <- function(coverage) {
  cat("\n--- Registry coverage ---\n")
  for (i in seq_len(nrow(coverage))) {
    row <- coverage[i, ]
    marker <- if (row$present_in_registry) "present" else "MISSING"
    cat(sprintf("  %-8s %-8s %4d tests\n", row$section, marker, row$n_tests))
  }

  missing_sections <- coverage$section[coverage$calls_register_fisher_test &
                                         !coverage$present_in_registry]
  if (length(missing_sections) > 0) {
    cat("\n  These sections register Fisher tests but wrote no registry row: ",
        paste(missing_sections, collapse = ", "), "\n", sep = "")
    cat("  Their tests are not validated in this run.\n")
  }

  unexpected <- coverage$section[!coverage$calls_register_fisher_test &
                                   coverage$present_in_registry]
  if (length(unexpected) > 0) {
    cat("\n  Registry rows come from sections this script did not expect: ",
        paste(unexpected, collapse = ", "), "\n", sep = "")
    cat("  Add them to EXPECTED_SECTIONS in this script.\n")
  }
}

print_concordance_summary <- function(summary_df) {
  cat("\n--- Fisher against label shuffle ---\n")
  for (level in CONCORDANCE_ORDER) {
    n <- sum(summary_df$concordance == level)
    cat(sprintf("  %-15s %3d / %3d (%.0f%%)\n", paste0(level, ":"), n,
                nrow(summary_df), 100 * n / nrow(summary_df)))
  }
}

print_weakened_tests <- function(summary_df) {
  weakened <- summary_df[summary_df$concordance == "Weakened", , drop = FALSE]
  cat("\n--- Tests significant by Fisher but not by the shuffle ---\n")

  if (nrow(weakened) == 0) {
    cat("  None.\n")
    return(invisible(NULL))
  }

  weakened <- weakened[order(weakened$empirical_p, decreasing = TRUE), ,
                       drop = FALSE]
  for (i in seq_len(nrow(weakened))) {
    row <- weakened[i, ]
    cat(sprintf("  %-44s Fisher %s | empirical p = %.4f | obs log2OR %.2f | null mean %.2f\n",
                row$key, fmt_p(row$fisher_p), row$empirical_p,
                row$observed_log2_or, row$null_mean_log2_or))
    cat(sprintf("      %s\n", row$description))
  }
}

print_degenerate_tests <- function(summary_df) {
  degenerate <- summary_df[summary_df$status != "Tested", , drop = FALSE]
  if (nrow(degenerate) == 0) return(invisible(NULL))

  cat("\n--- Tests with no usable null distribution ---\n")
  for (i in seq_len(nrow(degenerate))) {
    row <- degenerate[i, ]
    cat(sprintf("  %-44s %s (n = %s, row TRUE = %s, col TRUE = %s)\n",
                row$key, row$status, fmt_count(row$n_genes),
                fmt_count(row$n_row_true), fmt_count(row$n_col_true)))
  }
}

print_largest_corrections <- function(summary_df, n_show = 10) {
  tested <- summary_df[summary_df$status == "Tested", , drop = FALSE]
  if (nrow(tested) == 0) return(invisible(NULL))

  tested$shift <- abs(tested$null_mean_log2_or)
  tested <- tested[order(tested$shift, decreasing = TRUE), , drop = FALSE]
  n_show <- min(n_show, nrow(tested))

  cat("\n--- Largest chromosome-composition shifts in the null ---\n")
  for (i in seq_len(n_show)) {
    row <- tested[i, ]
    cat(sprintf("  %-44s null mean log2OR = %6.3f | expected both-TRUE stratified %.1f vs unstratified %.1f\n",
                row$key, row$null_mean_log2_or,
                row$expected_a_stratified, row$expected_a_unstratified))
  }
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  opt <- parse_section_options()
  out_dir <- opt$output_dir
  ntimes <- opt$ntimes
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n")
  cat("================================================================================\n")
  cat(sprintf("SECTION %s: LABEL-SHUFFLE VALIDATION OF EVERY GENE-LEVEL FISHER TEST\n",
              SECTION_ID))
  cat("================================================================================\n")
  cat(sprintf("Output directory: %s\n", out_dir))
  cat(sprintf("Registry:         %s\n", HANDOFF_PATHS$fisher_registry))
  cat(sprintf("Shuffles:         %s\n", format(ntimes, big.mark = ",")))
  cat(sprintf("Seed:             %d\n", PERM_SEED))
  cat(sprintf("regioneR:         %s\n", packageVersion("regioneR")))
  cat("\n")

  # --- Step 1: registry ------------------------------------------------------
  cat("STEP 1: Registry\n")
  registry <- load_registry(HANDOFF_PATHS$fisher_registry)
  coverage <- build_coverage_table(registry)
  print_coverage_summary(coverage)

  write_tsv_table(registry, out_dir, "40_04_registry_snapshot.tsv")
  write_tsv_table(coverage, out_dir, "40_04_registry_coverage.tsv")

  # --- Step 2: permutation ---------------------------------------------------
  cat("\nSTEP 2: Chromosome-stratified label shuffle\n")
  cache_path <- file.path(out_dir, CACHE_FILENAME)
  cached <- load_or_run_permutations(registry, ntimes, cache_path)
  results <- cached$results

  if (cached$ntimes != ntimes) {
    cat(sprintf("  The cache holds %s shuffles; --ntimes asked for %s. Reporting the cached run.\n",
                format(cached$ntimes, big.mark = ","),
                format(ntimes, big.mark = ",")))
    ntimes <- cached$ntimes
  }

  # --- Step 3: integrity check ----------------------------------------------
  cat("\nSTEP 3: Integrity check\n")
  check_against_registry(results, registry)

  # --- Step 4: result tables -------------------------------------------------
  cat("\nSTEP 4: Result tables\n")
  summary_df <- build_summary_table(results, ntimes)
  draws <- build_null_draws_table(results)
  composition <- build_composition_table(results)
  concordance_counts <- build_concordance_counts(summary_df)

  write_tsv_table(summary_df, out_dir, "40_04_permutation_summary.tsv")
  write_tsv_table(concordance_counts, out_dir, "40_04_concordance_counts.tsv")
  write_tsv_table(composition, out_dir, "40_04_chromosome_composition.tsv")
  write_tsv_table(draws, out_dir, "40_04_null_draws.tsv")

  print_concordance_summary(summary_df)
  print_weakened_tests(summary_df)
  print_degenerate_tests(summary_df)
  print_largest_corrections(summary_df)

  # --- Step 5: figures -------------------------------------------------------
  cat("\nSTEP 5: Figures\n")
  plot_forest(summary_df, ntimes, out_dir)
  plot_null_histograms(summary_df, draws, ntimes, out_dir)
  plot_null_violins(summary_df, draws, ntimes, out_dir)
  plot_analytic_vs_empirical(summary_df, ntimes, out_dir)
  plot_summary_table(summary_df, ntimes, out_dir)

  # --- Step 6: summary -------------------------------------------------------
  cat("\n")
  cat("================================================================================\n")
  cat(sprintf("SECTION %s SUMMARY\n", SECTION_ID))
  cat("================================================================================\n")
  cat(sprintf("Registered tests validated:  %d\n", nrow(summary_df)))
  cat(sprintf("Sections covered:            %d of %d that register tests\n",
              sum(coverage$calls_register_fisher_test &
                    coverage$present_in_registry),
              sum(coverage$calls_register_fisher_test)))
  cat(sprintf("Shuffles per test:           %s\n", format(ntimes, big.mark = ",")))
  cat(sprintf("Confirmed:                   %d\n",
              sum(summary_df$concordance == "Confirmed")))
  cat(sprintf("Weakened (Fisher only):      %d\n",
              sum(summary_df$concordance == "Weakened")))
  cat(sprintf("Strengthened (shuffle only): %d\n",
              sum(summary_df$concordance == "Strengthened")))
  cat(sprintf("Sign reversed:               %d\n",
              sum(summary_df$concordance == "Sign reversed")))
  cat(sprintf("Concordant, not significant: %d\n",
              sum(summary_df$concordance == "Concordant NS")))
  cat(sprintf("Degenerate:                  %d\n",
              sum(summary_df$concordance == "Degenerate")))
  cat(sprintf("Cache:                       %s\n", cache_path))
  cat(sprintf("Figures and tables written to: %s\n", out_dir))
  cat(sprintf("Section %s complete.\n", SECTION_ID))
}

main()
