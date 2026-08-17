# scripts/dmrseq/03_gather_results.R
#
# Phase 3 of the dmrseq pipeline: load per-chromosome DMR results (RDS files
# from Phase 2), combine into a single GRanges object, and write the final
# TSV table and summary file.
#
# Runs once per mode (local, blocks). Each invocation reads from
# results/dmrseq/.per_chr/{mode}/ and writes mch_dmrs_{mode}.tsv.
#
# Usage:
#   Rscript 03_gather_results.R \
#       --results-dir /path/to/results/dmrseq \
#       --mode local \
#       --split-dir /path/to/.dmrseq_splits \
#       --clean TRUE

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
})

CANONICAL_CHRS <- paste0("chr", c(1:19, "X", "Y"))

SAMPLE_ORDER <- c("ctrl_M1", "ctrl_M2", "ctrl_F1", "ctrl_F2",
                  "mut_M1", "mut_M2", "mut_F1", "mut_F2")

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

parse_cli_args <- function() {
  option_list <- list(
    make_option("--results-dir", type = "character", dest = "results_dir",
                default = NULL,
                help = "Results directory (results/dmrseq)"),
    make_option("--mode", type = "character", default = "local",
                help = "Which results to gather: 'local' or 'blocks' [default %default]"),
    make_option("--split-dir", type = "character", dest = "split_dir",
                default = NULL,
                help = "Phase 1 split directory (for cleanup)"),
    make_option("--clean", type = "logical", default = TRUE,
                help = "Delete split files after gathering [default %default]")
  )

  opt <- parse_args(OptionParser(option_list = option_list))
  for (req in c("results_dir")) {
    if (is.null(opt[[req]])) stop("Missing required argument: --", gsub("_", "-", req))
  }
  if (!opt$mode %in% c("local", "blocks")) {
    stop("--mode must be 'local' or 'blocks', got: ", opt$mode)
  }
  opt
}

# ---------------------------------------------------------------------------
# Write summary
# ---------------------------------------------------------------------------

write_summary <- function(all_dmrs, combined, mode, opt, lambda_dt) {
  summary_file <- file.path(opt$results_dir,
                            paste0("mch_dmr_summary_", mode, ".txt"))

  sink(summary_file)
  cat(sprintf("=== mCH DMR Discovery Summary (%s mode) ===\n\n", mode))

  if (is.null(combined) || length(combined) == 0) {
    cat("Result: No DMRs found on any chromosome.\n\n")
    if (mode == "local") {
      cat("This is a known possibility for non-CG data. CH sites lack the\n")
      cat("spatial correlation of CpG islands/shores that dmrseq was designed\n")
      cat("for. Absence of local DMRs does not necessarily mean absence of\n")
      cat("differential CH methylation — see gene-body results and block mode.\n\n")
    } else {
      cat("No large-scale methylation blocks detected.\n\n")
    }
  } else {
    n_dmrs <- length(combined)
    cat(sprintf("DMRs found: %d\n", n_dmrs))
    cat(sprintf("Chromosomes with DMRs: %s\n",
                paste(names(all_dmrs), collapse = ", ")))

    cat("\nPer-chromosome counts:\n")
    for (chr_name in names(all_dmrs)) {
      cat(sprintf("  %s: %d\n", chr_name, length(all_dmrs[[chr_name]])))
    }

    cat("\nWidth distribution:\n")
    w <- width(combined)
    cat(sprintf("  Min: %d  Median: %d  Max: %d  Mean: %.0f\n",
                min(w), median(w), max(w), mean(w)))

    if ("L" %in% names(mcols(combined))) {
      L <- combined$L
      cat(sprintf("\nSites per DMR (L):\n"))
      cat(sprintf("  Min: %d  Median: %d  Max: %d  Mean: %.0f\n",
                  min(L), median(L), max(L), mean(L)))
    }

    if ("qvalue" %in% names(mcols(combined))) {
      cat(sprintf("\nDMRs with q < 0.05: %d\n", sum(combined$qvalue < 0.05)))
      cat(sprintf("DMRs with q < 0.10: %d\n", sum(combined$qvalue < 0.10)))
    }

    if ("beta" %in% names(mcols(combined))) {
      cat("\nDirection:\n")
      cat(sprintf("  Hyper in mutant: %d\n", sum(combined$beta > 0)))
      cat(sprintf("  Hypo in mutant: %d\n", sum(combined$beta < 0)))
    }
  }

  cat(sprintf("\nMode: %s\n", mode))
  cat(sprintf("Detection type: %s\n",
              if (mode == "blocks") "large-scale methylation blocks"
              else "local differentially methylated regions"))

  if (!is.null(lambda_dt)) {
    cat("\nLambda rates:\n")
    for (i in seq_len(nrow(lambda_dt))) {
      cat(sprintf("  %s: %.4f\n", lambda_dt$sample[i], lambda_dt$lambda_rate[i]))
    }
  }
  sink()

  cat(sprintf("Summary written: %s\n", summary_file))
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  opt <- parse_cli_args()

  cat(sprintf("=== dmrseq Phase 3: Gather results (%s mode) ===\n\n", opt$mode))

  per_chr_dir <- file.path(opt$results_dir, ".per_chr", opt$mode)
  if (!dir.exists(per_chr_dir)) {
    stop("Per-chromosome results directory not found: ", per_chr_dir)
  }

  lambda_dt <- NULL
  lambda_file <- file.path(opt$results_dir, "lambda_rates.tsv")
  if (file.exists(lambda_file)) {
    lambda_dt <- fread(lambda_file)
  }

  all_dmrs <- list()
  n_processed <- 0L
  n_with_dmrs <- 0L

  for (chr_name in CANONICAL_CHRS) {
    rds_file <- file.path(per_chr_dir, paste0(chr_name, "_dmrs.rds"))
    if (!file.exists(rds_file)) {
      cat(sprintf("  %s: RDS file missing (skipping)\n", chr_name))
      next
    }

    chr_dmrs <- readRDS(rds_file)
    n_processed <- n_processed + 1L

    if (!is.null(chr_dmrs) && length(chr_dmrs) > 0) {
      all_dmrs[[chr_name]] <- chr_dmrs
      n_with_dmrs <- n_with_dmrs + 1L
      cat(sprintf("  %s: %d DMRs\n", chr_name, length(chr_dmrs)))
    } else {
      cat(sprintf("  %s: no DMRs\n", chr_name))
    }
  }

  cat(sprintf("\nProcessed: %d chromosomes, %d with DMRs\n",
              n_processed, n_with_dmrs))

  out_tsv <- file.path(opt$results_dir, paste0("mch_dmrs_", opt$mode, ".tsv"))

  if (length(all_dmrs) == 0) {
    cat("No DMRs found on any chromosome.\n")
    fwrite(data.table(note = paste0("No candidate regions passed the cutoff ",
                                    "on any chromosome (", opt$mode, " mode)")),
           out_tsv, sep = "\t")
  } else {
    combined <- do.call(c, unname(all_dmrs))
    n_dmrs <- length(combined)
    cat(sprintf("Total DMRs: %d\n", n_dmrs))

    dmr_df <- as.data.frame(combined)
    fwrite(as.data.table(dmr_df), out_tsv, sep = "\t")
  }

  cat(sprintf("Results written: %s\n", out_tsv))

  combined_for_summary <- if (length(all_dmrs) > 0) {
    do.call(c, unname(all_dmrs))
  } else {
    NULL
  }
  write_summary(all_dmrs, combined_for_summary, opt$mode, opt, lambda_dt)

  if (opt$clean && !is.null(opt$split_dir) && dir.exists(opt$split_dir)) {
    other_mode <- if (opt$mode == "local") "blocks" else "local"
    other_dir <- file.path(opt$results_dir, ".per_chr", other_mode)
    other_gathered <- file.exists(
      file.path(opt$results_dir, paste0("mch_dmrs_", other_mode, ".tsv"))
    )

    if (other_gathered) {
      cat("\nBoth modes gathered. Cleaning up split files...\n")
      unlink(opt$split_dir, recursive = TRUE)
      unlink(file.path(opt$results_dir, ".per_chr"), recursive = TRUE)
      cat("Cleanup complete.\n")
    } else {
      cat(sprintf("\nKeeping split files — %s mode not yet gathered.\n", other_mode))
    }
  }

  cat("\n=== Phase 3 complete ===\n")
}

main()
