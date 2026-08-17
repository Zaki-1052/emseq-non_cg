# scripts/dmrseq/01_split_samples.R
#
# Phase 1 of the dmrseq pipeline: read each sample's combined CH methylKit
# file, extract the lambda spike-in conversion rate, filter to canonical
# chromosomes, apply lambda correction to methylation counts, split by
# chromosome, and save per-chromosome RDS files to disk.
#
# Produces 168 RDS files (8 samples x 21 chromosomes) plus a lambda_rates.tsv
# metadata file. Phase 2 (02_dmrseq_chr.R) loads one chromosome at a time.
#
# Usage:
#   Rscript 01_split_samples.R \
#       --ch-dir /path/to/combined_ch \
#       --out-dir /path/to/results/dmrseq \
#       --split-dir /path/to/.dmrseq_splits \
#       --threads 16

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

CANONICAL_CHRS <- paste0("chr", c(1:19, "X", "Y"))

SAMPLE_ORDER <- c("ctrl_M1", "ctrl_M2", "ctrl_F1", "ctrl_F2",
                  "mut_M1", "mut_M2", "mut_F1", "mut_F2")

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

parse_cli_args <- function() {
  option_list <- list(
    make_option("--ch-dir", type = "character", dest = "ch_dir",
                default = NULL,
                help = "Directory with combined CH methylKit files"),
    make_option("--out-dir", type = "character", dest = "out_dir",
                default = NULL,
                help = "Output directory for final results (results/dmrseq)"),
    make_option("--split-dir", type = "character", dest = "split_dir",
                default = NULL,
                help = "Directory for per-chromosome RDS splits"),
    make_option("--threads", type = "integer", default = 16L,
                help = "data.table threads [default %default]")
  )

  opt <- parse_args(OptionParser(option_list = option_list))
  for (req in c("ch_dir", "out_dir", "split_dir")) {
    if (is.null(opt[[req]])) stop("Missing required argument: --", gsub("_", "-", req))
  }
  opt
}

# ---------------------------------------------------------------------------
# Split one sample to per-chromosome RDS files
# ---------------------------------------------------------------------------

split_sample_to_disk <- function(ch_file, sample_id, split_dir) {
  cat(sprintf("  Reading %s: %s\n", sample_id, ch_file))
  ch <- fread(ch_file, select = c("chr", "base", "coverage", "freqC"))
  cat(sprintf("    Total sites: %s\n", format(nrow(ch), big.mark = ",")))

  lambda <- ch[chr == "phage_lambda"]
  if (nrow(lambda) == 0) stop("No phage_lambda sites in ", ch_file)
  lambda_meth <- sum(as.integer(round(lambda$coverage * lambda$freqC / 100)))
  lambda_total <- sum(lambda$coverage)
  lambda_rate <- lambda_meth / lambda_total
  cat(sprintf("    Lambda rate: %.4f (%.2f%%)\n", lambda_rate, lambda_rate * 100))

  ch <- ch[chr %chin% CANONICAL_CHRS]
  cat(sprintf("    Canonical sites: %s\n", format(nrow(ch), big.mark = ",")))

  ch[, M := as.integer(round(coverage * freqC / 100))]
  ch[, M := pmax(0L, as.integer(round(M - lambda_rate * coverage)))]

  n_saved <- 0L
  for (chr_name in CANONICAL_CHRS) {
    chr_dt <- ch[chr == chr_name, .(base, M, Cov = coverage)]
    setkey(chr_dt, base)
    saveRDS(chr_dt, file.path(split_dir, paste0(sample_id, "_", chr_name, ".rds")))
    n_saved <- n_saved + nrow(chr_dt)
  }

  rm(ch, lambda); gc(verbose = FALSE)
  cat(sprintf("    Saved %s sites across %d chromosomes\n",
              format(n_saved, big.mark = ","), length(CANONICAL_CHRS)))

  lambda_rate
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  opt <- parse_cli_args()
  setDTthreads(opt$threads)

  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(opt$split_dir, recursive = TRUE, showWarnings = FALSE)

  cat("=== dmrseq Phase 1: Split samples by chromosome ===\n\n")
  cat("Parameters:\n")
  cat("  data.table threads:", getDTthreads(), "\n")
  cat("  CH dir:    ", opt$ch_dir, "\n")
  cat("  Split dir: ", opt$split_dir, "\n")
  cat("  Out dir:   ", opt$out_dir, "\n\n")

  lambda_rates <- numeric(length(SAMPLE_ORDER))
  names(lambda_rates) <- SAMPLE_ORDER

  for (i in seq_along(SAMPLE_ORDER)) {
    sid <- SAMPLE_ORDER[i]
    ch_file <- file.path(opt$ch_dir, paste0(sid, "_CH.methylKit.gz"))
    if (!file.exists(ch_file)) stop("File not found: ", ch_file)
    lambda_rates[i] <- split_sample_to_disk(ch_file, sid, opt$split_dir)
  }

  cat(sprintf("\nPhase 1 complete. %d samples split across %d chromosomes.\n",
              length(SAMPLE_ORDER), length(CANONICAL_CHRS)))
  cat("Lambda rates:\n")
  for (i in seq_along(SAMPLE_ORDER)) {
    cat(sprintf("  %s: %.4f (%.2f%%)\n",
                SAMPLE_ORDER[i], lambda_rates[i], lambda_rates[i] * 100))
  }

  lambda_dt <- data.table(
    sample = SAMPLE_ORDER,
    lambda_rate = lambda_rates
  )
  lambda_file <- file.path(opt$out_dir, "lambda_rates.tsv")
  fwrite(lambda_dt, lambda_file, sep = "\t")
  cat(sprintf("\nLambda rates saved to: %s\n", lambda_file))

  n_rds <- length(list.files(opt$split_dir, pattern = "\\.rds$"))
  cat(sprintf("RDS files created: %d\n", n_rds))

  cat("\n=== Phase 1 complete ===\n")
}

main()
