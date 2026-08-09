# scripts/06_mch_dmrseq.R
#
# Genome-wide CH DMR discovery with dmrseq.
#
# Processes chromosomes sequentially to bound memory usage. Phase 1 reads
# each sample once and saves per-chromosome splits as RDS to disk. Phase 2
# loads one chromosome at a time across all samples, merges to shared sites,
# builds BSseq, and runs dmrseq. Peak memory is bounded by the largest
# single chromosome (~5-10% of total data), not the full genome.
#
# Usage:
#   Rscript mch_dmrseq.R \
#       --ch-dir /path/to/combined_ch \
#       --out-dir /path/to/results/dmrseq

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(bsseq)
  library(dmrseq)
  library(GenomicRanges)
  library(BiocParallel)
})

CANONICAL_CHRS <- paste0("chr", c(1:19, "X", "Y"))

SAMPLE_ORDER <- c("ctrl_M1", "ctrl_M2", "ctrl_F1", "ctrl_F2",
                  "mut_M1", "mut_M2", "mut_F1", "mut_F2")

SAMPLE_GENO <- factor(c(rep("ctrl", 4), rep("mut", 4)),
                      levels = c("ctrl", "mut"))
SAMPLE_SEX <- factor(c("M", "M", "F", "F", "M", "M", "F", "F"))

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
                help = "Output directory"),
    make_option("--cutoff", type = "double", default = 0.005,
                help = "Methylation difference cutoff [default %default]"),
    make_option("--min-num-region", type = "integer", dest = "min_num_region",
                default = 5L,
                help = "Minimum CH sites per candidate region [default %default]"),
    make_option("--workers", type = "integer", default = 40L,
                help = "BiocParallel workers [default %default]"),
    make_option("--threads", type = "integer", default = getDTthreads(),
                help = "data.table threads [default %default]")
  )

  opt <- parse_args(OptionParser(option_list = option_list))
  for (req in c("ch_dir", "out_dir")) {
    if (is.null(opt[[req]])) stop("Missing required argument: --", gsub("_", "-", req))
  }
  opt
}

# ---------------------------------------------------------------------------
# Phase 1: Read each sample, extract lambda, split by chromosome, save to disk
# ---------------------------------------------------------------------------

split_sample_to_disk <- function(ch_file, sample_id, split_dir) {
  cat(sprintf("  Reading %s: %s\n", sample_id, ch_file))
  ch <- fread(ch_file, select = c("chr", "base", "coverage", "freqC"))
  cat(sprintf("    Total sites: %s\n", format(nrow(ch), big.mark = ",")))

  # Lambda extraction
  lambda <- ch[chr == "phage_lambda"]
  if (nrow(lambda) == 0) stop("No phage_lambda sites in ", ch_file)
  lambda_meth <- sum(as.integer(round(lambda$coverage * lambda$freqC / 100)))
  lambda_total <- sum(lambda$coverage)
  lambda_rate <- lambda_meth / lambda_total
  cat(sprintf("    Lambda rate: %.4f (%.2f%%)\n", lambda_rate, lambda_rate * 100))

  # Filter to canonical chromosomes
  ch <- ch[chr %chin% CANONICAL_CHRS]
  cat(sprintf("    Canonical sites: %s\n", format(nrow(ch), big.mark = ",")))

  # Compute M with lambda correction
  ch[, M := as.integer(round(coverage * freqC / 100))]
  ch[, M := pmax(0L, as.integer(round(M - lambda_rate * coverage)))]

  # Split by chromosome and save
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
# Phase 2: Process one chromosome — merge samples, build BSseq, run dmrseq
# ---------------------------------------------------------------------------

process_chromosome <- function(chr_name, split_dir, bp_param, cutoff, min_num_region) {
  cat(sprintf("\n--- %s ---\n", chr_name))

  # Load and sequentially merge all 8 samples
  accumulated <- NULL
  for (i in seq_along(SAMPLE_ORDER)) {
    sid <- SAMPLE_ORDER[i]
    rds_path <- file.path(split_dir, paste0(sid, "_", chr_name, ".rds"))
    dt <- readRDS(rds_path)

    m_col <- paste0("M_", i)
    cov_col <- paste0("Cov_", i)
    setnames(dt, c("M", "Cov"), c(m_col, cov_col))

    if (is.null(accumulated)) {
      accumulated <- dt
    } else {
      accumulated <- merge(accumulated, dt, by = "base")
    }
    rm(dt)
  }

  n_sites <- nrow(accumulated)
  cat(sprintf("  Shared sites: %s\n", format(n_sites, big.mark = ",")))

  if (n_sites < min_num_region) {
    cat("  Skipping — too few shared sites\n")
    rm(accumulated); gc(verbose = FALSE)
    return(NULL)
  }

  # Extract matrices
  m_cols <- paste0("M_", seq_along(SAMPLE_ORDER))
  cov_cols <- paste0("Cov_", seq_along(SAMPLE_ORDER))

  M_mat <- as.matrix(accumulated[, ..m_cols])
  Cov_mat <- as.matrix(accumulated[, ..cov_cols])
  colnames(M_mat) <- SAMPLE_ORDER
  colnames(Cov_mat) <- SAMPLE_ORDER
  storage.mode(M_mat) <- "integer"
  storage.mode(Cov_mat) <- "integer"

  pos_vec <- accumulated$base
  rm(accumulated); gc(verbose = FALSE)

  # Filter sites with zero coverage in any sample
  min_cov <- apply(Cov_mat, 1, min)
  keep <- min_cov >= 1L
  n_keep <- sum(keep)

  if (n_keep < min_num_region) {
    cat(sprintf("  Skipping — only %d sites with coverage in all samples\n", n_keep))
    rm(M_mat, Cov_mat, pos_vec); gc(verbose = FALSE)
    return(NULL)
  }

  if (n_keep < length(keep)) {
    M_mat <- M_mat[keep, , drop = FALSE]
    Cov_mat <- Cov_mat[keep, , drop = FALSE]
    pos_vec <- pos_vec[keep]
    cat(sprintf("  After coverage filter: %s\n", format(n_keep, big.mark = ",")))
  }

  # Build BSseq
  bs <- BSseq(chr = rep(chr_name, length(pos_vec)), pos = pos_vec,
              M = M_mat, Cov = Cov_mat,
              sampleNames = SAMPLE_ORDER)
  pData(bs)$genotype <- SAMPLE_GENO
  pData(bs)$sex <- SAMPLE_SEX
  bs <- sort(bs)

  rm(M_mat, Cov_mat, pos_vec); gc(verbose = FALSE)
  cat(sprintf("  BSseq: %s sites x %d samples\n",
              format(nrow(bs), big.mark = ","), ncol(bs)))

  # Run dmrseq
  chr_dmrs <- tryCatch(
    dmrseq(bs,
           testCovariate = "genotype",
           adjustCovariate = "sex",
           cutoff = cutoff,
           minNumRegion = min_num_region,
           BPPARAM = bp_param),
    error = function(e) {
      cat(sprintf("  dmrseq error: %s\n", conditionMessage(e)))
      NULL
    }
  )

  rm(bs); gc(verbose = FALSE)

  if (!is.null(chr_dmrs) && length(chr_dmrs) > 0) {
    cat(sprintf("  Found %d DMRs\n", length(chr_dmrs)))
  } else {
    cat("  No DMRs\n")
    chr_dmrs <- NULL
  }

  chr_dmrs
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  opt <- parse_cli_args()
  setDTthreads(opt$threads)
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

  split_dir <- file.path(opt$out_dir, ".chr_splits")
  dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)

  cat("=== Genome-wide mCH DMR Discovery (dmrseq, per-chromosome) ===\n\n")
  cat("Parameters:\n")
  cat("  cutoff:", opt$cutoff, "\n")
  cat("  minNumRegion:", opt$min_num_region, "\n")
  cat("  workers:", opt$workers, "\n")
  cat("  data.table threads:", getDTthreads(), "\n\n")

  # --- Phase 1: Read samples and split by chromosome ---
  cat("Phase 1: Reading samples and splitting by chromosome\n")
  lambda_rates <- numeric(length(SAMPLE_ORDER))
  names(lambda_rates) <- SAMPLE_ORDER

  for (i in seq_along(SAMPLE_ORDER)) {
    sid <- SAMPLE_ORDER[i]
    ch_file <- file.path(opt$ch_dir, paste0(sid, "_CH.methylKit.gz"))
    if (!file.exists(ch_file)) stop("File not found: ", ch_file)
    lambda_rates[i] <- split_sample_to_disk(ch_file, sid, split_dir)
  }

  cat(sprintf("\nPhase 1 complete. %d samples split across %d chromosomes.\n",
              length(SAMPLE_ORDER), length(CANONICAL_CHRS)))
  cat("Lambda rates:\n")
  for (i in seq_along(SAMPLE_ORDER)) {
    cat(sprintf("  %s: %.4f (%.2f%%)\n",
                SAMPLE_ORDER[i], lambda_rates[i], lambda_rates[i] * 100))
  }

  # --- Phase 2: Process each chromosome ---
  cat("\nPhase 2: Running dmrseq per chromosome\n")
  bp_param <- MulticoreParam(workers = opt$workers)
  all_dmrs <- list()
  total_shared <- 0L

  for (chr_name in CANONICAL_CHRS) {
    chr_dmrs <- process_chromosome(chr_name, split_dir, bp_param,
                                   opt$cutoff, opt$min_num_region)
    if (!is.null(chr_dmrs)) {
      all_dmrs[[chr_name]] <- chr_dmrs
    }
    gc(verbose = FALSE)
  }

  # --- Phase 3: Combine and write results ---
  cat("\n\nPhase 3: Writing results\n")

  if (length(all_dmrs) == 0) {
    cat("No DMRs found on any chromosome.\n")

    fwrite(data.table(note = "No candidate regions passed the cutoff on any chromosome"),
           file.path(opt$out_dir, "mch_dmrs.tsv"), sep = "\t")

    sink(file.path(opt$out_dir, "mch_dmr_summary.txt"))
    cat("=== mCH DMR Discovery Summary ===\n\n")
    cat("Result: No DMRs found on any chromosome.\n\n")
    cat("This is a known possibility for non-CG data. CH sites lack the\n")
    cat("spatial correlation of CpG islands/shores that dmrseq was designed\n")
    cat("for. Absence of DMRs does not necessarily mean absence of\n")
    cat("differential CH methylation — see gene-body results (Step 2b).\n\n")
    cat("Parameters:\n")
    cat(sprintf("  cutoff: %.3f\n", opt$cutoff))
    cat(sprintf("  minNumRegion: %d\n", opt$min_num_region))
    cat("\nLambda rates:\n")
    for (i in seq_along(SAMPLE_ORDER)) {
      cat(sprintf("  %s: %.4f\n", SAMPLE_ORDER[i], lambda_rates[i]))
    }
    sink()
  } else {
    combined <- do.call(c, unname(all_dmrs))
    n_dmrs <- length(combined)
    cat(sprintf("Total DMRs across all chromosomes: %d\n", n_dmrs))

    dmr_df <- as.data.frame(combined)
    fwrite(as.data.table(dmr_df),
           file.path(opt$out_dir, "mch_dmrs.tsv"), sep = "\t")

    sink(file.path(opt$out_dir, "mch_dmr_summary.txt"))
    cat("=== mCH DMR Discovery Summary ===\n\n")
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
    cat("\nParameters:\n")
    cat(sprintf("  cutoff: %.3f\n", opt$cutoff))
    cat(sprintf("  minNumRegion: %d\n", opt$min_num_region))
    cat("\nLambda rates:\n")
    for (i in seq_along(SAMPLE_ORDER)) {
      cat(sprintf("  %s: %.4f\n", SAMPLE_ORDER[i], lambda_rates[i]))
    }
    sink()
  }

  # Clean up split files
  cat("Cleaning up chromosome splits...\n")
  unlink(split_dir, recursive = TRUE)

  cat("\n=== Complete ===\n")
}

main()
