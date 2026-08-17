# scripts/dmrseq/02_dmrseq_chr.R
#
# Phase 2 of the dmrseq pipeline: load one chromosome's pre-split RDS files
# (all 8 samples), merge to shared sites, build BSseq, and run dmrseq.
#
# Runs as a SLURM array task — one chromosome per task.  Accepts --mode to
# switch between local DMR discovery (default) and large-scale block detection
# (block = TRUE in dmrseq).
#
# The fork/OpenMP fix: data.table threads are set to 1 at startup, before
# MulticoreParam forks the R process.  This prevents the "wrong args for
# environment subassignment" crash caused by forking with active OpenMP state.
#
# Usage:
#   Rscript 02_dmrseq_chr.R \
#       --chr chr1 \
#       --split-dir /path/to/.dmrseq_splits \
#       --out-dir /path/to/results/dmrseq \
#       --mode local \
#       --workers 28

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(bsseq)
  library(dmrseq)
  library(GenomicRanges)
  library(BiocParallel)
})

setDTthreads(1)

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
    make_option("--chr", type = "character", default = NULL,
                help = "Chromosome to process (e.g. chr1)"),
    make_option("--split-dir", type = "character", dest = "split_dir",
                default = NULL,
                help = "Directory with per-chromosome RDS splits from Phase 1"),
    make_option("--out-dir", type = "character", dest = "out_dir",
                default = NULL,
                help = "Output directory (results/dmrseq)"),
    make_option("--mode", type = "character", default = "local",
                help = "DMR detection mode: 'local' or 'blocks' [default %default]"),
    make_option("--cutoff", type = "double", default = 0.005,
                help = "Methylation difference cutoff [default %default]"),
    make_option("--min-num-region", type = "integer", dest = "min_num_region",
                default = 5L,
                help = "Minimum CH sites per candidate region [default %default]"),
    make_option("--workers", type = "integer", default = 28L,
                help = "BiocParallel workers [default %default]")
  )

  opt <- parse_args(OptionParser(option_list = option_list))
  for (req in c("chr", "split_dir", "out_dir")) {
    if (is.null(opt[[req]])) stop("Missing required argument: --", gsub("_", "-", req))
  }
  if (!opt$mode %in% c("local", "blocks")) {
    stop("--mode must be 'local' or 'blocks', got: ", opt$mode)
  }
  opt
}

# ---------------------------------------------------------------------------
# Process one chromosome
# ---------------------------------------------------------------------------

process_chromosome <- function(chr_name, split_dir, bp_param, cutoff,
                               min_num_region, block_mode) {
  cat(sprintf("\n--- %s (mode: %s) ---\n", chr_name,
              if (block_mode) "blocks" else "local"))

  accumulated <- NULL
  for (i in seq_along(SAMPLE_ORDER)) {
    sid <- SAMPLE_ORDER[i]
    rds_path <- file.path(split_dir, paste0(sid, "_", chr_name, ".rds"))
    if (!file.exists(rds_path)) stop("RDS file not found: ", rds_path)
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

  bs <- BSseq(chr = rep(chr_name, length(pos_vec)), pos = pos_vec,
              M = M_mat, Cov = Cov_mat,
              sampleNames = SAMPLE_ORDER)
  pData(bs)$genotype <- SAMPLE_GENO
  pData(bs)$sex <- SAMPLE_SEX
  bs <- sort(bs)

  rm(M_mat, Cov_mat, pos_vec); gc(verbose = FALSE)
  cat(sprintf("  BSseq: %s sites x %d samples\n",
              format(nrow(bs), big.mark = ","), ncol(bs)))

  chr_dmrs <- tryCatch(
    dmrseq(bs,
           testCovariate = "genotype",
           adjustCovariate = "sex",
           cutoff = cutoff,
           minNumRegion = min_num_region,
           block = block_mode,
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

  cat(sprintf("=== dmrseq Phase 2: %s (mode: %s) ===\n\n", opt$chr, opt$mode))
  cat("Threading:\n")
  cat(sprintf("  data.table threads: %d (fixed at 1 to prevent fork/OpenMP conflict)\n",
              getDTthreads()))
  cat(sprintf("  BiocParallel workers: %d\n", opt$workers))
  cat(sprintf("  cutoff: %s\n", opt$cutoff))
  cat(sprintf("  minNumRegion: %d\n", opt$min_num_region))
  cat(sprintf("  block mode: %s\n", opt$mode == "blocks"))
  cat("")

  bp_param <- MulticoreParam(workers = opt$workers)

  result_dir <- file.path(opt$out_dir, ".per_chr", opt$mode)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  block_mode <- opt$mode == "blocks"
  chr_dmrs <- process_chromosome(opt$chr, opt$split_dir, bp_param,
                                 opt$cutoff, opt$min_num_region, block_mode)

  result_file <- file.path(result_dir, paste0(opt$chr, "_dmrs.rds"))
  saveRDS(chr_dmrs, result_file)
  cat(sprintf("\nResult saved: %s\n", result_file))

  cat("\n=== Phase 2 complete ===\n")
}

main()
