# scripts/02b_mch_aggregate_features.R
#
# Per-sample sub-gene feature CH methylation aggregation.
#
# Reads one combined CHH+CHG methylKit file and aggregates CH site counts over
# sub-gene feature intervals (5' UTR, exon, splice donor, splice acceptor,
# intron, 3' UTR) instead of whole gene bodies. Section 50_01 compares
# methylation between feature types from these tables.
#
# Runs in parallel with step 02, not after it. Both read the step 01 output.
#
# The feature BEDs come from scripts/utils/generate_feature_beds.R. Each is
# BED6 with the name column "<gene_symbol>|<rank>". Feature type comes from
# the file name, so it is recorded per interval when the BED is read.
#
# Usage:
#   Rscript scripts/02b_mch_aggregate_features.R \
#       --sample-id ctrl_M1 \
#       --ch-file /path/to/combined_ch/ctrl_M1_CH.methylKit.gz \
#       --feature-dir /path/to/data/features \
#       --out-dir /path/to/results/02b_features

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

CANONICAL_CHRS <- paste0("chr", c(1:19, "X", "Y"))

# Feature type to BED file name. The order sets the factor level order used by
# section 50_01, running 5' to 3' along a transcript.
FEATURE_FILES <- c(
  "5UTR"                = "utr5_protein_coding.bed",
  "Exon"                = "exons_protein_coding.bed",
  "SpliceSite_Donor"    = "splice_donor_protein_coding.bed",
  "SpliceSite_Acceptor" = "splice_acceptor_protein_coding.bed",
  "Intron"              = "introns_protein_coding.bed",
  "3UTR"                = "utr3_protein_coding.bed"
)

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

parse_cli_args <- function() {
  option_list <- list(
    make_option("--sample-id", type = "character", dest = "sample_id",
                default = NULL,
                help = "Sample identifier (e.g. ctrl_M1)"),
    make_option("--ch-file", type = "character", dest = "ch_file",
                default = NULL,
                help = "Path to combined CH methylKit file (.gz)"),
    make_option("--feature-dir", type = "character", dest = "feature_dir",
                default = NULL,
                help = "Directory holding the feature BED files"),
    make_option("--out-dir", type = "character", dest = "out_dir",
                default = NULL,
                help = "Output directory for aggregated feature results"),
    make_option("--threads", type = "integer", dest = "threads",
                default = getDTthreads(),
                help = "Threads for data.table [default %default]")
  )

  opt <- parse_args(OptionParser(option_list = option_list))

  for (req in c("sample_id", "ch_file", "feature_dir", "out_dir")) {
    if (is.null(opt[[req]])) {
      stop("Missing required argument: --", gsub("_", "-", req))
    }
  }
  if (!file.exists(opt$ch_file)) stop("CH file not found: ", opt$ch_file)
  if (!dir.exists(opt$feature_dir)) {
    stop("Feature directory not found: ", opt$feature_dir)
  }

  for (fname in FEATURE_FILES) {
    fpath <- file.path(opt$feature_dir, fname)
    if (!file.exists(fpath)) {
      stop("Feature BED not found: ", fpath,
           "\nRun scripts/utils/generate_feature_beds.R first.")
    }
  }
  opt
}

# ---------------------------------------------------------------------------
# Feature intervals
# ---------------------------------------------------------------------------

#' Read every feature BED into one interval table.
#'
#' @param feature_dir Directory holding the BED files.
#' @return list with a data.table of interval metadata and a matching GRanges
load_feature_intervals <- function(feature_dir) {
  parts <- lapply(names(FEATURE_FILES), function(ftype) {
    fpath <- file.path(feature_dir, FEATURE_FILES[[ftype]])
    bed <- fread(fpath, header = FALSE, sep = "\t",
                 col.names = c("chr", "start", "end", "name", "score", "strand"))

    name_parts <- tstrsplit(bed$name, "|", fixed = TRUE)
    if (length(name_parts) != 2) {
      stop("Malformed name column in ", fpath,
           ". Expected '<gene_symbol>|<rank>'.")
    }

    dt <- data.table(
      gene_name = name_parts[[1]],
      feature_type = ftype,
      feature_rank = as.integer(name_parts[[2]]),
      chr = bed$chr,
      start = bed$start,
      end = bed$end,
      strand = bed$strand
    )
    cat(sprintf("  %-20s %8s intervals\n", ftype,
                format(nrow(dt), big.mark = ",")))
    dt
  })

  meta <- rbindlist(parts)
  meta[, feature_type := factor(feature_type, levels = names(FEATURE_FILES))]

  gr <- GRanges(
    seqnames = meta$chr,
    ranges = IRanges(start = meta$start + 1L, end = meta$end)
  )

  list(meta = meta, gr = gr)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  opt <- parse_cli_args()
  setDTthreads(opt$threads)
  cat("Sample:", opt$sample_id, "\n")
  cat("data.table threads:", getDTthreads(), "\n\n")

  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("Loading feature intervals from:", opt$feature_dir, "\n")
  features <- load_feature_intervals(opt$feature_dir)
  cat(sprintf("  Total: %s intervals across %s genes\n\n",
              format(nrow(features$meta), big.mark = ","),
              format(length(unique(features$meta$gene_name)), big.mark = ",")))

  # --- Read combined CH file ---
  cat("Reading:", opt$ch_file, "\n")
  ch <- fread(opt$ch_file, select = c("chr", "base", "coverage", "freqC"))
  cat(sprintf("  Total sites: %s\n", format(nrow(ch), big.mark = ",")))

  # --- Filter to canonical chromosomes ---
  ch <- ch[chr %chin% CANONICAL_CHRS]
  cat(sprintf("  Canonical chr sites: %s\n", format(nrow(ch), big.mark = ",")))

  # --- Compute methylated reads ---
  ch[, meth_reads := as.integer(round(coverage * freqC / 100))]

  # --- Convert to GRanges ---
  cat("Building GRanges...\n")
  ch_gr <- GRanges(
    seqnames = ch$chr,
    ranges = IRanges(start = ch$base, width = 1L)
  )
  ch_cov <- ch$coverage
  ch_meth <- ch$meth_reads
  rm(ch); gc(verbose = FALSE)

  # --- Find overlaps with feature intervals ---
  # A CH site can fall in intervals of more than one feature type, for example
  # a splice site window inside an exon. Each feature type counts it once.
  cat("Finding overlaps...\n")
  hits <- findOverlaps(ch_gr, features$gr)
  cat(sprintf("  Overlapping pairs: %s\n", format(length(hits), big.mark = ",")))
  rm(ch_gr); gc(verbose = FALSE)

  # --- Aggregate per interval ---
  cat("Aggregating per interval...\n")
  agg <- data.table(
    interval_idx = subjectHits(hits),
    coverage = ch_cov[queryHits(hits)],
    meth_reads = ch_meth[queryHits(hits)]
  )
  rm(hits, ch_cov, ch_meth); gc(verbose = FALSE)

  interval_stats <- agg[, .(
    n_ch_sites = .N,
    total_coverage = sum(as.numeric(coverage)),
    methylated_count = sum(as.numeric(meth_reads))
  ), by = interval_idx]
  rm(agg); gc(verbose = FALSE)

  # --- Build output for ALL intervals (zeros where no CH sites fall) ---
  result <- data.table(
    gene_name = features$meta$gene_name,
    feature_type = as.character(features$meta$feature_type),
    feature_rank = features$meta$feature_rank,
    chr = features$meta$chr,
    start = features$meta$start,
    end = features$meta$end,
    strand = features$meta$strand,
    feature_length = features$meta$end - features$meta$start,
    n_ch_sites = 0L,
    total_coverage = 0,
    methylated_count = 0
  )

  idx <- interval_stats$interval_idx
  result[idx, `:=`(
    n_ch_sites = as.integer(interval_stats$n_ch_sites),
    total_coverage = interval_stats$total_coverage,
    methylated_count = interval_stats$methylated_count
  )]

  result[, mch_rate := ifelse(total_coverage > 0,
                              methylated_count / total_coverage,
                              NA_real_)]
  result[, sample_id := opt$sample_id]

  setcolorder(result, c("sample_id", "gene_name", "feature_type", "feature_rank",
                        "chr", "start", "end", "strand", "feature_length",
                        "n_ch_sites", "total_coverage", "methylated_count",
                        "mch_rate"))

  out_file <- file.path(opt$out_dir, paste0(opt$sample_id, "_feature_mch.tsv"))
  fwrite(result, out_file, sep = "\t")
  cat("\nWrote", format(nrow(result), big.mark = ","),
      "intervals to:", out_file, "\n")

  # --- Summary per feature type ---
  cat("\n--- Summary by feature type ---\n")
  summary_dt <- result[, .(
    n_intervals = .N,
    n_with_sites = sum(n_ch_sites > 0),
    total_ch_sites = sum(n_ch_sites),
    mch_rate = sum(methylated_count) / sum(total_coverage)
  ), by = feature_type]
  summary_dt[, feature_type := factor(feature_type, levels = names(FEATURE_FILES))]
  setorder(summary_dt, feature_type)

  for (i in seq_len(nrow(summary_dt))) {
    cat(sprintf("  %-20s %8s intervals  %8s with sites  %12s CH sites  mCH %.4f (%.2f%%)\n",
                summary_dt$feature_type[i],
                format(summary_dt$n_intervals[i], big.mark = ","),
                format(summary_dt$n_with_sites[i], big.mark = ","),
                format(summary_dt$total_ch_sites[i], big.mark = ","),
                summary_dt$mch_rate[i],
                summary_dt$mch_rate[i] * 100))
  }
}

main()
