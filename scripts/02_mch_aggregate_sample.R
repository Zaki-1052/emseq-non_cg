# scripts/02_mch_aggregate_sample.R
#
# Per-sample gene-body CH methylation aggregation.
#
# Reads one combined CHH+CHG methylKit file, extracts the lambda spike-in
# conversion noise rate, filters to canonical chromosomes, and aggregates
# CH site counts over protein-coding gene bodies using GenomicRanges.
#
# Usage:
#   Rscript mch_aggregate_sample.R \
#       --sample-id ctrl_M1 \
#       --ch-file /path/to/combined_ch/ctrl_M1_CH.methylKit.gz \
#       --gene-bed /path/to/gene_bodies.protein_coding.bed \
#       --out-dir /path/to/results/02_aggregate/aggregated

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

CANONICAL_CHRS <- paste0("chr", c(1:19, "X", "Y"))

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
    make_option("--gene-bed", type = "character", dest = "gene_bed",
                default = NULL,
                help = "Path to gene body BED (gene_bodies.protein_coding.bed)"),
    make_option("--out-dir", type = "character", dest = "out_dir",
                default = NULL,
                help = "Output directory for aggregated results"),
    make_option("--threads", type = "integer", dest = "threads",
                default = getDTthreads(),
                help = "Threads for data.table [default %default]")
  )

  opt <- parse_args(OptionParser(option_list = option_list))

  for (req in c("sample_id", "ch_file", "gene_bed", "out_dir")) {
    if (is.null(opt[[req]])) {
      stop("Missing required argument: --", gsub("_", "-", req))
    }
  }
  if (!file.exists(opt$ch_file)) stop("CH file not found: ", opt$ch_file)
  if (!file.exists(opt$gene_bed)) stop("Gene BED not found: ", opt$gene_bed)
  opt
}

# ---------------------------------------------------------------------------
# Gene body GRanges
# ---------------------------------------------------------------------------

load_gene_bodies <- function(bed_path) {
  gene_bed <- fread(bed_path, header = FALSE,
                    col.names = c("chr", "start", "end", "gene_name", "gene_id",
                                  "strand", "gene_type", "gene_length"))
  # BED 0-based half-open -> GRanges 1-based closed: start+1, end unchanged
  gene_gr <- GRanges(
    seqnames = gene_bed$chr,
    ranges = IRanges(start = gene_bed$start + 1L, end = gene_bed$end)
  )
  list(bed = gene_bed, gr = gene_gr)
}

# ---------------------------------------------------------------------------
# Lambda spike-in extraction
# ---------------------------------------------------------------------------

extract_lambda <- function(ch, sample_id, out_dir) {
  lambda <- ch[chr == "phage_lambda"]
  if (nrow(lambda) == 0) {
    stop("No phage_lambda sites found — check chromosome naming in the input file")
  }

  lambda_meth <- sum(as.integer(round(lambda$coverage * lambda$freqC / 100)))
  lambda_total <- sum(lambda$coverage)
  lambda_rate <- lambda_meth / lambda_total

  cat(sprintf("  Lambda CH rate: %.4f (%.2f%%)\n", lambda_rate, lambda_rate * 100))
  cat(sprintf("  Lambda sites: %d, total coverage: %d\n",
              nrow(lambda), lambda_total))

  lambda_df <- data.table(
    sample_id = sample_id,
    lambda_ch_rate = lambda_rate,
    lambda_n_sites = nrow(lambda),
    lambda_total_cov = lambda_total,
    lambda_meth_reads = lambda_meth
  )
  lambda_file <- file.path(out_dir, paste0(sample_id, "_lambda.tsv"))
  fwrite(lambda_df, lambda_file, sep = "\t")
  cat("  Saved:", lambda_file, "\n")

  lambda_rate
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  opt <- parse_cli_args()
  setDTthreads(opt$threads)
  cat("Sample:", opt$sample_id, "\n")
  cat("data.table threads:", getDTthreads(), "\n\n")

  genes <- load_gene_bodies(opt$gene_bed)
  cat("Gene bodies loaded:", nrow(genes$bed), "genes\n\n")

  # --- Read combined CH file ---
  cat("Reading:", opt$ch_file, "\n")
  ch <- fread(opt$ch_file, select = c("chr", "base", "coverage", "freqC"))
  cat(sprintf("  Total sites: %s\n", format(nrow(ch), big.mark = ",")))

  # --- Lambda extraction (before chr filtering) ---
  cat("\nLambda spike-in:\n")
  lambda_rate <- extract_lambda(ch, opt$sample_id, opt$out_dir)

  # --- Filter to canonical chromosomes ---
  ch <- ch[chr %chin% CANONICAL_CHRS]
  cat(sprintf("\nCanonical chr sites: %s\n", format(nrow(ch), big.mark = ",")))

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

  # --- Find overlaps with gene bodies ---
  cat("Finding overlaps...\n")
  hits <- findOverlaps(ch_gr, genes$gr)
  cat(sprintf("  Overlapping pairs: %s\n", format(length(hits), big.mark = ",")))
  rm(ch_gr); gc(verbose = FALSE)

  # --- Aggregate per gene ---
  cat("Aggregating per gene...\n")
  agg <- data.table(
    gene_idx = subjectHits(hits),
    coverage = ch_cov[queryHits(hits)],
    meth_reads = ch_meth[queryHits(hits)]
  )
  rm(hits, ch_cov, ch_meth); gc(verbose = FALSE)

  gene_stats <- agg[, .(
    n_ch_sites = .N,
    total_coverage = sum(as.numeric(coverage)),
    meth_reads = sum(as.numeric(meth_reads))
  ), by = gene_idx]
  rm(agg); gc(verbose = FALSE)

  # --- Build output for ALL genes (zeros for genes without CH sites) ---
  result <- data.table(
    gene_name = genes$bed$gene_name,
    gene_id = genes$bed$gene_id,
    chr = genes$bed$chr,
    start = genes$bed$start,
    end = genes$bed$end,
    strand = genes$bed$strand,
    gene_type = genes$bed$gene_type,
    gene_length = genes$bed$gene_length,
    n_ch_sites = 0L,
    total_coverage = 0,
    meth_reads = 0
  )

  idx <- gene_stats$gene_idx
  result[idx, `:=`(
    n_ch_sites = as.integer(gene_stats$n_ch_sites),
    total_coverage = gene_stats$total_coverage,
    meth_reads = gene_stats$meth_reads
  )]

  out_file <- file.path(opt$out_dir, paste0(opt$sample_id, "_genebody_mch.tsv"))
  fwrite(result, out_file, sep = "\t")
  cat("\nWrote", nrow(result), "genes to:", out_file, "\n")

  # --- Summary ---
  has_sites <- result$n_ch_sites > 0
  genome_rate <- sum(result$meth_reads) / sum(result$total_coverage)

  cat("\n--- Summary ---\n")
  cat("Genes with CH sites:", sum(has_sites), "/", nrow(result), "\n")
  cat("Median CH sites/gene:", median(result$n_ch_sites[has_sites]), "\n")
  cat("Median total coverage/gene:", median(result$total_coverage[has_sites]), "\n")
  cat(sprintf("Gene-body mCH rate: %.4f (%.2f%%)\n", genome_rate, genome_rate * 100))
  cat(sprintf("Lambda CH rate: %.4f (%.2f%%)\n", lambda_rate, lambda_rate * 100))
}

main()
