# scripts/utils/generate_feature_beds.R
#
# Build sub-gene feature interval BEDs from a GENCODE GTF.
#
# Step 02b aggregates mCH over these intervals, and section 50_01 compares
# methylation between feature types. Run once during setup.
#
# Six feature types are written, one BED per type:
#   utr5_protein_coding.bed              5' UTR
#   exons_protein_coding.bed             exons
#   splice_donor_protein_coding.bed      donor windows at exon 3' ends
#   splice_acceptor_protein_coding.bed   acceptor windows at exon 5' ends
#   introns_protein_coding.bed           gene body minus exons
#   utr3_protein_coding.bed              3' UTR
#
# Each file is BED6. The name column is "<gene_symbol>|<rank>", where rank
# orders the intervals of one gene along the direction of transcription.
# Feature type comes from the file, so the name column stays unambiguous even
# for gene symbols that contain underscores.
#
# Intervals are reduced per gene, so overlapping transcript isoforms collapse
# into one set of non-overlapping intervals. A base counted once for a gene
# stays counted once.
#
# Usage:
#   Rscript scripts/utils/generate_feature_beds.R \
#       --gtf /path/to/gencode.vM25.annotation.gtf.gz \
#       --out-dir data/features \
#       --splice-window 10

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(data.table)
})

CANONICAL_CHRS <- paste0("chr", c(1:19, "X", "Y"))

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

parse_cli_args <- function() {
  option_list <- list(
    make_option("--gtf", type = "character", dest = "gtf", default = NULL,
                help = "Path to the GENCODE GTF (plain or gzipped)"),
    make_option("--out-dir", type = "character", dest = "out_dir", default = NULL,
                help = "Output directory for the feature BED files"),
    make_option("--splice-window", type = "integer", dest = "splice_window",
                default = 10L,
                help = "Half-width in bp for splice site windows [default %default]")
  )

  opt <- parse_args(OptionParser(option_list = option_list))

  for (req in c("gtf", "out_dir")) {
    if (is.null(opt[[req]])) {
      stop("Missing required argument: --", gsub("_", "-", req))
    }
  }
  if (!file.exists(opt$gtf)) stop("GTF not found: ", opt$gtf)
  if (opt$splice_window < 1L) stop("--splice-window must be at least 1")
  opt
}

# ---------------------------------------------------------------------------
# Protein-coding gene and transcript selection
# ---------------------------------------------------------------------------

#' Read the GTF and keep protein-coding records on canonical chromosomes.
load_gtf_protein_coding <- function(gtf_path) {
  cat("Reading GTF:", gtf_path, "\n")
  gtf <- rtracklayer::import(gtf_path)
  cat(sprintf("  %s records\n", format(length(gtf), big.mark = ",")))

  required <- c("type", "gene_type", "gene_name")
  missing <- setdiff(required, colnames(mcols(gtf)))
  if (length(missing) > 0) {
    stop("GTF is missing attributes: ", paste(missing, collapse = ", "))
  }

  gtf <- gtf[as.character(seqnames(gtf)) %in% CANONICAL_CHRS]
  gtf <- gtf[mcols(gtf)$gene_type == "protein_coding"]
  cat(sprintf("  %s protein-coding records on canonical chromosomes\n",
              format(length(gtf), big.mark = ",")))

  if (length(gtf) == 0) {
    stop("No protein-coding records survived filtering. Check the GTF build ",
         "and chromosome naming (expected UCSC style, for example chr1).")
  }
  gtf
}

#' Split combined UTR records into 5' and 3' UTRs by per-transcript CDS bounds.
#'
#' Older GENCODE releases emit a single "UTR" type. A UTR record is 5' when it
#' lies before the coding sequence in the direction of transcription, and 3'
#' when it lies after it. The comparison uses the coding bounds of the
#' transcript that owns the record. Isoforms of one gene can have different
#' coding bounds, so a gene-level CDS range would call a record 3' whenever
#' another isoform extends its coding region past that record.
#'
#' Records are returned unreduced, so that reduce_by_gene() runs after the
#' classification and the transcript identity is still available here.
#'
#' @param gtf GRanges of protein-coding GTF records with type and transcript_id.
#' @return list of two GRanges named five and three
split_utr_by_transcript <- function(gtf) {
  if (!"transcript_id" %in% colnames(mcols(gtf))) {
    stop("GTF uses a single UTR type but has no transcript_id attribute. ",
         "5' and 3' UTRs cannot be separated without it, because isoforms of ",
         "one gene can have different coding bounds.")
  }

  cds <- gtf[mcols(gtf)$type == "CDS"]
  if (length(cds) == 0) {
    stop("GTF has a single UTR type but no CDS records, so 5' and 3' UTRs ",
         "cannot be separated.")
  }
  utr <- gtf[mcols(gtf)$type == "UTR"]

  cds_tx <- as.character(mcols(cds)$transcript_id)
  if (anyNA(cds_tx)) {
    stop(sum(is.na(cds_tx)), " CDS records have no transcript_id.")
  }

  # Coding bounds of each transcript, in genomic coordinates.
  cds_bounds <- data.table(transcript_id = cds_tx,
                           cds_start = start(cds),
                           cds_end = end(cds))[
    , .(cds_start = min(cds_start), cds_end = max(cds_end)), by = transcript_id]

  utr_tx <- as.character(mcols(utr)$transcript_id)
  idx <- match(utr_tx, cds_bounds$transcript_id)
  if (anyNA(idx)) {
    orphans <- unique(utr_tx[is.na(idx)])
    stop(length(orphans), " transcripts carry UTR records but no CDS records, ",
         "so their UTRs cannot be classified. First few: ",
         paste(head(orphans, 5), collapse = ", "))
  }

  utr_strand <- as.character(strand(utr))
  if (any(utr_strand == "*")) {
    stop(sum(utr_strand == "*"), " UTR records have no strand.")
  }

  tx_cds_start <- cds_bounds$cds_start[idx]
  tx_cds_end   <- cds_bounds$cds_end[idx]
  minus <- utr_strand == "-"

  is_five  <- ifelse(minus, start(utr) > tx_cds_end, end(utr) < tx_cds_start)
  is_three <- ifelse(minus, end(utr) < tx_cds_start, start(utr) > tx_cds_end)

  # Every UTR record of a transcript lies wholly outside the coding bounds of
  # that same transcript. A record that is neither 5' nor 3' overlaps them.
  unclassified <- which(!is_five & !is_three)
  if (length(unclassified) > 0) {
    first <- unclassified[1]
    stop(length(unclassified), " UTR records overlap the coding bounds of ",
         "their own transcript. First: ", as.character(seqnames(utr))[first],
         ":", start(utr)[first], "-", end(utr)[first],
         " (", utr_tx[first], ")")
  }

  cat(sprintf("  %s UTR records: %s 5' and %s 3' across %s transcripts\n",
              format(length(utr), big.mark = ","),
              format(sum(is_five), big.mark = ","),
              format(sum(is_three), big.mark = ","),
              format(length(unique(utr_tx)), big.mark = ",")))

  list(five = utr[is_five], three = utr[is_three])
}

#' Reduce a GRanges to non-overlapping intervals within each gene.
#'
#' @param gr GRanges carrying a gene_name metadata column.
#' @return GRanges with gene_name, sorted by gene and coordinate
reduce_by_gene <- function(gr) {
  if (length(gr) == 0) return(gr)
  split_gr <- split(gr, mcols(gr)$gene_name)
  reduced <- GenomicRanges::reduce(split_gr)
  out <- unlist(reduced, use.names = TRUE)
  mcols(out)$gene_name <- names(out)
  names(out) <- NULL
  sort(out)
}

#' Rank the intervals of each gene along the direction of transcription.
#'
#' Rank 1 is the most 5' interval for the gene's strand.
#'
#' @param gr GRanges with a gene_name column.
#' @param gene_strand Named character vector mapping gene symbol to strand.
#' @return integer vector of ranks, parallel to gr
rank_within_gene <- function(gr, gene_strand) {
  dt <- data.table(
    idx = seq_along(gr),
    gene_name = mcols(gr)$gene_name,
    start = start(gr)
  )
  dt[, strand := gene_strand[gene_name]]
  if (anyNA(dt$strand)) {
    stop("Missing strand for genes: ",
         paste(head(unique(dt$gene_name[is.na(dt$strand)]), 5), collapse = ", "))
  }

  # Plus-strand genes rank by ascending start; minus-strand genes by descending.
  dt[, sort_key := ifelse(strand == "-", -start, start)]
  dt[, rank := frank(sort_key, ties.method = "first"), by = gene_name]
  setorder(dt, idx)
  as.integer(dt$rank)
}

# ---------------------------------------------------------------------------
# Splice site windows
# ---------------------------------------------------------------------------

#' Build splice donor and acceptor windows from exon intervals.
#'
#' The donor site sits at the 3' end of an exon and the acceptor at the 5' end,
#' both relative to the direction of transcription. Terminal boundaries of the
#' gene are dropped, because they are transcript ends rather than splice sites.
#'
#' @param exons GRanges of per-gene reduced exons with gene_name.
#' @param gene_strand Named character vector mapping gene symbol to strand.
#' @param half_width Half-width of each window in bp.
#' @return list of two GRanges named donor and acceptor
build_splice_windows <- function(exons, gene_strand, half_width) {
  dt <- data.table(
    gene_name = mcols(exons)$gene_name,
    chr = as.character(seqnames(exons)),
    start = start(exons),
    end = end(exons)
  )
  dt[, strand := gene_strand[gene_name]]
  if (anyNA(dt$strand)) {
    stop("build_splice_windows(): missing strand for genes: ",
         paste(head(unique(dt$gene_name[is.na(dt$strand)]), 5), collapse = ", "))
  }
  setorder(dt, gene_name, start)

  dt[, n_exons := .N, by = gene_name]
  dt[, exon_index := seq_len(.N), by = gene_name]

  # Single-exon genes have no splice sites.
  dt <- dt[n_exons > 1L]
  if (nrow(dt) == 0) {
    stop("No multi-exon genes found. Check the GTF exon records.")
  }

  # For a plus-strand gene the donor is the exon end and the acceptor is the
  # exon start. Minus-strand genes reverse the two.
  # The first exon in transcription order has no acceptor; the last has no donor.
  dt[, is_first_in_tx := ifelse(strand == "-", exon_index == n_exons, exon_index == 1L)]
  dt[, is_last_in_tx  := ifelse(strand == "-", exon_index == 1L, exon_index == n_exons)]

  dt[, donor_pos    := ifelse(strand == "-", start, end)]
  dt[, acceptor_pos := ifelse(strand == "-", end, start)]

  donor_dt <- dt[is_last_in_tx == FALSE]
  acceptor_dt <- dt[is_first_in_tx == FALSE]

  make_windows <- function(d, pos_col) {
    gr <- GRanges(
      seqnames = d$chr,
      ranges = IRanges(start = pmax(1L, d[[pos_col]] - half_width),
                       end = d[[pos_col]] + half_width),
      gene_name = d$gene_name
    )
    reduce_by_gene(gr)
  }

  list(
    donor = make_windows(donor_dt, "donor_pos"),
    acceptor = make_windows(acceptor_dt, "acceptor_pos")
  )
}

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

#' Write a feature GRanges as BED6.
#'
#' @param gr GRanges with gene_name.
#' @param gene_strand Named character vector mapping gene symbol to strand.
#' @param out_path Destination BED path.
#' @param label Feature name used in the log line.
write_feature_bed <- function(gr, gene_strand, out_path, label) {
  if (length(gr) == 0) {
    stop("No intervals generated for feature type: ", label)
  }

  ranks <- rank_within_gene(gr, gene_strand)

  bed <- data.table(
    chr = as.character(seqnames(gr)),
    start = start(gr) - 1L,          # BED is 0-based half-open
    end = end(gr),
    name = paste0(mcols(gr)$gene_name, "|", ranks),
    score = 0L,
    strand = gene_strand[mcols(gr)$gene_name]
  )
  setorder(bed, chr, start)

  fwrite(bed, out_path, sep = "\t", col.names = FALSE, quote = FALSE)

  n_genes <- length(unique(mcols(gr)$gene_name))
  total_bp <- sum(width(gr))
  cat(sprintf("  %-20s %8s intervals  %6s genes  %12s bp  -> %s\n",
              label,
              format(length(gr), big.mark = ","),
              format(n_genes, big.mark = ","),
              format(total_bp, big.mark = ","),
              basename(out_path)))
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  opt <- parse_cli_args()
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

  gtf <- load_gtf_protein_coding(opt$gtf)

  gene_records <- gtf[mcols(gtf)$type == "gene"]
  if (length(gene_records) == 0) {
    stop("GTF has no 'gene' records after filtering.")
  }

  gene_strand <- as.character(strand(gene_records))
  names(gene_strand) <- mcols(gene_records)$gene_name

  # One symbol can carry several gene records. Every interval set in this
  # script is grouped by symbol, so the symbol needs a single strand. Records
  # of one symbol on both strands would give wrong ranks and would swap the
  # donor and acceptor windows for part of that symbol.
  strand_by_name <- unique(data.frame(gene_name = names(gene_strand),
                                      strand = unname(gene_strand),
                                      stringsAsFactors = FALSE))
  conflicting <- strand_by_name$gene_name[duplicated(strand_by_name$gene_name)]
  if (length(conflicting) > 0) {
    stop(length(conflicting), " gene symbols have records on both strands: ",
         paste(head(conflicting, 5), collapse = ", "))
  }

  gene_strand <- gene_strand[!duplicated(names(gene_strand))]
  cat(sprintf("  %s protein-coding genes\n",
              format(length(gene_strand), big.mark = ",")))

  cat("\nBuilding feature intervals...\n")

  # --- Exons ---------------------------------------------------------------
  exon_records <- gtf[mcols(gtf)$type == "exon"]
  if (length(exon_records) == 0) stop("GTF has no 'exon' records.")
  exons <- reduce_by_gene(exon_records)

  # --- Introns: gene body minus exons --------------------------------------
  genes_reduced <- reduce_by_gene(gene_records)
  genes_by_gene <- split(genes_reduced, mcols(genes_reduced)$gene_name)
  exons_by_gene <- split(exons, mcols(exons)$gene_name)

  # setdiff() on two GRangesList objects pairs their elements by position.
  # Equal name vectors are what makes position mean the same gene in both.
  if (!identical(names(genes_by_gene), names(exons_by_gene))) {
    no_exons <- base::setdiff(names(genes_by_gene), names(exons_by_gene))
    no_gene  <- base::setdiff(names(exons_by_gene), names(genes_by_gene))
    stop("Gene and exon records cover different gene sets, so introns cannot ",
         "be computed. Genes with no exon record: ", length(no_exons), " (",
         paste(head(no_exons, 5), collapse = ", "), "). ",
         "Exon gene names with no gene record: ", length(no_gene), " (",
         paste(head(no_gene, 5), collapse = ", "), ").")
  }

  introns <- GenomicRanges::setdiff(genes_by_gene, exons_by_gene)
  introns <- unlist(introns, use.names = TRUE)
  mcols(introns)$gene_name <- names(introns)
  names(introns) <- NULL
  introns <- sort(introns)

  # --- UTRs ----------------------------------------------------------------
  record_types <- unique(as.character(mcols(gtf)$type))
  if ("five_prime_UTR" %in% record_types && "three_prime_UTR" %in% record_types) {
    utr5 <- reduce_by_gene(gtf[mcols(gtf)$type == "five_prime_UTR"])
    utr3 <- reduce_by_gene(gtf[mcols(gtf)$type == "three_prime_UTR"])
  } else if ("UTR" %in% record_types) {
    # Older GENCODE releases emit a single UTR type. Each record is classified
    # against the coding bounds of its own transcript, then reduced per gene.
    cat("  GTF uses a single UTR type; splitting by position relative to CDS.\n")
    utr_split <- split_utr_by_transcript(gtf)
    utr5 <- reduce_by_gene(utr_split$five)
    utr3 <- reduce_by_gene(utr_split$three)
  } else {
    stop("GTF has no UTR records. Expected five_prime_UTR/three_prime_UTR ",
         "or a combined UTR type.")
  }

  # --- Splice sites --------------------------------------------------------
  splice <- build_splice_windows(exons, gene_strand, opt$splice_window)

  cat("\nWriting BED files to:", opt$out_dir, "\n")
  write_feature_bed(utr5, gene_strand,
                    file.path(opt$out_dir, "utr5_protein_coding.bed"), "5UTR")
  write_feature_bed(exons, gene_strand,
                    file.path(opt$out_dir, "exons_protein_coding.bed"), "Exon")
  write_feature_bed(splice$donor, gene_strand,
                    file.path(opt$out_dir, "splice_donor_protein_coding.bed"),
                    "SpliceSite_Donor")
  write_feature_bed(splice$acceptor, gene_strand,
                    file.path(opt$out_dir, "splice_acceptor_protein_coding.bed"),
                    "SpliceSite_Acceptor")
  write_feature_bed(introns, gene_strand,
                    file.path(opt$out_dir, "introns_protein_coding.bed"), "Intron")
  write_feature_bed(utr3, gene_strand,
                    file.path(opt$out_dir, "utr3_protein_coding.bed"), "3UTR")

  cat("\nFeature BED generation complete.\n")
}

main()
