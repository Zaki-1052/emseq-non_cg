# scripts/sections/_shared_config.R
#
# Shared configuration for the mCH downstream section pipeline.
#
# Every section script sources this file as its first executable statement.
# It does three things:
#   1. Resolves all file paths from paths.yaml for the active environment.
#   2. Pre-loads the data objects that most sections need.
#   3. Defines shared constants, colour palettes, and helper functions.
#
# Section scripts locate this file with the standard bootstrap block:
#
#   .file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
#   if (length(.file_arg) != 1) stop("Cannot resolve script path. Run with Rscript.")
#   SECTIONS_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
#   source(file.path(SECTIONS_DIR, "_shared_config.R"))
#
# The active environment comes from the EMSEQ_ENV variable and defaults to
# "local". SLURM wrappers set EMSEQ_ENV=expanse.

# =============================================================================
# PACKAGES
# =============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(scales)
  library(ggrepel)
  library(pheatmap)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(ChIPseeker)
  library(data.table)
})

# =============================================================================
# PATH RESOLUTION
# =============================================================================

# Find the repository root that holds paths.yaml.
# Priority: EMSEQ_CODE_DIR variable, then the sourcing script's location,
# then the working directory.
resolve_code_dir <- function() {
  env_dir <- Sys.getenv("EMSEQ_CODE_DIR", "")
  if (nzchar(env_dir)) {
    if (!file.exists(file.path(env_dir, "paths.yaml"))) {
      stop("EMSEQ_CODE_DIR is set to '", env_dir,
           "' but that directory has no paths.yaml")
    }
    return(normalizePath(env_dir))
  }

  # SECTIONS_DIR is defined by the calling section script's bootstrap block.
  if (exists("SECTIONS_DIR", envir = parent.frame())) {
    candidate <- dirname(dirname(get("SECTIONS_DIR", envir = parent.frame())))
    if (file.exists(file.path(candidate, "paths.yaml"))) {
      return(normalizePath(candidate))
    }
  }

  # Walk up from the working directory.
  candidate <- normalizePath(getwd())
  repeat {
    if (file.exists(file.path(candidate, "paths.yaml"))) return(candidate)
    parent <- dirname(candidate)
    if (parent == candidate) break
    candidate <- parent
  }

  stop("Cannot find paths.yaml. Set EMSEQ_CODE_DIR to the repository root.")
}

CODE_DIR <- resolve_code_dir()
EMSEQ_ENV <- Sys.getenv("EMSEQ_ENV", "local")

.paths_file <- file.path(CODE_DIR, "paths.yaml")
.all_paths <- yaml::read_yaml(.paths_file)

if (!EMSEQ_ENV %in% names(.all_paths)) {
  stop("EMSEQ_ENV='", EMSEQ_ENV, "' has no block in ", .paths_file,
       ". Available: ", paste(names(.all_paths), collapse = ", "))
}

.roots <- .all_paths[[EMSEQ_ENV]]
DATA_DIR    <- .roots$data_dir
HIC_DIR     <- .roots$hic_dir
BIGWIGS_DIR <- .roots$bigwigs_dir

for (.root_name in c("data_dir", "hic_dir", "bigwigs_dir")) {
  if (is.null(.roots[[.root_name]])) {
    stop("paths.yaml block '", EMSEQ_ENV, "' is missing key: ", .root_name)
  }
}

# --- Core pipeline outputs ---------------------------------------------------
# Everything under data/ and results/ lives in CODE_DIR (the repo). DATA_DIR
# holds only the large raw inputs (combined CH methylKit files) that are too
# big to sync to the laptop.

DATA_PATHS <- list(
  mch_results   = file.path(CODE_DIR, "results", "03_differential", "mch_differential_results.tsv"),
  sample_matrix = file.path(CODE_DIR, "results", "03_differential", "mch_sample_matrix.tsv"),
  feature_dir   = file.path(CODE_DIR, "results", "02b_features"),
  gene_bed      = file.path(CODE_DIR, "data", "gencode.vM25.mouse.genes.annotation.bed")
)

# --- Histone mark consensus peaks (data/chip_peaks/) -------------------------

CHIP_PATHS <- list(
  ctcf     = file.path(CODE_DIR, "data", "chip_peaks", "CTCF.bed"),
  h3k27ac  = file.path(CODE_DIR, "data", "chip_peaks", "H3K27ac.bed"),
  h3k27me3 = file.path(CODE_DIR, "data", "chip_peaks", "H3K27me3.bed"),
  h3k4me1  = file.path(CODE_DIR, "data", "chip_peaks", "H3K4me1.bed"),
  h3k4me3  = file.path(CODE_DIR, "data", "chip_peaks", "H3K4me3.bed"),
  bivalent = file.path(CODE_DIR, "data", "chip_peaks", "Bivalent.bed")
)

# --- chromHMM emission state segmentations (data/chromatin/) -----------------

CHROMATIN_PATHS <- list(
  active_enhancer = file.path(CODE_DIR, "data", "chromatin", "activeenhancer.bed"),
  active_promoter = file.path(CODE_DIR, "data", "chromatin", "activepromoter.bed"),
  bivalent        = file.path(CODE_DIR, "data", "chromatin", "bivalent.bed")
)

# --- Quantitative differential binding (data/diffbind/ and data/) ------------

DIFFBIND_PATHS <- list(
  atac            = file.path(CODE_DIR, "data", "diffbind", "ATAC_diffbind.txt"),
  k27ac           = file.path(CODE_DIR, "data", "diffbind", "K27ac_diffbind.txt"),
  k27me3          = file.path(CODE_DIR, "data", "diffbind", "K27me3_diffbind.txt"),
  k119ub          = file.path(CODE_DIR, "data", "h2aub_diffbind.txt"),
  mecp2           = file.path(CODE_DIR, "data", "mecp2_diffbind.txt"),
  k119ub_consensus = file.path(CODE_DIR, "data", "K119ub_consensus_v3.bed"),
  mecp2_consensus  = file.path(CODE_DIR, "data", "MeCP2_adult_concensus_peakset_Conc4.txt"),
  atac_up         = file.path(CODE_DIR, "data", "diffbind", "ATAC_up.bed"),
  atac_down       = file.path(CODE_DIR, "data", "diffbind", "ATAC_down.bed"),
  atac_cons_ctrl  = file.path(CODE_DIR, "data", "diffbind", "ATAC_consensus_ctrl.bed"),
  atac_cons_mut   = file.path(CODE_DIR, "data", "diffbind", "ATAC_consensus_mut.bed"),
  k119ub_ctrl     = file.path(CODE_DIR, "data", "diffbind", "K119ub_ctrl.bed"),
  k119ub_mut      = file.path(CODE_DIR, "data", "diffbind", "K119ub_mut.bed"),
  k27ac_ctrl      = file.path(CODE_DIR, "data", "diffbind", "K27ac_ctrl.bed"),
  k27ac_mut       = file.path(CODE_DIR, "data", "diffbind", "K27ac_mut.bed"),
  k119ub_gene_signal = file.path(CODE_DIR, "data", "k119ub_gene_signal.tsv")
)

# --- MeCP2 peak-level and aging data (data/mecp2/) ---------------------------

MECP2_PATHS <- list(
  annotated  = file.path(CODE_DIR, "data", "mecp2", "MeCP2_annotated.txt"),
  up         = file.path(CODE_DIR, "data", "mecp2", "MeCP2_up.bed"),
  down       = file.path(CODE_DIR, "data", "mecp2", "MeCP2_down.bed"),
  ctrl_aging = file.path(CODE_DIR, "data", "mecp2", "MeCP2_ctrl_aging_diffbind.txt"),
  mut_aging  = file.path(CODE_DIR, "data", "mecp2", "MeCP2_mut_aging_diffbind.txt")
)

# --- Hi-C loops, compartments, subcompartments (data/hic/) -------------------

# compartments: HOMER getDiffExpression output, 25kb bins at 50kb step.
#   Column names contain spaces, so readers must use check.names = FALSE.
#   The six per-sample PC1 columns match "bedGraph avg over given bp" in
#   order ctrl_M1, ctrl_M2, ctrl_M3, mut_M1, mut_M2, mut_M3.
# subcompartments: CALDER2 labels at 100kb with columns chr, bin_start,
#   bin_end, ctrl_label, mut_label, continous_rank_ctrl, continous_rank_mut,
#   label_changed. Bin starts are 1-based.
HIC_PATHS <- list(
  loops           = file.path(CODE_DIR, "data", "hic", "characterized_loops.tsv"),
  compartments    = file.path(CODE_DIR, "data", "hic", "diffcompartments.txt"),
  subcompartments = file.path(CODE_DIR, "data", "hic", "calder2_subcompartments_100kb.tsv")
)

SUBCOMPARTMENT_ORDER <- c("A.1", "A.2", "B.1", "B.2")

SUBCOMPARTMENT_COLORS <- c(
  "A.1" = "#CB181D",
  "A.2" = "#FB6A4A",
  "B.1" = "#6BAED6",
  "B.2" = "#2171B5"
)

# --- Neuronal gene sets (data/neuronal/) -------------------------------------

GENESET_PATHS <- list(
  neuronal = file.path(CODE_DIR, "data", "neuronal", "neuronal_gene_set_go_derived.tsv"),
  synapse  = file.path(CODE_DIR, "data", "neuronal", "synapse_axon_gene_set.tsv")
)

# --- Sub-gene feature intervals (data/features/) -----------------------------

FEATURE_TYPES <- c("5UTR", "Exon", "SpliceSite_Donor", "SpliceSite_Acceptor",
                   "Intron", "3UTR")

FEATURE_PATHS <- list(
  `5UTR`                = file.path(CODE_DIR, "data", "features", "utr5_protein_coding.bed"),
  Exon                  = file.path(CODE_DIR, "data", "features", "exons_protein_coding.bed"),
  SpliceSite_Donor      = file.path(CODE_DIR, "data", "features", "splice_donor_protein_coding.bed"),
  SpliceSite_Acceptor   = file.path(CODE_DIR, "data", "features", "splice_acceptor_protein_coding.bed"),
  Intron                = file.path(CODE_DIR, "data", "features", "introns_protein_coding.bed"),
  `3UTR`                = file.path(CODE_DIR, "data", "features", "utr3_protein_coding.bed")
)

# --- BigWig signal tracks ----------------------------------------------------

BIGWIG_PATHS <- list(
  k27me3_ctrl = file.path(BIGWIGS_DIR, "H3K27me3Ctrl.bw"),
  k27me3_mut  = file.path(BIGWIGS_DIR, "H3K27me3Mut.bw"),
  k27ac_ctrl  = file.path(BIGWIGS_DIR, "H3K27acCtrl.bw"),
  k27ac_mut   = file.path(BIGWIGS_DIR, "H3K27acMut.bw"),
  k119ub_ctrl = file.path(BIGWIGS_DIR, "H2AK119ubCtrl.bw"),
  k119ub_mut  = file.path(BIGWIGS_DIR, "H2AK119ubMut.bw"),
  atac_ctrl   = file.path(BIGWIGS_DIR, "ATACctrl.bw"),
  atac_mut    = file.path(BIGWIGS_DIR, "ATACmut.bw"),
  k4me3_ctrl  = file.path(BIGWIGS_DIR, "H3K4me3Ctrl.bw"),
  k4me3_mut   = file.path(BIGWIGS_DIR, "H3K4me3Mut.bw"),
  mecp2_ctrl  = file.path(BIGWIGS_DIR, "MeCP2Ctrl.bw"),
  mecp2_mut   = file.path(BIGWIGS_DIR, "MeCP2Mut.bw")
)

# --- Section output roots ----------------------------------------------------

RESULTS_ROOT <- file.path(CODE_DIR, "results", "sections")

OUTPUT_PATHS <- list(
  chromatin  = file.path(RESULTS_ROOT, "10_chromatin"),
  chip       = file.path(RESULTS_ROOT, "20_chip_integration"),
  hic        = file.path(RESULTS_ROOT, "30_hic"),
  permutation = file.path(RESULTS_ROOT, "40_permutation"),
  features   = file.path(RESULTS_ROOT, "50_features"),
  mecp2      = file.path(RESULTS_ROOT, "60_mecp2"),
  neuronal        = file.path(RESULTS_ROOT, "70_neuronal"),
  cross_modality  = file.path(RESULTS_ROOT, "80_cross_modality")
)

# Cross-section handoff files. Sections that produce these write them here;
# sections that consume them read from here.
HANDOFF_PATHS <- list(
  gene_level_all_marks = file.path(OUTPUT_PATHS$chip, "gene_level_all_marks.tsv"),
  mecp2_no_meth_genes  = file.path(OUTPUT_PATHS$mecp2, "mecp2_no_meth_genes.tsv"),
  neuronal_gene_set    = file.path(OUTPUT_PATHS$neuronal, "neuronal_gene_set.tsv"),
  chromatin_state      = file.path(OUTPUT_PATHS$chromatin, "gene_chromatin_state.tsv"),
  fisher_registry      = file.path(RESULTS_ROOT, "fisher_registry")
)

# =============================================================================
# CONSTANTS
# =============================================================================

Q_THRESHOLD <- 0.05
TSS_THRESHOLD <- 2000
CANONICAL_CHRS <- paste0("chr", c(1:19, "X", "Y"))

SAMPLE_META <- data.frame(
  sample_id = c("ctrl_M1", "ctrl_M2", "ctrl_F1", "ctrl_F2",
                "mut_M1", "mut_M2", "mut_F1", "mut_F2"),
  genotype  = c("ctrl", "ctrl", "ctrl", "ctrl", "mut", "mut", "mut", "mut"),
  sex       = c("M", "M", "F", "F", "M", "M", "F", "F"),
  stringsAsFactors = FALSE
)

KEY_GENES <- c("Syt1", "Zbtb20", "Trpm3", "Epha3", "Mcu",
               "Cntnap2", "Lpp", "Dlgap1", "Arhgap26", "Cdh8")

# Chromatin state is recorded twice per gene: once for the promoter window and
# once for the gene body. A gene can carry an active promoter over a
# Polycomb-marked body, and a single label cannot express that.
#
# The Biomodal pipeline used one label and separated promoter states from distal
# states by distance to the nearest TSS. That works for the small DMR intervals
# it classified. A gene body contains its own TSS, so the distance is zero for
# every gene and the distal states never populate. Hence the split.

PROMOTER_STATE_ORDER <- c("Active_Promoter", "Bivalent_Promoter",
                          "Repressed_Promoter", "Unmarked_Promoter")

PROMOTER_STATE_COLORS <- c(
  "Active_Promoter"    = "#e41a1c",
  "Bivalent_Promoter"  = "#984ea3",
  "Repressed_Promoter" = "#756bb1",
  "Unmarked_Promoter"  = "#999999"
)

BODY_STATE_ORDER <- c("Active", "Mixed", "Polycomb",
                      "Enhancer_Marked", "Unmarked")

BODY_STATE_COLORS <- c(
  "Active"          = "#377eb8",
  "Mixed"           = "#ff7f00",
  "Polycomb"        = "#4daf4a",
  "Enhancer_Marked" = "#a6cee3",
  "Unmarked"        = "#999999"
)

COLORS <- list(
  direction = c("Hypermethylated" = "#D7191C", "Hypomethylated" = "#2C7BB6"),
  genotype  = c("ctrl" = "#2166AC", "mut" = "#B2182B"),
  sex       = c("F" = "#E377C2", "M" = "#17BECF"),
  significant = c("Significant" = "#E41A1C", "Not Significant" = "grey70"),
  mecp2     = c("MeCP2 Up" = "#D95F02", "MeCP2 Down" = "#7570B3",
                "Not Significant" = "grey70"),
  atac      = c("ATAC Up" = "#E6AB02", "ATAC Down" = "#66A61E",
                "Not Significant" = "grey70"),
  k119ub    = c("K119ub Gained" = "#756BB1", "K119ub Lost" = "#74C476",
                "Shared" = "grey70", "Not Significant" = "grey70"),
  h3k27ac   = c("H3K27ac Gained" = "#FF7F00", "H3K27ac Lost" = "#1F78B4",
                "Shared" = "grey70", "Not Significant" = "grey70"),
  h3k27me3  = c("H3K27me3 Gained" = "#6A3D9A", "H3K27me3 Lost" = "#B2DF8A",
                "Shared" = "grey70", "Not Significant" = "grey70"),
  compartment = c("A" = "#E41A1C", "B" = "#377EB8"),
  quadrant  = c("Q1" = "#D73027", "Q2" = "#FC8D59",
                "Q3" = "#4575B4", "Q4" = "#91BFDB")
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Consistent ggplot theme for all section figures.
theme_emseq <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, size = base_size),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold")
    )
}

#' Load a BED file of peaks as GRanges.
#'
#' Adds a "chr" prefix when the file uses bare chromosome names, so every
#' GRanges in the pipeline shares the UCSC naming style of mch_results.
#'
#' @param bed_path Path to a BED file (at least 3 columns, no header).
#' @param mark_name Display name used in the load message.
#' @return GRanges
load_chip_peaks <- function(bed_path, mark_name = "ChIP") {
  if (!file.exists(bed_path)) {
    stop(mark_name, " BED file not found: ", bed_path)
  }

  df <- data.table::fread(bed_path, header = FALSE, sep = "\t",
                          select = 1:3, col.names = c("chr", "start", "end"))
  df <- df[!grepl("^(track|browser|#)", chr)]
  df[, chr := as.character(chr)]
  df[!grepl("^chr", chr), chr := paste0("chr", chr)]

  gr <- GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start + 1L, end = df$end)
  )

  cat(sprintf("  %s: %d peaks\n", mark_name, length(gr)))
  gr
}

#' Load DiffBind results with a flexible column schema.
#'
#' Handles the summit-appended layout (Summit_Chr/Summit_Start/Summit_End),
#' the raw DiffBind layout (seqnames/start/end), and files that already use
#' Chr/Start/End. Standardises to Chr, Start, End, Fold, FDR, p.value and adds
#' a direction column classifying each peak as Gained, Lost, or Unchanged.
#'
#' @param filepath Path to the DiffBind table.
#' @param mark_name Display name used in the load message.
#' @param fdr_threshold FDR cutoff for the direction classification.
#' @return data.frame with standardised columns
load_diffbind_flex <- function(filepath, mark_name = "Mark", fdr_threshold = Q_THRESHOLD) {
  if (!file.exists(filepath)) {
    stop(mark_name, " DiffBind file not found: ", filepath)
  }

  df <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "", fill = TRUE)

  if ("Summit_Chr" %in% colnames(df)) {
    df$Chr   <- df$Summit_Chr
    df$Start <- df$Summit_Start
    df$End   <- df$Summit_End
  } else if ("seqnames" %in% colnames(df)) {
    df$Chr   <- df$seqnames
    df$Start <- df$start
    df$End   <- df$end
  } else if (!"Chr" %in% colnames(df)) {
    stop("Unrecognised DiffBind column schema in ", filepath,
         ". Expected Summit_Chr, seqnames, or Chr.")
  }

  if (!"p.value" %in% colnames(df) && "p-value" %in% colnames(df)) {
    df$p.value <- df[["p-value"]]
  }

  required <- c("Chr", "Start", "End", "Fold", "FDR")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("Missing required columns in ", filepath, ": ",
         paste(missing, collapse = ", "))
  }

  df$Chr <- as.character(df$Chr)

  na_coords <- is.na(df$Chr) | is.na(df$Start) | is.na(df$End)
  if (any(na_coords)) {
    cat(sprintf("  %s: dropping %d peaks with NA coordinates\n",
                mark_name, sum(na_coords)))
    df <- df[!na_coords, , drop = FALSE]
  }

  needs_prefix <- !grepl("^chr", df$Chr)
  df$Chr[needs_prefix] <- paste0("chr", df$Chr[needs_prefix])

  df$direction <- "Unchanged"
  df$direction[df$FDR < fdr_threshold & df$Fold > 0] <- "Gained"
  df$direction[df$FDR < fdr_threshold & df$Fold < 0] <- "Lost"

  cat(sprintf("  %s: %d peaks (%d gained, %d lost at FDR<%.2f)\n",
              mark_name, nrow(df),
              sum(df$direction == "Gained"), sum(df$direction == "Lost"),
              fdr_threshold))
  df
}

#' Convert a standardised DiffBind data.frame to GRanges.
#'
#' @param db data.frame from load_diffbind_flex()
#' @return GRanges carrying Fold, FDR, and direction metadata
diffbind_to_granges <- function(db) {
  GRanges(
    seqnames = db$Chr,
    ranges = IRanges(start = db$Start, end = db$End),
    Fold = db$Fold,
    FDR = db$FDR,
    direction = db$direction
  )
}

#' Compute histone mark overlap booleans for a query GRanges.
#'
#' @param query_gr GRanges to test.
#' @param chip_peaks_list Named list of GRanges. Names become column prefixes.
#' @return data.frame of logical columns named "<Mark>_overlap"
compute_chip_overlaps <- function(query_gr, chip_peaks_list) {
  overlaps <- lapply(chip_peaks_list, function(peaks) {
    countOverlaps(query_gr, peaks) > 0
  })
  names(overlaps) <- paste0(names(chip_peaks_list), "_overlap")
  as.data.frame(overlaps)
}

#' Check that an overlap table carries the five marks the classifiers read.
require_mark_overlaps <- function(overlaps, caller) {
  required <- c("h3k27ac_overlap", "h3k27me3_overlap", "h3k4me1_overlap",
                "h3k4me3_overlap", "bivalent_overlap")
  missing <- setdiff(required, colnames(overlaps))
  if (length(missing) > 0) {
    stop(caller, "() needs columns: ", paste(missing, collapse = ", "),
         "\nName the peak list elements ctcf, h3k27ac, h3k27me3, h3k4me1, ",
         "h3k4me3, bivalent when calling compute_chip_overlaps().")
  }
  invisible(TRUE)
}

#' Classify the promoter state of each gene.
#'
#' Call with overlaps computed over the promoter window (gene_promoters), not
#' over the gene body.
#'
#' Priority order:
#'   1. Bivalent_Promoter  a bivalent peak, or H3K4me3 together with H3K27me3
#'   2. Active_Promoter    H3K4me3
#'   3. Repressed_Promoter H3K27me3 without H3K27ac
#'   4. Unmarked_Promoter  none of the above
#'
#' @param overlaps data.frame from compute_chip_overlaps() over promoters.
#' @return factor with levels PROMOTER_STATE_ORDER
classify_promoter_state <- function(overlaps) {
  require_mark_overlaps(overlaps, "classify_promoter_state")

  h3k27ac  <- overlaps$h3k27ac_overlap
  h3k27me3 <- overlaps$h3k27me3_overlap
  h3k4me3  <- overlaps$h3k4me3_overlap
  bivalent <- overlaps$bivalent_overlap

  state <- rep("Unmarked_Promoter", nrow(overlaps))

  is_bivalent <- bivalent | (h3k4me3 & h3k27me3)
  state[is_bivalent] <- "Bivalent_Promoter"

  is_active <- !is_bivalent & h3k4me3
  state[is_active] <- "Active_Promoter"

  is_repressed <- !is_bivalent & !is_active & h3k27me3 & !h3k27ac
  state[is_repressed] <- "Repressed_Promoter"

  factor(state, levels = PROMOTER_STATE_ORDER)
}

#' Classify the chromatin state of each gene body.
#'
#' Call with overlaps computed over the gene body (gene_bodies).
#'
#' Priority order:
#'   1. Mixed            H3K27me3 together with H3K27ac
#'   2. Polycomb         H3K27me3 without H3K27ac
#'   3. Active           H3K27ac without H3K27me3
#'   4. Enhancer_Marked  H3K4me1 without H3K27ac or H3K27me3
#'   5. Unmarked         none of the above
#'
#' @param overlaps data.frame from compute_chip_overlaps() over gene bodies.
#' @return factor with levels BODY_STATE_ORDER
classify_body_state <- function(overlaps) {
  require_mark_overlaps(overlaps, "classify_body_state")

  h3k27ac  <- overlaps$h3k27ac_overlap
  h3k27me3 <- overlaps$h3k27me3_overlap
  h3k4me1  <- overlaps$h3k4me1_overlap

  state <- rep("Unmarked", nrow(overlaps))

  is_mixed <- h3k27me3 & h3k27ac
  state[is_mixed] <- "Mixed"

  is_polycomb <- !is_mixed & h3k27me3
  state[is_polycomb] <- "Polycomb"

  is_active <- !is_mixed & !is_polycomb & h3k27ac
  state[is_active] <- "Active"

  is_enhancer <- !is_mixed & !is_polycomb & !is_active & h3k4me1
  state[is_enhancer] <- "Enhancer_Marked"

  factor(state, levels = BODY_STATE_ORDER)
}

#' Assign each point to a quadrant from two effect-size vectors.
#'
#' Q1 both up, Q2 x down and y up, Q3 both down, Q4 x up and y down.
#'
#' @param dx Numeric vector for the x axis.
#' @param dy Numeric vector for the y axis.
#' @return character vector of quadrant labels
assign_quadrant <- function(dx, dy) {
  ifelse(dx > 0 & dy > 0, "Q1",
    ifelse(dx <= 0 & dy > 0, "Q2",
      ifelse(dx <= 0 & dy <= 0, "Q3", "Q4")))
}

#' Collapse ChIPseeker-annotated peaks to one row per gene.
#'
#' Two collapse rules:
#'   nearest_tss        keep the peak closest to the TSS
#'   median_significant median Fold of peaks passing fdr_threshold, falling
#'                      back to the median of all peaks when none pass
#'
#' @param annotated_peaks data.frame from ChIPseeker with SYMBOL, Fold, FDR,
#'   and (for nearest_tss) distanceToTSS columns.
#' @param method "nearest_tss" or "median_significant".
#' @param fdr_threshold FDR cutoff used by median_significant.
#' @param prefix Column-name prefix for the returned summary columns.
#' @return data.frame with one row per gene
aggregate_diffbind_by_gene <- function(annotated_peaks, method = "nearest_tss",
                                       fdr_threshold = Q_THRESHOLD,
                                       prefix = "mark") {
  if (!method %in% c("nearest_tss", "median_significant")) {
    stop("aggregate_diffbind_by_gene(): method must be 'nearest_tss' or ",
         "'median_significant', got '", method, "'")
  }

  required <- c("SYMBOL", "Fold", "FDR")
  if (method == "nearest_tss") required <- c(required, "distanceToTSS")
  missing <- setdiff(required, colnames(annotated_peaks))
  if (length(missing) > 0) {
    stop("aggregate_diffbind_by_gene() needs columns: ",
         paste(missing, collapse = ", "))
  }

  df <- annotated_peaks[!is.na(annotated_peaks$SYMBOL), , drop = FALSE]

  if (method == "nearest_tss") {
    out <- df %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::summarise(
        fold      = Fold[which.min(abs(distanceToTSS))],
        fdr       = FDR[which.min(abs(distanceToTSS))],
        min_fdr   = min(FDR, na.rm = TRUE),
        n_peaks   = dplyr::n(),
        n_sig     = sum(FDR < fdr_threshold, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    out <- df %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::summarise(
        fold = if (any(FDR < fdr_threshold, na.rm = TRUE)) {
          median(Fold[FDR < fdr_threshold], na.rm = TRUE)
        } else {
          median(Fold, na.rm = TRUE)
        },
        fdr     = min(FDR, na.rm = TRUE),
        min_fdr = min(FDR, na.rm = TRUE),
        n_peaks = dplyr::n(),
        n_sig   = sum(FDR < fdr_threshold, na.rm = TRUE),
        .groups = "drop"
      )
  }

  out$has_sig <- out$n_sig > 0
  colnames(out) <- c("gene_name",
                     paste0(prefix, c("_fold", "_fdr", "_min_fdr",
                                      "_n_peaks", "_n_sig", "_has_sig")))
  as.data.frame(out)
}

#' Annotate a DiffBind table to genes with ChIPseeker.
#'
#' @param db data.frame from load_diffbind_flex().
#' @param mark_name Display name used in messages.
#' @param tss_region Promoter window passed to annotatePeak().
#' @return data.frame of annotated peaks with SYMBOL and distanceToTSS
annotate_peaks_to_genes <- function(db, mark_name = "Mark",
                                    tss_region = c(-3000, 3000)) {
  gr <- diffbind_to_granges(db)
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

  cat(sprintf("  Annotating %s peaks to genes...\n", mark_name))
  anno <- ChIPseeker::annotatePeak(gr, TxDb = txdb, tssRegion = tss_region,
                                   annoDb = "org.Mm.eg.db", verbose = FALSE)
  as.data.frame(anno)
}

#' Run a gene-level Fisher's exact test and record it for permutation validation.
#'
#' This is the single entry point for every gene-level Fisher test in the
#' pipeline. It runs the test, writes the gene table behind it, and appends one
#' row to the shared registry. Section 40_04 reads that registry, reloads each
#' gene table, and validates the odds ratio by chromosome-stratified label
#' shuffle.
#'
#' Sections must call this rather than fisher.test() directly for gene-level
#' 2x2 tests, so that every such test reaches 40_04.
#'
#' @param section Section identifier, for example "10_01".
#' @param test_id Identifier for this test within the section.
#' @param description What the test asks, in one sentence.
#' @param gene_df data.frame with one row per gene. Needs gene_name, chr, and
#'   the two logical columns named by row_var and col_var.
#' @param row_var Name of the logical column forming the rows of the 2x2 table.
#' @param col_var Name of the logical column forming the columns.
#' @param output_dir Section output directory. The gene table is written to a
#'   fisher_tables subdirectory of it.
#' @param registry_path Registry file. Defaults to HANDOFF_PATHS$fisher_registry.
#' @return the htest object from fisher.test()
register_fisher_test <- function(section, test_id, description,
                                 gene_df, row_var, col_var, output_dir,
                                 registry_path = HANDOFF_PATHS$fisher_registry) {
  required <- c("gene_name", "chr", row_var, col_var)
  missing <- setdiff(required, colnames(gene_df))
  if (length(missing) > 0) {
    stop("register_fisher_test(): gene_df is missing columns: ",
         paste(missing, collapse = ", "))
  }

  keep <- !is.na(gene_df[[row_var]]) & !is.na(gene_df[[col_var]])
  gene_df <- gene_df[keep, c("gene_name", "chr", row_var, col_var), drop = FALSE]

  if (!is.logical(gene_df[[row_var]]) || !is.logical(gene_df[[col_var]])) {
    stop("register_fisher_test(): ", row_var, " and ", col_var,
         " must both be logical columns.")
  }
  if (nrow(gene_df) == 0) {
    stop("register_fisher_test(): no genes left after dropping NA for ",
         section, ":", test_id)
  }

  tab <- table(factor(gene_df[[row_var]], levels = c(TRUE, FALSE)),
               factor(gene_df[[col_var]], levels = c(TRUE, FALSE)))
  ft <- fisher.test(tab)

  table_dir <- file.path(output_dir, "fisher_tables")
  gene_table_path <- file.path(table_dir,
                               sprintf("%s_%s_genes.tsv", section, test_id))
  write_section_table(gene_df, gene_table_path)

  row <- data.frame(
    section = section,
    test_id = test_id,
    description = description,
    row_var = row_var,
    col_var = col_var,
    gene_table_path = gene_table_path,
    n_genes = nrow(gene_df),
    n_row_true = sum(gene_df[[row_var]]),
    n_col_true = sum(gene_df[[col_var]]),
    n_both_true = sum(gene_df[[row_var]] & gene_df[[col_var]]),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )

  # Each test writes its own shard so concurrent SLURM jobs never collide.
  # Section 40_04 reads every shard with list.files() and rbinds them.
  dir.create(registry_path, recursive = TRUE, showWarnings = FALSE)
  shard_path <- file.path(registry_path,
                          sprintf("%s_%s.tsv", section, test_id))
  write_section_table(row, shard_path)

  cat(sprintf("  Fisher %s:%s: OR=%.3f, p=%.3g (n=%s genes) [registered]\n",
              section, test_id, unname(ft$estimate), ft$p.value,
              format(nrow(gene_df), big.mark = ",")))
  ft
}

#' Patch regioneReloaded::chooseHclustMet for matrices with two or fewer rows.
#'
#' chooseHclustMet picks an hclust method by cophenetic correlation. With two
#' or fewer rows, dist() returns a single value, cor() returns NA, and the
#' method selection indexes with integer(0). This replacement uses the first
#' method directly in that case.
#'
#' Call after library(regioneReloaded) in the permutation sections.
patch_chooseHclustMet <- function() {
  safe_fn <- function(GM, scale = TRUE, vecMet = NULL, distHC = "euclidean") {
    if (scale == TRUE) GM <- scale(GM)
    if (is.null(vecMet)) {
      vecMet <- c("complete", "average", "single", "ward.D2",
                  "median", "centroid", "mcquitty")
    }
    mat_dist <- stats::dist(x = GM, method = distHC)

    if (nrow(GM) <= 2) {
      model <- stats::hclust(d = mat_dist, method = vecMet[1])
      methods::show(paste0("method for hclustering (<=2 elements, skipping cophenetic): ",
                           vecMet[1]))
      return(model)
    }

    resMetList <- lapply(seq_along(vecMet), FUN = function(i, mat_dist, vecMet) {
      stats::hclust(d = mat_dist, method = vecMet[[i]])
    }, mat_dist, vecMet)
    names(resMetList) <- vecMet
    resMetVec <- unlist(lapply(seq_along(resMetList), FUN = function(i, mat_dist, resMetList) {
      stats::cor(x = mat_dist, stats::cophenetic(resMetList[[i]]))
    }, mat_dist, resMetList))
    names(resMetVec) <- vecMet
    name_model <- vecMet[which(resMetVec == max(resMetVec))]
    if (length(name_model) > 1) name_model <- name_model[1]
    model <- resMetList[[name_model]]
    methods::show(paste0("method selected for hclustering: ", name_model))
    methods::show(resMetVec)
    model
  }
  assignInNamespace("chooseHclustMet", safe_fn, ns = "regioneReloaded")
  cat("  Patched regioneReloaded::chooseHclustMet for small-matrix safety.\n")
}

#' Load a one-column gene set file.
#'
#' @param path TSV with a header and a single "gene" column.
#' @param set_name Display name used in the load message.
#' @return character vector of unique gene symbols
load_gene_set <- function(path, set_name = "gene set") {
  if (!file.exists(path)) stop(set_name, " file not found: ", path)
  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "")
  if (!"gene" %in% colnames(df)) {
    stop(set_name, " file must have a 'gene' column: ", path)
  }
  genes <- unique(df$gene[!is.na(df$gene) & nzchar(df$gene)])
  cat(sprintf("  %s: %d genes\n", set_name, length(genes)))
  genes
}

#' Write a section result table as TSV.
#'
#' Stops when a character column holds a newline. Such a column is figure text,
#' not data: writing it unquoted splits each row across physical lines and the
#' file no longer parses. Build figure text at the plot with group_label() and
#' keep it out of the table.
#'
#' @param df data.frame to write.
#' @param path Destination TSV path.
#' @return the path, invisibly
write_section_table <- function(df, path) {
  df <- as.data.frame(df)

  has_newline <- vapply(df, function(col) {
    is.character(col) && any(grepl("\n", col, fixed = TRUE), na.rm = TRUE)
  }, logical(1))

  if (any(has_newline)) {
    stop("Cannot write ", basename(path), ": column(s) ",
         paste(names(df)[has_newline], collapse = ", "),
         " contain a newline, so each row would break across physical lines. ",
         "That column is figure text. Build it at the plot with group_label() ",
         "and drop it before writing.")
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  Wrote: %s (%s rows, %d cols)\n", basename(path),
              format(nrow(df), big.mark = ","), ncol(df)))
  invisible(path)
}

#' Summarise the distribution of one value column per group.
#'
#' Returns data only. The text drawn on a figure comes from group_label().
#' Keeping the two apart means a summary table is always writable as TSV.
#'
#' The first column is named after group_col, so a plot layer can join on it
#' by name.
#'
#' @param df data.frame holding the values.
#' @param group_col Name of the grouping column.
#' @param value_col Name of the numeric column.
#' @return data.frame with group_col, n, median, mean, q25, q75
summarise_groups <- function(df, group_col, value_col) {
  out <- df %>%
    dplyr::filter(!is.na(.data[[value_col]])) %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median = median(.data[[value_col]]),
      mean = mean(.data[[value_col]]),
      q25 = unname(quantile(.data[[value_col]], 0.25)),
      q75 = unname(quantile(.data[[value_col]], 0.75)),
      .groups = "drop"
    ) %>%
    as.data.frame()

  colnames(out)[1] <- group_col
  out
}

#' Build the on-plot annotation text for a summarise_groups() result.
#'
#' Every distribution figure in this pipeline shows n and the median, so the
#' numbers behind each group are readable without opening the table. Assign the
#' result to a column at the plot site, never in the table that gets written.
#'
#' @param summary_df data.frame from summarise_groups().
#' @param digits Digits for the median.
#' @param unit Text appended to the median, for example " kb".
#' @return character vector, one entry per group
group_label <- function(summary_df, digits = 3, unit = "") {
  required <- c("n", "median")
  missing <- setdiff(required, colnames(summary_df))
  if (length(missing) > 0) {
    stop("group_label() needs columns: ", paste(missing, collapse = ", "),
         ". Pass the result of summarise_groups().")
  }
  sprintf("n = %s\nmed = %s%s",
          format(summary_df$n, big.mark = ","),
          format(round(summary_df$median, digits), nsmall = digits),
          unit)
}

#' Build the annotation text for a single group from its raw values.
#'
#' For a one-off annotation where no summarise_groups() result exists.
#'
#' @param values Numeric vector for the group.
#' @param digits Digits for the median.
#' @return single string, for example "n = 412\nmed = 0.031"
group_annotation <- function(values, digits = 3) {
  values <- values[!is.na(values)]
  sprintf("n = %s\nmed = %s",
          format(length(values), big.mark = ","),
          format(round(median(values), digits), nsmall = digits))
}

# =============================================================================
# MULTI-FORMAT PLOT OUTPUT
# =============================================================================

source(file.path(CODE_DIR, "scripts", "utils", "multi_format_output.R"))

# =============================================================================
# DATA LOADING
# =============================================================================

cat("================================================================================\n")
cat("EM-seq mCH SECTION PIPELINE — SHARED CONFIG\n")
cat("================================================================================\n")
cat("Environment: ", EMSEQ_ENV, "\n", sep = "")
cat("Code dir:    ", CODE_DIR, "\n", sep = "")
cat("Data dir:    ", DATA_DIR, "\n", sep = "")
cat("Results root:", RESULTS_ROOT, "\n")
cat("\n")

# --- Gene-level mCH differential results -------------------------------------

if (!file.exists(DATA_PATHS$mch_results)) {
  stop("mCH differential results not found: ", DATA_PATHS$mch_results)
}

cat("Loading mCH differential results...\n")
mch_results <- read.table(DATA_PATHS$mch_results, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, quote = "")

.required_mch_cols <- c("gene_name", "gene_id", "chr", "start", "end", "strand",
                        "gene_length", "mch_ctrl", "mch_mut", "mch_diff",
                        "edger_logFC", "edger_fdr", "sig_fdr005")
.missing_mch <- setdiff(.required_mch_cols, colnames(mch_results))
if (length(.missing_mch) > 0) {
  stop("mch_differential_results.tsv is missing columns: ",
       paste(.missing_mch, collapse = ", "))
}

mch_results$mch_sig <- as.logical(mch_results$sig_fdr005)
# Direction comes from the model-adjusted logFC, not the raw mch_diff, so the
# direction label is consistent with the model that determined significance.
mch_results$mch_direction <- ifelse(mch_results$edger_logFC > 0,
                                    "Hypermethylated", "Hypomethylated")
mch_results$mch_hyper <- mch_results$mch_sig & mch_results$edger_logFC > 0
mch_results$mch_hypo  <- mch_results$mch_sig & mch_results$edger_logFC < 0
mch_results$neg_log10_fdr <- -log10(pmax(mch_results$edger_fdr, 1e-300))

cat(sprintf("  %s genes tested, %s significant (%s hyper, %s hypo) at FDR<%.2f\n",
            format(nrow(mch_results), big.mark = ","),
            format(sum(mch_results$mch_sig), big.mark = ","),
            format(sum(mch_results$mch_hyper), big.mark = ","),
            format(sum(mch_results$mch_hypo), big.mark = ","),
            Q_THRESHOLD))

# Gene body coordinates for the tested universe. Sections use this rather than
# the full GENCODE annotation so that overlap tests run on genes that have mCH
# measurements. DATA_PATHS$gene_bed holds the complete annotation when a
# section needs genes outside the tested set.
gene_bodies <- GRanges(
  seqnames = mch_results$chr,
  ranges = IRanges(start = mch_results$start + 1L, end = mch_results$end),
  strand = mch_results$strand,
  gene_name = mch_results$gene_name,
  gene_id = mch_results$gene_id,
  gene_length = mch_results$gene_length,
  mch_diff = mch_results$mch_diff,
  edger_logFC = mch_results$edger_logFC,
  edger_fdr = mch_results$edger_fdr,
  mch_sig = mch_results$mch_sig,
  mch_direction = mch_results$mch_direction
)

# Promoter windows, TSS +/- TSS_THRESHOLD, strand-aware. Sections classify
# promoter state over these and body state over gene_bodies. Both are parallel
# to mch_results, so a classification can be bound back by position.
gene_promoters <- promoters(gene_bodies,
                            upstream = TSS_THRESHOLD,
                            downstream = TSS_THRESHOLD)
# Genes near a chromosome start give a promoter window beginning below 1.
start(gene_promoters) <- pmax(1L, start(gene_promoters))

# --- Differential binding tables ---------------------------------------------

cat("\nLoading DiffBind results...\n")
k119ub_diffbind <- load_diffbind_flex(DIFFBIND_PATHS$k119ub, "H2AK119ub")
mecp2_diffbind  <- load_diffbind_flex(DIFFBIND_PATHS$mecp2,  "MeCP2")
atac_diffbind   <- load_diffbind_flex(DIFFBIND_PATHS$atac,   "ATAC")
k27ac_diffbind  <- load_diffbind_flex(DIFFBIND_PATHS$k27ac,  "H3K27ac")
k27me3_diffbind <- load_diffbind_flex(DIFFBIND_PATHS$k27me3, "H3K27me3")

DIFFBIND_TABLES <- list(
  atac   = atac_diffbind,
  k27ac  = k27ac_diffbind,
  k27me3 = k27me3_diffbind,
  k119ub = k119ub_diffbind
)

# --- Consensus peak sets ------------------------------------------------------

cat("\nLoading consensus peak sets...\n")
k119ub_consensus <- load_chip_peaks(DIFFBIND_PATHS$k119ub_consensus, "K119ub consensus")
mecp2_consensus  <- load_chip_peaks(DIFFBIND_PATHS$mecp2_consensus,  "MeCP2 consensus")

cat("\nShared config loaded.\n\n")
