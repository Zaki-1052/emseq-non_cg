# scripts/mca_gsea.R
# GO enrichment for gene-body mCA differential methylation results.
#
# Runs two complementary analyses on mca_differential_results.tsv:
#   GSEA  gseGO over every tested gene, ranked by log2_fc (no threshold)
#   ORA   enrichGO on the significant genes against the tested genes as universe,
#         run three ways: all significant, mCA-up only, mCA-down only
#
# Conventions (gseGO/enrichGO parameters, bitr mapping, ranked-list construction)
# follow section_61k_gsea_mecp2_k119ub.R and
# section_72_k119ub_neuronal_characterization.R in the biomodal downstream pipeline.
#
# Usage:
#   Rscript scripts/mca_gsea.R --results mca_results/mca_differential_results.tsv

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(svglite)
})

# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

option_list <- list(
  make_option("--results", type = "character", default = NULL,
              help = "Path to mca_differential_results.tsv [required]"),
  make_option("--outdir", type = "character", default = "mca_results/gsea",
              help = "Directory for output tables and figures [default: %default]"),
  make_option("--prefix", type = "character", default = "mca",
              help = "Output filename prefix [default: %default]"),
  make_option("--sig-basis", type = "character", default = "bonferroni", dest = "sig_basis",
              help = "Which correction defines the ORA gene set: bonferroni or fdr [default: %default]"),
  make_option("--ont", type = "character", default = "BP",
              help = "GO ontology: BP, MF, CC, or ALL [default: %default]"),
  make_option("--minsize", type = "integer", default = 15, dest = "minsize",
              help = "Minimum gene set size for GSEA [default: %default]"),
  make_option("--maxsize", type = "integer", default = 500, dest = "maxsize",
              help = "Maximum gene set size for GSEA [default: %default]"),
  make_option("--pcutoff", type = "double", default = 0.05, dest = "pcutoff",
              help = "p-value cutoff passed to gseGO/enrichGO; set to 1 to dump all terms [default: %default]"),
  make_option("--qcutoff", type = "double", default = 0.2, dest = "qcutoff",
              help = "q-value cutoff for enrichGO [default: %default]"),
  make_option("--showcat", type = "integer", default = 20, dest = "showcat",
              help = "Number of categories shown in each dotplot [default: %default]"),
  make_option("--neuronal-pattern", type = "character",
              default = "synap|neuron|axon|dendrit|nervous", dest = "neuronal_pattern",
              help = "Regex counted as a neuronal GO term in the summary [default: %default]"),
  make_option("--width", type = "double", default = 11,
              help = "Figure width, inches [default: %default]"),
  make_option("--height", type = "double", default = 9,
              help = "Figure height, inches [default: %default]"),
  make_option("--dpi", type = "integer", default = 300,
              help = "Raster resolution for PNG output [default: %default]")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  description = "GO enrichment (GSEA + ORA) for gene-body mCA differential results."
))

if (is.null(opt$results)) stop("--results is required")
if (!file.exists(opt$results)) stop("--results file not found: ", opt$results)
if (!opt$sig_basis %in% c("bonferroni", "fdr")) {
  stop("--sig-basis must be 'bonferroni' or 'fdr', got: ", opt$sig_basis)
}

# =============================================================================
# STYLE AND OUTPUT
# =============================================================================

# Mirrors theme_biomodal() from the biomodal downstream pipeline.
theme_mca <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, size = base_size),
      plot.caption = element_text(hjust = 0, size = base_size - 3, colour = "grey30"),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

save_multiformat <- function(plot, base_path, width, height, dpi) {
  dir.create(dirname(base_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(base_path, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(base_path, ".svg"), plot, width = width, height = height,
         device = svglite::svglite)
  ggsave(paste0(base_path, ".png"), plot, width = width, height = height,
         dpi = dpi, device = "png")
  cat(sprintf("  Saved: %s.{pdf,svg,png}\n", basename(base_path)))
  invisible(plot)
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
out <- function(suffix) file.path(opt$outdir, paste0(opt$prefix, "_", suffix))

# =============================================================================
# DATA
# =============================================================================

REQUIRED_COLS <- c("gene_name", "log2_fc", "mca_diff", "sig_bonferroni", "sig_fdr")

results <- read_tsv(opt$results, show_col_types = FALSE, progress = FALSE)
missing <- setdiff(REQUIRED_COLS, names(results))
if (length(missing) > 0) {
  stop("Results file is missing required columns: ", paste(missing, collapse = ", "))
}

sig_col <- if (opt$sig_basis == "bonferroni") "sig_bonferroni" else "sig_fdr"
results$is_sig <- results[[sig_col]]

basis_label <- if (opt$sig_basis == "bonferroni") "Bonferroni" else "BH FDR"

n_tested <- nrow(results)
sig_genes  <- results$gene_name[results$is_sig]
up_genes   <- results$gene_name[results$is_sig & results$mca_diff > 0]
down_genes <- results$gene_name[results$is_sig & results$mca_diff < 0]

cat(sprintf("Tested genes: %d\n", n_tested))
cat(sprintf("Significant (%s): %d  (%d up, %d down)\n",
            basis_label, length(sig_genes), length(up_genes), length(down_genes)))

# =============================================================================
# SYMBOL -> ENTREZ MAPPING
# =============================================================================

# Follows build_ranked_list() in section_61k: map to Entrez, keep the largest
# |value| per Entrez id when symbols collapse, then sort descending.
build_ranked_list <- function(symbols, values, label) {
  valid <- !is.na(symbols) & !is.na(values) & symbols != "" & is.finite(values)
  g <- symbols[valid]
  v <- values[valid]

  entrez <- suppressWarnings(
    bitr(g, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  )
  merged <- merge(
    data.frame(SYMBOL = g, value = v, stringsAsFactors = FALSE),
    entrez, by = "SYMBOL"
  )
  merged <- merged[order(-abs(merged$value)), ]
  merged <- merged[!duplicated(merged$ENTREZID), ]

  ranked <- setNames(merged$value, merged$ENTREZID)
  ranked <- sort(ranked, decreasing = TRUE)

  cat(sprintf("  %s: %d of %d genes mapped to Entrez (range %.3f to %.3f)\n",
              label, length(ranked), length(g), max(ranked), min(ranked)))
  ranked
}

map_to_entrez <- function(symbols, label) {
  entrez <- suppressWarnings(
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  )
  cat(sprintf("  %s: %d of %d genes mapped to Entrez\n",
              label, nrow(entrez), length(symbols)))
  unique(entrez$ENTREZID)
}

cat("\nMapping symbols to Entrez ids\n")
ranked <- build_ranked_list(results$gene_name, results$log2_fc, "Ranked list (log2_fc)")
universe_entrez <- map_to_entrez(results$gene_name, "ORA universe")
sig_entrez      <- map_to_entrez(sig_genes,  "ORA significant")
up_entrez       <- map_to_entrez(up_genes,   "ORA mCA-up")
down_entrez     <- map_to_entrez(down_genes, "ORA mCA-down")

# =============================================================================
# ENRICHMENT
# =============================================================================

# log2_fc is signed, so scoreType stays at the default "std" — the "pos" variant
# used in section_72 applies to unsigned signal magnitudes only.
run_gsea <- function(ranked_list, label) {
  cat(sprintf("\nRunning GSEA: %s\n", label))
  res <- gseGO(
    geneList      = ranked_list,
    OrgDb         = org.Mm.eg.db,
    ont           = opt$ont,
    minGSSize     = opt$minsize,
    maxGSSize     = opt$maxsize,
    pvalueCutoff  = opt$pcutoff,
    pAdjustMethod = "BH",
    verbose       = FALSE,
    seed          = TRUE,
    eps           = 0
  )
  res <- setReadable(res, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  cat(sprintf("  Terms returned: %d\n", nrow(res@result)))
  res
}

run_ora <- function(gene_entrez, label) {
  cat(sprintf("\nRunning ORA: %s (%d genes)\n", label, length(gene_entrez)))
  if (length(gene_entrez) == 0) {
    cat("  No genes; skipping\n")
    return(NULL)
  }
  res <- enrichGO(
    gene          = gene_entrez,
    universe      = universe_entrez,
    OrgDb         = org.Mm.eg.db,
    ont           = opt$ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = opt$pcutoff,
    qvalueCutoff  = opt$qcutoff,
    readable      = TRUE
  )
  cat(sprintf("  Terms returned: %d\n", nrow(res@result)))
  res
}

gsea_res <- run_gsea(ranked, sprintf("all tested genes ranked by log2_fc, GO %s", opt$ont))
ora_all  <- run_ora(sig_entrez,  sprintf("all significant (%s)", basis_label))
ora_up   <- run_ora(up_entrez,   "mCA higher in mutant")
ora_down <- run_ora(down_entrez, "mCA lower in mutant")

# =============================================================================
# TABLES
# =============================================================================

write_result <- function(res, suffix) {
  path <- paste0(out(suffix), ".tsv")
  df <- if (is.null(res)) data.frame() else as.data.frame(res@result)
  write_tsv(df, path)
  cat(sprintf("  Wrote %s (%d rows)\n", basename(path), nrow(df)))
  invisible(df)
}

cat("\nWriting tables\n")
gsea_df <- write_result(gsea_res, sprintf("gsea_go%s", tolower(opt$ont)))
all_df  <- write_result(ora_all,  sprintf("ora_all_go%s", tolower(opt$ont)))
up_df   <- write_result(ora_up,   sprintf("ora_up_go%s", tolower(opt$ont)))
down_df <- write_result(ora_down, sprintf("ora_down_go%s", tolower(opt$ont)))

# =============================================================================
# FIGURES
# =============================================================================

n_terms <- function(res) if (is.null(res)) 0 else nrow(res@result)

plot_gsea <- function(res, title) {
  # split/facet by NES sign so activated and suppressed sets read separately
  signs <- sign(res@result$NES)
  p <- if (length(unique(signs)) > 1) {
    dotplot(res, showCategory = opt$showcat, split = ".sign") +
      facet_grid(. ~ .sign)
  } else {
    dotplot(res, showCategory = opt$showcat)
  }
  p + labs(title = title,
           caption = sprintf("gseGO, GO %s, %d genes ranked by log2_fc, minGSSize %d, maxGSSize %d, BH",
                             opt$ont, length(ranked), opt$minsize, opt$maxsize)) +
    theme_mca()
}

plot_ora <- function(res, title, n_input) {
  dotplot(res, showCategory = opt$showcat) +
    labs(title = title,
         caption = sprintf("enrichGO, GO %s, %d genes vs %d-gene universe, BH, q < %g",
                           opt$ont, n_input, length(universe_entrez), opt$qcutoff)) +
    theme_mca()
}

cat("\nWriting figures\n")

if (n_terms(gsea_res) > 0) {
  save_multiformat(
    plot_gsea(gsea_res, sprintf("GSEA: mCA log2 fold change, GO %s", opt$ont)),
    out(sprintf("gsea_go%s_dotplot", tolower(opt$ont))),
    opt$width, opt$height, opt$dpi)
} else {
  cat("  GSEA returned no terms; no dotplot written\n")
}

ora_panels <- list(
  list(res = ora_all,  n = length(sig_entrez),  suffix = "ora_all",
       title = sprintf("ORA: all %s-significant genes, GO %s", basis_label, opt$ont)),
  list(res = ora_up,   n = length(up_entrez),   suffix = "ora_up",
       title = sprintf("ORA: mCA higher in mutant, GO %s", opt$ont)),
  list(res = ora_down, n = length(down_entrez), suffix = "ora_down",
       title = sprintf("ORA: mCA lower in mutant, GO %s", opt$ont))
)

for (panel in ora_panels) {
  if (n_terms(panel$res) > 0) {
    save_multiformat(
      plot_ora(panel$res, panel$title, panel$n),
      out(sprintf("%s_go%s_dotplot", panel$suffix, tolower(opt$ont))),
      opt$width, opt$height, opt$dpi)
  } else {
    cat(sprintf("  %s returned no terms; no dotplot written\n", panel$suffix))
  }
}

# =============================================================================
# SUMMARY
# =============================================================================

summarise_block <- function(df, label, con) {
  sig <- if (nrow(df) == 0) df else df[df$p.adjust < 0.05, , drop = FALSE]
  cat(sprintf("\n%s\n%s\n", label, strrep("-", nchar(label))), file = con)
  cat(sprintf("  Terms with q < 0.05: %d\n", nrow(sig)), file = con)
  if (nrow(sig) == 0) return(invisible(NULL))

  n_neuro <- sum(grepl(opt$neuronal_pattern, sig$Description, ignore.case = TRUE))
  cat(sprintf("  Matching /%s/: %d\n", opt$neuronal_pattern, n_neuro), file = con)

  if ("NES" %in% names(sig)) {
    cat(sprintf("  NES > 0: %d, NES < 0: %d\n",
                sum(sig$NES > 0), sum(sig$NES < 0)), file = con)
    ordered <- sig[order(-abs(sig$NES)), ]
    cat("  Top terms by |NES|:\n", file = con)
    for (i in seq_len(min(15, nrow(ordered)))) {
      cat(sprintf("    %-58s NES %+6.2f  q %.2e\n",
                  substr(ordered$Description[i], 1, 58),
                  ordered$NES[i], ordered$p.adjust[i]), file = con)
    }
  } else {
    ordered <- sig[order(sig$p.adjust), ]
    cat("  Top terms by q:\n", file = con)
    for (i in seq_len(min(15, nrow(ordered)))) {
      cat(sprintf("    %-58s %-9s q %.2e\n",
                  substr(ordered$Description[i], 1, 58),
                  ordered$GeneRatio[i], ordered$p.adjust[i]), file = con)
    }
  }
  invisible(NULL)
}

summary_path <- paste0(out("enrichment_summary"), ".txt")
con <- file(summary_path, open = "wt")
cat("Gene-body mCA GO enrichment\n===========================\n", file = con)
cat(sprintf("\nOntology:              GO %s\n", opt$ont), file = con)
cat(sprintf("Significance basis:    %s\n", basis_label), file = con)
cat(sprintf("Genes tested:          %d\n", n_tested), file = con)
cat(sprintf("Significant:           %d (%d up, %d down)\n",
            length(sig_genes), length(up_genes), length(down_genes)), file = con)
cat(sprintf("Ranked list (GSEA):    %d genes mapped to Entrez\n", length(ranked)), file = con)
cat(sprintf("ORA universe:          %d genes mapped to Entrez\n", length(universe_entrez)), file = con)

summarise_block(gsea_df, "GSEA (all genes ranked by log2_fc)", con)
summarise_block(all_df,  "ORA: all significant genes", con)
summarise_block(up_df,   "ORA: mCA higher in mutant", con)
summarise_block(down_df, "ORA: mCA lower in mutant", con)
close(con)

cat(sprintf("\nWrote %s\n", summary_path))
cat("\n")
writeLines(readLines(summary_path))
