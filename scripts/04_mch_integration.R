# scripts/04_mch_integration.R
#
# Integration of mCH differential results with K119ub, MeCP2, neuronal
# gene sets, gene length, and CG methylation/hydroxymethylation data.
#
# Sub-analyses:
#   4a. K119ub overlap and directional concordance
#   4b. MeCP2 overlap and directional concordance (parallel to 4a)
#   4c. Neuronal gene set enrichment and gene-length interaction
#   4d. CG mC/hmC x mCH cross-modality comparison
#
# Usage:
#   Rscript mch_integration.R \
#       --diff-results /path/to/mch_differential_results.tsv \
#       --sample-matrix /path/to/mch_sample_matrix.tsv \
#       --gene-bed /path/to/gene_bodies.protein_coding.bed \
#       --k119ub-consensus /path/to/K119ub_consensus_v3.bed \
#       --k119ub-diffbind /path/to/K119ub_diffbind_results_summit_appended_ap.txt \
#       --mecp2-consensus /path/to/MeCP2_adult_concensus_peakset_Conc4.txt \
#       --mecp2-diffbind /path/to/260804_MeCP2_adult_diffbind_results_filtered_Conc4.txt \
#       --neuronal-genes /path/to/72_neuronal_gene_set_go_derived.tsv \
#       --mc-dmr /path/to/DMR_mc_control__mutant.bed \
#       --hmc-dmr /path/to/DMR_hmc_control__mutant.bed \
#       --out-dir /path/to/results/integration \
#       --plot-dir /path/to/results/plots

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(clinfun)
  library(fgsea)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

script_dir <- if (interactive()) "scripts" else dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "utils", "multi_format_output.R"))

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

parse_cli_args <- function() {
  option_list <- list(
    make_option("--diff-results", type = "character", dest = "diff_results"),
    make_option("--sample-matrix", type = "character", dest = "sample_matrix"),
    make_option("--gene-bed", type = "character", dest = "gene_bed"),
    make_option("--k119ub-consensus", type = "character", dest = "k119ub_consensus"),
    make_option("--k119ub-diffbind", type = "character", dest = "k119ub_diffbind"),
    make_option("--mecp2-consensus", type = "character", dest = "mecp2_consensus"),
    make_option("--mecp2-diffbind", type = "character", dest = "mecp2_diffbind"),
    make_option("--neuronal-genes", type = "character", dest = "neuronal_genes"),
    make_option("--mc-dmr", type = "character", dest = "mc_dmr"),
    make_option("--hmc-dmr", type = "character", dest = "hmc_dmr"),
    make_option("--out-dir", type = "character", dest = "out_dir"),
    make_option("--plot-dir", type = "character", dest = "plot_dir"),
    make_option("--fdr-threshold", type = "double", dest = "fdr_threshold",
                default = 0.05,
                help = "FDR threshold for DiffBind classification [default %default]")
  )

  opt <- parse_args(OptionParser(option_list = option_list))

  required <- c("diff_results", "sample_matrix", "gene_bed",
                "k119ub_consensus", "k119ub_diffbind",
                "mecp2_consensus", "mecp2_diffbind",
                "neuronal_genes", "mc_dmr", "hmc_dmr",
                "out_dir", "plot_dir")
  for (req in required) {
    if (is.null(opt[[req]])) stop("Missing: --", gsub("_", "-", req))
  }
  opt
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

save_tsv <- function(dt, path) {
  fwrite(dt, path, sep = "\t")
  cat("  Saved:", path, "\n")
}

save_plot <- function(p, path, w = 7, h = 6) {
  save_multiformat_ggplot(p, path, width = w, height = h)
}

load_bed3 <- function(path) {
  df <- fread(path, header = FALSE, select = 1:3,
              col.names = c("chr", "start", "end"))
  GRanges(seqnames = df$chr, ranges = IRanges(start = df$start + 1L, end = df$end))
}

gene_bed_to_gr <- function(gene_bed) {
  GRanges(seqnames = gene_bed$chr,
          ranges = IRanges(start = gene_bed$start + 1L, end = gene_bed$end))
}

classify_diffbind <- function(db, fdr_thresh) {
  db$direction <- ifelse(db$FDR < fdr_thresh & db$Fold > 0, "gained",
                  ifelse(db$FDR < fdr_thresh & db$Fold < 0, "lost",
                         "unchanged"))
  db
}

# ---------------------------------------------------------------------------
# 4a. K119ub overlap
# ---------------------------------------------------------------------------

run_k119ub <- function(results, gene_bed, gene_gr, opt) {
  cat("\n=== 4a. K119ub Overlap ===\n")

  # Consensus overlap
  k119ub_gr <- load_bed3(opt$k119ub_consensus)
  cat("  K119ub consensus peaks:", length(k119ub_gr), "\n")

  overlap <- countOverlaps(gene_gr, k119ub_gr) > 0
  results$k119ub_overlap <- overlap[match(results$gene_name, gene_bed$gene_name)]

  sig <- !is.na(results$sig_fdr005) & results$sig_fdr005
  tbl <- table(significant = sig, k119ub = results$k119ub_overlap)
  ft <- fisher.test(tbl)
  cat(sprintf("  Fisher's exact: OR=%.3f, p=%.2e\n", ft$estimate, ft$p.value))

  fisher_dt <- data.table(
    test = "mCH_sig_x_K119ub_overlap",
    odds_ratio = as.numeric(ft$estimate),
    p_value = ft$p.value,
    ci_lower = ft$conf.int[1],
    ci_upper = ft$conf.int[2],
    n_sig_overlap = tbl["TRUE", "TRUE"],
    n_sig_no_overlap = tbl["TRUE", "FALSE"],
    n_nonsig_overlap = tbl["FALSE", "TRUE"],
    n_nonsig_no_overlap = tbl["FALSE", "FALSE"]
  )
  save_tsv(fisher_dt, file.path(opt$out_dir, "k119ub_overlap_fisher.tsv"))

  # DiffBind directional
  db <- fread(opt$k119ub_diffbind)
  db <- classify_diffbind(db, opt$fdr_threshold)
  cat("  K119ub DiffBind peaks:", nrow(db), "\n")
  cat("  Gained:", sum(db$direction == "gained"),
      " Lost:", sum(db$direction == "lost"),
      " Unchanged:", sum(db$direction == "unchanged"), "\n")

  # Use Peak coordinates for gene overlap
  db_gr <- GRanges(seqnames = db$Peak_Chr,
                   ranges = IRanges(start = db$Peak_Start + 1L, end = db$Peak_End))
  db_gr$direction <- db$direction

  gained_gr <- db_gr[db_gr$direction == "gained"]
  lost_gr <- db_gr[db_gr$direction == "lost"]

  k119ub_gained_all <- countOverlaps(gene_gr, gained_gr) > 0
  results$k119ub_gained <- k119ub_gained_all[match(results$gene_name, gene_bed$gene_name)]
  k119ub_lost_all <- countOverlaps(gene_gr, lost_gr) > 0
  results$k119ub_lost <- k119ub_lost_all[match(results$gene_name, gene_bed$gene_name)]
  results$k119ub_status <- ifelse(results$k119ub_gained, "gained",
                           ifelse(results$k119ub_lost, "lost", "unchanged"))

  # Wilcoxon tests
  wt_gv <- wilcox.test(results$edger_logFC[results$k119ub_status == "gained"],
                       results$edger_logFC[results$k119ub_status == "unchanged"])
  wt_lv <- wilcox.test(results$edger_logFC[results$k119ub_status == "lost"],
                       results$edger_logFC[results$k119ub_status == "unchanged"])
  cat(sprintf("  Wilcoxon gained vs unchanged: p=%.2e\n", wt_gv$p.value))
  cat(sprintf("  Wilcoxon lost vs unchanged: p=%.2e\n", wt_lv$p.value))

  wilcox_dt <- data.table(
    comparison = c("gained_vs_unchanged", "lost_vs_unchanged"),
    p_value = c(wt_gv$p.value, wt_lv$p.value),
    statistic = c(wt_gv$statistic, wt_lv$statistic),
    n_group1 = c(sum(results$k119ub_status == "gained"),
                 sum(results$k119ub_status == "lost")),
    n_group2 = c(sum(results$k119ub_status == "unchanged"),
                 sum(results$k119ub_status == "unchanged")),
    median_logFC_group1 = c(median(results$edger_logFC[results$k119ub_status == "gained"]),
                            median(results$edger_logFC[results$k119ub_status == "lost"])),
    median_logFC_group2 = c(median(results$edger_logFC[results$k119ub_status == "unchanged"]),
                            median(results$edger_logFC[results$k119ub_status == "unchanged"]))
  )
  save_tsv(wilcox_dt, file.path(opt$out_dir, "k119ub_direction_wilcoxon.tsv"))

  # Boxplot
  plot_df <- results[results$k119ub_status %in% c("gained", "unchanged", "lost"), ]
  plot_df$k119ub_status <- factor(plot_df$k119ub_status,
                                  levels = c("lost", "unchanged", "gained"))
  p <- ggplot(plot_df, aes(x = k119ub_status, y = edger_logFC, fill = k119ub_status)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    stat_summary(fun.data = function(x) data.frame(y = quantile(x, 0.75) + 1.5 * IQR(x), label = paste0("n=", format(length(x), big.mark = ","))), geom = "text", size = 2.8, vjust = -0.3) +
    stat_summary(fun = median, geom = "text", aes(label = sprintf("%.4f", after_stat(y))), vjust = 1.5, size = 2.5, colour = "grey20") +
    scale_fill_manual(values = c(lost = "#2166AC", unchanged = "grey70", gained = "#B2182B")) +
    labs(x = "K119ub Status in BAP1-KO", y = "mCH logFC (mut vs ctrl)",
         title = "mCH Change by K119ub Direction",
         subtitle = sprintf("Wilcoxon: gained vs unchanged p=%.1e, lost vs unchanged p=%.1e",
                            wt_gv$p.value, wt_lv$p.value)) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  save_plot(p, file.path(opt$plot_dir, "04_k119ub_direction_boxplot"))

  results
}

# ---------------------------------------------------------------------------
# 4b. MeCP2 binding
# ---------------------------------------------------------------------------

run_mecp2 <- function(results, gene_bed, gene_gr, opt) {
  cat("\n=== 4b. MeCP2 Binding ===\n")

  # Consensus overlap
  mecp2_gr <- load_bed3(opt$mecp2_consensus)
  cat("  MeCP2 consensus peaks:", length(mecp2_gr), "\n")

  overlap <- countOverlaps(gene_gr, mecp2_gr) > 0
  results$mecp2_overlap <- overlap[match(results$gene_name, gene_bed$gene_name)]

  sig <- !is.na(results$sig_fdr005) & results$sig_fdr005
  tbl <- table(significant = sig, mecp2 = results$mecp2_overlap)
  ft <- fisher.test(tbl)
  cat(sprintf("  Fisher's exact: OR=%.3f, p=%.2e\n", ft$estimate, ft$p.value))

  fisher_dt <- data.table(
    test = "mCH_sig_x_MeCP2_overlap",
    odds_ratio = as.numeric(ft$estimate),
    p_value = ft$p.value,
    ci_lower = ft$conf.int[1],
    ci_upper = ft$conf.int[2],
    n_sig_overlap = tbl["TRUE", "TRUE"],
    n_sig_no_overlap = tbl["TRUE", "FALSE"],
    n_nonsig_overlap = tbl["FALSE", "TRUE"],
    n_nonsig_no_overlap = tbl["FALSE", "FALSE"]
  )
  save_tsv(fisher_dt, file.path(opt$out_dir, "mecp2_overlap_fisher.tsv"))

  # DiffBind directional
  db <- fread(opt$mecp2_diffbind)
  db <- classify_diffbind(db, opt$fdr_threshold)
  cat("  MeCP2 DiffBind peaks:", nrow(db), "\n")
  cat("  Gained:", sum(db$direction == "gained"),
      " Lost:", sum(db$direction == "lost"),
      " Unchanged:", sum(db$direction == "unchanged"), "\n")

  # MeCP2 DiffBind uses seqnames/start/end columns
  db_gr <- GRanges(seqnames = db$seqnames,
                   ranges = IRanges(start = db$start + 1L, end = db$end))
  db_gr$direction <- db$direction

  gained_gr <- db_gr[db_gr$direction == "gained"]
  lost_gr <- db_gr[db_gr$direction == "lost"]

  mecp2_gained_all <- countOverlaps(gene_gr, gained_gr) > 0
  results$mecp2_gained <- mecp2_gained_all[match(results$gene_name, gene_bed$gene_name)]
  mecp2_lost_all <- countOverlaps(gene_gr, lost_gr) > 0
  results$mecp2_lost <- mecp2_lost_all[match(results$gene_name, gene_bed$gene_name)]
  results$mecp2_status <- ifelse(results$mecp2_gained, "gained",
                          ifelse(results$mecp2_lost, "lost", "unchanged"))

  # Wilcoxon tests
  wt_gv <- wilcox.test(results$edger_logFC[results$mecp2_status == "gained"],
                       results$edger_logFC[results$mecp2_status == "unchanged"])
  wt_lv <- wilcox.test(results$edger_logFC[results$mecp2_status == "lost"],
                       results$edger_logFC[results$mecp2_status == "unchanged"])
  cat(sprintf("  Wilcoxon gained vs unchanged: p=%.2e\n", wt_gv$p.value))
  cat(sprintf("  Wilcoxon lost vs unchanged: p=%.2e\n", wt_lv$p.value))

  wilcox_dt <- data.table(
    comparison = c("gained_vs_unchanged", "lost_vs_unchanged"),
    p_value = c(wt_gv$p.value, wt_lv$p.value),
    statistic = c(wt_gv$statistic, wt_lv$statistic),
    n_group1 = c(sum(results$mecp2_status == "gained"),
                 sum(results$mecp2_status == "lost")),
    n_group2 = c(sum(results$mecp2_status == "unchanged"),
                 sum(results$mecp2_status == "unchanged")),
    median_logFC_group1 = c(median(results$edger_logFC[results$mecp2_status == "gained"]),
                            median(results$edger_logFC[results$mecp2_status == "lost"])),
    median_logFC_group2 = c(median(results$edger_logFC[results$mecp2_status == "unchanged"]),
                            median(results$edger_logFC[results$mecp2_status == "unchanged"]))
  )
  save_tsv(wilcox_dt, file.path(opt$out_dir, "mecp2_direction_wilcoxon.tsv"))

  # Boxplot
  plot_df <- results[results$mecp2_status %in% c("gained", "unchanged", "lost"), ]
  plot_df$mecp2_status <- factor(plot_df$mecp2_status,
                                 levels = c("lost", "unchanged", "gained"))
  p <- ggplot(plot_df, aes(x = mecp2_status, y = edger_logFC, fill = mecp2_status)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    stat_summary(fun.data = function(x) data.frame(y = quantile(x, 0.75) + 1.5 * IQR(x), label = paste0("n=", format(length(x), big.mark = ","))), geom = "text", size = 2.8, vjust = -0.3) +
    stat_summary(fun = median, geom = "text", aes(label = sprintf("%.4f", after_stat(y))), vjust = 1.5, size = 2.5, colour = "grey20") +
    scale_fill_manual(values = c(lost = "#2166AC", unchanged = "grey70", gained = "#B2182B")) +
    labs(x = "MeCP2 Status in BAP1-KO", y = "mCH logFC (mut vs ctrl)",
         title = "mCH Change by MeCP2 Direction",
         subtitle = sprintf("Wilcoxon: gained vs unchanged p=%.1e, lost vs unchanged p=%.1e",
                            wt_gv$p.value, wt_lv$p.value)) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  save_plot(p, file.path(opt$plot_dir, "04_mecp2_direction_boxplot"))

  results
}

# ---------------------------------------------------------------------------
# 4c. Neuronal gene set and gene length
# ---------------------------------------------------------------------------

run_neuronal <- function(results, opt) {
  cat("\n=== 4c. Neuronal Gene Set / Gene Length ===\n")

  # Load neuronal genes
  neuronal <- fread(opt$neuronal_genes)
  neuronal_genes <- neuronal$gene
  cat("  Neuronal gene set:", length(neuronal_genes), "genes\n")

  results$neuronal <- results$gene_name %in% neuronal_genes
  n_match <- sum(results$neuronal)
  cat("  Matched in testable genes:", n_match, "/", length(neuronal_genes), "\n")

  # --- Wilcoxon rank-sum ---
  wt <- wilcox.test(results$edger_logFC[results$neuronal],
                    results$edger_logFC[!results$neuronal])
  cat(sprintf("  Wilcoxon neuronal vs non-neuronal: p=%.2e\n", wt$p.value))

  wilcox_dt <- data.table(
    comparison = "neuronal_vs_non_neuronal",
    p_value = wt$p.value,
    statistic = as.numeric(wt$statistic),
    n_neuronal = sum(results$neuronal),
    n_non_neuronal = sum(!results$neuronal),
    median_logFC_neuronal = median(results$edger_logFC[results$neuronal]),
    median_logFC_non_neuronal = median(results$edger_logFC[!results$neuronal])
  )
  save_tsv(wilcox_dt, file.path(opt$out_dir, "neuronal_wilcoxon.tsv"))

  # Boxplot
  p_box <- ggplot(results, aes(x = neuronal, y = edger_logFC, fill = neuronal)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    stat_summary(fun.data = function(x) data.frame(y = quantile(x, 0.75) + 1.5 * IQR(x), label = paste0("n=", format(length(x), big.mark = ","))), geom = "text", size = 2.8, vjust = -0.3) +
    stat_summary(fun = median, geom = "text", aes(label = sprintf("%.4f", after_stat(y))), vjust = 1.5, size = 2.5, colour = "grey20") +
    scale_fill_manual(values = c("TRUE" = "#E66101", "FALSE" = "grey70")) +
    scale_x_discrete(labels = c("Non-neuronal", "Neuronal")) +
    labs(x = NULL, y = "mCH logFC (mut vs ctrl)",
         title = "mCH Change: Neuronal vs Non-neuronal Genes",
         subtitle = sprintf("Wilcoxon p = %.1e", wt$p.value)) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  save_plot(p_box, file.path(opt$plot_dir, "04_neuronal_wilcoxon_boxplot"))

  # --- Linear interaction model ---
  cat("\n  Gene-length interaction model:\n")
  results$log_length <- log(results$gene_length)
  fit <- lm(edger_logFC ~ neuronal * log_length, data = results)
  coefs <- summary(fit)$coefficients
  cat("  Coefficients:\n")
  print(round(coefs, 6))

  interaction_dt <- as.data.table(coefs, keep.rownames = "term")
  setnames(interaction_dt, c("term", "estimate", "std_error", "t_value", "p_value"))
  save_tsv(interaction_dt, file.path(opt$out_dir, "genelength_interaction_model.tsv"))

  # --- Quintile stratification ---
  results$length_quintile <- cut(results$gene_length,
                                 breaks = quantile(results$gene_length, probs = 0:5/5),
                                 include.lowest = TRUE, labels = paste0("Q", 1:5))

  p_quint <- ggplot(results[!is.na(results$length_quintile), ],
                    aes(x = length_quintile, y = edger_logFC, fill = neuronal)) +
    geom_boxplot(outlier.size = 0.2, alpha = 0.7, position = position_dodge(0.8)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    stat_summary(fun.data = function(x) data.frame(y = quantile(x, 0.75) + 1.5 * IQR(x), label = paste0("n=", format(length(x), big.mark = ","))), geom = "text", size = 2.2, vjust = -0.3, position = position_dodge(0.8)) +
    scale_fill_manual(values = c("TRUE" = "#E66101", "FALSE" = "grey70"),
                      labels = c("Non-neuronal", "Neuronal")) +
    labs(x = "Gene Length Quintile (Q1=shortest, Q5=longest)",
         y = "mCH logFC (mut vs ctrl)",
         title = "mCH Change by Gene Length and Neuronal Status",
         fill = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  save_plot(p_quint, file.path(opt$plot_dir, "04_genelength_quintile_boxplot"), w = 9, h = 6)

  # --- Supplementary JT trend test ---
  quintile_effects <- results[!is.na(results$length_quintile),
                              .(effect = median(edger_logFC[neuronal]) -
                                         median(edger_logFC[!neuronal]),
                                n_neuronal = sum(neuronal),
                                n_non = sum(!neuronal)),
                              by = length_quintile]
  setorder(quintile_effects, length_quintile)
  cat("\n  Quintile effects (neuronal - non-neuronal median logFC):\n")
  print(quintile_effects)

  if (nrow(quintile_effects) >= 3) {
    jt <- jonckheere.test(quintile_effects$effect,
                          as.numeric(factor(quintile_effects$length_quintile)),
                          alternative = "increasing")
    cat(sprintf("  JT trend test: statistic=%.1f, p=%.4f\n",
                jt$statistic, jt$p.value))
  }

  # --- fgsea for custom neuronal gene set ---
  cat("\n  fgsea: neuronal gene set enrichment\n")
  stats_vec <- setNames(results$edger_logFC, results$gene_name)
  stats_vec <- sort(stats_vec, decreasing = TRUE)
  stats_vec <- stats_vec[is.finite(stats_vec)]

  fgsea_res <- fgsea(pathways = list(neuronal = neuronal_genes),
                     stats = stats_vec, minSize = 15)
  cat(sprintf("  NES=%.3f, p=%.2e, padj=%.2e\n",
              fgsea_res$NES, fgsea_res$pval, fgsea_res$padj))

  fgsea_dt <- as.data.table(fgsea_res)
  fgsea_dt$leadingEdge <- sapply(fgsea_dt$leadingEdge, paste, collapse = ";")
  save_tsv(fgsea_dt, file.path(opt$out_dir, "fgsea_neuronal.tsv"))

  p_enrich <- plotEnrichment(neuronal_genes, stats_vec) +
    labs(title = "fgsea: Neuronal Gene Set",
         subtitle = sprintf("NES=%.3f, p=%.2e", fgsea_res$NES, fgsea_res$pval)) +
    theme_bw(base_size = 12)
  save_plot(p_enrich, file.path(opt$plot_dir, "04_fgsea_neuronal_enrichment"))

  # --- Supplementary GO GSEA ---
  cat("\n  GO GSEA (biological process):\n")
  id_map <- tryCatch(
    bitr(names(stats_vec), fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Mm.eg.db),
    error = function(e) { cat("  bitr error:", conditionMessage(e), "\n"); NULL }
  )

  if (!is.null(id_map) && nrow(id_map) > 0) {
    gene_list <- stats_vec[id_map$SYMBOL]
    names(gene_list) <- id_map$ENTREZID
    gene_list <- sort(gene_list, decreasing = TRUE)

    gsea_go <- tryCatch(
      gseGO(geneList = gene_list, ont = "BP", OrgDb = org.Mm.eg.db,
            minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1,
            verbose = FALSE),
      error = function(e) { cat("  gseGO error:", conditionMessage(e), "\n"); NULL }
    )

    if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
      go_df <- as.data.table(as.data.frame(gsea_go))
      go_df$core_enrichment <- substr(go_df$core_enrichment, 1, 200)
      save_tsv(go_df, file.path(opt$out_dir, "gsea_go_results.tsv"))
      cat(sprintf("  GO terms at p<0.05: %d\n", sum(go_df$pvalue < 0.05)))
    } else {
      cat("  No GO GSEA results.\n")
      save_tsv(data.table(note = "No significant GO terms"),
               file.path(opt$out_dir, "gsea_go_results.tsv"))
    }
  }

  results
}

# ---------------------------------------------------------------------------
# 4d. CG x mCH cross-modality
# ---------------------------------------------------------------------------

run_cg_cross <- function(results, opt) {
  cat("\n=== 4d. CG x mCH Cross-Modality ===\n")

  # Load CG mC DMR
  mc <- fread(opt$mc_dmr)
  if (!grepl("^chr", mc[[1]][1])) mc[[1]] <- paste0("chr", mc[[1]])
  setnames(mc, c(names(mc)[1], "Name"), c("chr", "gene_name"), skip_absent = TRUE)
  if (!"gene_name" %in% names(mc) && ncol(mc) >= 14) {
    setnames(mc, names(mc)[14], "gene_name")
  }
  cat("  CG mC DMR genes:", nrow(mc), "\n")

  # Load CG hmC DMR
  hmc <- fread(opt$hmc_dmr)
  if (!grepl("^chr", hmc[[1]][1])) hmc[[1]] <- paste0("chr", hmc[[1]])
  if (!"gene_name" %in% names(hmc) && ncol(hmc) >= 14) {
    setnames(hmc, names(hmc)[14], "gene_name")
  }
  cat("  CG hmC DMR genes:", nrow(hmc), "\n")

  # Join with mCH results by gene name
  mc_sub <- mc[, .(gene_name, mc_mod_diff = mod_difference, mc_qvalue = dmr_qvalue)]
  hmc_sub <- hmc[, .(gene_name, hmc_mod_diff = mod_difference, hmc_qvalue = dmr_qvalue)]

  joined <- merge(results[, .(gene_name, edger_logFC, edger_fdr, mch_diff, gene_length)],
                  mc_sub, by = "gene_name", all.x = TRUE)
  joined <- merge(joined, hmc_sub, by = "gene_name", all.x = TRUE)

  has_mc <- !is.na(joined$mc_mod_diff)
  has_hmc <- !is.na(joined$hmc_mod_diff)
  cat("  Genes with CG mC data:", sum(has_mc), "\n")
  cat("  Genes with CG hmC data:", sum(has_hmc), "\n")

  # --- Spearman correlation ---
  cor_test <- cor.test(joined$mc_mod_diff[has_mc], joined$edger_logFC[has_mc],
                       method = "spearman")
  cat(sprintf("  Spearman (CG mC vs mCH): rho=%.3f, p=%.2e\n",
              cor_test$estimate, cor_test$p.value))

  cor_dt <- data.table(
    comparison = "CG_mC_vs_mCH",
    spearman_rho = as.numeric(cor_test$estimate),
    p_value = cor_test$p.value,
    n_genes = sum(has_mc)
  )
  save_tsv(cor_dt, file.path(opt$out_dir, "cg_mch_correlation.tsv"))

  # --- Quadrant analysis ---
  jmc <- joined[has_mc]
  jmc$mc_dir <- ifelse(jmc$mc_mod_diff > 0, "CG_mC_up", "CG_mC_down")
  jmc$mch_dir <- ifelse(jmc$edger_logFC > 0, "mCH_up", "mCH_down")
  jmc$quadrant <- paste(jmc$mc_dir, jmc$mch_dir, sep = " / ")

  quad_counts <- as.data.table(table(jmc$quadrant))
  setnames(quad_counts, c("quadrant", "n_genes"))
  save_tsv(quad_counts, file.path(opt$out_dir, "cg_mch_quadrant_counts.tsv"))
  cat("  Quadrant counts:\n")
  print(quad_counts)

  # Scatter plot
  p_scatter <- ggplot(jmc, aes(x = mc_mod_diff, y = edger_logFC)) +
    geom_point(size = 0.3, alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.5) +
    labs(x = "CG mC difference (mut - ctrl)",
         y = "mCH logFC (mut vs ctrl)",
         title = "CG Methylation vs Non-CG Methylation",
         subtitle = sprintf("Spearman rho=%.3f, p=%.2e",
                            cor_test$estimate, cor_test$p.value)) +
    theme_bw(base_size = 12)
  save_plot(p_scatter, file.path(opt$plot_dir, "04_cg_mch_scatter"))

  # --- Fisher's exact: CG-sig x mCH-sig ---
  mc_sig <- jmc$mc_qvalue < 0.05
  mch_sig <- jmc$edger_fdr < 0.05
  tbl <- table(mCH_sig = mch_sig, CG_mC_sig = mc_sig)
  ft <- fisher.test(tbl)
  cat(sprintf("  Fisher CG-sig x mCH-sig: OR=%.3f, p=%.2e\n",
              ft$estimate, ft$p.value))

  # --- Three-way: CG mC x CG hmC x mCH ---
  cat("\n  Three-way (mC / hmC / mCH):\n")
  j3 <- joined[has_mc & has_hmc]
  j3$coordinated <- j3$mc_mod_diff > 0 & j3$hmc_mod_diff < 0
  j3$mc_sig <- j3$mc_qvalue < 0.05
  j3$hmc_sig <- j3$hmc_qvalue < 0.05
  j3$concordant <- j3$mc_sig & j3$hmc_sig & j3$coordinated

  n_concordant <- sum(j3$concordant)
  cat(sprintf("  Concordant (mC up + hmC down, both sig): %d\n", n_concordant))

  if (n_concordant > 10) {
    wt_3way <- wilcox.test(j3$edger_logFC[j3$concordant],
                           j3$edger_logFC[!j3$concordant])
    cat(sprintf("  mCH logFC concordant vs rest: p=%.2e\n", wt_3way$p.value))
    cat(sprintf("    Concordant median logFC: %.4f\n",
                median(j3$edger_logFC[j3$concordant])))
    cat(sprintf("    Rest median logFC: %.4f\n",
                median(j3$edger_logFC[!j3$concordant])))
  }

  threeway_dt <- data.table(
    category = c("concordant_mCup_hmCdown", "other"),
    n_genes = c(n_concordant, nrow(j3) - n_concordant),
    median_mch_logFC = c(median(j3$edger_logFC[j3$concordant]),
                         median(j3$edger_logFC[!j3$concordant]))
  )
  save_tsv(threeway_dt, file.path(opt$out_dir, "cg_hmc_mch_threeway.tsv"))

  # Three-way scatter
  if (nrow(j3) > 0) {
    p3 <- ggplot(j3, aes(x = mc_mod_diff, y = edger_logFC,
                         color = hmc_mod_diff)) +
      geom_point(size = 0.5, alpha = 0.5) +
      scale_color_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B",
                            midpoint = 0, name = "hmC diff") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      labs(x = "CG mC difference", y = "mCH logFC",
           title = "Three-way: CG mC x CG hmC x mCH") +
      theme_bw(base_size = 12)
    save_plot(p3, file.path(opt$plot_dir, "04_cg_hmc_mch_threeway_scatter"))
  }
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

write_integration_summary <- function(opt) {
  result_files <- list.files(opt$out_dir, pattern = "\\.tsv$", full.names = TRUE)

  sink(file.path(opt$out_dir, "integration_summary.txt"))
  cat("=== mCH Integration Analysis Summary ===\n\n")
  cat("Output files:\n")
  for (f in result_files) {
    cat("  ", basename(f), "\n")
  }
  cat("\nPlot files:\n")
  plot_files <- list.files(opt$plot_dir, pattern = "\\.(pdf|png)$")
  for (f in plot_files) {
    cat("  ", f, "\n")
  }
  sink()
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  opt <- parse_cli_args()
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(opt$plot_dir, recursive = TRUE, showWarnings = FALSE)

  cat("=== mCH Integration Analyses ===\n\n")

  # Load differential results
  results <- fread(opt$diff_results)
  cat("Differential results loaded:", nrow(results), "genes\n")

  # Load gene body BED for GRanges
  gene_bed <- fread(opt$gene_bed, header = FALSE,
                    col.names = c("chr", "start", "end", "gene_name", "gene_id",
                                  "strand", "gene_type", "gene_length"))
  gene_gr <- gene_bed_to_gr(gene_bed)

  # Run sub-analyses
  results <- run_k119ub(results, gene_bed, gene_gr, opt)
  results <- run_mecp2(results, gene_bed, gene_gr, opt)
  results <- run_neuronal(results, opt)
  run_cg_cross(results, opt)

  write_integration_summary(opt)

  cat("\n=== Integration Complete ===\n")
}

main()
