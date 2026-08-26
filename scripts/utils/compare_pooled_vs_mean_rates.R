# scripts/utils/compare_pooled_vs_mean_rates.R
#
# Compare two gene-level mCH rate estimators:
#   1. rowMeans (current step 03): mean of per-sample rates, equal weight
#   2. Pooled (metagene-style): sum(meth) / sum(coverage), coverage-weighted
#
# Reads the per-sample aggregated counts and the differential results.
# Writes a comparison TSV and prints summary statistics.

agg_dir <- "results/02_aggregate/aggregated"
diff_path <- "results/03_differential/mch_differential_results.tsv"
out_path <- "results/03_differential/pooled_vs_mean_rate_comparison.tsv"

ctrl_ids <- c("ctrl_M1", "ctrl_M2", "ctrl_F1", "ctrl_F2")
mut_ids <- c("mut_M1", "mut_M2", "mut_F1", "mut_F2")

cat("Loading per-sample aggregated counts...\n")
load_sample <- function(sid) {
  path <- file.path(agg_dir, paste0(sid, "_genebody_mch.tsv"))
  d <- read.delim(path, stringsAsFactors = FALSE)
  data.frame(
    gene_name = d$gene_name,
    coverage = d$total_coverage,
    methylated = d$meth_reads,
    stringsAsFactors = FALSE
  )
}

ctrl_data <- lapply(ctrl_ids, load_sample)
mut_data <- lapply(mut_ids, load_sample)

genes <- ctrl_data[[1]]$gene_name
for (i in seq_along(ctrl_data)) {
  stopifnot(identical(ctrl_data[[i]]$gene_name, genes))
}
for (i in seq_along(mut_data)) {
  stopifnot(identical(mut_data[[i]]$gene_name, genes))
}

ctrl_cov <- do.call(cbind, lapply(ctrl_data, `[[`, "coverage"))
ctrl_met <- do.call(cbind, lapply(ctrl_data, `[[`, "methylated"))
mut_cov <- do.call(cbind, lapply(mut_data, `[[`, "coverage"))
mut_met <- do.call(cbind, lapply(mut_data, `[[`, "methylated"))

pooled_ctrl <- rowSums(ctrl_met) / rowSums(ctrl_cov)
pooled_mut <- rowSums(mut_met) / rowSums(mut_cov)
pooled_diff <- pooled_mut - pooled_ctrl

cat("Loading differential results for comparison...\n")
diff <- read.delim(diff_path, stringsAsFactors = FALSE)

m <- match(diff$gene_name, genes)
missing <- sum(is.na(m))
if (missing > 0) {
  cat(sprintf("  %d genes in diff results not found in aggregated data (filtered by coverage)\n", missing))
}
keep <- !is.na(m)
diff <- diff[keep, ]
m <- m[keep]

out <- data.frame(
  gene_name = diff$gene_name,
  mean_ctrl = diff$mch_ctrl,
  mean_mut = diff$mch_mut,
  mean_diff = diff$mch_diff,
  pooled_ctrl = pooled_ctrl[m],
  pooled_mut = pooled_mut[m],
  pooled_diff = pooled_diff[m],
  ctrl_delta = pooled_ctrl[m] - diff$mch_ctrl,
  mut_delta = pooled_mut[m] - diff$mch_mut,
  diff_delta = pooled_diff[m] - diff$mch_diff,
  edger_logFC = diff$edger_logFC,
  edger_fdr = diff$edger_fdr,
  sig_fdr005 = diff$sig_fdr005,
  stringsAsFactors = FALSE
)

write.table(out, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote %d genes to %s\n\n", nrow(out), out_path))

cat("=== Summary of pooled - mean differences ===\n")
for (col in c("ctrl_delta", "mut_delta", "diff_delta")) {
  v <- out[[col]]
  cat(sprintf("  %s:  median=%.2e  mean=%.2e  sd=%.2e  max|delta|=%.2e\n",
              col, median(v), mean(v), sd(v), max(abs(v))))
}

cat(sprintf("\nSpearman(mean_diff, pooled_diff): %.6f\n",
            cor(out$mean_diff, out$pooled_diff, method = "spearman")))
cat(sprintf("Spearman(mean_diff, edger_logFC): %.6f\n",
            cor(out$mean_diff, out$edger_logFC, method = "spearman")))
cat(sprintf("Spearman(pooled_diff, edger_logFC): %.6f\n",
            cor(out$pooled_diff, out$edger_logFC, method = "spearman")))

sig <- out[out$sig_fdr005 == TRUE, ]
mean_flipped <- sum(sig$edger_logFC < 0 & sig$mean_diff > 0) +
                sum(sig$edger_logFC > 0 & sig$mean_diff < 0)
pooled_flipped <- sum(sig$edger_logFC < 0 & sig$pooled_diff > 0) +
                  sum(sig$edger_logFC > 0 & sig$pooled_diff < 0)

cat(sprintf("\nAmong %d FDR<0.05 genes:\n", nrow(sig)))
cat(sprintf("  Direction disagreements with edgeR (mean_diff):   %d (%.1f%%)\n",
            mean_flipped, 100 * mean_flipped / nrow(sig)))
cat(sprintf("  Direction disagreements with edgeR (pooled_diff): %d (%.1f%%)\n",
            pooled_flipped, 100 * pooled_flipped / nrow(sig)))
