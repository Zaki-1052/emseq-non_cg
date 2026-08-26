# scripts/utils/metagene_equalweight.R
#
# Recompute the metagene profile using equal-weight-per-sample rates instead
# of pooled (coverage-weighted) rates. Compares against the existing pooled
# metagene to quantify how much the deep outlier samples shift the mutant line.

library(data.table)
library(ggplot2)

feature_dir <- "results/02b_features"
out_dir <- "results/sections/50_features/metagene_correction"

n_bins <- 50L
min_coverage <- 20L

samples <- data.table(
  sample_id = c("ctrl_M1", "ctrl_M2", "ctrl_F1", "ctrl_F2",
                "mut_M1", "mut_M2", "mut_F1", "mut_F2"),
  genotype = c(rep("ctrl", 4), rep("mut", 4))
)

body_features <- c("Exon", "Intron")

cat("Loading per-sample feature tables...\n")
all_samples <- list()
for (i in seq_len(nrow(samples))) {
  sid <- samples$sample_id[i]
  path <- file.path(feature_dir, paste0(sid, "_feature_mch.tsv"))
  dt <- fread(path)
  dt <- dt[feature_type %in% body_features]
  dt[, interval_key := paste(gene_name, feature_type, chr, start, end, sep = "|")]
  all_samples[[sid]] <- dt
  cat(sprintf("  %-10s %s intervals  %s total coverage\n",
              sid, format(nrow(dt), big.mark = ","),
              format(sum(dt$total_coverage), big.mark = ",")))
}

ref_keys <- all_samples[[1]]$interval_key
for (sid in samples$sample_id[-1]) {
  stopifnot(identical(all_samples[[sid]]$interval_key, ref_keys))
}

ref <- all_samples[[1]][, .(interval_key, gene_name, feature_type, chr, start, end,
                            strand, feature_length)]

cat("\nBuilding gene body coordinates...\n")
spans <- ref[, .(
  body_start = min(start),
  body_end = max(end),
  gene_strand = strand[1]
), by = gene_name]
spans <- spans[body_end > body_start]
cat(sprintf("  %s genes with valid body spans\n", format(nrow(spans), big.mark = ",")))

ref <- merge(ref, spans[, .(gene_name, body_start, body_end, gene_strand)],
             by = "gene_name")
ref[, body_length := body_end - body_start]
ref[, raw_start := (start - body_start) / body_length]
ref[, raw_end := (end - body_start) / body_length]
ref[, oriented_start := ifelse(gene_strand == "-", 1 - raw_end, raw_start)]
ref[, oriented_end := ifelse(gene_strand == "-", 1 - raw_start, raw_end)]
ref[, rel_start := pmax(0, pmin(1, oriented_start))]
ref[, rel_end := pmax(0, pmin(1, oriented_end))]
ref[, span := rel_end - rel_start]
ref <- ref[span > 0]

cat(sprintf("  %s intervals after filtering\n", format(nrow(ref), big.mark = ",")))

ref[, bin_first := as.integer(pmin(n_bins, floor(rel_start * n_bins) + 1))]
ref[, bin_last := as.integer(pmax(1, pmin(n_bins, ceiling(rel_end * n_bins))))]
ref[, n_bins_spanned := bin_last - bin_first + 1L]

surviving_keys <- ref$interval_key

cat("Expanding intervals into bins...\n")
slim <- ref[, .(interval_key, gene_name, rel_start, rel_end, span,
                bin_first, n_bins_spanned)]
expanded <- slim[rep(seq_len(nrow(slim)), slim$n_bins_spanned)]
expanded[, bin := bin_first + sequence(slim$n_bins_spanned) - 1L]
expanded[, bin_low := (bin - 1) / n_bins]
expanded[, bin_high := bin / n_bins]
expanded[, share := pmax(0, pmin(rel_end, bin_high) - pmax(rel_start, bin_low)) / span]

cat("Computing per-sample rates per gene x bin...\n")
gene_bin_rates <- list()

for (i in seq_len(nrow(samples))) {
  sid <- samples$sample_id[i]
  geno <- samples$genotype[i]
  dt <- all_samples[[sid]]

  dt_sub <- dt[match(surviving_keys, interval_key),
               .(interval_key, total_coverage, methylated_count)]

  exp_with_counts <- merge(expanded, dt_sub, by = "interval_key", sort = FALSE)

  gb <- exp_with_counts[, .(
    coverage = sum(total_coverage * share),
    methylated = sum(methylated_count * share)
  ), by = .(gene_name, bin)]

  gb <- gb[coverage >= min_coverage]
  gb[, rate := methylated / coverage]
  gb[, sample_id := sid]
  gb[, genotype := geno]
  gene_bin_rates[[sid]] <- gb
}

all_rates <- rbindlist(gene_bin_rates)
cat(sprintf("  %s gene x bin x sample records\n", format(nrow(all_rates), big.mark = ",")))

cat("Computing equal-weight and pooled profiles...\n")

sample_bin <- all_rates[, .(mean_rate = mean(rate)), by = .(genotype, sample_id, bin)]

equalweight <- sample_bin[, .(
  ew_rate = mean(mean_rate),
  ew_se = sd(mean_rate) / sqrt(.N)
), by = .(genotype, bin)]

ew_ctrl <- equalweight[genotype == "ctrl", .(bin, ew_ctrl = ew_rate, ew_ctrl_se = ew_se)]
ew_mut <- equalweight[genotype == "mut", .(bin, ew_mut = ew_rate, ew_mut_se = ew_se)]
ew_profile <- merge(ew_ctrl, ew_mut, by = "bin")
ew_profile[, ew_delta := ew_mut - ew_ctrl]

pooled <- all_rates[, .(
  total_cov = sum(coverage),
  total_met = sum(methylated)
), by = .(genotype, bin)]
pooled[, pooled_rate := total_met / total_cov]

pool_ctrl <- pooled[genotype == "ctrl", .(bin, pooled_ctrl = pooled_rate)]
pool_mut <- pooled[genotype == "mut", .(bin, pooled_mut = pooled_rate)]
pool_profile <- merge(pool_ctrl, pool_mut, by = "bin")
pool_profile[, pooled_delta := pooled_mut - pooled_ctrl]

profile <- merge(ew_profile, pool_profile, by = "bin")
profile[, mut_shift_pp := (pooled_mut - ew_mut) * 100]
profile[, ctrl_shift_pp := (pooled_ctrl - ew_ctrl) * 100]
profile[, delta_shift_pp := (pooled_delta - ew_delta) * 100]

out_tsv <- file.path(out_dir, "metagene_equalweight_vs_pooled.tsv")
write.table(profile, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Shift summary (pooled minus equal-weight, percentage points) ===\n")
cat(sprintf("  Ctrl line:  median=%.4f pp  mean=%.4f pp\n",
            median(profile$ctrl_shift_pp), mean(profile$ctrl_shift_pp)))
cat(sprintf("  Mut line:   median=%.4f pp  mean=%.4f pp\n",
            median(profile$mut_shift_pp), mean(profile$mut_shift_pp)))
cat(sprintf("  Delta:      median=%.4f pp  mean=%.4f pp\n",
            median(profile$delta_shift_pp), mean(profile$delta_shift_pp)))

cat(sprintf("\n  Equal-weight ctrl (mean across bins): %.5f (%.3f%%)\n",
            mean(ew_profile$ew_ctrl), mean(ew_profile$ew_ctrl) * 100))
cat(sprintf("  Equal-weight mut  (mean across bins): %.5f (%.3f%%)\n",
            mean(ew_profile$ew_mut), mean(ew_profile$ew_mut) * 100))
cat(sprintf("  Equal-weight delta:                   %.6f (%.4f%%)\n",
            mean(ew_profile$ew_delta), mean(ew_profile$ew_delta) * 100))
cat(sprintf("  Pooled ctrl:                          %.5f (%.3f%%)\n",
            mean(pool_profile$pooled_ctrl), mean(pool_profile$pooled_ctrl) * 100))
cat(sprintf("  Pooled mut:                           %.5f (%.3f%%)\n",
            mean(pool_profile$pooled_mut), mean(pool_profile$pooled_mut) * 100))
cat(sprintf("  Pooled delta:                         %.6f (%.4f%%)\n",
            mean(pool_profile$pooled_delta), mean(pool_profile$pooled_delta) * 100))

cat("\nPlotting...\n")

plot_data <- rbind(
  data.table(bin = profile$bin, rate = profile$ew_ctrl * 100,
             genotype = "Control", method = "Equal-weight"),
  data.table(bin = profile$bin, rate = profile$ew_mut * 100,
             genotype = "Mutant", method = "Equal-weight"),
  data.table(bin = profile$bin, rate = profile$pooled_ctrl * 100,
             genotype = "Control", method = "Pooled"),
  data.table(bin = profile$bin, rate = profile$pooled_mut * 100,
             genotype = "Mutant", method = "Pooled")
)

p <- ggplot(plot_data, aes(x = bin, y = rate, color = genotype, linetype = method)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("Control" = "#3A7CA5", "Mutant" = "#D4713B")) +
  scale_linetype_manual(values = c("Equal-weight" = "solid", "Pooled" = "dashed")) +
  labs(x = "Gene body position (5' to 3')",
       y = "mCH rate (%)",
       title = "Metagene mCH profile: equal-weight vs pooled",
       subtitle = "Dashed = pooled (coverage-weighted), Solid = equal-weight per sample",
       color = "Genotype", linetype = "Method") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

plot_path <- file.path(out_dir, "metagene_equalweight_vs_pooled.png")
ggsave(plot_path, p, width = 8, height = 5, dpi = 150)
cat(sprintf("Saved plot to %s\n", plot_path))

p_delta <- ggplot() +
  geom_line(data = data.table(bin = profile$bin, delta = profile$ew_delta * 100,
                              method = "Equal-weight"),
            aes(x = bin, y = delta, linetype = method), color = "#2E6B8A", linewidth = 0.8) +
  geom_line(data = data.table(bin = profile$bin, delta = profile$pooled_delta * 100,
                              method = "Pooled"),
            aes(x = bin, y = delta, linetype = method), color = "#2E6B8A", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  scale_linetype_manual(values = c("Equal-weight" = "solid", "Pooled" = "dashed")) +
  labs(x = "Gene body position (5' to 3')",
       y = "Delta mCH (mut - ctrl, percentage points)",
       title = "Metagene delta: equal-weight vs pooled",
       linetype = "Method") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

delta_path <- file.path(out_dir, "metagene_delta_equalweight_vs_pooled.png")
ggsave(delta_path, p_delta, width = 8, height = 4, dpi = 150)
cat(sprintf("Saved delta plot to %s\n", delta_path))

# =============================================================================
# LAMBDA-CORRECTED METAGENE
# =============================================================================

cat("\n=== Lambda-corrected metagene ===\n")

lambda_dir <- "results/02_aggregate/aggregated"
lambda_rates <- sapply(samples$sample_id, function(sid) {
  path <- file.path(lambda_dir, paste0(sid, "_lambda.tsv"))
  d <- read.delim(path, stringsAsFactors = FALSE)
  d$lambda_ch_rate[1]
})
names(lambda_rates) <- samples$sample_id

cat("Per-sample lambda rates:\n")
for (sid in samples$sample_id) {
  cat(sprintf("  %-12s %.4f%%\n", sid, lambda_rates[sid] * 100))
}
cat(sprintf("  Mean ctrl: %.4f%%   Mean mut: %.4f%%   Difference: %.4f pp\n",
            mean(lambda_rates[samples$sample_id[samples$genotype == "ctrl"]]) * 100,
            mean(lambda_rates[samples$sample_id[samples$genotype == "mut"]]) * 100,
            (mean(lambda_rates[samples$sample_id[samples$genotype == "mut"]]) -
             mean(lambda_rates[samples$sample_id[samples$genotype == "ctrl"]])) * 100))

all_rates[, lambda := lambda_rates[sample_id]]
all_rates[, rate_corrected := pmax(0, rate - lambda)]

sample_bin_corr <- all_rates[, .(
  mean_rate = mean(rate_corrected)
), by = .(genotype, sample_id, bin)]

ew_corr <- sample_bin_corr[, .(
  ew_rate = mean(mean_rate),
  ew_se = sd(mean_rate) / sqrt(.N)
), by = .(genotype, bin)]

ew_corr_ctrl <- ew_corr[genotype == "ctrl", .(bin, corr_ctrl = ew_rate, corr_ctrl_se = ew_se)]
ew_corr_mut <- ew_corr[genotype == "mut", .(bin, corr_mut = ew_rate, corr_mut_se = ew_se)]
corr_profile <- merge(ew_corr_ctrl, ew_corr_mut, by = "bin")
corr_profile[, corr_delta := corr_mut - corr_ctrl]

cat(sprintf("\n  Lambda-corrected equal-weight ctrl: %.5f (%.3f%%)\n",
            mean(corr_profile$corr_ctrl), mean(corr_profile$corr_ctrl) * 100))
cat(sprintf("  Lambda-corrected equal-weight mut:  %.5f (%.3f%%)\n",
            mean(corr_profile$corr_mut), mean(corr_profile$corr_mut) * 100))
cat(sprintf("  Lambda-corrected delta:             %.6f (%.4f%%)\n",
            mean(corr_profile$corr_delta), mean(corr_profile$corr_delta) * 100))
cat(sprintf("  Uncorrected equal-weight delta:     %.6f (%.4f%%)\n",
            mean(ew_profile$ew_delta), mean(ew_profile$ew_delta) * 100))

corr_tsv <- file.path(out_dir, "metagene_lambda_corrected.tsv")
full_corr <- merge(corr_profile, ew_profile, by = "bin")
write.table(full_corr, corr_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nPlotting lambda-corrected metagene...\n")

plot_corr <- rbind(
  data.table(bin = corr_profile$bin, rate = corr_profile$corr_ctrl * 100,
             se = corr_profile$corr_ctrl_se * 100,
             genotype = "Control", method = "Lambda-corrected"),
  data.table(bin = corr_profile$bin, rate = corr_profile$corr_mut * 100,
             se = corr_profile$corr_mut_se * 100,
             genotype = "Mutant", method = "Lambda-corrected"),
  data.table(bin = ew_profile$bin, rate = ew_profile$ew_ctrl * 100,
             se = ew_profile$ew_ctrl_se * 100,
             genotype = "Control", method = "Uncorrected"),
  data.table(bin = ew_profile$bin, rate = ew_profile$ew_mut * 100,
             se = ew_profile$ew_mut_se * 100,
             genotype = "Mutant", method = "Uncorrected")
)

p_corr <- ggplot(plot_corr, aes(x = bin, y = rate, color = genotype, linetype = method)) +
  geom_ribbon(aes(ymin = rate - se, ymax = rate + se, fill = genotype),
              alpha = 0.15, linetype = 0) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("Control" = "#3A7CA5", "Mutant" = "#D4713B")) +
  scale_fill_manual(values = c("Control" = "#3A7CA5", "Mutant" = "#D4713B")) +
  scale_linetype_manual(values = c("Lambda-corrected" = "solid", "Uncorrected" = "dashed")) +
  labs(x = "Gene body position (5' to 3')",
       y = "mCH rate (%)",
       title = "Metagene mCH: lambda-corrected vs uncorrected",
       subtitle = "Equal-weight per sample. Solid = lambda-subtracted, Dashed = raw",
       color = "Genotype", fill = "Genotype", linetype = "Method") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

corr_plot_path <- file.path(out_dir, "metagene_lambda_corrected.png")
ggsave(corr_plot_path, p_corr, width = 8, height = 5, dpi = 150)
cat(sprintf("Saved to %s\n", corr_plot_path))

p_corr_delta <- ggplot() +
  geom_ribbon(data = data.table(bin = corr_profile$bin,
                                delta = corr_profile$corr_delta * 100,
                                se = sqrt(corr_profile$corr_ctrl_se^2 +
                                          corr_profile$corr_mut_se^2) * 100),
              aes(x = bin, ymin = delta - se, ymax = delta + se),
              fill = "#2E6B8A", alpha = 0.2) +
  geom_line(data = data.table(bin = corr_profile$bin,
                              delta = corr_profile$corr_delta * 100,
                              method = "Lambda-corrected"),
            aes(x = bin, y = delta, linetype = method), color = "#2E6B8A", linewidth = 0.8) +
  geom_line(data = data.table(bin = ew_profile$bin,
                              delta = ew_profile$ew_delta * 100,
                              method = "Uncorrected"),
            aes(x = bin, y = delta, linetype = method), color = "#D4713B", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  scale_linetype_manual(values = c("Lambda-corrected" = "solid", "Uncorrected" = "dashed")) +
  labs(x = "Gene body position (5' to 3')",
       y = "Delta mCH (mut - ctrl, percentage points)",
       title = "Metagene delta: lambda-corrected vs uncorrected",
       subtitle = "Equal-weight per sample. SE ribbon on corrected.",
       linetype = "Method") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

corr_delta_path <- file.path(out_dir, "metagene_delta_lambda_corrected.png")
ggsave(corr_delta_path, p_corr_delta, width = 8, height = 4, dpi = 150)
cat(sprintf("Saved to %s\n", corr_delta_path))

cat("\nDone.\n")
