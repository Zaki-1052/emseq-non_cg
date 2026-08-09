# Load required libraries
library(GenomicRanges)
library(ggplot2)
library(dplyr)

# ---- 1. Read in your two DiffBind results ----
mecp2 <- read.table("/Users/janviramchandra/Desktop/mecp2.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
h2ak119ub <- read.table("/Users/janviramchandra/Desktop/h2aub.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ---- 2. Convert to GRanges ----
gr_mecp2 <- GRanges(seqnames = mecp2$seqnames,
                    ranges = IRanges(start = mecp2$start, end = mecp2$end))
mcols(gr_mecp2) <- mecp2[, c("Fold", "p.value", "FDR")]

gr_h2a <- GRanges(seqnames = h2ak119ub$seqnames,
                  ranges = IRanges(start = h2ak119ub$start, end = h2ak119ub$end))
mcols(gr_h2a) <- h2ak119ub[, c("Fold", "p.value", "FDR")]

# ---- 3. Find overlaps between the two peaksets ----
hits <- findOverlaps(gr_mecp2, gr_h2a)

n_common_peaks <- length(hits)
cat("Number of overlapping loci between MeCP2 and H2AK119Ub peaksets:", n_common_peaks, "\n")

# Build merged data frame from overlapping pairs
merged <- data.frame(
  mecp2_Fold = mcols(gr_mecp2)$Fold[queryHits(hits)],
  mecp2_FDR  = mcols(gr_mecp2)$FDR[queryHits(hits)],
  h2a_Fold   = mcols(gr_h2a)$Fold[subjectHits(hits)],
  h2a_FDR    = mcols(gr_h2a)$FDR[subjectHits(hits)]
)

# ---- 4. Compute overlap fraction of each MeCP2 peak covered by H2AUb ----
overlap_widths <- width(pintersect(gr_mecp2[queryHits(hits)], gr_h2a[subjectHits(hits)]))
mecp2_widths <- width(gr_mecp2[queryHits(hits)])
h2a_widths <- width(gr_h2a[queryHits(hits)])
merged$overlap_fraction_of_mecp2 <- overlap_widths / mecp2_widths
merged$overlap_fraction_of_h2a <- overlap_widths / h2a_widths

# Filter: require at least 80% of the MeCP2 peak to be covered by H2AUb
merged <- merged %>% filter(merged$overlap_fraction_of_mecp2 >= 0.5 & merged$overlap_fraction_of_h2a >= 0.5)

# Update n_common_peaks to reflect the filtered count
n_common_peaks <- nrow(merged)
cat("Number of overlapping loci after >=50% overlap filter:", n_common_peaks, "\n")

# ---- 5. Define significance ----
fdr_cutoff <- 0.05

# Option A: significant in EITHER dataset
# merged$sig <- ifelse(merged$mecp2_FDR < fdr_cutoff | merged$h2a_FDR < fdr_cutoff,
                 #    "Significant", "Not significant")

# Option B: significant in BOTH datasets (uncomment to use instead)
merged$sig <- ifelse(merged$mecp2_FDR < fdr_cutoff & merged$h2a_FDR < fdr_cutoff,
                      "Significant", "Not significant")

merged$sig <- factor(merged$sig, levels = c("Not significant", "Significant"))

# ---- 6. Assign quadrants ----
# Q1: MeCP2 up,   H2A up     (top right)
# Q2: MeCP2 down, H2A up     (top left)
# Q3: MeCP2 down, H2A down   (bottom left)
# Q4: MeCP2 up,   H2A down   (bottom right)
merged$quadrant <- case_when(
  merged$mecp2_Fold >= 0 & merged$h2a_Fold >= 0 ~ "Q1 (MeCP2 up / H2A up)",
  merged$mecp2_Fold <  0 & merged$h2a_Fold >= 0 ~ "Q2 (MeCP2 down / H2A up)",
  merged$mecp2_Fold <  0 & merged$h2a_Fold <  0 ~ "Q3 (MeCP2 down / H2A down)",
  merged$mecp2_Fold >= 0 & merged$h2a_Fold <  0 ~ "Q4 (MeCP2 up / H2A down)"
)

# ---- 7. Quadrant counts: ALL merged loci ----
quadrant_counts_all <- merged %>%
  count(quadrant) %>%
  mutate(percent = round(100 * n /n_common_peaks, 2))

cat("\nQuadrant breakdown — ALL overlapping loci:\n")
print(quadrant_counts_all)

# ---- 8. Quadrant counts: SIGNIFICANT loci only ----
library(tidyr)

quadrant_counts_sig <- merged %>%
  filter(sig == "Significant") %>%
  count(quadrant) %>%
  complete(
    quadrant = c(
      "Q1 (MeCP2 up / H2A up)",
      "Q2 (MeCP2 down / H2A up)",
      "Q3 (MeCP2 down / H2A down)",
      "Q4 (MeCP2 up / H2A down)"
    ),
    fill = list(n = 0)
  ) %>%
  mutate(
    percent = round(100 * n /n_common_peaks, 2),
    label = paste0("n = ", n, "\n", percent, "%")
  )

cat("\nQuadrant breakdown — SIGNIFICANT loci only:\n")
print(quadrant_counts_sig)

# ---- 9. Plot ----

merged <- merged %>% arrange(sig)

# Determine plot limits for placing labels
xmax <- max(abs(merged$mecp2_Fold), na.rm = TRUE)
ymax <- max(abs(merged$h2a_Fold), na.rm = TRUE)

quad_labels <- quadrant_counts_sig %>%
  mutate(
    x = c(
      0.95 * xmax,   # Q1
      -0.95 * xmax,   # Q2
      -0.92 * xmax,   # Q3
      0.95 * xmax    # Q4
    ),
    y = c(
      0.97 * ymax,   # Q1
      0.95 * ymax,   # Q2
      -0.95 * ymax,   # Q3
      -0.95 * ymax    # Q4
    )
  )

p <- ggplot(merged, aes(x = mecp2_Fold, y = h2a_Fold, color = sig)) +
  geom_point(size = 1.5, alpha = 0.8) +
  
  scale_color_manual(values = c(
    "Not significant" = "grey85",
    "Significant" = "deeppink"
  )) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  
  geom_text(
    data = quad_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  
  labs(
    x = "MeCP2 log2(Fold Change, mutant/control)",
    y = "H2AK119Ub log2(Fold Change, mutant/control)",
    color = "",
    title = paste0(
      "n = ",
      n_common_peaks,
      " overlapping loci (50% of MeCP2 peak and 50% of H2AUb peak must overlap)"
    )
  ) +
  
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top"
  )

print(p)

# ---- 10. Save ----

# ggsave("MeCP2_vs_H2AK119Ub_quadrant_plot.png",
#        p, width = 7, height = 7, dpi = 300)

# ggsave("MeCP2_vs_H2AK119Ub_quadrant_plot.svg",
#        p, width = 7, height = 7)

# ---- 11. Save quadrant summary ----

# write.csv(quadrant_counts_sig,
#           "quadrant_counts_significant_loci.csv",
#           row.names = FALSE)

# ---- 10. Save ----
#ggsave("MeCP2_vs_H2AK119Ub_quadrant_plot.png", p, width = 7, height = 7, dpi = 300)
#ggsave("MeCP2_vs_H2AK119Ub_quadrant_plot.svg", p, width = 7, height = 7)

# ---- 11. Save quadrant summary tables to file ----
#write.csv(quadrant_counts_all, "quadrant_counts_all_loci.csv", row.names = FALSE)
#write.csv(quadrant_counts_sig, "quadrant_counts_significant_loci.csv", row.names = FALSE)
