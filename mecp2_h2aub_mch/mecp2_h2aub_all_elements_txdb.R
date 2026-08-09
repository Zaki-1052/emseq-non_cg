############################################################
# MeCP2 vs H2AK119Ub gene-level analysis
#
# Flow:
# MeCP2 and H2AUb peaks
#      |
#      v
# mm10 gene annotation
#      |
#      v
# Multiple peaks per gene
#      |
#      v
# Peak weight = -log10(FDR)
#      |
#      v
# Weighted gene-level scores (weighted fold changes)

#NOTE: you can comment out the gene parts that you dont want to study.

############################################################

#libraries
library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(DOSE)

#upload diffbind .txt peak curated files (Conc = 4 or more)
mecp2 <- read.table(
  "/Users/janviramchandra/Desktop/mecp2.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

h2ak119ub <- read.table(
  "/Users/janviramchandra/Desktop/h2aub.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)


#convert peaks to GRanges
gr_mecp2 <- GRanges(
  seqnames = mecp2$seqnames,
  ranges = IRanges(
    start = mecp2$start,
    end = mecp2$end
  ),
  peak_id = paste0("mecp2_peak_", seq_len(nrow(mecp2))),
  Fold = mecp2$Fold,
  FDR = mecp2$FDR
)

gr_h2a <- GRanges(
  seqnames = h2ak119ub$seqnames,
  ranges = IRanges(
    start = h2ak119ub$start,
    end = h2ak119ub$end
  ),
  peak_id = paste0("h2a_peak_", seq_len(nrow(h2ak119ub))),
  Fold = h2ak119ub$Fold,
  FDR = h2ak119ub$FDR
)

#annotate MeCP2 and H2AUb peaks to genes (mm10)
peakAnno_mecp2 <- annotatePeak(
  gr_mecp2,
  TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
  annoDb = "org.Mm.eg.db",
  tssRegion = c(-2000, 2000)
)

anno_mecp2 <- as.data.frame(peakAnno_mecp2)

anno_mecp2$mecp2_Fold <- mcols(gr_mecp2)$Fold[
  match(
    anno_mecp2$peak_id,
    mcols(gr_mecp2)$peak_id
  )
]

anno_mecp2$mecp2_FDR <- mcols(gr_mecp2)$FDR[
  match(
    anno_mecp2$peak_id,
    mcols(gr_mecp2)$peak_id
  )
]

peakAnno_h2a <- annotatePeak(
  gr_h2a,
  TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
  annoDb = "org.Mm.eg.db",
  tssRegion = c(-2000, 2000),
)

anno_h2a <- as.data.frame(peakAnno_h2a)

anno_h2a$h2a_Fold <- mcols(gr_h2a)$Fold[
  match(
    anno_h2a$peak_id,
    mcols(gr_h2a)$peak_id
  )
]

anno_h2a$h2a_FDR <- mcols(gr_h2a)$FDR[
  match(
    anno_h2a$peak_id,
    mcols(gr_h2a)$peak_id
  )
]

#building annotation table
mecp2_df <- data.frame(
  
  SYMBOL =
    anno_mecp2$SYMBOL,
  
  geneID =
    anno_mecp2$geneId,
  
  annotation =
    anno_mecp2$annotation,
  
  distanceToTSS =
    anno_mecp2$distanceToTSS,
  
  mecp2_Fold =
    anno_mecp2$mecp2_Fold,
  
  mecp2_FDR =
    anno_mecp2$mecp2_FDR
  
)

h2a_df <- data.frame(
  
  SYMBOL =
    anno_h2a$SYMBOL,
  
  geneID =
    anno_h2a$geneId,
  
  annotation =
    anno_h2a$annotation,
  
  distanceToTSS =
    anno_h2a$distanceToTSS,
  
  h2a_Fold =
    anno_h2a$h2a_Fold,
  
  h2a_FDR =
    anno_h2a$h2a_FDR
  
)
table(anno_mecp2$annotation)

#filter for what you want
mecp2_df <- mecp2_df %>%
  mutate(
    feature = case_when(
 #     grepl("^Downstream", annotation) ~ "Downstream",
#      grepl("^Distal Intergenic", annotation) ~ "Distal Intergenic",
#      grepl("^Promoter", annotation) ~ "Promoter",
      grepl("^Exon", annotation) ~ "Exon",
      grepl("^Intron", annotation) ~ "Intron",
      grepl("^5' UTR", annotation) ~ "5' UTR",
      grepl("^3' UTR", annotation) ~ "3' UTR",
      TRUE ~ annotation
    )
  )

mecp2_df <- mecp2_df %>%
  filter(
    feature %in% c(
#      "Downstream",
#      "Distal Intergenic",
 #     "Promoter",
      "Exon",
      "Intron",
      "5' UTR",
      "3' UTR"
    ),
    !is.na(SYMBOL)
  )

h2a_df <- h2a_df %>%
  mutate(
    feature = case_when(
#      grepl("^Downstream", annotation) ~ "Downstream",
#      grepl("^Distal Intergenic", annotation) ~ "Distal Intergenic",
 #     grepl("^Promoter", annotation) ~ "Promoter",
      grepl("^Exon", annotation) ~ "Exon",
      grepl("^Intron", annotation) ~ "Intron",
      grepl("^5' UTR", annotation) ~ "5' UTR",
      grepl("^3' UTR", annotation) ~ "3' UTR",
      TRUE ~ annotation
    )
  )

h2a_df <- h2a_df %>%
  filter(
    feature %in% c(
 #     "Downstream",
#    "Distal Intergenic",
  #    "Promoter",
      "Exon",
      "Intron",
      "5' UTR",
      "3' UTR"
    ),
    !is.na(SYMBOL)
  )

#annotation category counts
# MeCP2 peak annotation counts
mecp2_annotation_counts <- anno_mecp2 %>%
  count(annotation) %>%
  arrange(desc(n))

print(mecp2_annotation_counts)

# H2AK119Ub peak annotation counts
h2a_annotation_counts <- anno_h2a %>%
  count(annotation) %>%
  arrange(desc(n))

print(h2a_annotation_counts)

# Extract annotation category only
annotation_summary_mecp2 <- anno_mecp2 %>%
  mutate(
    annotation_category = sub("\\s*\\(.*", "", annotation)
  ) %>%
  count(annotation_category, name = "peak_count") %>%
  arrange(desc(peak_count))

annotation_summary_h2a <- anno_h2a %>%
  mutate(
    annotation_category = sub("\\s*\\(.*", "", annotation)
  ) %>%
  count(annotation_category, name = "peak_count") %>%
  arrange(desc(peak_count))

# Print all annotation categories
cat("\n===== Annotation Categories =====\n")
print(annotation_summary_mecp2)
print(annotation_summary_h2a)

############################################################
#Calculate peak weights and collapse to gene level

# Peak Weight = -log(FDR)
# Gene FC = Sum(FC*peak weight of all peaks)/Sum(peak weight of all peaks)
############################################################

mecp2_df <- mecp2_df %>%
  mutate(
    mecp2_weight = -log10(pmax(mecp2_FDR,1e-300))
  )

h2a_df <- h2a_df %>%
  mutate(
    h2a_weight = -log10(pmax(h2a_FDR,1e-300))
  )

#collapse MeCP2 and H2AUb peaks to gene level
mecp2_gene <- mecp2_df %>%
  group_by(SYMBOL) %>%
  summarise(
    
    mecp2_Fold =
      sum(mecp2_Fold * mecp2_weight, na.rm=TRUE) /
      sum(mecp2_weight, na.rm=TRUE),
    
    mecp2_minFDR =
      min(mecp2_FDR, na.rm=TRUE),
    
    mecp2_sig_peaks =
      sum(mecp2_FDR < 0.05),
    
    n_mecp2_peaks=n(),
    
    .groups="drop"
    
  ) %>%
  mutate(
    mecp2_sig = mecp2_sig_peaks > 0
  )

h2a_gene <- h2a_df %>%
  group_by(SYMBOL) %>%
  summarise(
    
    h2a_Fold =
      sum(h2a_Fold * h2a_weight, na.rm=TRUE) /
      sum(h2a_weight, na.rm=TRUE),
    
    h2a_minFDR =
      min(h2a_FDR, na.rm=TRUE),
    
    h2a_sig_peaks =
      sum(h2a_FDR < 0.05),
    
    n_h2a_peaks =
      n(),
    
    .groups="drop"
    
  ) %>%
  mutate(
    h2a_sig = h2a_sig_peaks > 0
  )

############################################################
# Merge MeCP2 and H2AK119Ub gene-level tables
# Only genes containing BOTH datasets are retained
############################################################

gene_df <- inner_join(
  mecp2_gene,
  h2a_gene,
  by = "SYMBOL"
)

cat(
  "Number of genes with both MeCP2 and H2AK119Ub peaks:",
  nrow(gene_df),
  "\n"
)

############################################################
#Significance: if gene has at least one sig peak in both MeCP2 AND H2AUb
############################################################

fdr_cutoff <- 0.05

gene_df$sig <- ifelse(
  gene_df$mecp2_sig &
    gene_df$h2a_sig,
  "Significant",
  "Not significant"
)

gene_df$sig <- factor(
  gene_df$sig,
  levels = c(
    "Not significant",
    "Significant"
  )
  
)

############################################################
#Assign quadrants
############################################################

gene_df$quadrant <- case_when(
  
  gene_df$mecp2_Fold >= 0 &
    gene_df$h2a_Fold >= 0 ~
    "Q1 (MeCP2 up / H2A up)",
  
  gene_df$mecp2_Fold < 0 &
    gene_df$h2a_Fold >= 0 ~
    "Q2 (MeCP2 down / H2A up)",
  
  gene_df$mecp2_Fold < 0 &
    gene_df$h2a_Fold < 0 ~
    "Q3 (MeCP2 down / H2A down)",
  
  gene_df$mecp2_Fold >= 0 &
    gene_df$h2a_Fold < 0 ~
    "Q4 (MeCP2 up / H2A down)"
  
)

############################################################
# Spearman correlation (just needs fold changes)
############################################################

cor_test <- cor.test(
  gene_df$mecp2_Fold,
  gene_df$h2a_Fold,
  method = "spearman"
)

rho_label <- paste0(
  "Spearman \u03C1 = ",
  round(
    cor_test$estimate,
    2
  ),
  " | p = ",
  format.pval(
    cor_test$p.value,
    digits = 2,
    eps = 1e-10
  )
)

print(cor_test)

############################################################
#Quadrant stats
############################################################

total_genes <- nrow(gene_df)

total_sig_genes <- sum(
  gene_df$sig == "Significant"
)

#all genes per quadrant
quadrant_all <- gene_df %>%
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
  rename(
    n_total_genes = n
  ) %>%
  mutate(
    percent_total_genes =
      round(
        100 * n_total_genes / total_genes,
        1
      )
  )

#significant genes per quadrant
quadrant_sig <- gene_df %>%
  filter(
    sig == "Significant"
  ) %>%
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
  rename(
    n_sig_genes = n
  ) %>%
  mutate(
    percent_sig_genes =
      round(
        100 * n_sig_genes / total_sig_genes,
        1
      )
  )

#combine results
quadrant_counts <- left_join(
  quadrant_all,
  quadrant_sig,
  by = "quadrant"
)

#labels to appear in each quadrant with Ns and %s
quadrant_counts <- quadrant_counts %>%
  mutate(
    label =
      paste0(
        "Total: n=",
        n_total_genes,
        " (",
        percent_total_genes,
        "%)\n",
        "Sig: n=",
        n_sig_genes,
        " (",
        percent_sig_genes,
        "%)"
      )
  )

print(quadrant_counts)

#quadrant label positions
xmax <- max(
  abs(gene_df$mecp2_Fold),
  na.rm = TRUE
)

ymax <- max(
  abs(gene_df$h2a_Fold),
  na.rm = TRUE
)


quadrant_counts$x <- c(
  0.90*xmax,
  -0.30*xmax,
  -0.30*xmax,
  0.90*xmax
)

quadrant_counts$y <- c(
  0.10*ymax,
  0.95*ymax,
  -0.95*ymax,
  -0.95*ymax
)

############################################################
#label top 10 genes overall
#Label top 5 genes in quadrant 3
############################################################

#Top 10 strongest significant genes overall
#they all just happen to be in quadrant 1
#MeCP2 up / H2AK119Ub up
top_genes_all <- gene_df %>%
  filter(
    sig == "Significant"
  ) %>%
  
  mutate(
    label_strength =
      abs(mecp2_Fold) +
      abs(h2a_Fold)
  ) %>%
  
  arrange(
    desc(label_strength)
  ) %>%
  
  slice_head(
    n = 10
  )

#Top 5 Q3 genes specifically
#MeCP2 down / H2AK119Ub down
top_genes_Q3 <- gene_df %>%
  filter(
    sig == "Significant",
    quadrant ==
      "Q3 (MeCP2 down / H2A down)"
  ) %>%
  
  mutate(
    label_strength =
      abs(mecp2_Fold) +
      abs(h2a_Fold)
  ) %>%
  
  arrange(
    desc(label_strength)
  ) %>%
  
  slice_head(
    n = 5
  )

#don't duplicate labels if there are any overlaps
top_genes <- bind_rows(
  top_genes_all,
  top_genes_Q3
) %>%
  
  distinct(
    SYMBOL,
    .keep_all = TRUE
  )

############################################################
#generate scatter log2log2 plot
############################################################

p <- ggplot(
  gene_df,
  aes(
    x = mecp2_Fold,
    y = h2a_Fold,
    color = sig
  )
  
) +
  
  geom_point(
    size = 1.8,
    alpha = 0.8
    
  ) +
  
  scale_color_manual(
    values = c(
      "Not significant" = "grey85",
      "Significant" = "deeppink"
    )
    
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black"
    
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black"
    
  ) +
  
  #quadrant labels
  geom_text(
    data = quadrant_counts,
    aes(
      x = x,
      y = y,
      label = label
    ),
    
    inherit.aes = FALSE,
    size = 3,
    fontface = "bold",
    color = "black"
    
  ) +
  
  #gene labels
  geom_text_repel(
    data = top_genes,
    aes(
      x = mecp2_Fold,
      y = h2a_Fold,
      label = SYMBOL
    ),
    
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    fontface = "bold",
    color = "black"
    
  ) +
  
  #label axes
  labs(
    x = "MeCP2 log2(Fold Change, mutant/control)",
    y = "H2AK119Ub log2(Fold Change, mutant/control)",
    color = "",
    title = "Gene-level MeCP2 vs H2AK119Ub",
    subtitle =
      paste0(
        "n (total) = ",
        total_genes,
        " genes | ",
        "n (sig) = ",
        total_sig_genes,
        " genes ",
        "\n",
        rho_label
      )
    
  ) +
  
  theme_bw(
    base_size = 14
    
  ) +
  
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    plot.subtitle = element_text(
      size = 12,
      face = "bold"
    )
    
  )

print(p)

############################################################
#GO analyses
#ONLY significant genes in quadrant 1
# MeCP2 up / H2AK119Ub up
############################################################
universe_genes <- gene_df$SYMBOL

universe_entrez <- bitr(
  universe_genes,
  fromType="SYMBOL",
  toType="ENTREZID",
  OrgDb=org.Mm.eg.db
)

Q1_sig <- gene_df %>%
  filter(
    sig == "Significant",
    quadrant == "Q1 (MeCP2 up / H2A up)"
  )

Q1_symbols <- Q1_sig$SYMBOL

cat(
  "Number of significant Q1 genes:",
  length(Q1_symbols),
  "\n"
)

#convert gene symbols to Entrez IDs
Q1_entrez <- bitr(
  Q1_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

cat(
  "Genes converted to Entrez:",
  nrow(Q1_entrez),
  "\n"
)

#GO Biological Process
GO_BP <- enrichGO(
  gene = Q1_entrez$ENTREZID,
  universe = universe_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont="BP",
  pAdjustMethod="BH",
  readable=TRUE
)

head(GO_BP)

plot_GO_BP <- dotplot(
  GO_BP,
  showCategory = 20,
  title = "GO enrichment: Q1 MeCP2↑ / H2AK119Ub↑ genes"
)

#GO Molecular Function
GO_MF <- enrichGO(
  gene = Q1_entrez$ENTREZID,
  universe = universe_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont="MF",
  pAdjustMethod="BH",
  readable=TRUE
)

plot_GO_MF <-dotplot(
  GO_MF,
  showCategory = 20,
  title = "GO Molecular Function: Q1 MeCP2↑ / H2AK119Ub↑ genes"
)

#GO Cellular Component
GO_CC <- enrichGO(
  gene = Q1_entrez$ENTREZID,
  universe = universe_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont="CC",
  pAdjustMethod="BH",
  readable=TRUE
)

plot_GO_CC <-dotplot(
  GO_CC,
  showCategory = 20,
  title = "GO Cellular Component: Q1 MeCP2↑ / H2AK119Ub↑ genes"
)

############################################################
# 20. Save all plots
############################################################
setwd("/Users/janviramchandra/Desktop/adult_mecp2_vs_h2aub/just_genebody")

#save gene-level results as csv
write.csv(
  gene_df,
  "MeCP2_H2AK119Ub_gene_level_results.csv",
  row.names = FALSE
)

write.csv(
  quadrant_counts,
  "MeCP2_H2AK119Ub_gene_level_quadrants.csv",
  row.names = FALSE
)

write.csv(
  top_genes,
  "MeCP2_H2AK119Ub_top_labeled_genes.csv",
  row.names = FALSE
)
write.csv(
  anno_mecp2,
  "MeCP2_peak_annotation_table.csv",
  row.names = FALSE
)

write.csv(
  anno_h2a,
  "H2AK119Ub_peak_annotation_table.csv",
  row.names = FALSE
)

#save log2log2 quadrant plot
ggsave(
  "MeCP2_vs_H2AK119Ub_gene_level_plot.png",
  p,
  width = 10,
  height = 10,
  dpi = 300
)
ggsave(
  "MeCP2_vs_H2AK119Ub_gene_level_plot.svg",
  p,
  width = 10,
  height = 10
)

#save all GO analyses data as csv
write.csv(
  as.data.frame(GO_BP),
  "Q1_MeCP2_up_H2AK119Ub_up_GO_results.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(GO_MF),
  "Q1_MeCP2_up_H2AK119Ub_up_GO_MF_results.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(GO_CC),
  "Q1_MeCP2_up_H2AK119Ub_up_GO_CC_results.csv",
  row.names = FALSE
)

#save GO analyses plots
ggsave(
  "Q1_GO_BP_dotplot.png",
  plot_GO_BP,
  width = 8,
  height = 10,
  dpi = 300
)

ggsave(
  "Q1_GO_MF_dotplot.png",
  plot_GO_MF,
  width = 8,
  height = 10,
  dpi = 300
)

ggsave(
  "Q1_GO_CC_dotplot.png",
  plot_GO_CC,
  width = 8,
  height = 10,
  dpi = 300
)
