# scripts/install_r_packages.R
#
# Install R packages required for the mCH pipeline.
# Run once before submitting the pipeline if any packages are missing.
#
# Usage: Rscript install_r_packages.R

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c(
  # --- Core pipeline (steps 01-08) ---
  "edgeR", "limma", "GenomicRanges", "IRanges", "rtracklayer",
  "bsseq", "dmrseq", "BiocParallel",
  "clusterProfiler", "enrichplot", "DOSE", "fgsea", "goseq",
  "ChIPseeker", "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "org.Mm.eg.db", "GO.db", "AnnotationDbi",
  "data.table", "optparse",
  "here", "readr", "dplyr", "tidyr",
  "ggplot2", "patchwork", "ggrepel", "svglite",
  "clinfun",

  # --- Section pipeline: shared config ---
  "yaml",             # paths.yaml parsing
  "purrr",            # loaded by _shared_config.R
  "tibble",           # loaded by _shared_config.R
  "stringr",          # loaded by _shared_config.R
  "RColorBrewer",     # colour palettes
  "scales",           # axis formatting
  "pheatmap",         # correlation heatmaps (section 20_02)

  # --- Section pipeline: sub-gene features ---
  "GenomicFeatures",  # GTF parsing in generate_feature_beds.R
  "dunn.test",        # post-hoc pairwise tests (section 50_01)

  # --- Section pipeline: permutation tests (sections 40_01-40_04) ---
  "regioneR",
  "regioneReloaded",
  "BSgenome.Mmusculus.UCSC.mm10",

  # --- Section pipeline: models and set overlap ---
  "pROC",             # AUC for logistic models (section 20_02)
  "ggVennDiagram"     # gene set overlap (sections 60_01, 70_03)
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing:", pkg, "\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat("OK:", pkg, "(", as.character(packageVersion(pkg)), ")\n")
  }
}

cat("\nAll packages verified.\n")
