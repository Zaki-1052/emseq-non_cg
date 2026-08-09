# scripts/install_r_packages.R
#
# Install R packages required for the mCH pipeline.
# Run once before submitting the pipeline if any packages are missing.
#
# Usage: Rscript install_r_packages.R

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c(
  "edgeR", "limma", "GenomicRanges", "IRanges", "rtracklayer",
  "bsseq", "dmrseq", "BiocParallel",
  "clusterProfiler", "enrichplot", "DOSE", "fgsea", "goseq",
  "ChIPseeker", "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "org.Mm.eg.db", "GO.db", "AnnotationDbi",
  "data.table", "optparse",
  "here", "readr", "dplyr", "tidyr",
  "ggplot2", "patchwork", "ggrepel", "svglite",
  "clinfun"
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
