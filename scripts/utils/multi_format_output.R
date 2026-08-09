# scripts/utils/multi_format_output.R
# Utility functions for multi-format plot output (PDF + SVG + PNG + JPEG)
#
# Purpose:
#   Provides helper functions to output plots in four formats simultaneously:
#   - PDF (publication standard, vector)
#   - SVG (Illustrator-friendly, editable vector)
#   - PNG (raster, general use)
#   - JPEG (Google Slides, presentations)
#
# Usage:
#   source("scripts/utils/multi_format_output.R")
#
#   # For ggplot2 objects:
#   save_multiformat_ggplot(my_plot, "path/to/output/filename", width = 10, height = 8)
#
#   # For base R graphics:
#   save_multiformat_base(quote({ plot(x, y); lines(x, y) }), "path/to/output/filename", width = 10, height = 8)
#
# Note: Requires svglite package for high-quality SVG output
#   install.packages("svglite")


library(svglite)

save_multiformat_ggplot <- function(plot, base_path, width = 10, height = 8, dpi = 300, verbose = TRUE, use_subfolders = TRUE) {
  figure_name <- basename(base_path)
  parent_dir <- dirname(base_path)

  if (use_subfolders) {
    output_dir <- file.path(parent_dir, figure_name)
    file_prefix <- file.path(output_dir, figure_name)
  } else {
    output_dir <- parent_dir
    file_prefix <- base_path
  }

  if (!dir.exists(output_dir) && output_dir != ".") {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  pdf_path <- paste0(file_prefix, ".pdf")
  ggplot2::ggsave(pdf_path, plot, width = width, height = height)

  svg_path <- paste0(file_prefix, ".svg")
  ggplot2::ggsave(svg_path, plot, width = width, height = height, device = svglite::svglite)

  png_path <- paste0(file_prefix, ".png")
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = dpi, device = "png")

  jpg_path <- paste0(file_prefix, ".jpg")
  ggplot2::ggsave(jpg_path, plot, width = width, height = height, dpi = dpi, device = "jpeg")

  if (verbose) {
    if (use_subfolders) {
      cat(sprintf("  Saved: %s/{pdf,svg,png,jpg}\n", figure_name))
    } else {
      cat(sprintf("  Saved: %s.{pdf,svg,png,jpg}\n", figure_name))
    }
  }

  invisible(plot)
}


save_multiformat_base <- function(plot_expr, base_path, width = 10, height = 8, dpi = 300, verbose = TRUE, use_subfolders = TRUE) {
  figure_name <- basename(base_path)
  parent_dir <- dirname(base_path)

  if (use_subfolders) {
    output_dir <- file.path(parent_dir, figure_name)
    file_prefix <- file.path(output_dir, figure_name)
  } else {
    output_dir <- parent_dir
    file_prefix <- base_path
  }

  if (!dir.exists(output_dir) && output_dir != ".") {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  pdf_path <- paste0(file_prefix, ".pdf")
  pdf(pdf_path, width = width, height = height)
  tryCatch({
    eval(plot_expr)
  }, finally = {
    dev.off()
  })

  svg_path <- paste0(file_prefix, ".svg")
  svglite::svglite(svg_path, width = width, height = height)
  tryCatch({
    eval(plot_expr)
  }, finally = {
    dev.off()
  })

  png_path <- paste0(file_prefix, ".png")
  png(png_path, width = width * dpi, height = height * dpi, res = dpi)
  tryCatch({
    eval(plot_expr)
  }, finally = {
    dev.off()
  })

  jpg_path <- paste0(file_prefix, ".jpg")
  jpeg(jpg_path, width = width * dpi, height = height * dpi, res = dpi, quality = 95)
  tryCatch({
    eval(plot_expr)
  }, finally = {
    dev.off()
  })

  if (verbose) {
    if (use_subfolders) {
      cat(sprintf("  Saved: %s/{pdf,svg,png,jpg}\n", figure_name))
    } else {
      cat(sprintf("  Saved: %s.{pdf,svg,png,jpg}\n", figure_name))
    }
  }

  invisible(NULL)
}
