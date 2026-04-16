#!/usr/bin/env Rscript

wrapper_path <- sys.frame(1)$ofile
if (is.null(wrapper_path)) {
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
  wrapper_path <- if (length(file_arg)) file_arg[[1]] else "build_pathway_library.R"
}
project_root <- dirname(normalizePath(wrapper_path, winslash = "/", mustWork = TRUE))

source(file.path(project_root, "R", "load_project_code.R"), local = FALSE)

if (sys.nframe() == 0) {
  library(org.Hs.eg.db)

  result <- build_pathway_library(
    org_db = org.Hs.eg.db,
    species = "Homo sapiens",
    output_dir = "pathway_library"
  )
}
