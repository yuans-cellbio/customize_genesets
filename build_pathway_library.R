#!/usr/bin/env Rscript

wrapper_path <- sys.frame(1)$ofile
if (is.null(wrapper_path)) {
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
  wrapper_path <- if (length(file_arg)) file_arg[[1]] else "build_pathway_library.R"
}
project_root <- dirname(normalizePath(wrapper_path, winslash = "/", mustWork = TRUE))

source(file.path(project_root, "R", "load_project_code.R"), local = FALSE)
load_project_code(project_root)

if (sys.nframe() == 0) {
  result <- build_pathway_library(
    species = "human",
    output_dir = "pathway_library"
  )
}
