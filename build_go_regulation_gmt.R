#!/usr/bin/env Rscript

wrapper_path <- sys.frame(1)$ofile
if (is.null(wrapper_path)) {
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
  wrapper_path <- if (length(file_arg)) file_arg[[1]] else "build_go_regulation_gmt.R"
}
project_root <- dirname(normalizePath(wrapper_path, winslash = "/", mustWork = TRUE))

source(file.path(project_root, "R", "load_project_code.R"), local = FALSE)

if (sys.nframe() == 0) {
  library(org.Hs.eg.db)

  result_bp <- build_go_regulation_library(
    org_db = org.Hs.eg.db,
    include_IEA = TRUE,
    pos_keywords = DEFAULT_POS_KEYWORDS,
    neg_keywords = DEFAULT_NEG_KEYWORDS,
    output_dir = "go_dedup_gmt"
  )

  result_cc <- build_go_dedup_library(
    org_db = org.Hs.eg.db,
    ontology = "CC",
    include_IEA = TRUE,
    output_dir = "go_dedup_gmt"
  )

  result_mf <- build_go_dedup_library(
    org_db = org.Hs.eg.db,
    ontology = "MF",
    include_IEA = TRUE,
    output_dir = "go_dedup_gmt"
  )
}
