wrapper_path <- sys.frame(1)$ofile
if (is.null(wrapper_path)) {
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
  wrapper_path <- if (length(file_arg)) file_arg[[1]] else "build_human_pathway_library.R"
}
project_root <- dirname(normalizePath(wrapper_path, winslash = "/", mustWork = TRUE))

source(file.path(project_root, "R", "load_project_code.R"), local = FALSE)

library(org.Hs.eg.db)

pathways <- build_pathway_library(
  org_db = org.Hs.eg.db,
  species = "Homo sapiens",
  kegg_organism = "hsa",
  progeny_organism = "human",
  include_IEA = TRUE,
  pos_keywords = NULL,
  neg_keywords = NULL,
  min_size = 15,
  max_size = 500,
  dag_overlap_threshold = 0.8,
  jaccard_threshold = 0.7,
  msigdb_collections = c("H", "C2", "C7"),
  msigdb_subcollections = list(
    C2 = "CP:WIKIPATHWAYS",
    C7 = "IMMUNESIGDB"
  ),
  progeny_top_n = 100,
  collectri_organism = "human",
  collectri_min_targets = 0,
  collectri_max_targets = 2000,
  skip = character(0),
  output_dir = "pathway_library"
)
