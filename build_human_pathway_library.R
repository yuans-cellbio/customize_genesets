
packages <- c("build_pathway_library.R", "download_pathway_gmt.R", "build_go_regulation_gmt.R")
invisible(lapply(packages, source))

library(org.Hs.eg.db)

# # --- Build signed/unsigned BP library ---
# result_bp <- build_go_regulation_library(
#   org_db = org.Hs.eg.db,
#   include_IEA = TRUE,
#   min_size = 15,
#   max_size = 500,
#   output_dir = "go_dedup_gmt"
# )
# 
# # --- Build deduplicated CC library ---
# result_cc <- build_go_dedup_library(
#   org_db = org.Hs.eg.db,
#   ontology = "CC",
#   include_IEA = TRUE,
#   min_size = 15,
#   max_size = 500,
#   output_dir = "go_dedup_gmt"
# )
# 
# # --- Build deduplicated MF library ---
# result_mf <- build_go_dedup_library(
#   org_db = org.Hs.eg.db,
#   ontology = "MF",
#   include_IEA = TRUE,
#   min_size = 15,
#   max_size = 500,
#   output_dir = "go_dedup_gmt"
# )
# 
# # Download MSigDB, KEGG, Rectome, and PROGENy gene sets as GMT, 
# all_files <- download_all_pathway_gmt(
#   org_db = org.Hs.eg.db,
#   species = "Homo sapiens",
#   output_dir = "pathway_gmt"
# )

pathways <- build_pathway_library(
    org_db = org.Hs.eg.db,
    species = "Homo sapiens",
    kegg_organism = "hsa",
    progeny_organism = "human",
    include_IEA = TRUE,
    pos_keywords = NULL, # Use default
    neg_keywords = NULL, # Use default
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
    go_script = "build_go_regulation_gmt.R",
    pathway_script = "download_pathway_gmt.R",
    output_dir = "pathway_library"
    )

