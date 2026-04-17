internal_test_functions <- c(
  "build_directional_pair_map",
  "canonical_species_name",
  "dedup_by_jaccard",
  "enforce_paired_survival",
  "is_human_species",
  "is_mouse_species",
  "is_rat_species",
  "map_symbols_to_entrez",
  "normalize_collectri_direct_query",
  "organism_name_from_orgdb",
  "prepare_gsea_rankings",
  "resolve_ora_gene_input",
  "resolve_msigdb_db_species",
  "resolve_species_build_context",
  "select_ranked_genes_for_ora",
  "translate_msigdb_collection_code",
  "validate_signed_rank",
  "write_gmt_from_list"
)

if ("customizeGeneSets" %in% loadedNamespaces()) {
  pkg_ns <- asNamespace("customizeGeneSets")
  helper_env <- environment()

  for (function_name in internal_test_functions) {
    if (!exists(function_name, envir = helper_env, inherits = FALSE) &&
        exists(function_name, envir = pkg_ns, inherits = FALSE)) {
      assign(function_name, get(function_name, envir = pkg_ns), envir = helper_env)
    }
  }
}
