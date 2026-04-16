collect_versions <- function(org_db) {
  versions <- list()

  versions$R <- paste(R.version$major, R.version$minor, sep = ".")
  versions$platform <- R.version$platform

  versions$GO.db <- tryCatch({
    meta <- AnnotationDbi::metadata(GO.db::GO.db)
    list(
      package_version = safe_package_version("GO.db"),
      source_date = meta$value[meta$name == "GOSOURCEDATE"],
      source_url = meta$value[meta$name == "GOSOURCEURL"]
    )
  }, error = function(e) list(package_version = "unknown"))

  versions$org_db <- tryCatch({
    meta <- AnnotationDbi::metadata(org_db)
    package_name <- org_db$packageName
    list(
      package = package_name,
      version = safe_package_version(package_name),
      source = meta$value[meta$name == "DBSCHEMA"],
      built = meta$value[meta$name == "DBSCHEMAVERSION"]
    )
  }, error = function(e) list(package = "unknown", version = "unknown"))

  versions$msigdbr <- list(
    package_version = safe_package_version("msigdbr"),
    msigdb_version = tryCatch(
      sub("^(\\d+\\.\\d+).*", "\\1", safe_package_version("msigdbr")),
      error = function(e) "N/A"
    )
  )

  versions$KEGGREST <- list(package_version = safe_package_version("KEGGREST"))

  versions$reactome.db <- tryCatch({
    meta <- AnnotationDbi::metadata(reactome.db::reactome.db)
    list(
      package_version = safe_package_version("reactome.db"),
      db_date = meta$value[meta$name == "DBSCHEMAVERSION"]
    )
  }, error = function(e) list(package_version = "not installed", db_date = "N/A"))

  versions$progeny <- list(package_version = safe_package_version("progeny"))
  versions$decoupleR <- list(package_version = safe_package_version("decoupleR"))
  versions$Matrix <- safe_package_version("Matrix")

  versions
}

write_readme <- function(versions, params, file_inventory, output_dir) {
  lines <- character()
  add <- function(...) {
    lines <<- c(lines, paste0(...))
  }

  add("# Pathway Gene Set Library")
  add("")
  add("Auto-generated pathway library for enrichment analysis (ORA / GSEA).")
  add("")
  add("## Build Information")
  add("")
  add("| Field | Value |")
  add("|-------|-------|")
  add("| Date | ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), " |")
  add("| R version | ", versions$R, " |")
  add("| Platform | ", versions$platform, " |")
  add("")

  add("## Database Versions")
  add("")
  add("| Database | Package Version | Source Info |")
  add("|----------|-----------------|-------------|")
  add("| GO.db | ", versions$GO.db$package_version, " | Source date: ",
      ifelse(length(versions$GO.db$source_date), versions$GO.db$source_date, "N/A"), " |")
  add("| OrgDb (", versions$org_db$package, ") | ", versions$org_db$version, " | Schema: ",
      ifelse(length(versions$org_db$source), versions$org_db$source, "N/A"), " |")
  add("| msigdbr | ", versions$msigdbr$package_version, " | MSigDB v",
      ifelse(length(versions$msigdbr$msigdb_version), versions$msigdbr$msigdb_version, "N/A"), " |")
  add("| KEGGREST | ", versions$KEGGREST$package_version, " | Retrieved: ", format(Sys.Date(), "%Y-%m-%d"), " |")
  add("| reactome.db | ", versions$reactome.db$package_version, " | DB version: ",
      ifelse(length(versions$reactome.db$db_date), versions$reactome.db$db_date, "N/A"), " |")
  add("| progeny | ", versions$progeny$package_version, " | - |")
  add("| decoupleR (CollecTRI) | ", versions$decoupleR$package_version,
      " | Retrieved: ", format(Sys.Date(), "%Y-%m-%d"), " |")
  add("")

  add("## Methodology Notes")
  add("")
  add("- GO BP signed terms are classified from GO term names using explicit positive/negative keyword lists.")
  add("- Sibling pairing is conservative: only explicit text replacements (for example positive/negative regulation) are used when forcing pair retention.")
  add("- GO libraries are filtered by gene-set size, pruned by GO DAG ancestor overlap, and deduplicated by Jaccard overlap while preferring more specific GO terms.")
  add("- MSigDB uses mouse-native collections (`db_species = MM`) for mouse builds and human collections with ortholog mapping for rat builds.")
  add("- External libraries are exported as Entrez-based GMT files for direct use with clusterProfiler-style workflows.")
  add("")

  add("## Parameters")
  add("")
  add("```")
  for (param_name in names(params)) {
    value <- params[[param_name]]
    if (is.null(value) || !length(value)) {
      value <- "not used"
    } else if (is.character(value) && length(value) > 1L) {
      value <- paste(value, collapse = ", ")
    }
    add(param_name, " = ", value)
  }
  add("```")
  add("")

  add("## File Inventory")
  add("")
  for (section_name in names(file_inventory)) {
    add("### ", section_name)
    add("")
    for (filepath in file_inventory[[section_name]]) {
      if (!file.exists(filepath)) {
        add("- **", basename(filepath), "** - file not found")
        next
      }

      n_sets <- length(readLines(filepath, warn = FALSE))
      size_kb <- round(file.info(filepath)$size / 1024, 1)
      add("- **", basename(filepath), "** - ", n_sets, " gene sets (", size_kb, " KB)")
    }
    add("")
  }

  add("## Usage with clusterProfiler")
  add("")
  add("```r")
  add('library(clusterProfiler)')
  add('library(customizeGeneSets)')
  add("")
  add('gmt <- read.gmt("path/to/file.gmt")')
  add('gmt_names <- read_gmt_term2name("path/to/file.gmt")')
  add("")
  add("enricher(gene = my_genes, universe = bg, TERM2GENE = gmt, TERM2NAME = gmt_names)")
  add("GSEA(geneList = my_ranked_genes, TERM2GENE = gmt, TERM2NAME = gmt_names)")
  add("```")
  add("")

  add("## Gene ID Format")
  add("")
  add("All GMT files use **NCBI Entrez Gene IDs** for consistency.")
  add("Convert with `clusterProfiler::bitr()` if your data uses Ensembl or symbols.")

  readme_path <- file.path(output_dir, "README.md")
  writeLines(lines, readme_path)
  message(sprintf("Wrote %s", readme_path))
  readme_path
}

build_pathway_library <- function(
    org_db = NULL,
    species = "human",
    include_IEA = TRUE,
    pos_keywords = NULL,
    neg_keywords = NULL,
    pair_replacements = DEFAULT_DIRECTIONAL_PAIR_REPLACEMENTS,
    min_size = 15,
    max_size = 500,
    dag_overlap_threshold = 0.8,
    jaccard_threshold = 0.7,
    msigdb_collections = c("H", "C2"),
    msigdb_subcollections = NULL,
    progeny_top_n = 100,
    collectri_min_targets = 10,
    collectri_max_targets = 500,
    skip = character(0),
    go_script = "build_go_regulation_gmt.R",
    pathway_script = "download_pathway_gmt.R",
    output_dir = "pathway_library",
    ...
) {
  reject_database_species_overrides(list(...), fun_name = "build_pathway_library")

  if (!exists("build_go_regulation_library", mode = "function") && file.exists(go_script)) {
    source(go_script, local = FALSE)
  }
  if (!exists("download_msigdb_gmt", mode = "function") && file.exists(pathway_script)) {
    source(pathway_script, local = FALSE)
  }

  assert_positive_count(min_size, "min_size")
  assert_positive_count(max_size, "max_size")
  assert_min_max(min_size, max_size, "min_size", "max_size")
  assert_probability(dag_overlap_threshold, "dag_overlap_threshold")
  assert_probability(jaccard_threshold, "jaccard_threshold")
  assert_positive_count(progeny_top_n, "progeny_top_n")
  assert_positive_count(collectri_min_targets, "collectri_min_targets", allow_zero = TRUE)
  assert_positive_count(collectri_max_targets, "collectri_max_targets")
  assert_min_max(collectri_min_targets, collectri_max_targets, "collectri_min_targets", "collectri_max_targets")

  if (is.null(pos_keywords)) {
    pos_keywords <- DEFAULT_POS_KEYWORDS
  }
  if (is.null(neg_keywords)) {
    neg_keywords <- DEFAULT_NEG_KEYWORDS
  }

  harmonized <- resolve_species_build_context(
    species = species,
    org_db = org_db,
    skip = skip
  )
  org_db <- harmonized$org_db
  species <- harmonized$species
  msigdb_db_species <- harmonized$msigdb_db_species
  kegg_organism <- harmonized$kegg_organism
  progeny_organism <- harmonized$progeny_organism
  collectri_organism <- harmonized$collectri_organism
  skip <- harmonized$skip

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  file_inventory <- list()
  results <- list()

  if (!("go_bp" %in% skip)) {
    message("\n##########  GO Biological Process (signed + unsigned)  ##########\n")
    results$go_bp <- build_go_regulation_library(
      org_db = org_db,
      include_IEA = include_IEA,
      pos_keywords = pos_keywords,
      neg_keywords = neg_keywords,
      pair_replacements = pair_replacements,
      min_size = min_size,
      max_size = max_size,
      dag_overlap_threshold = dag_overlap_threshold,
      jaccard_threshold = jaccard_threshold,
      output_dir = output_dir
    )
    file_inventory[["GO Biological Process"]] <- c(
      results$go_bp$signed_gmt,
      results$go_bp$unsigned_gmt
    )
  }

  if (!("go_cc" %in% skip)) {
    message("\n##########  GO Cellular Component  ##########\n")
    results$go_cc <- build_go_dedup_library(
      org_db = org_db,
      ontology = "CC",
      include_IEA = include_IEA,
      min_size = min_size,
      max_size = max_size,
      dag_overlap_threshold = dag_overlap_threshold,
      jaccard_threshold = jaccard_threshold,
      output_dir = output_dir
    )
    file_inventory[["GO Cellular Component"]] <- results$go_cc$gmt_file
  }

  if (!("go_mf" %in% skip)) {
    message("\n##########  GO Molecular Function  ##########\n")
    results$go_mf <- build_go_dedup_library(
      org_db = org_db,
      ontology = "MF",
      include_IEA = include_IEA,
      min_size = min_size,
      max_size = max_size,
      dag_overlap_threshold = dag_overlap_threshold,
      jaccard_threshold = jaccard_threshold,
      output_dir = output_dir
    )
    file_inventory[["GO Molecular Function"]] <- results$go_mf$gmt_file
  }

  if (!("msigdb" %in% skip)) {
    message("\n##########  MSigDB  ##########\n")
    results$msigdb <- download_msigdb_gmt(
      species = species,
      db_species = msigdb_db_species,
      collections = msigdb_collections,
      subcollections = msigdb_subcollections,
      output_dir = output_dir
    )
    file_inventory[["MSigDB"]] <- unlist(results$msigdb, use.names = FALSE)
  }

  if (!("kegg" %in% skip)) {
    message("\n##########  KEGG  ##########\n")
    results$kegg <- download_kegg_gmt(
      organism = kegg_organism,
      output_dir = output_dir
    )
    file_inventory[["KEGG"]] <- results$kegg
  }

  if (!("reactome" %in% skip)) {
    message("\n##########  Reactome  ##########\n")
    results$reactome <- download_reactome_gmt(
      species = species,
      method = "reactome.db",
      output_dir = output_dir
    )
    file_inventory[["Reactome"]] <- results$reactome
  }

  if (!("progeny" %in% skip)) {
    message("\n##########  PROGENy  ##########\n")
    results$progeny <- download_progeny_gmt(
      org_db = org_db,
      organism = progeny_organism,
      top_n = progeny_top_n,
      output_dir = output_dir
    )
    file_inventory[["PROGENy"]] <- c(
      results$progeny$signed_gmt,
      results$progeny$unsigned_gmt
    )
  }

  if (!("collectri" %in% skip)) {
    message("\n##########  CollecTRI  ##########\n")
    results$collectri <- download_collectri_gmt(
      org_db = org_db,
      organism = collectri_organism,
      min_targets = collectri_min_targets,
      max_targets = collectri_max_targets,
      output_dir = output_dir
    )
    file_inventory[["CollecTRI"]] <- c(
      results$collectri$signed_gmt,
      results$collectri$unsigned_gmt
    )
  }

  message("\n##########  Writing README  ##########\n")
  versions <- collect_versions(org_db)
  params <- list(
    species = species,
    msigdb_db_species = msigdb_db_species,
    kegg_organism = kegg_organism,
    progeny_organism = progeny_organism,
    collectri_organism = collectri_organism,
    include_IEA = include_IEA,
    min_size = min_size,
    max_size = max_size,
    dag_overlap_threshold = dag_overlap_threshold,
    jaccard_threshold = jaccard_threshold,
    msigdb_collections = msigdb_collections,
    progeny_top_n = progeny_top_n,
    collectri_min_targets = collectri_min_targets,
    collectri_max_targets = collectri_max_targets,
    pos_keywords = pos_keywords,
    neg_keywords = neg_keywords,
    skipped = if (length(skip)) skip else "none"
  )
  readme_path <- write_readme(versions, params, file_inventory, output_dir)

  total_files <- sum(vapply(file_inventory, length, integer(1)))
  message(sprintf(
    "\n========== COMPLETE ==========\n%d GMT files written to %s\nREADME: %s",
    total_files,
    output_dir,
    readme_path
  ))

  results$file_inventory <- file_inventory
  results$versions <- versions
  results$readme <- readme_path
  invisible(results)
}
