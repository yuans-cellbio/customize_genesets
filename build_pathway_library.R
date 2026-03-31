#!/usr/bin/env Rscript
# =============================================================================
# build_pathway_library.R
#
# Master pipeline: builds a complete, version-tracked pathway library by
# combining GO-based gene sets (signed/unsigned/CC/MF) with external databases
# (MSigDB, KEGG, Reactome, PROGENy).
#
# Outputs:
#   - GMT files for all databases
#   - README.md with full provenance: database versions, package versions,
#     retrieval date, parameters, and file inventory
#
# Usage:
#   source("build_pathway_library.R")
#   result <- build_pathway_library(org_db = org.Hs.eg.db)
#
# Dependencies:
#   source("build_go_regulation_gmt.R")
#   source("download_pathway_gmt.R")
#   BiocManager::install(c("GO.db", "org.Hs.eg.db", "AnnotationDbi",
#                           "reactome.db", "progeny", "KEGGREST"))
#   install.packages(c("msigdbr", "Matrix"))
# =============================================================================

# =============================================================================
# VERSION INFO COLLECTORS
# =============================================================================

#' Collect version and metadata from installed packages and databases.
#' @param org_db OrgDb object
#' @return named list of version strings
collect_versions <- function(org_db) {
  versions <- list()
  
  # R
  versions$R <- paste(R.version$major, R.version$minor, sep = ".")
  versions$platform <- R.version$platform
  
  # GO.db
  versions$GO.db <- tryCatch({
    meta <- AnnotationDbi::metadata(GO.db)
    list(
      package_version = as.character(packageVersion("GO.db")),
      source_date = meta$value[meta$name == "GOSOURCEDATE"],
      source_url  = meta$value[meta$name == "GOSOURCEURL"]
    )
  }, error = function(e) list(package_version = "unknown"))
  
  # OrgDb
  versions$org_db <- tryCatch({
    meta <- AnnotationDbi::metadata(org_db)
    pkgName <- org_db$packageName
    list(
      package   = pkgName,
      version   = sessionInfo(pkgName)$otherPkgs[[1]]$Version,
      source    = meta$value[meta$name == "DBSCHEMA"],
      built     = meta$value[meta$name == "DBSCHEMAVERSION"]
    )
  }, error = function(e) list(package = "unknown"))
  
  # msigdbr
  versions$msigdbr <- tryCatch({
    list(
      package_version = as.character(packageVersion("msigdbr")),
      msigdb_version  = sub("^(\\d+\\.\\d+).*", "\\1", as.character(packageVersion("msigdbr")))
    )
  }, error = function(e) list(package_version = "not installed"))
  
  # KEGGREST
  versions$KEGGREST <- tryCatch({
    list(package_version = as.character(packageVersion("KEGGREST")))
  }, error = function(e) list(package_version = "not installed"))
  
  # reactome.db
  versions$reactome.db <- tryCatch({
    meta <- AnnotationDbi::metadata(reactome.db::reactome.db)
    list(
      package_version = as.character(packageVersion("reactome.db")),
      db_date         = meta$value[meta$name == "DBSCHEMAVERSION"]
    )
  }, error = function(e) list(package_version = "not installed"))
  
  # progeny
  versions$progeny <- tryCatch({
    list(package_version = as.character(packageVersion("progeny")))
  }, error = function(e) list(package_version = "not installed"))
  
  # decoupleR
  versions$decoupleR <- tryCatch({
    list(package_version = as.character(packageVersion("decoupleR")))
  }, error = function(e) list(package_version = "not installed"))
  
  # Matrix
  versions$Matrix <- tryCatch({
    as.character(packageVersion("Matrix"))
  }, error = function(e) "not installed")
  
  versions
}

# =============================================================================
# README WRITER
# =============================================================================

#' Write a README.md file documenting the library build.
#'
#' @param versions output from collect_versions()
#' @param params named list of parameters used
#' @param file_inventory named list: logical section -> character vector of file paths
#' @param output_dir directory containing the GMT files
write_readme <- function(versions, params, file_inventory, output_dir) {
  lines <- character()
  add <- function(...) lines <<- c(lines, paste0(...))
  
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
  add("|----------|----------------|-------------|")
  
  # GO.db
  go <- versions$GO.db
  add("| GO.db | ", go$package_version, " | ",
      "Source date: ", ifelse(length(go$source_date), go$source_date, "N/A"), " |")
  
  # OrgDb
  org <- versions$org_db
  add("| OrgDb (", org$package, ") | ", org$version, " | ",
      "Schema: ", ifelse(length(org$source), org$source, "N/A"), " |")
  
  # msigdbr
  ms <- versions$msigdbr
  add("| msigdbr | ", ms$package_version, " | ",
      "MSigDB v", ifelse(length(ms$msigdb_version), ms$msigdb_version, "N/A"), " |")
  
  # KEGGREST
  add("| KEGGREST | ", versions$KEGGREST$package_version, " | ",
      "Retrieved: ", format(Sys.Date(), "%Y-%m-%d"), " |")
  
  # reactome.db
  rx <- versions$reactome.db
  add("| reactome.db | ", rx$package_version, " | ",
      "DB version: ", ifelse(length(rx$db_date), rx$db_date, "N/A"), " |")
  
  # progeny
  add("| progeny | ", versions$progeny$package_version, " | - |")
  
  # decoupleR / CollecTRI
  add("| decoupleR (CollecTRI) | ", versions$decoupleR$package_version,
      " | Retrieved: ", format(Sys.Date(), "%Y-%m-%d"), " |")
  
  add("")
  add("## Parameters")
  add("")
  add("```")
  for (nm in names(params)) {
    val <- params[[nm]]
    if (is.character(val) && length(val) > 1) {
      val <- paste(val, collapse = ", ")
    }
    add(nm, " = ", val)
  }
  add("```")
  add("")
  
  add("## File Inventory")
  add("")
  
  for (section in names(file_inventory)) {
    add("### ", section)
    add("")
    files <- file_inventory[[section]]
    for (f in files) {
      if (file.exists(f)) {
        fname <- basename(f)
        # Count lines (= number of gene sets)
        n_sets <- length(readLines(f, warn = FALSE))
        size_kb <- round(file.info(f)$size / 1024, 1)
        add("- **", fname, "** — ", n_sets, " gene sets (", size_kb, " KB)")
      } else {
        add("- **", basename(f), "** — file not found")
      }
    }
    add("")
  }
  
  add("## Usage with clusterProfiler")
  add("")
  add("```r")
  add('library(clusterProfiler)')
  add('source("build_go_regulation_gmt.R")  # for read_gmt_term2name()')
  add("")
  add("# Load any GMT file")
  add('gmt       <- read.gmt("path/to/file.gmt")')
  add('gmt_names <- read_gmt_term2name("path/to/file.gmt")')
  add("")
  add("# ORA")
  add("enricher(gene = my_genes, universe = bg, TERM2GENE = gmt, TERM2NAME = gmt_names)")
  add("")
  add("# GSEA")
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

# =============================================================================
# MASTER PIPELINE
# =============================================================================

#' Build a complete, version-tracked pathway gene set library.
#'
#' Runs GO regulation pipeline (signed + unsigned BP, CC, MF) and downloads
#' KEGG, Reactome, MSigDB Hallmark, and PROGENy. Writes a README.md with
#' full provenance for reproducibility.
#'
#' @param org_db OrgDb annotation object (e.g. org.Hs.eg.db)
#' @param species species name (default "Homo sapiens")
#' @param kegg_organism KEGG organism code (default "hsa")
#' @param progeny_organism PROGENy organism (default "human")
#' @param include_IEA include electronically inferred GO annotations? (default TRUE)
#' @param pos_keywords positive-direction GO keywords
#' @param neg_keywords negative-direction GO keywords
#' @param min_size minimum gene set size (default 15)
#' @param max_size maximum gene set size (default 500)
#' @param dag_overlap_threshold DAG pruning threshold (default 0.8)
#' @param jaccard_threshold Jaccard dedup threshold (default 0.7)
#' @param msigdb_collections MSigDB collections to download (default: H and C2)
#' @param msigdb_subcollections optional subcollection filters (named list)
#' @param progeny_top_n PROGENy top genes per pathway (default 100)
#' @param collectri_organism CollecTRI organism (default "human")
#' @param collectri_min_targets CollecTRI minimum targets per TF (default 10)
#' @param collectri_max_targets CollecTRI maximum targets per TF (default 500)
#' @param skip character vector of databases to skip.
#'   Options: "go_bp", "go_cc", "go_mf", "msigdb", "kegg", "reactome", "progeny", "collectri"
#' @param go_script path to build_go_regulation_gmt.R
#' @param pathway_script path to download_pathway_gmt.R
#' @param output_dir output directory (default "pathway_library")
#' @return invisible list of all results and file paths
build_pathway_library <- function(
    org_db,
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
    msigdb_collections = c("H", "C2"),
    msigdb_subcollections = NULL,
    progeny_top_n = 100,
    collectri_organism = "human",
    collectri_min_targets = 10,
    collectri_max_targets = 500,
    skip = character(0),
    go_script = "build_go_regulation_gmt.R",
    pathway_script = "download_pathway_gmt.R",
    output_dir = "pathway_library"
) {
  # --- Source dependency scripts ---
  source(go_script, local = FALSE)
  source(pathway_script, local = FALSE)
  
  # Use defaults from the GO script if not provided
  if (is.null(pos_keywords)) pos_keywords <- DEFAULT_POS_KEYWORDS
  if (is.null(neg_keywords)) neg_keywords <- DEFAULT_NEG_KEYWORDS
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Track all output files
  file_inventory <- list()
  results <- list()
  
  # =========================================================================
  # GO BP: signed + unsigned
  # =========================================================================
  if (!"go_bp" %in% skip) {
    message("\n##########  GO Biological Process (signed + unsigned)  ##########\n")
    results$go_bp <- build_go_regulation_library(
      org_db = org_db,
      include_IEA = include_IEA,
      pos_keywords = pos_keywords,
      neg_keywords = neg_keywords,
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
  
  # =========================================================================
  # GO CC
  # =========================================================================
  if (!"go_cc" %in% skip) {
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
  
  # =========================================================================
  # GO MF
  # =========================================================================
  if (!"go_mf" %in% skip) {
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
  
  # =========================================================================
  # MSigDB
  # =========================================================================
  if (!"msigdb" %in% skip) {
    message("\n##########  MSigDB  ##########\n")
    results$msigdb <- download_msigdb_gmt(
      species = species,
      collections = msigdb_collections,
      subcollections = msigdb_subcollections,
      output_dir = output_dir
    )
    file_inventory[["MSigDB"]] <- unlist(results$msigdb)
  }
  
  # =========================================================================
  # KEGG
  # =========================================================================
  if (!"kegg" %in% skip) {
    message("\n##########  KEGG  ##########\n")
    results$kegg <- download_kegg_gmt(
      organism = kegg_organism,
      output_dir = output_dir
    )
    file_inventory[["KEGG"]] <- results$kegg
  }
  
  # =========================================================================
  # Reactome
  # =========================================================================
  if (!"reactome" %in% skip) {
    message("\n##########  Reactome  ##########\n")
    results$reactome <- download_reactome_gmt(
      species = species,
      method = "reactome.db",
      output_dir = output_dir
    )
    file_inventory[["Reactome"]] <- results$reactome
  }
  
  # =========================================================================
  # PROGENy
  # =========================================================================
  if (!"progeny" %in% skip) {
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
  
  # =========================================================================
  # CollecTRI
  # =========================================================================
  if (!"collectri" %in% skip) {
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
  
  # =========================================================================
  # VERSION INFO & README
  # =========================================================================
  message("\n##########  Writing README  ##########\n")
  
  versions <- collect_versions(org_db)
  
  params <- list(
    species = species,
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
    skipped = if (length(skip) > 0) skip else "none"
  )
  
  readme_path <- write_readme(versions, params, file_inventory, output_dir)
  
  # =========================================================================
  # SUMMARY
  # =========================================================================
  total_files <- sum(vapply(file_inventory, length, integer(1)))
  message(sprintf(
    "\n========== COMPLETE ==========\n%d GMT files written to %s\nREADME: %s",
    total_files, output_dir, readme_path
  ))
  
  results$file_inventory <- file_inventory
  results$versions <- versions
  results$readme <- readme_path
  
  invisible(results)
}

# =============================================================================
# RUN EXAMPLE
# =============================================================================

if (sys.nframe() == 0) {
  library(org.Hs.eg.db)
  
  result <- build_pathway_library(
    org_db = org.Hs.eg.db,
    species = "Homo sapiens",
    output_dir = "pathway_library"
  )
  
  # Skip slow databases for quick testing:
  # result <- build_pathway_library(
  #   org_db = org.Hs.eg.db,
  #   skip = c("kegg"),  # KEGG requires internet and is slow
  #   output_dir = "pathway_library"
  # )
  
  # Mouse example:
  # library(org.Mm.eg.db)
  # result_mouse <- build_pathway_library(
  #   org_db = org.Mm.eg.db,
  #   species = "Mus musculus",
  #   kegg_organism = "mmu",
  #   progeny_organism = "mouse",
  #   output_dir = "pathway_library_mouse"
  # )
}