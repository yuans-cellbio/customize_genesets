#!/usr/bin/env Rscript
# =============================================================================
# download_pathway_gmt.R
#
# Download pathway databases as local GMT files (Entrez IDs) for use with
# clusterProfiler::enricher() / GSEA() and fgsea.
#
# Supported databases:
#   - MSigDB (all collections or specific ones, via msigdbr)
#   - KEGG (via KEGGREST, requires internet)
#   - Reactome (via reactome.db or msigdbr)
#   - PROGENy (via progeny, signed weights -> split into pos/neg GMT)
#
# Dependencies:
#   install.packages(c("msigdbr", "KEGGREST"))
#   BiocManager::install(c("reactome.db", "AnnotationDbi", "progeny"))
#   # org.Hs.eg.db or equivalent for KEGG symbol-to-entrez conversion
#
# All outputs use Entrez IDs for consistency with GO GMT files.
# =============================================================================

# =============================================================================
# SHARED HELPERS
# =============================================================================

#' Write a named list of gene sets to GMT format
#' @param gene_sets named list: set name -> character vector of Entrez IDs
#' @param descriptions named character vector: set name -> description (optional)
#' @param filepath output file path
write_gmt_from_list <- function(gene_sets, descriptions = NULL, filepath) {
  lines <- vapply(names(gene_sets), function(name) {
    desc <- if (!is.null(descriptions) && name %in% names(descriptions)) {
      descriptions[[name]]
    } else {
      "NA"
    }
    genes <- gene_sets[[name]]
    paste(c(name, desc, genes), collapse = "\t")
  }, character(1))

  writeLines(lines, filepath)
  message(sprintf("Wrote %d gene sets to %s", length(gene_sets), filepath))
}

#' Read a TERM2NAME mapping from a GMT file (same as in build_go_regulation_gmt.R)
#' @param filepath path to GMT file
#' @return data.frame with columns: ID, Name
read_gmt_term2name <- function(filepath) {
  lines <- readLines(filepath)
  parts <- strsplit(lines, "\t", fixed = TRUE)
  data.frame(
    ID   = vapply(parts, `[`, character(1), 1),
    Name = vapply(parts, `[`, character(1), 2),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 1. MSigDB
# =============================================================================

#' Download MSigDB gene sets as GMT files.
#'
#' Uses the msigdbr package (data is bundled, no internet required).
#' Can export all collections or specific ones.
#'
#' @param species species name (default "Homo sapiens"). See msigdbr::msigdbr_species()
#' @param collections character vector of collection codes to export.
#'   Use "all" for everything, or specific codes like c("H", "C2", "C5").
#'   Common collections:
#'     H  - Hallmark (50 sets, high-level biological states)
#'     C2 - Curated (includes KEGG, Reactome, WikiPathways, BioCarta, CGP)
#'     C3 - Regulatory target (TF targets, miRNA targets)
#'     C5 - Ontology (GO BP, CC, MF, HPO)
#'     C6 - Oncogenic signatures
#'     C7 - Immunologic signatures
#'     C8 - Cell type signatures
#' @param subcollections optional named list mapping collection -> subcollection filter.
#'   e.g. list(C2 = "CP:KEGG_MEDICUS") to get only KEGG from C2.
#'   See msigdbr::msigdbr_collections() for available subcollections.
#' @param output_dir directory for output files
#' @return invisible list of file paths
download_msigdb_gmt <- function(
  species = "Homo sapiens",
  collections = "all",
  subcollections = NULL,
  output_dir = "."
) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Install msigdbr: install.packages('msigdbr')")
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  species_tag <- gsub(" ", "_", tolower(species))

  # Get available collections
  avail <- msigdbr::msigdbr_collections()

  if (identical(collections, "all")) {
    collections <- unique(avail$gs_collection)
  }

  files <- list()

  for (coll in collections) {
    subcoll_filter <- subcollections[[coll]]

    if (!is.null(subcoll_filter)) {
      m_df <- msigdbr::msigdbr(species = species, collection = coll, subcollection = subcoll_filter)
      file_tag <- paste0(tolower(coll), "_", tolower(gsub("[: ]", "_", subcoll_filter)))
    } else {
      m_df <- msigdbr::msigdbr(species = species, collection = coll)
      file_tag <- tolower(coll)
    }

    if (nrow(m_df) == 0) {
      message(sprintf("  Skipping %s: no gene sets found", coll))
      next
    }

    # Convert to named list using Entrez IDs
    # Remove NAs and empty strings
    m_df <- m_df[!is.na(m_df$ncbi_gene) & m_df$ncbi_gene != "", ]

    gene_sets <- split(m_df$ncbi_gene, m_df$gs_name)
    gene_sets <- lapply(gene_sets, unique)

    # Build descriptions
    desc_df <- unique(m_df[, c("gs_name", "gs_description")])
    descriptions <- setNames(desc_df$gs_description, desc_df$gs_name)

    filepath <- file.path(output_dir, sprintf("msigdb_%s_%s.gmt", file_tag, species_tag))
    write_gmt_from_list(gene_sets, descriptions, filepath)
    files[[coll]] <- filepath
  }

  invisible(files)
}

# =============================================================================
# 2. KEGG
# =============================================================================

#' Download KEGG pathway gene sets as a GMT file.
#'
#' Uses KEGGREST to fetch current pathway-gene mappings. Requires internet.
#'
#' @param organism KEGG organism code (default "hsa" for human).
#'   Common codes: "hsa" (human), "mmu" (mouse), "rno" (rat).
#' @param output_dir directory for output file
#' @return invisible file path
download_kegg_gmt <- function(
  organism = "hsa",
  output_dir = "."
) {
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    stop("Install KEGGREST: BiocManager::install('KEGGREST')")
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Get list of pathways
  message("Fetching KEGG pathway list...")
  pathways <- KEGGREST::keggList("pathway", organism)
  pathway_ids <- sub(paste0("^path:", organism), "", names(pathways))
  pathway_names <- sub(paste0(" - .*$"), "", pathways)

  message(sprintf("  Found %d pathways", length(pathway_ids)))

  # Fetch gene members for each pathway
  message("Fetching gene members (this may take a few minutes)...")
  gene_sets <- list()
  descriptions <- character()

  for (i in seq_along(pathway_ids)) {
    pid <- pathway_ids[i]

    tryCatch({
      genes <- KEGGREST::keggGet(paste0("path:", pid))[[1]]$GENE

      if (!is.null(genes)) {
        # KEGG returns alternating: entrez_id, "gene_symbol; description"
        # Extract only the Entrez IDs (odd positions)
        entrez_ids <- genes[seq(1, length(genes), by = 2)]
        gene_sets[[pid]] <- unique(entrez_ids)
        descriptions[[pid]] <- pathway_names[i]
      }
    }, error = function(e) {
      message(sprintf("    Warning: failed to fetch %s", pid))
    })

    # Rate limiting — KEGG requires reasonable request rates
    if (i %% 10 == 0) {
      message(sprintf("    Processed %d / %d pathways", i, length(pathway_ids)))
      Sys.sleep(1)
    }
  }

  # Remove empty sets
  gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) > 0]
  message(sprintf("  %d pathways with gene annotations", length(gene_sets)))

  filepath <- file.path(output_dir, sprintf("kegg_%s.gmt", organism))
  write_gmt_from_list(gene_sets, descriptions, filepath)

  invisible(filepath)
}

# =============================================================================
# 3. REACTOME
# =============================================================================

#' Download Reactome pathway gene sets as a GMT file.
#'
#' Two methods available:
#'   method = "reactome.db" — uses the Bioconductor reactome.db package (offline)
#'   method = "msigdbr"     — extracts Reactome from MSigDB C2:CP:REACTOME
#'
#' @param species species name for msigdbr method (default "Homo sapiens")
#' @param method one of "reactome.db" or "msigdbr" (default "reactome.db")
#' @param output_dir directory for output file
#' @return invisible file path
download_reactome_gmt <- function(
  species = "Homo sapiens",
  method = "reactome.db",
  output_dir = "."
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (method == "msigdbr") {
    # --- Method: extract from MSigDB ---
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
      stop("Install msigdbr: install.packages('msigdbr')")
    }

    message("Fetching Reactome from MSigDB C2:CP:REACTOME...")
    m_df <- msigdbr::msigdbr(species = species, collection = "C2", subcollection = "CP:REACTOME")
    m_df <- m_df[!is.na(m_df$ncbi_gene) & m_df$ncbi_gene != "", ]

    gene_sets <- split(m_df$ncbi_gene, m_df$gs_name)
    gene_sets <- lapply(gene_sets, unique)

    desc_df <- unique(m_df[, c("gs_name", "gs_description")])
    descriptions <- setNames(desc_df$gs_description, desc_df$gs_name)

    species_tag <- gsub(" ", "_", tolower(species))
    filepath <- file.path(output_dir, sprintf("reactome_msigdb_%s.gmt", species_tag))

  } else if (method == "reactome.db") {
    # --- Method: reactome.db ---
    if (!requireNamespace("reactome.db", quietly = TRUE)) {
      stop("Install reactome.db: BiocManager::install('reactome.db')")
    }

    message("Fetching Reactome pathways from reactome.db...")

    # Get pathway-to-Entrez mappings
    pathway_genes <- AnnotationDbi::select(reactome.db::reactome.db,
      keys = AnnotationDbi::keys(reactome.db::reactome.db, keytype = "PATHID"),
      keytype = "PATHID",
      columns = c("PATHID", "ENTREZID")
    )
    pathway_genes <- pathway_genes[!is.na(pathway_genes$ENTREZID), ]

    # Get pathway names
    pathway_names_df <- AnnotationDbi::select(reactome.db::reactome.db,
      keys = unique(pathway_genes$PATHID),
      keytype = "PATHID",
      columns = c("PATHID", "PATHNAME")
    )

    # Filter to species of interest
    species_prefix <- sub("^(\\S+).*", "\\1", species)
    pathway_names_df <- pathway_names_df[
      grepl(paste0("^", species_prefix), pathway_names_df$PATHNAME, ignore.case = TRUE), ]

    # Keep only pathways for this species
    keep_ids <- unique(pathway_names_df$PATHID)
    pathway_genes <- pathway_genes[pathway_genes$PATHID %in% keep_ids, ]

    gene_sets <- split(pathway_genes$ENTREZID, pathway_genes$PATHID)
    gene_sets <- lapply(gene_sets, unique)

    descriptions <- setNames(pathway_names_df$PATHNAME, pathway_names_df$PATHID)
    # Remove species prefix from names for cleaner output
    descriptions <- sub(paste0("^", species, ": "), "", descriptions)

    species_tag <- gsub(" ", "_", tolower(species))
    filepath <- file.path(output_dir, sprintf("reactome_%s.gmt", species_tag))

  } else {
    stop("method must be one of: 'reactome.db', 'msigdbr'")
  }

  message(sprintf("  %d Reactome pathways", length(gene_sets)))
  write_gmt_from_list(gene_sets, descriptions, filepath)

  invisible(filepath)
}

# =============================================================================
# 4. PROGENy (SIGNED PATHWAY SIGNATURES)
# =============================================================================

#' Export PROGENy pathway signatures as GMT files.
#'
#' PROGENy assigns signed weights to genes per pathway. This function exports
#' two GMT files:
#'   1. Separate signed sets: "PATHWAY_UP" and "PATHWAY_DN" for each pathway
#'   2. Combined unsigned sets: all footprint genes per pathway regardless of sign
#'
#' Genes are filtered to the top N most significant per pathway.
#'
#' @param org_db an OrgDb object for symbol-to-Entrez conversion
#' @param organism one of "human" or "mouse" (default "human")
#' @param top_n number of top significant genes per pathway (default 100)
#' @param output_dir directory for output files
#' @return invisible list with signed_gmt and unsigned_gmt file paths
download_progeny_gmt <- function(
    org_db,
    organism = "human",
    top_n = 100,
    output_dir = "."
) {
  if (!requireNamespace("progeny", quietly = TRUE)) {
    stop("Install progeny: BiocManager::install('progeny')")
  }
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get full PROGENy model
  model <- if (organism == "human") {
    progeny::model_human_full
  } else if (organism == "mouse") {
    progeny::model_mouse_full
  } else {
    stop("organism must be 'human' or 'mouse'")
  }
  
  message(sprintf("PROGENy %s model: %d genes x %d pathways",
                  organism, length(unique(model$gene)), length(unique(model$pathway))
  ))
  
  # Filter to top N most significant genes per pathway
  model_top <- do.call(rbind, lapply(split(model, model$pathway), function(df) {
    rownames(df) <- NULL
    df <- df[order(df$p.value), ]
    head(df, top_n)
  }))
  rownames(model_top) <- NULL
  
  # Guard against factor columns from do.call(rbind, ...) or the model itself
  model_top$gene    <- as.character(model_top$gene)
  model_top$pathway <- as.character(model_top$pathway)
  
  message(sprintf("  Using top %d genes per pathway (%d total gene-pathway pairs)",
                  top_n, nrow(model_top)
  ))
  
  # Convert gene symbols to Entrez IDs
  symbols <- unique(model_top$gene)
  message(sprintf("  Converting %d unique gene symbols to Entrez IDs...", length(symbols)))
  message(sprintf("  Example symbols: %s", paste(head(symbols, 10), collapse = ", ")))
  
  # mapIds handles missing/invalid keys gracefully (returns NA per key)
  sym_to_entrez <- AnnotationDbi::mapIds(
    org_db,
    keys = symbols,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  n_mapped <- sum(!is.na(sym_to_entrez))
  if (n_mapped == 0) {
    # Fall back: try ALIAS keytype (some PROGENy genes use older aliases)
    message("  No matches via SYMBOL, trying ALIAS...")
    sym_to_entrez <- AnnotationDbi::mapIds(
      org_db,
      keys = symbols,
      keytype = "ALIAS",
      column = "ENTREZID",
      multiVals = "first"
    )
    n_mapped <- sum(!is.na(sym_to_entrez))
  }
  
  message(sprintf("  Mapped %d / %d gene symbols to Entrez IDs", n_mapped, length(symbols)))
  
  # Build signed gene sets
  signed_sets <- list()
  unsigned_sets <- list()
  
  for (pw in unique(model_top$pathway)) {
    pw_genes <- model_top[model_top$pathway == pw, ]
    
    # Ensure character indexing into the named vector (factors use integer codes)
    pos_symbols <- as.character(pw_genes$gene[pw_genes$weight > 0])
    neg_symbols <- as.character(pw_genes$gene[pw_genes$weight < 0])
    
    pos_entrez <- na.omit(sym_to_entrez[pos_symbols])
    neg_entrez <- na.omit(sym_to_entrez[neg_symbols])
    all_entrez <- na.omit(sym_to_entrez[as.character(pw_genes$gene)])
    
    if (length(pos_entrez) > 0) {
      signed_sets[[paste0(pw, "_UP")]] <- unname(pos_entrez)
    }
    if (length(neg_entrez) > 0) {
      signed_sets[[paste0(pw, "_DN")]] <- unname(neg_entrez)
    }
    if (length(all_entrez) > 0) {
      unsigned_sets[[pw]] <- unname(all_entrez)
    }
  }
  
  # Write signed GMT
  signed_desc <- setNames(
    sub("_UP$", " upregulated genes",
        sub("_DN$", " downregulated genes", names(signed_sets))),
    names(signed_sets)
  )
  signed_file <- file.path(output_dir, sprintf("progeny_signed_%s_top%d.gmt", organism, top_n))
  write_gmt_from_list(signed_sets, signed_desc, signed_file)
  
  # Write unsigned GMT
  unsigned_desc <- setNames(
    paste(names(unsigned_sets), "regulated genes"),
    names(unsigned_sets)
  )
  unsigned_file <- file.path(output_dir, sprintf("progeny_unsigned_%s_top%d.gmt", organism, top_n))
  write_gmt_from_list(unsigned_sets, unsigned_desc, unsigned_file)
  
  message(sprintf(
    "\n=== PROGENy Summary ===\nOrganism:  %s\nTop genes: %d per pathway\nSigned:    %d sets in %s\nUnsigned:  %d sets in %s",
    organism, top_n,
    length(signed_sets), signed_file,
    length(unsigned_sets), unsigned_file
  ))
  
  invisible(list(
    signed_gmt = signed_file,
    unsigned_gmt = unsigned_file
  ))
}

# =============================================================================
# 5. CollecTRI (SIGNED TF REGULONS)
# =============================================================================

#' Export CollecTRI transcription factor regulons as GMT files.
#'
#' CollecTRI (via decoupleR/OmnipathR) provides signed TF-target interactions.
#' This function exports:
#'   1. Signed GMT: separate "TF_ACT" (activated targets) and "TF_REP"
#'      (repressed targets) sets per TF
#'   2. Unsigned GMT: all targets per TF regardless of mode
#'
#' Gene symbols are converted to Entrez IDs for consistency.
#'
#' @param org_db an OrgDb object for symbol-to-Entrez conversion
#' @param organism one of "human", "mouse", or "rat" (default "human")
#' @param split_complexes logical; split TF complexes into subunits? (default FALSE)
#' @param min_targets minimum number of targets per TF to include (default 10)
#' @param max_targets maximum number of targets per TF to include (default 500)
#' @param output_dir directory for output files
#' @return invisible list with signed_gmt and unsigned_gmt file paths
download_collectri_gmt <- function(
    org_db,
    organism = "human",
    split_complexes = FALSE,
    min_targets = 10,
    max_targets = 500,
    output_dir = "."
) {
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop("Install decoupleR: BiocManager::install('decoupleR')")
  }
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Retrieve CollecTRI network
  message(sprintf("Fetching CollecTRI network (organism: %s)...", organism))
  net <- tryCatch(
    decoupleR::get_collectri(
      organism = organism,
      split_complexes = split_complexes
    ),
    error = function(e) {
      stop(
        "Failed to retrieve CollecTRI. This requires internet access to ",
        "query OmniPath.\nError: ", e$message
      )
    }
  )
  
  n_tf <- length(unique(net$source))
  n_interactions <- nrow(net)
  message(sprintf("  %d TFs, %d interactions", n_tf, n_interactions))
  
  # Convert target gene symbols to Entrez IDs
  symbols <- unique(net$target)
  message(sprintf("  Converting %d unique target symbols to Entrez IDs...", length(symbols)))
  
  sym_to_entrez <- AnnotationDbi::mapIds(
    org_db,
    keys = symbols,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  n_mapped <- sum(!is.na(sym_to_entrez))
  if (n_mapped == 0) {
    message("  No matches via SYMBOL, trying ALIAS...")
    sym_to_entrez <- AnnotationDbi::mapIds(
      org_db,
      keys = symbols,
      keytype = "ALIAS",
      column = "ENTREZID",
      multiVals = "first"
    )
    n_mapped <- sum(!is.na(sym_to_entrez))
  }
  
  message(sprintf("  Mapped %d / %d target symbols to Entrez IDs", n_mapped, length(symbols)))
  
  if (n_mapped == 0) {
    warning("No gene symbols could be mapped. Check org_db matches the organism.")
    return(invisible(list(signed_gmt = NULL, unsigned_gmt = NULL)))
  }
  
  # Add Entrez IDs to network
  net$entrez <- sym_to_entrez[net$target]
  net <- net[!is.na(net$entrez), ]
  
  # Also convert TF source symbols for the description
  tf_symbols <- unique(net$source)
  tf_to_entrez <- AnnotationDbi::mapIds(
    org_db,
    keys = tf_symbols,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  # Build signed gene sets
  signed_sets <- list()
  signed_desc <- character()
  unsigned_sets <- list()
  unsigned_desc <- character()
  
  for (tf in unique(net$source)) {
    tf_net <- net[net$source == tf, ]
    
    act_targets <- unique(tf_net$entrez[tf_net$mor > 0])
    rep_targets <- unique(tf_net$entrez[tf_net$mor < 0])
    all_targets <- unique(tf_net$entrez)
    
    # TF Entrez ID for the GMT ID field (fall back to symbol if unmapped)
    tf_entrez <- tf_to_entrez[tf]
    tf_id <- ifelse(is.na(tf_entrez), tf, tf_entrez)
    
    # Signed sets: use TF_ACT / TF_REP naming
    if (length(act_targets) >= min_targets && length(act_targets) <= max_targets) {
      set_id <- paste0(tf, "_ACT")
      signed_sets[[set_id]] <- act_targets
      signed_desc[[set_id]] <- sprintf(
        "Genes activated by %s", tf
      )
    }
    if (length(rep_targets) >= min_targets && length(rep_targets) <= max_targets) {
      set_id <- paste0(tf, "_REP")
      signed_sets[[set_id]] <- rep_targets
      signed_desc[[set_id]] <- sprintf(
        "Genes repressed by %s", tf
      )
    }
    
    # Unsigned set: all targets
    if (length(all_targets) >= min_targets && length(all_targets) <= max_targets) {
      unsigned_sets[[tf]] <- all_targets
      unsigned_desc[[tf]] <- sprintf(
        "Genes targeted by %s", tf
      )
    }
  }
  
  # Write signed GMT
  signed_file <- file.path(output_dir, sprintf("collectri_signed_%s.gmt", organism))
  if (length(signed_sets) > 0) {
    write_gmt_from_list(signed_sets, signed_desc, signed_file)
  } else {
    writeLines(character(0), signed_file)
    message(sprintf("Wrote 0 gene sets to %s (no TFs passed min_targets filter)", signed_file))
  }
  
  # Write unsigned GMT
  unsigned_file <- file.path(output_dir, sprintf("collectri_unsigned_%s.gmt", organism))
  if (length(unsigned_sets) > 0) {
    write_gmt_from_list(unsigned_sets, unsigned_desc, unsigned_file)
  } else {
    writeLines(character(0), unsigned_file)
    message(sprintf("Wrote 0 gene sets to %s", unsigned_file))
  }
  
  # Summary
  n_act <- sum(grepl("_ACT$", names(signed_sets)))
  n_rep <- sum(grepl("_REP$", names(signed_sets)))
  
  message(sprintf(
    paste0(
      "\n=== CollecTRI Summary ===",
      "\nOrganism:        %s",
      "\nTarget filter:   %d-%d",
      "\nSigned:          %d sets (%d ACT, %d REP) in %s",
      "\nUnsigned:        %d sets in %s",
      "\nTotal TFs:       %d (from %d in full network)",
      "\nSplit complexes: %s"
    ),
    organism, min_targets, max_targets,
    length(signed_sets), n_act, n_rep, signed_file,
    length(unsigned_sets), unsigned_file,
    length(unsigned_sets), n_tf,
    ifelse(split_complexes, "yes", "no")
  ))
  
  invisible(list(
    signed_gmt = signed_file,
    unsigned_gmt = unsigned_file
  ))
}

# =============================================================================
# 6. DOWNLOAD ALL
# =============================================================================

#' Download all supported pathway databases as GMT files.
#'
#' @param org_db an OrgDb object (for PROGENy/CollecTRI symbol conversion)
#' @param species species name (default "Homo sapiens")
#' @param kegg_organism KEGG organism code (default "hsa")
#' @param progeny_organism PROGENy organism (default "human")
#' @param collectri_organism CollecTRI organism (default "human")
#' @param collectri_min_targets minimum targets per TF for CollecTRI (default 10)
#' @param collectri_max_targets maximum targets per TF for CollecTRI (default 500)
#' @param msigdb_collections which MSigDB collections to download
#'   (default: Hallmark + C2 curated pathways)
#' @param output_dir directory for all output files
#' @return invisible list of all file paths
download_all_pathway_gmt <- function(
    org_db,
    species = "Homo sapiens",
    kegg_organism = "hsa",
    progeny_organism = "human",
    collectri_organism = "human",
    collectri_min_targets = 10,
    collectri_max_targets = 500,
    msigdb_collections = c("H", "C2"),
    output_dir = "pathway_gmt"
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  files <- list()
  
  # MSigDB
  message("\n========== MSigDB ==========")
  files$msigdb <- download_msigdb_gmt(
    species = species,
    collections = msigdb_collections,
    output_dir = output_dir
  )
  
  # KEGG
  message("\n========== KEGG ==========")
  files$kegg <- download_kegg_gmt(
    organism = kegg_organism,
    output_dir = output_dir
  )
  
  # Reactome
  message("\n========== Reactome ==========")
  files$reactome <- download_reactome_gmt(
    species = species,
    method = "reactome.db",
    output_dir = output_dir
  )
  
  # PROGENy
  message("\n========== PROGENy ==========")
  files$progeny <- download_progeny_gmt(
    org_db = org_db,
    organism = progeny_organism,
    output_dir = output_dir
  )
  
  # CollecTRI
  message("\n========== CollecTRI ==========")
  files$collectri <- download_collectri_gmt(
    org_db = org_db,
    organism = collectri_organism,
    min_targets = collectri_min_targets,
    max_targets = collectri_max_targets,
    output_dir = output_dir
  )
  
  message("\n========== Done ==========")
  message(sprintf("All GMT files saved to: %s", output_dir))
  
  invisible(files)
}

# =============================================================================
# RUN EXAMPLE
# =============================================================================

if (sys.nframe() == 0) {
  library(org.Hs.eg.db)
  
  # Download everything
  all_files <- download_all_pathway_gmt(
    org_db = org.Hs.eg.db,
    species = "Homo sapiens",
    output_dir = "pathway_gmt"
  )
  
  # Or download individually:
  #
  # # Just Hallmark
  # download_msigdb_gmt(
  #   collections = "H",
  #   output_dir = "pathway_gmt"
  # )
  #
  # # Just KEGG Medicus from MSigDB (curated human disease pathways)
  # download_msigdb_gmt(
  #   collections = "C2",
  #   subcollections = list(C2 = "CP:KEGG_MEDICUS"),
  #   output_dir = "pathway_gmt"
  # )
  #
  # # Reactome via msigdbr (alternative to reactome.db)
  # download_reactome_gmt(method = "msigdbr", output_dir = "pathway_gmt")
  #
  # # PROGENy with more genes per pathway
  # download_progeny_gmt(
  #   org_db = org.Hs.eg.db,
  #   top_n = 500,
  #   output_dir = "pathway_gmt"
  # )
  #
  # # CollecTRI TF regulons with stricter filter
  # download_collectri_gmt(
  #   org_db = org.Hs.eg.db,
  #   organism = "human",
  #   min_targets = 20,
  #   output_dir = "pathway_gmt"
  # )
  
  # --- Usage with clusterProfiler ---
  # library(clusterProfiler)
  # source("build_go_regulation_gmt.R")  # for read_gmt_term2name()
  #
  # hallmark_gmt   <- read.gmt("pathway_gmt/msigdb_h_homo_sapiens.gmt")
  # hallmark_names <- read_gmt_term2name("pathway_gmt/msigdb_h_homo_sapiens.gmt")
  #
  # result <- enricher(
  #   gene      = my_genes,
  #   universe  = background,
  #   TERM2GENE = hallmark_gmt,
  #   TERM2NAME = hallmark_names
  # )
}
