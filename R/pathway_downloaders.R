resolve_msigdb_db_species <- function(species, db_species = NULL) {
  if (!is.null(db_species)) {
    resolved <- toupper(trimws(as.character(db_species)))
    if (!resolved %in% c("HS", "MM")) {
      warning(
        sprintf(
          "msigdb_db_species = '%s' is not one of the validated values ('HS', 'MM'). Proceeding as requested.",
          resolved
        ),
        call. = FALSE
      )
    }
    return(resolved)
  }

  if (is_mouse_species(species)) {
    message("Using mouse-native MSigDB with db_species = 'MM'.")
    return("MM")
  }

  if (is_rat_species(species)) {
    warning(
      "Rat builds use human MSigDB (db_species = 'HS') with ortholog mapping because msigdbr does not provide a rat-native MSigDB database.",
      call. = FALSE
    )
    return("HS")
  }

  "HS"
}

translate_msigdb_collection_code <- function(collection, db_species = "HS") {
  collection <- toupper(trimws(as.character(collection)))

  if (db_species != "MM") {
    return(collection)
  }

  mouse_collection_map <- c(
    H = "MH",
    C1 = "M1",
    C2 = "M2",
    C3 = "M3",
    C5 = "M5",
    C7 = "M7",
    C8 = "M8"
  )

  mouse_collection_map[[collection]] %||% collection
}

download_msigdb_gmt <- function(
  species = "Homo sapiens",
  db_species = NULL,
  collections = "all",
  subcollections = NULL,
  output_dir = "."
) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Install msigdbr: install.packages('msigdbr')", call. = FALSE)
  }

  species <- canonical_species_name(species)
  db_species <- resolve_msigdb_db_species(species, db_species = db_species)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  species_tag <- gsub(" ", "_", tolower(species))
  available <- msigdbr::msigdbr_collections(db_species = db_species)

  if (identical(collections, "all")) {
    collections <- unique(available$gs_collection)
  }

  files <- list()
  for (requested_collection in collections) {
    collection <- translate_msigdb_collection_code(requested_collection, db_species = db_species)
    if (!(collection %in% available$gs_collection)) {
      stop(
        sprintf(
          "Collection '%s' is not available for msigdb_db_species = '%s'. For mouse builds, set msigdb_db_species = 'HS' if you want ortholog-mapped human collections instead.",
          requested_collection,
          db_species
        ),
        call. = FALSE
      )
    }

    if (db_species == "MM" && collection != toupper(trimws(as.character(requested_collection)))) {
      message(sprintf(
        "Using mouse-native MSigDB collection '%s' for requested collection '%s'.",
        collection,
        requested_collection
      ))
    }

    subcollection <- NULL
    if (!is.null(subcollections)) {
      subcollection <- subcollections[[requested_collection]]
      if (is.null(subcollection) && requested_collection != collection) {
        subcollection <- subcollections[[collection]]
      }
    }

    if (is.null(subcollection)) {
      msigdb_df <- msigdbr::msigdbr(
        db_species = db_species,
        species = species,
        collection = collection
      )
      file_tag <- tolower(collection)
    } else {
      msigdb_df <- msigdbr::msigdbr(
        db_species = db_species,
        species = species,
        collection = collection,
        subcollection = subcollection
      )
      file_tag <- paste0(tolower(collection), "_", tolower(gsub("[: ]", "_", subcollection)))
    }

    if (!nrow(msigdb_df)) {
      message(sprintf("  Skipping %s: no gene sets found", collection))
      next
    }

    msigdb_df <- msigdb_df[!is.na(msigdb_df$ncbi_gene) & msigdb_df$ncbi_gene != "", , drop = FALSE]
    gene_sets <- lapply(split(msigdb_df$ncbi_gene, msigdb_df$gs_name), unique)
    description_df <- unique(msigdb_df[, c("gs_name", "gs_description"), drop = FALSE])
    descriptions <- stats::setNames(description_df$gs_description, description_df$gs_name)

    filepath <- file.path(output_dir, sprintf("msigdb_%s_%s.gmt", file_tag, species_tag))
    write_gmt_from_list(gene_sets, descriptions, filepath)
    files[[requested_collection]] <- filepath
  }

  invisible(files)
}

download_kegg_gmt <- function(
  organism = "hsa",
  output_dir = "."
) {
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    stop("Install KEGGREST: BiocManager::install('KEGGREST')", call. = FALSE)
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  message("Fetching KEGG pathway list...")
  pathways <- KEGGREST::keggList("pathway", organism)
  raw_ids <- sub("^path:", "", names(pathways))
  pathway_ids <- ifelse(grepl(paste0("^", organism), raw_ids), raw_ids, paste0(organism, raw_ids))
  pathway_names <- sub(" - .*$", "", unname(pathways))

  message(sprintf("  Found %d pathways", length(pathway_ids)))
  message("Fetching gene members (this may take a few minutes)...")

  gene_sets <- list()
  descriptions <- character()

  for (idx in seq_along(pathway_ids)) {
    pathway_id <- pathway_ids[[idx]]
    tryCatch({
      pathway_record <- KEGGREST::keggGet(paste0("path:", pathway_id))[[1]]
      genes <- pathway_record$GENE

      if (!is.null(genes)) {
        entrez_ids <- genes[seq(1, length(genes), by = 2)]
        gene_sets[[pathway_id]] <- unique(as.character(entrez_ids))
        descriptions[[pathway_id]] <- pathway_names[[idx]]
      }
    }, error = function(e) {
      message(sprintf("    Warning: failed to fetch %s", pathway_id))
    })

    if (idx %% 10 == 0) {
      message(sprintf("    Processed %d / %d pathways", idx, length(pathway_ids)))
      Sys.sleep(1)
    }
  }

  gene_sets <- gene_sets[lengths(gene_sets) > 0]
  message(sprintf("  %d pathways with gene annotations", length(gene_sets)))

  filepath <- file.path(output_dir, sprintf("kegg_%s.gmt", organism))
  write_gmt_from_list(gene_sets, descriptions, filepath)
  invisible(filepath)
}

download_reactome_gmt <- function(
  species = "Homo sapiens",
  method = "reactome.db",
  output_dir = "."
) {
  species <- canonical_species_name(species)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (method == "msigdbr") {
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
      stop("Install msigdbr: install.packages('msigdbr')", call. = FALSE)
    }

    db_species <- resolve_msigdb_db_species(species)
    reactome_collection <- translate_msigdb_collection_code("C2", db_species = db_species)
    message("Fetching Reactome from MSigDB C2:CP:REACTOME...")
    reactome_df <- msigdbr::msigdbr(
      db_species = db_species,
      species = species,
      collection = reactome_collection,
      subcollection = "CP:REACTOME"
    )
    reactome_df <- reactome_df[!is.na(reactome_df$ncbi_gene) & reactome_df$ncbi_gene != "", , drop = FALSE]

    gene_sets <- lapply(split(reactome_df$ncbi_gene, reactome_df$gs_name), unique)
    description_df <- unique(reactome_df[, c("gs_name", "gs_description"), drop = FALSE])
    descriptions <- stats::setNames(description_df$gs_description, description_df$gs_name)
    species_tag <- gsub(" ", "_", tolower(species))
    filepath <- file.path(output_dir, sprintf("reactome_msigdb_%s.gmt", species_tag))
  } else if (method == "reactome.db") {
    if (!requireNamespace("reactome.db", quietly = TRUE)) {
      stop("Install reactome.db: BiocManager::install('reactome.db')", call. = FALSE)
    }

    message("Fetching Reactome pathways from reactome.db...")

    pathway_genes <- AnnotationDbi::select(
      reactome.db::reactome.db,
      keys = AnnotationDbi::keys(reactome.db::reactome.db, keytype = "PATHID"),
      keytype = "PATHID",
      columns = c("PATHID", "ENTREZID")
    )
    pathway_genes <- pathway_genes[!is.na(pathway_genes$ENTREZID), , drop = FALSE]

    pathway_names_df <- AnnotationDbi::select(
      reactome.db::reactome.db,
      keys = unique(pathway_genes$PATHID),
      keytype = "PATHID",
      columns = c("PATHID", "PATHNAME")
    )

    species_prefix <- paste0(species, ": ")
    prefixed <- startsWith(pathway_names_df$PATHNAME, species_prefix)
    if (any(prefixed)) {
      pathway_names_df <- pathway_names_df[prefixed, , drop = FALSE]
      pathway_names_df$PATHNAME <- sub(species_prefix, "", pathway_names_df$PATHNAME, fixed = TRUE)
    }

    keep_ids <- unique(pathway_names_df$PATHID)
    pathway_genes <- pathway_genes[pathway_genes$PATHID %in% keep_ids, , drop = FALSE]

    gene_sets <- lapply(split(pathway_genes$ENTREZID, pathway_genes$PATHID), unique)
    descriptions <- stats::setNames(pathway_names_df$PATHNAME, pathway_names_df$PATHID)
    species_tag <- gsub(" ", "_", tolower(species))
    filepath <- file.path(output_dir, sprintf("reactome_%s.gmt", species_tag))
  } else {
    stop("method must be one of: 'reactome.db', 'msigdbr'", call. = FALSE)
  }

  message(sprintf("  %d Reactome pathways", length(gene_sets)))
  write_gmt_from_list(gene_sets, descriptions, filepath)
  invisible(filepath)
}

download_progeny_gmt <- function(
    org_db,
    organism = "human",
    top_n = 100,
    output_dir = "."
) {
  if (!requireNamespace("progeny", quietly = TRUE)) {
    stop("Install progeny: BiocManager::install('progeny')", call. = FALSE)
  }

  assert_positive_count(top_n, "top_n")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  model <- switch(
    organism,
    "human" = progeny::model_human_full,
    "mouse" = progeny::model_mouse_full,
    stop("organism must be 'human' or 'mouse'", call. = FALSE)
  )

  message(sprintf(
    "PROGENy %s model: %d genes x %d pathways",
    organism,
    length(unique(model$gene)),
    length(unique(model$pathway))
  ))

  model_top <- do.call(rbind, lapply(split(model, model$pathway), function(pathway_df) {
    rownames(pathway_df) <- NULL
    pathway_df <- pathway_df[order(pathway_df$p.value), , drop = FALSE]
    head(pathway_df, top_n)
  }))
  rownames(model_top) <- NULL
  model_top$gene <- as.character(model_top$gene)
  model_top$pathway <- as.character(model_top$pathway)

  message(sprintf(
    "  Using top %d genes per pathway (%d total gene-pathway pairs)",
    top_n,
    nrow(model_top)
  ))

  symbols <- unique(model_top$gene)
  message(sprintf("  Converting %d unique gene symbols to Entrez IDs...", length(symbols)))
  message(sprintf("  Example symbols: %s", paste(head(symbols, 10), collapse = ", ")))

  symbol_to_entrez <- map_symbols_to_entrez(org_db, symbols)
  n_mapped <- sum(!is.na(symbol_to_entrez))
  message(sprintf("  Mapped %d / %d gene symbols to Entrez IDs", n_mapped, length(symbols)))

  signed_sets <- list()
  unsigned_sets <- list()

  for (pathway in unique(model_top$pathway)) {
    pathway_df <- model_top[model_top$pathway == pathway, , drop = FALSE]
    pos_entrez <- unique(na.omit(symbol_to_entrez[as.character(pathway_df$gene[pathway_df$weight > 0])]))
    neg_entrez <- unique(na.omit(symbol_to_entrez[as.character(pathway_df$gene[pathway_df$weight < 0])]))
    all_entrez <- unique(na.omit(symbol_to_entrez[as.character(pathway_df$gene)]))

    if (length(pos_entrez)) {
      signed_sets[[paste0(pathway, "_UP")]] <- unname(pos_entrez)
    }
    if (length(neg_entrez)) {
      signed_sets[[paste0(pathway, "_DN")]] <- unname(neg_entrez)
    }
    if (length(all_entrez)) {
      unsigned_sets[[pathway]] <- unname(all_entrez)
    }
  }

  signed_descriptions <- stats::setNames(
    sub("_UP$", " upregulated genes", sub("_DN$", " downregulated genes", names(signed_sets))),
    names(signed_sets)
  )
  signed_file <- file.path(output_dir, sprintf("progeny_signed_%s_top%d.gmt", organism, top_n))
  write_gmt_from_list(signed_sets, signed_descriptions, signed_file)

  unsigned_descriptions <- stats::setNames(
    paste(names(unsigned_sets), "regulated genes"),
    names(unsigned_sets)
  )
  unsigned_file <- file.path(output_dir, sprintf("progeny_unsigned_%s_top%d.gmt", organism, top_n))
  write_gmt_from_list(unsigned_sets, unsigned_descriptions, unsigned_file)

  message(sprintf(
    paste0(
      "\n=== PROGENy Summary ===",
      "\nOrganism:  %s",
      "\nTop genes: %d per pathway",
      "\nSigned:    %d sets in %s",
      "\nUnsigned:  %d sets in %s"
    ),
    organism,
    top_n,
    length(signed_sets),
    signed_file,
    length(unsigned_sets),
    unsigned_file
  ))

  invisible(list(
    signed_gmt = signed_file,
    unsigned_gmt = unsigned_file
  ))
}

collectri_taxonomy_id <- function(organism) {
  canonical <- canonical_species_name(organism)

  switch(
    canonical,
    "Homo sapiens" = 9606L,
    "Mus musculus" = 10090L,
    "Rattus norvegicus" = 10116L,
    stop(sprintf("Unsupported CollecTRI organism: %s", organism), call. = FALSE)
  )
}

normalize_collectri_direct_query <- function(collectri, split_complexes = FALSE) {
  required_cols <- c(
    "source",
    "source_genesymbol",
    "target_genesymbol",
    "is_stimulation",
    "is_inhibition"
  )

  missing_cols <- setdiff(required_cols, colnames(collectri))
  if (length(missing_cols)) {
    stop(
      sprintf(
        "Direct CollecTRI query is missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  cols <- c("source", "source_genesymbol", "target_genesymbol", "is_stimulation", "is_inhibition")
  is_complex <- grepl("COMPLEX", collectri$source, fixed = TRUE)
  collectri_interactions <- collectri[!is_complex, cols, drop = FALSE]
  collectri_complex <- collectri[is_complex, cols, drop = FALSE]

  if (!split_complexes && nrow(collectri_complex)) {
    complex_symbols <- as.character(collectri_complex$source_genesymbol)
    complex_symbols[grepl("JUN", complex_symbols) | grepl("FOS", complex_symbols)] <- "AP1"
    complex_symbols[grepl("REL", complex_symbols) | grepl("NFKB", complex_symbols)] <- "NFKB"
    collectri_complex$source_genesymbol <- complex_symbols
  }

  collapsed <- unique(rbind(collectri_interactions, collectri_complex))
  collapsed <- collapsed[
    !is.na(collapsed$source_genesymbol) &
      !is.na(collapsed$target_genesymbol),
    ,
    drop = FALSE
  ]
  collapsed$mor <- ifelse(collapsed$is_stimulation == 1, 1, -1)

  data.frame(
    source = tools::toTitleCase(tolower(as.character(collapsed$source_genesymbol))),
    target = tools::toTitleCase(tolower(as.character(collapsed$target_genesymbol))),
    mor = as.numeric(collapsed$mor),
    stringsAsFactors = FALSE
  )
}

fetch_collectri_network <- function(organism = "human", split_complexes = FALSE) {
  primary_result <- tryCatch(
    decoupleR::get_collectri(
      organism = organism,
      split_complexes = split_complexes
    ),
    error = identity
  )

  if (!inherits(primary_result, "error")) {
    return(primary_result)
  }

  if (!is_rat_species(organism)) {
    stop(primary_result$message, call. = FALSE)
  }

  warning(
    paste(
      "decoupleR::get_collectri failed for rat; trying an OmnipathR fallback based on a direct CollecTRI interactions query.",
      "This path still depends on OmnipathR and Ensembl species resolution."
    ),
    call. = FALSE
  )

  taxonomy_id <- collectri_taxonomy_id(organism)
  fallback_attempts <- list(
    list(
      label = "OmnipathR::collectri(query_type = 'interactions', organism = 'rat')",
      run = function() OmnipathR::collectri(
        query_type = "interactions",
        organism = "rat",
        genesymbol = TRUE,
        loops = TRUE
      )
    ),
    list(
      label = "OmnipathR::collectri(query_type = 'interactions', organism = 10116L)",
      run = function() OmnipathR::collectri(
        query_type = "interactions",
        organism = taxonomy_id,
        genesymbol = TRUE,
        loops = TRUE
      )
    ),
    list(
      label = "OmnipathR::static_table(query = 'interactions', resource = 'collectri', organism = 10116L)",
      run = function() OmnipathR::static_table(
        query = "interactions",
        resource = "collectri",
        organism = taxonomy_id
      )
    )
  )

  fallback_errors <- character(0)
  for (attempt in fallback_attempts) {
    raw_result <- tryCatch(attempt$run(), error = identity)
    if (inherits(raw_result, "error")) {
      fallback_errors <- c(fallback_errors, sprintf("%s: %s", attempt$label, raw_result$message))
      next
    }

    normalized <- tryCatch(
      normalize_collectri_direct_query(raw_result, split_complexes = split_complexes),
      error = identity
    )
    if (!inherits(normalized, "error")) {
      return(normalized)
    }

    fallback_errors <- c(fallback_errors, sprintf("%s normalization: %s", attempt$label, normalized$message))
  }

  stop(
    paste(
      "Failed to retrieve rat CollecTRI through decoupleR and all OmnipathR fallbacks.",
      sprintf("Primary error: %s", primary_result$message),
      sprintf("Fallback errors: %s", paste(fallback_errors, collapse = " | "))
    ),
    call. = FALSE
  )
}

download_collectri_gmt <- function(
    org_db,
    organism = "human",
    split_complexes = FALSE,
    min_targets = 10,
    max_targets = 500,
    output_dir = "."
) {
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop("Install decoupleR: BiocManager::install('decoupleR')", call. = FALSE)
  }

  assert_positive_count(min_targets, "min_targets", allow_zero = TRUE)
  assert_positive_count(max_targets, "max_targets")
  assert_min_max(min_targets, max_targets, "min_targets", "max_targets")

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  message(sprintf("Fetching CollecTRI network (organism: %s)...", organism))
  network <- tryCatch(
    fetch_collectri_network(
      organism = organism,
      split_complexes = split_complexes
    ),
    error = function(e) {
      stop(
        "Failed to retrieve CollecTRI. This step depends on OmniPath/OmnipathR and may require network access or cached static tables.\nError: ",
        e$message,
        call. = FALSE
      )
    }
  )

  n_tf <- length(unique(network$source))
  n_interactions <- nrow(network)
  message(sprintf("  %d TFs, %d interactions", n_tf, n_interactions))

  target_symbols <- unique(network$target)
  message(sprintf("  Converting %d unique target symbols to Entrez IDs...", length(target_symbols)))
  target_to_entrez <- map_symbols_to_entrez(org_db, target_symbols)
  n_mapped <- sum(!is.na(target_to_entrez))
  message(sprintf("  Mapped %d / %d target symbols to Entrez IDs", n_mapped, length(target_symbols)))

  if (!n_mapped) {
    warning("No gene symbols could be mapped. Check org_db matches the organism.")
    return(invisible(list(signed_gmt = NULL, unsigned_gmt = NULL)))
  }

  network$entrez <- target_to_entrez[network$target]
  network <- network[!is.na(network$entrez), , drop = FALSE]

  signed_sets <- list()
  signed_descriptions <- character()
  unsigned_sets <- list()
  unsigned_descriptions <- character()

  for (tf in unique(network$source)) {
    tf_network <- network[network$source == tf, , drop = FALSE]

    activated_targets <- unique(tf_network$entrez[tf_network$mor > 0])
    repressed_targets <- unique(tf_network$entrez[tf_network$mor < 0])
    all_targets <- unique(tf_network$entrez)

    if (length(activated_targets) > 0 &&
        length(activated_targets) >= min_targets &&
        length(activated_targets) <= max_targets) {
      set_id <- paste0(tf, "_ACT")
      signed_sets[[set_id]] <- activated_targets
      signed_descriptions[[set_id]] <- sprintf("Genes activated by %s", tf)
    }

    if (length(repressed_targets) > 0 &&
        length(repressed_targets) >= min_targets &&
        length(repressed_targets) <= max_targets) {
      set_id <- paste0(tf, "_REP")
      signed_sets[[set_id]] <- repressed_targets
      signed_descriptions[[set_id]] <- sprintf("Genes repressed by %s", tf)
    }

    if (length(all_targets) > 0 &&
        length(all_targets) >= min_targets &&
        length(all_targets) <= max_targets) {
      unsigned_sets[[tf]] <- all_targets
      unsigned_descriptions[[tf]] <- sprintf("Genes targeted by %s", tf)
    }
  }

  signed_file <- file.path(output_dir, sprintf("collectri_signed_%s.gmt", organism))
  if (length(signed_sets)) {
    write_gmt_from_list(signed_sets, signed_descriptions, signed_file)
  } else {
    writeLines(character(0), signed_file)
    message(sprintf("Wrote 0 gene sets to %s (no TFs passed the target-size filter)", signed_file))
  }

  unsigned_file <- file.path(output_dir, sprintf("collectri_unsigned_%s.gmt", organism))
  if (length(unsigned_sets)) {
    write_gmt_from_list(unsigned_sets, unsigned_descriptions, unsigned_file)
  } else {
    writeLines(character(0), unsigned_file)
    message(sprintf("Wrote 0 gene sets to %s", unsigned_file))
  }

  n_activated <- sum(grepl("_ACT$", names(signed_sets)))
  n_repressed <- sum(grepl("_REP$", names(signed_sets)))

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
    organism,
    min_targets,
    max_targets,
    length(signed_sets),
    n_activated,
    n_repressed,
    signed_file,
    length(unsigned_sets),
    unsigned_file,
    length(unsigned_sets),
    n_tf,
    if (split_complexes) "yes" else "no"
  ))

  invisible(list(
    signed_gmt = signed_file,
    unsigned_gmt = unsigned_file
  ))
}

download_all_pathway_gmt <- function(
    org_db,
    species = "Homo sapiens",
    msigdb_db_species = NULL,
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

  message("\n========== MSigDB ==========")
  files$msigdb <- download_msigdb_gmt(
    species = species,
    db_species = msigdb_db_species,
    collections = msigdb_collections,
    output_dir = output_dir
  )

  message("\n========== KEGG ==========")
  files$kegg <- download_kegg_gmt(
    organism = kegg_organism,
    output_dir = output_dir
  )

  message("\n========== Reactome ==========")
  files$reactome <- download_reactome_gmt(
    species = species,
    method = "reactome.db",
    output_dir = output_dir
  )

  message("\n========== PROGENy ==========")
  files$progeny <- download_progeny_gmt(
    org_db = org_db,
    organism = progeny_organism,
    output_dir = output_dir
  )

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
