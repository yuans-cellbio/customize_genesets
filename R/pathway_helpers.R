suppressPackageStartupMessages({
  library(AnnotationDbi)
})

`%||%` <- function(x, y) {
  if (is.null(x)) {
    y
  } else {
    x
  }
}

assert_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 || x > 1) {
    stop(sprintf("%s must be a single numeric value between 0 and 1", name), call. = FALSE)
  }
}

assert_positive_count <- function(x, name, allow_zero = FALSE) {
  lower_bound <- if (allow_zero) 0 else 1
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < lower_bound || x != as.integer(x)) {
    qualifier <- if (allow_zero) "a non-negative integer" else "a positive integer"
    stop(sprintf("%s must be %s", name, qualifier), call. = FALSE)
  }
}

assert_min_max <- function(min_value, max_value, min_name, max_name) {
  if (min_value > max_value) {
    stop(sprintf("%s must be less than or equal to %s", min_name, max_name), call. = FALSE)
  }
}

organism_name_from_orgdb <- function(org_db) {
  meta <- tryCatch(AnnotationDbi::metadata(org_db), error = function(e) NULL)
  if (is.null(meta)) {
    return("unknown")
  }

  organism <- meta$value[meta$name == "ORGANISM"]
  if (!length(organism)) {
    return("unknown")
  }

  organism[[1]]
}

organism_slug_from_orgdb <- function(org_db) {
  gsub(" ", "_", tolower(organism_name_from_orgdb(org_db)))
}

normalize_species_label <- function(species) {
  if (!is.character(species) || length(species) != 1L || is.na(species) || !nzchar(trimws(species))) {
    stop("species must be a single non-empty character string", call. = FALSE)
  }

  tolower(trimws(species))
}

canonical_species_name <- function(species) {
  normalized <- normalize_species_label(species)

  switch(
    normalized,
    "human" = "Homo sapiens",
    "homo sapiens" = "Homo sapiens",
    "mouse" = "Mus musculus",
    "mus musculus" = "Mus musculus",
    "rat" = "Rattus norvegicus",
    "rattus norvegicus" = "Rattus norvegicus",
    species
  )
}

is_human_species <- function(species) {
  canonical_species_name(species) == "Homo sapiens"
}

is_mouse_species <- function(species) {
  canonical_species_name(species) == "Mus musculus"
}

is_rat_species <- function(species) {
  canonical_species_name(species) == "Rattus norvegicus"
}

safe_package_version <- function(package_name) {
  tryCatch(
    as.character(packageVersion(package_name)),
    error = function(e) "not installed"
  )
}

compact_gene_vector <- function(genes) {
  genes <- as.character(genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]
  unique(genes)
}

compact_gene_sets <- function(gene_sets) {
  stats::setNames(
    lapply(gene_sets, compact_gene_vector),
    names(gene_sets)
  )
}

write_gmt_from_list <- function(gene_sets, descriptions = NULL, filepath) {
  gene_sets <- compact_gene_sets(gene_sets)
  gene_sets <- gene_sets[lengths(gene_sets) > 0]

  lines <- vapply(names(gene_sets), function(name) {
    desc <- if (!is.null(descriptions) && name %in% names(descriptions)) {
      descriptions[[name]]
    } else {
      "NA"
    }

    paste(c(name, desc, gene_sets[[name]]), collapse = "\t")
  }, character(1))

  writeLines(lines, filepath)
  message(sprintf("Wrote %d gene sets to %s", length(gene_sets), filepath))
  invisible(filepath)
}

write_gmt_from_terms <- function(gene_sets, term_df, filepath) {
  gene_sets <- compact_gene_sets(gene_sets)
  gene_sets <- gene_sets[lengths(gene_sets) > 0]
  term_lookup <- stats::setNames(term_df$TERM, term_df$GOID)

  lines <- vapply(names(gene_sets), function(go_id) {
    term_name <- term_lookup[[go_id]] %||% go_id
    paste(c(go_id, term_name, gene_sets[[go_id]]), collapse = "\t")
  }, character(1))

  writeLines(lines, filepath)
  message(sprintf("Wrote %d gene sets to %s", length(gene_sets), filepath))
  invisible(filepath)
}

read_gmt_term2name <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(sprintf("GMT file not found: %s", filepath), call. = FALSE)
  }

  lines <- readLines(filepath, warn = FALSE)
  if (!length(lines)) {
    return(data.frame(ID = character(0), Name = character(0), stringsAsFactors = FALSE))
  }

  parts <- strsplit(lines, "\t", fixed = TRUE)
  data.frame(
    ID = vapply(parts, `[`, character(1), 1),
    Name = vapply(parts, `[`, character(1), 2),
    stringsAsFactors = FALSE
  )
}

map_symbols_to_entrez <- function(
    org_db,
    symbols,
    keytypes = c("SYMBOL", "ALIAS"),
    multiVals = "first"
) {
  symbols <- compact_gene_vector(symbols)
  if (!length(symbols)) {
    return(stats::setNames(character(0), character(0)))
  }

  mappings <- stats::setNames(rep(NA_character_, length(symbols)), symbols)

  for (keytype in keytypes) {
    attempted <- tryCatch(
      AnnotationDbi::mapIds(
        org_db,
        keys = symbols,
        keytype = keytype,
        column = "ENTREZID",
        multiVals = multiVals
      ),
      error = function(e) NULL
    )

    if (is.null(attempted)) {
      next
    }

    fillable <- is.na(mappings) & !is.na(attempted[names(mappings)])
    mappings[fillable] <- as.character(attempted[names(mappings)][fillable])
  }

  mappings
}

source_project_modules <- function(project_root = ".") {
  module_paths <- file.path(
    project_root,
    "R",
    c(
      "go_library.R",
      "pathway_downloaders.R",
      "pathway_library_builder.R",
      "pathway_analysis.R"
    )
  )

  invisible(lapply(module_paths, source, local = FALSE))
}
