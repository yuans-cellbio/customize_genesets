DEFAULT_POS_KEYWORDS <- c(
  "positive regulation", "activation of", "upregulation of",
  "induction of", "stimulation of", "promotion of",
  "enhancement of", "potentiation of"
)

DEFAULT_NEG_KEYWORDS <- c(
  "negative regulation", "inhibition of", "downregulation of",
  "suppression of", "repression of", "attenuation of",
  "blockade of", "sequestering of", "desensitization of",
  "deactivation of", "inactivation of", "degradation of"
)

DEFAULT_DIRECTIONAL_PAIR_REPLACEMENTS <- list(
  "positive regulation" = c("negative regulation"),
  "activation of" = c("inhibition of", "deactivation of", "inactivation of"),
  "upregulation of" = c("downregulation of"),
  "up-regulation of" = c("down-regulation of"),
  "up regulation of" = c("down regulation of"),
  "stimulation of" = c("suppression of")
)

discover_directional_keywords <- function(
  pos_keywords = DEFAULT_POS_KEYWORDS,
  neg_keywords = DEFAULT_NEG_KEYWORDS,
  candidate_words = NULL
) {
  bp_terms <- AnnotationDbi::select(
    GO.db::GO.db,
    keys = "BP",
    keytype = "ONTOLOGY",
    columns = c("GOID", "TERM")
  )

  terms_lower <- tolower(bp_terms$TERM)

  if (is.null(candidate_words)) {
    candidate_words <- c(
      "positive regulation", "activation of", "upregulation of",
      "induction of", "stimulation of", "promotion of",
      "enhancement of", "augmentation of", "potentiation of",
      "facilitation of", "up-regulation of", "up regulation of",
      "negative regulation", "inhibition of", "downregulation of",
      "suppression of", "repression of", "attenuation of",
      "blockade of", "inactivation of", "desensitization of",
      "deactivation of", "silencing of", "down-regulation of",
      "down regulation of", "sequestering of", "destabilization of",
      "degradation of"
    )
  }

  known <- c(pos_keywords, neg_keywords)
  results <- data.frame(
    keyword = character(0),
    n_terms = integer(0),
    already_known = logical(0),
    example_terms = character(0),
    stringsAsFactors = FALSE
  )

  for (keyword in candidate_words) {
    matches <- grep(tolower(keyword), terms_lower, fixed = TRUE)
    if (!length(matches)) {
      next
    }

    examples <- utils::head(bp_terms$TERM[matches], 3)
    results <- rbind(
      results,
      data.frame(
        keyword = keyword,
        n_terms = length(matches),
        already_known = tolower(keyword) %in% tolower(known),
        example_terms = paste(examples, collapse = " | "),
        stringsAsFactors = FALSE
      )
    )
  }

  results <- results[order(-results$n_terms), , drop = FALSE]
  rownames(results) <- NULL

  message(sprintf(
    "Found %d candidate keywords with GO BP matches (%d already in your lists)",
    nrow(results),
    sum(results$already_known)
  ))

  unseen <- results[!results$already_known, , drop = FALSE]
  if (nrow(unseen)) {
    message("\nKeywords not yet in your lists:")
    for (idx in seq_len(nrow(unseen))) {
      message(sprintf(
        '  "%s" (%d terms) e.g.: %s',
        unseen$keyword[idx],
        unseen$n_terms[idx],
        unseen$example_terms[idx]
      ))
    }
  }

  results
}

earliest_match_pos <- function(term, keywords) {
  positions <- vapply(keywords, function(keyword) {
    match_pos <- regexpr(keyword, term, fixed = TRUE)
    if (match_pos == -1L) {
      Inf
    } else {
      as.double(match_pos)
    }
  }, double(1))

  min(positions)
}

normalize_directional_pair_replacements <- function(pair_replacements) {
  normalized <- lapply(pair_replacements, function(values) unique(tolower(as.character(values))))
  names(normalized) <- tolower(names(pair_replacements))

  reversed <- list()
  for (from in names(normalized)) {
    for (to in normalized[[from]]) {
      reversed[[to]] <- unique(c(reversed[[to]], from))
    }
  }

  merged <- utils::modifyList(normalized, reversed, keep.null = TRUE)
  merged[order(names(merged))]
}

build_directional_pair_map <- function(
    term_df,
    pair_replacements = DEFAULT_DIRECTIONAL_PAIR_REPLACEMENTS
) {
  replacement_map <- normalize_directional_pair_replacements(pair_replacements)
  term_lookup <- stats::setNames(tolower(term_df$TERM), term_df$GOID)
  name_to_id <- stats::setNames(term_df$GOID, tolower(term_df$TERM))
  pair_map <- as.list(stats::setNames(rep(NA_character_, nrow(term_df)), term_df$GOID))

  for (go_id in names(pair_map)) {
    term <- term_lookup[[go_id]]
    found <- NA_character_

    for (from in names(replacement_map)) {
      if (!grepl(from, term, fixed = TRUE)) {
        next
      }

      for (candidate_keyword in replacement_map[[from]]) {
        candidate_name <- sub(from, candidate_keyword, term, fixed = TRUE)
        candidate_id <- name_to_id[candidate_name]
        if (length(candidate_id) && !is.na(candidate_id[[1]])) {
          found <- unname(candidate_id)
          break
        }
      }

      if (!is.na(found)) {
        break
      }
    }

    pair_map[[go_id]] <- found
  }

  pair_map
}

classify_regulation_terms <- function(
  pos_keywords = DEFAULT_POS_KEYWORDS,
  neg_keywords = DEFAULT_NEG_KEYWORDS,
  pair_replacements = DEFAULT_DIRECTIONAL_PAIR_REPLACEMENTS
) {
  bp_terms <- AnnotationDbi::select(
    GO.db::GO.db,
    keys = "BP",
    keytype = "ONTOLOGY",
    columns = c("GOID", "TERM")
  )

  terms_lower <- tolower(bp_terms$TERM)
  pos_keywords_lower <- tolower(pos_keywords)
  neg_keywords_lower <- tolower(neg_keywords)
  direction <- character(nrow(bp_terms))

  for (idx in seq_along(terms_lower)) {
    pos_pos <- earliest_match_pos(terms_lower[[idx]], pos_keywords_lower)
    neg_pos <- earliest_match_pos(terms_lower[[idx]], neg_keywords_lower)

    if (is.infinite(pos_pos) && is.infinite(neg_pos)) {
      direction[[idx]] <- "unsigned"
    } else if (pos_pos < neg_pos) {
      direction[[idx]] <- "positive"
    } else if (neg_pos < pos_pos) {
      direction[[idx]] <- "negative"
    } else {
      direction[[idx]] <- "unsigned"
      message(sprintf("  Tie at position %s for term: %s", pos_pos, bp_terms$TERM[[idx]]))
    }
  }

  list(
    positive = bp_terms$GOID[direction == "positive"],
    negative = bp_terms$GOID[direction == "negative"],
    unsigned = bp_terms$GOID[direction == "unsigned"],
    pair_map = build_directional_pair_map(bp_terms, pair_replacements = pair_replacements),
    term_df = bp_terms
  )
}

get_go_genes <- function(go_ids, org_db, include_IEA = TRUE) {
  mapping <- AnnotationDbi::select(
    org_db,
    keys = go_ids,
    keytype = "GOALL",
    columns = c("ENTREZID", "EVIDENCEALL")
  )

  if (!include_IEA) {
    mapping <- mapping[mapping$EVIDENCEALL != "IEA", , drop = FALSE]
  }

  mapping <- unique(mapping[, c("GOALL", "ENTREZID"), drop = FALSE])
  mapping <- mapping[!is.na(mapping$ENTREZID), , drop = FALSE]
  split(as.character(mapping$ENTREZID), mapping$GOALL)
}

filter_by_size <- function(gene_sets, min_size = 15, max_size = 500) {
  assert_positive_count(min_size, "min_size")
  assert_positive_count(max_size, "max_size")
  assert_min_max(min_size, max_size, "min_size", "max_size")

  gene_sets <- compact_gene_sets(gene_sets)
  sizes <- vapply(gene_sets, length, integer(1))
  gene_sets[sizes >= min_size & sizes <= max_size]
}

go_ancestor_lookup <- function(ontology = "BP") {
  switch(
    toupper(ontology),
    "BP" = as.list(GO.db::GOBPANCESTOR),
    "CC" = as.list(GO.db::GOCCANCESTOR),
    "MF" = as.list(GO.db::GOMFANCESTOR),
    stop("ontology must be one of: BP, CC, MF", call. = FALSE)
  )
}

get_go_term_specificity <- function(go_ids, ontology = "BP") {
  ancestor_lookup <- go_ancestor_lookup(ontology)
  stats::setNames(vapply(go_ids, function(go_id) {
    ancestors <- ancestor_lookup[[go_id]]
    length(unique(ancestors[!is.na(ancestors)]))
  }, integer(1)), go_ids)
}

prune_by_dag <- function(gene_sets, ontology = "BP", overlap_threshold = 0.8) {
  assert_probability(overlap_threshold, "overlap_threshold")

  go_ids <- names(gene_sets)
  if (length(go_ids) < 2L) {
    return(gene_sets)
  }

  ancestor_lookup <- go_ancestor_lookup(ontology)
  to_drop <- character(0)

  for (child_id in go_ids) {
    ancestors <- ancestor_lookup[[child_id]]
    if (is.null(ancestors)) {
      next
    }

    for (ancestor_id in intersect(ancestors, go_ids)) {
      overlap <- length(intersect(gene_sets[[child_id]], gene_sets[[ancestor_id]]))
      frac <- overlap / length(gene_sets[[ancestor_id]])
      if (frac >= overlap_threshold) {
        to_drop <- c(to_drop, ancestor_id)
      }
    }
  }

  to_drop <- unique(to_drop)
  message(sprintf("  DAG pruning: dropping %d redundant ancestor terms", length(to_drop)))
  gene_sets[setdiff(go_ids, to_drop)]
}

pick_cluster_representative <- function(cluster_indices, set_ids, set_sizes, representative_scores = NULL) {
  if (is.null(representative_scores)) {
    return(cluster_indices[[which.max(set_sizes[cluster_indices])]])
  }

  scores <- representative_scores[set_ids[cluster_indices]]
  scores[is.na(scores)] <- -Inf
  best <- cluster_indices[scores == max(scores)]

  if (length(best) > 1L) {
    best_sizes <- set_sizes[best]
    best <- best[best_sizes == min(best_sizes)]
  }

  if (length(best) > 1L) {
    best_names <- set_ids[best]
    best <- best[order(best_names)][1]
  }

  best[[1]]
}

dedup_by_jaccard <- function(gene_sets, jaccard_threshold = 0.7, representative_scores = NULL) {
  assert_probability(jaccard_threshold, "jaccard_threshold")

  gene_sets <- compact_gene_sets(gene_sets)
  set_ids <- names(gene_sets)
  n_sets <- length(set_ids)
  if (n_sets < 2L) {
    return(gene_sets)
  }

  all_genes <- unique(unlist(gene_sets, use.names = FALSE))
  gene_index <- stats::setNames(seq_along(all_genes), all_genes)

  row_indices <- integer(0)
  col_indices <- integer(0)
  for (col_idx in seq_along(set_ids)) {
    genes <- gene_sets[[set_ids[[col_idx]]]]
    row_indices <- c(row_indices, gene_index[genes])
    col_indices <- c(col_indices, rep(col_idx, length(genes)))
  }

  incidence_matrix <- Matrix::sparseMatrix(
    i = row_indices,
    j = col_indices,
    x = 1,
    dims = c(length(all_genes), n_sets),
    dimnames = list(all_genes, set_ids)
  )

  message(sprintf(
    "  Jaccard dedup: built %d x %d sparse matrix (%.1f%% sparse)",
    nrow(incidence_matrix),
    ncol(incidence_matrix),
    (1 - Matrix::nnzero(incidence_matrix) / (nrow(incidence_matrix) * ncol(incidence_matrix))) * 100
  ))

  intersection_matrix <- crossprod(incidence_matrix)
  set_sizes <- diff(incidence_matrix@p)
  tri_summary <- Matrix::summary(Matrix::triu(intersection_matrix, k = 1))

  if (!nrow(tri_summary)) {
    message("  Jaccard dedup: no overlapping pairs found")
    return(gene_sets)
  }

  message(sprintf("  Jaccard dedup: evaluating %d overlapping pairs", nrow(tri_summary)))
  union_sizes <- set_sizes[tri_summary$i] + set_sizes[tri_summary$j] - tri_summary$x
  jaccard_values <- tri_summary$x / union_sizes
  keep_pairs <- which(jaccard_values >= jaccard_threshold)

  if (!length(keep_pairs)) {
    message("  Jaccard dedup: no pairs above threshold")
    return(gene_sets)
  }

  message(sprintf(
    "  Jaccard dedup: %d pairs above threshold %.2f",
    length(keep_pairs),
    jaccard_threshold
  ))

  parent <- seq_len(n_sets)
  find_root <- function(x) {
    while (parent[[x]] != x) {
      parent[[x]] <<- parent[[parent[[x]]]]
      x <- parent[[x]]
    }
    x
  }
  unite <- function(a, b) {
    root_a <- find_root(a)
    root_b <- find_root(b)
    if (root_a != root_b) {
      parent[[root_a]] <<- root_b
    }
  }

  for (pair_idx in keep_pairs) {
    unite(tri_summary$i[[pair_idx]], tri_summary$j[[pair_idx]])
  }

  roots <- vapply(seq_len(n_sets), find_root, integer(1))
  clusters <- split(seq_len(n_sets), roots)
  keep_indices <- vapply(
    clusters,
    pick_cluster_representative,
    integer(1),
    set_ids = set_ids,
    set_sizes = set_sizes,
    representative_scores = representative_scores
  )

  n_dropped <- n_sets - length(keep_indices)
  message(sprintf("  Jaccard dedup: dropped %d redundant terms", n_dropped))
  gene_sets[set_ids[sort(keep_indices)]]
}

enforce_paired_survival <- function(gene_sets, all_gene_sets, pair_map) {
  survivors <- names(gene_sets)
  restored <- character(0)

  for (go_id in survivors) {
    sibling <- pair_map[[go_id]]
    if (!is.na(sibling) && !(sibling %in% survivors) && sibling %in% names(all_gene_sets)) {
      gene_sets[[sibling]] <- all_gene_sets[[sibling]]
      survivors <- c(survivors, sibling)
      restored <- c(restored, sibling)
    }
  }

  restored <- unique(restored)
  if (length(restored)) {
    message(sprintf("  Paired survival: restored %d dropped siblings", length(restored)))
  } else {
    message("  Paired survival: no siblings needed restoring")
  }

  gene_sets
}

build_go_regulation_library <- function(
  org_db,
  include_IEA = TRUE,
  pos_keywords = DEFAULT_POS_KEYWORDS,
  neg_keywords = DEFAULT_NEG_KEYWORDS,
  pair_replacements = DEFAULT_DIRECTIONAL_PAIR_REPLACEMENTS,
  min_size = 15,
  max_size = 500,
  dag_overlap_threshold = 0.8,
  jaccard_threshold = 0.7,
  output_dir = "."
) {
  assert_positive_count(min_size, "min_size")
  assert_positive_count(max_size, "max_size")
  assert_min_max(min_size, max_size, "min_size", "max_size")
  assert_probability(dag_overlap_threshold, "dag_overlap_threshold")
  assert_probability(jaccard_threshold, "jaccard_threshold")

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  evidence_tag <- if (include_IEA) "all_evidence" else "curated_only"
  organism <- organism_name_from_orgdb(org_db)
  organism_slug <- organism_slug_from_orgdb(org_db)

  message(sprintf("Organism: %s", organism))
  message("Classifying GO BP terms by directional keywords...")
  message(sprintf(
    "  Positive keywords: %s\n  Negative keywords: %s",
    paste(pos_keywords, collapse = ", "),
    paste(neg_keywords, collapse = ", ")
  ))

  terms <- classify_regulation_terms(
    pos_keywords = pos_keywords,
    neg_keywords = neg_keywords,
    pair_replacements = pair_replacements
  )

  message(sprintf(
    "  Found: %d positive, %d negative, %d unsigned",
    length(terms$positive),
    length(terms$negative),
    length(terms$unsigned)
  ))

  n_pairs <- sum(!is.na(unlist(terms$pair_map[terms$positive], use.names = FALSE)))
  message(sprintf("  Explicit sibling pairs found: %d", n_pairs))

  message("\n=== Processing SIGNED library ===")
  signed_ids <- c(terms$positive, terms$negative)
  message("Fetching gene annotations...")
  signed_genes_all <- get_go_genes(signed_ids, org_db, include_IEA = include_IEA)
  message(sprintf("  Annotated %d / %d terms", length(signed_genes_all), length(signed_ids)))

  message("Size filtering...")
  signed_genes_all <- filter_by_size(signed_genes_all, min_size, max_size)
  message(sprintf("  %d terms after size filter", length(signed_genes_all)))

  pos_gene_sets <- signed_genes_all[intersect(names(signed_genes_all), terms$positive)]
  neg_gene_sets <- signed_genes_all[intersect(names(signed_genes_all), terms$negative)]

  message("Processing positive direction terms...")
  pos_gene_sets <- prune_by_dag(pos_gene_sets, ontology = "BP", overlap_threshold = dag_overlap_threshold)
  pos_gene_sets <- dedup_by_jaccard(
    pos_gene_sets,
    jaccard_threshold = jaccard_threshold,
    representative_scores = get_go_term_specificity(names(pos_gene_sets), ontology = "BP")
  )
  message(sprintf("  %d positive terms after dedup", length(pos_gene_sets)))

  message("Processing negative direction terms...")
  neg_gene_sets <- prune_by_dag(neg_gene_sets, ontology = "BP", overlap_threshold = dag_overlap_threshold)
  neg_gene_sets <- dedup_by_jaccard(
    neg_gene_sets,
    jaccard_threshold = jaccard_threshold,
    representative_scores = get_go_term_specificity(names(neg_gene_sets), ontology = "BP")
  )
  message(sprintf("  %d negative terms after dedup", length(neg_gene_sets)))

  signed_gene_sets <- c(pos_gene_sets, neg_gene_sets)
  message("Enforcing paired survival...")
  signed_gene_sets <- enforce_paired_survival(signed_gene_sets, signed_genes_all, terms$pair_map)

  signed_file <- file.path(
    output_dir,
    sprintf("go_%s_bp_dedup_signed_%s.gmt", organism_slug, evidence_tag)
  )
  write_gmt_from_terms(signed_gene_sets, terms$term_df, signed_file)

  message("\n=== Processing UNSIGNED library ===")
  message("Fetching gene annotations...")
  unsigned_gene_sets <- get_go_genes(terms$unsigned, org_db, include_IEA = include_IEA)
  message(sprintf("  Annotated %d / %d terms", length(unsigned_gene_sets), length(terms$unsigned)))

  message("Size filtering...")
  unsigned_gene_sets <- filter_by_size(unsigned_gene_sets, min_size, max_size)
  message(sprintf("  %d terms after size filter", length(unsigned_gene_sets)))

  message("DAG pruning...")
  unsigned_gene_sets <- prune_by_dag(unsigned_gene_sets, ontology = "BP", overlap_threshold = dag_overlap_threshold)
  message(sprintf("  %d terms after DAG pruning", length(unsigned_gene_sets)))

  message("Jaccard deduplication...")
  unsigned_gene_sets <- dedup_by_jaccard(
    unsigned_gene_sets,
    jaccard_threshold = jaccard_threshold,
    representative_scores = get_go_term_specificity(names(unsigned_gene_sets), ontology = "BP")
  )
  message(sprintf("  %d terms in final unsigned library", length(unsigned_gene_sets)))

  unsigned_file <- file.path(
    output_dir,
    sprintf("go_%s_bp_dedup_unsigned_%s.gmt", organism_slug, evidence_tag)
  )
  write_gmt_from_terms(unsigned_gene_sets, terms$term_df, unsigned_file)

  n_positive <- sum(names(signed_gene_sets) %in% terms$positive)
  n_negative <- sum(names(signed_gene_sets) %in% terms$negative)

  paired_count <- 0L
  counted <- character(0)
  for (go_id in names(signed_gene_sets)) {
    if (go_id %in% counted) {
      next
    }

    sibling <- terms$pair_map[[go_id]]
    if (!is.na(sibling) && sibling %in% names(signed_gene_sets)) {
      paired_count <- paired_count + 1L
      counted <- c(counted, go_id, sibling)
    }
  }

  orphan_count <- length(signed_gene_sets) - length(counted)
  message(sprintf(
    paste0(
      "\n=== Summary ===",
      "\nOrganism: %s",
      "\nIEA:      %s",
      "\nSigned:   %d terms (%d positive, %d negative)",
      "\n          %d complete pairs, %d orphans",
      "\nUnsigned: %d terms",
      "\nOutput:   %s"
    ),
    organism,
    if (include_IEA) "included" else "excluded",
    length(signed_gene_sets),
    n_positive,
    n_negative,
    paired_count,
    orphan_count,
    length(unsigned_gene_sets),
    output_dir
  ))

  invisible(list(
    signed = signed_gene_sets,
    unsigned = unsigned_gene_sets,
    pair_map = terms$pair_map,
    term_df = terms$term_df,
    signed_gmt = signed_file,
    unsigned_gmt = unsigned_file
  ))
}

build_go_dedup_library <- function(
  org_db,
  ontology = "CC",
  include_IEA = TRUE,
  min_size = 15,
  max_size = 500,
  dag_overlap_threshold = 0.8,
  jaccard_threshold = 0.7,
  output_dir = "."
) {
  ontology <- toupper(ontology)
  if (!(ontology %in% c("BP", "CC", "MF"))) {
    stop("ontology must be one of: BP, CC, MF", call. = FALSE)
  }

  assert_positive_count(min_size, "min_size")
  assert_positive_count(max_size, "max_size")
  assert_min_max(min_size, max_size, "min_size", "max_size")
  assert_probability(dag_overlap_threshold, "dag_overlap_threshold")
  assert_probability(jaccard_threshold, "jaccard_threshold")

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  evidence_tag <- if (include_IEA) "all_evidence" else "curated_only"
  organism <- organism_name_from_orgdb(org_db)
  organism_slug <- organism_slug_from_orgdb(org_db)
  ontology_name <- switch(
    ontology,
    "BP" = "Biological Process",
    "CC" = "Cellular Component",
    "MF" = "Molecular Function"
  )

  message(sprintf("Organism: %s", organism))
  message(sprintf("Ontology: %s (%s)", ontology, ontology_name))
  message("Fetching GO terms...")

  all_terms <- AnnotationDbi::select(
    GO.db::GO.db,
    keys = ontology,
    keytype = "ONTOLOGY",
    columns = c("GOID", "TERM")
  )
  message(sprintf("  %d total %s terms", nrow(all_terms), ontology))

  message("Fetching gene annotations...")
  gene_sets <- get_go_genes(all_terms$GOID, org_db, include_IEA = include_IEA)
  message(sprintf("  Annotated %d / %d terms", length(gene_sets), nrow(all_terms)))

  message("Size filtering...")
  gene_sets <- filter_by_size(gene_sets, min_size, max_size)
  message(sprintf("  %d terms after size filter", length(gene_sets)))

  message("DAG pruning...")
  gene_sets <- prune_by_dag(gene_sets, ontology = ontology, overlap_threshold = dag_overlap_threshold)
  message(sprintf("  %d terms after DAG pruning", length(gene_sets)))

  message("Jaccard deduplication...")
  gene_sets <- dedup_by_jaccard(
    gene_sets,
    jaccard_threshold = jaccard_threshold,
    representative_scores = get_go_term_specificity(names(gene_sets), ontology = ontology)
  )
  message(sprintf("  %d terms in final library", length(gene_sets)))

  gmt_file <- file.path(
    output_dir,
    sprintf("go_%s_%s_dedup_%s.gmt", organism_slug, tolower(ontology), evidence_tag)
  )
  write_gmt_from_terms(gene_sets, all_terms, gmt_file)

  message(sprintf(
    paste0(
      "\n=== Summary ===",
      "\nOrganism:  %s",
      "\nOntology:  %s (%s)",
      "\nIEA:       %s",
      "\nGene sets: %d",
      "\nOutput:    %s"
    ),
    organism,
    ontology,
    ontology_name,
    if (include_IEA) "included" else "excluded",
    length(gene_sets),
    gmt_file
  ))

  invisible(list(
    gene_sets = gene_sets,
    term_df = all_terms,
    gmt_file = gmt_file
  ))
}
