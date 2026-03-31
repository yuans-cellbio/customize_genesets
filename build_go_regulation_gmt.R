#!/usr/bin/env Rscript
# =============================================================================
# build_go_regulation_gmt.R
#
# Build deduplicated gene set libraries from Gene Ontology.
#
# Two main functions:
#   1. build_go_regulation_library() — For BP: produces signed (directional)
#      and unsigned GMT files with paired survival enforcement.
#   2. build_go_dedup_library() — For CC, MF, or BP without direction:
#      produces a single deduplicated GMT file.
#
# Deduplication pipeline: size filter -> DAG pruning -> sparse Jaccard dedup
#
# Dependencies:
#   BiocManager::install(c("GO.db", "AnnotationDbi"))
#   install.packages("Matrix")  # ships with base R
#   Plus an OrgDb package, e.g. org.Hs.eg.db (human) or org.Mm.eg.db (mouse)
# =============================================================================

suppressPackageStartupMessages({
  library(GO.db)
  library(AnnotationDbi)
  library(Matrix)
})

# =============================================================================
# 1. DEFAULT DIRECTIONAL KEYWORDS
# =============================================================================

# Keywords indicating positive/activating direction
DEFAULT_POS_KEYWORDS <- c(
  "positive regulation", "activation of", "upregulation of",
  "induction of", "stimulation of", "promotion of",
  "enhancement of", "potentiation of"
)

# Keywords indicating negative/inhibiting direction
DEFAULT_NEG_KEYWORDS <- c(
  "negative regulation", "inhibition of", "downregulation of",
  "suppression of", "repression of", "attenuation of",
  "blockade of", "sequestering of", "desensitization of",
  "deactivation of", "inactivation of", "degradation of"
)

# =============================================================================
# 2. KEYWORD DISCOVERY HELPER
# =============================================================================

#' Scan all GO BP terms for potential directional keywords not yet in your lists.
#'
#' @param pos_keywords character vector of known positive keywords
#' @param neg_keywords character vector of known negative keywords
#' @param candidate_words character vector of words to check. If NULL, uses built-in list.
#' @return data.frame with columns: keyword, n_terms, already_known, example_terms
discover_directional_keywords <- function(
  pos_keywords = DEFAULT_POS_KEYWORDS,
  neg_keywords = DEFAULT_NEG_KEYWORDS,
  candidate_words = NULL
) {
  bp_terms <- AnnotationDbi::select(GO.db,
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

  for (kw in candidate_words) {
    matches <- grep(tolower(kw), terms_lower, fixed = TRUE)
    if (length(matches) > 0) {
      examples <- head(bp_terms$TERM[matches], 3)
      results <- rbind(results, data.frame(
        keyword = kw,
        n_terms = length(matches),
        already_known = tolower(kw) %in% tolower(known),
        example_terms = paste(examples, collapse = " | "),
        stringsAsFactors = FALSE
      ))
    }
  }

  results <- results[order(-results$n_terms), ]
  rownames(results) <- NULL

  message(sprintf(
    "Found %d candidate keywords with GO BP matches (%d already in your lists)",
    nrow(results), sum(results$already_known)
  ))

  new_kw <- results[!results$already_known, ]
  if (nrow(new_kw) > 0) {
    message("\nKeywords NOT yet in your lists:")
    for (i in seq_len(nrow(new_kw))) {
      message(sprintf(
        '  "%s" (%d terms) e.g.: %s',
        new_kw$keyword[i], new_kw$n_terms[i], new_kw$example_terms[i]
      ))
    }
  }

  results
}

# =============================================================================
# 3. CLASSIFY GO TERMS
# =============================================================================

#' Find the position of the earliest (leftmost) keyword match in a string.
#' @param term character(1) lowercase term name
#' @param keywords character vector of keywords
#' @return integer position of earliest match, or Inf
earliest_match_pos <- function(term, keywords) {
  positions <- vapply(keywords, function(kw) {
    m <- regexpr(kw, term, fixed = TRUE)
    if (m == -1L) Inf else as.double(m)
  }, double(1))
  min(positions)
}

#' Classify GO BP terms into positive, negative, or unsigned.
#'
#' When a term contains both positive and negative keywords,
#' the leftmost keyword determines the direction.
#'
#' @param pos_keywords character vector of positive-direction keywords
#' @param neg_keywords character vector of negative-direction keywords
#' @return list with: positive, negative, unsigned (GO ID vectors),
#'   pair_map (named list: signed GO ID -> sibling GO ID or NA),
#'   term_df (data.frame of all BP terms)
classify_regulation_terms <- function(
  pos_keywords = DEFAULT_POS_KEYWORDS,
  neg_keywords = DEFAULT_NEG_KEYWORDS
) {
  bp_terms <- AnnotationDbi::select(GO.db,
    keys = "BP",
    keytype = "ONTOLOGY",
    columns = c("GOID", "TERM")
  )

  terms_lower <- tolower(bp_terms$TERM)
  pos_keywords_lower <- tolower(pos_keywords)
  neg_keywords_lower <- tolower(neg_keywords)

  direction <- character(nrow(bp_terms))

  for (i in seq_along(terms_lower)) {
    pos_pos <- earliest_match_pos(terms_lower[i], pos_keywords_lower)
    neg_pos <- earliest_match_pos(terms_lower[i], neg_keywords_lower)

    if (is.infinite(pos_pos) && is.infinite(neg_pos)) {
      direction[i] <- "unsigned"
    } else if (pos_pos < neg_pos) {
      direction[i] <- "positive"
    } else if (neg_pos < pos_pos) {
      direction[i] <- "negative"
    } else {
      direction[i] <- "unsigned"
      message(sprintf("  Tie at position %d for term: %s", pos_pos, bp_terms$TERM[i]))
    }
  }

  pos_ids <- bp_terms$GOID[direction == "positive"]
  neg_ids <- bp_terms$GOID[direction == "negative"]
  unsigned_ids <- bp_terms$GOID[direction == "unsigned"]

  # --- Build sibling pair map ---
  term_lookup <- setNames(bp_terms$TERM, bp_terms$GOID)
  name_to_id <- setNames(bp_terms$GOID, tolower(bp_terms$TERM))
  pair_map <- list()

  n_pairs <- min(length(pos_keywords), length(neg_keywords))
  swap_pairs <- lapply(seq_len(n_pairs), function(i) {
    list(pos = pos_keywords_lower[i], neg = neg_keywords_lower[i])
  })

  find_sibling <- function(go_id, current_direction) {
    term <- tolower(term_lookup[[go_id]])
    for (sp in swap_pairs) {
      if (current_direction == "positive" && grepl(sp$pos, term, fixed = TRUE)) {
        candidate_name <- sub(sp$pos, sp$neg, term, fixed = TRUE)
        candidate_id <- name_to_id[candidate_name]
        if (!is.na(candidate_id)) return(unname(candidate_id))
      } else if (current_direction == "negative" && grepl(sp$neg, term, fixed = TRUE)) {
        candidate_name <- sub(sp$neg, sp$pos, term, fixed = TRUE)
        candidate_id <- name_to_id[candidate_name]
        if (!is.na(candidate_id)) return(unname(candidate_id))
      }
    }
    return(NA_character_)
  }

  for (pid in pos_ids) pair_map[[pid]] <- find_sibling(pid, "positive")
  for (nid in neg_ids) pair_map[[nid]] <- find_sibling(nid, "negative")

  list(
    positive = pos_ids,
    negative = neg_ids,
    unsigned = unsigned_ids,
    pair_map = pair_map,
    term_df = bp_terms
  )
}

# =============================================================================
# 4. GET GENE ANNOTATIONS
# =============================================================================

#' Map GO term IDs to Entrez gene IDs
#' @param go_ids character vector of GO IDs
#' @param org_db an OrgDb object
#' @param include_IEA logical; if FALSE, exclude IEA evidence codes
#' @return named list: GO ID -> character vector of Entrez IDs
get_go_genes <- function(go_ids, org_db, include_IEA = TRUE) {
  mapping <- AnnotationDbi::select(org_db,
    keys = go_ids,
    keytype = "GOALL",
    columns = "ENTREZID"
  )

  if (!include_IEA) {
    mapping <- mapping[mapping$EVIDENCEALL != "IEA", ]
  }

  mapping <- unique(mapping[, c("GOALL", "ENTREZID")])
  mapping <- mapping[!is.na(mapping$ENTREZID), ]

  split(mapping$ENTREZID, mapping$GOALL)
}

# =============================================================================
# 5. SIZE FILTER
# =============================================================================

#' Filter gene sets by size
#' @param gene_sets named list of character vectors
#' @param min_size minimum number of genes (default 15)
#' @param max_size maximum number of genes (default 500)
#' @return filtered named list
filter_by_size <- function(gene_sets, min_size = 15, max_size = 500) {
  sizes <- vapply(gene_sets, length, integer(1))
  gene_sets[sizes >= min_size & sizes <= max_size]
}

# =============================================================================
# 6. DAG-BASED PRUNING
# =============================================================================

#' Remove redundant ancestor terms using GO DAG structure.
#' @param gene_sets named list of character vectors (GO ID -> genes)
#' @param ontology one of "BP", "CC", "MF" (default "BP")
#' @param overlap_threshold proportion threshold (default 0.8)
#' @return pruned named list
prune_by_dag <- function(gene_sets, ontology = "BP", overlap_threshold = 0.8) {
  go_ids <- names(gene_sets)
  if (length(go_ids) < 2) return(gene_sets)

  anc_obj <- switch(ontology,
    "BP" = GOBPANCESTOR,
    "CC" = GOCCANCESTOR,
    "MF" = GOMFANCESTOR,
    stop("ontology must be one of: BP, CC, MF")
  )
  anc_env <- as.list(anc_obj)
  to_drop <- character(0)

  for (child_id in go_ids) {
    ancestors <- anc_env[[child_id]]
    if (is.null(ancestors)) next

    shared_ancestors <- intersect(ancestors, go_ids)

    for (anc_id in shared_ancestors) {
      overlap <- length(intersect(gene_sets[[child_id]], gene_sets[[anc_id]]))
      frac <- overlap / length(gene_sets[[anc_id]])

      if (frac >= overlap_threshold) {
        to_drop <- c(to_drop, anc_id)
      }
    }
  }

  to_drop <- unique(to_drop)
  message(sprintf("  DAG pruning: dropping %d redundant ancestor terms", length(to_drop)))
  gene_sets[setdiff(go_ids, to_drop)]
}

# =============================================================================
# 7. JACCARD DEDUPLICATION (SPARSE MATRIX)
# =============================================================================

#' Deduplicate gene sets by Jaccard similarity using sparse matrix computation.
#'
#' Builds a sparse binary genes x terms matrix, computes all pairwise
#' intersection sizes via a single matrix multiplication, derives Jaccard,
#' then clusters terms above the threshold and keeps the largest per cluster.
#'
#' @param gene_sets named list of character vectors
#' @param jaccard_threshold threshold above which to merge (default 0.7)
#' @return deduplicated named list
dedup_by_jaccard <- function(gene_sets, jaccard_threshold = 0.7) {
  go_ids <- names(gene_sets)
  n_sets <- length(go_ids)
  if (n_sets < 2) return(gene_sets)

  # Build sparse binary matrix: rows = genes, columns = GO terms
  all_genes <- unique(unlist(gene_sets, use.names = FALSE))
  gene_idx <- setNames(seq_along(all_genes), all_genes)

  # Construct triplets for sparse matrix
  col_indices <- integer(0)
  row_indices <- integer(0)

  for (j in seq_along(go_ids)) {
    genes <- gene_sets[[go_ids[j]]]
    rows <- gene_idx[genes]
    row_indices <- c(row_indices, rows)
    col_indices <- c(col_indices, rep(j, length(rows)))
  }

  M <- sparseMatrix(
    i = row_indices,
    j = col_indices,
    x = 1,
    dims = c(length(all_genes), n_sets),
    dimnames = list(all_genes, go_ids)
  )

  message(sprintf(
    "  Jaccard dedup: built %d x %d sparse matrix (%.1f%% sparse)",
    nrow(M), ncol(M), (1 - nnzero(M) / (nrow(M) * ncol(M))) * 100
  ))

  # Intersection matrix: t(M) %*% M gives pairwise intersection sizes
  # This is the key operation — handled by optimized sparse algebra
  intersection_mat <- crossprod(M)  # equivalent to t(M) %*% M but faster

  # Set sizes
  set_sizes <- diff(M@p)  # column counts for dgCMatrix (fast)

  # Convert intersection to Jaccard (only for upper triangle)
  # Union(A,B) = |A| + |B| - |A ∩ B|
  # Jaccard = |A ∩ B| / (|A| + |B| - |A ∩ B|)

  # Extract upper triangle nonzero entries (pairs with nonzero intersection)
  inter_tri <- triu(intersection_mat, k = 1)  # strict upper triangle

  # Get the (i, j, x) triplets from the upper triangle
  tri_summary <- summary(inter_tri)

  if (nrow(tri_summary) == 0) {
    message("  Jaccard dedup: no overlapping pairs found")
    return(gene_sets)
  }

  message(sprintf("  Jaccard dedup: evaluating %d overlapping pairs", nrow(tri_summary)))

  # Compute Jaccard for each pair
  union_sizes <- set_sizes[tri_summary$i] + set_sizes[tri_summary$j] - tri_summary$x
  jaccard_vals <- tri_summary$x / union_sizes

  # Filter to pairs above threshold
  above <- which(jaccard_vals >= jaccard_threshold)

  if (length(above) == 0) {
    message("  Jaccard dedup: no pairs above threshold")
    return(gene_sets)
  }

  message(sprintf("  Jaccard dedup: %d pairs above threshold %.2f", length(above), jaccard_threshold))

  # Union-find clustering on the pairs above threshold
  merge_i <- tri_summary$i[above]
  merge_j <- tri_summary$j[above]

  parent <- seq_len(n_sets)
  find <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  unite <- function(a, b) {
    ra <- find(a)
    rb <- find(b)
    if (ra != rb) parent[ra] <<- rb
  }

  for (k in seq_along(merge_i)) {
    unite(merge_i[k], merge_j[k])
  }

  # Find root for each node
  roots <- vapply(seq_len(n_sets), find, integer(1))
  clusters <- split(seq_len(n_sets), roots)

  # Keep largest set per cluster
  keep_idx <- vapply(clusters, function(cl) {
    cl[which.max(set_sizes[cl])]
  }, integer(1))

  n_dropped <- n_sets - length(keep_idx)
  message(sprintf("  Jaccard dedup: dropped %d redundant terms", n_dropped))
  gene_sets[go_ids[sort(keep_idx)]]
}

# =============================================================================
# 8. ENFORCE PAIRED SURVIVAL
# =============================================================================

#' After dedup, force-include siblings of surviving signed terms.
#'
#' @param gene_sets named list of surviving signed gene sets
#' @param all_gene_sets named list of ALL signed gene sets (pre-dedup, post size filter)
#' @param pair_map named list mapping each GO ID to its sibling GO ID (or NA)
#' @return gene_sets with missing siblings restored
enforce_paired_survival <- function(gene_sets, all_gene_sets, pair_map) {
  survivors <- names(gene_sets)
  restored <- character(0)

  for (id in survivors) {
    sibling <- pair_map[[id]]
    if (!is.na(sibling) && !(sibling %in% survivors) && sibling %in% names(all_gene_sets)) {
      gene_sets[[sibling]] <- all_gene_sets[[sibling]]
      survivors <- c(survivors, sibling)
      restored <- c(restored, sibling)
    }
  }

  restored <- unique(restored)
  if (length(restored) > 0) {
    message(sprintf("  Paired survival: restored %d dropped siblings", length(restored)))
  } else {
    message("  Paired survival: no siblings needed restoring")
  }

  gene_sets
}

# =============================================================================
# 9. GMT EXPORT AND IMPORT HELPERS
# =============================================================================

#' Write gene sets to GMT format
#' Column order: GO_ID (term identifier), term_name (description), genes...
#' This ensures clusterProfiler::read.gmt() uses GO IDs as term identifiers.
#'
#' @param gene_sets named list of character vectors (GO ID -> genes)
#' @param term_df data.frame with GOID and TERM columns
#' @param filepath output file path
write_gmt <- function(gene_sets, term_df, filepath) {
  term_lookup <- setNames(term_df$TERM, term_df$GOID)

  lines <- vapply(names(gene_sets), function(id) {
    term_name <- term_lookup[[id]]
    if (is.na(term_name)) term_name <- id
    genes <- gene_sets[[id]]
    paste(c(id, term_name, genes), collapse = "\t")
  }, character(1))

  writeLines(lines, filepath)
  message(sprintf("Wrote %d gene sets to %s", length(gene_sets), filepath))
}

#' Read a TERM2NAME mapping from a GMT file.
#' Extracts column 1 (GO ID) and column 2 (term name) for use with
#' clusterProfiler::enricher() and GSEA().
#'
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
# 10. MAIN PIPELINE
# =============================================================================

#' Build signed and unsigned GO regulation gene set libraries
#'
#' @param org_db an OrgDb annotation object (e.g. org.Hs.eg.db, org.Mm.eg.db)
#' @param include_IEA logical; include electronically inferred annotations? (default TRUE)
#' @param pos_keywords character vector of positive-direction keywords
#' @param neg_keywords character vector of negative-direction keywords
#' @param min_size minimum gene set size (default 15)
#' @param max_size maximum gene set size (default 500)
#' @param dag_overlap_threshold DAG pruning overlap threshold (default 0.8)
#' @param jaccard_threshold Jaccard dedup threshold (default 0.7)
#' @param output_dir directory for GMT output files
#' @return invisible list of the gene set collections and metadata
build_go_regulation_library <- function(
  org_db,
  include_IEA = TRUE,
  pos_keywords = DEFAULT_POS_KEYWORDS,
  neg_keywords = DEFAULT_NEG_KEYWORDS,
  min_size = 15,
  max_size = 500,
  dag_overlap_threshold = 0.8,
  jaccard_threshold = 0.7,
  output_dir = "."
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  iea_tag <- ifelse(include_IEA, "all_evidence", "curated_only")

  # Detect organism
  org_name <- tryCatch(
    AnnotationDbi::metadata(org_db)[
      AnnotationDbi::metadata(org_db)$name == "ORGANISM", "value"
    ],
    error = function(e) "unknown"
  )
  message(sprintf("Organism: %s", org_name))

  # --- Classify terms ---
  message("Classifying GO BP terms by directional keywords...")
  message(sprintf(
    "  Positive keywords: %s\n  Negative keywords: %s",
    paste(pos_keywords, collapse = ", "),
    paste(neg_keywords, collapse = ", ")
  ))

  terms <- classify_regulation_terms(pos_keywords, neg_keywords)
  message(sprintf(
    "  Found: %d positive, %d negative, %d unsigned",
    length(terms$positive), length(terms$negative), length(terms$unsigned)
  ))

  n_paired_raw <- sum(!is.na(unlist(terms$pair_map[terms$positive])))
  message(sprintf(
    "  Paired terms (have a sibling in opposite direction): %d",
    n_paired_raw
  ))

  # =========================================================================
  # SIGNED LIBRARY
  # =========================================================================
  message("\n=== Processing SIGNED library ===")
  signed_ids <- c(terms$positive, terms$negative)

  message("Fetching gene annotations...")
  signed_genes_all <- get_go_genes(signed_ids, org_db, include_IEA = include_IEA)
  message(sprintf("  Annotated %d / %d terms", length(signed_genes_all), length(signed_ids)))

  message("Size filtering...")
  signed_genes_all <- filter_by_size(signed_genes_all, min_size, max_size)
  message(sprintf("  %d terms after size filter", length(signed_genes_all)))

  pos_ids_surviving <- intersect(names(signed_genes_all), terms$positive)
  neg_ids_surviving <- intersect(names(signed_genes_all), terms$negative)

  message("Processing positive direction terms...")
  pos_genes <- signed_genes_all[pos_ids_surviving]
  pos_genes <- prune_by_dag(pos_genes, ontology = "BP", overlap_threshold = dag_overlap_threshold)
  pos_genes <- dedup_by_jaccard(pos_genes, jaccard_threshold)
  message(sprintf("  %d positive terms after dedup", length(pos_genes)))

  message("Processing negative direction terms...")
  neg_genes <- signed_genes_all[neg_ids_surviving]
  neg_genes <- prune_by_dag(neg_genes, ontology = "BP", overlap_threshold = dag_overlap_threshold)
  neg_genes <- dedup_by_jaccard(neg_genes, jaccard_threshold)
  message(sprintf("  %d negative terms after dedup", length(neg_genes)))

  signed_genes <- c(pos_genes, neg_genes)
  message("Enforcing paired survival...")
  signed_genes <- enforce_paired_survival(
    signed_genes, signed_genes_all, terms$pair_map
  )
  
  signed_file <- file.path(output_dir, sprintf("go_%s_bp_dedup_signed_%s.gmt", tolower(gsub(" ", "_", org_name)), iea_tag))
  write_gmt(signed_genes, terms$term_df, signed_file)

  # =========================================================================
  # UNSIGNED LIBRARY
  # =========================================================================
  message("\n=== Processing UNSIGNED library ===")

  message("Fetching gene annotations...")
  unsigned_genes <- get_go_genes(terms$unsigned, org_db, include_IEA = include_IEA)
  message(sprintf("  Annotated %d / %d terms", length(unsigned_genes), length(terms$unsigned)))

  message("Size filtering...")
  unsigned_genes <- filter_by_size(unsigned_genes, min_size, max_size)
  message(sprintf("  %d terms after size filter", length(unsigned_genes)))

  message("DAG pruning...")
  unsigned_genes <- prune_by_dag(unsigned_genes, ontology = "BP", overlap_threshold = dag_overlap_threshold)
  message(sprintf("  %d terms after DAG pruning", length(unsigned_genes)))

  message("Jaccard deduplication...")
  unsigned_genes <- dedup_by_jaccard(unsigned_genes, jaccard_threshold)
  message(sprintf("  %d terms in final unsigned library", length(unsigned_genes)))
  
  unsigned_file <- file.path(output_dir, sprintf("go_%s_bp_dedup_unsigned_%s.gmt", tolower(gsub(" ", "_", org_name)), iea_tag))
  write_gmt(unsigned_genes, terms$term_df, unsigned_file)

  # =========================================================================
  # SUMMARY
  # =========================================================================
  n_pos <- sum(names(signed_genes) %in% terms$positive)
  n_neg <- sum(names(signed_genes) %in% terms$negative)

  n_paired <- 0
  counted <- character(0)
  for (id in names(signed_genes)) {
    if (id %in% counted) next
    sib <- terms$pair_map[[id]]
    if (!is.na(sib) && sib %in% names(signed_genes)) {
      n_paired <- n_paired + 1
      counted <- c(counted, id, sib)
    }
  }
  n_orphan <- length(signed_genes) - length(counted)

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
    org_name,
    ifelse(include_IEA, "included", "excluded"),
    length(signed_genes), n_pos, n_neg,
    n_paired, n_orphan,
    length(unsigned_genes),
    output_dir
  ))

  invisible(list(
    signed = signed_genes,
    unsigned = unsigned_genes,
    pair_map = terms$pair_map,
    term_df = terms$term_df,
    signed_gmt = signed_file,
    unsigned_gmt = unsigned_file
  ))
}

# =============================================================================
# 11. GENERAL GO DEDUP LIBRARY (CC, MF, or BP without direction)
# =============================================================================

#' Build a deduplicated gene set library from any GO ontology.
#'
#' This is a simplified pipeline without directional classification or paired
#' survival — suitable for Cellular Component, Molecular Function, or any
#' GO ontology where direction is not applicable.
#'
#' Pipeline: get all terms -> annotate genes -> size filter -> DAG prune ->
#'           Jaccard dedup -> write GMT
#'
#' @param org_db an OrgDb annotation object (e.g. org.Hs.eg.db, org.Mm.eg.db)
#' @param ontology one of "CC", "MF", or "BP" (default "CC")
#' @param include_IEA logical; include electronically inferred annotations? (default TRUE)
#' @param min_size minimum gene set size (default 15)
#' @param max_size maximum gene set size (default 500)
#' @param dag_overlap_threshold DAG pruning overlap threshold (default 0.8)
#' @param jaccard_threshold Jaccard dedup threshold (default 0.7)
#' @param output_dir directory for GMT output files
#' @return invisible list with gene_sets, term_df, and gmt_file path
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
  stopifnot(ontology %in% c("BP", "CC", "MF"))

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  iea_tag <- ifelse(include_IEA, "all_evidence", "curated_only")

  ontology_name <- switch(ontology,
    "BP" = "Biological Process",
    "CC" = "Cellular Component",
    "MF" = "Molecular Function"
  )

  # Detect organism
  org_name <- tryCatch(
    AnnotationDbi::metadata(org_db)[
      AnnotationDbi::metadata(org_db)$name == "ORGANISM", "value"
    ],
    error = function(e) "unknown"
  )
  message(sprintf("Organism: %s", org_name))
  message(sprintf("Ontology: %s (%s)", ontology, ontology_name))

  # --- Get all terms for this ontology ---
  message("Fetching GO terms...")
  all_terms <- AnnotationDbi::select(GO.db,
    keys = ontology,
    keytype = "ONTOLOGY",
    columns = c("GOID", "TERM")
  )
  message(sprintf("  %d total %s terms", nrow(all_terms), ontology))

  go_ids <- all_terms$GOID

  # --- Annotate genes ---
  message("Fetching gene annotations...")
  gene_sets <- get_go_genes(go_ids, org_db, include_IEA = include_IEA)
  message(sprintf("  Annotated %d / %d terms", length(gene_sets), length(go_ids)))

  # --- Size filter ---
  message("Size filtering...")
  gene_sets <- filter_by_size(gene_sets, min_size, max_size)
  message(sprintf("  %d terms after size filter", length(gene_sets)))

  # --- DAG pruning ---
  message("DAG pruning...")
  gene_sets <- prune_by_dag(gene_sets, ontology = ontology, overlap_threshold = dag_overlap_threshold)
  message(sprintf("  %d terms after DAG pruning", length(gene_sets)))

  # --- Jaccard dedup ---
  message("Jaccard deduplication...")
  gene_sets <- dedup_by_jaccard(gene_sets, jaccard_threshold)
  message(sprintf("  %d terms in final library", length(gene_sets)))

  # --- Write GMT ---
  gmt_file <- file.path(
    output_dir,
    sprintf("go_%s_%s_dedup_%s.gmt", tolower(gsub(" ", "_", org_name)), tolower(ontology), iea_tag)
  )
  write_gmt(gene_sets, all_terms, gmt_file)

  # --- Summary ---
  message(sprintf(
    paste0(
      "\n=== Summary ===",
      "\nOrganism:  %s",
      "\nOntology:  %s (%s)",
      "\nIEA:       %s",
      "\nGene sets: %d",
      "\nOutput:    %s"
    ),
    org_name, ontology, ontology_name,
    ifelse(include_IEA, "included", "excluded"),
    length(gene_sets),
    gmt_file
  ))

  invisible(list(
    gene_sets = gene_sets,
    term_df = all_terms,
    gmt_file = gmt_file
  ))
}

# =============================================================================
# RUN EXAMPLE
# =============================================================================

if (sys.nframe() == 0) {
  # --- Optional: discover additional directional keywords ---
  # discovery <- discover_directional_keywords()
  # print(discovery)

  library(org.Hs.eg.db)

  # --- Build signed/unsigned BP library ---
  result_bp <- build_go_regulation_library(
    org_db = org.Hs.eg.db,
    include_IEA = TRUE,
    pos_keywords = DEFAULT_POS_KEYWORDS,
    neg_keywords = DEFAULT_NEG_KEYWORDS,
    min_size = 15,
    max_size = 500,
    output_dir = "go_dedup_gmt"
  )

  # --- Build deduplicated CC library ---
  result_cc <- build_go_dedup_library(
    org_db = org.Hs.eg.db,
    ontology = "CC",
    include_IEA = TRUE,
    min_size = 15,
    max_size = 500,
    output_dir = "go_dedup_gmt"
  )

  # --- Build deduplicated MF library ---
  result_mf <- build_go_dedup_library(
    org_db = org.Hs.eg.db,
    ontology = "MF",
    include_IEA = TRUE,
    min_size = 15,
    max_size = 500,
    output_dir = "go_dedup_gmt"
  )

  # --- Usage with clusterProfiler ---
  # library(clusterProfiler)
  #
  # # Load gene sets and term names
  # signed_gmt  <- read.gmt(result_bp$signed_gmt)
  # signed_names <- read_gmt_term2name(result_bp$signed_gmt)
  #
  # cc_gmt   <- read.gmt(result_cc$gmt_file)
  # cc_names <- read_gmt_term2name(result_cc$gmt_file)
  #
  # # ORA
  # ora_result <- enricher(
  #   gene      = my_genes,
  #   universe  = background,
  #   TERM2GENE = signed_gmt,
  #   TERM2NAME = signed_names
  # )
  #
  # # GSEA
  # gsea_result <- GSEA(
  #   geneList  = my_ranked_genes,
  #   TERM2GENE = signed_gmt,
  #   TERM2NAME = signed_names
  # )

  # Mouse example (uncomment to use):
  # library(org.Mm.eg.db)
  # result_mouse <- build_go_regulation_library(
  #   org_db = org.Mm.eg.db,
  #   include_IEA = FALSE,
  #   output_dir = "go_regulation_gmt_mouse"
  # )
  # result_mouse_cc <- build_go_dedup_library(
  #   org_db = org.Mm.eg.db,
  #   ontology = "CC",
  #   output_dir = "go_regulation_gmt_mouse"
  # )
}
