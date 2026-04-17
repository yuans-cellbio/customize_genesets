load_pathway_libraries <- function(pathway_dir = "pathway_library") {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Install clusterProfiler to load GMT files", call. = FALSE)
  }

  signed_paths <- list.files(pathway_dir, pattern = ".*_signed.*\\.gmt$", full.names = TRUE)
  unsigned_paths <- setdiff(
    list.files(pathway_dir, pattern = "\\.gmt$", full.names = TRUE),
    signed_paths
  )

  load_single_library <- function(filepath) {
    list(
      lib_name = sub("\\.gmt$", "", basename(filepath)),
      term2gene = clusterProfiler::read.gmt(filepath),
      term2name = read_gmt_term2name(filepath)
    )
  }

  list(
    signed = lapply(signed_paths, load_single_library),
    unsigned = lapply(unsigned_paths, load_single_library)
  )
}

validate_signed_rank <- function(signed_rank, arg_name = "signed_rank") {
  if (!is.numeric(signed_rank)) {
    stop(sprintf("%s must be a numeric vector", arg_name), call. = FALSE)
  }
  if (is.null(names(signed_rank))) {
    stop(sprintf("%s must be a named vector of Entrez IDs", arg_name), call. = FALSE)
  }

  keep <- is.finite(signed_rank) & !is.na(names(signed_rank)) & nzchar(names(signed_rank))
  signed_rank <- signed_rank[keep]
  if (!length(signed_rank)) {
    stop(sprintf("%s does not contain any usable genes", arg_name), call. = FALSE)
  }

  signed_rank[!duplicated(names(signed_rank))]
}

validate_gene_ids <- function(gene_ids, arg_name) {
  if (!is.atomic(gene_ids)) {
    stop(sprintf("%s must be an atomic vector of gene IDs", arg_name), call. = FALSE)
  }

  gene_ids <- as.character(gene_ids)
  gene_ids <- gene_ids[!is.na(gene_ids) & nzchar(gene_ids)]
  gene_ids <- unique(gene_ids)
  if (!length(gene_ids)) {
    stop(sprintf("%s does not contain any usable genes", arg_name), call. = FALSE)
  }

  gene_ids
}

select_ranked_genes_for_ora <- function(signed_rank, gene_fraction = 0.1) {
  signed_rank <- validate_signed_rank(signed_rank)
  assert_probability(gene_fraction, "gene_fraction")
  if (gene_fraction <= 0) {
    stop("gene_fraction must be greater than 0", call. = FALSE)
  }

  abs_rank <- sort(abs(signed_rank), decreasing = TRUE)
  top_n <- min(length(abs_rank), max(1L, as.integer(ceiling(length(abs_rank) * gene_fraction))))
  selected <- names(abs_rank)[seq_len(top_n)]

  list(
    selected = selected,
    up = names(signed_rank)[signed_rank > 0 & names(signed_rank) %in% selected],
    down = names(signed_rank)[signed_rank < 0 & names(signed_rank) %in% selected]
  )
}

resolve_ora_gene_input <- function(
    signed_rank = NULL,
    selected_genes = NULL,
    gene_fraction = NULL,
    universe = NULL
) {
  using_ranked_input <- !is.null(signed_rank)
  using_selected_input <- !is.null(selected_genes)

  if (using_ranked_input == using_selected_input) {
    stop("Provide exactly one of signed_rank or selected_genes", call. = FALSE)
  }

  if (using_ranked_input) {
    signed_rank <- validate_signed_rank(signed_rank, arg_name = "signed_rank")
    gene_fraction <- gene_fraction %||% 0.1
    selected <- select_ranked_genes_for_ora(signed_rank, gene_fraction = gene_fraction)
    universe <- if (is.null(universe)) names(signed_rank) else validate_gene_ids(universe, "universe")
    if (!all(selected$selected %in% universe)) {
      stop("universe must include all selected genes", call. = FALSE)
    }

    return(list(
      mode = "ranked",
      selected = selected,
      universe = universe
    ))
  }

  selected_genes <- validate_signed_rank(selected_genes, arg_name = "selected_genes")
  if (!is.null(gene_fraction)) {
    assert_probability(gene_fraction, "gene_fraction")
    if (!isTRUE(all.equal(gene_fraction, 1))) {
      stop("gene_fraction must be 1 when selected_genes is supplied", call. = FALSE)
    }
  }
  if (is.null(universe)) {
    stop("universe must be supplied when selected_genes is used", call. = FALSE)
  }

  universe <- validate_gene_ids(universe, "universe")
  selected <- list(
    selected = names(selected_genes),
    up = names(selected_genes)[selected_genes > 0],
    down = names(selected_genes)[selected_genes < 0]
  )
  if (!all(selected$selected %in% universe)) {
    stop("universe must include all selected genes", call. = FALSE)
  }

  list(
    mode = "selected",
    selected = selected,
    universe = universe
  )
}

prepare_gsea_rankings <- function(signed_rank) {
  signed_rank <- validate_signed_rank(signed_rank)
  list(
    signed = sort(signed_rank, decreasing = TRUE),
    unsigned_magnitude = sort(abs(signed_rank), decreasing = TRUE)
  )
}

format_enrichment_result <- function(result, org_db, direction, library_name) {
  if (is.null(result)) {
    return(NULL)
  }

  result_df <- as.data.frame(result)
  if (!nrow(result_df)) {
    return(NULL)
  }

  readable_result <- tryCatch(
    clusterProfiler::setReadable(result, org_db, "ENTREZID"),
    error = function(e) result
  )
  result_df <- as.data.frame(readable_result)
  if (!nrow(result_df)) {
    return(NULL)
  }

  result_df$Direction <- direction
  result_df$Library <- library_name
  result_df[, c("Library", "Direction", setdiff(names(result_df), c("Library", "Direction"))), drop = FALSE]
}

adjust_enrichment_results <- function(result_df, p_adjust_method = "BH") {
  if (is.null(result_df) || !nrow(result_df)) {
    return(result_df)
  }

  group_key <- interaction(result_df$Library, result_df$Direction, drop = TRUE)
  result_df$p.adjust <- stats::ave(
    result_df$pvalue,
    group_key,
    FUN = function(x) stats::p.adjust(x, method = p_adjust_method)
  )
  result_df
}

run_ora <- function(
    signed_rank = NULL,
    selected_genes = NULL,
    gene_fraction = NULL,
    org_db,
    universe = NULL,
    signed_library,
    unsigned_library,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 15,
    maxGSSize = 500
) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Install clusterProfiler to run ORA", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Install dplyr to combine ORA results", call. = FALSE)
  }

  ora_input <- resolve_ora_gene_input(
    signed_rank = signed_rank,
    selected_genes = selected_genes,
    gene_fraction = gene_fraction,
    universe = universe
  )
  selected <- ora_input$selected
  universe <- ora_input$universe

  run_single_ora <- function(library_entry, gene_ids, direction_label) {
    if (!length(gene_ids)) {
      return(NULL)
    }

    result <- clusterProfiler::enricher(
      gene = gene_ids,
      universe = universe,
      TERM2GENE = library_entry$term2gene,
      TERM2NAME = library_entry$term2name,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      pAdjustMethod = "none",
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    )

    format_enrichment_result(
      result,
      org_db = org_db,
      direction = direction_label,
      library_name = library_entry$lib_name
    )
  }

  signed_results <- c(
    lapply(signed_library, run_single_ora, gene_ids = selected$up, direction_label = "up"),
    lapply(signed_library, run_single_ora, gene_ids = selected$down, direction_label = "down")
  )
  unsigned_results <- lapply(
    unsigned_library,
    run_single_ora,
    gene_ids = selected$selected,
    direction_label = "unsigned"
  )

  combined <- dplyr::bind_rows(c(signed_results, unsigned_results))
  adjust_enrichment_results(combined, p_adjust_method = pAdjustMethod)
}

run_gsea <- function(
    signed_rank,
    org_db,
    signed_library,
    unsigned_library,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 15,
    maxGSSize = 500,
    nPermSimple = 10000
) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Install clusterProfiler to run GSEA", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Install dplyr to combine GSEA results", call. = FALSE)
  }

  rankings <- prepare_gsea_rankings(signed_rank)

  run_single_gsea <- function(library_entry, gene_list, direction_label, score_type) {
    result <- clusterProfiler::GSEA(
      geneList = gene_list,
      TERM2GENE = library_entry$term2gene,
      TERM2NAME = library_entry$term2name,
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = "none",
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      eps = 0,
      nPermSimple = nPermSimple,
      scoreType = score_type,
      seed = TRUE
    )

    format_enrichment_result(
      result,
      org_db = org_db,
      direction = direction_label,
      library_name = library_entry$lib_name
    )
  }

  signed_results <- lapply(
    signed_library,
    run_single_gsea,
    gene_list = rankings$signed,
    direction_label = "signed",
    score_type = "std"
  )
  unsigned_results <- lapply(
    unsigned_library,
    run_single_gsea,
    gene_list = rankings$unsigned_magnitude,
    direction_label = "unsigned_magnitude",
    score_type = "pos"
  )

  combined <- dplyr::bind_rows(c(signed_results, unsigned_results))
  adjust_enrichment_results(combined, p_adjust_method = pAdjustMethod)
}
