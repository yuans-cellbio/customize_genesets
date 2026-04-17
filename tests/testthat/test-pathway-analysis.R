test_that("ORA selection uses an exact top-n fraction", {
  signed_rank <- c(A = 5, B = -4, C = 3, D = -2, E = 1)

  selected <- select_ranked_genes_for_ora(signed_rank, gene_fraction = 0.4)

  expect_setequal(selected$selected, c("A", "B"))
  expect_equal(selected$up, "A")
  expect_equal(selected$down, "B")
})

test_that("ORA input resolution defaults to ranked-selection mode", {
  signed_rank <- c(A = 5, B = -4, C = 3, D = -2, E = 1)

  resolved <- resolve_ora_gene_input(signed_rank = signed_rank, gene_fraction = 0.4)

  expect_equal(resolved$mode, "ranked")
  expect_setequal(resolved$selected$selected, c("A", "B"))
  expect_equal(resolved$selected$up, "A")
  expect_equal(resolved$selected$down, "B")
  expect_equal(resolved$universe, names(signed_rank))
})

test_that("ORA input resolution supports explicitly selected genes", {
  selected_genes <- c(A = 2, B = -1, C = 0)

  resolved <- resolve_ora_gene_input(
    selected_genes = selected_genes,
    gene_fraction = 1,
    universe = c("A", "B", "C", "D")
  )

  expect_equal(resolved$mode, "selected")
  expect_setequal(resolved$selected$selected, c("A", "B", "C"))
  expect_equal(resolved$selected$up, "A")
  expect_equal(resolved$selected$down, "B")
  expect_equal(resolved$universe, c("A", "B", "C", "D"))
})

test_that("ORA input resolution enforces compatible mode arguments", {
  signed_values <- c(A = 2, B = -1)

  expect_error(
    resolve_ora_gene_input(),
    "Provide exactly one of signed_rank or selected_genes"
  )
  expect_error(
    resolve_ora_gene_input(signed_rank = signed_values, selected_genes = signed_values),
    "Provide exactly one of signed_rank or selected_genes"
  )
  expect_error(
    resolve_ora_gene_input(selected_genes = signed_values),
    "universe must be supplied when selected_genes is used"
  )
  expect_error(
    resolve_ora_gene_input(
      selected_genes = signed_values,
      gene_fraction = 0.5,
      universe = c("A", "B", "C")
    ),
    "gene_fraction must be 1 when selected_genes is supplied"
  )
  expect_error(
    resolve_ora_gene_input(
      selected_genes = signed_values,
      universe = c("A", "C")
    ),
    "universe must include all selected genes"
  )
})

test_that("rank validation removes duplicate and unusable names", {
  signed_rank <- setNames(c(1, 2, 3, NA_real_, -1), c("101", "101", "", "102", "103"))

  validated <- validate_signed_rank(signed_rank)

  expect_equal(names(validated), c("101", "103"))
  expect_equal(unname(validated), c(1, -1))
})

test_that("GSEA rankings keep sign for signed sets and magnitude for unsigned sets", {
  signed_rank <- c(A = -2, B = 5, C = 1)

  rankings <- prepare_gsea_rankings(signed_rank)

  expect_equal(names(rankings$signed), c("B", "C", "A"))
  expect_equal(names(rankings$unsigned_magnitude), c("B", "A", "C"))
  expect_equal(unname(rankings$unsigned_magnitude), c(5, 2, 1))
})
