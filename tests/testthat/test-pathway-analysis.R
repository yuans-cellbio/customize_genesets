test_that("ORA selection uses an exact top-n fraction", {
  signed_rank <- c(A = 5, B = -4, C = 3, D = -2, E = 1)

  selected <- select_ranked_genes_for_ora(signed_rank, gene_fraction = 0.4)

  expect_setequal(selected$selected, c("A", "B"))
  expect_equal(selected$up, "A")
  expect_equal(selected$down, "B")
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
