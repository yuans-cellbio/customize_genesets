test_that("directional pair map uses explicit text replacements", {
  term_df <- data.frame(
    GOID = c("GO:1", "GO:2", "GO:3"),
    TERM = c(
      "positive regulation of foo",
      "negative regulation of foo",
      "promotion of foo"
    ),
    stringsAsFactors = FALSE
  )

  pair_map <- build_directional_pair_map(term_df)

  expect_equal(pair_map[["GO:1"]], "GO:2")
  expect_equal(pair_map[["GO:2"]], "GO:1")
  expect_true(is.na(pair_map[["GO:3"]]))
})

test_that("jaccard dedup can prefer the more specific representative", {
  gene_sets <- list(
    GO_A = c("1", "2", "3", "4"),
    GO_B = c("1", "2", "3", "4"),
    GO_C = c("7", "8")
  )

  deduped <- dedup_by_jaccard(
    gene_sets,
    jaccard_threshold = 0.9,
    representative_scores = c(GO_A = 1, GO_B = 5, GO_C = 1)
  )

  expect_setequal(names(deduped), c("GO_B", "GO_C"))
})

test_that("paired survival restores missing siblings after deduplication", {
  surviving <- list("GO:1" = c("1", "2"))
  all_sets <- list(
    "GO:1" = c("1", "2"),
    "GO:2" = c("3", "4")
  )
  pair_map <- list("GO:1" = "GO:2", "GO:2" = "GO:1")

  restored <- enforce_paired_survival(surviving, all_sets, pair_map)

  expect_setequal(names(restored), c("GO:1", "GO:2"))
})
