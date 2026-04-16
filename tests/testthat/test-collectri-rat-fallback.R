test_that("direct rat CollecTRI normalization collapses complexes and signs edges", {
  raw_collectri <- data.frame(
    source = c("STAT1", "COMPLEX_AP1", "COMPLEX_NFKB"),
    source_genesymbol = c("STAT1", "JUN_FOS", "REL_NFKB1"),
    target_genesymbol = c("IL6", "TNF", "CCL2"),
    is_stimulation = c(1, 1, 0),
    is_inhibition = c(0, 0, 1),
    stringsAsFactors = FALSE
  )

  normalized <- normalize_collectri_direct_query(raw_collectri, split_complexes = FALSE)

  expect_equal(colnames(normalized), c("source", "target", "mor"))
  expect_equal(normalized$source, c("Stat1", "Ap1", "Nfkb"))
  expect_equal(normalized$target, c("Il6", "Tnf", "Ccl2"))
  expect_equal(normalized$mor, c(1, 1, -1))
})
