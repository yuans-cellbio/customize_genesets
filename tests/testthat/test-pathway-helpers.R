test_that("GMT helpers round-trip names and clean duplicate genes", {
  tmp_gmt <- tempfile(fileext = ".gmt")

  write_gmt_from_list(
    gene_sets = list(SetA = c("1", "2", "2", NA_character_)),
    descriptions = c(SetA = "Example set"),
    filepath = tmp_gmt
  )

  term2name <- read_gmt_term2name(tmp_gmt)
  expect_equal(term2name$ID, "SetA")
  expect_equal(term2name$Name, "Example set")

  parts <- strsplit(readLines(tmp_gmt), "\t", fixed = TRUE)[[1]]
  expect_equal(parts, c("SetA", "Example set", "1", "2"))
})

test_that("symbol mapping keeps the first successful keytype", {
  skip_if_not_installed("org.Hs.eg.db")
  library(org.Hs.eg.db)

  mapping <- map_symbols_to_entrez(org.Hs.eg.db, c("TP53", "P53"))

  expect_false(is.na(mapping[["TP53"]]))
  expect_false(is.na(mapping[["P53"]]))
})
