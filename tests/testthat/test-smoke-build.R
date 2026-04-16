test_that("offline build smoke test completes", {
  skip_if_not_installed("GO.db")
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("reactome.db")
  skip_if_not_installed("progeny")

  library(org.Hs.eg.db)

  out_dir <- tempfile("pathway_smoke_")
  result <- build_pathway_library(
    org_db = org.Hs.eg.db,
    species = "Homo sapiens",
    msigdb_collections = c("H"),
    skip = c("kegg", "collectri"),
    output_dir = out_dir
  )

  expect_true(file.exists(result$go_bp$signed_gmt))
  expect_true(file.exists(result$go_cc$gmt_file))
  expect_true(file.exists(result$go_mf$gmt_file))
  expect_true(file.exists(result$reactome))
  expect_true(file.exists(result$progeny$signed_gmt))
  expect_true(file.exists(result$readme))
})
