test_that("species helpers canonicalize common names", {
  expect_equal(canonical_species_name("mouse"), "Mus musculus")
  expect_equal(canonical_species_name("rat"), "Rattus norvegicus")
  expect_true(is_mouse_species("Mus musculus"))
  expect_true(is_rat_species("rat"))
  expect_true(is_human_species("human"))
})

test_that("MSigDB database selection uses mouse-native sets and rat ortholog mapping", {
  expect_message(
    expect_equal(resolve_msigdb_db_species("Mus musculus"), "MM"),
    "mouse-native MSigDB"
  )

  expect_warning(
    expect_equal(resolve_msigdb_db_species("Rattus norvegicus"), "HS"),
    "ortholog mapping"
  )

  expect_equal(resolve_msigdb_db_species("Homo sapiens"), "HS")
  expect_equal(resolve_msigdb_db_species("Mus musculus", db_species = "HS"), "HS")
})

test_that("mouse-native MSigDB collection codes are translated correctly", {
  expect_equal(translate_msigdb_collection_code("H", db_species = "MM"), "MH")
  expect_equal(translate_msigdb_collection_code("C2", db_species = "MM"), "M2")
  expect_equal(translate_msigdb_collection_code("C7", db_species = "MM"), "M7")
  expect_equal(translate_msigdb_collection_code("H", db_species = "HS"), "H")
})

test_that("builder harmonization aligns mouse settings and guards rat PROGENy", {
  mouse_cfg <- suppressWarnings(harmonize_species_build_config(
    species = "mouse",
    kegg_organism = "hsa",
    progeny_organism = "human",
    collectri_organism = "human",
    skip = character(0)
  ))

  expect_equal(mouse_cfg$species, "Mus musculus")
  expect_equal(mouse_cfg$kegg_organism, "mmu")
  expect_equal(mouse_cfg$progeny_organism, "mouse")
  expect_equal(mouse_cfg$collectri_organism, "mouse")
  expect_false("progeny" %in% mouse_cfg$skip)

  expect_warning(
    rat_cfg <- harmonize_species_build_config(
      species = "Rattus norvegicus",
      kegg_organism = "hsa",
      progeny_organism = "human",
      collectri_organism = "human",
      skip = character(0)
    ),
    "Rat build requested"
  )

  expect_equal(rat_cfg$species, "Rattus norvegicus")
  expect_equal(rat_cfg$kegg_organism, "rno")
  expect_equal(rat_cfg$collectri_organism, "rat")
  expect_true("progeny" %in% rat_cfg$skip)
})
