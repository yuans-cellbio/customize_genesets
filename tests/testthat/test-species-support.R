test_that("species helpers canonicalize supported labels", {
  expect_equal(canonical_species_name("human"), "Homo sapiens")
  expect_equal(canonical_species_name("mouse"), "Mus musculus")
  expect_equal(canonical_species_name("mosue"), "Mus musculus")
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

  expect_message(
    expect_equal(resolve_msigdb_db_species("Rattus norvegicus"), "HS"),
    "ortholog-mapped human MSigDB"
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

test_that("species build context infers matching OrgDb objects and internal keys", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("org.Rn.eg.db")

  human_cfg <- resolve_species_build_context("human")
  expect_equal(human_cfg$species, "Homo sapiens")
  expect_equal(human_cfg$org_db_package, "org.Hs.eg.db")
  expect_equal(organism_name_from_orgdb(human_cfg$org_db), "Homo sapiens")
  expect_equal(human_cfg$kegg_organism, "hsa")
  expect_equal(human_cfg$progeny_organism, "human")
  expect_equal(human_cfg$collectri_organism, "human")
  expect_false("progeny" %in% human_cfg$skip)

  mouse_cfg <- resolve_species_build_context("mouse")
  expect_equal(mouse_cfg$species, "Mus musculus")
  expect_equal(mouse_cfg$org_db_package, "org.Mm.eg.db")
  expect_equal(organism_name_from_orgdb(mouse_cfg$org_db), "Mus musculus")
  expect_equal(mouse_cfg$msigdb_db_species, "MM")
  expect_equal(mouse_cfg$kegg_organism, "mmu")
  expect_equal(mouse_cfg$progeny_organism, "mouse")
  expect_equal(mouse_cfg$collectri_organism, "mouse")

  expect_message(
    rat_cfg <- resolve_species_build_context("rat"),
    "Rat builds omit PROGENy"
  )
  expect_equal(rat_cfg$species, "Rattus norvegicus")
  expect_equal(rat_cfg$org_db_package, "org.Rn.eg.db")
  expect_equal(organism_name_from_orgdb(rat_cfg$org_db), "Rattus norvegicus")
  expect_equal(rat_cfg$msigdb_db_species, "HS")
  expect_equal(rat_cfg$kegg_organism, "rno")
  expect_null(rat_cfg$progeny_organism)
  expect_equal(rat_cfg$collectri_organism, "rat")
  expect_true("progeny" %in% rat_cfg$skip)
})

test_that("species build context rejects mismatched OrgDb inputs", {
  skip_if_not_installed("org.Hs.eg.db")

  library(org.Hs.eg.db)

  expect_error(
    resolve_species_build_context("mouse", org_db = org.Hs.eg.db),
    "supplied org_db is built for"
  )
})

test_that("build_pathway_library rejects retired per-database species overrides", {
  expect_error(
    do.call(build_pathway_library, list(species = "human", kegg_organism = "hsa")),
    "no longer accepts per-database species overrides"
  )
})
