library(testthat)

if (file.exists(file.path("R", "load_project_code.R"))) {
  source(file.path("R", "load_project_code.R"), local = FALSE)
  load_project_code()
  test_dir(file.path("tests", "testthat"))
} else {
  library(customizeGeneSets)
  test_check("customizeGeneSets")
}
