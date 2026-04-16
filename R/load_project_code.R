load_project_code <- function(project_root = NULL) {
  if (is.null(project_root)) {
    loader_path <- NULL
    for (frame in rev(sys.frames())) {
      if (!is.null(frame$ofile)) {
        loader_path <- frame$ofile
        break
      }
    }

    project_root <- if (is.null(loader_path)) {
      normalizePath(".", winslash = "/", mustWork = TRUE)
    } else {
      dirname(dirname(normalizePath(loader_path, winslash = "/", mustWork = TRUE)))
    }
  }

  source(file.path(project_root, "R", "pathway_helpers.R"), local = FALSE)
  source_project_modules(project_root)
  invisible(project_root)
}
