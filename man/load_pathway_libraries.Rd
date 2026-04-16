\name{load_pathway_libraries}
\alias{load_pathway_libraries}
\title{Load Signed and Unsigned GMT Libraries}
\description{
Read GMT files from a directory and split them into signed and unsigned
libraries suitable for \code{run_ora()} and \code{run_gsea()}.
}
\usage{
load_pathway_libraries(pathway_dir = "pathway_library")
}
\arguments{
  \item{pathway_dir}{Directory containing GMT files. Files whose names contain
  \code{"_signed"} are treated as signed libraries; the rest are treated as
  unsigned libraries.}
}
\value{
A list with \code{signed} and \code{unsigned} elements. Each element is a list
of library entries containing \code{lib_name}, \code{term2gene}, and
  \code{term2name}.
}
