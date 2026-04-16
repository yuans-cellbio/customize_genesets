\name{download_msigdb_gmt}
\alias{download_msigdb_gmt}
\title{Export MSigDB Collections as GMT Files}
\description{
Export one or more MSigDB collections as GMT files using \pkg{msigdbr}. Mouse
builds use the mouse-native MSigDB database by default; rat builds use
ortholog-mapped human MSigDB because a rat-native MSigDB database is not
available in \pkg{msigdbr}.
}
\usage{
download_msigdb_gmt(
  species = "Homo sapiens",
  db_species = NULL,
  collections = "all",
  subcollections = NULL,
  output_dir = "."
)
}
\arguments{
  \item{species}{Target species name.}
  \item{db_species}{Optional MSigDB database selector such as \code{"HS"} or
  \code{"MM"}. Leave as \code{NULL} to use the package defaults.}
  \item{collections}{Collection code or vector of collection codes. Use
  \code{"all"} to export every collection available for the chosen database.}
  \item{subcollections}{Optional named list mapping collection codes to
  subcollections.}
  \item{output_dir}{Directory where GMT files are written.}
}
\value{
An invisible named list mapping requested collection codes to written GMT file
paths.
}
