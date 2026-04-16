\name{download_kegg_gmt}
\alias{download_kegg_gmt}
\title{Download KEGG Pathway GMT Files}
\description{
Retrieve KEGG pathway gene memberships with \pkg{KEGGREST} and write them as a
GMT file.
}
\usage{
download_kegg_gmt(
  organism = "hsa",
  output_dir = "."
)
}
\arguments{
  \item{organism}{KEGG organism code such as \code{"hsa"}, \code{"mmu"}, or
  \code{"rno"}.}
  \item{output_dir}{Directory where the GMT file is written.}
}
\value{
The written GMT filepath, returned invisibly.
}
