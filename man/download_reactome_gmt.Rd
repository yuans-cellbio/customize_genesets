\name{download_reactome_gmt}
\alias{download_reactome_gmt}
\title{Export Reactome Pathways as GMT Files}
\description{
Export Reactome pathways either from \pkg{reactome.db} or from the Reactome
subset of MSigDB.
}
\usage{
download_reactome_gmt(
  species = "Homo sapiens",
  method = "reactome.db",
  output_dir = "."
)
}
\arguments{
  \item{species}{Target species name.}
  \item{method}{Reactome source. Supported values are \code{"reactome.db"} and
  \code{"msigdbr"}.}
  \item{output_dir}{Directory where the GMT file is written.}
}
\value{
The written GMT filepath, returned invisibly.
}
