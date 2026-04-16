\name{download_progeny_gmt}
\alias{download_progeny_gmt}
\title{Export PROGENy Gene Sets as GMT Files}
\description{
Export signed and unsigned PROGENy pathway footprints as Entrez-based GMT
files.
}
\usage{
download_progeny_gmt(
    org_db,
    organism = "human",
    top_n = 100,
    output_dir = "."
)
}
\arguments{
  \item{org_db}{An \code{OrgDb} object for the target species.}
  \item{organism}{PROGENy model selector. Supported values are \code{"human"}
  and \code{"mouse"}.}
  \item{top_n}{Number of lowest-\code{p.value} genes retained per PROGENy
  pathway.}
  \item{output_dir}{Directory where the GMT files are written.}
}
\value{
An invisible list with \code{signed_gmt} and \code{unsigned_gmt} paths.
}
