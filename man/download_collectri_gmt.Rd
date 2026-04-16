\name{download_collectri_gmt}
\alias{download_collectri_gmt}
\title{Export CollecTRI Regulator Gene Sets as GMT Files}
\description{
Export signed and unsigned CollecTRI transcription factor target sets as
Entrez-based GMT files. Human and mouse builds use \pkg{decoupleR}; rat builds
use the direct \pkg{OmnipathR} CollecTRI interactions route and then normalize
the interaction table into the same output format.
}
\usage{
download_collectri_gmt(
    org_db,
    organism = "human",
    split_complexes = FALSE,
    min_targets = 10,
    max_targets = 500,
    output_dir = "."
)
}
\arguments{
  \item{org_db}{An \code{OrgDb} object for the target species.}
  \item{organism}{Organism selector. Supported values are human, mouse, and
  rat.}
  \item{split_complexes}{Logical; keep CollecTRI complexes split into their
  reported components instead of collapsing AP1 and NFKB-style complexes.}
  \item{min_targets}{Minimum number of mapped targets required for a TF gene
  set.}
  \item{max_targets}{Maximum number of mapped targets allowed for a TF gene
  set.}
  \item{output_dir}{Directory where the GMT files are written.}
}
\value{
An invisible list with \code{signed_gmt} and \code{unsigned_gmt} paths.
}
