\name{build_go_dedup_library}
\alias{build_go_dedup_library}
\title{Build a Deduplicated GO Ontology Library}
\description{
Build a deduplicated GO library for one ontology (\code{BP}, \code{CC}, or
\code{MF}) after size filtering, DAG pruning, and Jaccard deduplication.
}
\usage{
build_go_dedup_library(
  org_db,
  ontology = "CC",
  include_IEA = TRUE,
  min_size = 15,
  max_size = 500,
  dag_overlap_threshold = 0.8,
  jaccard_threshold = 0.7,
  output_dir = "."
)
}
\arguments{
  \item{org_db}{An \code{OrgDb} object for the target species.}
  \item{ontology}{GO ontology code: \code{"BP"}, \code{"CC"}, or
  \code{"MF"}.}
  \item{include_IEA}{Logical; include electronically inferred GO annotations.}
  \item{min_size}{Minimum allowed gene-set size.}
  \item{max_size}{Maximum allowed gene-set size.}
  \item{dag_overlap_threshold}{Ancestor-overlap threshold for DAG pruning.}
  \item{jaccard_threshold}{Jaccard-overlap threshold for deduplication.}
  \item{output_dir}{Directory where the GMT file is written.}
}
\value{
An invisible list with the final gene sets, the GO term table, and the path to
the written GMT file.
}
\examples{
\dontrun{
library(org.Hs.eg.db)
build_go_dedup_library(org.Hs.eg.db, ontology = "CC", output_dir = "go_dedup_gmt")
}
}
