\name{build_go_regulation_library}
\alias{build_go_regulation_library}
\title{Build Signed and Unsigned GO Biological Process Libraries}
\description{
Classify GO Biological Process terms into positive, negative, and unsigned
groups using explicit directional keywords, then export signed and unsigned GMT
libraries after size filtering, DAG pruning, and Jaccard deduplication.
}
\usage{
build_go_regulation_library(
  org_db,
  include_IEA = TRUE,
  pos_keywords = DEFAULT_POS_KEYWORDS,
  neg_keywords = DEFAULT_NEG_KEYWORDS,
  pair_replacements = DEFAULT_DIRECTIONAL_PAIR_REPLACEMENTS,
  min_size = 15,
  max_size = 500,
  dag_overlap_threshold = 0.8,
  jaccard_threshold = 0.7,
  output_dir = "."
)
}
\arguments{
  \item{org_db}{An \code{OrgDb} object for the target species.}
  \item{include_IEA}{Logical; include electronically inferred GO annotations.}
  \item{pos_keywords}{Character vector of positive-direction GO BP keywords.}
  \item{neg_keywords}{Character vector of negative-direction GO BP keywords.}
  \item{pair_replacements}{Named list of explicit text replacements used to
  identify positive/negative sibling GO terms.}
  \item{min_size}{Minimum allowed gene-set size.}
  \item{max_size}{Maximum allowed gene-set size.}
  \item{dag_overlap_threshold}{Ancestor-overlap threshold for pruning GO DAG
  ancestors.}
  \item{jaccard_threshold}{Jaccard-overlap threshold for deduplicating GO
  terms.}
  \item{output_dir}{Directory where the signed and unsigned GMT files are
  written.}
}
\value{
An invisible list with the signed and unsigned gene sets, the GO term table,
the directional pair map, and the paths to the written GMT files.
}
\examples{
\dontrun{
library(org.Hs.eg.db)
build_go_regulation_library(org.Hs.eg.db, output_dir = "go_dedup_gmt")
}
}
