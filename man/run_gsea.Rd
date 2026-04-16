\name{run_gsea}
\alias{run_gsea}
\title{Run GSEA on Signed and Unsigned Libraries}
\description{
Run GSEA against signed and unsigned pathway libraries using a named numeric
vector of signed gene-level statistics.
}
\usage{
run_gsea(
    signed_rank,
    org_db,
    signed_library,
    unsigned_library,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 15,
    maxGSSize = 500,
    nPermSimple = 10000
)
}
\arguments{
  \item{signed_rank}{Named numeric vector of gene-level signed statistics using
  Entrez Gene IDs as names.}
  \item{org_db}{An \code{OrgDb} object used to make the enrichment results
  human-readable.}
  \item{signed_library}{List of signed library entries from
  \code{load_pathway_libraries()}.}
  \item{unsigned_library}{List of unsigned library entries from
  \code{load_pathway_libraries()}.}
  \item{pvalueCutoff}{Raw \code{p.value} cutoff passed to
  \code{clusterProfiler::GSEA()}.}
  \item{pAdjustMethod}{Multiple-testing correction method applied after
  combining results within each library-direction group.}
  \item{minGSSize}{Minimum gene-set size passed to \code{clusterProfiler::GSEA()}.}
  \item{maxGSSize}{Maximum gene-set size passed to \code{clusterProfiler::GSEA()}.}
  \item{nPermSimple}{Number of simple permutations passed to
  \code{clusterProfiler::GSEA()}.}
}
\value{
A data frame of combined GSEA results, or \code{NULL} when no enrichments are
returned.
}
\details{
Signed libraries use the signed statistic, while unsigned libraries use the
absolute statistic and a positive-only score type.
}
