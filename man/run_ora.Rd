\name{run_ora}
\alias{run_ora}
\title{Run Over-Representation Analysis on Signed and Unsigned Libraries}
\description{
Run ORA against signed and unsigned pathway libraries using a named numeric
vector of signed gene-level statistics.
}
\usage{
run_ora(
    signed_rank,
    gene_fraction = NULL,
    top_frac_threshold = NULL,
    org_db,
    universe = NULL,
    signed_library,
    unsigned_library,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 15,
    maxGSSize = 500
)
}
\arguments{
  \item{signed_rank}{Named numeric vector of gene-level signed statistics using
  Entrez Gene IDs as names.}
  \item{gene_fraction}{Fraction of genes to carry into ORA based on the
  absolute statistic.}
  \item{top_frac_threshold}{Backward-compatible alias for
  \code{gene_fraction}.}
  \item{org_db}{An \code{OrgDb} object used to make the enrichment results
  human-readable.}
  \item{universe}{Optional background gene universe. Defaults to the names of
  \code{signed_rank}.}
  \item{signed_library}{List of signed library entries from
  \code{load_pathway_libraries()}.}
  \item{unsigned_library}{List of unsigned library entries from
  \code{load_pathway_libraries()}.}
  \item{pvalueCutoff}{Raw \code{p.value} cutoff passed to
  \code{clusterProfiler::enricher()}.}
  \item{qvalueCutoff}{\code{q.value} cutoff passed to
  \code{clusterProfiler::enricher()}.}
  \item{pAdjustMethod}{Multiple-testing correction method applied after
  combining results within each library-direction group.}
  \item{minGSSize}{Minimum gene-set size passed to
  \code{clusterProfiler::enricher()}.}
  \item{maxGSSize}{Maximum gene-set size passed to
  \code{clusterProfiler::enricher()}.}
}
\value{
A data frame of combined ORA results, or \code{NULL} when no enrichments are
returned.
}
\details{
Signed libraries are tested separately against the up- and down-ranked subsets.
Unsigned libraries are tested against the union of genes selected by absolute
rank magnitude.
}
