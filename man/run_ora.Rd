\name{run_ora}
\alias{run_ora}
\title{Run Over-Representation Analysis on Signed and Unsigned Libraries}
\description{
Run ORA against signed and unsigned pathway libraries using either a full
ranked vector or an explicitly selected gene foreground.
}
\usage{
run_ora(
    signed_rank = NULL,
    selected_genes = NULL,
    gene_fraction = NULL,
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
  \item{signed_rank}{Optional named numeric vector of gene-level signed
  statistics using Entrez Gene IDs as names. Supply this to use ranked-selection
  mode.}
  \item{selected_genes}{Optional named numeric vector of selected gene-level
  signed statistics using Entrez Gene IDs as names. Supply this to use
  explicit-hit mode, where the provided genes are used directly as the ORA
  foreground.}
  \item{gene_fraction}{Fraction of genes to carry into ORA based on the
  absolute statistic when \code{signed_rank} is supplied. Leave \code{NULL} or
  set it to \code{1} when using \code{selected_genes}.}
  \item{org_db}{An \code{OrgDb} object used to make the enrichment results
  human-readable.}
  \item{universe}{Background gene universe. Defaults to the names of
  \code{signed_rank} in ranked-selection mode and must be supplied in
  explicit-hit mode.}
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
In ranked-selection mode, signed libraries are tested separately against the up-
and down-ranked subsets chosen by \code{gene_fraction}, while unsigned
libraries are tested against their union.

In explicit-hit mode, the supplied \code{selected_genes} are used directly as
the ORA foreground. Their numeric values are only used to split signed
libraries into up- and down-selected subsets; unsigned libraries use gene
membership only.
}
