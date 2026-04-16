\name{customizeGeneSets-package}
\alias{customizeGeneSets-package}
\alias{customizeGeneSets}
\docType{package}
\title{Build Pathway Gene Set Libraries for Enrichment Analysis}
\description{
Build reproducible pathway gene set libraries from GO, MSigDB, KEGG,
Reactome, PROGENy, and CollecTRI, then use those libraries in ORA or GSEA
workflows with consistent signed versus unsigned conventions.
}
\details{
The package is organized around three tasks:
\enumerate{
  \item Build complete or source-specific GMT libraries.
  \item Load GMT libraries into \code{TERM2GENE}/\code{TERM2NAME} structures.
  \item Run ORA or GSEA on signed and unsigned libraries in a consistent way.
}

The top-level \code{build_pathway_library()} API accepts \code{species} as the
public organism switch and infers the matching \code{OrgDb} package plus the
database-specific organism keys internally.
}
\seealso{
\code{\link{build_pathway_library}},
\code{\link{download_all_pathway_gmt}},
\code{\link{load_pathway_libraries}},
\code{\link{run_ora}},
\code{\link{run_gsea}}
}
\keyword{package}
