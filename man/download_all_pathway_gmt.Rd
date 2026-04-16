\name{download_all_pathway_gmt}
\alias{download_all_pathway_gmt}
\title{Download or Export Non-GO Pathway GMT Libraries}
\description{
Run the non-GO pathway exporters for MSigDB, KEGG, Reactome, PROGENy, and
CollecTRI using the species-aware defaults of the package.
}
\usage{
download_all_pathway_gmt(
    org_db = NULL,
    species = "human",
    collectri_min_targets = 10,
    collectri_max_targets = 500,
    msigdb_collections = c("H", "C2"),
    output_dir = "pathway_gmt",
    ...
)
}
\arguments{
  \item{org_db}{Optional \code{OrgDb} object or package name. Leave as
  \code{NULL} to infer the matching \code{OrgDb} from \code{species}.}
  \item{species}{Species selector. Supported values are human, mouse, and rat,
  along with their canonical Latin names.}
  \item{collectri_min_targets}{Minimum number of mapped targets required for a
  CollecTRI TF gene set.}
  \item{collectri_max_targets}{Maximum number of mapped targets allowed for a
  CollecTRI TF gene set.}
  \item{msigdb_collections}{Character vector of MSigDB collections to export.}
  \item{output_dir}{Directory where GMT files are written.}
  \item{...}{Reserved for future use. Per-database species overrides are not
  accepted by this API.}
}
\value{
An invisible list containing the output paths returned by each completed
exporter.
}
\details{
Rat builds omit PROGENy automatically because this package does not provide a
rat PROGENy model.
}
