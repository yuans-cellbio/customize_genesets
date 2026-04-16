\name{build_pathway_library}
\alias{build_pathway_library}
\title{Build a Full Pathway Library}
\description{
Build a species-aware pathway gene set library and write GMT files plus a
provenance README to an output directory.
}
\usage{
build_pathway_library(
    org_db = NULL,
    species = "human",
    include_IEA = TRUE,
    pos_keywords = NULL,
    neg_keywords = NULL,
    pair_replacements = DEFAULT_DIRECTIONAL_PAIR_REPLACEMENTS,
    min_size = 15,
    max_size = 500,
    dag_overlap_threshold = 0.8,
    jaccard_threshold = 0.7,
    msigdb_collections = c("H", "C2"),
    msigdb_subcollections = NULL,
    progeny_top_n = 100,
    collectri_min_targets = 10,
    collectri_max_targets = 500,
    skip = character(0),
    go_script = "build_go_regulation_gmt.R",
    pathway_script = "download_pathway_gmt.R",
    output_dir = "pathway_library",
    ...
)
}
\arguments{
  \item{org_db}{Optional \code{OrgDb} object or package name. Leave as
  \code{NULL} to infer the matching \code{OrgDb} from \code{species}.}
  \item{species}{Species selector. Supported values are human, mouse, and rat,
  along with their canonical Latin names.}
  \item{include_IEA}{Logical; include electronically inferred GO annotations.}
  \item{pos_keywords}{Optional character vector of positive-direction GO BP
  keywords. Defaults to the package keyword list.}
  \item{neg_keywords}{Optional character vector of negative-direction GO BP
  keywords. Defaults to the package keyword list.}
  \item{pair_replacements}{Named list of explicit directional text
  replacements used when preserving positive/negative sibling GO terms.}
  \item{min_size}{Minimum allowed gene-set size for GO and regulator-derived
  sets.}
  \item{max_size}{Maximum allowed gene-set size for GO and regulator-derived
  sets.}
  \item{dag_overlap_threshold}{Ancestor-overlap threshold used for GO DAG
  pruning.}
  \item{jaccard_threshold}{Jaccard-overlap threshold used for GO
  deduplication.}
  \item{msigdb_collections}{Character vector of MSigDB collections to export.}
  \item{msigdb_subcollections}{Optional named list mapping collection codes to
  subcollections.}
  \item{progeny_top_n}{Number of lowest-\code{p.value} genes retained per
  PROGENy pathway.}
  \item{collectri_min_targets}{Minimum number of mapped targets required for a
  CollecTRI TF gene set.}
  \item{collectri_max_targets}{Maximum number of mapped targets allowed for a
  CollecTRI TF gene set.}
  \item{skip}{Character vector of build stages to skip. Supported stage names
  include \code{"go_bp"}, \code{"go_cc"}, \code{"go_mf"}, \code{"msigdb"},
  \code{"kegg"}, \code{"reactome"}, \code{"progeny"}, and
  \code{"collectri"}.}
  \item{go_script}{Compatibility argument for source-tree usage. Usually leave
  unchanged.}
  \item{pathway_script}{Compatibility argument for source-tree usage. Usually
  leave unchanged.}
  \item{output_dir}{Directory where GMT files and the build README are
  written.}
  \item{...}{Reserved for future use. Per-database species overrides are not
  accepted by this API.}
}
\value{
An invisible named list containing per-source results, a file inventory, a
version summary, and the path to the generated README.
}
\details{
This is the recommended entry point for most users. The function canonicalizes
\code{species}, infers the matching \code{OrgDb} when needed, applies
species-specific defaults for KEGG, MSigDB, PROGENy, and CollecTRI, and then
runs the requested build stages in sequence.
}
\examples{
\dontrun{
build_pathway_library(
  species = "human",
  msigdb_collections = c("H", "C2"),
  skip = c("kegg", "collectri"),
  output_dir = "pathway_library"
)
}
}
