# Customize Gene Sets

Build reproducible pathway GMT libraries for enrichment analysis, with a focus on defensible GO processing, explicit provenance, and code that is easier to troubleshoot than a single monolithic script.

This repo does two related jobs:

1. It builds pathway libraries from GO, MSigDB, KEGG, Reactome, PROGENy, and CollecTRI.
2. It provides helper functions for running ORA and GSEA against those libraries in a consistent way.

The generated human library lives in [`pathway_library/README.md`](./pathway_library/README.md), which records the exact build date, package versions, parameters, and file inventory for the checked-in artifacts.

## Why This Repo Exists

Many pathway workflows mix several pain points:

- GMT files come from different sources and use different identifier systems.
- GO collections are often highly redundant and hard to interpret.
- Signed and unsigned libraries are easily conflated in downstream analysis.
- Older scripts tend to accumulate hidden assumptions that are hard to audit.

This project tries to make those assumptions visible and testable.

## Core Methodology

- All exported GMT files use NCBI Entrez Gene IDs for consistency across sources.
- GO Biological Process signed terms are classified from explicit positive and negative keyword lists.
- Signed GO sibling retention is conservative: only explicit text replacements such as positive versus negative regulation are used when keeping paired terms together.
- GO libraries are filtered by gene-set size, pruned using GO DAG ancestor overlap, and then deduplicated by Jaccard overlap.
- When Jaccard-deduplicating GO sets, the pipeline prefers the more specific term rather than simply the largest one.
- Mouse builds use native mouse MSigDB collections through `msigdbr` with `db_species = "MM"`.
- Rat builds use human MSigDB collections with ortholog mapping because `msigdbr` does not provide a rat-native MSigDB database.
- External resources are exported into the same GMT-style interface so they can be used side by side in `clusterProfiler` workflows.

## Species Support

- Human: fully supported by the checked-in build and all main workflows.
- Mouse: supported for GO, MSigDB, Reactome, KEGG, PROGENy, and CollecTRI.
- Rat: supported for GO, MSigDB, Reactome, and KEGG.
- Rat CollecTRI: the pipeline now tries a rat-specific OmnipathR fallback based on direct CollecTRI interactions, but this path is still fragile and depends on OmniPath/Ensembl behavior.
- Rat PROGENy: not supported in this pipeline. The top-level builder now warns and skips that step automatically.

## Common Use Cases

- Build a defensible human pathway library for transcriptomics or proteomics enrichment analysis.
- Replace inherited GMT files with a versioned build that records package versions and retrieval dates.
- Generate a signed GO BP library for analyses where directionality matters.
- Build mouse-native MSigDB pathway libraries without falling back to human ortholog mapping.
- Build a rat pathway library while preserving explicit warnings about unsupported steps.
- Combine broad pathway databases like Reactome and KEGG with footprint-style resources like PROGENy and CollecTRI.
- Customize the GO keyword lists or pathway-source parameters for a specific biological question.
- Troubleshoot one data source at a time instead of rerunning the full pipeline.

## Project Layout

- `R/pathway_helpers.R`: shared validation, GMT IO, ID mapping, species helpers, and module-loading helpers.
- `R/go_library.R`: GO keyword classification, sibling pairing, DAG pruning, Jaccard deduplication, and GO GMT builders.
- `R/pathway_downloaders.R`: MSigDB, KEGG, Reactome, PROGENy, and CollecTRI download/export functions.
- `R/pathway_library_builder.R`: top-level orchestration, species-aware defaults, and provenance README generation.
- `R/pathway_analysis.R`: helpers for loading libraries and running ORA or GSEA consistently.
- `build_go_regulation_gmt.R`: lightweight wrapper for GO-focused builds.
- `download_pathway_gmt.R`: lightweight wrapper for downloading/exporting non-GO sources.
- `build_pathway_library.R`: lightweight wrapper for the default full build.
- `build_human_pathway_library.R`: checked-in human build entry point used to regenerate `pathway_library/`.
- `pathway_analysis.Rmd`: example downstream analysis notebook using the shared helpers.
- `tests/testthat/`: unit tests plus a smoke test for the offline-capable build path.

## Typical Workflows

### 1. Use the checked-in human library

If you just want the current human GMT files, start in `pathway_library/`.

```r
library(clusterProfiler)
source("R/load_project_code.R")

gmt <- read.gmt("pathway_library/reactome_homo_sapiens.gmt")
gmt_names <- read_gmt_term2name("pathway_library/reactome_homo_sapiens.gmt")
```

### 2. Rebuild the default human library

```r
source("R/load_project_code.R")
library(org.Hs.eg.db)

build_pathway_library(
  org_db = org.Hs.eg.db,
  species = "Homo sapiens",
  output_dir = "pathway_library"
)
```

### 3. Rebuild with custom parameters

Use this when you want to change GO filtering thresholds, select different MSigDB collections, or skip network-backed resources during development.

```r
source("R/load_project_code.R")
library(org.Hs.eg.db)

build_pathway_library(
  org_db = org.Hs.eg.db,
  species = "Homo sapiens",
  include_IEA = FALSE,
  min_size = 20,
  max_size = 300,
  msigdb_collections = c("H", "C2"),
  skip = c("kegg", "collectri"),
  output_dir = "pathway_library_custom"
)
```

### 4. Build a mouse library

The top-level builder now aligns the species-specific settings automatically for mouse builds:

- `species` is canonicalized to `Mus musculus`
- `msigdb_db_species` defaults to `MM`
- requested human-style collection codes such as `H`, `C2`, and `C7` are translated to mouse-native codes such as `MH`, `M2`, and `M7`
- `kegg_organism` is set to `mmu`
- `progeny_organism` is set to `mouse`
- `collectri_organism` is set to `mouse`

```r
source("R/load_project_code.R")
library(org.Mm.eg.db)

build_pathway_library(
  org_db = org.Mm.eg.db,
  species = "Mus musculus",
  output_dir = "pathway_library_mouse"
)
```

If you want to be completely explicit:

```r
source("R/load_project_code.R")
library(org.Mm.eg.db)

build_pathway_library(
  org_db = org.Mm.eg.db,
  species = "Mus musculus",
  msigdb_db_species = "MM",
  kegg_organism = "mmu",
  progeny_organism = "mouse",
  collectri_organism = "mouse",
  output_dir = "pathway_library_mouse"
)
```

### 5. Build a rat library

Rat builds also align the species-specific settings:

- `species` is canonicalized to `Rattus norvegicus`
- `msigdb_db_species` defaults to `HS`
- rat MSigDB uses ortholog mapping from the human MSigDB database
- rat MSigDB filenames therefore keep the human-style collection codes such as `msigdb_h_rattus_norvegicus.gmt`
- `kegg_organism` is set to `rno`
- `collectri_organism` is set to `rat`
- PROGENy is skipped with a warning
- CollecTRI emits a warning because the rat path uses a fallback workaround and may still fail

```r
source("R/load_project_code.R")
library(org.Rn.eg.db)

build_pathway_library(
  org_db = org.Rn.eg.db,
  species = "Rattus norvegicus",
  output_dir = "pathway_library_rat"
)
```

If you want to be explicit about the rat-safe path:

```r
source("R/load_project_code.R")
library(org.Rn.eg.db)

build_pathway_library(
  org_db = org.Rn.eg.db,
  species = "Rattus norvegicus",
  msigdb_db_species = "HS",
  kegg_organism = "rno",
  collectri_organism = "rat",
  skip = "progeny",
  output_dir = "pathway_library_rat"
)
```

### 6. Build only GO libraries

Useful when you are tuning GO-specific methodology and do not want network traffic or non-GO dependencies on the critical path.

```r
source("R/load_project_code.R")
library(org.Hs.eg.db)

build_go_regulation_library(
  org_db = org.Hs.eg.db,
  output_dir = "go_dedup_gmt"
)

build_go_dedup_library(
  org_db = org.Hs.eg.db,
  ontology = "CC",
  output_dir = "go_dedup_gmt"
)
```

### 7. Run ORA and GSEA with the helper layer

```r
source("R/load_project_code.R")
library(org.Hs.eg.db)

libs <- load_pathway_libraries("pathway_library")

signed_rank <- c("7157" = 5.2, "1956" = -4.0, "5290" = 3.1)

ora_result <- run_ora(
  signed_rank = signed_rank,
  gene_fraction = 0.10,
  org_db = org.Hs.eg.db,
  signed_library = libs$signed,
  unsigned_library = libs$unsigned
)

gsea_result <- run_gsea(
  signed_rank = signed_rank,
  org_db = org.Hs.eg.db,
  signed_library = libs$signed,
  unsigned_library = libs$unsigned
)
```

## ORA And GSEA Interpretation Notes

- Signed libraries are intended for directional interpretation.
- Unsigned libraries are intended for magnitude-based interpretation.
- In `run_ora()`, the ORA input list is selected explicitly by `gene_fraction` from the ranked vector.
- In `run_gsea()`, signed libraries use the signed statistic, while unsigned libraries use the absolute statistic.

That distinction is important: unsigned enrichments should be interpreted as strong perturbation, not necessarily upregulation.

## Dependencies

Common packages used by the pipeline include:

- `GO.db`
- `AnnotationDbi`
- `Matrix`
- `org.Hs.eg.db`
- `org.Mm.eg.db`
- `org.Rn.eg.db`
- `msigdbr`
- `reactome.db`
- `progeny`
- `decoupleR`
- `clusterProfiler`
- `testthat`

## Offline Versus Network-Backed Steps

- Offline-capable: GO, MSigDB via `msigdbr`, Reactome via `reactome.db`, PROGENy, and most unit tests.
- Network-backed: KEGG and CollecTRI.

That split matters when debugging in restricted environments: the smoke test intentionally exercises the offline-capable path first.

## Testing

Run the full test suite with:

```r
source("R/load_project_code.R")
testthat::test_dir("tests/testthat")
```

What the tests cover:

- GMT read/write helpers and ID-mapping behavior.
- Species canonicalization and MSigDB database selection.
- GO sibling pairing and Jaccard-dedup selection behavior.
- ORA and GSEA input-preparation helpers.
- An end-to-end smoke build for the offline-capable pipeline.

## Troubleshooting

- If a build fails in KEGG or CollecTRI, retry with `skip = c("kegg", "collectri")` to isolate the offline path.
- If a rat build prints a PROGENy warning, that is expected: the builder is protecting you from an unsupported rat PROGENy configuration.
- If a rat build prints a CollecTRI warning, that is also expected: the pipeline is trying a fallback around the known `get_collectri(organism = "rat")` failure mode.
- If you need mouse-native MSigDB collections, keep `msigdb_db_species = "MM"` or let the mouse defaults choose it automatically.
- If you need rat MSigDB, expect human MSigDB ortholog mapping through `msigdbr`.
- If a downstream analysis fails because of gene identifiers, confirm your inputs are Entrez IDs or convert them before enrichment.
- If you want to inspect the exact settings used for the checked-in library, read [`pathway_library/README.md`](./pathway_library/README.md).
- If you want to change biological assumptions, start in `R/go_library.R` and `R/pathway_downloaders.R` instead of editing the wrapper scripts.

## Recommended Entry Points

- Use `build_human_pathway_library.R` to regenerate the checked-in human artifacts.
- Use `build_pathway_library()` for a normal full rebuild from code.
- Use `build_go_regulation_library()` and `build_go_dedup_library()` when working specifically on GO methodology.
- Use `run_ora()` and `run_gsea()` when you want downstream analyses to follow the same signed-versus-unsigned conventions as the repo.
