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
- `build_pathway_library()` and `download_all_pathway_gmt()` now treat `species` as the public species switch and infer the matching `OrgDb` plus per-database organism keys internally.
- External resources are exported into the same GMT-style interface so they can be used side by side in `clusterProfiler` workflows.

## Species Support

- Human: fully supported by the checked-in build and all main workflows.
- Mouse: supported for GO, MSigDB, Reactome, KEGG, PROGENy, and CollecTRI.
- Rat: supported for GO, MSigDB, Reactome, KEGG, and CollecTRI.
- Rat CollecTRI: built from the direct `OmnipathR::collectri(query_type = "interactions", organism = "rat")` route, then normalized into the same signed and unsigned GMT interface used elsewhere in the repo.
- Rat PROGENy: not supported in this pipeline. Rat builds omit that step by design and report it as an informational message.

## Common Use Cases

- Build a defensible human pathway library for transcriptomics or proteomics enrichment analysis.
- Replace inherited GMT files with a versioned build that records package versions and retrieval dates.
- Generate a signed GO BP library for analyses where directionality matters.
- Build mouse-native MSigDB pathway libraries without falling back to human ortholog mapping.
- Build a rat pathway library with the supported rat sources wired automatically and unsupported sources omitted cleanly.
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

## Installation

Install the repo as a local source package:

```r
install.packages("path/to/customize_genesets", repos = NULL, type = "source")
library(customizeGeneSets)
```

If you are working directly from a source checkout without installing the package:

```r
source("R/load_project_code.R")
load_project_code()
```

Once the package is installed, the user-facing functions also have package help pages such as `?build_pathway_library`, `?run_ora`, and `?download_collectri_gmt`.

## Typical Workflows

### 1. Use the checked-in human library

If you just want the current human GMT files, start in `pathway_library/`.

```r
library(clusterProfiler)
library(customizeGeneSets)

gmt <- read.gmt("pathway_library/reactome_homo_sapiens.gmt")
gmt_names <- read_gmt_term2name("pathway_library/reactome_homo_sapiens.gmt")
```

### 2. Rebuild the default human library

```r
library(customizeGeneSets)

build_pathway_library(
  species = "human",
  output_dir = "pathway_library"
)
```

### 3. Rebuild with custom parameters

Use this when you want to change GO filtering thresholds, select different MSigDB collections, or skip network-backed resources during development.

```r
library(customizeGeneSets)

build_pathway_library(
  species = "human",
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
- the matching `OrgDb` is inferred as `org.Mm.eg.db`
- mouse-native MSigDB is selected with `db_species = "MM"`
- requested human-style collection codes such as `H`, `C2`, and `C7` are translated to mouse-native codes such as `MH`, `M2`, and `M7`
- KEGG, PROGENy, and CollecTRI receive the correct mouse organism keys internally

```r
library(customizeGeneSets)

build_pathway_library(
  species = "mouse",
  output_dir = "pathway_library_mouse"
)
```

If you want to use a custom `OrgDb` object anyway, you still can:

```r
library(customizeGeneSets)
library(org.Mm.eg.db)

build_pathway_library(
  species = "mouse",
  org_db = org.Mm.eg.db,
  output_dir = "pathway_library_mouse"
)
```

### 5. Build a rat library

Rat builds also align the species-specific settings:

- `species` is canonicalized to `Rattus norvegicus`
- the matching `OrgDb` is inferred as `org.Rn.eg.db`
- rat CollecTRI uses the direct OmnipathR interactions query and the same normalization logic as the standalone workaround
- rat PROGENy is omitted automatically because there is no rat PROGENy model in this pipeline
- rat MSigDB uses ortholog mapping from the human MSigDB database
- rat MSigDB filenames therefore keep the human-style collection codes such as `msigdb_h_rattus_norvegicus.gmt`
- KEGG and CollecTRI receive the correct rat organism keys internally

```r
library(customizeGeneSets)

build_pathway_library(
  species = "rat",
  output_dir = "pathway_library_rat"
)
```

If you want to pass a custom `OrgDb`, it must match the declared species:

```r
library(customizeGeneSets)
library(org.Rn.eg.db)

build_pathway_library(
  species = "rat",
  org_db = org.Rn.eg.db,
  output_dir = "pathway_library_rat"
)
```

If you need to control database-specific organism routes manually, call the lower-level `download_*()` function for that source instead of the top-level builder.

### 6. Build only GO libraries

Useful when you are tuning GO-specific methodology and do not want network traffic or non-GO dependencies on the critical path.

```r
library(customizeGeneSets)
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
library(customizeGeneSets)
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

## Function Reference

The package exposes a small set of user-facing functions. The sections below are the README version of the package help pages.

### `build_pathway_library()`

Purpose:
Build the full pathway library for one species, including GO, MSigDB, KEGG, Reactome, PROGENy, and CollecTRI unless stages are skipped.

Arguments:
- `org_db`: optional `OrgDb` object or package name. Leave `NULL` to let the function infer the correct `OrgDb` from `species`.
- `species`: one of `human`, `mouse`, or `rat`, or the matching Latin species name.
- `include_IEA`: include or exclude electronically inferred GO annotations.
- `pos_keywords`, `neg_keywords`, `pair_replacements`: control how GO BP signed terms are detected and paired.
- `min_size`, `max_size`, `dag_overlap_threshold`, `jaccard_threshold`: control GO filtering and deduplication behavior.
- `msigdb_collections`, `msigdb_subcollections`: control which MSigDB collections are exported.
- `progeny_top_n`: number of top PROGENy footprint genes retained per pathway.
- `collectri_min_targets`, `collectri_max_targets`: size filter for CollecTRI transcription-factor target sets.
- `skip`: stages to omit, such as `kegg`, `collectri`, or `progeny`.
- `output_dir`: output directory for GMT files and the provenance README.

Returns:
- An invisible named list with per-stage outputs such as `go_bp`, `go_cc`, `msigdb`, `reactome`, `progeny`, and `collectri` when those stages are run.
- `file_inventory`, `versions`, and `readme` entries describing what was built and where the provenance README was written.

### `build_go_regulation_library()`

Purpose:
Build signed and unsigned GO Biological Process libraries using explicit directional keywords.

Arguments:
- `org_db`: required `OrgDb` object for the target species.
- `include_IEA`: include or exclude electronically inferred GO annotations.
- `pos_keywords`, `neg_keywords`, `pair_replacements`: control directional classification and sibling retention.
- `min_size`, `max_size`, `dag_overlap_threshold`, `jaccard_threshold`: control GO filtering and deduplication.
- `output_dir`: output directory for the signed and unsigned GO BP GMT files.

Returns:
- An invisible list containing `signed`, `unsigned`, `pair_map`, `term_df`, `signed_gmt`, and `unsigned_gmt`.

### `build_go_dedup_library()`

Purpose:
Build a deduplicated GO library for one ontology (`BP`, `CC`, or `MF`).

Arguments:
- `org_db`: required `OrgDb` object for the target species.
- `ontology`: GO ontology to export.
- `include_IEA`, `min_size`, `max_size`, `dag_overlap_threshold`, `jaccard_threshold`, `output_dir`: same roles as in the signed GO builder.

Returns:
- An invisible list with `gene_sets`, `term_df`, and `gmt_file`.

### `discover_directional_keywords()`

Purpose:
Inspect GO BP term names for directional phrases that are not yet in the current positive and negative keyword lists.

Arguments:
- `pos_keywords`, `neg_keywords`: current keyword lists.
- `candidate_words`: optional phrases to scan instead of the built-in candidate list.

Returns:
- A data frame describing matched candidate keywords and example GO terms.

### `download_all_pathway_gmt()`

Purpose:
Export the non-GO sources only: MSigDB, KEGG, Reactome, PROGENy, and CollecTRI.

Arguments:
- `org_db`, `species`: same species-resolution behavior as `build_pathway_library()`.
- `collectri_min_targets`, `collectri_max_targets`: target-size filter for CollecTRI.
- `msigdb_collections`: MSigDB collections to export.
- `output_dir`: output directory for all written GMT files.

Returns:
- An invisible named list of output paths from the component download functions.

### Source-specific GMT exporters

`download_msigdb_gmt()`:
- Inputs: `species`, optional `db_species`, `collections`, `subcollections`, `output_dir`.
- Output: invisible named list from requested collection code to GMT filepath.

`download_kegg_gmt()`:
- Inputs: KEGG organism code such as `hsa`, `mmu`, or `rno`, plus `output_dir`.
- Output: invisible GMT filepath.

`download_reactome_gmt()`:
- Inputs: `species`, `method` (`reactome.db` or `msigdbr`), `output_dir`.
- Output: invisible GMT filepath.

`download_progeny_gmt()`:
- Inputs: `org_db`, `organism` (`human` or `mouse`), `top_n`, `output_dir`.
- Output: invisible list with `signed_gmt` and `unsigned_gmt`.

`download_collectri_gmt()`:
- Inputs: `org_db`, `organism`, `split_complexes`, `min_targets`, `max_targets`, `output_dir`.
- Output: invisible list with `signed_gmt` and `unsigned_gmt`.

### `load_pathway_libraries()`

Purpose:
Load a directory of GMT files into the `TERM2GENE` and `TERM2NAME` structures expected by the downstream ORA and GSEA helpers.

Arguments:
- `pathway_dir`: directory containing GMT files.

Returns:
- A list with `signed` and `unsigned` elements, each containing library entries with `lib_name`, `term2gene`, and `term2name`.

### `run_ora()`

Purpose:
Run ORA on signed and unsigned libraries using a named numeric vector of signed gene-level statistics.

Arguments:
- `signed_rank`: named numeric vector of Entrez-based signed statistics.
- `gene_fraction`: fraction of top absolute-ranked genes used for ORA.
- `top_frac_threshold`: backward-compatible alias for `gene_fraction`.
- `org_db`: `OrgDb` object used to make results readable.
- `universe`: optional background gene universe.
- `signed_library`, `unsigned_library`: outputs from `load_pathway_libraries()`.
- `pvalueCutoff`, `qvalueCutoff`, `pAdjustMethod`, `minGSSize`, `maxGSSize`: parameters forwarded to the enrichment workflow.

Returns:
- A combined ORA result data frame, or `NULL` if no enrichment results are returned.

### `run_gsea()`

Purpose:
Run GSEA on signed and unsigned libraries with consistent directional semantics.

Arguments:
- `signed_rank`: named numeric vector of Entrez-based signed statistics.
- `org_db`: `OrgDb` object used to make results readable.
- `signed_library`, `unsigned_library`: outputs from `load_pathway_libraries()`.
- `pvalueCutoff`, `pAdjustMethod`, `minGSSize`, `maxGSSize`, `nPermSimple`: parameters forwarded to the GSEA workflow.

Returns:
- A combined GSEA result data frame, or `NULL` if no enrichment results are returned.

### `read_gmt_term2name()`

Purpose:
Read the term identifier and description columns from a GMT file for use as `TERM2NAME`.

Arguments:
- `filepath`: path to a GMT file.

Returns:
- A two-column data frame with `ID` and `Name`.

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
- `OmnipathR`
- `clusterProfiler`
- `testthat`

## Offline Versus Network-Backed Steps

- Offline-capable: GO, MSigDB via `msigdbr`, Reactome via `reactome.db`, PROGENy, and most unit tests.
- Network-backed: KEGG and CollecTRI.

That split matters when debugging in restricted environments: the smoke test intentionally exercises the offline-capable path first.

## Testing

Run the full test suite with:

```r
system("Rscript tests/testthat.R")
```

What the tests cover:

- GMT read/write helpers and ID-mapping behavior.
- Species canonicalization and MSigDB database selection.
- GO sibling pairing and Jaccard-dedup selection behavior.
- ORA and GSEA input-preparation helpers.
- An end-to-end smoke build for the offline-capable pipeline.

## Troubleshooting

- If a build fails in KEGG or CollecTRI, retry with `skip = c("kegg", "collectri")` to isolate the offline path.
- If a rat build reports that PROGENy was omitted, that is expected and does not invalidate the rest of the rat library.
- If rat CollecTRI fails, treat that as a real network or OmniPath problem rather than an expected warning path.
- If you need mouse-native MSigDB collections, let the mouse defaults choose them automatically.
- If you need rat MSigDB, expect human MSigDB ortholog mapping through `msigdbr`.
- If a downstream analysis fails because of gene identifiers, confirm your inputs are Entrez IDs or convert them before enrichment.
- If you want to inspect the exact settings used for the checked-in library, read [`pathway_library/README.md`](./pathway_library/README.md).
- If you want to change biological assumptions, start in `R/go_library.R` and `R/pathway_downloaders.R` instead of editing the wrapper scripts.

## Recommended Entry Points

- Use `build_human_pathway_library.R` to regenerate the checked-in human artifacts.
- Use `build_pathway_library()` for a normal full rebuild from code.
- Use `build_go_regulation_library()` and `build_go_dedup_library()` when working specifically on GO methodology.
- Use `run_ora()` and `run_gsea()` when you want downstream analyses to follow the same signed-versus-unsigned conventions as the repo.
