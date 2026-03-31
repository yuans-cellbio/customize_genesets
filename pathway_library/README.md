# Pathway Gene Set Library

Auto-generated pathway library for enrichment analysis (ORA / GSEA).

## Build Information

| Field | Value |
|-------|-------|
| Date | 2026-03-26 14:23:52 EDT |
| R version | 4.5.1 |
| Platform | x86_64-w64-mingw32 |

## Database Versions

| Database | Package Version | Source Info |
|----------|----------------|-------------|
| GO.db | 3.21.0 | Source date: 2025-02-06 |
| OrgDb (org.Hs.eg.db) | 3.21.0 | Schema: HUMAN_DB |
| msigdbr | 26.1.0 | MSigDB v26.1 |
| KEGGREST | 1.48.1 | Retrieved: 2026-03-26 |
| reactome.db | 1.92.0 | DB version: 92 |
| progeny | 1.30.0 | - |
| decoupleR (CollecTRI) | 2.14.0 | Retrieved: 2026-03-26 |

## Parameters

```
species = Homo sapiens
kegg_organism = hsa
progeny_organism = human
collectri_organism = human
include_IEA = TRUE
min_size = 15
max_size = 500
dag_overlap_threshold = 0.8
jaccard_threshold = 0.7
msigdb_collections = H, C2, C7
progeny_top_n = 100
collectri_min_targets = 0
collectri_max_targets = 2000
pos_keywords = positive regulation, activation of, upregulation of, induction of, stimulation of, promotion of, enhancement of, potentiation of
neg_keywords = negative regulation, inhibition of, downregulation of, suppression of, repression of, attenuation of, blockade of, sequestering of, desensitization of, deactivation of, inactivation of, degradation of
skipped = none
```

## File Inventory

### GO Biological Process

- **go_homo_sapiens_bp_dedup_signed_all_evidence.gmt** — 813 gene sets (340.8 KB)
- **go_homo_sapiens_bp_dedup_unsigned_all_evidence.gmt** — 2763 gene sets (1320.9 KB)

### GO Cellular Component

- **go_homo_sapiens_cc_dedup_all_evidence.gmt** — 461 gene sets (248.8 KB)

### GO Molecular Function

- **go_homo_sapiens_mf_dedup_all_evidence.gmt** — 744 gene sets (315.8 KB)

### MSigDB

- **msigdb_h_homo_sapiens.gmt** — 50 gene sets (41.4 KB)
- **msigdb_c2_cp_wikipathways_homo_sapiens.gmt** — 925 gene sets (283.5 KB)
- **msigdb_c7_immunesigdb_homo_sapiens.gmt** — 4872 gene sets (5962.1 KB)

### KEGG

- **kegg_hsa.gmt** — 357 gene sets (204.9 KB)

### Reactome

- **reactome_homo_sapiens.gmt** — 2746 gene sets (839.9 KB)

### PROGENy

- **progeny_signed_human_top100.gmt** — 25 gene sets (8.1 KB)
- **progeny_unsigned_human_top100.gmt** — 14 gene sets (7.6 KB)

### CollecTRI

- **collectri_signed_human.gmt** — 2372 gene sets (295.6 KB)
- **collectri_unsigned_human.gmt** — 1186 gene sets (248.4 KB)

## Usage with clusterProfiler

```r
library(clusterProfiler)
source("build_go_regulation_gmt.R")  # for read_gmt_term2name()

# Load any GMT file
gmt       <- read.gmt("path/to/file.gmt")
gmt_names <- read_gmt_term2name("path/to/file.gmt")

# ORA
enricher(gene = my_genes, universe = bg, TERM2GENE = gmt, TERM2NAME = gmt_names)

# GSEA
GSEA(geneList = my_ranked_genes, TERM2GENE = gmt, TERM2NAME = gmt_names)
```

## Gene ID Format

All GMT files use **NCBI Entrez Gene IDs** for consistency.
Convert with `clusterProfiler::bitr()` if your data uses Ensembl or symbols.
