# Pathway Analysis Guide

This document provides guidance on performing pathway enrichment analysis using KEGG, GO, and ConsensusPathDB.

---

## Overview

Pathway analysis identifies biological processes, molecular functions, and signaling pathways enriched in a set of proteins (typically differentially expressed).

| Database | Type | Organisms | Access |
|----------|------|-----------|--------|
| KEGG | Pathways | 500+ species | clusterProfiler, web |
| GO | Ontology (BP, MF, CC) | All | clusterProfiler, web |
| ConsensusPathDB | Integrated | Human, mouse, yeast | Python client, web |
| Reactome | Pathways | Human, mouse | ReactomePA, web |

---

## Input Preparation

### Filter Significant Proteins

```r
# Filter by significance thresholds
myQval <- 0.01
myFC <- 0.58

significant <- all_comparisons[
  all_comparisons$Qvalue < myQval &
  abs(all_comparisons$AVG.Log2.Ratio) >= myFC,
]

# Separate up and down regulated
up_regulated <- significant[significant$AVG.Log2.Ratio > 0, ]
down_regulated <- significant[significant$AVG.Log2.Ratio < 0, ]

# Get gene lists
all_genes <- significant$Genes
up_genes <- up_regulated$Genes
down_genes <- down_regulated$Genes
```

### Convert Gene IDs

```r
library(clusterProfiler)
library(org.Hs.eg.db)  # Human
# library(org.Mm.eg.db)  # Mouse

# Convert gene symbols to Entrez IDs
gene_ids <- bitr(
  all_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# For GO/KEGG with gene names
gene_ids_symbol <- bitr(
  all_genes,
  fromType = "SYMBOL",
  toType = c("ENTREZID", "ENSEMBL", "UNIPROT"),
  OrgDb = org.Hs.eg.db
)
```

---

## KEGG Pathway Analysis

```r
library(clusterProfiler)

# Enrichment analysis
kegg_result <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = "hsa",  # "mmu" for mouse
  keyType = "kegg",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  minGSSize = 10,
  maxGSSize = 500
)

# View results
head(as.data.frame(kegg_result))

# Visualizations
dotplot(kegg_result, showCategory = 20, title = "KEGG Pathway Enrichment")
barplot(kegg_result, showCategory = 20)
cnetplot(kegg_result, categorySize = "pvalue", foldChange = NULL)
emapplot(kegg_result)

# Save results
write.csv(as.data.frame(kegg_result), "kegg_results.csv")
ggsave("kegg_dotplot.png", dotplot(kegg_result, showCategory = 20),
       width = 10, height = 8, dpi = 300)
```

### KEGG with Fold Change

```r
# Create named vector of fold changes
gene_fc <- significant$AVG.Log2.Ratio
names(gene_fc) <- gene_ids$ENTREZID[match(significant$Genes, gene_ids$SYMBOL)]
gene_fc <- sort(gene_fc, decreasing = TRUE)

# Gene Set Enrichment Analysis (GSEA)
kegg_gsea <- gseKEGG(
  geneList = gene_fc,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05
)

gseaplot2(kegg_gsea, geneSetID = 1:3)
```

---

## GO Enrichment Analysis

```r
library(clusterProfiler)

# Biological Process
go_bp <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)

# Cellular Component
go_cc <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)

# Molecular Function
go_mf <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)

# All three combined
go_all <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pvalueCutoff = 0.05,
  readable = TRUE
)

# Visualizations
dotplot(go_bp, showCategory = 20, title = "GO Biological Process")
dotplot(go_cc, showCategory = 20, title = "GO Cellular Component")
dotplot(go_mf, showCategory = 20, title = "GO Molecular Function")
```

---

## ConsensusPathDB

ConsensusPathDB integrates multiple pathway databases. Use the Python client in `python/client.py`.

### Web Interface

1. Go to http://cpdb.molgen.mpg.de/
2. Upload gene list (one per line or UniProt IDs)
3. Select analysis type (over-representation, enrichment)
4. Choose databases (KEGG, Reactome, BioCarta, etc.)
5. Download results

### Python Client (SOAP)

```python
# Located at: ~/.claude/Skills/Proteomics/python/client.py

from cpdb_services import cpdbServicesPortTypeService
from cpdb_services_types import *

# Initialize client
client = cpdbServicesPortTypeService()

# Map gene identifiers to CPDB IDs
mapped = client.mapAccessionNumbers(
    accessionNumbers=["P12345", "Q67890"],
    accType="uniprot"
)

# Over-representation analysis
result = client.overRepresentationAnalysis(
    entityType="genes",
    fsetType="P",  # P = pathways
    cpdbIds=mapped.cpdbIds,
    accType="uniprot",
    pThreshold=0.01
)
```

---

## Reactome Pathway Analysis

```r
library(ReactomePA)

reactome_result <- enrichPathway(
  gene = gene_ids$ENTREZID,
  organism = "human",
  pvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(reactome_result, showCategory = 20)
viewPathway("R-HSA-109582", foldChange = gene_fc, readable = TRUE)
```

---

## Interpretation Guidelines

### Significance Metrics

| Metric | Description | Threshold |
|--------|-------------|-----------|
| p-value | Enrichment significance | < 0.05 |
| q-value (padj) | FDR-adjusted p-value | < 0.1 or 0.05 |
| Count | Number of genes in pathway | Higher is better |
| GeneRatio | Proportion of pathway genes | Higher = more specific |

### Common Pitfalls

1. **Multiple testing:** Always use adjusted p-values
2. **Background set:** Use expressed genes as background, not all genes
3. **Redundant terms:** GO terms are hierarchical; use simplify()
4. **Pathway size:** Very large pathways (>500 genes) are less informative

### Reducing Redundancy

```r
# Simplify GO results (remove redundant terms)
go_bp_simplified <- simplify(
  go_bp,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

dotplot(go_bp_simplified, showCategory = 20)
```

---

## Required R Packages

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db",      # Human annotation
  "org.Mm.eg.db",      # Mouse annotation
  "DOSE",              # Disease ontology
  "enrichplot",        # Visualization
  "ReactomePA",        # Reactome pathways
  "pathview"           # KEGG pathway visualization
))
```

---

## Output Files

| File | Contents |
|------|----------|
| `kegg_results.csv` | KEGG enrichment results table |
| `go_bp_results.csv` | GO Biological Process results |
| `kegg_dotplot.png` | KEGG dot plot visualization |
| `go_bp_dotplot.png` | GO BP dot plot |
| `pathway_network.png` | Enrichment map (pathway similarity) |
