# PathwayAnalysis Workflow

**Trigger:** "pathway analysis", "KEGG enrichment", "ConsensusPathDB", "GO enrichment"

---

## Overview

Perform pathway enrichment analysis on differentially expressed proteins using KEGG, GO, or ConsensusPathDB.

---

## Available Methods

| Method | Tool | Best For |
|--------|------|----------|
| KEGG | clusterProfiler | Metabolic and signaling pathways |
| GO | clusterProfiler | Gene Ontology (BP, MF, CC) |
| ConsensusPathDB | Python client | Multiple database integration |

---

## Execution Steps

### Step 1: Prepare Input

Filter proteins to those meeting significance thresholds:
- Default q-value threshold: 0.01
- Default fold-change threshold: 0.58 (1.5x)

### Step 2: Option A - ConsensusPathDB Script

```bash
cd [PROJECT_DIR]
Rscript ~/.claude/Skills/Proteomics/rscripts/ConsensusPathDB_23_0411_v03.R
```

### Step 3: Option B - clusterProfiler (R)

**KEGG Pathway Analysis:**
```r
library(clusterProfiler)
library(org.Hs.eg.db)  # or org.Mm.eg.db for mouse
library(ggplot2)

# Read significant proteins
sig_proteins <- read.csv("data/significant_proteins.csv")
gene_list <- sig_proteins$Genes

# Convert gene symbols to Entrez IDs
gene_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# KEGG enrichment
kegg_result <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = "hsa",  # "mmu" for mouse
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

# Dot plot
dotplot(kegg_result, showCategory = 20, title = "KEGG Pathway Enrichment")
ggsave("output/kegg_dotplot.png", width = 10, height = 8, dpi = 300)

# Save results
write.csv(as.data.frame(kegg_result), "output/kegg_results.csv")
```

**GO Enrichment Analysis:**
```r
# Biological Process
go_bp <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)

dotplot(go_bp, showCategory = 20, title = "GO Biological Process")
ggsave("output/go_bp_dotplot.png", width = 10, height = 8, dpi = 300)

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
```

### Step 4: Option C - ConsensusPathDB (Python)

Using the Python client in `python/`:

```python
# Located at: ~/.claude/Skills/Proteomics/python/client.py
# See client.py for SOAP web service usage

# Basic usage:
# 1. Map accession numbers to CPDB IDs
# 2. Run over-representation analysis
# 3. Retrieve pathway results
```

---

## Required R Packages

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db",    # Human
  "org.Mm.eg.db",    # Mouse
  "DOSE",
  "enrichplot"
))
```

---

## Visualization Options

| Plot Type | Function | Description |
|-----------|----------|-------------|
| Dot plot | `dotplot()` | Size = gene count, color = p-value |
| Bar plot | `barplot()` | Horizontal bars by enrichment |
| Cnet plot | `cnetplot()` | Network of pathways and genes |
| Heatplot | `heatplot()` | Pathway-gene matrix |
| Enrichment map | `emapplot()` | Pathway similarity network |

---

## Example Output

```
User: "Run KEGG pathway analysis on my significant proteins"

Claude: Running the **PathwayAnalysis** workflow from the **Proteomics** skill...

        Found 523 significantly altered proteins
        - Up-regulated: 287
        - Down-regulated: 236

        KEGG Pathway Analysis Results:

        Top 5 Enriched Pathways:
        1. Ribosome (hsa03010) - 45 genes, q = 1.2e-15
        2. Oxidative phosphorylation (hsa00190) - 32 genes, q = 3.4e-10
        3. Metabolic pathways (hsa01100) - 89 genes, q = 5.6e-08
        4. Carbon metabolism (hsa01200) - 28 genes, q = 2.1e-07
        5. Protein processing in ER (hsa04141) - 24 genes, q = 8.9e-06

        Output files:
        - output/kegg_dotplot.png
        - output/kegg_results.csv
```
