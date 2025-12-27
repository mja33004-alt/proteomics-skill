# SaspAnalysis Workflow

**Trigger:** "SASP analysis", "senescence factors", "core SASP"

---

## Overview

Analyze Senescence-Associated Secretory Phenotype (SASP) factors in your proteomics dataset. SASP proteins are secreted by senescent cells and include cytokines, chemokines, growth factors, and proteases.

---

## Core SASP Reference Data

Located in `~/.claude/Skills/Proteomics/data/`:
- `Human_Core_SASP.csv` - Human SASP factors with scores
- `Mouse_Core_SASP.csv` - Mouse SASP orthologs

**Score columns:**
- **IR** - Ionizing Radiation-induced senescence
- **RAS** - Oncogene-induced senescence (RAS)
- **ATV** - Atazanavir-induced senescence

---

## Key SASP Factors

| Category | Examples |
|----------|----------|
| Cytokines | IL1A, IL1B, IL6, IL8 (CXCL8) |
| Chemokines | CXCL1, CXCL2, CCL2 (MCP1) |
| Growth factors | GDF15, IGFBP3, IGFBP7 |
| Proteases | MMP1, MMP3, MMP9, MMP12 |
| Other | SERPINE1 (PAI-1), FN1, TIMP1 |

---

## Execution Steps

### Step 1: Gather Information

Ask user for:
1. **Protein intensity file** - For protein-level analysis
2. **Comparison file** - For fold change analysis
3. **Organism** - "human" or "mouse"
4. **Q-value threshold** - Default: 0.01

### Step 2: Option A - Use Plot_Workup Script

The Core_SASP function is included in `Plot_Workup_V10.R`:

```bash
cd [PROJECT_DIR]
Rscript ~/.claude/Skills/Proteomics/rscripts/Plot_Workup_V10.R
```

### Step 3: Option B - Generate Custom Code

```r
library(tidyverse)
library(gplots)

# Parameters
organism <- "human"
myQval <- 0.01
myFC <- 0.58
skill_data <- "~/.claude/Skills/Proteomics/data"

# Load SASP reference
if (organism == "human") {
  SASP <- read.csv(file.path(skill_data, "Human_Core_SASP.csv"), stringsAsFactors = FALSE)
} else {
  SASP <- read.csv(file.path(skill_data, "Mouse_Core_SASP.csv"), stringsAsFactors = FALSE)
}

# Load protein intensity data (for heatmap)
protein <- read.csv("data/protein_report.csv", stringsAsFactors = FALSE)
condition <- read.csv("data/ConditionSetup.csv", stringsAsFactors = FALSE)

# Prepare intensity matrix
iMat <- as.matrix(protein[, grep("PG.Quantity", names(protein))])
rownames(iMat) <- protein$PG.Genes

# Select SASP proteins in dataset
my_SASP <- iMat[rownames(iMat) %in% SASP$Genes, ]
cat(paste("Found", nrow(my_SASP), "Core SASP proteins in dataset\n"))

# Heatmap of SASP protein intensities
myRamp <- colorRampPalette(colors = c("#0571b0", "#f7f7f7", "#ca0020"))
myheatcol <- condition$Color

tiff("output/SASP_proteins_heatmap.tiff", res = 300, height = 8, width = 8, units = "in")
heatmap.2(
  t(scale(t(log10(my_SASP)))),
  col = myRamp,
  trace = "none",
  labRow = rownames(my_SASP),
  ColSideColors = myheatcol,
  margins = c(8, 8)
)
dev.off()

# Save protein intensities
write.csv(my_SASP, "output/SASP_protein_intensities.csv")

# Load comparison data
all_comparisons <- read.csv("data/candidates.csv", stringsAsFactors = FALSE)

# Select SASP proteins in comparison data
my_SASP_comp <- all_comparisons[all_comparisons$Genes %in% SASP$Genes, ]

# Add significance coloring
my_SASP_comp$Color <- ifelse(
  my_SASP_comp$AVG.Log2.Ratio >= myFC & my_SASP_comp$Qvalue < myQval, "Red",
  ifelse(my_SASP_comp$AVG.Log2.Ratio <= -myFC & my_SASP_comp$Qvalue < myQval, "Blue", "Grey")
)

# Pivot for multiple comparisons
SASP_wide <- my_SASP_comp %>%
  select(Genes, Comparison..group1.group2., AVG.Log2.Ratio, Qvalue, Color) %>%
  pivot_wider(
    id_cols = Genes,
    names_from = Comparison..group1.group2.,
    values_from = c(AVG.Log2.Ratio, Qvalue, Color)
  )

# Save comparison data
write.csv(SASP_wide, "output/SASP_comparison.csv")

# Comparison heatmap
temp <- my_SASP_comp %>%
  select(Genes, Comparison..group1.group2., AVG.Log2.Ratio) %>%
  pivot_wider(id_cols = Genes, names_from = Comparison..group1.group2., values_from = AVG.Log2.Ratio)

temp_mat <- as.matrix(temp[, -1])
rownames(temp_mat) <- temp$Genes

tiff("output/SASP_comparison_heatmap.tiff", res = 300, height = 8, width = 8, units = "in")
heatmap.2(
  temp_mat,
  col = myRamp,
  trace = "none",
  labRow = rownames(temp_mat),
  margins = c(12, 8)
)
dev.off()

# Summary
cat("\nSASP Analysis Summary:\n")
cat(paste("Total Core SASP proteins:", nrow(SASP), "\n"))
cat(paste("Quantified in dataset:", nrow(my_SASP), "\n"))
cat(paste("Significantly altered:", sum(my_SASP_comp$Color != "Grey"), "\n"))
```

---

## Output Files

| File | Description |
|------|-------------|
| `SASP_protein_intensities.csv` | Intensity values for SASP proteins |
| `SASP_proteins_heatmap.tiff` | Heatmap of SASP protein intensities |
| `SASP_comparison.csv` | Fold change and significance for each comparison |
| `SASP_comparison_heatmap.tiff` | Heatmap of log2 fold changes |

---

## Interpretation

**Strong SASP signature indicators:**
- Multiple cytokines elevated (IL1B, IL6, IL8)
- MMP family members upregulated
- GDF15 strongly induced
- SERPINE1 (PAI-1) elevated

**Context considerations:**
- SASP can be beneficial (wound healing, immune clearance)
- SASP can be detrimental (chronic inflammation, cancer promotion)
- Tissue and disease context matters

---

## Example Output

```
User: "Run SASP analysis on my aging kidney samples"

Claude: Running the **SaspAnalysis** workflow from the **Proteomics** skill...

        SASP Analysis Summary
        =====================
        Core SASP proteins in reference: 175
        Quantified in your dataset: 47 (26.9%)
        Significantly altered (q < 0.01): 23

        Top upregulated SASP factors:
        | Gene | Log2FC | q-value | IR score | RAS score |
        |------|--------|---------|----------|-----------|
        | GDF15 | 2.34  | 0.0001  | 2.54     | 5.07      |
        | MMP1  | 1.89  | 0.0012  | 6.87     | 7.85      |
        | CXCL1 | 1.67  | 0.0023  | 8.65     | 9.23      |
        | IL1B  | 1.45  | 0.0034  | 7.23     | 8.12      |

        Top downregulated SASP factors:
        | Gene  | Log2FC | q-value |
        |-------|--------|---------|
        | TIMP2 | -1.23  | 0.0056  |
        | IGFBP5| -0.98  | 0.0078  |

        Interpretation:
        Strong SASP signature detected with elevated cytokines (IL1B),
        chemokines (CXCL1), and proteases (MMP1). Consider senescence-related
        mechanisms in your experimental context.

        Output files:
        - output/SASP_proteins_heatmap.tiff
        - output/SASP_comparison_heatmap.tiff
        - output/SASP_comparison.csv
```
