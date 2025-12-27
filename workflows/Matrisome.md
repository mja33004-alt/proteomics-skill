# Matrisome Workflow

**Trigger:** "matrisome analysis", "ECM proteins", "extracellular matrix"

---

## Overview

Analyze extracellular matrix (ECM) proteins in your dataset using the MatrisomeDB reference database. Identifies collagens, glycoproteins, proteoglycans, and matrisome-associated proteins.

---

## Matrisome Categories

| Division | Category | Examples |
|----------|----------|----------|
| **Core matrisome** | Collagens | COL1A1, COL4A1, COL6A1 |
| | ECM Glycoproteins | FN1, LAMA1, FBLN1 |
| | Proteoglycans | DCN, BGN, VCAN |
| **Matrisome-associated** | ECM-affiliated Proteins | ANXA2, GPC1 |
| | ECM Regulators | MMP2, TIMP1, LOX |
| | Secreted Factors | TGFB1, VEGFA, BMP2 |

---

## Reference Data

Located in `~/.claude/Skills/Proteomics/data/`:
- `matrisome_hs_masterlist.csv` - Human matrisome (1,000+ proteins)
- `matrisome_mm_masterlist.csv` - Mouse matrisome

---

## Execution Steps

### Step 1: Gather Information

Ask user for:
1. **Comparison data file** - CSV with fold change/q-value
2. **Organism** - "human" or "mouse"
3. **Output directory**

### Step 2: Option A - Use Plot_Workup Script

The Matrisome function is included in `Plot_Workup_V10.R`:

```bash
cd [PROJECT_DIR]
Rscript ~/.claude/Skills/Proteomics/rscripts/Plot_Workup_V10.R
```

### Step 3: Option B - Generate Custom Code

```r
library(tidyverse)
library(ggplot2)

# Parameters
organism <- "human"  # or "mouse"
skill_data <- "~/.claude/Skills/Proteomics/data"

# Load matrisome reference
if (organism == "mouse") {
  ECM <- read.csv(file.path(skill_data, "matrisome_mm_masterlist.csv"), stringsAsFactors = FALSE)
} else {
  ECM <- read.csv(file.path(skill_data, "matrisome_hs_masterlist.csv"), stringsAsFactors = FALSE)
}

# Load comparison data
all_comparisons <- read.csv("data/candidates.csv", stringsAsFactors = FALSE)

# Join with matrisome data
my_ECM <- right_join(all_comparisons, ECM, by = c("Genes" = "Gene.Symbol"))

# Select relevant columns
my_ECM <- my_ECM[, c("Genes", "ProteinDescriptions", "Division", "Category",
                     "Comparison..group1.group2.", "AVG.Log2.Ratio", "Qvalue",
                     "GO.Cellular.Component")]

# Pivot for multiple comparisons
my_ECM_wide <- my_ECM %>%
  pivot_wider(
    id_cols = c(Genes, ProteinDescriptions, Division, Category, GO.Cellular.Component),
    names_from = Comparison..group1.group2.,
    values_from = c(AVG.Log2.Ratio, Qvalue)
  )

# Save full results
write.csv(my_ECM_wide, "output/matrisome_analysis.csv")

# Filter to quantified proteins only
my_ECM_quant <- my_ECM_wide %>% filter(!is.na(ProteinDescriptions))

# Set category order
my_ECM_quant$Category <- factor(my_ECM_quant$Category,
  levels = c("Collagens", "ECM Glycoproteins", "Proteoglycans",
             "ECM-affiliated Proteins", "ECM Regulators", "Secreted Factors"))

# Generate histogram
tiff("output/matrisome_histogram.tiff", res = 300, height = 5, width = 7, units = "in")
my_ECM_quant %>%
  group_by(Category) %>%
  count() %>%
  ggplot(aes(x = Category, y = n, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c(
    alpha('#FD6467', 1), alpha('#FD6467', 0.5), alpha('#FD6467', 0.2),
    alpha('#5B1A18', 1), alpha('#5B1A18', 0.5), alpha('#5B1A18', 0.2)
  )) +
  scale_y_continuous(name = "Number of proteins", expand = expansion(mult = c(0, 0.1))) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 16, color = "black"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    legend.position = "none"
  )
dev.off()

# Summary
cat("\nMatrisome Analysis Summary:\n")
print(table(my_ECM_quant$Division, my_ECM_quant$Category))
```

---

## Output Files

| File | Description |
|------|-------------|
| `matrisome_analysis.csv` | All matrisome proteins with fold change/q-values |
| `matrisome_histogram.tiff` | Bar chart by category |

---

## Example Output

```
User: "Run matrisome analysis on my kidney data"

Claude: Running the **Matrisome** workflow from the **Proteomics** skill...

        Matrisome Analysis Summary
        ==========================
        Total matrisome proteins quantified: 187

        Core Matrisome:
        - Collagens: 23
        - ECM Glycoproteins: 67
        - Proteoglycans: 12

        Matrisome-associated:
        - ECM-affiliated Proteins: 34
        - ECM Regulators: 28
        - Secreted Factors: 23

        Top significantly altered (q < 0.01):
        - COL1A1: log2FC = 2.34 (Up in Diseased)
        - FN1: log2FC = 1.89 (Up in Diseased)
        - MMP9: log2FC = 1.67 (Up in Diseased)

        Output files:
        - output/matrisome_analysis.csv
        - output/matrisome_histogram.tiff
```
