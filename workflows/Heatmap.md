# Heatmap Workflow

**Trigger:** "heatmap", "PCA", "correlation plot", "sample clustering"

---

## Overview

Generate sample-level visualizations including PCA plots, heatmaps, and correlation matrices to assess data quality and sample relationships.

---

## Execution Steps

### Step 1: Gather Information

Ask user for:
1. **Protein intensity file** - CSV with samples as columns
2. **Condition setup file** - CSV mapping samples to conditions/colors
3. **Sample column pattern** - Regex to extract sample IDs
4. **Output directory**

### Step 2: Option A - Use Full Pipeline

The `Plot_Workup_V10.R` script includes PCA, heatmap, and correlation in one run:

```bash
cd [PROJECT_DIR]
Rscript ~/.claude/Skills/Proteomics/rscripts/Plot_Workup_V10.R
```

### Step 3: Option B - Generate Custom Code

**PCA Plot:**
```r
library(ggplot2)
library(ggrepel)

# Read and prepare data
protein <- read.csv("data/protein_report.csv", stringsAsFactors = FALSE)
condition <- read.csv("data/ConditionSetup.csv", stringsAsFactors = FALSE)

# Extract intensity matrix
iMat <- as.matrix(protein[, grep("PG.Quantity", names(protein))])
rownames(iMat) <- protein$PG.Genes

# Median normalize
meds <- apply(iMat, 2, median, na.rm = TRUE)
nMat <- sweep(iMat, 2, meds/mean(meds), FUN = "/")

# Remove incomplete cases and zeros
pcMat <- nMat[complete.cases(nMat), ]
pcMat[pcMat == 0] <- 1

# Run PCA
pcRes <- prcomp(t(log2(pcMat)), center = TRUE, scale. = TRUE)
pcSum <- summary(pcRes)

# Create labels
PC1label <- paste0("PC1, ", round(100 * pcSum$importance["Proportion of Variance", "PC1"], 1), "% of variance")
PC2label <- paste0("PC2, ", round(100 * pcSum$importance["Proportion of Variance", "PC2"], 1), "% of variance")

# Plot
pcPlotFrame <- data.frame(
  treatment = condition$Condition,
  sample = colnames(nMat),
  pcRes$x[, 1:5],
  color = condition$Color
)

ggplot(pcPlotFrame, aes(PC1, PC2, color = color, shape = treatment)) +
  geom_point(size = 3) +
  scale_x_continuous(name = PC1label) +
  scale_y_continuous(name = PC2label) +
  scale_color_identity() +
  stat_ellipse(aes(color = color)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")

ggsave("output/pca_plot.tiff", dpi = 300, height = 8, width = 8, units = "in")
```

**Heatmap:**
```r
library(gplots)

# Using pcMat from PCA code above
myRamp <- colorRampPalette(colors = c("#0571b0", "#f7f7f7", "#ca0020"))
myheatcol <- condition$Color

tiff("output/heatmap.tiff", res = 300, height = 8, width = 8, units = "in")
heatmap.2(
  t(scale(t(log10(pcMat)))),
  col = myRamp,
  trace = "none",
  labRow = FALSE,
  ColSideColors = myheatcol
)
dev.off()
```

**Correlation Plot:**
```r
library(corrplot)

# Calculate correlation matrix
M <- cor(iMat, use = "complete.obs")

# Color palette
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

tiff("output/correlation.tiff", res = 300, height = 8, width = 8, units = "in")
corrplot(M, method = "color", col = col(200), tl.col = "black")
dev.off()

# Save correlation matrix
write.csv(M, "output/correlation_matrix.csv")
```

---

## Condition Setup File Format

The condition setup file should have:

| Column | Description |
|--------|-------------|
| Sample | Sample ID matching column names |
| Condition | Treatment/group name |
| Color | Hex color code (e.g., "#FF0000") |

Example:
```csv
Sample,Condition,Color
RCCN2-01,Cortex,#E41A1C
RCCN2-02,Cortex,#E41A1C
RCCN2-03,Medulla,#377EB8
RCCN2-04,Medulla,#377EB8
```

---

## Required R Packages

```r
install.packages(c("ggplot2", "ggrepel", "gplots", "corrplot", "RColorBrewer"))
```

---

## Output Files

| File | Description |
|------|-------------|
| `[batch]_pca.tiff` | PCA scatter plot with ellipses |
| `[batch]_heatmap.tiff` | Hierarchical clustering heatmap |
| `[batch]_correlation.tiff` | Sample correlation matrix |
| `[batch]_correlation.csv` | Correlation values |
