# Normalize Workflow

**Trigger:** "normalize data", "apply normalization", "median/quantile/loess normalize"

---

## Overview

Apply normalization methods to proteomics intensity data to correct for technical variation between samples.

## Available Methods

| Method | Use Case | Description |
|--------|----------|-------------|
| **Median** | General use, robust | Scales each sample so all medians are equal |
| **Quantile** | Cross-sample comparison | Forces identical intensity distributions |
| **Cyclic Loess** | Paired samples | Iterative local regression between pairs |
| **Lowess** | General smoothing | Local weighted regression |
| **VSN** | High variance data | Variance stabilization using arsinh transform |
| **Sample Loading** | TMT/labeling experiments | Scales based on total protein loading |

---

## Execution Steps

### Step 1: Gather Information

Ask user for:
1. **Input file path** - CSV with protein intensities (columns = samples, rows = proteins)
2. **Sample column pattern** - Regex to identify sample columns (e.g., `"JB\\d_\\d+"`)
3. **Desired normalization method(s)**

### Step 2: Option A - Use Existing Script

For TMT multi-batch experiments or comprehensive normalization:

```bash
cd [PROJECT_DIR]
Rscript ~/.claude/Skills/Proteomics/rscripts/normalization/Step_1_Normalization.R
```

**Note:** Script may need modification for specific input format.

### Step 3: Option B - Generate Custom R Code

For custom normalization, generate R code following these patterns:

**Median Normalization:**
```r
library(tidyverse)

# Read data
iMat <- as.matrix(read.csv("data/protein_intensities.csv", row.names=1))

# Median normalization
meds <- apply(iMat, 2, median, na.rm=TRUE)
MnMat <- sweep(iMat, 2, meds/mean(meds), FUN="/")

# Write output
write.csv(MnMat, "output/median_normalized.csv")
```

**Quantile Normalization:**
```r
library(preprocessCore)

iMat <- as.matrix(read.csv("data/protein_intensities.csv", row.names=1))
QnMat <- 2^normalize.quantiles(log2(iMat))
colnames(QnMat) <- colnames(iMat)
rownames(QnMat) <- rownames(iMat)

write.csv(QnMat, "output/quantile_normalized.csv")
```

**Cyclic Loess Normalization:**
```r
library(limma)

iMat <- as.matrix(read.csv("data/protein_intensities.csv", row.names=1))
CLnMat <- 2^normalizeCyclicLoess(data.frame(log2(iMat)), method = "pairs")

write.csv(CLnMat, "output/cyclic_loess_normalized.csv")
```

**VSN (Variance Stabilizing):**
```r
library(vsn)

iMat <- as.matrix(read.csv("data/protein_intensities.csv", row.names=1))
VSNnMat <- 2^justvsn(iMat)

write.csv(VSNnMat, "output/vsn_normalized.csv")
```

### Step 4: Verify Output

After normalization, verify:
- Output file exists and is readable
- No NA/Inf values introduced
- Column names (samples) preserved
- Row names (protein IDs) preserved
- Boxplot of log2 intensities should be aligned

---

## Required R Packages

```r
install.packages(c("tidyverse"))
BiocManager::install(c("limma", "preprocessCore", "vsn"))
```

---

## Example Execution

```
User: "Normalize my proteomics data using quantile normalization"

Claude: Running the **Normalize** workflow from the **Proteomics** skill...

        What is the path to your protein intensity CSV file?

User: ~/Desktop/RCCN_Kidney/data/protein_report.csv

Claude: [Generates R script]
        [Executes: Rscript normalize_quantile.R]

        Normalization complete. Output saved to:
        ~/Desktop/RCCN_Kidney/output/quantile_normalized.csv

        Summary:
        - 5,432 proteins processed
        - 12 samples normalized
        - No missing values introduced
```
