# ExcelWorkup Workflow

**Trigger:** "create Excel report", "filter by q-value", "generate data tables"

---

## Overview

Generate formatted Excel workbooks containing protein identifications, comparison data, and filtered results by significance thresholds.

---

## Execution Steps

### Step 1: Gather Information

Ask user for:
1. **Protein report file** - CSV with identified proteins
2. **Comparison file** - CSV with fold change/q-value data
3. **Output filename** - Where to save Excel file
4. **Batch name** - For sheet labeling
5. **Significance thresholds** - Which q-value filters (0.05, 0.01, 0.001)
6. **Fold-change cutoff** - Default: 0.58

### Step 2: Option A - Use Existing Script

```bash
cd [PROJECT_DIR]
Rscript ~/.claude/Skills/Proteomics/rscripts/Excel_Workup_v05.R
```

**Required modifications in script:**
- Line 10: `setwd("~/Desktop/[YOUR_PROJECT]/")`
- Line 12: `myoutput <- "output/Data_Tables/[OUTPUT_FILE].xlsx"`
- Line 14: `batch <- "[YOUR_BATCH]"`
- Lines 19-21: Set q-value filter flags
- Line 24: `myFC <- 0.58`
- Lines 26-28: Update input file paths

### Step 3: Option B - Generate Custom Code

```r
library(openxlsx)
library(tidyverse)
library(scales)
library(stringr)

# Parameters
myoutput <- "output/Data_Tables/analysis.xlsx"
batch <- "MyBatch"
myFC <- 0.58

# Read data
protein <- read.csv("data/protein_report.csv", stringsAsFactors = FALSE)
all_comparisons <- read.csv("data/candidates.csv", stringsAsFactors = FALSE)

# Format protein count
proIDs <- paste(comma(nrow(protein)), "Protein Groups Identified with >= 2 Unique Peptides")

# Initialize workbook
wb <- createWorkbook()

# Sheet 1: All Identified Proteins
addWorksheet(wb, sheetName = "Identified Protein Groups")
writeData(wb, sheet = 1, x = paste(batch, "- Protein Group Identifications"), startRow = 1)
writeData(wb, sheet = 1, x = proIDs, startRow = 2)
writeData(wb, sheet = 1, x = protein, startRow = 4)

# Sheet 2: All Comparisons (unfiltered)
addWorksheet(wb, sheetName = "All Comparisons")
writeData(wb, sheet = 2, x = paste(batch, "- All Comparisons"), startRow = 1)
writeData(wb, sheet = 2, x = proIDs, startRow = 2)
writeData(wb, sheet = 2, x = "Altered Protein Groups with No Filter", startRow = 3)
writeData(wb, sheet = 2, x = all_comparisons, startRow = 5)

# Function to add filtered sheets
add_filtered_sheet <- function(myQval) {
  filtered <- all_comparisons[
    all_comparisons$Absolute.AVG.Log2.Ratio >= myFC &
    all_comparisons$Qvalue <= myQval,
  ]

  sheet_name <- paste("q <=", myQval)
  addWorksheet(wb, sheetName = sheet_name)

  writeData(wb, sheet = sheet_name,
            x = paste(batch, "- Significantly Altered"),
            startRow = 1)
  writeData(wb, sheet = sheet_name,
            x = paste(comma(nrow(filtered)),
                     "proteins with |log2(FC)| >=", myFC,
                     "& q-value <=", myQval),
            startRow = 2)
  writeData(wb, sheet = sheet_name, x = filtered, startRow = 4)
}

# Add sheets for each threshold
add_filtered_sheet(0.05)
add_filtered_sheet(0.01)
add_filtered_sheet(0.001)

# Save workbook
saveWorkbook(wb, myoutput, overwrite = TRUE)

# Optionally open
# openXL(myoutput)
```

---

## Expected Input Format

**Protein Report (from Spectronaut/DIA-NN):**
- `PG.Genes` - Gene symbols
- `PG.ProteinDescriptions` - Protein names
- `PG.Qvalue` - Global q-value
- `PG.UniProtIds` - UniProt accessions
- `PG.Quantity.*` columns - Sample intensities

**Comparison/Candidates File:**
- `Comparison..group1.group2.` - Comparison name
- `Genes` - Gene symbols
- `AVG.Log2.Ratio` - Log2 fold change
- `Qvalue` - Statistical significance
- `Absolute.AVG.Log2.Ratio` - Absolute fold change

---

## Output Structure

| Sheet | Contents |
|-------|----------|
| Identified Protein Groups | All proteins with >= 2 unique peptides |
| All Comparisons | Unfiltered comparison data |
| q <= 0.05 | Filtered: \|log2FC\| >= 0.58 AND q <= 0.05 |
| q <= 0.01 | Filtered: \|log2FC\| >= 0.58 AND q <= 0.01 |
| q <= 0.001 | Filtered: \|log2FC\| >= 0.58 AND q <= 0.001 |

---

## Required R Packages

```r
install.packages(c("openxlsx", "tidyverse", "scales", "stringr"))
```

---

## Example Output

```
User: "Generate an Excel report with q-value filtering"

Claude: Running the **ExcelWorkup** workflow from the **Proteomics** skill...

        Created Excel workbook: output/Data_Tables/RCCN2_analysis.xlsx

        Sheets created:
        1. Identified Protein Groups (5,432 proteins)
        2. All Comparisons (16,296 entries)
        3. q <= 0.05 (1,245 significant proteins)
        4. q <= 0.01 (876 significant proteins)
        5. q <= 0.001 (423 significant proteins)

        Parameters used:
        - Fold change threshold: |log2FC| >= 0.58 (1.5x)
        - Batch: RCCN2
```
