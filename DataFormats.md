# Data Formats Reference

This document describes the expected input and output file formats for the Proteomics skill.

---

## Input File Formats

### Protein Report (from Spectronaut/DIA-NN)

**Expected filename pattern:** `*_Protein_Report_2pep.csv`

| Column | Description | Example |
|--------|-------------|---------|
| `PG.Genes` | Gene symbol(s) | "ACTB" or "ACTB;ACTG1" |
| `PG.ProteinDescriptions` | Protein name(s) | "Actin, cytoplasmic 1" |
| `PG.Qvalue` | Global q-value | 0.0001 |
| `PG.ProteinNames` | Short protein names | "ACTB_HUMAN" |
| `PG.UniProtIds` | UniProt accession(s) | "P60709" |
| `PG.BiologicalProcess` | GO BP terms | "cytoskeleton organization" |
| `PG.CellularComponent` | GO CC terms | "cytoplasm" |
| `PG.MolecularFunction` | GO MF terms | "ATP binding" |
| `PG.Quantity.[Sample]` | Intensity per sample | 1234567.89 |

### Comparison/Candidates File

**Expected filename pattern:** `*_candidates_2pep*.csv`

| Column | Description | Example |
|--------|-------------|---------|
| `Comparison..group1.group2.` | Comparison name | "Treatment / Control" |
| `Genes` | Gene symbol | "ACTB" |
| `ProteinDescriptions` | Protein name | "Actin, cytoplasmic 1" |
| `ProteinNames` | Short name | "ACTB_HUMAN" |
| `UniProtIds` | UniProt accession | "P60709" |
| `ProteinGroups` | Protein group ID | "PG.12345" |
| `Group` | Sample group | "Treatment" |
| `AVG.Log2.Ratio` | Log2 fold change | 1.234 |
| `Qvalue` | Statistical significance | 0.001 |
| `Absolute.AVG.Log2.Ratio` | Absolute fold change | 1.234 |
| `Pvalue` | Unadjusted p-value | 0.0001 |
| `X..of.Ratios` | Number of ratios | 6 |

### Condition Setup File

**For batch/condition information and visualization colors**

| Column | Description | Example |
|--------|-------------|---------|
| `Sample` | Sample ID (matches column names) | "RCCN2-01" |
| `Condition` | Treatment/group name | "Cortex" |
| `Color` | Hex color code | "#E41A1C" |

**Example:**
```csv
Sample,Condition,Color
RCCN2-01,Cortex,#E41A1C
RCCN2-02,Cortex,#E41A1C
RCCN2-03,Medulla,#377EB8
RCCN2-04,Medulla,#377EB8
```

---

## Output File Formats

### Generated Excel Workbooks

| Sheet | Contents |
|-------|----------|
| Identified Protein Groups | All proteins with >= 2 unique peptides |
| All Comparisons | Unfiltered comparison data |
| Signif Altered - q <= X | Filtered by q-value threshold |
| [Group1] vs. [Group2] - q <= X | Per-comparison filtered data |

### Generated Plots

| Type | Format | Dimensions | DPI |
|------|--------|------------|-----|
| Volcano | TIFF | 5x5 inches | 300 |
| Heatmap | TIFF | 8x8 inches | 300 |
| PCA | TIFF | 8x8 inches | 300 |
| Correlation | TIFF | variable | 300 |
| LOPIT | TIFF | 8.5x5.74 inches | 300 |
| KEGG dotplot | PNG | 10x8 inches | 300 |
| Matrisome histogram | TIFF | 7x5 inches | 300 |

---

## Reference Data Formats

### MISEV2018_EV_Markers.txt

Tab-separated values with columns:
- Category (1-5)
- Protein (gene symbol)
- Accession (UniProt ID)
- Type/Description

### Core SASP Files (CSV)

| Column | Description |
|--------|-------------|
| Uniprot | UniProt accession |
| Genes | Gene symbol |
| IR | Ionizing radiation score |
| RAS | RAS oncogene score |
| ATV | Atazanavir score |

### Matrisome Files (CSV)

| Column | Description |
|--------|-------------|
| Division | "Core matrisome" or "Matrisome-associated" |
| Category | Subcategory (Collagens, ECM Glycoproteins, etc.) |
| Gene.Symbol | Gene symbol |
| [Additional annotation columns] | UniProt, description, etc. |

---

## Project Directory Structure

For running the full analysis scripts, organize data as:

```
[PROJECT_DIR]/
├── data/
│   ├── [date]_[batch]_Protein_Report_2pep.csv
│   ├── [date]_[batch]_candidates_2pep_v1.csv
│   └── [date]_[batch]_ConditionSetup.csv
└── output/
    ├── Data_Tables/
    │   └── [Excel reports]
    └── [plot files]
```
