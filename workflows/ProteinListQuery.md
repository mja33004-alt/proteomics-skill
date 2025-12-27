# ProteinListQuery Workflow

**Trigger:** "check EV markers", "MISEV proteins", "exosome markers", "blood contaminants", "cross-reference protein list"

---

## Overview

Cross-reference a list of proteins against curated reference databases including MISEV2018 EV markers, exosome markers, blood contaminants, and more.

---

## Available Reference Lists

Located in `~/.claude/Skills/Proteomics/data/`:

| List | File | Format | Key Info |
|------|------|--------|----------|
| MISEV2018 EV Markers | `MISEV2018_EV_Markers.txt` | TSV | 500+ proteins, Categories 1-5 |
| EV Categories | `MISEV2018_EV_Categories.txt` | Text | Category definitions |
| Exosome Markers | `Exosome_Protein_Markers.txt` | TSV | CD63, CD81, CD9, TSG101, etc. |
| Blood Proteins | `Top_10_Blood_Proteins.txt` | Text | Most abundant blood contaminants |
| Apolipoproteins | `Apolipoproteins.txt` | Text | APOA1, APOB, etc. |

---

## MISEV2018 Categories

| Category | Description | Examples |
|----------|-------------|----------|
| 1 | All EVs (non-tissue specific) | Tetraspanins (CD63, CD81, CD9), Integrins |
| 2 | All EVs (membrane-binding ability) | Annexins, HSPs, ESCRT proteins |
| 3 | Purity control markers | Lipoproteins, ribosomal proteins (non-EV) |
| 4 | EV subtypes/pathologic state | Disease-specific markers |
| 5 | Functional EV components | Cytokines, growth factors |

---

## Execution Steps

### Step 1: Identify Query Type

Ask user:
1. What reference list to query against?
2. User's protein list format (gene names, UniProt IDs, or accessions)
3. Whether to return categories/annotations

### Step 2: Load and Cross-Reference

**R Code for MISEV2018 Query:**
```r
# Load reference data
skill_data <- "~/.claude/Skills/Proteomics/data"
misev <- read.delim(file.path(skill_data, "MISEV2018_EV_Markers.txt"), stringsAsFactors = FALSE)

# Load user protein list
user_proteins <- read.csv("my_proteins.csv")$Genes  # or however formatted

# Cross-reference
matched <- misev[misev$Protein %in% user_proteins, ]

# Summary by category
table(matched$Category)

# Detailed output
write.csv(matched, "output/misev_matches.csv")
```

**R Code for Blood Contaminant Check:**
```r
blood_proteins <- readLines(file.path(skill_data, "Top_10_Blood_Proteins.txt"))
blood_proteins <- blood_proteins[blood_proteins != ""]

# Find contaminants in dataset
contaminants <- user_proteins[user_proteins %in% blood_proteins]

if (length(contaminants) > 0) {
  warning(paste("Found", length(contaminants), "potential blood contaminants:", paste(contaminants, collapse = ", ")))
}
```

**R Code for Exosome Marker Check:**
```r
exosome <- read.delim(file.path(skill_data, "Exosome_Protein_Markers.txt"), stringsAsFactors = FALSE)

# Canonical exosome markers
canonical <- c("CD63", "CD81", "CD9", "TSG101", "ALIX", "SDCBP")

matched_canonical <- user_proteins[user_proteins %in% canonical]
cat(paste("Canonical exosome markers found:", paste(matched_canonical, collapse = ", "), "\n"))
```

### Step 3: Generate Report

Output a summary including:
- Total proteins queried
- Number of matches found
- Breakdown by category (for MISEV2018)
- List of matched proteins with annotations

---

## Quick Query Mode

For simple yes/no queries:

```
User: "Is CD63 an EV marker?"

Claude: Yes, CD63 is a Category 1 EV marker (non-tissue specific tetraspanin).
        UniProt: P08962
        Common use: Standard exosome marker (CD63, CD81, CD9 trinity)
        Notes: Lysosomal-associated membrane protein 3 (LAMP-3)
```

---

## Example Execution

```
User: "Check which of my proteins are MISEV2018 EV markers"

Claude: Running the **ProteinListQuery** workflow from the **Proteomics** skill...

        Please provide your protein list file path.

User: ~/Desktop/project/significant_proteins.csv

Claude: MISEV2018 EV Marker Analysis
        ============================
        Total proteins queried: 523
        Matches found: 47 (9.0%)

        By Category:
        - Category 1 (All EVs): 12 proteins
          - Tetraspanins: CD63, CD81, CD9
          - Integrins: ITGB1, ITGA5

        - Category 2 (Membrane-binding): 18 proteins
          - Annexins: ANXA1, ANXA2, ANXA5, ANXA6
          - ESCRT: TSG101, CHMP2A

        - Category 3 (Purity controls): 8 proteins
          - Apolipoproteins: APOA1, APOE (consider contamination)

        - Category 5 (Functional): 9 proteins
          - Growth factors: TGFB1, EGF

        Full results saved to: output/misev_matches.csv
```
