# Proteomics Skill

A Claude Code skill for label-free quantitative proteomics analysis. Provides workflows for data normalization, visualization (volcano plots, heatmaps, PCA), pathway enrichment analysis (KEGG, ConsensusPathDB), and protein list cross-referencing.

## Installation

Clone this repository directly to your Claude Code Skills directory:

```bash
git clone https://github.com/jobburt-labs/proteomics-skill.git ~/.claude/Skills/Proteomics
```

## R Dependencies

This skill requires R with the following packages installed:

### CRAN Packages

```r
install.packages(c(
  "tidyverse",
  "ggplot2",
  "openxlsx",
  "pheatmap",
  "gplots",
  "corrplot",
  "RColorBrewer",
  "VennDiagram",
  "eulerr",
  "scales",
  "stringr"
))
```

### Bioconductor Packages

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "limma",
  "preprocessCore",
  "vsn",
  "clusterProfiler",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "DOSE",
  "enrichplot",
  "ReactomePA",
  "pathview",
  "pRoloc",
  "pRolocdata"
))
```

### One-liner Installation

```r
# Install all dependencies at once
install.packages(c("tidyverse", "ggplot2", "openxlsx", "pheatmap", "gplots", "corrplot", "RColorBrewer", "VennDiagram", "eulerr", "scales", "stringr"))

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "preprocessCore", "vsn", "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "DOSE", "enrichplot", "ReactomePA", "pathview", "pRoloc", "pRolocdata"))
```

## Python Dependencies (Optional)

For ConsensusPathDB SOAP client:

```bash
pip install zeep  # SOAP client for Python 3
```

## Directory Structure

```
Proteomics/
├── SKILL.md                    # Main skill definition
├── README.md                   # This file
├── DataFormats.md              # Input/output file specifications
├── NormalizationMethods.md     # Statistical normalization reference
├── VisualizationPatterns.md    # ggplot2 plotting patterns
├── PathwayAnalysisGuide.md     # Pathway enrichment methodology
├── workflows/
│   ├── Normalize.md            # Data normalization
│   ├── VolcanoPlot.md          # Volcano plot generation
│   ├── Heatmap.md              # Heatmap/PCA/correlation
│   ├── PathwayAnalysis.md      # KEGG/ConsensusPathDB
│   ├── ProteinListQuery.md     # Protein list cross-reference
│   ├── ExcelWorkup.md          # Excel report generation
│   ├── Matrisome.md            # ECM/Matrisome analysis
│   └── SaspAnalysis.md         # SASP factor analysis
├── rscripts/
│   ├── Plot_Workup_V10.R       # Visualization pipeline
│   ├── Excel_Workup_v05.R      # Excel report generation
│   ├── ConsensusPathDB_23_0411_v03.R
│   ├── toolkit.R               # Library loading
│   ├── barplots.R              # Bar plot utility
│   └── normalization/
│       ├── Step_1_Normalization.R
│       └── 2201_Label_Free_Functions.R
├── python/
│   ├── client.py               # ConsensusPathDB SOAP client
│   ├── cpdb_services.py
│   └── cpdb_services_types.py
├── data/
│   ├── MISEV2018_EV_Markers.txt
│   ├── MISEV2018_EV_Categories.txt
│   ├── Exosome_Protein_Markers.txt
│   ├── Top_10_Blood_Proteins.txt
│   ├── Apolipoproteins.txt
│   ├── Human_Core_SASP.csv
│   ├── Mouse_Core_SASP.csv
│   ├── matrisome_hs_masterlist.csv
│   └── matrisome_mm_masterlist.csv
└── tools/
    └── .gitkeep
```

## Workflows

| Workflow | Trigger | Description |
|----------|---------|-------------|
| Normalize | "normalize", "median normalization" | Apply normalization methods |
| VolcanoPlot | "volcano plot", "fold change plot" | Generate volcano plots |
| Heatmap | "heatmap", "PCA", "correlation" | PCA, heatmaps, correlation plots |
| PathwayAnalysis | "pathway", "KEGG", "enrichment" | Pathway enrichment analysis |
| ProteinListQuery | "EV markers", "check against" | Cross-reference protein lists |
| ExcelWorkup | "Excel report", "filter proteins" | Generate filtered Excel output |
| Matrisome | "matrisome", "ECM proteins" | ECM/Matrisome analysis |
| SaspAnalysis | "SASP", "senescence" | Core SASP factor analysis |

## Usage Examples

```
"Create a volcano plot for my proteomics data"
"Normalize my protein intensities using quantile normalization"
"Which of my proteins are MISEV2018 EV markers?"
"Run KEGG pathway analysis on significantly altered proteins"
"Check for SASP factors in my aging samples"
"Generate a heatmap of my top differentially expressed proteins"
```

## Reference Data

| Dataset | Description | Organisms |
|---------|-------------|-----------|
| MISEV2018_EV_Markers | Extracellular vesicle marker proteins | Human |
| Core_SASP | Senescence-associated secretory phenotype factors | Human, Mouse |
| MatrisomeDB | Extracellular matrix proteins | Human, Mouse |
| Apolipoproteins | Blood contamination markers | Human |

## License

This skill is provided as-is for proteomics data analysis.
