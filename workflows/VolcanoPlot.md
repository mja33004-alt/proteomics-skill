# VolcanoPlot Workflow

**Trigger:** "volcano plot", "create volcano", "visualize fold change"

---

## Overview

Generate volcano plots to visualize differential protein expression, showing statistical significance (q-value) vs. fold change.

## Required Input

A comparison data file with columns:
- `Genes` or `PG.Genes` - Gene symbols
- `AVG.Log2.Ratio` or `Log2FC` - Log2 fold change
- `Qvalue` - Statistical significance (adjusted p-value)
- `Comparison..group1.group2.` - Comparison identifier

---

## Execution Steps

### Step 1: Gather Parameters

Ask user for:
1. **Comparison data file** - CSV with fold change and q-value data
2. **Q-value threshold** - Default: 0.01
3. **Fold-change threshold** - Default: 0.58 (1.5x)
4. **Which comparison(s)** - Specific or all comparisons

### Step 2: Option A - Use Full Pipeline Script

For complete analysis with PCA, heatmaps, and volcano plots:

```bash
cd [PROJECT_DIR]
Rscript ~/.claude/Skills/Proteomics/rscripts/Plot_Workup_V10.R
```

**Required modifications in script:**
- Line 11: `setwd("~/Desktop/[YOUR_PROJECT]/")`
- Line 14: `batch <- "[YOUR_BATCH]"`
- Line 17: `organism <- "human"` or `"mouse"`
- Line 33: `mypattern <- "[YOUR_PATTERN]"`
- Lines 36-38: Set q-value filter flags
- Line 41: `myFC <- 0.58`

### Step 3: Option B - Generate Custom Code

For standalone volcano plot generation:

```r
library(ggplot2)
library(stringr)

# Parameters
myQval <- 0.01
myFC <- 0.58
mycolors <- c("Blue" = "Blue", "Red" = "Red", "Gray" = "Gray")

# Read data
all_comparisons <- read.csv("data/candidates.csv", stringsAsFactors = FALSE)
all_comparisons$minuslogqval <- -1 * log10(all_comparisons$Qvalue)

# Volcano plot function
volcano <- function(comparison_name) {
  df <- all_comparisons[all_comparisons$Comparison..group1.group2. == comparison_name, ]

  # Assign colors based on significance
  df$Color <- ifelse(df$AVG.Log2.Ratio >= myFC & df$Qvalue < myQval, "Red", "Gray")
  df$Color <- ifelse(df$AVG.Log2.Ratio <= -myFC & df$Qvalue < myQval, "Blue", df$Color)

  # Count up/down regulated
  down <- sum(df$Color == "Blue")
  up <- sum(df$Color == "Red")

  # Create plot
  p <- ggplot(df, aes(x = AVG.Log2.Ratio, y = minuslogqval, col = Color)) +
    geom_point() +
    geom_vline(xintercept = c(-myFC, myFC), linetype = "dashed") +
    geom_hline(yintercept = -log10(myQval), linetype = "dashed") +
    scale_color_manual(values = mycolors) +
    ylab("-Log10(q-value)") +
    scale_x_continuous(name = "Log2(fold change)", limits = c(-5, 5)) +
    ylim(c(0, 100)) +
    theme_classic() +
    theme(axis.title = element_text(size = 20, color = "black"),
          axis.text = element_text(size = 18, color = "black"),
          legend.position = "none")

  # Save plot
  filename <- paste0("output/", gsub(" / ", "_vs_", comparison_name), "_volcano.tiff")
  ggsave(filename, plot = p, dpi = 300, width = 5, height = 5, units = "in")

  cat(paste0(comparison_name, ": ", up, " up, ", down, " down\n"))
}

# Apply to all comparisons
sapply(unique(all_comparisons$Comparison..group1.group2.), volcano)
```

### Step 4: Report Results

For each comparison, report:
- Number of up-regulated proteins (Red)
- Number of down-regulated proteins (Blue)
- File saved location

---

## Customization Options

| Option | How to Modify |
|--------|---------------|
| X-axis range | `limits = c(-X, X)` in `scale_x_continuous` |
| Y-axis range | `ylim(c(0, Y))` |
| Add gene labels | Add `geom_text_repel(aes(label = Genes), data = subset(...))` |
| Change colors | Modify `mycolors` vector |
| Adjust thresholds | Change `myFC` and `myQval` values |
| Output format | Change device in `ggsave()` to "pdf", "png", etc. |

---

## Example Output

```
User: "Create volcano plots for my kidney comparison data"

Claude: Running the **VolcanoPlot** workflow from the **Proteomics** skill...

        Generated 3 volcano plots:

        1. Cortex_vs_Medulla_volcano.tiff
           - 245 Up-regulated (Red)
           - 189 Down-regulated (Blue)
           - Thresholds: q < 0.01, |log2FC| > 0.58

        2. Diseased_vs_Control_volcano.tiff
           - 523 Up-regulated
           - 412 Down-regulated

        All files saved to: ~/Desktop/RCCN_Kidney/output/
```
