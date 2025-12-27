# Visualization Patterns Reference

This document contains ggplot2 and base R plotting patterns extracted from the Proteomics skill R scripts.

---

## Volcano Plot

```r
library(ggplot2)
library(stringr)

# Parameters
myQval <- 0.01
myFC <- 0.58
mycolors <- c("Blue" = "Blue", "Red" = "Red", "Gray" = "Gray")

volcano <- function(comparison_name, df) {
  # Filter to comparison
  plot_df <- df[df$Comparison..group1.group2. == comparison_name, ]
  plot_df$minuslogqval <- -1 * log10(plot_df$Qvalue)

  # Assign colors
  plot_df$Color <- ifelse(plot_df$AVG.Log2.Ratio >= myFC & plot_df$Qvalue < myQval, "Red", "Gray")
  plot_df$Color <- ifelse(plot_df$AVG.Log2.Ratio <= -myFC & plot_df$Qvalue < myQval, "Blue", plot_df$Color)

  # Create plot
  ggplot(plot_df, aes(x = AVG.Log2.Ratio, y = minuslogqval, col = Color)) +
    geom_point() +
    geom_vline(xintercept = c(-myFC, myFC), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(myQval), linetype = "dashed", color = "black") +
    scale_color_manual(values = mycolors) +
    ylab("-Log10(q-value)") +
    scale_x_continuous(name = "Log2(fold change)", limits = c(-5, 5)) +
    ylim(c(0, 100)) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 20, color = "black"),
      axis.text = element_text(size = 18, color = "black"),
      legend.position = "none"
    )
}
```

---

## PCA Plot

```r
library(ggplot2)

pca_plot <- function(pcRes, condition_df) {
  pcSum <- summary(pcRes)

  PC1label <- paste0("PC1, ", round(100 * pcSum$importance["Proportion of Variance", "PC1"], 1), "% of variance")
  PC2label <- paste0("PC2, ", round(100 * pcSum$importance["Proportion of Variance", "PC2"], 1), "% of variance")

  pcPlotFrame <- data.frame(
    treatment = condition_df$Condition,
    sample = rownames(pcRes$x),
    pcRes$x[, 1:5],
    color = condition_df$Color
  )

  ggplot(pcPlotFrame, aes(PC1, PC2, color = color, shape = treatment)) +
    geom_point(size = 3) +
    scale_x_continuous(name = PC1label) +
    scale_y_continuous(name = PC2label) +
    scale_color_identity() +
    stat_ellipse(aes(color = color)) +
    theme_bw(base_size = 14) +
    theme(legend.position = "right")
}
```

---

## Heatmap

```r
library(gplots)

heatmap_plot <- function(mat, color_sidebar) {
  myRamp <- colorRampPalette(colors = c("#0571b0", "#f7f7f7", "#ca0020"))

  heatmap.2(
    t(scale(t(log10(mat)))),
    col = myRamp,
    trace = "none",
    labRow = FALSE,
    ColSideColors = color_sidebar,
    key = TRUE,
    density.info = "none"
  )
}
```

---

## Correlation Plot

```r
library(corrplot)

correlation_plot <- function(mat) {
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  M <- cor(mat, use = "complete.obs")

  corrplot(
    M,
    method = "color",
    col = col(200),
    tl.col = "black",
    tl.srt = 45,
    addCoef.col = "black",
    number.cex = 0.7
  )

  return(M)
}
```

---

## Violin Plot

```r
library(ggplot2)

violin_plot <- function(df, group_col, value_col, fill_colors) {
  ggplot(df, aes(x = factor(.data[[group_col]]), y = .data[[value_col]],
                  fill = .data[[group_col]], color = .data[[group_col]])) +
    geom_violin(trim = TRUE, scale = "width") +
    stat_summary(fun = mean, geom = "point", shape = 23, fill = "black", size = 3) +
    geom_jitter(shape = 16, position = position_jitter(0.1), color = alpha("black", 0.2)) +
    scale_fill_manual(values = alpha(fill_colors, 0.7)) +
    scale_color_manual(values = fill_colors) +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = 18, color = "black"),
      axis.text.y = element_text(size = 16, color = "black"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "bottom"
    )
}
```

---

## Matrisome Histogram

```r
library(ggplot2)

matrisome_histogram <- function(df) {
  category_order <- c("Collagens", "ECM Glycoproteins", "Proteoglycans",
                      "ECM-affiliated Proteins", "ECM Regulators", "Secreted Factors")

  df$Category <- factor(df$Category, levels = category_order)

  df %>%
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
}
```

---

## LOPIT Plot

```r
library(ggplot2)
library(pRoloc)
library(pRolocdata)

lopit_plot <- function(df, reference_data, comparison_name) {
  # Filter to significant proteins
  myCan <- df[!is.na(df$Absolute.AVG.Log2.Ratio) & df$Color != "Grey", ]

  ggplot(reference_data, aes(x = Dimension.1, y = Dimension.2)) +
    geom_point(alpha = 0.1) +
    geom_point(
      data = myCan,
      aes(x = Dimension.1, y = Dimension.2, size = Absolute.AVG.Log2.Ratio),
      color = myCan$Color,
      alpha = 0.3
    ) +
    labs(
      title = comparison_name,
      x = "t-SNE Dim. 1",
      y = "t-SNE Dim. 2"
    ) +
    guides(size = guide_legend(title = "|log2(FC)| Ratio")) +
    theme_bw(base_size = 12)
}
```

---

## Pathway Dotplot

```r
library(ggplot2)

pathway_dotplot <- function(enrichment_result, top_n = 20) {
  # For clusterProfiler results
  dotplot(enrichment_result, showCategory = top_n) +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10)
    )
}

# Custom dotplot for other enrichment results
custom_pathway_dotplot <- function(df) {
  ggplot(df, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(x = "Gene Ratio", y = "") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    )
}
```

---

## Color Palettes

```r
# Blue-White-Red diverging (for fold change)
myRamp <- colorRampPalette(colors = c("#0571b0", "#f7f7f7", "#ca0020"))

# Red-Yellow-Blue diverging
col_ramp <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

# Categorical (paired)
library(RColorBrewer)
paired_colors <- brewer.pal(12, "Paired")

# Matrisome colors
matrisome_colors <- c(
  alpha('#FD6467', 1), alpha('#FD6467', 0.5), alpha('#FD6467', 0.2),
  alpha('#5B1A18', 1), alpha('#5B1A18', 0.5), alpha('#5B1A18', 0.2)
)
```

---

## Saving Plots

```r
# TIFF (publication quality)
ggsave("plot.tiff", dpi = 300, width = 5, height = 5, units = "in")

# PDF (vector)
ggsave("plot.pdf", width = 5, height = 5, units = "in")

# PNG (web/presentation)
ggsave("plot.png", dpi = 300, width = 5, height = 5, units = "in")

# Base R TIFF
tiff("plot.tiff", res = 300, height = 5, width = 5, units = "in")
# ... plotting code ...
dev.off()
```
