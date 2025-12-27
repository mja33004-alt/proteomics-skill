# Normalization Methods Reference

This document describes the normalization methods available in the Proteomics skill and when to use each.

---

## Method Overview

| Method | Best For | Assumptions |
|--------|----------|-------------|
| **Median** | General use | Most proteins unchanged between samples |
| **Quantile** | Cross-experiment comparison | Samples have similar intensity distributions |
| **Cyclic Loess** | Technical replicates | Paired samples with systematic bias |
| **Lowess** | General smoothing | Local intensity-dependent bias |
| **VSN** | High variance data | Heteroscedasticity (variance depends on mean) |
| **Sample Loading** | TMT/labeling experiments | Known or estimable loading amounts |

---

## Median Normalization

**Description:** Scales each sample so all sample medians are equal to the global median.

**Assumptions:**
- Most proteins are unchanged between samples
- Systematic technical bias affects all proteins equally

**Implementation:**
```r
meds <- apply(iMat, 2, median, na.rm = TRUE)
MnMat <- sweep(iMat, 2, meds / mean(meds), FUN = "/")
```

**When to use:**
- Default choice for most proteomics experiments
- Robust to outliers
- Works well when ~50% or more proteins are unchanged

---

## Quantile Normalization

**Description:** Forces all samples to have identical intensity distributions by ranking and averaging.

**Assumptions:**
- Samples should have the same overall distribution
- Differences are primarily technical, not biological

**Implementation:**
```r
library(preprocessCore)
QnMat <- 2^normalize.quantiles(log2(iMat))
```

**When to use:**
- Cross-experiment comparisons
- When you expect similar global protein expression
- Microarray-style analysis

**Caution:**
- May remove true biological differences if samples are very different
- Not recommended when comparing very different conditions

---

## Cyclic Loess Normalization

**Description:** Iteratively applies local regression (loess) between pairs of samples to remove intensity-dependent bias.

**Assumptions:**
- Bias between samples is intensity-dependent
- Most proteins are unchanged

**Implementation:**
```r
library(limma)
CLnMat <- 2^normalizeCyclicLoess(data.frame(log2(iMat)), method = "pairs")
```

**When to use:**
- When MA plots show intensity-dependent bias (banana shape)
- Paired experimental designs
- Technical replicates

---

## Lowess Normalization

**Description:** Local weighted scatterplot smoothing for removing intensity-dependent effects.

**Implementation:**
```r
library(limma)
LFnMat <- 2^normalizeCyclicLoess(data.frame(log2(iMat)))
```

**When to use:**
- Similar to Cyclic Loess but less aggressive
- When loess fits show consistent patterns

---

## VSN (Variance Stabilizing Normalization)

**Description:** Uses arsinh transformation to stabilize variance across the intensity range.

**Assumptions:**
- Variance is heteroscedastic (changes with mean)
- Data follows a specific noise model

**Implementation:**
```r
library(vsn)
VSNnMat <- 2^justvsn(iMat)
```

**When to use:**
- When variance increases with intensity (common in proteomics)
- For subsequent statistical tests that assume homoscedasticity
- Use `meanSdPlot(justvsn(iMat))` to verify stabilization

---

## Sample Loading Normalization

**Description:** Scales based on total protein loading, often used for TMT/iTRAQ labeled experiments.

**Assumptions:**
- Total protein per sample is known or can be estimated
- Differences in total signal are technical, not biological

**Implementation:**
```r
# Calculate target (global average)
target <- mean(colSums(iMat, na.rm = TRUE))

# Calculate normalization factors
norm_facs <- target / colSums(iMat, na.rm = TRUE)

# Apply
SLnMat <- sweep(iMat, 2, norm_facs, FUN = "*")
```

**For multi-batch TMT:**
```r
# Separate by experiment/batch
exp1 <- select(dMat, ends_with(".1"))
exp2 <- select(dMat, ends_with(".2"))

# Normalize each batch to same target
target <- mean(c(colSums(exp1, na.rm = TRUE), colSums(exp2, na.rm = TRUE)))

norm_facs <- target / colSums(exp1, na.rm = TRUE)
exp1 <- sweep(exp1, 2, norm_facs, FUN = "*")

# Repeat for each batch, then combine
```

---

## Choosing a Method

**Decision Tree:**

1. **Is this TMT/labeling data?**
   - Yes → Start with Sample Loading, then consider Median or Quantile
   - No → Continue

2. **Do MA plots show intensity-dependent bias?**
   - Yes → Use Cyclic Loess or Lowess
   - No → Continue

3. **Do you need to compare across experiments?**
   - Yes → Consider Quantile
   - No → Continue

4. **Is variance heteroscedastic?**
   - Yes → Consider VSN
   - No → Use Median (default)

---

## Quality Control

After normalization, verify with:

1. **Boxplots:** Log2 intensities should be aligned
   ```r
   boxplot(log2(nMat), las = 2)
   ```

2. **MA plots:** Should center around y = 0
   ```r
   # Compare sample 1 vs sample 2
   M <- log2(nMat[,1]) - log2(nMat[,2])
   A <- 0.5 * (log2(nMat[,1]) + log2(nMat[,2]))
   plot(A, M)
   abline(h = 0, col = "red")
   ```

3. **No NA/Inf introduced:**
   ```r
   sum(is.na(nMat))
   sum(is.infinite(nMat))
   ```

4. **Sample names preserved:**
   ```r
   identical(colnames(nMat), colnames(iMat))
   ```

---

## Required R Packages

```r
# CRAN
install.packages("tidyverse")

# Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("limma", "preprocessCore", "vsn"))
```
