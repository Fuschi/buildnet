# buildnet

**buildnet** provides lightweight tools to build **correlation-based networks** from numeric matrices,  
such as microbial abundance or feature data.  

It includes functions to compute correlations, derive p-values, apply multiple thresholding strategies,  
and easily convert adjacency matrices into `igraph` network objects.

---

## 🚀 Installation

```r
# Install from GitHub
# install.packages("remotes")
remotes::install_github("Fuschi/buildnet")
```

---

## 🧩 Core functions

| Function | Description |
|-----------|-------------|
| `calculate_correlation()` | Compute pairwise correlations (Pearson, Spearman, Kendall). |
| `calculate_pvalue()` | Compute pairwise correlation p-values. |
| `threshold_absolute()` | Keep correlations above a minimum absolute value. |
| `threshold_density()` | Keep strongest correlations up to a target edge density. |
| `threshold_pvalue()` | Keep significant correlations after p-value adjustment. |
| `adjacency_to_graph()` | Convert an adjacency matrix into an `igraph` object. |
| `build_corr_net()` | High-level wrapper combining all steps into one call. |

---

## 🧠 Example

```r
library(buildnet)

# Simulated data (rows = samples, cols = features)
set.seed(1)
x <- matrix(rnorm(1000), nrow = 100, ncol = 10)
colnames(x) <- paste0("taxa_", 1:10)

# Build a correlation network using p-value thresholding
g <- build_corr_net(
  x,
  cor_method    = "pearson",
  thresh_method = "p-value",
  thresh_value  = 0.05,
  adjust        = "BH",
  output        = "graph"
)

# Inspect resulting network
g
```

---

## ⚙️ Design

- Explicit, modular steps (correlation → threshold → graph)  
- Compatible with any numeric matrix (e.g., microbial abundance, gene expression)  
- Integrates with `omicrat`

---
