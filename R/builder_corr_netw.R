# =========================
# Correlation & P-values
# =========================

#' @title Calculate correlation matrix
#' @description
#' Compute a correlation matrix from a numeric matrix where **columns are variables**
#' (features/taxa) and rows are samples. Handles zero-variance columns by returning NA.
#'
#' @param x Numeric matrix or data.frame (columns = variables to correlate).
#' @param method Correlation method: one of \code{"pearson"}, \code{"spearman"}, \code{"kendall"}.
#'
#' @return A numeric symmetric correlation matrix with column/row names from \code{x}.
#' @export
calculate_correlation <- function(x, method = c("pearson", "spearman", "kendall")) {
  method <- match.arg(method)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x)) cli::cli_abort("`x` must be a numeric matrix/data.frame.")
  if (is.null(colnames(x))) cli::cli_abort("`x` must have column names (features/taxa).")
  
  p <- ncol(x)
  cor_mat <- matrix(NA_real_, p, p)
  colnames(cor_mat) <- rownames(cor_mat) <- colnames(x)
  
  for (i in seq_len(p)) {
    for (j in i:p) {
      xi <- x[, i]
      xj <- x[, j]
      # guard against degenerate columns
      if (all(is.na(xi)) || all(is.na(xj)) ||
          sd(xi, na.rm = TRUE) == 0 || sd(xj, na.rm = TRUE) == 0) {
        cor_val <- NA_real_
      } else {
        tst <- suppressWarnings(stats::cor.test(xi, xj, method = method))
        cor_val <- unname(tst$estimate)
      }
      cor_mat[i, j] <- cor_mat[j, i] <- cor_val
    }
  }
  cor_mat
}


#' @title Calculate pairwise correlation p-values
#' @description
#' Compute a matrix of p-values from pairwise correlation tests for a numeric matrix
#' where **columns are variables** and rows are samples. Uses the same method options
#' as \code{calculate_correlation()}.
#'
#' @param x Numeric matrix or data.frame (columns = variables to correlate).
#' @param method Correlation method: \code{"pearson"}, \code{"spearman"}, or \code{"kendall"}.
#'
#' @return A symmetric matrix of p-values with the same dimnames as \code{x}.
#' @export
calculate_pvalue <- function(x, method = c("pearson", "spearman", "kendall")) {
  method <- match.arg(method)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x)) cli::cli_abort("`x` must be a numeric matrix/data.frame.")
  if (is.null(colnames(x))) cli::cli_abort("`x` must have column names (features/taxa).")
  
  p <- ncol(x)
  pval_mat <- matrix(NA_real_, p, p)
  colnames(pval_mat) <- rownames(pval_mat) <- colnames(x)
  
  for (i in seq_len(p)) {
    for (j in i:p) {
      xi <- x[, i]
      xj <- x[, j]
      if (all(is.na(xi)) || all(is.na(xj)) ||
          sd(xi, na.rm = TRUE) == 0 || sd(xj, na.rm = TRUE) == 0) {
        p_val <- NA_real_
      } else {
        tst <- suppressWarnings(stats::cor.test(xi, xj, method = method))
        p_val <- unname(tst$p.value)
      }
      pval_mat[i, j] <- pval_mat[j, i] <- p_val
    }
  }
  pval_mat
}


# =========================
# Thresholds
# =========================

#' @title Threshold correlation matrix by absolute value
#' @description
#' Build an adjacency matrix by applying a minimum absolute-correlation cutoff.
#'
#' @param cor_mat Correlation matrix (numeric, symmetric).
#' @param min_cor Minimum absolute correlation to include as an edge. In `[0, 1]`.
#'
#' @return Numeric symmetric adjacency matrix. Weights retain the (signed) correlation
#' wherever \code{|r| >= min_cor}; zeros elsewhere. Diagonal is set to 0.
#' @export
threshold_absolute <- function(cor_mat, min_cor = 0.5) {
  if (!is.matrix(cor_mat)) cli::cli_abort("`cor_mat` must be a matrix.")
  if (!is.numeric(min_cor) || min_cor < 0 || min_cor > 1)
    cli::cli_abort("`min_cor` must be a number in [0, 1].")
  
  adj <- (abs(cor_mat) >= min_cor) * cor_mat
  diag(adj) <- 0
  adj
}


#' @title Threshold correlation matrix to achieve target edge density
#' @description
#' Keep the largest absolute correlations until the desired undirected edge density.
#'
#' @param cor_mat Correlation matrix (numeric, symmetric).
#' @param density Target density in (0, 1]. For p nodes, target edges ~ \code{density * p*(p-1)/2}.
#'
#' @return Numeric symmetric adjacency matrix with retained correlations as weights.
#' @export
threshold_density <- function(cor_mat, density = 0.05) {
  if (!is.matrix(cor_mat)) cli::cli_abort("`cor_mat` must be a matrix.")
  if (!is.numeric(density) || density <= 0 || density > 1)
    cli::cli_abort("`density` must be in (0, 1].")
  
  ut <- upper.tri(cor_mat)
  cor_vals <- abs(cor_mat[ut])
  
  if (length(cor_vals) == 0 || all(is.na(cor_vals))) {
    adj <- matrix(0, nrow(cor_mat), ncol(cor_mat))
    dimnames(adj) <- dimnames(cor_mat)
    return(adj)
  }
  
  cutoff <- stats::quantile(cor_vals, probs = 1 - density, na.rm = TRUE, names = FALSE)
  adj <- (abs(cor_mat) >= cutoff) * cor_mat
  diag(adj) <- 0
  adj
}


#' @title Threshold by p-value significance
#' @description
#' Build a (0/1) adjacency mask from a p-value matrix after multiple testing correction.
#' This function decides *which* edges are significant; carry weights from your
#' correlation matrix downstream.
#'
#' @param pval_mat Square matrix of p-values (symmetric).
#' @param alpha Significance level in (0, 1].
#' @param adjust P-value adjustment: \code{"none"}, \code{"holm"}, \code{"hochberg"}, \code{"hommel"},
#'   \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}.
#'
#' @return Symmetric 0/1 numeric matrix where 1 marks a significant edge.
#' @export
threshold_pvalue <- function(pval_mat, alpha = 0.05,
                             adjust = c("none", "holm", "hochberg", "hommel",
                                        "bonferroni", "BH", "BY", "fdr")) {
  adjust <- match.arg(adjust)
  
  if (!is.matrix(pval_mat)) cli::cli_abort("`pval_mat` must be a matrix.")
  if (!is.numeric(alpha) || alpha <= 0 || alpha > 1)
    cli::cli_abort("`alpha` must be in (0, 1].")
  
  p <- nrow(pval_mat)
  if (p != ncol(pval_mat)) cli::cli_abort("`pval_mat` must be square.")
  if (is.null(colnames(pval_mat))) colnames(pval_mat) <- seq_len(p)
  if (is.null(rownames(pval_mat))) rownames(pval_mat) <- colnames(pval_mat)
  
  ut <- upper.tri(pval_mat)
  pv_ut <- pval_mat[ut]
  if (length(pv_ut) == 0 || all(is.na(pv_ut))) {
    adj <- matrix(0, p, p); dimnames(adj) <- dimnames(pval_mat); return(adj)
  }
  
  pv_adj <- stats::p.adjust(pv_ut, method = adjust)
  sig_ut <- as.numeric(pv_adj <= alpha)
  
  adj <- matrix(0, p, p)
  adj[ut] <- sig_ut
  adj <- adj + t(adj)
  diag(adj) <- 0
  dimnames(adj) <- dimnames(pval_mat)
  adj
}


# =========================
# Graph conversion
# =========================

#' @title Convert adjacency matrix to igraph
#' @description
#' Convert a symmetric adjacency matrix into an undirected \code{igraph} object.
#' If all edges are zero, return an empty graph with named vertices.
#'
#' @param adj Numeric symmetric adjacency matrix.
#' @param weighted Logical; when \code{TRUE}, carries non-zero entries as \code{weight}.
#'
#' @return An \code{igraph} object.
#' @export
adjacency_to_graph <- function(adj, weighted = TRUE) {
  if (!is.matrix(adj)) cli::cli_abort("`adj` must be a matrix.")
  if (nrow(adj) != ncol(adj)) cli::cli_abort("`adj` must be square.")
  if (is.null(colnames(adj))) colnames(adj) <- seq_len(ncol(adj))
  if (is.null(rownames(adj))) rownames(adj) <- colnames(adj)
  
  if (all(adj == 0, na.rm = TRUE)) {
    g <- igraph::make_empty_graph(n = ncol(adj), directed = FALSE)
    igraph::V(g)$name <- colnames(adj)
    return(g)
  }
  
  igraph::graph_from_adjacency_matrix(
    adjmatrix = adj,
    mode      = "undirected",
    weighted  = weighted,
    diag      = FALSE
  )
}


# =========================
# High-level wrapper
# =========================

#' @title Build a correlation network in one call
#' @description
#' Orchestrates correlation, thresholding, and graph conversion. Internally uses
#' \code{calculate_correlation()}, \code{calculate_pvalue()}, \code{threshold_*()},
#' and \code{adjacency_to_graph()}. Designed to plug into \code{set_omic()}.
#'
#' @param x Numeric matrix or data.frame with **columns = variables/features/taxa** and rows = samples.
#' @param cor_method \code{"pearson"}, \code{"spearman"}, or \code{"kendall"}.
#' @param thresh_method \code{"absolute"}, \code{"density"}, or \code{"p-value"}.
#' @param thresh_value Numeric threshold depending on \code{thresh_method}:
#'   \itemize{
#'     \item \code{"absolute"}: minimum absolute correlation in `[0, 1]`
#'     \item \code{"density"}: target edge density in (0, 1]
#'     \item \code{"p-value"}: alpha significance level in (0, 1]
#'   }
#' @param adjust Multiple-testing correction for \code{"p-value"}: \code{"none"}, \code{"holm"},
#'   \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}.
#' @param output \code{"graph"}, \code{"adjacency"}, or \code{"both"}.
#'
#' @return As requested by \code{output}: an \code{igraph}, an adjacency matrix, or a list with both.
#' @export
build_corr_net <- function(x,
                           cor_method    = c("pearson", "spearman", "kendall"),
                           thresh_method = c("absolute", "density", "p-value"),
                           thresh_value,
                           adjust        = c("none", "holm", "hochberg", "hommel",
                                             "bonferroni", "BH", "BY", "fdr"),
                           use_abs       = TRUE,
                           output        = c("graph", "adjacency", "both")) {
  cor_method    <- match.arg(cor_method)
  thresh_method <- match.arg(thresh_method)
  output        <- match.arg(output)
  adjust        <- match.arg(adjust)
  
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x)) cli::cli_abort("`x` must be a numeric matrix/data.frame.")
  if (is.null(colnames(x))) cli::cli_abort("`x` must have column names (features/taxa).")
  if (missing(thresh_value))
    cli::cli_abort("`thresh_value` is required: absolute=[0,1], density=(0,1], p-value=(0,1].")
  
  # 1) Correlation
  cor_mat <- calculate_correlation(x, method = cor_method)
  
  # 2) Threshold
  if (identical(thresh_method, "absolute")) {
    if (!is.numeric(thresh_value) || thresh_value < 0 || thresh_value > 1)
      cli::cli_abort("For {.val absolute}, `thresh_value` must be in [0,1].")
    adj <- threshold_absolute(cor_mat, min_cor = thresh_value)
    
  } else if (identical(thresh_method, "density")) {
    if (!is.numeric(thresh_value) || thresh_value <= 0 || thresh_value > 1)
      cli::cli_abort("For {.val density}, `thresh_value` must be in (0,1].")
    adj <- threshold_density(cor_mat, density = thresh_value)
    
  } else { # "p-value"
    if (!is.numeric(thresh_value) || thresh_value <= 0 || thresh_value > 1)
      cli::cli_abort("For {.val p-value}, `thresh_value` (alpha) must be in (0,1].")
    
    pval_mat <- calculate_pvalue(x, method = cor_method)
    sig_mask <- threshold_pvalue(pval_mat, alpha = thresh_value, adjust = adjust)
    adj <- (sig_mask != 0) * cor_mat
    diag(adj) <- 0
  }
  
  # 3) Output
  if (identical(output, "adjacency")) {
    return(adj)
  } else if (identical(output, "graph")) {
    return(adjacency_to_graph(adj, weighted = TRUE))
  } else {
    return(list(adjacency = adj, graph = adjacency_to_graph(adj, weighted = TRUE)))
  }
}


