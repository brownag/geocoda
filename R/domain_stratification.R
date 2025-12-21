#' Identify Spatial Strata from PCA Clustering
#'
#' Automatically detect spatial domains or strata based on clustering of Principal
#' Component scores from ILR values. This suggests natural compositional boundaries
#' within the study area, useful for domain stratification when non-stationarity is detected.
#'
#' @param data A data frame or sf object with columns `x`, `y` for spatial coordinates
#'   and columns `ilr1`, `ilr2`, etc., with ILR values.
#' @param n_strata Numeric. Number of strata to identify (default 2). Can be a vector
#'   of candidate values for automatic selection via silhouette analysis.
#' @param method Character string specifying clustering approach: `"kmeans"` (default)
#'   uses K-means clustering; `"hierarchical"` uses hierarchical clustering with
#'   automatic cut height detection.
#' @param plot Logical. If `TRUE` (default), visualize strata spatially and via PCA.
#'
#' @return A list containing:
#'   - `strata`: Factor vector indicating stratum membership for each observation
#'   - `n_strata`: Number of identified strata
#'   - `cluster_centers`: Cluster centers in PCA space
#'   - `silhouette_widths`: Silhouette widths for each observation (quality of assignment)
#'   - `recommendation`: Character string with interpretation and suggested actions
#'   - `summary`: Data frame with stratum-wise statistics
#'
#' @details
#' **Methodology**:
#' 1. Perform PCA on ILR values
#' 2. Apply K-means clustering on PC1-PC2 scores
#' 3. Assign each spatial location to a cluster/stratum
#' 4. Optionally optimize cluster count via silhouette analysis
#'
#' **Interpretation**:
#' - Silhouette width: -1 (worst assignment) to +1 (best assignment)
#' - Mean silhouette > 0.5 indicates strong clustering
#' - Strata with low silhouette widths may need adjustment (merge or increase n_strata)
#'
#' **Usage in Workflow**:
#' After detecting non-stationarity via `gc_assess_stationarity()`, use this function
#' to suggest how to partition the domain. Then apply the modeling workflow independently
#' to each stratum for improved modeling.
#'
#' @examples
#' \dontrun{
#' # Create sample data with distinct compositional domains
#' n_points <- 100
#' set.seed(42)
#' data <- data.frame(
#'   x = runif(n_points, 0, 100),
#'   y = runif(n_points, 0, 100),
#'   ilr1 = c(rnorm(50, 0.3, 0.6), rnorm(50, -0.3, 0.6)),
#'   ilr2 = c(rnorm(50, 0.1, 0.5), rnorm(50, 0.2, 0.5))
#' )
#'
#' # Identify spatial strata
#' strata_result <- gc_identify_strata(data, n_strata = 2, plot = TRUE)
#' print(strata_result$recommendation)
#' print(strata_result$summary)
#'
#' # Use stratum information for domain-specific modeling
#' stratum_1 <- data[strata_result$strata == 1, ]
#' stratum_2 <- data[strata_result$strata == 2, ]
#' # ... apply gc_ilr_model(), gc_sim_composition() separately
#' }
#'
#' @importFrom stats prcomp kmeans dist hclust cutree
#' @importFrom cluster silhouette
#' @export
gc_identify_strata <- function(data,
                               n_strata = 2,
                               method = "kmeans",
                               plot = TRUE) {
  if (methods::is(data, "data.frame")) {
    if (!all(c("x", "y") %in% colnames(data))) {
      stop("data must contain columns 'x' and 'y'")
    }
    data_df <- data
  } else if (methods::is(data, "sf")) {
    coords <- sf::st_coordinates(data)
    data_dropped <- sf::st_drop_geometry(data)
    data_df <- data.frame(
      x = coords[, 1],
      y = coords[, 2],
      data_dropped
    )
  } else {
    stop("data must be a data frame or sf object")
  }

  ilr_cols <- grep("^ilr\\d+$", colnames(data_df), value = TRUE)
  if (length(ilr_cols) < 2) {
    stop("data must contain at least 2 ILR columns (ilr1, ilr2, ...)")
  }

  method <- match.arg(method, c("kmeans", "hierarchical"))

  if (length(n_strata) > 1) {
    # Auto-select optimal number of strata via silhouette analysis
    result <- .identify_strata_auto(data_df, ilr_cols, n_strata, method, plot)
  } else {
    # Use specified number of strata
    if (method == "kmeans") {
      result <- .identify_strata_kmeans(data_df, ilr_cols, n_strata, plot)
    } else {
      result <- .identify_strata_hierarchical(data_df, ilr_cols, n_strata, plot)
    }
  }

  return(result)
}

.identify_strata_kmeans <- function(data_df, ilr_cols, n_strata, plot) {
  ilr_matrix <- as.matrix(data_df[, ilr_cols])

  # PCA on ILR values
  pca_result <- stats::prcomp(ilr_matrix, scale = TRUE)
  pca_scores <- pca_result$x[, 1:2]

  # K-means clustering
  km_result <- stats::kmeans(pca_scores, centers = n_strata, nstart = 10)

  # Compute silhouette widths
  sil <- cluster::silhouette(km_result$cluster, stats::dist(pca_scores))
  sil_widths <- sil[, "sil_width"]

  mean_sil <- mean(sil_widths)
  quality <- if (mean_sil > 0.5) {
    "strong"
  } else if (mean_sil > 0.25) {
    "moderate"
  } else {
    "weak"
  }

  recommendation <- sprintf(
    "Clustering quality: %s (mean silhouette = %.3f). %d strata identified.",
    quality, mean_sil,
    n_strata
  )

  if (quality == "weak") {
    recommendation <- paste(
      recommendation,
      "Consider merging strata or checking for bimodal distributions."
    )
  }

  # Spatial visualization
  if (plot) {
    graphics::par(mfrow = c(1, 2))

    # PCA space
    graphics::plot(
      pca_scores,
      col = km_result$cluster,
      main = "Strata in PCA Space",
      xlab = sprintf("PC1 (%.1f%%)", pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100),
      ylab = sprintf("PC2 (%.1f%%)", pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100),
      pch = 16,
      cex = 0.8
    )
    graphics::points(km_result$centers, col = 1:n_strata, pch = 3, cex = 2, lwd = 2)

    # Spatial layout
    colors <- c("red", "blue", "green", "orange", "purple", "brown", "pink", "gray")[1:n_strata]
    graphics::plot(
      data_df$x, data_df$y,
      col = colors[km_result$cluster],
      main = "Strata in Space",
      xlab = "X Coordinate",
      ylab = "Y Coordinate",
      pch = 16,
      cex = 0.8
    )

    graphics::par(mfrow = c(1, 1))
  }

  # Summary statistics
  summary_df <- data.frame(
    Stratum = 1:n_strata,
    N_Observations = as.numeric(table(km_result$cluster)),
    Mean_Silhouette = by(sil_widths, km_result$cluster, mean),
    Quality = ifelse(by(sil_widths, km_result$cluster, mean) > 0.5, "Good", "Fair")
  )
  rownames(summary_df) <- NULL

  list(
    strata = factor(km_result$cluster),
    n_strata = n_strata,
    cluster_centers = km_result$centers,
    silhouette_widths = sil_widths,
    pca_loadings = pca_result$rotation,
    pca_scores = pca_result$x,
    recommendation = recommendation,
    summary = summary_df
  )
}

.identify_strata_hierarchical <- function(data_df, ilr_cols, n_strata, plot) {
  ilr_matrix <- as.matrix(data_df[, ilr_cols])

  # PCA on ILR values
  pca_result <- stats::prcomp(ilr_matrix, scale = TRUE)
  pca_scores <- pca_result$x[, 1:2]

  # Hierarchical clustering
  hc_result <- stats::hclust(stats::dist(pca_scores), method = "ward.D2")

  # Cut dendrogram at specified height
  strata_assignment <- stats::cutree(hc_result, k = n_strata)

  # Compute silhouette widths
  sil <- cluster::silhouette(strata_assignment, stats::dist(pca_scores))
  sil_widths <- sil[, "sil_width"]

  mean_sil <- mean(sil_widths)
  quality <- if (mean_sil > 0.5) {
    "strong"
  } else if (mean_sil > 0.25) {
    "moderate"
  } else {
    "weak"
  }

  recommendation <- sprintf(
    "Hierarchical clustering quality: %s (mean silhouette = %.3f). %d strata identified.",
    quality, mean_sil,
    n_strata
  )

  if (quality == "weak") {
    recommendation <- paste(
      recommendation,
      "Consider alternative cut height or checking for continuous gradients."
    )
  }

  # Spatial visualization
  if (plot) {
    graphics::par(mfrow = c(1, 2))

    # PCA space
    graphics::plot(
      pca_scores,
      col = strata_assignment,
      main = "Strata in PCA Space (Hierarchical)",
      xlab = sprintf("PC1 (%.1f%%)", pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100),
      ylab = sprintf("PC2 (%.1f%%)", pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100),
      pch = 16,
      cex = 0.8
    )

    # Spatial layout
    colors <- c("red", "blue", "green", "orange", "purple", "brown", "pink", "gray")[1:n_strata]
    graphics::plot(
      data_df$x, data_df$y,
      col = colors[strata_assignment],
      main = "Strata in Space",
      xlab = "X Coordinate",
      ylab = "Y Coordinate",
      pch = 16,
      cex = 0.8
    )

    graphics::par(mfrow = c(1, 1))
  }

  # Summary statistics
  summary_df <- data.frame(
    Stratum = 1:n_strata,
    N_Observations = as.numeric(table(strata_assignment)),
    Mean_Silhouette = by(sil_widths, strata_assignment, mean),
    Quality = ifelse(by(sil_widths, strata_assignment, mean) > 0.5, "Good", "Fair")
  )
  rownames(summary_df) <- NULL

  list(
    strata = factor(strata_assignment),
    n_strata = n_strata,
    silhouette_widths = sil_widths,
    pca_loadings = pca_result$rotation,
    pca_scores = pca_result$x,
    recommendation = recommendation,
    summary = summary_df
  )
}

.identify_strata_auto <- function(data_df, ilr_cols, n_strata_candidates, method, plot) {
  ilr_matrix <- as.matrix(data_df[, ilr_cols])
  pca_result <- stats::prcomp(ilr_matrix, scale = TRUE)
  pca_scores <- pca_result$x[, 1:2]

  best_sil <- -Inf
  best_n <- n_strata_candidates[1]
  sil_results <- list()

  for (n in n_strata_candidates) {
    if (method == "kmeans") {
      km_result <- stats::kmeans(pca_scores, centers = n, nstart = 10)
      assignment <- km_result$cluster
    } else {
      hc_result <- stats::hclust(stats::dist(pca_scores), method = "ward.D2")
      assignment <- stats::cutree(hc_result, k = n)
    }

    sil <- cluster::silhouette(assignment, stats::dist(pca_scores))
    mean_sil <- mean(sil[, "sil_width"])
    sil_results[[as.character(n)]] <- mean_sil

    if (mean_sil > best_sil) {
      best_sil <- mean_sil
      best_n <- n
    }
  }

  # Use best result
  if (method == "kmeans") {
    km_result <- stats::kmeans(pca_scores, centers = best_n, nstart = 10)
    assignment <- km_result$cluster
  } else {
    hc_result <- stats::hclust(stats::dist(pca_scores), method = "ward.D2")
    assignment <- stats::cutree(hc_result, k = best_n)
  }

  sil <- cluster::silhouette(assignment, stats::dist(pca_scores))
  sil_widths <- sil[, "sil_width"]

  recommendation <- sprintf(
    "Optimal strata count: %d (mean silhouette = %.3f) selected from candidates %s.",
    best_n, best_sil,
    paste(n_strata_candidates, collapse = ", ")
  )

  if (plot) {
    graphics::par(mfrow = c(1, 3))

    # Silhouette scores
    graphics::plot(
      as.numeric(names(sil_results)),
      unlist(sil_results),
      type = "b",
      main = "Silhouette vs Number of Strata",
      xlab = "Number of Strata",
      ylab = "Mean Silhouette Width",
      ylim = c(0, 1)
    )
    graphics::abline(v = best_n, col = "red", lty = 2)

    # PCA space
    colors <- c("red", "blue", "green", "orange", "purple", "brown", "pink", "gray")[1:best_n]
    graphics::plot(
      pca_scores,
      col = colors[assignment],
      main = "Optimal Strata in PCA Space",
      xlab = sprintf("PC1 (%.1f%%)", pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100),
      ylab = sprintf("PC2 (%.1f%%)", pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100),
      pch = 16,
      cex = 0.8
    )

    # Spatial layout
    graphics::plot(
      data_df$x, data_df$y,
      col = colors[assignment],
      main = "Optimal Strata in Space",
      xlab = "X Coordinate",
      ylab = "Y Coordinate",
      pch = 16,
      cex = 0.8
    )

    graphics::par(mfrow = c(1, 1))
  }

  # Summary statistics
  summary_df <- data.frame(
    Stratum = 1:best_n,
    N_Observations = as.numeric(table(assignment)),
    Mean_Silhouette = by(sil_widths, assignment, mean),
    Quality = ifelse(by(sil_widths, assignment, mean) > 0.5, "Good", "Fair")
  )
  rownames(summary_df) <- NULL

  list(
    strata = factor(assignment),
    n_strata = best_n,
    silhouette_widths = sil_widths,
    pca_loadings = pca_result$rotation,
    pca_scores = pca_result$x,
    recommendation = recommendation,
    summary = summary_df
  )
}
