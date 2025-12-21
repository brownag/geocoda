#' Assess Stationarity of Compositional Spatial Data
#'
#' Evaluate whether the spatial correlation structure is constant (stationary) across
#' the study area using PCA biplot visualization and spatial clustering diagnostics.
#' Non-stationarity suggests the domain should be stratified or that a global Linear
#' Model of Coregionalization may be inappropriate.
#'
#' @param data A data frame or sf object with columns `x`, `y` for spatial coordinates
#'   and columns `ilr1`, `ilr2`, etc., with ILR values.
#' @param method Character string specifying diagnostic approach: `"biplot"` (default)
#'   uses PCA and spatial coloring; `"local"` computes local covariance within
#'   spatial windows.
#' @param plot Logical. If `TRUE` (default), display biplot colored by spatial location
#'   to visually inspect for clusters or gradients.
#' @param window_size Numeric. Size of local window for `method = "local"` (default `NULL`
#'   uses automatic sizing based on spatial extent).
#'
#' @return A list containing:
#'   - `stationary`: Logical indicating whether global stationarity assumption is supported
#'   - `method_used`: Character string indicating which diagnostic was used
#'   - `pca_loadings`: PCA loadings matrix (for biplot method)
#'   - `pca_scores`: Principal component scores
#'   - `spatial_variance`: Estimated variance of spatial means across local windows
#'   - `recommendation`: Character string with interpretation and recommended action
#'   - `summary`: Data frame with diagnostic statistics
#'
#' @details
#' **Biplot Method** (default):
#' Performs PCA on ILR values and creates a biplot with points colored by spatial
#' location (x, y). If the biplot shows distinct clusters that spatially correspond
#' to different regions of the study area, this suggests non-stationarity (e.g.,
#' distinct geological units with different covariance structures).
#'
#' **Local Method**:
#' Divides the spatial domain into regular windows and computes covariance separately
#' in each window. Large variation in covariance between windows indicates non-stationarity.
#'
#' **Interpretation**:
#' - Stationary: Points are randomly scattered in biplot, no spatial clustering
#' - Non-Stationary: Distinct biplot clusters correspond to distinct spatial regions
#' - Action: If non-stationary, stratify the domain and apply the kriging workflow
#'   independently to each stratum.
#'
#' @examples
#' \dontrun{
#' # Create sample data with spatial coordinates and ILR values
#' n_points <- 100
#' set.seed(42)
#' data <- data.frame(
#'   x = runif(n_points, 0, 100),
#'   y = runif(n_points, 0, 100),
#'   ilr1 = rnorm(n_points, mean = 0.5, sd = 0.8),
#'   ilr2 = rnorm(n_points, mean = -0.2, sd = 0.6)
#' )
#'
#' # Assess stationarity
#' stationarity <- gc_assess_stationarity(data, method = "biplot", plot = TRUE)
#' print(stationarity$recommendation)
#'
#' # Or use local window method
#' stationarity_local <- gc_assess_stationarity(data, method = "local")
#' print(stationarity_local$summary)
#' }
#'
#' @importFrom stats prcomp cov
#' @importFrom methods is
#' @export
gc_assess_stationarity <- function(data,
                                   method = "biplot",
                                   plot = TRUE,
                                   window_size = NULL) {
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

  method <- match.arg(method, c("biplot", "local"))

  if (method == "biplot") {
    result <- .assess_stationarity_biplot(data_df, ilr_cols, plot)
  } else {
    result <- .assess_stationarity_local(data_df, ilr_cols, window_size)
  }

  return(result)
}

.assess_stationarity_biplot <- function(data_df, ilr_cols, plot) {
  ilr_matrix <- as.matrix(data_df[, ilr_cols])

  pca_result <- stats::prcomp(ilr_matrix, scale = TRUE)

  pc1_scores <- pca_result$x[, 1]
  pc2_scores <- pca_result$x[, 2]

  spatial_pc1_var <- stats::var(tapply(pc1_scores, cut(data_df$x, breaks = 5), mean), na.rm = TRUE)
  spatial_pc2_var <- stats::var(tapply(pc2_scores, cut(data_df$y, breaks = 5), mean), na.rm = TRUE)

  spatial_variance <- (spatial_pc1_var + spatial_pc2_var) / 2

  loadings_scale <- sqrt(colSums(pca_result$rotation^2))
  is_stationary <- spatial_variance < 0.1

  if (plot) {
    # Create color palette from blue to red scaled by X coordinate
    color_ramp_func <- grDevices::colorRamp(c("blue", "red"))
    normalized_x <- (data_df$x - min(data_df$x)) / (max(data_df$x) - min(data_df$x))
    color_matrix <- color_ramp_func(normalized_x)
    colors <- grDevices::rgb(color_matrix[, 1], color_matrix[, 2], color_matrix[, 3], maxColorValue = 255)
    
    graphics::plot(
      pc1_scores, pc2_scores,
      col = colors,
      main = "PCA Biplot: ILR Space (colored by X coordinate)",
      xlab = sprintf("PC1 (%.1f%%)", pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100),
      ylab = sprintf("PC2 (%.1f%%)", pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100),
      pch = 16,
      cex = 0.8
    )
    graphics::abline(h = 0, v = 0, lty = 2, col = "lightgray")
  }

  recommendation <- if (is_stationary) {
    "Stationarity supported: PCA scores show random spatial distribution. Global LMC model is appropriate."
  } else {
    "Potential non-stationarity detected: PCA scores cluster spatially. Consider stratifying domain into geological/soil units."
  }

  list(
    stationary = is_stationary,
    method_used = "biplot",
    pca_loadings = pca_result$rotation,
    pca_scores = pca_result$x,
    spatial_variance = spatial_variance,
    recommendation = recommendation,
    summary = data.frame(
      Method = "Biplot (PCA)",
      SpatialVariance = round(spatial_variance, 4),
      Stationary = is_stationary,
      Interpretation = recommendation,
      stringsAsFactors = FALSE
    )
  )
}

.assess_stationarity_local <- function(data_df, ilr_cols, window_size) {
  if (is.null(window_size)) {
    x_range <- diff(range(data_df$x))
    y_range <- diff(range(data_df$y))
    window_size <- (x_range + y_range) / 2 / 4
  }

  x_breaks <- seq(min(data_df$x), max(data_df$x), by = window_size)
  y_breaks <- seq(min(data_df$y), max(data_df$y), by = window_size)

  window_covs <- list()
  window_means <- list()
  window_n <- 0

  for (i in seq_len(length(x_breaks) - 1)) {
    for (j in seq_len(length(y_breaks) - 1)) {
      x_in <- data_df$x >= x_breaks[i] & data_df$x < x_breaks[i + 1]
      y_in <- data_df$y >= y_breaks[j] & data_df$y < y_breaks[j + 1]
      window_idx <- x_in & y_in

      if (sum(window_idx) >= 3) {
        window_n <- window_n + 1
        window_data <- as.matrix(data_df[window_idx, ilr_cols])

        window_cov <- stats::cov(window_data)
        window_mean <- colMeans(window_data)

        window_covs[[window_n]] <- window_cov
        window_means[[window_n]] <- window_mean
      }
    }
  }

  if (window_n < 2) {
    warning("Not enough data in windows. Returning inconclusive result.")
    return(list(
      stationary = NA,
      method_used = "local",
      spatial_variance = NA,
      recommendation = "Insufficient data for local window analysis",
      summary = data.frame(
        Method = "Local Windows",
        WindowsAnalyzed = window_n,
        Stationary = NA,
        stringsAsFactors = FALSE
      )
    ))
  }

  global_cov <- stats::cov(as.matrix(data_df[, ilr_cols]))

  covariance_diffs <- sapply(window_covs, function(w) {
    sqrt(sum((w - global_cov)^2))
  })

  mean_cov_diff <- mean(covariance_diffs)
  global_cov_scale <- sqrt(sum(global_cov^2))

  normalized_diff <- mean_cov_diff / global_cov_scale

  is_stationary <- normalized_diff < 0.15

  recommendation <- if (is_stationary) {
    sprintf(
      "Stationarity supported: Local covariances vary minimally (%.3f) from global. Global LMC model is appropriate.",
      normalized_diff
    )
  } else {
    sprintf(
      "Potential non-stationarity detected: Local covariances vary substantially (%.3f). Domain stratification recommended.",
      normalized_diff
    )
  }

  spatial_variance <- stats::var(unlist(window_means))

  list(
    stationary = is_stationary,
    method_used = "local",
    spatial_variance = spatial_variance,
    covariance_variation = mean_cov_diff,
    recommendation = recommendation,
    summary = data.frame(
      Method = "Local Windows",
      WindowsAnalyzed = window_n,
      NormalizedCovDiff = round(normalized_diff, 4),
      Stationary = is_stationary,
      Interpretation = recommendation,
      stringsAsFactors = FALSE
    )
  )
}

#' Assess Gaussianity of ILR Values
#'
#' Test whether the ILR-transformed compositional data follow a multivariate normal
#' distribution. Sequential Gaussian Simulation assumes normality; departures from
#' normality (e.g., bimodal or skewed distributions) may compromise simulation accuracy.
#'
#' @param ilr_values A matrix or data frame with ILR values (each column is one dimension).
#' @param method Character string specifying normality test: `"anderson"` (default)
#'   uses Anderson-Darling test; `"shapiro"` uses Shapiro-Wilk (univariate only);
#'   `"mardia"` uses Mardia's multivariate skewness/kurtosis.
#' @param plot Logical. If `TRUE` (default), display Q-Q plots and histograms.
#' @param alpha Significance level for normality test (default 0.05).
#'
#' @return A list containing:
#'   - `gaussian`: Logical indicating whether multivariate normality is supported
#'   - `method_used`: Character string indicating which test was used
#'   - `p_values`: P-values from normality tests for each dimension (Anderson) or multivariate (Mardia)
#'   - `skewness`: Multivariate skewness (if method = "mardia")
#'   - `kurtosis`: Multivariate excess kurtosis (if method = "mardia")
#'   - `recommendation`: Character string with interpretation and suggested remedies
#'   - `summary`: Data frame with dimension-wise test statistics
#'
#' @details
#' **Anderson-Darling Test** (univariate, recommended for 2-3D):
#' Tests each ILR dimension separately against normality. More sensitive to tail
#' deviations than Shapiro-Wilk. Dimensions failing the test indicate non-Gaussianity
#' in that specific balance (e.g., a bimodal distribution in sand vs silt balance).
#'
#' **Shapiro-Wilk Test** (univariate):
#' Limited to D ≤ 5000 samples. Less sensitive to tail behavior than Anderson-Darling.
#'
#' **Mardia's Test** (multivariate):
#' Tests for joint departures from normality via multivariate skewness and kurtosis.
#' More holistic but less specific about which dimensions violate normality.
#'
#' **Interpretation and Remedies**:
#' - Gaussian: Proceed with Sequential Gaussian Simulation (SGS)
#' - Non-Gaussian: Consider:
#'   1. **Anamorphosis**: Apply normal-score transform before simulation, back-transform after
#'   2. **Domain stratification**: Bimodal distributions suggest distinct geological units
#'   3. **t-location-scale model**: Alternative to Gaussian for heavier tails
#'   4. Accept and document uncertainty: Some non-Gaussianity is common in real soil data
#'
#' @examples
#' \dontrun{
#' # Test ILR values
#' samples <- data.frame(
#'   sand = c(20, 25, 30, 22, 18, 35),
#'   silt = c(60, 55, 50, 58, 62, 45),
#'   clay = c(20, 20, 20, 20, 20, 20)
#' )
#'
#' ilr_vals <- compositions::ilr(compositions::acomp(samples))
#' colnames(ilr_vals) <- c("ilr1", "ilr2")
#'
#' # Assess Gaussianity
#' gaussianity <- gc_assess_gaussianity(ilr_vals, method = "anderson", plot = TRUE)
#' print(gaussianity$recommendation)
#' }
#'
#' @importFrom stats qqnorm qqline pnorm shapiro.test
#' @export
gc_assess_gaussianity <- function(ilr_values,
                                   method = "anderson",
                                   plot = TRUE,
                                   alpha = 0.05) {
  if (is.data.frame(ilr_values)) {
    ilr_matrix <- as.matrix(ilr_values)
  } else if (is.matrix(ilr_values)) {
    ilr_matrix <- ilr_values
  } else {
    stop("ilr_values must be a matrix or data frame")
  }

  method <- match.arg(method, c("anderson", "shapiro", "mardia"))

  if (method == "anderson") {
    result <- .assess_gaussianity_anderson(ilr_matrix, plot, alpha)
  } else if (method == "shapiro") {
    result <- .assess_gaussianity_shapiro(ilr_matrix, plot, alpha)
  } else {
    result <- .assess_gaussianity_mardia(ilr_matrix, plot, alpha)
  }

  return(result)
}

.anderson_darling_test <- function(x, alpha) {
  n <- length(x)
  x_sorted <- sort(x)

  u <- stats::pnorm(x_sorted, mean = mean(x), sd = stats::sd(x))

  u_clipped <- pmax(pmin(u, 1 - 1e-10), 1e-10)

  a_squared <- -n - (1 / n) * sum((2 * seq_len(n) - 1) * (log(u_clipped) + log(1 - rev(u_clipped))))

  a_squared_adj <- a_squared * (1 + 0.75 / n + 2.25 / n^2)

  if (a_squared_adj < 0.576) {
    p_value <- 1 - exp(-13.436 + 101.14 * a_squared_adj - 223.73 * a_squared_adj^2)
  } else if (a_squared_adj < 0.656) {
    p_value <- 1 - exp(-8.318 + 42.796 * a_squared_adj - 59.938 * a_squared_adj^2)
  } else if (a_squared_adj < 2.491) {
    p_value <- exp(-0.9177 + 4.618 * a_squared_adj - 12.03 * a_squared_adj^2)
  } else if (a_squared_adj < 5.054) {
    p_value <- exp(0.02979 - 0.71313 * (a_squared_adj - 3))
  } else {
    p_value <- exp(0.04153 - 1.111 * (a_squared_adj - 5))
  }

  list(statistic = a_squared_adj, p_value = p_value)
}

.assess_gaussianity_anderson <- function(ilr_matrix, plot, alpha) {
  p_values <- numeric(ncol(ilr_matrix))
  colnames(ilr_matrix) <- colnames(ilr_matrix) %||% paste0("ILR", seq_len(ncol(ilr_matrix)))

  for (i in seq_len(ncol(ilr_matrix))) {
    test_result <- .anderson_darling_test(ilr_matrix[, i], alpha)
    p_values[i] <- test_result$p_value
  }

  is_gaussian <- all(p_values > alpha)

  if (plot) {
    old_par <- graphics::par(mfrow = c(ncol(ilr_matrix), 2), mar = c(4, 4, 2, 1))
    on.exit(graphics::par(old_par))

    for (i in seq_len(ncol(ilr_matrix))) {
      stats::qqnorm(ilr_matrix[, i], main = sprintf("%s Q-Q Plot", colnames(ilr_matrix)[i]))
      stats::qqline(ilr_matrix[, i], col = 2)

      graphics::hist(ilr_matrix[, i], main = sprintf("%s Histogram", colnames(ilr_matrix)[i]),
                     xlab = colnames(ilr_matrix)[i], probability = TRUE)
      x_seq <- seq(min(ilr_matrix[, i]), max(ilr_matrix[, i]), length.out = 100)
      graphics::lines(x_seq, stats::dnorm(x_seq, mean = mean(ilr_matrix[, i]), sd = stats::sd(ilr_matrix[, i])),
                      col = 2, lwd = 2)
    }
  }

  recommendation <- if (is_gaussian) {
    "Multivariate normality supported: Proceed with Sequential Gaussian Simulation (SGS)."
  } else {
    failed_dims <- which(p_values <= alpha)
    sprintf(
      "Non-Gaussianity detected in dimension(s) %s. Recommend anamorphosis (normal-score transform) or domain stratification.",
      paste(failed_dims, collapse = ", ")
    )
  }

  summary_df <- data.frame(
    ILR_Dimension = colnames(ilr_matrix),
    PValue = round(p_values, 4),
    Gaussian = p_values > alpha,
    stringsAsFactors = FALSE
  )

  list(
    gaussian = is_gaussian,
    method_used = "anderson",
    p_values = p_values,
    recommendation = recommendation,
    summary = summary_df
  )
}

.assess_gaussianity_shapiro <- function(ilr_matrix, plot, alpha) {
  p_values <- numeric(ncol(ilr_matrix))
  colnames(ilr_matrix) <- colnames(ilr_matrix) %||% paste0("ILR", seq_len(ncol(ilr_matrix)))

  for (i in seq_len(ncol(ilr_matrix))) {
    test_result <- stats::shapiro.test(ilr_matrix[, i])
    p_values[i] <- test_result$p.value
  }

  is_gaussian <- all(p_values > alpha)

  if (plot) {
    old_par <- graphics::par(mfrow = c(ncol(ilr_matrix), 2), mar = c(4, 4, 2, 1))
    on.exit(graphics::par(old_par))

    for (i in seq_len(ncol(ilr_matrix))) {
      stats::qqnorm(ilr_matrix[, i], main = sprintf("%s Q-Q Plot", colnames(ilr_matrix)[i]))
      stats::qqline(ilr_matrix[, i], col = 2)

      graphics::hist(ilr_matrix[, i], main = sprintf("%s Histogram", colnames(ilr_matrix)[i]),
                     xlab = colnames(ilr_matrix)[i], probability = TRUE)
      x_seq <- seq(min(ilr_matrix[, i]), max(ilr_matrix[, i]), length.out = 100)
      graphics::lines(x_seq, stats::dnorm(x_seq, mean = mean(ilr_matrix[, i]), sd = stats::sd(ilr_matrix[, i])),
                      col = 2, lwd = 2)
    }
  }

  recommendation <- if (is_gaussian) {
    "Multivariate normality supported: Proceed with Sequential Gaussian Simulation (SGS)."
  } else {
    failed_dims <- which(p_values <= alpha)
    sprintf(
      "Non-Gaussianity detected in dimension(s) %s. Recommend anamorphosis (normal-score transform) or domain stratification.",
      paste(failed_dims, collapse = ", ")
    )
  }

  summary_df <- data.frame(
    ILR_Dimension = colnames(ilr_matrix),
    PValue = round(p_values, 4),
    Gaussian = p_values > alpha,
    stringsAsFactors = FALSE
  )

  list(
    gaussian = is_gaussian,
    method_used = "shapiro",
    p_values = p_values,
    recommendation = recommendation,
    summary = summary_df
  )
}

.assess_gaussianity_mardia <- function(ilr_matrix, plot, alpha) {
  n <- nrow(ilr_matrix)
  p <- ncol(ilr_matrix)

  mean_vec <- colMeans(ilr_matrix)
  cov_mat <- stats::cov(ilr_matrix)

  centered <- sweep(ilr_matrix, 2, mean_vec)

  cov_inv <- tryCatch(
    solve(cov_mat),
    error = function(e) {
      warning("Covariance matrix is singular. Using pseudo-inverse.")
      MASS::ginv(cov_mat)
    }
  )

  mahal_dist_sq <- rowSums((centered %*% cov_inv) * centered)

  g1p_sum <- 0
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j) {
        term <- sum((centered[i, ] %*% cov_inv) * centered[j, ])^3
        g1p_sum <- g1p_sum + term
      }
    }
  }
  b1p <- g1p_sum / (n^2)

  b2p <- mean((mahal_dist_sq - p) * (mahal_dist_sq - p - 2))

  skewness <- (n * b1p) / 6
  kurtosis <- (b2p - p * (p + 2)) / sqrt(8 * p * (p + 2))

  skewness_pval <- 2 * (1 - stats::pnorm(abs(skewness)))
  kurtosis_pval <- 2 * (1 - stats::pnorm(abs(kurtosis)))

  is_gaussian <- (skewness_pval > alpha) && (kurtosis_pval > alpha)

  if (plot) {
    graphics::hist(mahal_dist_sq, breaks = 20, probability = TRUE,
                   main = "Mahalanobis Distance Distribution",
                   xlab = "D²", ylab = "Density")
    x_chisq <- seq(0, max(mahal_dist_sq), length.out = 100)
    graphics::lines(x_chisq, stats::dchisq(x_chisq, df = p), col = 2, lwd = 2,
                    label = "Chi-square reference")
  }

  recommendation <- if (is_gaussian) {
    "Multivariate normality supported via Mardia's test: Proceed with Sequential Gaussian Simulation (SGS)."
  } else {
    issues <- c()
    if (skewness_pval <= alpha) issues <- c(issues, "excess skewness")
    if (kurtosis_pval <= alpha) issues <- c(issues, "excess kurtosis")
    sprintf(
      "Non-Gaussianity detected (%s). Recommend anamorphosis (normal-score transform) or domain stratification.",
      paste(issues, collapse = " and ")
    )
  }

  summary_df <- data.frame(
    Test = c("Skewness", "Kurtosis"),
    Statistic = c(round(skewness, 4), round(kurtosis, 4)),
    PValue = c(round(skewness_pval, 4), round(kurtosis_pval, 4)),
    Gaussian = c(skewness_pval > alpha, kurtosis_pval > alpha),
    stringsAsFactors = FALSE
  )

  list(
    gaussian = is_gaussian,
    method_used = "mardia",
    skewness = skewness,
    kurtosis = kurtosis,
    p_values = c(skewness_pval, kurtosis_pval),
    recommendation = recommendation,
    summary = summary_df
  )
}

#' Apply Normal-Score (Anamorphosis) Transform
#'
#' Transform simulated values back from normal space to data space when the original
#' ILR values deviate from Gaussianity. This preserves the statistical properties of
#' the original data distribution while maintaining spatial correlation structure.
#'
#' @param sim_ilr_values Matrix of simulated ILR values in normal space (from SGS).
#' @param ref_ilr_values Matrix of reference (original) ILR values used to estimate
#'   the anamorphosis function.
#' @param despike Logical. If `TRUE` (default), apply empirical despiking to remove
#'   extreme values that may arise from the tail of the transform.
#' @param despike_threshold Numeric. Threshold for despiking (default 1.5 times the
#'   range of reference data).
#'
#' @return A matrix of back-transformed values in the same space as the reference data.
#'
#' @details
#' **Methodology**:
#' The anamorphosis transform is implemented via empirical quantile matching:
#' 1. Compute quantiles of reference data (original)
#' 2. Compute quantiles of simulated data (normal space)
#' 3. Match simulated quantiles to reference quantiles (inverse transform)
#' 4. Interpolate for values between quantiles
#' 5. Optionally despike extreme values
#'
#' **Despiking**:
#' When simulations produce extreme values in the tail of the normal distribution,
#' the back-transform may produce unrealistic values outside the range of the
#' reference data. Despiking caps these values at the reference range bounds.
#'
#' @examples
#' \dontrun{
#' # Simulate ILR values (assume from SGS)
#' sim_values <- matrix(rnorm(100, 0, 1), ncol = 2)
#'
#' # Reference data (original ILR values)
#' ref_values <- matrix(rnorm(50, 0.1, 0.9), ncol = 2)
#'
#' # Apply anamorphosis back-transform
#' back_transformed <- gc_apply_anamorphosis(sim_values, ref_values)
#' print(range(back_transformed))
#' }
#'
#' @importFrom stats quantile approx
#' @export
gc_apply_anamorphosis <- function(sim_ilr_values,
                                   ref_ilr_values,
                                   despike = TRUE,
                                   despike_threshold = NULL) {
  if (ncol(sim_ilr_values) != ncol(ref_ilr_values)) {
    stop("sim_ilr_values and ref_ilr_values must have the same number of dimensions")
  }

  n_sim <- nrow(sim_ilr_values)
  n_dims <- ncol(sim_ilr_values)

  if (is.null(despike_threshold)) {
    despike_threshold <- 1.5 * mean(apply(ref_ilr_values, 2, function(x) diff(range(x))))
  }

  back_transformed <- matrix(0, nrow = n_sim, ncol = n_dims)

  for (dim in seq_len(n_dims)) {
    ref_vals <- ref_ilr_values[, dim]
    sim_vals <- sim_ilr_values[, dim]

    ref_quantiles <- stats::quantile(ref_vals, probs = seq(0, 1, by = 0.1), names = FALSE)
    sim_quantiles <- stats::quantile(sim_vals, probs = seq(0, 1, by = 0.1), names = FALSE)

    back_transform_func <- function(x) {
      stats::approx(sim_quantiles, ref_quantiles, xout = x, method = "linear", rule = 2)$y
    }

    back_transformed[, dim] <- sapply(sim_vals, back_transform_func)

    if (despike) {
      ref_range <- range(ref_vals)
      # Ensure all values are within reference range, replacing NAs with min/max
      values <- back_transformed[, dim]
      values[is.na(values)] <- ref_range[1]
      back_transformed[, dim] <- pmax(pmin(values, ref_range[2]), ref_range[1])
    }
  }

  colnames(back_transformed) <- colnames(sim_ilr_values)

  return(back_transformed)
}
