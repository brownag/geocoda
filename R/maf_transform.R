#' Compute Minimum/Maximum Autocorrelation Factors (MAF) Rotation
#'
#' Compute the MAF rotation matrix from ILR covariance and spatial structure.
#' MAF decorrelates multivariate data at a specified spatial lag, ordering factors
#' by decreasing spatial autocorrelation (high autocorrelation → low autocorrelation).
#'
#' @param ilr_params List from [gc_ilr_params()] containing mean, cov, names
#' @param variogram_models List of fitted gstat vgm objects (one per ILR dimension)
#' @param lag_h Distance for MAF computation. If NULL, defaults to mean_range/3
#' @param method Character: "ilr" (apply to ILR coordinates, default) or reserved for future use
#'
#' @return List with components:
#'   - `rotation_matrix`: D×D rotation matrix for transforming to MAF space
#'   - `eigenvalues`: Vector of eigenvalues ordered by decreasing autocorrelation
#'   - `lag_h`: Distance used for computation
#'   - `method`: Method used ("ilr")
#'   - `pca_matrix`: Lag-0 PCA transformation matrix (internal use)
#'
#' @details
#' **Two-stage decorrelation process:**
#'
#' 1. **Stage 1 - PCA at lag 0**: Decorrelate collocated values
#'    - Compute eigendecomposition of ILR covariance Σ(0)
#'    - Standardize: Y = Σ(0)^(-1/2) (X - μ)
#'
#' 2. **Stage 2 - MAF at lag h**: Diagonalize spatial structure
#'    - Build covariance matrix at lag h: C(h)
#'    - Compute semivariogram matrix: Γ(h) = Σ(0) - C(h)
#'    - Diagonalize in standardized space: Γ_Y(h) = Σ(0)^(-1/2) Γ(h) Σ(0)^(-1/2)
#'    - Extract eigenvectors as MAF rotation
#'
#' **Result**: MAF factors Z = R * Y are ordered by eigenvalues
#' (high eigenvalue → high spatial variance → signal),
#' (low eigenvalue → low spatial variance → noise).
#'
#' @examples
#' \dontrun{
#' # Synthetic compositional data
#' set.seed(42)
#' sand <- rnorm(50, 60, 15)
#' silt <- rnorm(50, 25, 10)
#' clay <- 100 - sand - silt
#' data <- data.frame(x = rnorm(50), y = rnorm(50),
#'                    sand = pmax(0, sand), silt = pmax(0, silt),
#'                    clay = pmax(0, clay))
#'
#' # ILR parameters and variogram
#' ilr_params <- gc_ilr_params(data[, c("sand", "silt", "clay")])
#' vgm_fit <- gc_fit_vgm(ilr_params, data, aggregate = TRUE)
#'
#' # Compute MAF
#' maf_obj <- gc_compute_maf(ilr_params, list(vgm_fit, vgm_fit), lag_h = 10)
#'
#' # MAF rotation matrix
#' print(maf_obj$rotation_matrix)
#' }
#'
#' @export
gc_compute_maf <- function(ilr_params, variogram_models, lag_h = NULL, method = "ilr") {
  stopifnot(is.list(ilr_params), !is.null(ilr_params$cov))
  stopifnot(is.list(variogram_models), length(variogram_models) > 0)
  stopifnot(method %in% c("ilr"))

  n_ilr <- ncol(ilr_params$cov)

  # Default lag_h: mean range / 3
  if (is.null(lag_h)) {
    ranges <- sapply(variogram_models, function(v) {
      if (nrow(v) > 0 && "range" %in% names(v)) {
        v$range[nrow(v)]  # Last row typically has the range
      } else {
        NA_real_
      }
    })
    ranges <- ranges[!is.na(ranges)]
    if (length(ranges) > 0) {
      lag_h <- mean(ranges) / 3
    } else {
      # Fallback if ranges not available
      lag_h <- 1.0
      warning("Could not extract ranges from variograms; using lag_h = 1.0")
    }
  }

  # Stage 1: PCA at lag 0 (decorrelate collocated)
  eig_0 <- eigen(ilr_params$cov, symmetric = TRUE)
  Lambda_0_inv_sqrt <- diag(1 / sqrt(pmax(eig_0$values, 1e-10)))  # Avoid division by zero
  Sigma_inv_sqrt <- eig_0$vectors %*% Lambda_0_inv_sqrt %*% t(eig_0$vectors)

  # Stage 2: Build covariance at lag h
  C_h <- gc_covariance_matrix(ilr_params, variogram_models, lag_h = lag_h)

  # Stage 3: Semivariogram matrix at lag h
  Gamma_h <- ilr_params$cov - C_h

  # Stage 4: Transform to standardized space and diagonalize
  Gamma_Y <- Sigma_inv_sqrt %*% Gamma_h %*% Sigma_inv_sqrt

  # Eigendecomposition of semivariogram in standardized space
  eig_gamma <- eigen(Gamma_Y, symmetric = TRUE)

  # Stage 5: MAF rotation matrix
  # Columns of V are eigenvectors (orthonormal); transpose for row multiplication
  V <- eig_gamma$vectors
  rotation_matrix <- t(V)  # Orthogonal matrix: R @ R^T = I

  # Order by eigenvalue magnitude (high → low, signal → noise)
  ord <- order(eig_gamma$values, decreasing = TRUE)
  rotation_matrix <- rotation_matrix[ord, ]
  eigenvalues <- eig_gamma$values[ord]

  # Return standardized list
  result <- list(
    rotation_matrix = rotation_matrix,
    eigenvalues = eigenvalues,
    lag_h = lag_h,
    method = method,
    pca_matrix = Sigma_inv_sqrt
  )

  class(result) <- "gc_maf_object"
  return(result)
}

#' Apply MAF Transformation to Data
#'
#' Transform ILR coordinates to MAF space using a computed rotation matrix.
#'
#' @param data Data frame containing ILR columns (ilr1, ilr2, ...)
#' @param maf_object List from [gc_compute_maf()] containing rotation_matrix
#' @param ilr_mean Vector of ILR means (optional; used for centering if provided)
#'
#' @return Data frame with new columns maf1, maf2, ... appended to original data
#'
#' @details
#' Applies the transformation: Z = R * (Y - mean(Y))
#' where Y are ILR coordinates and R is the MAF rotation matrix.
#'
#' @examples
#' \dontrun{
#' # After gc_compute_maf()
#' data_ilr <- data.frame(ilr1 = rnorm(20), ilr2 = rnorm(20))
#' data_ilr$ilr_mean <- c(0.1, -0.2)
#'
#' data_maf <- gc_transform_maf(data_ilr, maf_obj)
#' head(data_maf)  # Now has maf1, maf2 columns
#' }
#'
#' @export
gc_transform_maf <- function(data, maf_object, ilr_mean = NULL) {
  stopifnot(is.data.frame(data))
  stopifnot(class(maf_object) == "gc_maf_object")

  # Extract ILR columns
  ilr_cols <- grep("^ilr\\d+$", names(data), value = TRUE)
  if (length(ilr_cols) == 0) {
    stop("No ILR columns found in data")
  }

  ilr_matrix <- as.matrix(data[, ilr_cols])

  # Center: use provided mean or compute from data
  if (is.null(ilr_mean)) {
    ilr_centered <- scale(ilr_matrix, center = TRUE, scale = FALSE)
  } else {
    stopifnot(length(ilr_mean) == ncol(ilr_matrix))
    ilr_centered <- sweep(ilr_matrix, 2, ilr_mean, "-")
  }

  # Two-stage transformation:
  # Stage 1: Standardize (PCA at lag 0)
  Y_standardized <- ilr_centered %*% t(maf_object$pca_matrix)

  # Stage 2: Apply MAF rotation: Z = Y_standardized * V^T
  maf_matrix <- Y_standardized %*% t(maf_object$rotation_matrix)

  # Create output dataframe
  maf_df <- as.data.frame(maf_matrix)
  names(maf_df) <- paste0("maf", seq_len(ncol(maf_matrix)))

  # Append to original data
  result <- cbind(data, maf_df)
  return(result)
}

#' Inverse MAF Transformation
#'
#' Back-transform from MAF space to ILR coordinates.
#'
#' @param maf_data Data frame containing MAF columns (maf1, maf2, ...)
#' @param maf_object List from [gc_compute_maf()] containing rotation_matrix
#' @param ilr_mean Vector of ILR means (optional; used for centering if provided)
#'
#' @return Data frame with new columns ilr1, ilr2, ... appended to original data
#'
#' @details
#' Applies: Y = R^(-1) * Z + μ
#' where Z are MAF factors and R is the MAF rotation matrix.
#'
#' The inverse exists if the rotation matrix is invertible (which it always is
#' for orthogonal MAF rotations).
#'
#' @examples
#' \dontrun{
#' # Back-transform simulated MAF factors to ILR
#' maf_sims <- data.frame(maf1 = rnorm(20), maf2 = rnorm(20))
#' ilr_sims <- gc_inverse_maf(maf_sims, maf_obj, ilr_mean = c(0.1, -0.2))
#' head(ilr_sims)  # Now has ilr1, ilr2 columns
#' }
#'
#' @export
gc_inverse_maf <- function(maf_data, maf_object, ilr_mean = NULL) {
  stopifnot(is.data.frame(maf_data))
  stopifnot(class(maf_object) == "gc_maf_object")

  # Extract MAF columns
  maf_cols <- grep("^maf\\d+$", names(maf_data), value = TRUE)
  if (length(maf_cols) == 0) {
    stop("No MAF columns found in data")
  }

  maf_matrix <- as.matrix(maf_data[, maf_cols])

  # Inverse two-stage transformation:
  # Stage 1: Inverse MAF rotation (Z @ V^T since V is orthonormal)
  rotation_inv <- solve(maf_object$rotation_matrix)
  Y_standardized <- maf_matrix %*% t(rotation_inv)  # Orthogonal: t(solve(t(V))) = V

  # Stage 2: Inverse standardization
  pca_inv <- solve(maf_object$pca_matrix)
  ilr_centered <- Y_standardized %*% t(pca_inv)

  # Add mean if provided
  if (is.null(ilr_mean)) {
    ilr_matrix <- ilr_centered
  } else {
    stopifnot(length(ilr_mean) == ncol(ilr_centered))
    ilr_matrix <- sweep(ilr_centered, 2, ilr_mean, "+")
  }

  # Create output dataframe
  ilr_df <- as.data.frame(ilr_matrix)
  names(ilr_df) <- paste0("ilr", seq_len(ncol(ilr_matrix)))

  # Append to original data
  result <- cbind(maf_data, ilr_df)
  return(result)
}
