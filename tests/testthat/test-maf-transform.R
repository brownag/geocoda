test_that("gc_vgm_to_cov converts variogram to covariance correctly", {
  # Simple variogram: sill=1, nugget=0.1, range=50
  vgm_fit <- gstat::vgm(psill = 0.9, model = "Exp", range = 50, nugget = 0.1)

  # At distance 0, C(0) should equal total sill
  C_0 <- gc_vgm_to_cov(vgm_fit, distance = 0)
  expect_equal(C_0, 1.0, tolerance = 0.01)

  # At large distance, C(h) should approach 0 (uncorrelated)
  C_large <- gc_vgm_to_cov(vgm_fit, distance = 500)
  expect_true(C_large < 0.01)  # Near zero, well below nugget
})

test_that("gc_covariance_matrix builds symmetric matrix", {
  # Create simple ILR parameters
  set.seed(42)
  ilr_cov <- matrix(c(1.0, 0.3, 0.3, 1.0), nrow = 2)
  ilr_params <- list(
    mean = c(0, 0),
    cov = ilr_cov,
    names = c("sand", "silt", "clay")
  )

  # Two variogram models (one per dimension)
  vgm1 <- gstat::vgm(psill = 1.0, model = "Exp", range = 100, nugget = 0.1)
  vgm2 <- gstat::vgm(psill = 1.0, model = "Exp", range = 100, nugget = 0.1)
  vgm_list <- list(vgm1, vgm2)

  # Covariance matrix at lag 50
  C_h <- gc_covariance_matrix(ilr_params, vgm_list, lag_h = 50)

  # Check symmetry
  expect_equal(C_h, t(C_h), tolerance = 1e-10)

  # Check dimensions
  expect_equal(dim(C_h), c(2, 2))

  # Check diagonal is positive
  expect_true(all(diag(C_h) > 0))
})

test_that("gc_compute_maf produces orthogonal rotation matrix", {
  # Create synthetic ILR parameters
  set.seed(42)
  ilr_cov <- matrix(c(1.0, 0.4, 0.4, 1.0), nrow = 2)
  ilr_params <- list(
    mean = c(0, 0),
    cov = ilr_cov,
    names = c("sand", "silt", "clay")
  )

  vgm1 <- gstat::vgm(psill = 1.0, model = "Exp", range = 100, nugget = 0.1)
  vgm2 <- gstat::vgm(psill = 1.0, model = "Exp", range = 100, nugget = 0.1)
  vgm_list <- list(vgm1, vgm2)

  # Compute MAF
  maf_obj <- gc_compute_maf(ilr_params, vgm_list, lag_h = 50)

  # Check rotation matrix is approximately orthogonal: R * R^T ≈ I
  R <- maf_obj$rotation_matrix
  identity_approx <- R %*% t(R)
  expected_identity <- diag(2)

  expect_equal(identity_approx, expected_identity, tolerance = 1e-10)
})

test_that("gc_compute_maf eigenvalues ordered correctly", {
  set.seed(42)
  ilr_cov <- matrix(c(1.0, 0.2, 0.2, 1.0), nrow = 2)
  ilr_params <- list(
    mean = c(0, 0),
    cov = ilr_cov,
    names = c("sand", "silt", "clay")
  )

  vgm1 <- gstat::vgm(psill = 1.0, model = "Exp", range = 100, nugget = 0.1)
  vgm2 <- gstat::vgm(psill = 1.0, model = "Exp", range = 100, nugget = 0.1)
  vgm_list <- list(vgm1, vgm2)

  maf_obj <- gc_compute_maf(ilr_params, vgm_list, lag_h = 50)

  # Eigenvalues should be in decreasing order
  eigenvals <- maf_obj$eigenvalues
  expect_true(eigenvals[1] >= eigenvals[2])
})

test_that("gc_transform_maf produces MAF columns", {
  # Create data with ILR columns
  set.seed(42)
  data <- data.frame(
    x = rnorm(20),
    y = rnorm(20),
    ilr1 = rnorm(20),
    ilr2 = rnorm(20)
  )

  # Create MAF object
  ilr_cov <- matrix(c(1.0, 0.3, 0.3, 1.0), nrow = 2)
  ilr_params <- list(
    mean = c(0, 0),
    cov = ilr_cov,
    names = c("sand", "silt", "clay")
  )

  vgm1 <- gstat::vgm(psill = 1.0, model = "Exp", range = 100, nugget = 0.1)
  vgm2 <- gstat::vgm(psill = 1.0, model = "Exp", range = 100, nugget = 0.1)
  vgm_list <- list(vgm1, vgm2)

  maf_obj <- gc_compute_maf(ilr_params, vgm_list, lag_h = 50)

  # Transform data
  data_maf <- gc_transform_maf(data, maf_obj)

  # Check MAF columns exist
  expect_true("maf1" %in% names(data_maf))
  expect_true("maf2" %in% names(data_maf))

  # Check original columns still present
  expect_true("x" %in% names(data_maf))
  expect_true("y" %in% names(data_maf))
  expect_true("ilr1" %in% names(data_maf))
})

test_that("gc_inverse_maf recovers original ILR values (round-trip)", {
  # Create data with ILR columns
  set.seed(42)
  ilr_original <- data.frame(
    ilr1 = rnorm(20),
    ilr2 = rnorm(20)
  )

  # Create MAF object
  ilr_cov <- matrix(c(1.0, 0.3, 0.3, 1.0), nrow = 2)
  ilr_params <- list(
    mean = c(0.1, -0.2),
    cov = ilr_cov,
    names = c("sand", "silt", "clay")
  )

  vgm1 <- gstat::vgm(psill = 1.0, model = "Exp", range = 100, nugget = 0.1)
  vgm2 <- gstat::vgm(psill = 1.0, model = "Exp", range = 100, nugget = 0.1)
  vgm_list <- list(vgm1, vgm2)

  maf_obj <- gc_compute_maf(ilr_params, vgm_list, lag_h = 50)

  # Forward and backward transform
  data_with_coords <- cbind(x = rnorm(20), y = rnorm(20), ilr_original)
  data_maf <- gc_transform_maf(data_with_coords, maf_obj, ilr_mean = ilr_params$mean)
  data_recovered <- gc_inverse_maf(data_maf, maf_obj, ilr_mean = ilr_params$mean)

  # Extract recovered ILR values
  ilr_recovered <- data_recovered[, c("ilr1", "ilr2")]

  # Should be close to original (within numerical precision)
  expect_equal(as.matrix(ilr_recovered), as.matrix(ilr_original), tolerance = 1e-10)
})

test_that("MAF factors are uncorrelated after transformation", {
  # Create synthetic correlated ILR data (without mvtnorm dependency)
  set.seed(42)
  n <- 50
  x <- runif(n, 0, 100)
  y <- runif(n, 0, 100)

  # Create correlated ILR data using simple multivariate normal
  mu <- c(0, 0)
  sigma <- matrix(c(1.0, 0.5, 0.5, 1.0), nrow = 2)

  # Cholesky decomposition for correlation
  L <- chol(sigma)
  z <- matrix(rnorm(2 * n), nrow = n)
  ilr_data <- z %*% L

  # Create data frame
  data <- data.frame(
    x = x,
    y = y,
    ilr1 = ilr_data[, 1],
    ilr2 = ilr_data[, 2]
  )

  # Compute ILR parameters from data
  ilr_cov <- cov(data[, c("ilr1", "ilr2")])
  ilr_params <- list(
    mean = colMeans(data[, c("ilr1", "ilr2")]),
    cov = ilr_cov,
    names = c("sand", "silt", "clay")
  )

  # Create variogram models
  vgm1 <- gstat::vgm(psill = 0.9, model = "Exp", range = 50, nugget = 0.1)
  vgm2 <- gstat::vgm(psill = 0.9, model = "Exp", range = 50, nugget = 0.1)
  vgm_list <- list(vgm1, vgm2)

  # Compute MAF
  maf_obj <- gc_compute_maf(ilr_params, vgm_list, lag_h = 25)

  # Transform to MAF space
  data_maf <- gc_transform_maf(data, maf_obj, ilr_mean = ilr_params$mean)

  # Extract MAF columns
  maf_data <- data_maf[, c("maf1", "maf2")]

  # Compute correlation between MAF factors
  maf_cor <- cor(maf_data)

  # Off-diagonal correlations should be small (nearly uncorrelated)
  expect_true(abs(maf_cor[1, 2]) < 0.3)  # Relaxed tolerance for synthetic data
})

test_that("gc_ilr_model with use_maf=TRUE stores MAF object", {
  skip_if_not_installed("gstat")

  # Create simple ILR parameters
  ilr_cov <- matrix(c(1.0, 0.3, 0.3, 1.0), nrow = 2)
  ilr_params <- list(
    mean = c(0, 0),
    cov = ilr_cov,
    names = c("sand", "silt", "clay")
  )

  # Variogram model
  vgm_fit <- gstat::vgm(psill = 1.0, model = "Exp", range = 50, nugget = 0.1)

  # Build model with MAF
  model_maf <- gc_ilr_model(ilr_params, vgm_fit, use_maf = TRUE)

  # Check attributes
  expect_true(!is.null(attr(model_maf, "maf_object")))
  expect_true(!is.null(attr(model_maf, "ilr_params")))
  expect_true(!is.null(attr(model_maf, "ilr_mean_original")))
})

test_that("gc_ilr_model with use_maf=FALSE does not store MAF object", {
  skip_if_not_installed("gstat")

  # Create simple ILR parameters
  ilr_cov <- matrix(c(1.0, 0.3, 0.3, 1.0), nrow = 2)
  ilr_params <- list(
    mean = c(0, 0),
    cov = ilr_cov,
    names = c("sand", "silt", "clay")
  )

  # Variogram model
  vgm_fit <- gstat::vgm(psill = 1.0, model = "Exp", range = 50, nugget = 0.1)

  # Build model without MAF (default)
  model_no_maf <- gc_ilr_model(ilr_params, vgm_fit, use_maf = FALSE)

  # Check attributes
  expect_true(is.null(attr(model_no_maf, "maf_object")))
  expect_true(!is.null(attr(model_no_maf, "ilr_params")))
})
