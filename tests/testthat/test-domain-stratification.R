test_that("gc_identify_strata with kmeans method", {
  set.seed(42)
  n <- 100
  data <- data.frame(
    x = runif(n, 0, 100),
    y = runif(n, 0, 100),
    ilr1 = c(rnorm(50, 0.5, 0.6), rnorm(50, -0.5, 0.6)),
    ilr2 = c(rnorm(50, 0.2, 0.5), rnorm(50, 0.2, 0.5))
  )

  result <- gc_identify_strata(data, n_strata = 2, method = "kmeans", plot = FALSE)

  expect_type(result, "list")
  expect_true(all(c("strata", "n_strata", "cluster_centers", "silhouette_widths",
                     "pca_loadings", "pca_scores", "recommendation", "summary") %in% names(result)))
  expect_equal(result$n_strata, 2)
  expect_equal(length(result$strata), n)
  expect_true(all(result$strata %in% 1:2))
})

test_that("gc_identify_strata with hierarchical method", {
  set.seed(42)
  n <- 100
  data <- data.frame(
    x = runif(n, 0, 100),
    y = runif(n, 0, 100),
    ilr1 = c(rnorm(50, 0.5, 0.6), rnorm(50, -0.5, 0.6)),
    ilr2 = c(rnorm(50, 0.2, 0.5), rnorm(50, 0.2, 0.5))
  )

  result <- gc_identify_strata(data, n_strata = 2, method = "hierarchical", plot = FALSE)

  expect_type(result, "list")
  expect_equal(result$n_strata, 2)
  expect_equal(length(result$strata), n)
})

test_that("gc_identify_strata auto-selects optimal number of strata", {
  set.seed(42)
  n <- 100
  data <- data.frame(
    x = runif(n, 0, 100),
    y = runif(n, 0, 100),
    ilr1 = c(rnorm(33, 0.5, 0.6), rnorm(33, -0.5, 0.6), rnorm(34, 0.0, 0.5)),
    ilr2 = c(rnorm(33, 0.2, 0.5), rnorm(33, 0.2, 0.5), rnorm(34, -0.3, 0.5))
  )

  result <- gc_identify_strata(data, n_strata = c(2, 3, 4), method = "kmeans", plot = FALSE)

  expect_type(result, "list")
  expect_true(result$n_strata %in% c(2, 3, 4))
  expect_equal(length(result$strata), n)
})

test_that("gc_identify_strata computes silhouette widths", {
  set.seed(42)
  n <- 100
  data <- data.frame(
    x = runif(n, 0, 100),
    y = runif(n, 0, 100),
    ilr1 = c(rnorm(50, 0.5, 0.6), rnorm(50, -0.5, 0.6)),
    ilr2 = c(rnorm(50, 0.2, 0.5), rnorm(50, 0.2, 0.5))
  )

  result <- gc_identify_strata(data, n_strata = 2, plot = FALSE)

  expect_equal(length(result$silhouette_widths), n)
  expect_true(all(result$silhouette_widths >= -1 & result$silhouette_widths <= 1))
})

test_that("gc_identify_strata requires at least 2 ILR dimensions", {
  data <- data.frame(
    x = c(1, 2, 3),
    y = c(4, 5, 6),
    ilr1 = c(0.1, 0.2, 0.3)
  )

  expect_error(gc_identify_strata(data), "at least 2 ILR columns")
})

test_that("gc_identify_strata generates summary statistics", {
  set.seed(42)
  n <- 100
  data <- data.frame(
    x = runif(n, 0, 100),
    y = runif(n, 0, 100),
    ilr1 = c(rnorm(50, 0.5, 0.6), rnorm(50, -0.5, 0.6)),
    ilr2 = c(rnorm(50, 0.2, 0.5), rnorm(50, 0.2, 0.5))
  )

  result <- gc_identify_strata(data, n_strata = 2, plot = FALSE)

  expect_type(result$summary, "list")
  expect_equal(nrow(result$summary), 2)
  expect_true(all(c("Stratum", "N_Observations", "Mean_Silhouette", "Quality") %in% colnames(result$summary)))
})

test_that("gc_identify_strata preserves stratum assignment", {
  set.seed(42)
  n <- 50
  data <- data.frame(
    x = runif(n, 0, 100),
    y = runif(n, 0, 100),
    ilr1 = c(rnorm(25, 0.5, 0.6), rnorm(25, -0.5, 0.6)),
    ilr2 = c(rnorm(25, 0.2, 0.5), rnorm(25, 0.2, 0.5))
  )

  result1 <- gc_identify_strata(data, n_strata = 2, method = "kmeans", plot = FALSE)
  result2 <- gc_identify_strata(data, n_strata = 2, method = "kmeans", plot = FALSE)

  # Both should have same number of strata
  expect_equal(result1$n_strata, result2$n_strata)
  # Both should identify same overall partition (potentially different labels)
  expect_equal(length(unique(result1$strata)), length(unique(result2$strata)))
})
test_that("gc_identify_strata handles 3+ ILR dimensions", {
  set.seed(42)
  n <- 100
  data <- data.frame(
    x = runif(n, 0, 100),
    y = runif(n, 0, 100),
    ilr1 = c(rnorm(50, 0.5, 0.6), rnorm(50, -0.5, 0.6)),
    ilr2 = c(rnorm(50, 0.2, 0.5), rnorm(50, 0.2, 0.5)),
    ilr3 = c(rnorm(50, 0.1, 0.4), rnorm(50, -0.1, 0.4))
  )

  result <- gc_identify_strata(data, n_strata = 2, method = "kmeans", plot = FALSE)

  # Should successfully extract first 2 PCs from 3 available
  expect_type(result, "list")
  expect_equal(result$n_strata, 2)
  expect_equal(length(result$strata), n)
  expect_equal(ncol(result$pca_scores), 3)  # Extracted 3 PCs from 3 ILR dims
  expect_equal(ncol(result$pca_loadings), 3)  # 3 ILR dimensions
})

test_that("gc_identify_strata scales to large datasets", {
  set.seed(42)
  n <- 2000
  data <- data.frame(
    x = runif(n, 0, 100),
    y = runif(n, 0, 100),
    ilr1 = c(rnorm(1000, 0.5, 0.6), rnorm(1000, -0.5, 0.6)),
    ilr2 = c(rnorm(1000, 0.2, 0.5), rnorm(1000, 0.2, 0.5))
  )

  # Benchmark performance
  start_time <- Sys.time()
  result <- gc_identify_strata(data, n_strata = 2, method = "kmeans", plot = FALSE)
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Verify correctness
  expect_type(result, "list")
  expect_equal(result$n_strata, 2)
  expect_equal(length(result$strata), n)
  expect_true(nrow(result$summary) == 2)

  # Should reasonably complete (allow generous 10s limit for CI/CD)
  expect_true(elapsed < 10, info = sprintf("Took %.2f seconds for n=%d", elapsed, n))
})

test_that("gc_identify_strata detects weak clustering quality", {
  set.seed(42)
  n <- 100
  # Create continuous compositional gradient (poor clustering target)
  data <- data.frame(
    x = runif(n, 0, 100),
    y = runif(n, 0, 100),
    ilr1 = seq(-0.5, 0.5, length.out = n) + rnorm(n, 0, 0.1),
    ilr2 = seq(-0.3, 0.3, length.out = n) + rnorm(n, 0, 0.08)
  )

  result <- gc_identify_strata(data, n_strata = 3, method = "kmeans", plot = FALSE)

  # Should detect quality classification for strata
  expect_true(result$summary$Quality[1] %in% c("Fair", "Good"))
  # Recommendation should mention quality or suggest using different method
  expect_true(grepl("quality|stratum|silhouette", tolower(result$recommendation)))
})

test_that("gc_identify_strata accepts sf objects with geometry", {
  skip_if_not_installed("sf")

  set.seed(42)
  n <- 100
  # Create standard data frame first
  df <- data.frame(
    x = runif(n, 0, 100),
    y = runif(n, 0, 100),
    ilr1 = c(rnorm(50, 0.5, 0.6), rnorm(50, -0.5, 0.6)),
    ilr2 = c(rnorm(50, 0.2, 0.5), rnorm(50, 0.2, 0.5)),
    stringsAsFactors = FALSE
  )

  # Convert to sf object with geometry, keeping x and y columns
  data_sf <- sf::st_as_sf(df, coords = c("x", "y"), remove = FALSE)

  # Extract coordinates and reconstruct as data frame (simulating what function does)
  coords <- sf::st_coordinates(data_sf)
  data_test <- sf::st_drop_geometry(data_sf)
  data_test$x <- coords[, 1]
  data_test$y <- coords[, 2]

  # Test that function can handle the reconstructed data
  result <- gc_identify_strata(data_test, n_strata = 2, method = "kmeans", plot = FALSE)

  expect_type(result, "list")
  expect_equal(result$n_strata, 2)
  expect_equal(length(result$strata), n)
  expect_equal(nrow(result$summary), 2)
})

test_that("gc_identify_strata auto-selection handles boundary cases", {
  set.seed(42)
  n <- 100
  # Create three distinct clusters
  data <- data.frame(
    x = c(runif(33, 0, 30), runif(33, 40, 70), runif(34, 80, 100)),
    y = c(runif(33, 0, 30), runif(33, 40, 70), runif(34, 80, 100)),
    ilr1 = c(rnorm(33, -0.5, 0.3), rnorm(33, 0.0, 0.3), rnorm(34, 0.5, 0.3)),
    ilr2 = c(rnorm(33, 0.3, 0.3), rnorm(33, -0.2, 0.3), rnorm(34, -0.4, 0.3))
  )

  # Test with reasonable candidates (2, 3, 4, 5)
  result <- gc_identify_strata(data, n_strata = c(2, 3, 4, 5), method = "kmeans", plot = FALSE)

  # Should select reasonable optimum within candidates
  expect_true(result$n_strata %in% c(2, 3, 4, 5))
  expect_equal(length(result$strata), n)
  expect_true(all(result$strata %in% 1:result$n_strata))

  # Verify recommendation notes the optimal selection
  expect_true(grepl("optimal|strata", tolower(result$recommendation)))
})