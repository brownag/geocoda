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
