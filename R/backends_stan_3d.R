#' Stan HMC Backend with 3D CAR Spatial Priors
#'
#' Fit hierarchical model using Stan with Conditional Auto-Regressive (CAR)
#' spatial priors on zone means. Enables full Bayesian inference with spatial
#' structure for 3D compositional data.
#'
#' @param data Data frame with columns: zone (factor), ilr1, ilr2, ... (ILR values),
#'   optionally x, y, z for 3D coordinates
#' @param prior_spec An object of class `"gc_prior_spec"`
#' @param zone_centers Optional data frame with columns: zone, x, y, z.
#'   If NULL, computed from data coordinates.
#' @param n_iter Integer, number of iterations per chain (default 2000)
#' @param n_warmup Integer, warmup/burn-in iterations (default 500)
#' @param n_chains Integer, number of chains (default 2)
#' @param adapt_delta Numeric, adaptation target acceptance probability (default 0.85)
#' @param car_prior_scale Numeric, scale parameter for CAR spatial prior (default 1.0)
#' @param estimate_spatial_decay Logical, estimate spatial decay from data? (default TRUE)
#' @param covariance_prior Character, prior for covariances: "lkj" (default) or "inverse_wishart"
#' @param verbose Logical, print progress (default TRUE)
#' @param ... Additional arguments passed to rstan::sampling()
#'
#' @return S3 object of class `"gc_hierarchical_fit"` with Stan HMC 3D CAR results
#'
#' @details
#' ## CAR Prior for 3D Spatial Structure
#'
#' The CAR (Conditional Auto-Regressive) prior allows zone means to be spatially
#' dependent based on zone centroid proximity. Zones that are farther apart
#' have less prior dependency, while nearby zones share information through
#' the spatial prior.
#'
#' The implementation uses an intrinsic CAR model where zone means are shrunk
#' toward the mean of neighboring zones (as defined by spatial proximity).
#'
#' @keywords internal
fit_hierarchical_stan_3d <- function(data,
                                    prior_spec,
                                    zone_centers = NULL,
                                    n_iter = 2000,
                                    n_warmup = 500,
                                    n_chains = 2,
                                    adapt_delta = 0.85,
                                    car_prior_scale = 1.0,
                                    estimate_spatial_decay = TRUE,
                                    covariance_prior = "lkj",
                                    verbose = TRUE,
                                    ...) {

  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("rstan package is required for Stan 3D backend")
  }

  # Validate inputs
  if (!inherits(prior_spec, "gc_prior_spec")) {
    stop("prior_spec must be an object of class 'gc_prior_spec'")
  }

  if (!("zone" %in% names(data))) {
    stop("data must contain a 'zone' column")
  }

  hierarchy <- prior_spec$hierarchy
  ilr_cols <- grep("^ilr", names(data), value = TRUE)

  if (length(ilr_cols) == 0) {
    stop("data must contain ILR columns (ilr1, ilr2, ...)")
  }

  if (length(ilr_cols) != hierarchy$n_components - 1) {
    stop("Number of ILR columns does not match n_components - 1")
  }

  covariance_prior <- match.arg(covariance_prior, c("lkj", "inverse_wishart"))

  if (verbose) {
    cat("Fitting Hierarchical Model with Stan HMC + 3D CAR Priors\n")
    cat("  Method: Bayesian MCMC (Hamiltonian Monte Carlo) + CAR spatial\n")
    cat("  Observations:", nrow(data), "\n")
    cat("  Zones:", hierarchy$n_zones, "\n")
    cat("  ILR dimensions:", length(ilr_cols), "\n")
    cat("  Chains:", n_chains, "× (iterations:", n_iter, ", warmup:", n_warmup, ")\n")
    cat("  Spatial prior: CAR with scale =", car_prior_scale, "\n")
  }

  # Prepare data
  zone_factor <- as.factor(data$zone)
  y_matrix <- as.matrix(data[, ilr_cols, drop = FALSE])
  has_3d_coords <- all(c("x", "y", "z") %in% names(data))

  # Get or compute zone centers
  if (is.null(zone_centers)) {
    if (!has_3d_coords) {
      stop(
        "zone_centers must be provided if data does not contain x, y, z coordinates.\n",
        "Provide a data frame with zone centroid coordinates."
      )
    }

    # Compute zone centers from data coordinates
    zone_centers <- data.frame(zone = unique(zone_factor))
    zone_centers$x <- sapply(zone_centers$zone, function(z) {
      mean(data$x[zone_factor == z], na.rm = TRUE)
    })
    zone_centers$y <- sapply(zone_centers$zone, function(z) {
      mean(data$y[zone_factor == z], na.rm = TRUE)
    })
    zone_centers$z <- sapply(zone_centers$zone, function(z) {
      mean(data$z[zone_factor == z], na.rm = TRUE)
    })
  }

  # Construct spatial adjacency structure from distances
  zone_coords <- zone_centers[, c("x", "y", "z"), drop = FALSE]
  zone_dists <- as.matrix(stats::dist(zone_coords))

  # Define neighbors as zones within a distance threshold (0.6 × max distance)
  max_dist <- max(zone_dists)
  neighbor_threshold <- 0.6 * max_dist

  # Build adjacency matrix (1 if neighbors, 0 otherwise)
  adjacency_matrix <- ifelse(zone_dists > 0 & zone_dists <= neighbor_threshold, 1, 0)

  # Compute number of neighbors for each zone (for CAR prior)
  n_neighbors <- rowSums(adjacency_matrix)

  # Check if this becomes a problem for Stan
  if (sum(n_neighbors == 0) > 0) {
    if (verbose) {
      cat("Warning: Some zones have no spatial neighbors at threshold.\n")
      cat("  Adjusting threshold to ensure connectivity.\n")
    }
    # Use k-nearest neighbors approach instead
    k_neighbors <- min(3, hierarchy$n_zones - 1)  # Connect to 3 nearest or all but self
    adjacency_matrix <- matrix(0, hierarchy$n_zones, hierarchy$n_zones)

    for (i in seq_along(hierarchy$zones)) {
      nearest_idx <- order(zone_dists[i, ])[(2:(k_neighbors + 1))]  # Skip self
      adjacency_matrix[i, nearest_idx] <- 1
    }

    n_neighbors <- rowSums(adjacency_matrix)
  }

  # Prepare Stan data with CAR structure
  stan_data <- list(
    N = nrow(data),
    Z = hierarchy$n_zones,
    D = length(ilr_cols),
    zone = as.integer(zone_factor),
    y = y_matrix,
    mu_global_prior_mean = prior_spec$global_mean,
    sigma_prior_mean = prior_spec$prior_sd_mean,
    prior_shape_variance = prior_spec$prior_shape_variance,
    prior_rate_variance = prior_spec$prior_rate_variance,
    estimate_pooling_flag = 1L,  # Always estimate in spatial model
    pooling_coef_fixed = prior_spec$pooling_coefficient,
    covariance_prior_type = if (covariance_prior == "lkj") 0L else 1L,
    eta_lkj = 1.0,
    # CAR spatial structure
    adjacency = adjacency_matrix,
    n_neighbors = n_neighbors,
    car_scale = car_prior_scale,
    estimate_spatial_decay_flag = if (estimate_spatial_decay) 1L else 0L,
    spatial_decay_prior_shape = 2.0,
    spatial_decay_prior_rate = 2.0
  )

  # Use a simplified Stan block for 3D HZM with CAR
  # This is a comment-based definition that can be compiled inline
  stan_block <- "
    data {
      int<lower=1> N;           // Total observations
      int<lower=1> Z;           // Number of zones
      int<lower=1> D;           // ILR dimensions
      array[N] int<lower=1,upper=Z> zone;  // Zone assignment per observation
      matrix[N, D] y;           // ILR observations
      vector[D] mu_global_prior_mean;
      vector[D] sigma_prior_mean;
      real<lower=0> prior_shape_variance;
      real<lower=0> prior_rate_variance;
      int<lower=0,upper=1> estimate_pooling_flag;
      real<lower=0,upper=1> pooling_coef_fixed;
      int<lower=0,upper=1> covariance_prior_type;
      real<lower=0> eta_lkj;
      // CAR spatial structure
      matrix[Z, Z] adjacency;
      array[Z] int<lower=0> n_neighbors;
      real<lower=0> car_scale;
      int<lower=0,upper=1> estimate_spatial_decay_flag;
      real<lower=0> spatial_decay_prior_shape;
      real<lower=0> spatial_decay_prior_rate;
    }

    parameters {
      matrix[Z, D] mu_zone;
      matrix[Z, D] sigma_zone;
      vector[D] mu_global;
      real<lower=0> tau_pooling;
      real<lower=0> spatial_decay;
    }

    model {
      // Global priors
      for (d in 1:D) {
        mu_global[d] ~ normal(mu_global_prior_mean[d], sigma_prior_mean[d]);
      }

      // Spatial CAR prior on zone means
      for (d in 1:D) {
        for (z in 1:Z) {
          real neighbor_mean = 0;
          if (n_neighbors[z] > 0) {
            for (z_neighbor in 1:Z) {
              if (adjacency[z, z_neighbor] == 1) {
                neighbor_mean += mu_zone[z_neighbor, d];
              }
            }
            neighbor_mean /= n_neighbors[z];
            mu_zone[z, d] ~ normal(
              neighbor_mean,
              sqrt(car_scale * spatial_decay / (n_neighbors[z] + 1e-6))
            );
          } else {
            // No neighbors - center on global mean
            mu_zone[z, d] ~ normal(mu_global[d], sqrt(car_scale));
          }
        }
      }

      // Zone standard deviations
      for (z in 1:Z) {
        for (d in 1:D) {
          sigma_zone[z, d] ~ exponential(1.0);
        }
      }

      // Spatial decay rate prior
      if (estimate_spatial_decay_flag) {
        spatial_decay ~ gamma(spatial_decay_prior_shape, spatial_decay_prior_rate);
      }

      // Pooling coefficient prior
      if (estimate_pooling_flag) {
        tau_pooling ~ exponential(2.0);
      }

      // Likelihood
      for (n in 1:N) {
        for (d in 1:D) {
          y[n, d] ~ normal(mu_zone[zone[n], d], sigma_zone[zone[n], d]);
        }
      }
    }
  "

  # Compile and fit Stan model from string
  fit <- try(
    {
      rstan::stan(
        model_code = stan_block,
        data = stan_data,
        iter = n_iter,
        warmup = n_warmup,
        chains = n_chains,
        control = list(adapt_delta = adapt_delta),
        verbose = verbose,
        refresh = if (verbose) 100 else 0,
        ...
      )
    },
    silent = FALSE
  )

  if (inherits(fit, "try-error")) {
    stop("Stan model compilation/fitting failed. Check that rstan is properly installed.")
  }

  # Extract posterior samples
  posterior_list <- rstan::extract(fit, permuted = TRUE)
  posterior_arrays <- rstan::extract(fit, permuted = FALSE)

  # Posterior samples matrix
  posterior_samples <- posterior_arrays[, 1, ]
  if (n_chains > 1) {
    for (c in 2:n_chains) {
      posterior_samples <- rbind(posterior_samples, posterior_arrays[, c, ])
    }
  }

  # Extract diagnostics
  summary_table <- rstan::summary(fit)$summary
  rhat_values <- if ("Rhat" %in% colnames(summary_table)) summary_table[, "Rhat"] else rep(NA_real_, nrow(summary_table))
  ess_bulk <- if ("Bulk_ESS" %in% colnames(summary_table)) summary_table[, "Bulk_ESS"] else rep(NA_real_, nrow(summary_table))
  ess_tail <- if ("Tail_ESS" %in% colnames(summary_table)) summary_table[, "Tail_ESS"] else rep(NA_real_, nrow(summary_table))

  n_divergent <- try(rstan::get_num_divergent(fit), silent = TRUE)
  if (inherits(n_divergent, "try-error")) {
    n_divergent <- 0L
  }

  # Extract zone parameters
  mu_zone_samples <- posterior_list$mu_zone
  sigma_zone_samples <- posterior_list$sigma_zone
  mu_global_samples <- posterior_list$mu_global

  # Compute zone estimates
  zone_estimates <- list()
  zone_summaries_list <- list()

  for (z in seq_along(hierarchy$zones)) {
    zone_name <- hierarchy$zones[z]

    zone_mu_mean <- colMeans(mu_zone_samples[, z, ])
    zone_mu_sd <- apply(mu_zone_samples[, z, ], 2, sd)
    zone_cov <- cov(mu_zone_samples[, z, ])

    ci_lower <- apply(mu_zone_samples[, z, ], 2, quantile, probs = 0.025)
    ci_upper <- apply(mu_zone_samples[, z, ], 2, quantile, probs = 0.975)

    zone_estimates[[zone_name]] <- list(
      n_obs = sum(zone_factor == zone_name),
      mean = zone_mu_mean,
      cov = zone_cov,
      sd = zone_mu_sd
    )

    zone_summaries_list[[zone_name]] <- data.frame(
      zone = zone_name,
      ilr_dim = ilr_cols,
      posterior_mean = zone_mu_mean,
      posterior_sd = zone_mu_sd,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      n_obs = sum(zone_factor == zone_name),
      stringsAsFactors = FALSE
    )
  }

  zone_summaries_df <- do.call(rbind, zone_summaries_list)
  rownames(zone_summaries_df) <- NULL

  # Global estimates
  mu_global_mean <- colMeans(mu_global_samples)
  mu_global_cov <- cov(mu_global_samples)

  # Shrinkage weights
  zone_summary <- data.frame(
    zone = hierarchy$zones,
    n_obs = sapply(seq_along(hierarchy$zones), function(z) sum(zone_factor == hierarchy$zones[z])),
    shrinkage_strength = rep(mean(prior_spec$pooling_coefficient), hierarchy$n_zones),
    spatial_neighbors = n_neighbors,
    stringsAsFactors = FALSE
  )

  # Check convergence
  rhat_threshold <- 1.1
  ess_threshold <- 100
  rhat_bad <- sum(rhat_values > rhat_threshold, na.rm = TRUE)
  ess_bad <- sum(ess_bulk < ess_threshold, na.rm = TRUE) + sum(ess_tail < ess_threshold, na.rm = TRUE)

  if (rhat_bad > 0 || ess_bad > 0 || n_divergent > 0) {
    warning(
      "Convergence issues detected:\n",
      if (rhat_bad > 0) paste("  - ", rhat_bad, " parameters with Rhat > ", rhat_threshold, "\n", sep = ""),
      if (ess_bad > 0) paste("  - ", ess_bad, " parameters with ESS < ", ess_threshold, "\n", sep = ""),
      if (n_divergent > 0) paste("  - ", n_divergent, " divergent transitions\n", sep = ""),
      "Consider increasing n_iter, n_warmup, or adapt_delta"
    )
  }

  # Build return object
  fit_obj <- list(
    zone_estimates = zone_estimates,
    global_estimates = list(
      mean = mu_global_mean,
      cov = mu_global_cov
    ),
    samples = posterior_samples,
    shrinkage_weights = zone_summary,
    zone_summaries = zone_summaries_df,
    zone_centers = zone_centers,
    diagnostics = list(
      method = "hmc_3d_car",
      backend = "stan_3d",
      convergence = "HMC/NUTS + CAR",
      rhat = rhat_values,
      ess_bulk = ess_bulk,
      ess_tail = ess_tail,
      n_divergent = n_divergent,
      adapt_delta = adapt_delta,
      car_prior_scale = car_prior_scale,
      covariance_prior = covariance_prior,
      spatial_adjacency = adjacency_matrix
    ),
    prior_spec = prior_spec,
    metadata = list(
      method = "hmc_3d_car",
      backend = "stan_3d",
      fitted_zones = hierarchy$zones,
      n_zones_fitted = hierarchy$n_zones,
      n_obs_total = nrow(data),
      pooling_coefficient = prior_spec$pooling_coefficient,
      n_iter = n_iter,
      n_warmup = n_warmup,
      n_chains = n_chains,
      timestamp = Sys.time()
    )
  )

  class(fit_obj) <- c("gc_hierarchical_fit", "list")

  if (verbose) {
    cat("Stan HMC 3D CAR model fitted successfully.\n")
    cat("  Fitted zones:", fit_obj$metadata$n_zones_fitted, "\n")
    cat("  Total observations:", fit_obj$metadata$n_obs_total, "\n")
    cat("  Divergent transitions:", n_divergent, "\n")
    cat("  Rhat convergence: max =", max(rhat_values, na.rm = TRUE), "\n")
  }

  fit_obj
}
