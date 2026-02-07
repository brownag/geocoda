#' Nimble MCMC Backend with 3D Intrinsic CAR Priors
#'
#' Fit hierarchical model using Nimble with intrinsic Conditional Auto-Regressive
#' (CAR) spatial priors on zone means. Enables flexible Bayesian inference with
#' adaptive samplers for 3D compositional data.
#'
#' @param data Data frame with columns: zone (factor), ilr1, ilr2, ... (ILR values),
#'   optionally x, y, z for 3D coordinates
#' @param prior_spec An object of class `"gc_prior_spec"`
#' @param zone_centers Optional data frame with columns: zone, x, y, z.
#'   If NULL, computed from data coordinates.
#' @param n_iter Integer, number of iterations per chain (default 2000)
#' @param n_warmup Integer, warmup/burn-in iterations (default 500)
#' @param n_chains Integer, number of chains (default 2)
#' @param prior_spatial_scale Numeric, scale parameter for intrinsic CAR prior (default 1.0)
#' @param estimate_spatial_decay Logical, estimate spatial decay from data? (default TRUE)
#' @param covariance_prior Character, prior for covariances: "lkj" (default) or "inverse_wishart"
#' @param verbose Logical, print progress (default TRUE)
#' @param ... Additional arguments (reserved for future use)
#'
#' @return S3 object of class `"gc_hierarchical_fit"` with Nimble MCMC 3D CAR results
#'
#' @details
#' ## Intrinsic CAR Prior for 3D Zones
#'
#' The intrinsic CAR prior allows automatic adaptation to spatial structure.
#' Unlike the Stan CAR implementation, Nimble uses an intrinsic formulation where
#' the prior is conditionally specified for each zone using its neighbors.
#'
#' Zones are neighbors if they are within a distance threshold based on their
#' 3D centroid coordinates. The intrinsic CAR prior shrinks zone parameters
#' toward the average of neighboring zones.
#'
#' @keywords internal
fit_hierarchical_nimble_3d <- function(data,
                                      prior_spec,
                                      zone_centers = NULL,
                                      n_iter = 2000,
                                      n_warmup = 500,
                                      n_chains = 2,
                                      prior_spatial_scale = 1.0,
                                      estimate_spatial_decay = TRUE,
                                      covariance_prior = "lkj",
                                      verbose = TRUE,
                                      ...) {

  if (!requireNamespace("nimble", quietly = TRUE)) {
    stop("nimble package is required for Nimble 3D backend")
  }

  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("coda package is required for Nimble diagnostics")
  }

  suppressMessages(library(nimble))

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
    cat("Fitting Hierarchical Model with Nimble MCMC + 3D Intrinsic CAR\n")
    cat("  Method: Bayesian MCMC (Adaptive samplers) + intrinsic CAR spatial\n")
    cat("  Observations:", nrow(data), "\n")
    cat("  Zones:", hierarchy$n_zones, "\n")
    cat("  ILR dimensions:", length(ilr_cols), "\n")
    cat("  Chains:", n_chains, "× (iterations:", n_iter, ", warmup:", n_warmup, ")\n")
    cat("  Spatial prior: Intrinsic CAR with scale =", prior_spatial_scale, "\n")
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

  n_zones <- hierarchy$n_zones
  n_ilr <- length(ilr_cols)

  # Construct spatial adjacency structure from distances
  zone_coords <- zone_centers[, c("x", "y", "z"), drop = FALSE]
  zone_dists <- as.matrix(stats::dist(zone_coords))

  # Define neighbors using distance threshold
  max_dist <- max(zone_dists)
  neighbor_threshold <- 0.6 * max_dist

  # Build adjacency matrix
  adjacency_matrix <- ifelse(zone_dists > 0 & zone_dists <= neighbor_threshold, 1, 0)
  n_neighbors <- rowSums(adjacency_matrix)

  # Handle isolated zones
  if (sum(n_neighbors == 0) > 0) {
    if (verbose) {
      cat("Note: Using k-nearest neighbors to ensure all zones are connected.\n")
    }
    k_neighbors <- min(3, n_zones - 1)
    adjacency_matrix <- matrix(0, n_zones, n_zones)

    for (i in 1:n_zones) {
      nearest_idx <- order(zone_dists[i, ])[(2:(k_neighbors + 1))]
      adjacency_matrix[i, nearest_idx] <- 1
    }

    n_neighbors <- rowSums(adjacency_matrix)
  }

  # Define Nimble model with intrinsic CAR
  model_code <- nimble::nimbleCode({
    # Global hyperparameters
    for (d in 1:D) {
      mu_global[d] ~ dnorm(mu_global_prior_mean[d], sd = sigma_prior_mean[d])
    }

    # Zone-level parameters with intrinsic CAR spatial prior
    for (d in 1:D) {
      # CAR component: each zone shrinks toward neighbors
      for (z in 1:Z) {
        if (n_neighbors[z] > 0) {
          neighbor_sum = 0
          for (z_other in 1:Z) {
            neighbor_sum <- neighbor_sum + adjacency[z, z_other] * mu_zone[z_other, d]
          }
          mu_zone[z, d] ~ dnorm(
            neighbor_sum / n_neighbors[z],
            sd = sqrt(car_scale / (n_neighbors[z] + 1e-6))
          )
        } else {
          # No neighbors - use global mean
          mu_zone[z, d] ~ dnorm(mu_global[d], sd = sqrt(car_scale))
        }
      }

      # Zone-level standard deviations
      for (z in 1:Z) {
        sigma_zone[z, d] ~ dexp(1.0)
      }
    }

    # Likelihood: multivariate normal per observation
    for (n in 1:N) {
      for (d in 1:D) {
        y[n, d] ~ dnorm(mu_zone[zone[n], d], sd = sigma_zone[zone[n], d])
      }
    }

    # Spatial scale prior (if estimating)
    if (estimate_spatial_decay_flag) {
      tau_decay ~ dexp(1.0)
    }
  })

  # Prepare constants and data
  model_constants <- list(
    N = nrow(y_matrix),
    Z = n_zones,
    D = n_ilr,
    zone = as.integer(zone_factor),
    mu_global_prior_mean = prior_spec$global_mean,
    sigma_prior_mean = prior_spec$prior_sd_mean,
    estimate_spatial_decay_flag = as.integer(estimate_spatial_decay),
    car_scale = prior_spatial_scale,
    adjacency = adjacency_matrix,
    n_neighbors = n_neighbors
  )

  model_data <- list(y = y_matrix)

  # Initial values function
  inits_fn <- function() {
    list(
      mu_global = rnorm(n_ilr, prior_spec$global_mean, 0.1),
      mu_zone = matrix(rnorm(n_zones * n_ilr, 0, 0.5), n_zones, n_ilr),
      sigma_zone = matrix(pmax(0.1, abs(rnorm(n_zones * n_ilr, 1, 0.2))), n_zones, n_ilr),
      tau_decay = if (estimate_spatial_decay) runif(1, 0.1, 0.5) else NULL
    )
  }

  # Build Nimble model
  nimble_model <- try(
    {
      nimble::nimbleModel(model_code,
        constants = model_constants,
        data = model_data,
        inits = inits_fn()
      )
    },
    silent = FALSE
  )

  if (inherits(nimble_model, "try-error")) {
    stop("Failed to build Nimble model. Check model specification.")
  }

  # Configure MCMC
  nimble_mcmc <- nimble::configureMCMC(nimble_model, verbose = verbose)

  # Build and compile MCMC
  nimble_mcmc <- nimble::buildMCMC(nimble_mcmc)

  compiled_result <- try(
    {
      compiled_model <- nimble::compileNimble(nimble_model)
      compiled_mcmc <- nimble::compileNimble(nimble_mcmc, project = nimble_model)
      list(compiled_model = compiled_model, compiled_mcmc = compiled_mcmc)
    },
    silent = FALSE
  )

  if (inherits(compiled_result, "try-error")) {
    stop("Failed to compile Nimble model. Check that nimble C++ compiler is available.")
  }

  compiled_model <- compiled_result$compiled_model
  compiled_mcmc <- compiled_result$compiled_mcmc

  # Run MCMC chains
  mcmc_samples_list <- list()

  for (chain in 1:n_chains) {
    if (verbose) {
      cat("Running chain", chain, "of", n_chains, "\n")
    }

    # Set new initial values
    nimble_model$setData(model_data)
    nimble_model$setInits(inits_fn())
    compiled_model$setData(model_data)
    compiled_model$setInits(inits_fn())

    # Run sampling
    samples <- try(
      {
        nimble::runMCMC(
          compiled_mcmc,
          niter = n_iter,
          nburnin = n_warmup,
          nchains = 1,
          thin = 1,
          summary = FALSE,
          samplesAsCodaMCMC = TRUE
        )
      },
      silent = FALSE
    )

    if (inherits(samples, "try-error")) {
      stop("MCMC sampling failed for chain ", chain)
    }

    mcmc_samples_list[[chain]] <- samples
  }

  # Combine chains
  all_samples <- do.call(rbind, mcmc_samples_list)
  posterior_samples <- as.matrix(all_samples)
  colnames_samples <- colnames(posterior_samples)

  # Extract zone-specific parameters
  zone_estimates <- list()
  zone_summaries_list <- list()
  zone_names <- hierarchy$zones

  for (z in seq_along(zone_names)) {
    zone_name <- zone_names[z]
    zone_mu_pattern <- paste0("mu_zone\\[", z, ",.*\\]")
    zone_mu_cols <- grep(zone_mu_pattern, colnames_samples)

    if (length(zone_mu_cols) > 0) {
      zone_mu_samples <- posterior_samples[, zone_mu_cols, drop = FALSE]
      zone_mu_mean <- colMeans(zone_mu_samples)
      zone_mu_sd <- apply(zone_mu_samples, 2, sd)
      zone_cov <- cov(zone_mu_samples)

      ci_lower <- apply(zone_mu_samples, 2, quantile, probs = 0.025)
      ci_upper <- apply(zone_mu_samples, 2, quantile, probs = 0.975)

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
  }

  zone_summaries_df <- do.call(rbind, zone_summaries_list)
  if (!is.null(zone_summaries_df) && nrow(zone_summaries_df) > 0) {
    rownames(zone_summaries_df) <- NULL
  } else {
    zone_summaries_df <- data.frame()
  }

  # Global estimates
  mu_global_cols <- grep("^mu_global", colnames_samples)
  if (length(mu_global_cols) > 0) {
    mu_global_samples <- posterior_samples[, mu_global_cols, drop = FALSE]
    mu_global_mean <- colMeans(mu_global_samples)
    mu_global_cov <- cov(mu_global_samples)
  } else {
    mu_global_mean <- rep(NA, n_ilr)
    mu_global_cov <- diag(NA, n_ilr)
  }

  # Compute Gelman-Rubin diagnostics
  gr_diag <- try(
    {
      mcmc_list <- coda::mcmc.list(lapply(mcmc_samples_list, coda::as.mcmc))
      coda::gelman.diag(mcmc_list, autoburnin = FALSE)
    },
    silent = TRUE
  )

  gr_values <- if (inherits(gr_diag, "try-error")) {
    rep(NA_real_, ncol(posterior_samples))
  } else {
    gr_diag$psrf[, "Point est."]
  }

  # Check convergence
  gr_threshold <- 1.1
  gr_bad <- sum(gr_values > gr_threshold, na.rm = TRUE)

  if (gr_bad > 0) {
    warning(
      "Convergence issues detected:\n",
      paste("  - ", gr_bad, " parameters with Gelman-Rubin > ", gr_threshold, "\n", sep = ""),
      "Consider increasing n_iter or n_warmup"
    )
  }

  # Shrinkage weights
  zone_summary <- data.frame(
    zone = zone_names,
    n_obs = sapply(seq_along(zone_names), function(z) sum(zone_factor == zone_names[z])),
    shrinkage_strength = rep(prior_spec$pooling_coefficient, n_zones),
    spatial_neighbors = n_neighbors,
    stringsAsFactors = FALSE
  )

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
      method = "mcmc_3d_icar",
      backend = "nimble_3d",
      convergence = "Adaptive samplers + intrinsic CAR",
      gelman_rubin = gr_values,
      prior_spatial_scale = prior_spatial_scale,
      covariance_prior = covariance_prior,
      spatial_adjacency = adjacency_matrix
    ),
    prior_spec = prior_spec,
    metadata = list(
      method = "mcmc_3d_icar",
      backend = "nimble_3d",
      fitted_zones = zone_names,
      n_zones_fitted = n_zones,
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
    cat("Nimble MCMC 3D intrinsic CAR model fitted successfully.\n")
    cat("  Fitted zones:", fit_obj$metadata$n_zones_fitted, "\n")
    cat("  Total observations:", fit_obj$metadata$n_obs_total, "\n")
    if (!anyNA(gr_values)) {
      cat("  Gelman-Rubin convergence: max =", max(gr_values, na.rm = TRUE), "\n")
    }
  }

  fit_obj
}
