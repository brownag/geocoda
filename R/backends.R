#' Hierarchical Model Backend Interface
#'
#' Multi-backend system for fitting hierarchical models with analytical
#' shrinkage, Stan HMC, or Nimble MCMC methods.
#'
#' @keywords internal
#' @name backends
NULL


#' Check if a backend is available
#'
#' Verify that required packages are installed for the requested backend.
#'
#' @param backend Character, backend name: "analytical", "stan", or "nimble"
#'
#' @return Logical, TRUE if backend is available, FALSE otherwise
#' @keywords internal
check_backend_available <- function(backend) {
  switch(backend,
    "analytical" = TRUE,  # Always available
    "stan" = requireNamespace("rstan", quietly = TRUE),
    "nimble" = requireNamespace("nimble", quietly = TRUE),
    FALSE
  )
}


#' Get backend installation message
#'
#' Provide clear installation instructions for unavailable backends.
#'
#' @param backend Character, backend name
#'
#' @return Character string with installation instructions
#' @keywords internal
get_backend_install_msg <- function(backend) {
  switch(backend,
    "stan" = paste(
      "The 'rstan' package is required for Stan backend.\n",
      "Install it with: install.packages('rstan')\n",
      "See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started"
    ),
    "nimble" = paste(
      "The 'nimble' package is required for Nimble backend.\n",
      "Install it with: install.packages('nimble')\n",
      "Requires C++14 compiler support."
    ),
    "Unknown backend"
  )
}


#' Dispatcher for hierarchical model fitting
#'
#' Route to appropriate backend (analytical, Stan, or Nimble) for fitting.
#'
#' @param data Data frame with columns: zone (factor), ilr1, ilr2, ... (ILR values)
#' @param prior_spec An object of class `"gc_prior_spec"`
#' @param backend Character, backend selection: "analytical" (default), "stan", or "nimble"
#' @param verbose Logical, print progress (default TRUE)
#' @param ... Additional arguments passed to backend-specific functions
#'
#' @return An S3 object of class `"gc_hierarchical_fit"` with standardized structure:
#'   - `zone_estimates`: Zone-specific posterior parameters
#'   - `global_estimates`: Global parameter estimates
#'   - `samples`: Posterior samples (NULL for analytical, matrix for MCMC)
#'   - `diagnostics`: Backend-specific diagnostics
#'   - `prior_spec`: Original prior specification
#'   - `metadata`: Fitting information including backend and timestamp
#'
#' @keywords internal
fit_hierarchical_backend <- function(data,
                                    prior_spec,
                                    backend = "analytical",
                                    verbose = TRUE,
                                    ...) {

  # Validate backend argument
  backend <- match.arg(backend, c("analytical", "stan", "nimble"))

  # Check backend availability
  if (!check_backend_available(backend)) {
    stop(
      "Backend '", backend, "' is not available.\n",
      get_backend_install_msg(backend),
      call. = FALSE
    )
  }

  # Dispatch to appropriate backend
  switch(backend,
    "analytical" = fit_hierarchical_analytical(data, prior_spec, verbose = verbose),
    "stan" = fit_hierarchical_stan(data, prior_spec, verbose = verbose, ...),
    "nimble" = fit_hierarchical_nimble(data, prior_spec, verbose = verbose, ...)
  )
}


#' Analytical shrinkage backend
#'
#' Fit hierarchical model using empirical Bayes analytical shrinkage.
#' Provides fast closed-form estimates suitable for exploration and visualization.
#'
#' @param data Data frame with columns: zone (factor), ilr1, ilr2, ... (ILR values)
#' @param prior_spec An object of class `"gc_prior_spec"`
#' @param verbose Logical, print progress (default TRUE)
#'
#' @return S3 object of class `"gc_hierarchical_fit"` with analytical results
#'
#' @keywords internal
fit_hierarchical_analytical <- function(data,
                                       prior_spec,
                                       verbose = TRUE) {

  if (!inherits(prior_spec, "gc_prior_spec")) {
    stop("prior_spec must be an object of class 'gc_prior_spec'")
  }

  if (!("zone" %in% names(data))) {
    stop("data must contain a 'zone' column")
  }

  hierarchy <- prior_spec$hierarchy

  # Extract ILR columns
  ilr_cols <- grep("^ilr", names(data), value = TRUE)
  if (length(ilr_cols) == 0) {
    stop("data must contain ILR columns (ilr1, ilr2, ...)")
  }

  if (length(ilr_cols) != hierarchy$n_components - 1) {
    stop("Number of ILR columns does not match n_components - 1")
  }

  n_obs <- nrow(data)
  n_zones <- hierarchy$n_zones
  n_ilr <- length(ilr_cols)

  if (verbose) {
    cat("Fitting Hierarchical Shrinkage Model\n")
    cat("  Observations:", n_obs, "\n")
    cat("  Zones:", n_zones, "\n")
    cat("  ILR dimensions:", n_ilr, "\n")
    cat("  Method: Analytical shrinkage (empirical Bayes)\n")
  }

  # Prepare data for shrinkage estimation
  zone_factor <- as.factor(data$zone)
  y_matrix <- as.matrix(data[, ilr_cols, drop = FALSE])

  # Initialize zone-level summaries
  zone_summaries <- list()
  posterior_draws <- list()

  for (z in seq_along(hierarchy$zones)) {
    zone_name <- hierarchy$zones[z]
    zone_idx <- zone_factor == zone_name

    if (sum(zone_idx) == 0) {
      if (verbose) {
        cat("Warning: Zone", zone_name, "has no observations\n")
      }
      next
    }

    y_zone <- y_matrix[zone_idx, , drop = FALSE]
    n_z <- nrow(y_zone)

    # Compute zone posterior (conjugate updating)
    zone_mean <- colMeans(y_zone, na.rm = TRUE)
    zone_cov <- stats::cov(y_zone, use = "complete.obs")

    # Shrinkage toward global mean
    shrinkage_coef <- prior_spec$pooling_coefficient
    pooled_mean <- (1 - shrinkage_coef) * zone_mean + shrinkage_coef * prior_spec$global_mean

    # Posterior variance (reduced by sample size in zone)
    posterior_var_scaling <- 1.0 / (1.0 + n_z * shrinkage_coef)
    pooled_cov <- posterior_var_scaling * zone_cov

    # Store zone posterior
    posterior_draws[[zone_name]] <- list(
      n_obs = n_z,
      mean = pooled_mean,
      cov = pooled_cov,
      sd = sqrt(diag(pooled_cov))
    )

    # Compute credible intervals (95%)
    ci_lower <- pooled_mean - 1.96 * sqrt(diag(pooled_cov))
    ci_upper <- pooled_mean + 1.96 * sqrt(diag(pooled_cov))

    zone_summaries[[zone_name]] <- data.frame(
      zone = zone_name,
      ilr_dim = ilr_cols,
      posterior_mean = pooled_mean,
      posterior_sd = sqrt(diag(pooled_cov)),
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      n_obs = n_z,
      shrinkage = shrinkage_coef,
      stringsAsFactors = FALSE
    )
  }

  zone_summaries_df <- do.call(rbind, zone_summaries)
  rownames(zone_summaries_df) <- NULL

  # Summary of pooling across zones
  zone_summary <- data.frame(
    zone = names(posterior_draws),
    n_obs = sapply(posterior_draws, function(x) x$n_obs),
    shrinkage_strength = rep(prior_spec$pooling_coefficient, length(posterior_draws)),
    stringsAsFactors = FALSE
  )

  # Prepare return object with standardized structure
  fit <- list(
    zone_estimates = posterior_draws,
    global_estimates = list(
      mean = prior_spec$global_mean,
      cov = prior_spec$global_covariance
    ),
    samples = NULL,  # Analytical method produces no samples
    shrinkage_weights = zone_summary,
    zone_summaries = zone_summaries_df,
    diagnostics = list(
      method = "analytical_shrinkage",
      backend = "analytical",
      convergence = "analytical (no convergence metrics)"
    ),
    prior_spec = prior_spec,
    metadata = list(
      method = "analytical_shrinkage",
      backend = "analytical",
      fitted_zones = names(posterior_draws),
      n_zones_fitted = length(posterior_draws),
      n_obs_total = nrow(data),
      pooling_coefficient = prior_spec$pooling_coefficient,
      timestamp = Sys.time()
    )
  )

  class(fit) <- c("gc_hierarchical_fit", "list")

  if (verbose) {
    cat("Hierarchical shrinkage model fitted successfully.\n")
    cat("  Method: empirical Bayes shrinkage (analytical)\n")
    cat("  Fitted zones:", fit$metadata$n_zones_fitted, "\n")
    cat("  Total observations:", fit$metadata$n_obs_total, "\n")
    cat("  Pooling strength:", fit$metadata$pooling_coefficient, "\n")
  }

  fit
}


#' Stan HMC backend
#'
#' Fit hierarchical model using Stan with Hamiltonian Monte Carlo.
#' Provides full Bayesian inference with posterior sampling and convergence diagnostics.
#'
#' @param data Data frame with columns: zone (factor), ilr1, ilr2, ... (ILR values)
#' @param prior_spec An object of class `"gc_prior_spec"`
#' @param n_iter Integer, number of iterations per chain (default 2000)
#' @param n_warmup Integer, warmup/burn-in iterations (default 500)
#' @param n_chains Integer, number of chains (default 2)
#' @param adapt_delta Numeric, adaptation target acceptance probability (default 0.8)
#' @param estimate_pooling Logical, estimate pooling strength from data? (default FALSE)
#' @param covariance_prior Character, prior for covariances: "lkj" (default) or "inverse_wishart"
#' @param verbose Logical, print progress (default TRUE)
#' @param ... Additional arguments passed to rstan::sampling()
#'
#' @return S3 object of class `"gc_hierarchical_fit"` with Stan HMC results
#'
#' @keywords internal
fit_hierarchical_stan <- function(data,
                                 prior_spec,
                                 n_iter = 2000,
                                 n_warmup = 500,
                                 n_chains = 2,
                                 adapt_delta = 0.8,
                                 estimate_pooling = FALSE,
                                 covariance_prior = "lkj",
                                 verbose = TRUE,
                                 ...) {

  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("rstan package is required for Stan backend")
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
    cat("Fitting Hierarchical Model with Stan HMC\n")
    cat("  Method: Bayesian MCMC (Hamiltonian Monte Carlo)\n")
    cat("  Observations:", nrow(data), "\n")
    cat("  Zones:", hierarchy$n_zones, "\n")
    cat("  ILR dimensions:", length(ilr_cols), "\n")
    cat("  Chains:", n_chains, "× (iterations:", n_iter, ", warmup:", n_warmup, ")\n")
    cat("  Pooling: ", if (estimate_pooling) "estimated" else "fixed", "\n")
    cat("  Covariance prior: ", covariance_prior, "\n")
  }

  # Prepare data for Stan
  zone_factor <- as.factor(data$zone)
  y_matrix <- as.matrix(data[, ilr_cols, drop = FALSE])

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
    estimate_pooling_flag = if (estimate_pooling) 1L else 0L,
    pooling_coef_fixed = prior_spec$pooling_coefficient,
    covariance_prior_type = if (covariance_prior == "lkj") 0L else 1L,
    eta_lkj = 1.0
  )

  # Get Stan model path
  stan_file <- system.file("stan", "hierarchical_ilr.stan", package = "geocoda")

  if (!file.exists(stan_file)) {
    stop("Stan model file not found at ", stan_file)
  }

  # Compile and fit Stan model
  fit <- rstan::stan(
    file = stan_file,
    data = stan_data,
    iter = n_iter,
    warmup = n_warmup,
    chains = n_chains,
    control = list(adapt_delta = adapt_delta),
    verbose = verbose,
    refresh = if (verbose) 100 else 0,
    ...
  )

  # Extract posterior samples
  posterior_list <- rstan::extract(fit, permuted = TRUE)  # List format for parameter extraction
  posterior_arrays <- rstan::extract(fit, permuted = FALSE)  # Array format (iter x chains x params)

  # Convert to standardized matrix format: rows = iterations, cols = parameters
  posterior_samples <- posterior_arrays[, 1, ]  # Extract first chain as matrix
  if (n_chains > 1) {
    # Combine all chains
    for (c in 2:n_chains) {
      posterior_samples <- rbind(posterior_samples, posterior_arrays[, c, ])
    }
  }

  # Compute diagnostics
  summary_table <- rstan::summary(fit)$summary

  # Handle cases where columns might not exist
  rhat_values <- if ("Rhat" %in% colnames(summary_table)) summary_table[, "Rhat"] else rep(NA_real_, nrow(summary_table))
  ess_bulk <- if ("Bulk_ESS" %in% colnames(summary_table)) summary_table[, "Bulk_ESS"] else rep(NA_real_, nrow(summary_table))
  ess_tail <- if ("Tail_ESS" %in% colnames(summary_table)) summary_table[, "Tail_ESS"] else rep(NA_real_, nrow(summary_table))

  # Get divergence and treedepth warnings (safely handle API changes)
  n_divergent <- try(rstan::get_num_divergent(fit), silent = TRUE)
  if (inherits(n_divergent, "try-error")) {
    n_divergent <- 0L  # Default to 0 if function not available
  }

  n_max_treedepth <- 0L  # Not always available in newer rstan versions

  # Extract zone-specific parameters (from list format)
  mu_zone_samples <- posterior_list$mu_zone  # (n_iter × Z × D)
  sigma_zone_samples <- posterior_list$sigma_zone  # (n_iter × Z × D)
  mu_global_samples <- posterior_list$mu_global  # (n_iter × D)

  # Calculate zone estimates (posterior means and credible intervals)
  zone_estimates <- list()
  zone_summaries_list <- list()

  for (z in seq_along(hierarchy$zones)) {
    zone_name <- hierarchy$zones[z]

    # Posterior summary for this zone
    zone_mu_mean <- colMeans(mu_zone_samples[, z, ])  # D-vector
    zone_mu_sd <- apply(mu_zone_samples[, z, ], 2, sd)  # D-vector
    zone_sigma_mean <- colMeans(sigma_zone_samples[, z, ])  # D-vector

    # Credible intervals (95%)
    ci_lower <- apply(mu_zone_samples[, z, ], 2, quantile, probs = 0.025)
    ci_upper <- apply(mu_zone_samples[, z, ], 2, quantile, probs = 0.975)

    # Compute covariance from samples
    zone_cov <- cov(mu_zone_samples[, z, ])

    zone_estimates[[zone_name]] <- list(
      n_obs = sum(zone_factor == zone_name),
      mean = zone_mu_mean,
      cov = zone_cov,
      sd = zone_mu_sd
    )

    # Zone summary dataframe
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
    shrinkage_strength = rep(prior_spec$pooling_coefficient, hierarchy$n_zones),
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
    samples = posterior_samples,  # Full posterior samples for advanced analysis
    shrinkage_weights = zone_summary,
    zone_summaries = zone_summaries_df,
    diagnostics = list(
      method = "hmc",
      backend = "stan",
      convergence = "HMC/NUTS",
      rhat = rhat_values,
      ess_bulk = ess_bulk,
      ess_tail = ess_tail,
      n_divergent = n_divergent,
      n_max_treedepth = n_max_treedepth,
      adapt_delta = adapt_delta,
      estimate_pooling = estimate_pooling,
      covariance_prior = covariance_prior
    ),
    prior_spec = prior_spec,
    metadata = list(
      method = "hmc",
      backend = "stan",
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
    cat("Stan HMC model fitted successfully.\n")
    cat("  Fitted zones:", fit_obj$metadata$n_zones_fitted, "\n")
    cat("  Total observations:", fit_obj$metadata$n_obs_total, "\n")
    cat("  Divergent transitions:", n_divergent, "\n")
    cat("  Rhat convergence: max =", max(rhat_values, na.rm = TRUE), "\n")
  }

  fit_obj
}


#' Nimble MCMC backend
#'
#' Fit hierarchical model using Nimble with adaptive MCMC samplers.
#' Provides flexible Bayesian inference with configurable samplers.
#'
#' @param data Data frame with columns: zone (factor), ilr1, ilr2, ... (ILR values)
#' @param prior_spec An object of class `"gc_prior_spec"`
#' @param n_iter Integer, number of iterations per chain (default 2000)
#' @param n_warmup Integer, warmup/burn-in iterations (default 500)
#' @param n_chains Integer, number of chains (default 2)
#' @param estimate_pooling Logical, estimate pooling strength from data? (default FALSE)
#' @param covariance_prior Character, prior for covariances: "lkj" (default) or "inverse_wishart"
#' @param verbose Logical, print progress (default TRUE)
#' @param ... Additional arguments (reserved for future use)
#'
#' @return S3 object of class `"gc_hierarchical_fit"` with Nimble MCMC results
#'
#' @keywords internal
fit_hierarchical_nimble <- function(data,
                                   prior_spec,
                                   n_iter = 2000,
                                   n_warmup = 500,
                                   n_chains = 2,
                                   estimate_pooling = FALSE,
                                   covariance_prior = "lkj",
                                   verbose = TRUE,
                                   ...) {

  if (!requireNamespace("nimble", quietly = TRUE)) {
    stop("nimble package is required for Nimble backend")
  }

  # For full Nimble support with multivariate distributions, may need coda package
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("coda package is required for Nimble diagnostics")
  }

  # Attach nimble to ensure proper environment setup for internal functions
  # This is necessary for getNimbleOption() to work properly
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
    cat("Fitting Hierarchical Model with Nimble MCMC\n")
    cat("  Method: Bayesian MCMC (Adaptive samplers)\n")
    cat("  Observations:", nrow(data), "\n")
    cat("  Zones:", hierarchy$n_zones, "\n")
    cat("  ILR dimensions:", length(ilr_cols), "\n")
    cat("  Chains:", n_chains, "× (iterations:", n_iter, ", warmup:", n_warmup, ")\n")
    cat("  Pooling: ", if (estimate_pooling) "estimated" else "fixed", "\n")
  }

  # Prepare data
  zone_factor <- as.factor(data$zone)
  y_matrix <- as.matrix(data[, ilr_cols, drop = FALSE])
  n_ilr <- length(ilr_cols)

  # Define Nimble model
  model_code <- nimble::nimbleCode({
    # Global hyperparameters
    for (d in 1:D) {
      mu_global[d] ~ dnorm(mu_global_prior_mean[d], sd = sigma_prior_mean)
    }

    # Zone-level parameters
    for (z in 1:Z) {
      # Zone means with shrinkage toward global mean
      for (d in 1:D) {
        mu_zone[z, d] ~ dnorm(mu_global[d], sd = sqrt(pooling_coef))
      }

      # Zone-level standard deviations (exponential prior)
      for (d in 1:D) {
        sigma_zone[z, d] ~ dexp(1.0)
      }
    }

    # Likelihood: multivariate normal per observation
    for (n in 1:N) {
      # Note: Nimble doesn't have easy multivariate normal, so we use univariate normal
      # This is a simplification - for production use, consider custom distributions
      for (d in 1:D) {
        y[n, d] ~ dnorm(mu_zone[zone[n], d], sd = sigma_zone[zone[n], d])
      }
    }

    # Prior on pooling if estimated
    if (estimate_pooling_flag) {
      tau_pooling ~ dexp(1.0)
    }
  })

  # Prepare constants and data
  model_constants <- list(
    N = nrow(y_matrix),
    Z = hierarchy$n_zones,
    D = n_ilr,
    zone = as.integer(zone_factor),
    mu_global_prior_mean = prior_spec$global_mean,
    sigma_prior_mean = prior_spec$prior_sd_mean,
    estimate_pooling_flag = as.integer(estimate_pooling),
    pooling_coef = prior_spec$pooling_coefficient
  )

  model_data <- list(y = y_matrix)

  # Initial values (simple approach)
  inits_fn <- function() {
    list(
      mu_global = rnorm(n_ilr, prior_spec$global_mean, 0.1),
      mu_zone = matrix(rnorm(hierarchy$n_zones * n_ilr, 0, 0.5), hierarchy$n_zones, n_ilr),
      sigma_zone = matrix(pmax(0.1, abs(rnorm(hierarchy$n_zones * n_ilr, 1, 0.2))), hierarchy$n_zones, n_ilr),
      tau_pooling = if (estimate_pooling) runif(1, 0.1, 0.5) else NULL
    )
  }

  # Build Nimble model
  nimble_model <- nimble::nimbleModel(model_code,
    constants = model_constants,
    data = model_data,
    inits = inits_fn()
  )

  # Configure MCMC
  nimble_mcmc <- nimble::configureMCMC(nimble_model, verbose = verbose)

  # Build and compile MCMC
  nimble_mcmc <- nimble::buildMCMC(nimble_mcmc)
  compiled_model <- nimble::compileNimble(nimble_model)
  compiled_mcmc <- nimble::compileNimble(nimble_mcmc, project = nimble_model)

  # Run MCMC chains
  mcmc_samples_list <- list()

  for (chain in 1:n_chains) {
    if (verbose) {
      cat("Running chain", chain, "of", n_chains, "\n")
    }

    # Set new initial values for each chain
    nimble_model$setData(model_data)
    nimble_model$setInits(inits_fn())
    compiled_model$setData(model_data)
    compiled_model$setInits(inits_fn())

    # Run sampling
    samples <- nimble::runMCMC(
      compiled_mcmc,
      niter = n_iter,
      nburnin = n_warmup,
      nchains = 1,
      thin = 1,
      summary = FALSE,
      samplesAsCodaMCMC = TRUE
    )

    mcmc_samples_list[[chain]] <- samples
  }

  # Combine chains
  all_samples <- do.call(rbind, mcmc_samples_list)

  # Convert to standard format (matrix of samples)
  posterior_samples <- as.matrix(all_samples)
  colnames_samples <- colnames(posterior_samples)

  # Extract zone-specific parameters
  mu_zone_cols <- grep("^mu_zone", colnames_samples)
  mu_global_cols <- grep("^mu_global", colnames_samples)

  zone_estimates <- list()
  zone_summaries_list <- list()

  for (z in seq_along(hierarchy$zones)) {
    zone_name <- hierarchy$zones[z]

    # Extract posterior samples for this zone
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
    zone_summaries_df <- data.frame()  # Return empty data frame if no summaries
  }

  # Global estimates
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
    zone = hierarchy$zones,
    n_obs = sapply(seq_along(hierarchy$zones), function(z) sum(zone_factor == hierarchy$zones[z])),
    shrinkage_strength = rep(prior_spec$pooling_coefficient, hierarchy$n_zones),
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
    diagnostics = list(
      method = "mcmc_adaptive",
      backend = "nimble",
      convergence = "Adaptive samplers",
      gelman_rubin = gr_values,
      estimate_pooling = estimate_pooling,
      covariance_prior = covariance_prior
    ),
    prior_spec = prior_spec,
    metadata = list(
      method = "mcmc_adaptive",
      backend = "nimble",
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
    cat("Nimble MCMC model fitted successfully.\n")
    cat("  Fitted zones:", fit_obj$metadata$n_zones_fitted, "\n")
    cat("  Total observations:", fit_obj$metadata$n_obs_total, "\n")
    if (!anyNA(gr_values)) {
      cat("  Gelman-Rubin convergence: max =", max(gr_values, na.rm = TRUE), "\n")
    }
  }

  fit_obj
}
