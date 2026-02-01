// Hierarchical Bayesian Model for ILR-Transformed Compositional Data
//
// This model implements a hierarchical structure for zone-specific estimates
// of ILR-transformed compositional data with global shrinkage priors.
//
// Key parameters:
// - estimate_pooling_flag: 0 = fixed pooling_coef, 1 = estimate tau_pooling
// - covariance_prior_type: 0 = LKJ (default), 1 = Inverse-Wishart

data {
  // Dimensions
  int<lower=1> N;                      // Total number of observations
  int<lower=1> Z;                      // Number of zones
  int<lower=1> D;                      // ILR dimensions (n_components - 1)

  // Data
  int<lower=1,upper=Z> zone[N];        // Zone assignment for each observation
  matrix[N, D] y;                      // ILR observations (N x D)

  // Prior specification (from prior_spec object)
  vector[D] mu_global_prior_mean;      // Prior mean for global mean
  real<lower=0> sigma_prior_mean;      // Prior SD for mean parameters

  real<lower=0> prior_shape_variance;  // Gamma prior shape for variances
  real<lower=0> prior_rate_variance;   // Gamma prior rate for variances

  // Sampling options
  int<lower=0,upper=1> estimate_pooling_flag;     // 1 = estimate tau_pooling
  real<lower=0,upper=1> pooling_coef_fixed;       // Used if estimate_pooling_flag=0

  int<lower=0,upper=1> covariance_prior_type;     // 0 = LKJ, 1 = Inverse-Wishart

  // Regularization
  real<lower=0> eta_lkj;               // LKJ correlation concentration (1 = uniform)
}

parameters {
  // Global hyperparameters
  vector[D] mu_global;                 // Global mean for ILR coordinates

  // Zone-specific parameters
  matrix[Z, D] mu_zone;                // Zone means (Z x D)

  // Covariance structure
  vector<lower=0>[D] sigma_zone[Z];    // Zone-level standard deviations (array)

  // Conditional: either estimate pooling strength or use fixed
  real<lower=0.001> tau_pooling[estimate_pooling_flag];  // Pooling variance scaling

  // Correlation matrices (for LKJ prior on covariances)
  corr_matrix[D] Rho_zone[Z];          // Correlation matrices per zone
}

transformed parameters {
  // Zone covariance matrices (from SD and correlation)
  cov_matrix[D] Sigma_zone[Z];

  for (z in 1:Z) {
    // Construct covariance from SD and correlation: Sigma = diag(sigma) * Rho * diag(sigma)
    Sigma_zone[z] = diag_pre_multiply(sigma_zone[z], Rho_zone[z]) *
                    diag_post_multiply(Rho_zone[z], sigma_zone[z]);
  }

  // Pooling coefficient (either fixed or estimated)
  real pooling_coef = estimate_pooling_flag ? tau_pooling[1] : pooling_coef_fixed;
}

model {
  // Priors for global hyperparameters
  mu_global ~ normal(mu_global_prior_mean, sigma_prior_mean);

  // Zone-level priors with hierarchical structure
  for (z in 1:Z) {
    // Pooling toward global mean: zone means shrink toward global mean
    mu_zone[z] ~ normal(mu_global', pooling_coef);

    // Standard deviation priors (weak, half-normal)
    sigma_zone[z] ~ exponential(1.0);

    // Correlation prior
    if (covariance_prior_type == 0) {
      // LKJ prior for correlations
      Rho_zone[z] ~ lkj_corr(eta_lkj);
    } else if (covariance_prior_type == 1) {
      // For Inverse-Wishart equivalent, use LKJ with higher eta
      Rho_zone[z] ~ lkj_corr(2.0);
    }
  }

  // Prior on pooling parameter (if estimated)
  if (estimate_pooling_flag) {
    tau_pooling[1] ~ exponential(1.0);
  }

  // Likelihood: observations from zone-specific distributions
  for (n in 1:N) {
    y[n] ~ multi_normal(mu_zone[zone[n]]', Sigma_zone[zone[n]]);
  }
}

generated quantities {
  // For posterior predictive checks and diagnostics
  vector[N] log_lik;

  for (n in 1:N) {
    log_lik[n] = multi_normal_lpdf(y[n] | mu_zone[zone[n]]', Sigma_zone[zone[n]]);
  }
}
