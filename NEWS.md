# geocoda (development version)

## New Features - Phase 2: Multi-Backend Hierarchical Models

### MCMC Backends for Hierarchical Models

**Implemented Stan HMC Backend:**

- Hamiltonian Monte Carlo with NUTS (No-U-Turn Sampler) via rstan
- Full Bayesian posterior inference with real convergence diagnostics
- Diagnostic metrics: Rhat (potential scale reduction), ESS (effective sample size), divergences
- Flexible pooling strength: fixed (from prior) or estimated from data
- Flexible covariance priors: LKJ (correlation + variance decomposition) or Inverse-Wishart
- Posterior samples for advanced downstream analysis

**Implemented Nimble MCMC Backend:**

- Adaptive MCMC samplers with automatic configuration via nimble
- Gelman-Rubin convergence diagnostics (potential scale reduction factor)
- Multi-chain execution with parallel computation support
- Same flexible pooling and covariance prior options as Stan
- Posterior samples matrix for full Bayesian analysis

**Backend Interface:**

- New `backend` parameter in `gc_fit_hierarchical_model()`: `"analytical"` (default), `"stan"`, or `"nimble"`
- Default behavior unchanged: `backend="analytical"` preserves v0.2.0 analytical shrinkage method
- All backends return standardized `gc_hierarchical_fit` objects with consistent structure:
  - `zone_estimates`: posterior means, SDs, credible intervals per zone
  - `global_estimates`: global hyperparameter estimates
  - `samples`: NULL for analytical, posterior samples matrix for MCMC backends
  - `diagnostics`: backend-specific convergence metrics
  - `metadata`: fitting information and timestamps

**New Parameters for MCMC Backends:**

- `n_iter`: Total number of MCMC iterations (default 2000)
- `n_warmup`: Burn-in iterations (default 500)
- `n_chains`: Number of parallel chains (default 2)
- `estimate_pooling`: Estimate pooling strength from data (default FALSE for fixed)
- `covariance_prior`: "lkj" (default) or "inverse_wishart"
- Stan-specific: `adapt_delta` for sampler adaptation (default 0.8)

### Breaking Changes

- **Removed deprecated MCMC parameters** from `gc_fit_hierarchical_model()`
  - Removed: `n_iter`, `n_burnin`, `n_chains`, `adapt_delta` from direct arguments to main function
  - These parameters are now backend-specific (only used with `backend="stan"` or `backend="nimble"`)
  - Old code using these parameters will error with clear instruction to specify backend

### Documentation & Vignettes

- **New vignette**: "Hierarchical Model Backends" - comprehensive guide to all three backends
  - When to use each backend (decision tree and recommendations)
  - Performance comparison and timing benchmarks
  - Detailed examples for each backend
  - Advanced topics: estimated pooling, covariance priors, parallel computation
  - Troubleshooting common issues
- **Updated function documentation**: `gc_fit_hierarchical_model()` includes:
  - Backend selection examples
  - Parameter descriptions for MCMC backends
  - Diagnostic interpretation guidance
  - Cross-backend comparison examples

### Dependencies

**New Suggests (optional):**

- `rstan (>= 2.26.0)` - Stan HMC backend
- `nimble (>= 0.13.0)` - Nimble MCMC backend
- `coda (>= 0.19-4)` - Convergence diagnostics for Nimble

MCMC backends are optional. Analytical backend always available without additional dependencies.

### Backward Compatibility

- ✓ Default behavior unchanged: `gc_fit_hierarchical_model(data, priors)` uses analytical backend
- ✓ All existing tests continue to pass without modification
- ✓ Existing code using analytical method requires no changes
- ✗ Code using old MCMC parameters must be updated to specify backend or remove parameters

### Testing

- 11 new tests for Stan and Nimble backends in `test-hierarchical-backends.R`
- Tests properly skip when optional dependencies unavailable
- Integration test script: `tests/manual/test-backends-integration.R`
- All 95 automated tests passing, 10 skipped when MCMC packages unavailable

# geocoda 0.2.0

## Documentation Corrections

### Hierarchical Zone Methodology Clarification

- **IMPORTANT**: Corrected misleading documentation in hierarchical zone methods
- `gc_fit_hierarchical_model()` uses empirical Bayes **shrinkage estimation**, NOT full Bayesian MCMC
- Updated function documentation to accurately describe the analytical shrinkage method
- MCMC parameters (`n_iter`, `n_burnin`, `n_chains`, `adapt_delta`) now trigger deprecation warnings
  - These parameters are reserved for future MCMC implementation in v0.3.0
  - Current implementation is analytical (fast) shrinkage-based estimation
- Removed misleading "convergence diagnostics" (hardcoded Rhat/n_eff values that had no meaning)
- Return object clarification:
  - `zone_estimates` contains analytical estimates (not MCMC posterior samples)
  - `shrinkage_weights` shows pooling applied to each zone
  - `metadata$method` explicitly set to "analytical_shrinkage"
- Updated validation function to check zone coverage instead of MCMC convergence
- Added FAQ entry: "Is hierarchical modeling true Bayesian MCMC?" with roadmap to v0.3.0

## Major Features

### Documentation

- Enhanced README with workflow diagrams and examples
- New Getting Started vignette for beginners
- 5 workflow vignettes covering core functionality
- Real-world case study and FAQ vignettes
- Function examples added to key functions

### Diagnostic Functions
- Cross-validation framework (LOO + K-fold)
- Uncertainty quantification and entropy metrics
- Stationarity testing and quality reporting

### SDA/soilDB Integration
- Live SSURGO queries via soilDB package
- Component aggregation and weighting
- Tile-based optimization for large areas
- Hierarchical data preparation

### Advanced Examples
- Depth-stratified 3D soil mapping
- Multi-zone hierarchical analysis

## Technical Improvements
- Parallel query execution
- Smart caching (24-hour TTL)
- Multi-source data fusion
- Comprehensive validation framework

## Dependencies
- New: `digest`, `parallel` (imports)
- New: `soilDB` (suggests)

## Getting Started
- `vignette("Getting Started")` for beginners
- `vignette("SSURGO Integration")` for production use

