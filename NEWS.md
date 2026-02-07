# geocoda (development version)

## Per-Zone Risk Assessment & Carbon Compliance (Phase 3.0)

### New Functions

**Per-Zone Risk Assessment:**
- `gc_risk_assessment_per_zone()` - Expected loss under asymmetric costs, stratified by zone
  - Asymmetric cost modeling (high_under vs high_over penalties)
  - Multi-component support (assess several soil properties simultaneously)
  - Zone-aware decision support for resource allocation
  - Output: Risk metrics with probability breakdown per zone

**Per-Zone Carbon Auditing:**
- `gc_carbon_audit_per_zone()` - Carbon stock estimation with compliance verification
  - SoilGrids empirical carbon model and custom carbon fraction options
  - Bootstrap-based confidence intervals
  - Compliance status: Verified / Uncertain / At Risk
  - Credit issuance with conservative adjustments

### Features

- Zone-stratified decision support following Phase 2 ensemble integration
- Flexible zone specification: zone_col, zone_definition (sf/factor/character)
- Full 39-test comprehensive test suite validating all aspects
- Examples: 05_per_zone_risk_assessment_decision_support.Rmd with complete workflow
- 800+ lines of production-ready example code

### Documentation

- Roxygen2-generated comprehensive man pages
- Example demonstrates cost matrix definition, decision matrices, sensitivity analysis

---

## Advanced Diagnostics & Uncertainty Quantification (Phase 3.1)

### New Functions

**Cross-Validation:**
- `gc_cross_validate_per_zone()` - K-fold CV metrics stratified by zone
  - Identifies under-validated zones needing additional data
  - Separate RMSE/MAE/R² per zone and fold
  - Supports k-fold or leave-one-out strategies
  - Custom model function support

**Compositional Entropy:**
- `gc_compute_entropy_per_zone()` - Shannon/Simpson entropy per zone
  - Identifies compositional uncertainty hotspots
  - Normalized to [0,1] for cross-zone comparison
  - Dominant component identification
  - Uncertainty classification (Low/Medium/High)

**Bootstrap Parameter Uncertainty:**
- `gc_bootstrap_uncertainty_per_zone()` - Parameter CI from resampling
  - Zone-specific bootstrap-based confidence intervals
  - Assesses parameter estimate stability
  - Configurable confidence levels and replications
  - Detects whether zones are meaningfully distinct

### Features

- Comprehensive 16-test scenario coverage (input validation, computation, edge cases)
- Example: 06_per_zone_model_diagnostics.Rmd with integrated diagnostic dashboard
- Composite confidence scoring combining all three metrics
- Sampling recommendation framework based on diagnostics

### Documentation

- Vignette 09: Advanced Diagnostics & Uncertainty Quantification (600+ lines)
- Diagnostic best practices and interpretation guidelines
- Usage patterns for iterative model refinement

---

## 3D Stratified Geostatistics (Phase 3.2)

### New Functions

**3D Variogram Fitting:**
- `gc_fit_vgm_3d_per_zone()` - Zone-specific 3D variogram models
  - Lateral vs vertical anisotropy quantification
  - 3D distance computation respecting domain bounds
  - Per-zone variogram summaries with range estimates
  - Anisotropy ratio calculation (lateral/vertical)

**Nested Hierarchy Definition:**
- `gc_define_nested_hierarchy()` - Multi-scale domain × depth structures
  - Lateral domains (geological/topographic) × depth strata (pedogenic)
  - Configurable depth breaks for flexible stratification
  - Creates 2D grid of hierarchical units for modeling

**3D Hierarchical Simulation:**
- `gc_sim_hierarchical_3d_per_zone()` - 3D ensemble respecting zone/depth structure
  - Generates compositional realizations on 3D prediction grid
  - Depth-dependent variability modeling (compaction effect)
  - Zone × depth hierarchical parameter support
  - Full ensemble with spatial coordinates

### Features

- 23-test comprehensive validation of 3D functions
- Examples: 
  - 07_stratified_3d_geostatistics.Rmd (theory + 3D variogram patterns)
  - 08_complete_3d_stratified_workflow.Rmd (end-to-end risk and carbon application)
- Depth profile extraction and visualization templates
- Spatial variability by depth assessment

### Documentation

- Man pages with detailed usage and conceptual explanations
- Examples include depth-stratified composition profiles and uncertainty mapping
- 3D-specific best practices (anisotropy interpretation, depth modeling)

---

## Overall Phase 3 Summary

### Total Additions

- **6 new functions** (2 risk + 3 diagnostics + 3 3D spatial)
- **4 comprehensive examples** (05-08_*.Rmd: 2200+ lines total)
- **3 test suites** (120+ tests across diagnostics, 3D, and integration)
- **2 documentation files** (Vignette 09, Complete Stratification Guide)
- **README enhancement** with stratified modeling overview

### Backward Compatibility

All Phase 3 additions are new functions; no breaking changes to existing API.
Works seamlessly with Phase 1-2 infrastructure (gc_identify_strata, HZM, ensemble).

### Integration Pattern

Phase 3.0-3.2 follows consistent per-zone suffix pattern:
```
gc_ensemble_per_zone() [Phase 2]
gc_percentile_map_per_zone() [Phase 2]
gc_probability_map_per_zone() [Phase 2]
gc_ensemble_quality_report_per_zone() [Phase 2]
gc_risk_assessment_per_zone() [Phase 3.0]
gc_carbon_audit_per_zone() [Phase 3.0]
gc_cross_validate_per_zone() [Phase 3.1]
gc_compute_entropy_per_zone() [Phase 3.1]
gc_bootstrap_uncertainty_per_zone() [Phase 3.1]
gc_fit_vgm_3d_per_zone() [Phase 3.2]
gc_define_nested_hierarchy() [Phase 3.2]
gc_sim_hierarchical_3d_per_zone() [Phase 3.2]
```

---

## Risk Assessment & Decision Support Legacy (Phase 3.0 Precursor)

### Risk Assessment Functions

**Probability Mapping:**
- `gc_probability_map()` - Calculate probability of exceeding threshold at each location
  - Supports multiple operators: >, >=, <, <=, ==
  - NA handling for robust computation
  - Output: SpatRaster with probability values in [0, 1]

**Percentile Mapping:**
- `gc_percentile_map()` - Compute percentiles (e.g., P10, P50, P90) from realizations
  - Standard uncertainty quantification framework
  - Custom percentile lists and naming formats
  - Conservative/optimistic bounds for decision-making

**Risk Assessment with Loss Functions:**
- `gc_risk_assessment()` - Integrate asymmetric costs into decision framework
  - Linear (symmetric), asymmetric (two-sided), step (binary) loss functions
  - Custom loss function support for advanced users
  - Expected loss maps for optimal decision boundaries
  - Cost parameters: a (false positive cost), b (false negative cost)

**Carbon Credit Auditing:**
- `gc_carbon_audit()` - Specialized wrapper for carbon stock verification
  - Conservative, expected, optimistic audit types
  - Confidence-based decision thresholds
  - Human-readable audit reports with summary statistics
  - Works with spatial (SpatRaster) and non-spatial (data.frame) data
  - Integration with carbon credit programs (payment-for-ecosystem-services)

### Documentation & Vignettes

Vignette coverage includes:
- Ensemble Analysis & Risk Assessment (vignette 02)
- Complex workflow and advanced examples (vignette 00, 10)
- Parameter selection guidelines (vignette 06)
- FAQ and troubleshooting (vignette 07)

### Features

- Post-processes simulation realizations for decision support
- No changes to existing simulation functions or output structures
- Works seamlessly with `gc_sim_composition()`, `gc_simulate_zones()`, `gc_sim_hierarchical()`
- Compatible with all hierarchical backends (analytical, Stan, Nimble)
- No new external dependencies added

### Backward Compatibility

- [Yes] Fully backward compatible
- [Yes] All existing tests continue to pass
- [Yes] Risk assessment functions are opt-in (new functions, no modifications to existing)
- [Yes] Existing code requires zero modifications

### Testing

- 144 automated tests for risk assessment functions
- Unit tests: probability mapping, percentile calculation, loss functions
- Integration tests: with existing geocoda functions
- Validation tests: against theoretical probability distributions
- Edge case handling: missing values, extreme thresholds, custom loss functions

---

## Multi-Backend Hierarchical Models

### MCMC Backends

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

- Vignette 08: "Hierarchical Model Backends: Deep Dive"
  - When to use each backend (analytical, Stan, Nimble)
  - Performance benchmarks and examples
  - Troubleshooting backend-specific issues
- Function documentation: `gc_fit_hierarchical_model()` with backend selection examples

### Dependencies

**New Suggests (optional):**

- `rstan (>= 2.26.0)` - Stan HMC backend
- `nimble (>= 0.13.0)` - Nimble MCMC backend
- `coda (>= 0.19-4)` - Convergence diagnostics for Nimble

MCMC backends are optional. Analytical backend always available without additional dependencies.

### Backward Compatibility

- [Yes] Default behavior unchanged: `gc_fit_hierarchical_model(data, priors)` uses analytical backend
- [Yes] All existing tests continue to pass without modification
- [Yes] Existing code using analytical method requires no changes
- [No] Code using old MCMC parameters must be updated to specify backend or remove parameters

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

- Comprehensive vignette collection covering:
  - Soil texture workflow and complete examples
  - Ensemble analysis and risk assessment
  - Hierarchical and 3D spatial modeling
  - SSURGO integration  
  - Parameter selection and troubleshooting
- Function-level documentation with examples

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

See vignettes for detailed examples:
- Soil Texture Workflow (vignette 00) - Complete beginner workflow
- SSURGO Integration (vignette 05) - Production data integration
- FAQ & Troubleshooting (vignette 07) - Common issues and solutions

