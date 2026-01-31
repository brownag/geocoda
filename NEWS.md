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

