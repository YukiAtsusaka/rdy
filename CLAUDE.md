# Development Tasks for `rdy` Package

*A Practical Guide to Ranking Data Analysis in the Social Sciences*

**Reference Document:** `doc/overview.pdf`
**Last Updated:** 2026-01-28

---

## Package Overview

The `rdy` package provides ready-to-use utilities for ranking data analysis in social science research. This document outlines development tasks organized by the 9 key areas of ranking data analysis workflow.

**Companion Package:** This package works alongside `rankingQ` (https://sysilviakim.com/rankingQ) by Seo-young Silvia Kim. See integration notes throughout this document.

---

## 1. Data Collection

### High Priority
- [ ] Implement bias correction functions for measurement error
  - [ ] Plug-in estimator for anchor question method
  - [ ] Inverse probability weighting estimator
  - [ ] Add `correct_bias()` function that takes anchor column and main ranking columns
  - [ ] Documentation with examples from Atsusaka & Kim (2024)
  - [ ] Ensure compatibility with `rankingQ::imprr_direct()` and `rankingQ::imprr_weight()`

### Medium Priority
- [ ] Create helper functions for data validation
  - [ ] `validate_rankings()` - Check for valid ranking format
  - [ ] `check_ties()` - Detect and report tied rankings
  - [ ] `detect_random()` - Flag potential random responses using anchor questions

### Low Priority
- [ ] Add data collection guidelines vignette
  - [ ] Best practices for number of items (< 7)
  - [ ] Ranking format recommendations (drag-and-drop vs matrix)
  - [ ] Item order randomization strategies

---

## 2. Descriptive Statistics

### Current Status
✅ **Available in `rankingQ` package:**
- `rankingQ::imprr_direct()` - Modal rankings and empirical distributions
- `rankingQ::imprr_weight()` - Weighted ranking distributions
- These handle the univariate/bivariate descriptive analysis shown in overview (pages 8-11)

### Integration Tasks
- [ ] Document workflow between `rankingQ` and `rdy`
  - [ ] Create vignette showing how to use `rankingQ` for descriptives, then `rdy` for inference
  - [ ] Examples using `rankingQ::imprr_direct()` → `rdy::diff_in_means()`
  - [ ] Examples using `rankingQ::imprr_weight()` → `rdy::viz_diff()`

- [ ] Ensure data format compatibility
  - [ ] Test that outputs from `rankingQ` functions work seamlessly with `rdy` functions
  - [ ] Document any required data transformations
  - [ ] Add helper functions if needed for format conversion

### Unique to `rdy` (not duplicating `rankingQ`)
- [ ] Extend `viz_diff()` to accept `rankingQ` output
  - [ ] Allow visualization of descriptive stats from `imprr_direct()`
  - [ ] Overlay treatment vs control distributions

- [ ] Summary statistics specifically for experiments
  - [ ] Balance checks for ranking experiments
  - [ ] Covariate balance tables by treatment group

---

## 3. Rank Distances

### Current Status
✅ `rank_dist()` implemented with:
- Kendall distance (discordant pairs)
- Footrule distance (absolute differences)
- Rho distance (squared differences)
- Bootstrap and analytic SE methods

### Enhancement Tasks
- [ ] Add more distance metrics
  - [ ] Weighted Kendall tau
  - [ ] Cayley distance (minimum transpositions)
  - [ ] Hamming distance (position-based)

- [ ] Improve visualization
  - [ ] `plot.rank_gap()` S3 method for rank_gap objects
  - [ ] Comparison plot for multiple distance metrics

- [ ] Performance optimization
  - [ ] Consider Rcpp implementation for large datasets
  - [ ] Parallel bootstrap option

- [ ] Integration with `rankingQ`
  - [ ] Accept output from `imprr_direct()` or `imprr_weight()`
  - [ ] Compute distances between modal rankings from different groups

---

## 4. Modeling

### High Priority
- [ ] Integration with Plackett-Luce models
  - [ ] Wrapper functions for `mlogit::mlogit()`
  - [ ] Helper for data transformation (wide to long)
  - [ ] `fit_pl_model()` - Streamlined interface

- [ ] Support for mixed effects models
  - [ ] Integration with `lme4` or similar
  - [ ] Random effects for respondent-level variation

### Medium Priority
- [ ] Alternative modeling frameworks
  - [ ] Bradley-Terry models for pairwise comparisons
  - [ ] Thurstonian models
  - [ ] Distance-based models

### Low Priority
- [ ] Model diagnostics
  - [ ] Goodness-of-fit tests
  - [ ] Residual analysis for ranking models
  - [ ] Model comparison tools (AIC, BIC, etc.)

---

## 5. Drawing Quantities of Interest

### High Priority (NOTED AS "FORTHCOMING" IN OVERVIEW)
- [ ] **Implement clarify-style simulation-based inference** (Page 14-15)
  - [ ] `sim_quantities()` - Main function for simulation-based inference
  - [ ] Extract coefficients from fitted Plackett-Luce models
  - [ ] Vary target independent variable, fix others to typical values
  - [ ] Draw rankings from PL model with estimated parameters (use `rpluce()`)
  - [ ] Compute desired quantities (average ranks, probabilities, etc.)
  - [ ] Repeat 1000 times for 95% CI

- [ ] Support for various quantities of interest
  - [ ] Average ranks by covariate values
  - [ ] Predicted probabilities for pairwise comparisons
  - [ ] Top-k probabilities by covariate profiles
  - [ ] First differences and marginal effects

### Medium Priority
- [ ] Visualization of quantities of interest
  - [ ] `plot_quantities()` - Coefficient plots with uncertainties
  - [ ] Predicted probability plots
  - [ ] First difference plots

---

## 6. Sampling and Simulations

### Current Status
✅ `rpluce()` implemented for single probability vector (Page 15-16)

### Enhancement Tasks
- [ ] **Extend `rpluce()` for conditional means** (Page 16-17)
  - [ ] Support unit-specific probability vectors
  - [ ] `rpluce_conditional()` - Draw from i.i.d. PL distributions with varying parameters
  - [ ] Accept data frame with covariates and coefficients
  - [ ] Generate rankings conditional on unit characteristics

- [ ] **Support for treatment assignment** (Page 17-18)
  - [ ] `generate_experiment()` - Complete DGP for ranking experiments
  - [ ] Generate Y(0) under control probabilities
  - [ ] Define individual causal effects on item utilities
  - [ ] Generate Y(1) under treatment probabilities
  - [ ] Assign treatment: D ~ Bernoulli(p)
  - [ ] Compute observed outcomes: Y_obs = Y(1)*D + Y(0)*(1-D)

- [ ] Alternative sampling methods
  - [ ] Mallows model sampler
  - [ ] Uniform random rankings
  - [ ] Permutation-based sampling

---

## 7. Power Calculation

### Current Status
✅ `sim_power()` implemented with Monte Carlo approach (Page 19-20)

### Enhancement Tasks
- [ ] Improve `sim_power()` function
  - [ ] Add progress bar for long simulations
  - [ ] Parallel processing option
  - [ ] More flexible effect size specifications
  - [ ] Support for different experimental designs (block randomization, etc.)

- [ ] Analytical power calculations
  - [ ] Formula-based power for simple designs
  - [ ] Comparison with simulation results

- [ ] **Handle degenerate outcomes in control group** (Page 20-21)
  - [ ] Document the problem: Var(Y|D=0) = 0
  - [ ] Warning messages when SD(Y|D=0) = 0
  - [ ] Explain implications: no heterogeneity, small sample, floor/ceiling effects
  - [ ] Alternative standardization approaches
  - [ ] Sensitivity analysis for small samples

---

## 8. Missing Data

### High Priority (MAJOR GAP IN CURRENT PACKAGE)
- [ ] **Implement non-parametric bias correction** (Atsusaka & Liu, 2025) (Page 21-22)
  - [ ] `handle_partial_rankings()` - Main function for partial rankings
  - [ ] Detect partial vs complete rankings
  - [ ] Implement bias correction estimators
  - [ ] Bootstrap inference for corrected estimates
  - [ ] Document why listwise deletion fails even under MCAR

- [ ] Diagnostic functions for missing data
  - [ ] `diagnose_missing()` - Summary of missingness patterns by item and rank
  - [ ] Create tables like the example on page 21 (Alaska RCV data)
  - [ ] Visualization of missing data structure
  - [ ] Tests for MCAR assumption

- [ ] Multiple imputation approaches
  - [ ] Impute missing ranks using observed patterns
  - [ ] Integration with `mice` package
  - [ ] Proper uncertainty propagation

### Medium Priority
- [ ] Documentation and vignettes
  - [ ] Explain behavioral vs design-based missingness
  - [ ] Guide for choosing missing data method
  - [ ] Examples from real datasets (Alaska RCV data)

---

## 9. Extensions

### High Priority
- [ ] Support for weighted data
  - [ ] Survey weights in `diff_in_means()`
  - [ ] Weighted rank distances
  - [ ] Integration with `survey` package
  - [ ] Coordinate with `rankingQ::imprr_weight()`

- [ ] Regression adjustment
  - [ ] Covariate adjustment in treatment effect estimation
  - [ ] Lin's estimator for ranking outcomes
  - [ ] Doubly robust estimation

### Medium Priority
- [ ] Heterogeneous treatment effects
  - [ ] Subgroup analysis
  - [ ] Interaction terms
  - [ ] CATE estimation

- [ ] Sensitivity analysis
  - [ ] Robustness to measurement error
  - [ ] Sensitivity to missing data assumptions
  - [ ] Placebo tests and falsification

---

## Package Infrastructure

### Documentation
- [ ] Create comprehensive vignettes
  - [ ] **"Introduction to rdy and rankingQ"** - Joint workflow, division of labor
  - [ ] "Treatment Effect Estimation for Rankings" - Deep dive on `diff_in_means()`
  - [ ] "Power Analysis for Ranking Experiments" - Guide to `sim_power()`
  - [ ] "Handling Missing Rankings" - Missing data methods
  - [ ] "Advanced: Plackett-Luce Models" - Integration with modeling

- [ ] Improve function documentation
  - [ ] Add more examples to all exported functions
  - [ ] Cross-reference `rankingQ` functions where appropriate
  - [ ] Include references to methodological papers

### Testing
- [ ] Create comprehensive test suite
  - [ ] Unit tests for all exported functions
  - [ ] Integration tests for workflows
  - [ ] Tests for edge cases (degenerate outcomes, etc.)
  - [ ] Tests for compatibility with `rankingQ` output
  - [ ] Use `testthat` framework

- [ ] Add example datasets
  - [ ] Congressional representation data (from overview pages 4-5)
  - [ ] Alaska RCV data (with missing ranks, page 21)
  - [ ] Simulated datasets for testing

### Package Quality
- [ ] Set up continuous integration
  - [ ] GitHub Actions for R CMD check
  - [ ] Code coverage with `covr`
  - [ ] Automated testing on multiple platforms

- [ ] Improve code organization
  - [ ] Review and refactor existing functions
  - [ ] Consistent coding style (tidyverse style guide)
  - [ ] Optimize performance bottlenecks

- [ ] Prepare for CRAN submission
  - [ ] Fix any R CMD check warnings/notes
  - [ ] Write NEWS.md
  - [ ] Update README with badges
  - [ ] Create pkgdown website

---

## Integration with `rankingQ`

### Division of Labor
**`rankingQ` focuses on:**
- Descriptive statistics (`imprr_direct`, `imprr_weight`)
- Modal rankings
- Empirical distributions
- Visualization of ranking patterns

**`rdy` focuses on:**
- Treatment effect estimation (`diff_in_means`)
- Statistical inference (standard errors, p-values)
- Power analysis (`sim_power`)
- Sampling/simulation (`rpluce`)
- Rank distances (`rank_dist`)
- Missing data methods (forthcoming)
- Bias correction (forthcoming)

### Coordination Tasks
- [ ] Joint vignette with Silvia Kim
  - [ ] Show complete workflow from data → description → inference
  - [ ] Example: `rankingQ::imprr_direct()` → `rdy::diff_in_means()` → `rdy::viz_diff()`

- [ ] Ensure data format compatibility
  - [ ] Test that ranking data works seamlessly across both packages
  - [ ] Document standard format (wide vs long)
  - [ ] Create conversion helpers if needed

- [ ] Cross-package citation
  - [ ] Each package should reference the other in documentation
  - [ ] Joint paper or technical report?

---

## Methodological Papers to Implement

1. **Atsusaka (2025)** - Ballot order effects (Political Analysis 33(1):64-72)
   - ✅ Already supported via `diff_in_means()`

2. **Atsusaka & Liu (2025)** - Missing data approach (Working Paper)
   - ❌ Not yet implemented - **HIGH PRIORITY**

3. **Atsusaka & Kim (2024)** - Measurement error correction (Political Analysis pp. 1-22)
   - ❌ Not yet implemented - **HIGH PRIORITY**
   - Joint work with Silvia Kim, ensure coordination

4. **Atsusaka & Hayes (2025)** - Representation preferences
   - Use as motivating example (pages 3-5, 11-14 of overview)

---

## Priorities Summary

### Immediate (Next 2 weeks)
1. **Clarify-style simulation inference** (`sim_quantities()`) - noted as "forthcoming" on page 14
2. **Extend `rpluce()` for conditional means** - needed for simulation inference
3. Document integration with `rankingQ`
4. Create joint workflow vignette

### Short-term (Next 1-2 months)
1. **Missing data handling** (partial rankings) - major gap
2. **Measurement error correction** (anchor questions) - major gap
3. Test suite development with `rankingQ` compatibility tests
4. Example datasets from overview

### Long-term (Next 3-6 months)
1. CRAN submission preparation
2. Advanced modeling integration (Plackett-Luce wrappers)
3. Performance optimization (Rcpp for rank distances)
4. Comprehensive documentation website

---

## Notes

- Version 0.0.0.9000 indicates development package
- Current functions are working but need comprehensive testing
- Package structure is solid, main gaps are in missing features
- Strong collaboration with `rankingQ` - avoid duplication, ensure integration
- Overview document (22 pages) provides excellent roadmap for development

---

## Contact

**Package Author:** Yuki Atsusaka (yuki.atsusaka@gmail.com)
- GitHub: https://github.com/YukiAtsusaka/rdy
- Personal Website: https://atsusaka.org/

**Collaboration:** Seo-young Silvia Kim
- rankingQ: https://sysilviakim.com/rankingQ
