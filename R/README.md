# rdtp -- RD Tipping-Point Cutoff Search

**Version 1.0.0** | R CMD check: 0 errors, 0 warnings | 128 tests passing

An R package that searches over candidate cutoff values on a discrete forcing
variable to find, for each unit, the cutoff that maximizes R-squared from a
linear-spline regression discontinuity (RD) regression.

Based on the threshold-detection method in:

> McEachin, A., T. Domina, and A. Penner. 2020. "Heterogeneous Effects of
> Early Algebra across California Middle Schools." *Journal of Policy Analysis
> and Management* 39(3): 772--800.
> [doi:10.1002/pam.22202](https://doi.org/10.1002/pam.22202)

## Installation

```r
# From GitHub:
remotes::install_github("amceachin-code/rdtp", subdir = "R")

# Or from local source:
install.packages("path/to/rd_tipping_point/R", repos = NULL, type = "source")
```

**Dependencies:** Only base R (`stats`). Robust/cluster SEs optionally use
`sandwich` and `lmtest` (listed in `Suggests`).

## Quick Example

```r
library(rdtp)

# Generate synthetic data: 5 schools with known cutoffs
set.seed(12345)
n <- 500
dat <- data.frame(
  school     = rep(1:5, each = 100),
  test_score = 200 + sample(0:300, n, replace = TRUE)
)
true_cuts <- c(300, 350, 400, 325, 375)
dat$above <- as.integer(dat$test_score >= true_cuts[dat$school])
dat$outcome <- 50 + 0.1 * dat$test_score + 20 * dat$above + rnorm(n, 0, 3)

# Run the search
fit <- rdtp(dat, depvar = "outcome", forcing = "test_score",
            by = "school", searchrange = c(250, 450))
print(fit)
summary(fit)
plot(fit, type = "cutoffs")
```

## API

```r
rdtp(data, depvar, forcing, by,
     bandwidth = 75, searchrange = NULL, minobs = 10,
     vce = "ols", cluster = NULL, subset = NULL,
     saving = NULL, replace = FALSE,
     verbose = FALSE, level = 0.95)
```

| Argument | Default | Description |
|----------|---------|-------------|
| `data` | (required) | A data.frame containing all variables |
| `depvar` | (required) | Name of the dependent (outcome) variable |
| `forcing` | (required) | Name of the discrete forcing (running) variable |
| `by` | (required) | Name of the unit identifier variable |
| `bandwidth` | 75 | Symmetric bandwidth around each candidate cutoff |
| `searchrange` | `NULL` | Numeric vector `c(min, max)` to restrict candidates |
| `minobs` | 10 | Minimum observations required on each side |
| `vce` | `"ols"` | Variance estimator: `"ols"`, `"robust"` (HC1), or `"cluster"` |
| `cluster` | `NULL` | Name of the cluster variable (required when `vce = "cluster"`) |
| `subset` | `NULL` | Logical/numeric vector to subset rows |
| `saving` | `NULL` | File path to save results (`.rds` or `.csv`) |
| `replace` | `FALSE` | Overwrite existing file when using `saving` |
| `verbose` | `FALSE` | Display per-unit progress |
| `level` | 0.95 | Confidence level (0--1 scale) |

## Return Value

An S3 object of class `"rdtp"` with components:

| Component | Description |
|-----------|-------------|
| `results` | data.frame with columns: `unit`, `cutoff`, `r2`, `beta`, `se`, `tstat`, `n_left`, `n_right`, `n_total`, `pred_left`, `pred_right` |
| `n_units` | Number of units searched |
| `n_found` | Number of units with a valid cutoff |
| `n_skipped` | Number of units with no valid cutoff |
| `mean_r2` | Mean R-squared across found units |
| `mean_beta` | Mean discontinuity estimate |
| `mean_tstat` | Mean t-statistic |

S3 methods: `print()`, `summary()`, `plot()`.

## Comparison with Stata Version

The original Stata implementation (`rdtp.ado`) is archived in `archive/stata/`.
The R and Stata versions implement the identical algorithm -- OLS linear-spline
RD regression with R-squared maximization -- and produce numerically identical
results on the same data (within floating-point tolerance). See
`tests/cross-validation/` for the comparison script.

Key differences in interface:

| Feature | R | Stata |
|---------|---|-------|
| VCE syntax | `vce = "robust"` | `vce(robust)` |
| Cluster | `vce = "cluster", cluster = "var"` | `vce(cluster var)` |
| Confidence level | `level = 0.95` (0--1) | `level(95)` (0--100) |
| Verbose output | `verbose = TRUE` | `noisily` |
| Save results | `saving = "file.rds"` | `saving(file)` |

## Project Structure

```
rd_tipping_point/
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── README.md
├── .Rbuildignore
├── .gitignore
├── R/
│   ├── rdtp-package.R          # Package-level docs
│   ├── rdtp.R                  # Main function + S3 methods
│   ├── rdtp_search_unit.R      # Per-unit search logic (internal)
│   └── utils.R                 # Validation + VCE helpers (internal)
├── man/
│   ├── rdtp.Rd
│   ├── print.rdtp.Rd
│   ├── summary.rdtp.Rd
│   ├── print.summary.rdtp.Rd
│   ├── plot.rdtp.Rd
│   └── rdtp-package.Rd
├── tests/
│   ├── testthat.R
│   ├── testthat/
│   │   ├── test-rdtp.R               # Core functionality tests
│   │   ├── test-input-validation.R   # Edge cases & input validation
│   │   ├── test-vce.R                # Robust/cluster SE tests
│   │   └── test-exact-ols.R          # Exact OLS coefficient verification
│   └── cross-validation/
│       └── cross_validate.R          # R vs Stata comparison
├── vignettes/
│   └── getting-started.Rmd
└── archive/stata/                    # Original Stata implementation
    ├── rdtp.ado
    ├── rdtp.sthlp
    ├── README.md
    ├── PROGRESS.md
    ├── RD_algorithm_para_school_linear_spline_RR1.do
    └── tests/
        └── test_rdtp.do
```

## Limitations

- **Discrete forcing variable.** Designed for forcing variables with a finite
  number of distinct values (e.g., integer test scores).
- **R-squared is a heuristic.** Not a formal statistical test for threshold
  existence.
- **No second-stage estimation.** Identifies cutoffs only; treatment-effect
  estimation should be conducted separately.
- **Symmetric bandwidth only.** Asymmetric bandwidths are not supported.
- **No multiple-testing correction.** The "best" cutoff may be selected by
  chance when searching over many candidates.

## References

- Hansen, B. E. 2000. "Sample Splitting and Threshold Estimation."
  *Econometrica* 68(3): 575--603.
  [doi:10.1111/1468-0262.00124](https://doi.org/10.1111/1468-0262.00124)
- Lee, D. S., and T. Lemieux. 2010. "Regression Discontinuity Designs in
  Economics." *Journal of Economic Literature* 48(2): 281--355.
  [doi:10.1257/jel.48.2.281](https://doi.org/10.1257/jel.48.2.281)
- McEachin, A., T. Domina, and A. Penner. 2020. "Heterogeneous Effects of
  Early Algebra across California Middle Schools." *Journal of Policy Analysis
  and Management* 39(3): 772--800.
  [doi:10.1002/pam.22202](https://doi.org/10.1002/pam.22202)

## Author

Andrew McEachin
