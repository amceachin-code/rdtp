# rdtp -- RD Tipping-Point Cutoff Search

A Stata command that searches over candidate cutoff values on a discrete forcing variable to find, for each unit, the cutoff that maximizes R-squared from a linear-spline regression discontinuity (RD) regression.

Based on the threshold-detection method in:

> McEachin, A., T. Domina, and A. Penner. 2020. "Heterogeneous Effects of Early Algebra across California Middle Schools." *Journal of Policy Analysis and Management* 39(3): 772--800. [doi:10.1002/pam.22202](https://doi.org/10.1002/pam.22202)

## Overview

In many policy settings the true RD cutoff is unknown and varies across units (schools, districts, hospitals, etc.). For example, different schools may set different minimum test-score thresholds for course placement, but these cutoffs are not recorded in administrative data. `rdtp` automates the search: for each unit it sweeps over every unique value of the forcing variable (within an optional range), fits a linear-spline RD regression at each candidate cutoff, and selects the one with the highest R-squared.

The model estimated at each candidate cutoff *c* is:

```
depvar = b0 + b1*cut + b2*force + b3*(cut*force) + e
```

where `cut = 1(forcing >= c)` and `force = forcing - c`, using only observations within a symmetric bandwidth of *c*.

`rdtp` performs **only** the cutoff search. Second-stage treatment-effect estimation, bandwidth selection, and multiple-testing adjustments should be conducted separately once cutoffs have been identified.

## Requirements

- **Stata 16.0** or later

## Files

| File | Lines | Description |
|------|------:|-------------|
| `rdtp.ado` | 612 | Main Stata command |
| `rdtp.sthlp` | 389 | Stata help file (viewable with `help rdtp`) |
| `README.md` | -- | This file |

## Installation

Copy `rdtp.ado` and `rdtp.sthlp` to any directory on your Stata adopath, or point Stata to this folder:

```stata
adopath + "/path/to/rd_tipping_point"
```

To verify the installation:

```stata
help rdtp
```

## Syntax

```
rdtp depvar forcingvar [if] [in], by(varname) [options]
```

### Required

| Argument | Description |
|----------|-------------|
| `depvar` | Dependent (outcome) variable |
| `forcingvar` | Discrete forcing (running) variable |
| `by(varname)` | Unit identifier -- the search is performed separately for each unique value (numeric or string) |

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `bandwidth(#)` | 75 | Symmetric bandwidth around each candidate cutoff |
| `searchrange(# #)` | all unique values | Restrict candidate cutoffs to [min, max] |
| `minobs(#)` | 10 | Minimum observations required on each side of a candidate cutoff |
| `vce(vcetype)` | OLS | Variance estimator passed to `regress` (e.g., `robust`, `cluster varname`) |
| `saving(filename)` | -- | Save per-unit results to a `.dta` file |
| `replace` | -- | Allow `saving()` to overwrite an existing file |
| `noisily` | -- | Display per-unit progress during the search |
| `level(#)` | 95 | Confidence level for internal regressions |

## Quick Examples

### Basic usage

```stata
rdtp math_score test_score, by(school_id)
```

### Restrict the search range and use cluster-robust standard errors

```stata
rdtp math_score test_score, by(school_id) ///
    searchrange(225 500) bandwidth(50) vce(cluster test_score)
```

### Save results to disk with per-unit progress output

```stata
rdtp math_score test_score, by(school_id) ///
    saving(rd_results) replace noisily
```

### Post-processing the saved results

```stata
use rd_results, clear
histogram cutoff if !missing(cutoff), width(5) title("Detected Cutoffs")
summarize r2 beta tstat if !missing(cutoff), detail
```

### Access stored results programmatically

```stata
rdtp math_score test_score, by(school_id) searchrange(225 500)
display "Found cutoffs for " r(n_found) " of " r(n_units) " units"
matrix list r(results)
```

## Stored Results

`rdtp` is an r-class command. After execution it stores:

### Scalars

| Name | Description |
|------|-------------|
| `r(n_units)` | Number of units searched |
| `r(n_found)` | Number of units with a valid cutoff |
| `r(n_skipped)` | Number of units with no valid cutoff |
| `r(bandwidth)` | Bandwidth used |
| `r(minobs)` | Minimum observations per side |
| `r(mean_r2)` | Mean R-squared across units with valid cutoffs |
| `r(mean_beta)` | Mean discontinuity estimate across units with valid cutoffs |
| `r(mean_tstat)` | Mean t-statistic across units with valid cutoffs |
| `r(level)` | Confidence level |

### Macros

| Name | Description |
|------|-------------|
| `r(cmd)` | `rdtp` |
| `r(cmdline)` | Full command as typed |
| `r(depvar)` | Dependent variable name |
| `r(forcing)` | Forcing variable name |
| `r(byvar)` | Unit identifier variable name |
| `r(vce)` | VCE specification (if specified) |
| `r(searchrange)` | Search range (if specified) |
| `r(saving)` | Path to saved results file (if specified) |

### Matrix

| Name | Dimensions | Description |
|------|-----------|-------------|
| `r(results)` | n_units x 10 | Per-unit results |

Columns of `r(results)`: `cutoff`, `r2`, `beta`, `se`, `tstat`, `n_left`, `n_right`, `n_total`, `pred_left`, `pred_right`. Units for which no valid cutoff was found have missing values in all columns.

## Bugs Fixed from the Original Do-File

The original implementation (`RD_algorithm_para_school_linear_spline_RR1.do`) contained six bugs that `rdtp` corrects:

1. **-99 sentinels instead of missing values.** The original stored `-99` when no valid cutoff was found for a unit. This is indistinguishable from a legitimate data value and silently corrupts downstream summary statistics. `rdtp` uses Stata-native missing values (`.`).

2. **Uncaptured `bob` variable.** The original created a temporary variable named `bob` without `capture` or `tempvar`, causing the program to fail if a variable named `bob` already existed. `rdtp` uses `tempvar` and `capture` throughout.

3. **pctL/pctR label swap.** The original computed predicted values at the cutoff from the left and right but assigned them to the wrong variable names. `rdtp` computes `pred_left = _b[_cons]` (predicted value approaching from below) and `pred_right = _b[_cons] + _b[cut]` (predicted value approaching from above), with correct labels.

4. **Disk re-reads on every iteration.** The original re-read the full dataset from disk (`use ... , clear`) inside the unit loop, which is slow and fragile (depends on a hard-coded file path). `rdtp` uses `preserve`/`restore` to hold the dataset in memory.

5. **Incremental append instead of postfile.** The original opened, appended to, and re-saved a results file on every unit iteration, an O(n^2) pattern that becomes very slow with many units. `rdtp` uses Stata's `postfile` facility to write results in a single pass.

6. **Hard-coded paths and magic numbers.** The original hard-coded file paths, bandwidth, and search-range values. `rdtp` exposes all of these as configurable options with sensible defaults.

## Limitations

- **Discrete forcing variable.** Designed for forcing variables with a finite number of distinct values (e.g., integer test scores). With continuous forcing variables the candidate set may be extremely large.
- **R-squared is a heuristic.** Selecting the cutoff that maximizes R-squared is not a formal statistical test for the existence of a threshold.
- **No second-stage estimation.** `rdtp` identifies candidate cutoffs but does not estimate causal treatment effects.
- **Symmetric bandwidth only.** Asymmetric bandwidths are not supported.
- **No multiple-testing correction.** When searching over many candidates, the "best" cutoff may be selected by chance. Researchers should assess robustness using domain knowledge and sensitivity analyses.

## References

- Hansen, B. E. 2000. "Sample Splitting and Threshold Estimation." *Econometrica* 68(3): 575--603. [doi:10.1111/1468-0262.00124](https://doi.org/10.1111/1468-0262.00124)
- Lee, D. S., and T. Lemieux. 2010. "Regression Discontinuity Designs in Economics." *Journal of Economic Literature* 48(2): 281--355. [doi:10.1257/jel.48.2.281](https://doi.org/10.1257/jel.48.2.281)
- McEachin, A., T. Domina, and A. Penner. 2020. "Heterogeneous Effects of Early Algebra across California Middle Schools." *Journal of Policy Analysis and Management* 39(3): 772--800. [doi:10.1002/pam.22202](https://doi.org/10.1002/pam.22202)

## Author

Andrew McEachin
