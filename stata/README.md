# rdtp (Stata)

RD Tipping-Point Cutoff Search -- Stata command

## Requirements

- **Stata 16.0** or later

## Installation

Copy `rdtp.ado` and `rdtp.sthlp` to any directory on your Stata adopath, or
point Stata to this folder:

```stata
adopath + "/path/to/rd_tipping_point/stata"
```

To verify:

```stata
help rdtp
```

## Files

| File | Description |
|------|-------------|
| `rdtp.ado` | Main Stata command (612 lines) |
| `rdtp.sthlp` | Stata help file (viewable with `help rdtp`) |
| `tests/test_rdtp.do` | Test suite (18 tests) |

## Syntax

```
rdtp depvar forcingvar [if] [in], by(varname) [options]
```

### Required

| Argument | Description |
|----------|-------------|
| `depvar` | Dependent (outcome) variable |
| `forcingvar` | Discrete forcing (running) variable |
| `by(varname)` | Unit identifier -- search is performed separately for each unique value |

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `bandwidth(#)` | 75 | Symmetric bandwidth around each candidate cutoff |
| `searchrange(# #)` | all unique values | Restrict candidate cutoffs to [min, max] |
| `minobs(#)` | 10 | Minimum observations required on each side |
| `vce(vcetype)` | OLS | Variance estimator (e.g., `robust`, `cluster varname`) |
| `saving(filename)` | -- | Save per-unit results to a `.dta` file |
| `replace` | -- | Allow `saving()` to overwrite an existing file |
| `noisily` | -- | Display per-unit progress during the search |
| `level(#)` | 95 | Confidence level for internal regressions |

## Quick Examples

```stata
* Basic usage
rdtp math_score test_score, by(school_id)

* Restrict search range + cluster-robust SEs
rdtp math_score test_score, by(school_id) ///
    searchrange(225 500) bandwidth(50) vce(cluster test_score)

* Save results
rdtp math_score test_score, by(school_id) ///
    saving(rd_results) replace noisily
```

## Stored Results

`rdtp` is an r-class command. Key stored results:

| Name | Description |
|------|-------------|
| `r(n_units)` | Number of units searched |
| `r(n_found)` | Number of units with a valid cutoff |
| `r(mean_r2)` | Mean R-squared across units with valid cutoffs |
| `r(mean_beta)` | Mean discontinuity estimate |
| `r(results)` | n_units x 10 matrix of per-unit results |

See `help rdtp` for the complete list.

## References

- McEachin, A., T. Domina, and A. Penner. 2020. "Heterogeneous Effects of
  Early Algebra across California Middle Schools." *Journal of Policy Analysis
  and Management* 39(3): 772--800.
