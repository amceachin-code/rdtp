# rdtp

RD Tipping-Point Cutoff Search for Discrete Forcing Variables

Searches over candidate cutoff values on a discrete forcing variable to find,
for each unit, the cutoff that maximizes R-squared from a linear-spline
regression discontinuity (RD) regression. Based on the threshold-detection
method in McEachin, Domina & Penner (2020, JPAM).

## Implementations

| Language | Folder | Package Name | Install |
|----------|--------|-------------|---------|
| **R** | [`R/`](R/) | `rdtp` | `remotes::install_github("amceachin-code/rdtp", subdir = "R")` |
| **Stata** | [`stata/`](stata/) | `rdtp` | Copy files to adopath |

Both implementations use the same algorithm and produce numerically equivalent
results on identical data.

## Method

For each unit (identified by a grouping variable), `rdtp` sweeps over every
unique value of the forcing variable within an optional search range. At each
candidate cutoff *c*, it fits:

```
depvar = b0 + b1*cut + b2*force + b3*(cut*force) + e
```

where `cut = 1(forcing >= c)` and `force = forcing - c`, using observations
within a symmetric bandwidth of *c*. The cutoff with the highest R-squared is
selected.

## Limitations

- Designed for discrete forcing variables with a finite number of distinct values.
- R-squared maximization is a heuristic, not a formal statistical test.
- Performs only the cutoff search -- second-stage treatment-effect estimation
  should be conducted separately.
- Symmetric bandwidth only.

## References

- Hansen, B. E. 2000. "Sample Splitting and Threshold Estimation."
  *Econometrica* 68(3): 575--603.
- Lee, D. S., and T. Lemieux. 2010. "Regression Discontinuity Designs in
  Economics." *Journal of Economic Literature* 48(2): 281--355.
- McEachin, A., T. Domina, and A. Penner. 2020. "Heterogeneous Effects of
  Early Algebra across California Middle Schools." *Journal of Policy Analysis
  and Management* 39(3): 772--800.

## License

GPL (>= 3)
