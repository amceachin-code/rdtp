# ===============================================================================
# rdtp_search_unit.R -- Internal per-unit cutoff search
# ===============================================================================
#
# For a single unit's data, sweeps over candidate cutoffs on the discrete
# forcing variable, fits a linear-spline RD regression at each, and selects
# the cutoff with the maximum R-squared.
#
# Model at each candidate cutoff c:
#   depvar = beta0 + beta1 * cut + beta2 * force + beta3 * (cut * force) + eps
#   where cut = 1(forcing >= c), force = forcing - c
#
# This is an INTERNAL function -- not exported.
# ===============================================================================


# -----------------------------------------------------------------------------
# .rdtp_search_unit()
#
# Arguments:
#   unit_data   -- data.frame containing depvar, forcing, and optionally
#                 the cluster variable for this unit only
#   depvar      -- character: name of the dependent variable column
#   forcing     -- character: name of the forcing variable column
#   bandwidth   -- numeric: symmetric bandwidth around each candidate
#   searchrange -- numeric(2) or NULL: restrict candidates to [min, max]
#   minobs      -- integer: minimum observations required on each side
#   vce         -- character: "ols", "robust", or "cluster"
#   cluster     -- character or NULL: name of the cluster variable column
#   level       -- numeric: confidence level (0, 1)
#
# Returns:
#   A named list with components:
#     found      -- logical: was a valid cutoff found?
#     cutoff     -- numeric: best cutoff (NA if not found)
#     r2         -- numeric: R-squared at best cutoff (NA if not found)
#     beta       -- numeric: discontinuity estimate (coef on "cut")
#     se         -- numeric: standard error of beta
#     tstat      -- numeric: t-statistic (beta / se)
#     n_left     -- integer: obs on left side at best cutoff
#     n_right    -- integer: obs on right side at best cutoff
#     n_total    -- integer: total obs in best regression
#     pred_left  -- numeric: predicted outcome at cutoff from left
#     pred_right -- numeric: predicted outcome at cutoff from right
# -----------------------------------------------------------------------------

# Note: `level` is accepted for forward compatibility (e.g., future CI
# computation) but is not currently used in the search algorithm.
.rdtp_search_unit <- function(unit_data, depvar, forcing, bandwidth,
                              searchrange, minobs, vce, cluster, level) {

    # -- Default return for "no valid cutoff found" -----------------------
    empty_result <- list(
        found      = FALSE,
        cutoff     = NA_real_,
        r2         = NA_real_,
        beta       = NA_real_,
        se         = NA_real_,
        tstat      = NA_real_,
        n_left     = NA_integer_,
        n_right    = NA_integer_,
        n_total    = NA_integer_,
        pred_left  = NA_real_,
        pred_right = NA_real_
    )

    # -- Extract vectors --------------------------------------------------
    y <- unit_data[[depvar]]
    x <- unit_data[[forcing]]

    # -- Identify candidate cutoff values ---------------------------------
    # Use sorted unique values of the forcing variable
    candidates <- sort(unique(x))

    # Filter by searchrange if specified
    if (!is.null(searchrange)) {
        candidates <- candidates[candidates >= searchrange[1] &
                                 candidates <= searchrange[2]]
    }

    # No candidates available
    if (length(candidates) == 0) {
        return(empty_result)
    }

    # -- Search over candidate cutoffs ------------------------------------
    # Track the best (highest R-squared) cutoff for this unit
    best_r2     <- -1
    best_cutoff <- NA_real_
    best_beta   <- NA_real_
    best_se     <- NA_real_
    best_tstat  <- NA_real_
    best_nleft  <- NA_integer_
    best_nright <- NA_integer_
    best_ntotal <- NA_integer_
    best_predL  <- NA_real_
    best_predR  <- NA_real_
    found_any   <- FALSE

    for (c_val in candidates) {

        # .. Construct the RD variables ..
        cut   <- as.integer(x >= c_val)
        force <- x - c_val
        int   <- cut * force

        # .. Identify observations within bandwidth ..
        in_bw <- abs(force) <= bandwidth

        # .. Count observations on each side within bandwidth ..
        nl <- sum(cut == 0L & in_bw)
        nr <- sum(cut == 1L & in_bw)

        # Skip if either side has too few observations
        if (nl < minobs || nr < minobs) next

        # .. Build regression data (within bandwidth only) ..
        reg_idx <- which(in_bw)
        reg_df  <- data.frame(
            y     = y[reg_idx],
            cut   = cut[reg_idx],
            force = force[reg_idx],
            int   = int[reg_idx]
        )

        # .. Fit the linear-spline RD regression ..
        #    y = beta0 + beta1 * cut + beta2 * force + beta3 * int + eps
        fit <- tryCatch(
            lm(y ~ cut + force + int, data = reg_df),
            error = function(e) NULL
        )

        # Skip if regression failed
        if (is.null(fit)) next

        # .. Extract key quantities ..
        cf <- coef(fit)

        # Skip if cut coefficient is NA (collinearity) or exactly 0
        # (Stata drops collinear variables and sets coef to 0)
        if (is.na(cf["cut"]) || cf["cut"] == 0) next

        # Cache summary(fit) -- avoid recomputing in the hot loop
        fit_summary <- summary(fit)

        # Get F-statistic; skip if NA (model not identified)
        fstat <- fit_summary$fstatistic
        if (is.null(fstat) || is.na(fstat[1])) next

        # .. Get R-squared (always from OLS, regardless of VCE) ..
        r2 <- fit_summary$r.squared

        # .. Update best if this R-squared is the new maximum ..
        #    Tie-breaking: strict > means the first candidate (lowest
        #    forcing value, since candidates are sorted) wins ties.
        if (r2 > best_r2) {

            # Get SE for the cut coefficient using specified VCE
            se_cut <- tryCatch({
                if (vce == "cluster") {
                    clust_vec <- unit_data[[cluster]][reg_idx]
                    .rdtp_get_se(fit, vce = vce, cluster_vec = clust_vec)
                } else {
                    .rdtp_get_se(fit, vce = vce)
                }
            }, error = function(e) {
                # Fall back to OLS SE if VCE computation fails for this cutoff
                fit_summary$coefficients["cut", "Std. Error"]
            })

            best_r2     <- r2
            best_cutoff <- c_val
            best_beta   <- cf["cut"]
            best_se     <- se_cut
            best_tstat  <- cf["cut"] / se_cut
            best_nleft  <- nl
            best_nright <- nr
            best_ntotal <- nrow(reg_df)

            # Predicted outcome at cutoff from the LEFT (cut=0, force=0):
            #   y_hat_left = beta0 = intercept
            best_predL <- cf["(Intercept)"]

            # Predicted outcome at cutoff from the RIGHT (cut=1, force=0):
            #   y_hat_right = beta0 + beta1 = intercept + coef(cut)
            best_predR <- cf["(Intercept)"] + cf["cut"]

            found_any <- TRUE
        }
    }

    # -- Return results ---------------------------------------------------
    if (!found_any) {
        return(empty_result)
    }

    list(
        found      = TRUE,
        cutoff     = unname(best_cutoff),
        r2         = unname(best_r2),
        beta       = unname(best_beta),
        se         = unname(best_se),
        tstat      = unname(best_tstat),
        n_left     = as.integer(best_nleft),
        n_right    = as.integer(best_nright),
        n_total    = as.integer(best_ntotal),
        pred_left  = unname(best_predL),
        pred_right = unname(best_predR)
    )
}
