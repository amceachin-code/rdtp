# ═══════════════════════════════════════════════════════════════════════════════
# test-exact-ols.R — Verify rdtp's OLS results match manual lm() calls
# ═══════════════════════════════════════════════════════════════════════════════
#
# For each unit, verify that rdtp's beta, se, r2, and predicted values
# match a manual lm() call on the same data at the same cutoff.
# This ensures the R implementation's internal OLS is correct independent
# of the Stata comparison.
# ═══════════════════════════════════════════════════════════════════════════════


test_that("rdtp results match manual lm() for each unit", {
    # ── Generate data with known cutoffs ─────────────────────────────────
    set.seed(54321)
    n <- 500
    dat <- data.frame(
        school     = rep(1:5, each = 100),
        test_score = 200 + sample(0:300, n, replace = TRUE)
    )
    true_cuts <- c(300, 350, 400, 325, 375)
    dat$true_cutoff <- true_cuts[dat$school]
    dat$above   <- as.integer(dat$test_score >= dat$true_cutoff)
    dat$outcome <- 50 + 0.1 * dat$test_score + 20 * dat$above + rnorm(n, 0, 3)

    # ── Run rdtp ─────────────────────────────────────────────────────────
    fit <- rdtp(dat, depvar = "outcome", forcing = "test_score",
                by = "school", searchrange = c(250, 450), bandwidth = 75)

    # ── Verify each unit against manual lm() ─────────────────────────────
    for (i in seq_len(nrow(fit$results))) {
        row <- fit$results[i, ]

        # Skip units with no valid cutoff
        if (is.na(row$cutoff)) next

        school_id <- as.integer(row$unit)
        c_val     <- row$cutoff

        # Subset to this school
        sdat <- dat[dat$school == school_id, ]

        # Construct RD variables at the detected cutoff
        sdat$cut   <- as.integer(sdat$test_score >= c_val)
        sdat$force <- sdat$test_score - c_val
        sdat$int   <- sdat$cut * sdat$force

        # Restrict to bandwidth
        sdat_bw <- sdat[abs(sdat$force) <= 75, ]

        # Fit the manual lm()
        manual_fit <- lm(outcome ~ cut + force + int, data = sdat_bw)
        manual_cf  <- coef(manual_fit)
        manual_se  <- summary(manual_fit)$coefficients["cut", "Std. Error"]
        manual_r2  <- summary(manual_fit)$r.squared

        # ── Assert exact match (within floating-point tolerance) ─────────

        # Beta (discontinuity estimate)
        expect_equal(
            row$beta, unname(manual_cf["cut"]),
            tolerance = 1e-10,
            label = sprintf("school %d: beta", school_id)
        )

        # Standard error
        expect_equal(
            row$se, manual_se,
            tolerance = 1e-10,
            label = sprintf("school %d: se", school_id)
        )

        # R-squared
        expect_equal(
            row$r2, manual_r2,
            tolerance = 1e-10,
            label = sprintf("school %d: r2", school_id)
        )

        # Predicted values at cutoff
        # pred_left = intercept (cut=0, force=0)
        expect_equal(
            row$pred_left, unname(manual_cf["(Intercept)"]),
            tolerance = 1e-10,
            label = sprintf("school %d: pred_left", school_id)
        )

        # pred_right = intercept + beta_cut (cut=1, force=0)
        expect_equal(
            row$pred_right,
            unname(manual_cf["(Intercept)"] + manual_cf["cut"]),
            tolerance = 1e-10,
            label = sprintf("school %d: pred_right", school_id)
        )

        # N total
        expect_equal(
            row$n_total, nrow(sdat_bw),
            label = sprintf("school %d: n_total", school_id)
        )

        # t-statistic = beta / se
        expect_equal(
            row$tstat, unname(manual_cf["cut"]) / manual_se,
            tolerance = 1e-10,
            label = sprintf("school %d: tstat", school_id)
        )
    }
})


test_that("predicted values are consistent with beta sign", {
    set.seed(12345)
    n <- 500
    dat <- data.frame(
        school     = rep(1:5, each = 100),
        test_score = 200 + sample(0:300, n, replace = TRUE)
    )
    true_cuts <- c(300, 350, 400, 325, 375)
    dat$above   <- as.integer(dat$test_score >= true_cuts[dat$school])
    dat$outcome <- 50 + 0.1 * dat$test_score + 20 * dat$above + rnorm(n, 0, 3)

    fit <- rdtp(dat, depvar = "outcome", forcing = "test_score",
                by = "school", searchrange = c(250, 450))

    for (i in seq_len(nrow(fit$results))) {
        row <- fit$results[i, ]
        if (is.na(row$cutoff)) next

        # pred_right - pred_left should equal beta
        expect_equal(
            row$pred_right - row$pred_left,
            row$beta,
            tolerance = 1e-10,
            label = sprintf("unit %s: pred_right - pred_left = beta", row$unit)
        )
    }
})
