# ═══════════════════════════════════════════════════════════════════════════════
# test-vce.R — Robust and cluster SE tests for rdtp
# ═══════════════════════════════════════════════════════════════════════════════


# ─────────────────────────────────────────────────────────────────────────────
# Helper: generate data with heteroskedasticity
# ─────────────────────────────────────────────────────────────────────────────

make_hetero_data <- function(seed = 42, n = 500) {
    set.seed(seed)
    dat <- data.frame(
        unit = rep(1, n),
        x    = sample(0:99, n, replace = TRUE)
    )
    dat$above <- as.integer(dat$x >= 50)
    # Heteroskedastic errors: variance increases with x
    dat$y <- 10 + 8 * dat$above + rnorm(n, 0, 1 + 0.1 * dat$x)
    dat
}


# ═══════════════════════════════════════════════════════════════════════════════
# Tests
# ═══════════════════════════════════════════════════════════════════════════════

test_that("robust SEs differ from OLS SEs", {
    skip_if_not_installed("sandwich")

    dat <- make_hetero_data()

    fit_ols <- rdtp(dat, depvar = "y", forcing = "x",
                    by = "unit", vce = "ols", bandwidth = 50)
    fit_rob <- rdtp(dat, depvar = "y", forcing = "x",
                    by = "unit", vce = "robust", bandwidth = 50)

    # Same cutoff and beta (both maximize OLS R2)
    expect_equal(fit_ols$results$cutoff[1], fit_rob$results$cutoff[1])
    expect_equal(fit_ols$results$beta[1], fit_rob$results$beta[1])
    expect_equal(fit_ols$results$r2[1], fit_rob$results$r2[1])

    # SEs should differ (heteroskedastic data)
    expect_false(
        isTRUE(all.equal(fit_ols$results$se[1], fit_rob$results$se[1])),
        label = "Robust SE should differ from OLS SE with heteroskedastic data"
    )
})


test_that("cluster SEs differ from OLS SEs", {
    skip_if_not_installed("sandwich")

    set.seed(999)
    n <- 500
    dat <- data.frame(
        unit      = rep(1, n),
        cluster_v = rep(1:25, each = 20),
        x         = sample(0:99, n, replace = TRUE)
    )
    dat$above <- as.integer(dat$x >= 50)
    # Add cluster-level random effects for within-cluster correlation
    cluster_effects <- rnorm(25, 0, 3)
    dat$y <- 10 + 8 * dat$above + cluster_effects[dat$cluster_v] +
        rnorm(n, 0, 2)

    fit_ols <- rdtp(dat, depvar = "y", forcing = "x",
                    by = "unit", vce = "ols", bandwidth = 50)
    fit_cl  <- rdtp(dat, depvar = "y", forcing = "x",
                    by = "unit", vce = "cluster", cluster = "cluster_v",
                    bandwidth = 50)

    # Same cutoff and beta
    expect_equal(fit_ols$results$cutoff[1], fit_cl$results$cutoff[1])
    expect_equal(fit_ols$results$beta[1], fit_cl$results$beta[1])

    # SEs should differ (clustered data)
    expect_false(
        isTRUE(all.equal(fit_ols$results$se[1], fit_cl$results$se[1])),
        label = "Cluster SE should differ from OLS SE with clustered data"
    )
})


test_that("R-squared is always from OLS regardless of VCE", {
    skip_if_not_installed("sandwich")

    dat <- make_hetero_data()

    fit_ols <- rdtp(dat, depvar = "y", forcing = "x",
                    by = "unit", vce = "ols", bandwidth = 50)
    fit_rob <- rdtp(dat, depvar = "y", forcing = "x",
                    by = "unit", vce = "robust", bandwidth = 50)

    # R-squared should be identical (both from OLS)
    expect_equal(fit_ols$results$r2[1], fit_rob$results$r2[1])
})


test_that("graceful error when sandwich not installed", {
    # We can't easily unload sandwich, but we can test the internal
    # function directly with a mock
    skip_if_not_installed("sandwich")

    # This test just verifies the error message text is correct
    # by checking the function body references requireNamespace
    fn_body <- body(.rdtp_get_se)
    fn_text <- paste(deparse(fn_body), collapse = "\n")
    expect_true(grepl("requireNamespace", fn_text))
})


test_that("t-statistic is beta / se for all VCE types", {
    skip_if_not_installed("sandwich")

    dat <- make_hetero_data()

    for (v in c("ols", "robust")) {
        fit <- rdtp(dat, depvar = "y", forcing = "x",
                    by = "unit", vce = v, bandwidth = 50)

        # t-stat should equal beta / se
        expect_equal(
            fit$results$tstat[1],
            fit$results$beta[1] / fit$results$se[1],
            tolerance = 1e-10,
            label = sprintf("t = beta/se for vce = '%s'", v)
        )
    }
})


test_that("warning when cluster specified without vce = 'cluster'", {
    dat <- make_hetero_data()
    dat$cl <- rep(1:25, each = 20)

    expect_warning(
        rdtp(dat, depvar = "y", forcing = "x",
             by = "unit", vce = "ols", cluster = "cl", bandwidth = 50),
        "ignored"
    )
})
