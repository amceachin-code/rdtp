# ===============================================================================
# test-rdtp.R -- Core test suite for rdtp, mirroring the 18 Stata tests
# ===============================================================================
#
# Each test generates synthetic data with known properties, runs rdtp(),
# and verifies that return values match expectations.
#
# Shared helpers (make_test_data, make_minimal_data) are in helper-data.R,
# which testthat auto-sources before running tests.
# ===============================================================================


# ===============================================================================
# TEST 1: Basic functionality — 5 units with known cutoffs
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 1: basic functionality with 5 units", {
    dat <- make_test_data()
    fit <- rdtp(dat, depvar = "outcome", forcing = "test_score",
                by = "school", searchrange = c(250, 450))

    expect_s3_class(fit, "rdtp")
    expect_equal(fit$n_units, 5)
    expect_equal(fit$n_found, 5)
    expect_equal(fit$n_skipped, 0)
    expect_equal(fit$bandwidth, 75)
    expect_equal(fit$minobs, 10)
    expect_true(fit$mean_r2 > 0)
    expect_true(fit$mean_r2 < 1)
    expect_true(fit$mean_beta > 0)  # discontinuity is positive (+20)
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 2: Cutoff accuracy — detected cutoffs near true values
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 2: cutoff accuracy within +/- 25", {
    dat <- make_test_data()
    fit <- rdtp(dat, depvar = "outcome", forcing = "test_score",
                by = "school", searchrange = c(250, 450))

    true_cuts <- c(300, 350, 400, 325, 375)

    # Results are sorted by unit, so school 1-5 map to rows 1-5
    for (i in 1:5) {
        expect_true(
            abs(fit$results$cutoff[i] - true_cuts[i]) <= 25,
            label = sprintf("school %d: detected=%g, true=%d",
                            i, fit$results$cutoff[i], true_cuts[i])
        )
    }
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 3: String by-variable
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 3: string by-variable works", {
    set.seed(99999)
    n <- 300
    dat <- data.frame(
        school_name = rep(c("Alpha", "Beta", "Gamma"), each = 100),
        score       = 1 + sample(0:99, n, replace = TRUE),
        stringsAsFactors = FALSE
    )
    true_c <- c(Alpha = 40, Beta = 60, Gamma = 50)
    dat$above <- as.integer(dat$score >= true_c[dat$school_name])
    dat$y <- 10 + 5 * dat$above + rnorm(n, 0, 2)

    fit <- rdtp(dat, depvar = "y", forcing = "score",
                by = "school_name", bandwidth = 50)

    expect_equal(fit$n_units, 3)
    expect_equal(fit$n_found, 3)
    expect_equal(fit$byvar, "school_name")
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 4: No valid cutoffs (minobs too high)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 4: no valid cutoffs when minobs is too high", {
    set.seed(11111)
    dat <- data.frame(
        unit = rep(1, 20),
        x    = 1:20,
        y    = rnorm(20)
    )

    # With only 20 obs and minobs=50, no candidate can qualify
    fit <- rdtp(dat, depvar = "y", forcing = "x",
                by = "unit", minobs = 50)

    expect_equal(fit$n_units, 1)
    expect_equal(fit$n_found, 0)
    expect_equal(fit$n_skipped, 1)
    expect_true(is.na(fit$mean_r2))
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 5: Single unit
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 5: single unit works", {
    set.seed(22222)
    n <- 200
    dat <- data.frame(
        unit = rep(1, n),
        x    = sample(0:99, n, replace = TRUE)
    )
    dat$above <- as.integer(dat$x >= 50)
    dat$y <- 3 + 8 * dat$above + rnorm(n, 0, 2)

    fit <- rdtp(dat, depvar = "y", forcing = "x",
                by = "unit", bandwidth = 50)

    expect_equal(fit$n_units, 1)
    expect_equal(fit$n_found, 1)
    expect_true(fit$mean_beta > 0)
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 6: searchrange() restricts candidates
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 6: searchrange restricts candidates", {
    set.seed(33333)
    n <- 200
    dat <- data.frame(
        unit = rep(1, n),
        x    = sample(0:99, n, replace = TRUE)
    )
    dat$above <- as.integer(dat$x >= 50)
    dat$y <- 3 + 8 * dat$above + rnorm(n, 0, 2)

    # Search only 60-90: true cutoff at 50 is excluded
    fit <- rdtp(dat, depvar = "y", forcing = "x",
                by = "unit", searchrange = c(60, 90), bandwidth = 50)

    expect_equal(fit$n_found, 1)
    expect_true(fit$results$cutoff[1] >= 60)
    expect_true(fit$results$cutoff[1] <= 90)
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 7: Non-integer forcing variable (GPA-like values)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 7: non-integer forcing variable", {
    set.seed(44444)
    n <- 300
    dat <- data.frame(
        unit = rep(1, n),
        gpa  = 1.0 + sample(0:29, n, replace = TRUE) / 10
    )
    dat$above <- as.integer(dat$gpa >= 2.5)
    dat$y <- 10 + 4 * dat$above + rnorm(n, 0, 2)

    fit <- rdtp(dat, depvar = "y", forcing = "gpa",
                by = "unit", searchrange = c(1.5, 3.5), bandwidth = 2)

    expect_equal(fit$n_found, 1)
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 8: vce = "robust" works
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 8: vce = 'robust' works", {
    skip_if_not_installed("sandwich")

    set.seed(55555)
    n <- 200
    dat <- data.frame(
        unit = rep(1, n),
        x    = sample(0:99, n, replace = TRUE)
    )
    dat$above <- as.integer(dat$x >= 50)
    dat$y <- 3 + 8 * dat$above + rnorm(n, 0, 2)

    fit <- rdtp(dat, depvar = "y", forcing = "x",
                by = "unit", vce = "robust", bandwidth = 50)

    expect_equal(fit$n_found, 1)
    expect_equal(fit$vce, "robust")
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 9: vce = "cluster" works
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 9: vce = 'cluster' works", {
    skip_if_not_installed("sandwich")

    set.seed(66666)
    n <- 500
    dat <- data.frame(
        school    = rep(1:5, each = 100),
        classroom = rep(1:20, each = 25),
        x         = 200 + sample(0:300, n, replace = TRUE)
    )
    true_cuts <- c(300, 350, 400, 325, 375)
    dat$true_c <- true_cuts[dat$school]
    dat$above  <- as.integer(dat$x >= dat$true_c)
    dat$y      <- 50 + 0.1 * dat$x + 15 * dat$above + rnorm(n, 0, 5)

    fit <- rdtp(dat, depvar = "y", forcing = "x",
                by = "school", vce = "cluster", cluster = "classroom",
                searchrange = c(250, 450))

    expect_equal(fit$n_found, 5)
    expect_equal(fit$vce, "cluster")
    expect_equal(fit$cluster_var, "classroom")
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 10: VCE pre-validation catches invalid cluster variable
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 10: VCE pre-validation catches bad cluster variable", {
    set.seed(77777)
    dat <- data.frame(
        unit = rep(1, 200),
        x    = sample(0:99, 200, replace = TRUE),
        y    = rnorm(200)
    )

    # Should error immediately — nonexistent cluster variable
    expect_error(
        rdtp(dat, depvar = "y", forcing = "x",
             by = "unit", vce = "cluster", cluster = "nonexistent_var"),
        "not found in data"
    )
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 11: saving() produces a valid file
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 11: saving() produces valid output", {
    set.seed(88888)
    n <- 200
    dat <- data.frame(
        school = rep(1:2, each = 100),
        x      = sample(0:99, n, replace = TRUE)
    )
    true_c <- c(40, 60)
    dat$above <- as.integer(dat$x >= true_c[dat$school])
    dat$y <- 5 + 6 * dat$above + rnorm(n, 0, 2)

    tmpfile <- tempfile(fileext = ".rds")
    on.exit(unlink(tmpfile), add = TRUE)

    fit <- rdtp(dat, depvar = "y", forcing = "x",
                by = "school", saving = tmpfile, bandwidth = 50)

    # Verify the saved file exists and has the right structure
    expect_true(file.exists(tmpfile))
    saved <- readRDS(tmpfile)
    expect_equal(nrow(saved), 2)
    expect_true(all(c("unit", "cutoff", "r2", "beta", "se", "tstat",
                      "n_left", "n_right", "n_total",
                      "pred_left", "pred_right") %in% names(saved)))
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 12: Results structure — correct dimensions and column names
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 12: results structure is correct", {
    set.seed(10101)
    n <- 300
    dat <- data.frame(
        unit = rep(1:3, each = 100),
        x    = sample(0:99, n, replace = TRUE)
    )
    dat$above <- as.integer(dat$x >= 50)
    dat$y <- 3 + 5 * dat$above + rnorm(n, 0, 2)

    fit <- rdtp(dat, depvar = "y", forcing = "x",
                by = "unit", bandwidth = 50)

    expect_equal(nrow(fit$results), 3)
    expect_equal(ncol(fit$results), 11)
    expected_cols <- c("unit", "cutoff", "r2", "beta", "se", "tstat",
                       "n_left", "n_right", "n_total",
                       "pred_left", "pred_right")
    expect_equal(names(fit$results), expected_cols)
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 13: Empty estimation sample errors gracefully
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 13: empty estimation sample errors", {
    set.seed(12121)
    dat <- data.frame(
        unit = rep(1, 100),
        x    = sample(0:99, 100, replace = TRUE),
        y    = rnorm(100)
    )

    # Use subset that selects zero observations
    expect_error(
        rdtp(dat, depvar = "y", forcing = "x",
             by = "unit", subset = dat$x > 999),
        "no observations"
    )
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 14: verbose option runs without error
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 14: verbose option works", {
    set.seed(14141)
    n <- 200
    dat <- data.frame(
        unit = rep(1:2, each = 100),
        x    = sample(0:99, n, replace = TRUE)
    )
    dat$above <- as.integer(dat$x >= 50)
    dat$y <- 3 + 5 * dat$above + rnorm(n, 0, 2)

    expect_message(
        rdtp(dat, depvar = "y", forcing = "x",
             by = "unit", bandwidth = 50, verbose = TRUE)
    )
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 15: pred_left < pred_right for positive discontinuity
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 15: pred_left < pred_right for positive jump", {
    set.seed(15151)
    n <- 500
    dat <- data.frame(
        unit = rep(1, n),
        x    = sample(0:99, n, replace = TRUE)
    )
    dat$above <- as.integer(dat$x >= 50)
    # Large positive jump: outcome increases by 20 at cutoff
    dat$y <- 10 + 20 * dat$above + rnorm(n, 0, 1)

    fit <- rdtp(dat, depvar = "y", forcing = "x",
                by = "unit", bandwidth = 50)

    expect_true(fit$results$pred_left[1] < fit$results$pred_right[1])
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 16: saving() with replace
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 16: saving with and without replace", {
    set.seed(16161)
    n <- 200
    dat <- data.frame(
        unit = rep(1, n),
        x    = sample(0:99, n, replace = TRUE)
    )
    dat$above <- as.integer(dat$x >= 50)
    dat$y <- 5 + 6 * dat$above + rnorm(n, 0, 2)

    tmpfile <- tempfile(fileext = ".rds")
    on.exit(unlink(tmpfile), add = TRUE)

    # First save — should succeed
    fit1 <- rdtp(dat, depvar = "y", forcing = "x",
                 by = "unit", saving = tmpfile, bandwidth = 50)
    expect_true(file.exists(tmpfile))

    # Second save without replace — should error
    expect_error(
        rdtp(dat, depvar = "y", forcing = "x",
             by = "unit", saving = tmpfile, bandwidth = 50),
        "already exists"
    )

    # Second save with replace — should succeed
    fit2 <- rdtp(dat, depvar = "y", forcing = "x",
                 by = "unit", saving = tmpfile, replace = TRUE,
                 bandwidth = 50)
    expect_true(file.exists(tmpfile))
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 17: subset argument restricts to selected units
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 17: subset argument works", {
    set.seed(17171)
    n <- 400
    dat <- data.frame(
        school = rep(1:4, each = 100),
        x      = sample(0:99, n, replace = TRUE)
    )
    dat$above <- as.integer(dat$x >= 50)
    dat$y <- 5 + 6 * dat$above + rnorm(n, 0, 2)

    # Restrict to schools 1 and 2 only
    fit <- rdtp(dat, depvar = "y", forcing = "x",
                by = "school", bandwidth = 50,
                subset = dat$school <= 2)

    expect_equal(fit$n_units, 2)
    expect_equal(fit$n_found, 2)
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 18: level option
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Test 18: level option is stored correctly", {
    set.seed(18181)
    n <- 200
    dat <- data.frame(
        unit = rep(1, n),
        x    = sample(0:99, n, replace = TRUE)
    )
    dat$above <- as.integer(dat$x >= 50)
    dat$y <- 5 + 6 * dat$above + rnorm(n, 0, 2)

    fit <- rdtp(dat, depvar = "y", forcing = "x",
                by = "unit", bandwidth = 50, level = 0.90)

    expect_equal(fit$n_found, 1)
    expect_equal(fit$level, 0.90)
})
