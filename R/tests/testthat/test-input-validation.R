# ===============================================================================
# test-input-validation.R -- Edge cases and error messages for rdtp
# ===============================================================================
#
# Shared helpers (make_test_data, make_minimal_data) are in helper-data.R.
# ===============================================================================


# ===============================================================================
# Tests
# ═══════════════════════════════════════════════════════════════════════════════

test_that("data must be a data.frame", {
    expect_error(
        rdtp(as.matrix(make_minimal_data()), "y", "x", "unit"),
        "must be a data.frame"
    )
})

test_that("data must have rows", {
    dat <- make_minimal_data()
    expect_error(
        rdtp(dat[0, ], "y", "x", "unit"),
        "zero rows"
    )
})

test_that("depvar must name an existing numeric column", {
    dat <- make_minimal_data()
    expect_error(rdtp(dat, "nonexistent", "x", "unit"), "not found")
    dat$s <- letters[1:50]
    expect_error(rdtp(dat, "s", "x", "unit"), "must be numeric")
})

test_that("forcing must name an existing numeric column", {
    dat <- make_minimal_data()
    expect_error(rdtp(dat, "y", "nonexistent", "unit"), "not found")
    dat$s <- letters[1:50]
    expect_error(rdtp(dat, "y", "s", "unit"), "must be numeric")
})

test_that("by must name an existing column", {
    dat <- make_minimal_data()
    expect_error(rdtp(dat, "y", "x", "nonexistent"), "not found")
})

test_that("by cannot be same as depvar or forcing", {
    dat <- make_minimal_data()
    expect_error(rdtp(dat, "y", "x", "y"), "same as the dependent")
    expect_error(rdtp(dat, "y", "x", "x"), "same as the forcing")
})

test_that("bandwidth must be positive", {
    dat <- make_minimal_data()
    expect_error(rdtp(dat, "y", "x", "unit", bandwidth = 0), "positive")
    expect_error(rdtp(dat, "y", "x", "unit", bandwidth = -5), "positive")
})

test_that("minobs must be >= 1", {
    dat <- make_minimal_data()
    expect_error(rdtp(dat, "y", "x", "unit", minobs = 0), "at least 1")
})

test_that("searchrange must be valid", {
    dat <- make_minimal_data()
    expect_error(
        rdtp(dat, "y", "x", "unit", searchrange = c(50, 10)),
        "less than maximum"
    )
    expect_error(
        rdtp(dat, "y", "x", "unit", searchrange = c(50, 50)),
        "less than maximum"
    )
    expect_error(
        rdtp(dat, "y", "x", "unit", searchrange = 50),
        "length 2"
    )
})

test_that("vce must be valid", {
    dat <- make_minimal_data()
    expect_error(
        rdtp(dat, "y", "x", "unit", vce = "hc3"),
        "must be one of"
    )
})

test_that("cluster must be specified with vce = 'cluster'", {
    dat <- make_minimal_data()
    expect_error(
        rdtp(dat, "y", "x", "unit", vce = "cluster"),
        "must be specified"
    )
})

test_that("cluster variable must exist in data", {
    dat <- make_minimal_data()
    expect_error(
        rdtp(dat, "y", "x", "unit", vce = "cluster",
             cluster = "nonexistent"),
        "not found in data"
    )
})

test_that("level must be in (0, 1)", {
    dat <- make_minimal_data()
    expect_error(rdtp(dat, "y", "x", "unit", level = 0), "in \\(0, 1\\)")
    expect_error(rdtp(dat, "y", "x", "unit", level = 1), "in \\(0, 1\\)")
    expect_error(rdtp(dat, "y", "x", "unit", level = 95), "in \\(0, 1\\)")
})

test_that("saving file extension must be .rds or .csv", {
    dat <- make_minimal_data()
    expect_error(
        rdtp(dat, "y", "x", "unit", saving = "results.xlsx"),
        "rds or .csv"
    )
})

test_that("NAs in key variables are handled", {
    dat <- make_minimal_data()
    dat$y[1:5] <- NA
    dat$x[6:10] <- NA

    # Should run without error (NAs are dropped)
    fit <- rdtp(dat, depvar = "y", forcing = "x", by = "unit")
    expect_s3_class(fit, "rdtp")
})

test_that("saving to CSV works", {
    set.seed(123)
    dat <- data.frame(
        unit = rep(1, 100),
        x    = sample(0:99, 100, replace = TRUE),
        y    = rnorm(100)
    )
    dat$above <- as.integer(dat$x >= 50)
    dat$y <- dat$y + 5 * dat$above

    tmpfile <- tempfile(fileext = ".csv")
    on.exit(unlink(tmpfile), add = TRUE)

    fit <- rdtp(dat, depvar = "y", forcing = "x",
                by = "unit", saving = tmpfile, bandwidth = 50)

    expect_true(file.exists(tmpfile))
    saved <- utils::read.csv(tmpfile, stringsAsFactors = FALSE)
    expect_equal(nrow(saved), 1)
    expect_true("cutoff" %in% names(saved))
})
