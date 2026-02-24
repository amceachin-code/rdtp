# ===============================================================================
# helper-data.R -- Shared test data generators (auto-sourced by testthat)
# ===============================================================================


# -----------------------------------------------------------------------------
# make_test_data()
#
# Generates a data.frame with 5 schools (100 obs each), discrete test_score
# in [200, 500], and a positive discontinuity of +20 at each school's true
# cutoff. True cutoffs: 300, 350, 400, 325, 375.
# -----------------------------------------------------------------------------

make_test_data <- function(seed = 12345) {
    set.seed(seed)
    n <- 500
    dat <- data.frame(
        school     = rep(1:5, each = 100),
        test_score = 200 + sample(0:300, n, replace = TRUE)
    )
    true_cuts <- c(300, 350, 400, 325, 375)
    dat$true_cutoff <- true_cuts[dat$school]
    dat$above   <- as.integer(dat$test_score >= dat$true_cutoff)
    dat$outcome <- 50 + 0.1 * dat$test_score + 20 * dat$above + rnorm(n, 0, 3)
    dat
}


# -----------------------------------------------------------------------------
# make_minimal_data()
#
# Minimal valid data.frame for quick input-validation tests.
# -----------------------------------------------------------------------------

make_minimal_data <- function() {
    data.frame(
        unit = rep(1, 50),
        x    = 1:50,
        y    = rnorm(50)
    )
}
