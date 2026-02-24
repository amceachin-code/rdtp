# ===============================================================================
# rdtp.R -- Main function and S3 methods for the rdtp package
# ===============================================================================
#
# RD Tipping-Point Search for Discrete Forcing Variables
#
# Searches over candidate cutoff values on a discrete running variable to
# find, for each unit, the cutoff that maximizes R-squared from a
# linear-spline regression discontinuity (RD) regression.
#
# Based on McEachin, Domina, & Penner (2020, JPAM)
# doi:10.1002/pam.22202
# ===============================================================================


# -----------------------------------------------------------------------------
# rdtp() -- Main exported function
# -----------------------------------------------------------------------------

#' RD Tipping-Point Cutoff Search for Discrete Forcing Variables
#'
#' Searches over candidate cutoff values on a discrete running variable to
#' find, for each unit, the cutoff that maximizes R-squared from a
#' linear-spline regression discontinuity (RD) regression.
#'
#' @param data A data.frame containing all variables.
#' @param depvar Character string naming the dependent (outcome) variable.
#' @param forcing Character string naming the discrete forcing (running) variable.
#' @param by Character string naming the unit identifier variable. The search
#'   is performed separately for each unique value. May be numeric or character.
#' @param bandwidth Numeric. Symmetric bandwidth around each candidate cutoff.
#'   Only observations with \code{|forcing - c| <= bandwidth} are included in
#'   the regression. Default is 75.
#' @param searchrange Numeric vector of length 2 \code{c(min, max)}, or
#'   \code{NULL}. If specified, restricts candidate cutoffs to values of the
#'   forcing variable in \code{[min, max]}. Default is \code{NULL} (all unique
#'   values).
#' @param minobs Integer. Minimum number of observations required on each side
#'   of a candidate cutoff (within the bandwidth) for the regression to be
#'   attempted. Default is 10.
#' @param vce Character string specifying the variance-covariance estimator.
#'   One of \code{"ols"} (default), \code{"robust"} (HC1, matching Stata), or
#'   \code{"cluster"}.
#' @param cluster Character string naming the cluster variable (required when
#'   \code{vce = "cluster"}, ignored otherwise).
#' @param subset Logical or numeric vector for subsetting rows of \code{data}
#'   before the search. Default is \code{NULL} (use all rows).
#' @param saving Character string. If specified, saves the per-unit results to
#'   the given file path. Supported extensions: \code{.rds} (default) and
#'   \code{.csv}.
#' @param replace Logical. If \code{TRUE}, overwrite an existing file when
#'   using \code{saving}. Default is \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, display per-unit progress during
#'   the search. Default is \code{FALSE}.
#' @param level Numeric in \code{(0, 1)}. Confidence level for reporting.
#'   Default is 0.95.
#'
#' @return An object of class \code{"rdtp"}, which is a list containing:
#' \describe{
#'   \item{results}{A data.frame with one row per unit and columns: \code{unit},
#'     \code{cutoff}, \code{r2}, \code{beta}, \code{se}, \code{tstat},
#'     \code{n_left}, \code{n_right}, \code{n_total}, \code{pred_left},
#'     \code{pred_right}. Units with no valid cutoff have \code{NA} values.}
#'   \item{n_units}{Number of units searched.}
#'   \item{n_found}{Number of units with a valid cutoff.}
#'   \item{n_skipped}{Number of units with no valid cutoff.}
#'   \item{mean_r2}{Mean R-squared across found units (NA if n_found == 0).}
#'   \item{mean_beta}{Mean discontinuity estimate across found units.}
#'   \item{mean_tstat}{Mean t-statistic across found units.}
#'   \item{bandwidth}{Bandwidth used.}
#'   \item{minobs}{Minimum observations per side.}
#'   \item{level}{Confidence level.}
#'   \item{vce}{VCE type used.}
#'   \item{cluster_var}{Cluster variable name (or NULL).}
#'   \item{depvar}{Dependent variable name.}
#'   \item{forcing}{Forcing variable name.}
#'   \item{byvar}{Unit identifier variable name.}
#'   \item{searchrange}{Search range (or NULL).}
#'   \item{saving}{File path where results were saved (or NULL).}
#'   \item{call}{The matched call.}
#' }
#'
#' @details
#' For each unit identified by \code{by}, the algorithm:
#' \enumerate{
#'   \item Obtains the set of unique values of the forcing variable within
#'     \code{searchrange} (or all unique values if unspecified).
#'   \item For each candidate cutoff \eqn{c}, constructs:
#'     \deqn{cut = 1(forcing \ge c), \quad force = forcing - c, \quad
#'       interact = cut \times force}
#'   \item Counts observations on each side within \code{bandwidth}. Skips if
#'     either side has fewer than \code{minobs} observations.
#'   \item Estimates:
#'     \deqn{depvar = \beta_0 + \beta_1 \cdot cut + \beta_2 \cdot force +
#'       \beta_3 \cdot interact + \varepsilon}
#'   \item Selects the cutoff \eqn{c^*} that maximizes \eqn{R^2}.
#' }
#'
#' Predicted values at the cutoff:
#' \itemize{
#'   \item \code{pred_left} = \eqn{\hat{\beta}_0} (intercept)
#'   \item \code{pred_right} = \eqn{\hat{\beta}_0 + \hat{\beta}_1}
#' }
#'
#' @references
#' McEachin, A., T. Domina, and A. Penner. 2020. "Heterogeneous effects of
#' early algebra across California middle schools." \emph{Journal of Policy
#' Analysis and Management} 39(3): 772-800.
#' \doi{10.1002/pam.22202}
#'
#' Hansen, B. E. 2000. "Sample splitting and threshold estimation."
#' \emph{Econometrica} 68(3): 575-603.
#' \doi{10.1111/1468-0262.00124}
#'
#' @examples
#' # Generate synthetic data: 5 schools with known cutoffs
#' set.seed(12345)
#' n <- 500
#' dat <- data.frame(
#'   school     = rep(1:5, each = 100),
#'   test_score = 200 + sample(0:300, n, replace = TRUE)
#' )
#' true_cuts <- c(300, 350, 400, 325, 375)
#' dat$above <- as.integer(dat$test_score >= true_cuts[dat$school])
#' dat$outcome <- 50 + 0.1 * dat$test_score + 20 * dat$above + rnorm(n, 0, 3)
#'
#' # Run the tipping-point search
#' fit <- rdtp(dat, depvar = "outcome", forcing = "test_score",
#'             by = "school", searchrange = c(250, 450))
#' print(fit)
#' summary(fit)
#'
#' @export
rdtp <- function(data, depvar, forcing, by,
                 bandwidth = 75, searchrange = NULL, minobs = 10,
                 vce = "ols", cluster = NULL, subset = NULL,
                 saving = NULL, replace = FALSE,
                 verbose = FALSE, level = 0.95) {

    # ---------------------------------------------------------------------
    # 1. CAPTURE CALL AND VALIDATE INPUTS
    # ---------------------------------------------------------------------
    cl <- match.call()

    saving <- .rdtp_validate_inputs(
        data = data, depvar = depvar, forcing = forcing, by = by,
        bandwidth = bandwidth, searchrange = searchrange, minobs = minobs,
        vce = vce, cluster = cluster, subset = subset, saving = saving,
        replace = replace, level = level
    )

    # ---------------------------------------------------------------------
    # 2. APPLY SUBSET AND REMOVE NAs
    # ---------------------------------------------------------------------

    # Apply subset if provided
    if (!is.null(subset)) {
        data <- data[subset, , drop = FALSE]
        if (nrow(data) == 0) {
            stop("no observations remaining after applying 'subset'",
                 call. = FALSE)
        }
    }

    # Remove rows with NAs in depvar, forcing, or by
    keep_cols <- c(depvar, forcing, by)
    if (vce == "cluster" && !is.null(cluster)) {
        keep_cols <- c(keep_cols, cluster)
    }
    complete <- complete.cases(data[, keep_cols, drop = FALSE])
    data <- data[complete, , drop = FALSE]

    if (nrow(data) == 0) {
        stop("no observations in estimation sample (after removing NAs)",
             call. = FALSE)
    }

    # ---------------------------------------------------------------------
    # 3. ENUMERATE UNITS AND PRE-VALIDATE VCE
    # ---------------------------------------------------------------------

    units <- sort(unique(data[[by]]))
    n_units <- length(units)

    if (n_units == 0) {
        stop("no units found in 'by' variable", call. = FALSE)
    }

    # VCE pre-validation: run a test to catch invalid specifications early
    if (vce != "ols") {
        first_unit_data <- data[data[[by]] == units[1], , drop = FALSE]
        .rdtp_vce_prevalidation(vce, cluster, first_unit_data)
    }

    # ---------------------------------------------------------------------
    # 4. MAIN LOOP -- iterate over units, search over cutoffs
    # ---------------------------------------------------------------------

    # Pre-allocate results list
    results_list <- vector("list", n_units)
    n_found   <- 0L
    n_skipped <- 0L

    for (i in seq_along(units)) {
        u <- units[i]

        # Subset to this unit
        unit_mask <- data[[by]] == u
        unit_data <- data[unit_mask, , drop = FALSE]

        if (nrow(unit_data) == 0) {
            # No observations for this unit -- record as skipped
            results_list[[i]] <- list(
                unit = as.character(u), found = FALSE,
                cutoff = NA_real_, r2 = NA_real_, beta = NA_real_,
                se = NA_real_, tstat = NA_real_,
                n_left = NA_integer_, n_right = NA_integer_,
                n_total = NA_integer_,
                pred_left = NA_real_, pred_right = NA_real_
            )
            n_skipped <- n_skipped + 1L
            if (verbose) {
                message(sprintf("  [%d/%d] %s: no observations -- skipped",
                                i, n_units, as.character(u)))
            }
            next
        }

        # Run the per-unit search
        res <- .rdtp_search_unit(
            unit_data   = unit_data,
            depvar      = depvar,
            forcing     = forcing,
            bandwidth   = bandwidth,
            searchrange = searchrange,
            minobs      = minobs,
            vce         = vce,
            cluster     = cluster,
            level       = level
        )

        # Record results
        res$unit <- as.character(u)
        results_list[[i]] <- res

        if (res$found) {
            n_found <- n_found + 1L
            if (verbose) {
                message(sprintf(
                    "  [%d/%d] %s: cutoff = %g, R2 = %.4f, beta = %.3f",
                    i, n_units, as.character(u),
                    res$cutoff, res$r2, res$beta
                ))
            }
        } else {
            n_skipped <- n_skipped + 1L
            if (verbose) {
                message(sprintf("  [%d/%d] %s: no valid cutoff -- skipped",
                                i, n_units, as.character(u)))
            }
        }

        # Progress display (every 50 units, unless verbose is on)
        if (!verbose && i %% 50 == 0 && i < n_units) {
            message(sprintf("  ... processed %d of %d units", i, n_units))
        }
    }

    # ---------------------------------------------------------------------
    # 5. ASSEMBLE RESULTS DATA.FRAME
    # ---------------------------------------------------------------------

    results <- data.frame(
        unit       = vapply(results_list, function(r) r$unit, character(1)),
        cutoff     = vapply(results_list, function(r) r$cutoff, numeric(1)),
        r2         = vapply(results_list, function(r) r$r2, numeric(1)),
        beta       = vapply(results_list, function(r) r$beta, numeric(1)),
        se         = vapply(results_list, function(r) r$se, numeric(1)),
        tstat      = vapply(results_list, function(r) r$tstat, numeric(1)),
        n_left     = vapply(results_list, function(r) r$n_left, integer(1)),
        n_right    = vapply(results_list, function(r) r$n_right, integer(1)),
        n_total    = vapply(results_list, function(r) r$n_total, integer(1)),
        pred_left  = vapply(results_list, function(r) r$pred_left, numeric(1)),
        pred_right = vapply(results_list, function(r) r$pred_right, numeric(1)),
        stringsAsFactors = FALSE
    )

    # ---------------------------------------------------------------------
    # 6. SUMMARY STATISTICS
    # ---------------------------------------------------------------------

    if (n_found > 0) {
        found_mask <- !is.na(results$cutoff)
        mean_r2    <- mean(results$r2[found_mask])
        mean_beta  <- mean(results$beta[found_mask])
        mean_tstat <- mean(results$tstat[found_mask])
    } else {
        mean_r2    <- NA_real_
        mean_beta  <- NA_real_
        mean_tstat <- NA_real_
    }

    # ---------------------------------------------------------------------
    # 7. SAVE RESULTS IF REQUESTED
    # ---------------------------------------------------------------------

    if (!is.null(saving)) {
        ext <- tolower(tools::file_ext(saving))
        if (ext == "csv") {
            utils::write.csv(results, file = saving, row.names = FALSE)
        } else {
            saveRDS(results, file = saving)
        }
    }

    # ---------------------------------------------------------------------
    # 8. BUILD AND RETURN S3 OBJECT
    # ---------------------------------------------------------------------

    structure(
        list(
            results     = results,
            n_units     = n_units,
            n_found     = n_found,
            n_skipped   = n_skipped,
            mean_r2     = mean_r2,
            mean_beta   = mean_beta,
            mean_tstat  = mean_tstat,
            bandwidth   = bandwidth,
            minobs      = minobs,
            level       = level,
            vce         = vce,
            cluster_var = cluster,
            depvar      = depvar,
            forcing     = forcing,
            byvar       = by,
            searchrange = searchrange,
            saving      = saving,
            call        = cl
        ),
        class = "rdtp"
    )
}


# ===============================================================================
# S3 METHODS
# ===============================================================================


# -----------------------------------------------------------------------------
# print.rdtp()
# -----------------------------------------------------------------------------

#' Print method for rdtp objects
#'
#' Displays a summary header, per-unit results table (if <= 50 units), and
#' aggregate statistics.
#'
#' @param x An object of class \code{"rdtp"}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#' @export
print.rdtp <- function(x, ...) {

    # -- Header -----------------------------------------------------------
    cat("\n")
    cat("RD Tipping-Point Search\n")
    cat(strrep("-", 60), "\n")
    cat(sprintf("Dependent variable:  %s\n", x$depvar))
    cat(sprintf("Forcing variable:    %s\n", x$forcing))
    cat(sprintf("Unit variable:       %s\n", x$byvar))
    cat(sprintf("Bandwidth:           %g\n", x$bandwidth))
    cat(sprintf("Min obs per side:    %d\n", x$minobs))
    if (x$vce == "ols") {
        cat("VCE:                 OLS (default)\n")
    } else if (x$vce == "cluster") {
        cat(sprintf("VCE:                 cluster (%s)\n", x$cluster_var))
    } else {
        cat(sprintf("VCE:                 %s\n", x$vce))
    }
    if (!is.null(x$searchrange)) {
        cat(sprintf("Search range:        [%g, %g]\n",
                    x$searchrange[1], x$searchrange[2]))
    } else {
        cat("Search range:        (all unique values)\n")
    }
    cat(sprintf("Units to search:     %d\n", x$n_units))
    cat(strrep("-", 60), "\n")

    # -- Per-unit table (if <= 50 units and at least one found) -----------
    if (x$n_units <= 50 && x$n_found > 0) {
        cat("\n")
        cat(strrep("-", 77), "\n")
        cat(sprintf("%18s | %10s %8s %11s %10s %9s %8s\n",
                    "Unit", "Cutoff", "R2", "Beta", "SE", "t", "N"))
        cat(strrep("-", 77), "\n")

        for (i in seq_len(nrow(x$results))) {
            u_lbl <- x$results$unit[i]
            # Truncate long names for display
            if (nchar(u_lbl) > 17) {
                u_lbl <- paste0(substr(u_lbl, 1, 14), "...")
            }

            if (!is.na(x$results$cutoff[i])) {
                cat(sprintf("%18s | %10g %8.4f %11.3f %10.3f %9.3f %8d\n",
                            u_lbl,
                            x$results$cutoff[i],
                            x$results$r2[i],
                            x$results$beta[i],
                            x$results$se[i],
                            x$results$tstat[i],
                            x$results$n_total[i]))
            } else {
                cat(sprintf("%18s | %10s %8s %11s %10s %9s %8s\n",
                            u_lbl, ".", ".", ".", ".", ".", "."))
            }
        }
        cat(strrep("-", 77), "\n")
    }

    # -- Summary ----------------------------------------------------------
    cat("\n")
    cat(strrep("-", 60), "\n")
    cat("Results\n")
    cat(strrep("-", 60), "\n")
    cat(sprintf("Units searched:       %d\n", x$n_units))
    cat(sprintf("Cutoffs found:        %d\n", x$n_found))
    cat(sprintf("Units skipped/failed: %d\n", x$n_skipped))

    if (x$n_found > 0) {
        cat("\nSummary (units with valid cutoffs):\n")
        cat(sprintf("  Mean R2:            %9.4f\n", x$mean_r2))
        cat(sprintf("  Mean discontinuity: %9.4f\n", x$mean_beta))
        cat(sprintf("  Mean t-statistic:   %9.4f\n", x$mean_tstat))
    }

    if (x$n_units > 50 && x$n_found > 0) {
        cat(sprintf(
            "\n(Per-unit table suppressed: %d units > 50. Use saving() to inspect individual results.)\n",
            x$n_units))
    }

    if (!is.null(x$saving)) {
        cat(sprintf("\nResults saved to: %s\n", x$saving))
    }

    cat("\n")
    invisible(x)
}


# -----------------------------------------------------------------------------
# summary.rdtp()
# -----------------------------------------------------------------------------

#' Summary method for rdtp objects
#'
#' Produces a summary object with distributional statistics for cutoffs,
#' R-squared, and beta across units.
#'
#' @param object An object of class \code{"rdtp"}.
#' @param ... Additional arguments (ignored).
#'
#' @return An object of class \code{"summary.rdtp"}.
#' @export
summary.rdtp <- function(object, ...) {

    found_mask <- !is.na(object$results$cutoff)
    found_results <- object$results[found_mask, , drop = FALSE]

    # Compute distributional summaries for key columns
    if (nrow(found_results) > 0) {
        cutoff_dist <- c(
            summary(found_results$cutoff),
            sd = stats::sd(found_results$cutoff)
        )
        r2_dist <- c(
            summary(found_results$r2),
            sd = stats::sd(found_results$r2)
        )
        beta_dist <- c(
            summary(found_results$beta),
            sd = stats::sd(found_results$beta)
        )
        tstat_dist <- c(
            summary(found_results$tstat),
            sd = stats::sd(found_results$tstat)
        )
    } else {
        cutoff_dist <- r2_dist <- beta_dist <- tstat_dist <- NULL
    }

    structure(
        list(
            n_units       = object$n_units,
            n_found       = object$n_found,
            n_skipped     = object$n_skipped,
            mean_r2       = object$mean_r2,
            mean_beta     = object$mean_beta,
            mean_tstat    = object$mean_tstat,
            bandwidth     = object$bandwidth,
            vce           = object$vce,
            cutoff_dist   = cutoff_dist,
            r2_dist       = r2_dist,
            beta_dist     = beta_dist,
            tstat_dist    = tstat_dist,
            found_results = found_results,
            call          = object$call
        ),
        class = "summary.rdtp"
    )
}


# -----------------------------------------------------------------------------
# print.summary.rdtp()
# -----------------------------------------------------------------------------

#' Print method for summary.rdtp objects
#'
#' @param x An object of class \code{"summary.rdtp"}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#' @export
print.summary.rdtp <- function(x, ...) {

    cat("\nRD Tipping-Point Search -- Summary\n")
    cat(strrep("=", 60), "\n\n")
    cat(sprintf("Call: %s\n\n", deparse(x$call, width.cutoff = 500L)))
    cat(sprintf("Units searched: %d  |  Found: %d  |  Skipped: %d\n",
                x$n_units, x$n_found, x$n_skipped))
    cat(sprintf("Bandwidth: %g  |  VCE: %s\n\n",
                x$bandwidth, x$vce))

    if (x$n_found > 0) {
        cat("Distribution of results (units with valid cutoffs):\n")
        cat(strrep("-", 60), "\n")

        # Format and display each distribution
        .print_dist <- function(label, dist) {
            cat(sprintf("\n  %s:\n", label))
            cat(sprintf("    Min:    %10.4f\n", dist["Min."]))
            cat(sprintf("    1st Qu: %10.4f\n", dist["1st Qu."]))
            cat(sprintf("    Median: %10.4f\n", dist["Median"]))
            cat(sprintf("    Mean:   %10.4f\n", dist["Mean"]))
            cat(sprintf("    3rd Qu: %10.4f\n", dist["3rd Qu."]))
            cat(sprintf("    Max:    %10.4f\n", dist["Max."]))
            cat(sprintf("    SD:     %10.4f\n", dist["sd"]))
        }

        .print_dist("Cutoff", x$cutoff_dist)
        .print_dist("R-squared", x$r2_dist)
        .print_dist("Beta (discontinuity)", x$beta_dist)
        .print_dist("t-statistic", x$tstat_dist)
        cat("\n")
    } else {
        cat("No units had valid cutoffs.\n\n")
    }

    invisible(x)
}


# -----------------------------------------------------------------------------
# plot.rdtp()
# -----------------------------------------------------------------------------

#' Plot method for rdtp objects
#'
#' Produces base R histograms of detected cutoffs, R-squared values, or
#' discontinuity estimates across units.
#'
#' @param x An object of class \code{"rdtp"}.
#' @param type Character string. What to plot: \code{"cutoffs"} (default),
#'   \code{"r2"}, or \code{"beta"}.
#' @param ... Additional arguments passed to \code{\link[graphics]{hist}}.
#'
#' @return Invisibly returns the histogram object (from \code{\link[graphics]{hist}}).
#' @export
plot.rdtp <- function(x, type = "cutoffs", ...) {

    valid_types <- c("cutoffs", "r2", "beta")
    if (!(type %in% valid_types)) {
        stop(sprintf("'type' must be one of: %s",
                     paste(shQuote(valid_types), collapse = ", ")),
             call. = FALSE)
    }

    found_mask <- !is.na(x$results$cutoff)
    if (sum(found_mask) == 0) {
        stop("no units with valid cutoffs to plot", call. = FALSE)
    }

    found <- x$results[found_mask, , drop = FALSE]

    h <- switch(type,
        cutoffs = graphics::hist(found$cutoff,
                                 main = "Distribution of Detected Cutoffs",
                                 xlab = "Cutoff",
                                 col  = "steelblue",
                                 border = "white",
                                 ...),
        r2      = graphics::hist(found$r2,
                                 main = expression("Distribution of " * R^2),
                                 xlab = expression(R^2),
                                 col  = "steelblue",
                                 border = "white",
                                 ...),
        beta    = graphics::hist(found$beta,
                                 main = "Distribution of Discontinuity Estimates",
                                 xlab = expression(hat(beta)[1] ~ "(discontinuity)"),
                                 col  = "steelblue",
                                 border = "white",
                                 ...)
    )

    invisible(h)
}
