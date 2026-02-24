# ===============================================================================
# utils.R -- Internal validation and VCE helper functions for rdtp
# ===============================================================================


# -----------------------------------------------------------------------------
# .rdtp_validate_inputs()
#
# Validates all user-supplied arguments before entering the main loop.
# Errors early with informative messages. Called at the top of rdtp().
# -----------------------------------------------------------------------------

.rdtp_validate_inputs <- function(data, depvar, forcing, by,
                                  bandwidth, searchrange, minobs,
                                  vce, cluster, subset, saving,
                                  replace, level) {

    # -- data must be a data.frame ------------------------------------------
    if (!is.data.frame(data)) {
        stop("'data' must be a data.frame", call. = FALSE)
    }
    if (nrow(data) == 0) {
        stop("'data' has zero rows", call. = FALSE)
    }

    # -- depvar, forcing, by must be character strings naming columns --------
    if (!is.character(depvar) || length(depvar) != 1) {
        stop("'depvar' must be a single character string naming a column in data",
             call. = FALSE)
    }
    if (!is.character(forcing) || length(forcing) != 1) {
        stop("'forcing' must be a single character string naming a column in data",
             call. = FALSE)
    }
    if (!is.character(by) || length(by) != 1) {
        stop("'by' must be a single character string naming a column in data",
             call. = FALSE)
    }

    # -- columns must exist in data -----------------------------------------
    if (!(depvar %in% names(data))) {
        stop(sprintf("column '%s' (depvar) not found in data", depvar),
             call. = FALSE)
    }
    if (!(forcing %in% names(data))) {
        stop(sprintf("column '%s' (forcing) not found in data", forcing),
             call. = FALSE)
    }
    if (!(by %in% names(data))) {
        stop(sprintf("column '%s' (by) not found in data", by),
             call. = FALSE)
    }

    # -- depvar and forcing must be numeric ---------------------------------
    if (!is.numeric(data[[depvar]])) {
        stop(sprintf("'%s' (depvar) must be numeric", depvar), call. = FALSE)
    }
    if (!is.numeric(data[[forcing]])) {
        stop(sprintf("'%s' (forcing) must be numeric", forcing), call. = FALSE)
    }

    # -- by() cannot duplicate depvar or forcing ----------------------------
    if (by == depvar) {
        stop("'by' variable cannot be the same as the dependent variable",
             call. = FALSE)
    }
    if (by == forcing) {
        stop("'by' variable cannot be the same as the forcing variable",
             call. = FALSE)
    }

    # -- bandwidth must be strictly positive --------------------------------
    if (!is.numeric(bandwidth) || length(bandwidth) != 1 ||
        is.na(bandwidth) || bandwidth <= 0) {
        stop("'bandwidth' must be a single positive number", call. = FALSE)
    }

    # -- minobs must be at least 1 ------------------------------------------
    if (!is.numeric(minobs) || length(minobs) != 1 ||
        is.na(minobs) || minobs < 1) {
        stop("'minobs' must be at least 1", call. = FALSE)
    }

    # -- searchrange validation ---------------------------------------------
    if (!is.null(searchrange)) {
        if (!is.numeric(searchrange) || length(searchrange) != 2) {
            stop("'searchrange' must be a numeric vector of length 2 (min, max)",
                 call. = FALSE)
        }
        if (any(is.na(searchrange))) {
            stop("'searchrange' cannot contain NA values", call. = FALSE)
        }
        if (searchrange[1] >= searchrange[2]) {
            stop(sprintf(
                "searchrange minimum (%.4g) must be less than maximum (%.4g)",
                searchrange[1], searchrange[2]),
                call. = FALSE)
        }
    }

    # -- vce validation -----------------------------------------------------
    valid_vce <- c("ols", "robust", "cluster")
    if (!is.character(vce) || length(vce) != 1 || !(vce %in% valid_vce)) {
        stop(sprintf("'vce' must be one of: %s",
                     paste(shQuote(valid_vce), collapse = ", ")),
             call. = FALSE)
    }

    # -- cluster variable validation ----------------------------------------
    if (vce == "cluster") {
        if (is.null(cluster)) {
            stop("'cluster' must be specified when vce = \"cluster\"",
                 call. = FALSE)
        }
        if (!is.character(cluster) || length(cluster) != 1) {
            stop("'cluster' must be a single character string naming a column in data",
                 call. = FALSE)
        }
        if (!(cluster %in% names(data))) {
            stop(sprintf("cluster variable '%s' not found in data", cluster),
                 call. = FALSE)
        }
    }
    if (vce != "cluster" && !is.null(cluster)) {
        warning("'cluster' is ignored when vce != \"cluster\"", call. = FALSE)
    }

    # -- subset validation --------------------------------------------------
    if (!is.null(subset)) {
        if (!is.logical(subset) && !is.numeric(subset)) {
            stop("'subset' must be a logical or numeric vector", call. = FALSE)
        }
    }

    # -- saving validation --------------------------------------------------
    if (!is.null(saving)) {
        if (!is.character(saving) || length(saving) != 1) {
            stop("'saving' must be a single character string (file path)",
                 call. = FALSE)
        }
        # Check file extension -- support .rds and .csv
        ext <- tolower(tools::file_ext(saving))
        if (!(ext %in% c("rds", "csv", ""))) {
            stop("'saving' file must have .rds or .csv extension (or no extension for .rds default)",
                 call. = FALSE)
        }
        # If no extension, default to .rds
        if (ext == "") {
            saving <- paste0(saving, ".rds")
        }
        # Check if file exists (unless replace = TRUE)
        if (!replace && file.exists(saving)) {
            stop(sprintf("file '%s' already exists. Use replace = TRUE to overwrite.",
                         saving),
                 call. = FALSE)
        }
    }

    # -- level validation ---------------------------------------------------
    if (!is.numeric(level) || length(level) != 1 ||
        is.na(level) || level <= 0 || level >= 1) {
        stop("'level' must be a single number in (0, 1)", call. = FALSE)
    }

    # Return the (possibly modified) saving path
    invisible(saving)
}


# -----------------------------------------------------------------------------
# .rdtp_get_se()
#
# Extracts the standard error for the "cut" coefficient from a fitted lm
# object, using the specified VCE method.
#
# Arguments:
#   fit      -- a fitted lm() object
#   vce      -- "ols", "robust", or "cluster"
#   cluster  -- cluster vector (same length as model data), or NULL
#
# Returns: scalar SE for the "cut" coefficient (second coefficient).
# -----------------------------------------------------------------------------

.rdtp_get_se <- function(fit, vce = "ols", cluster_vec = NULL) {

    if (vce == "ols") {
        # Standard OLS standard errors
        return(summary(fit)$coefficients["cut", "Std. Error"])
    }

    # sandwich is needed for both robust and cluster
    if (!requireNamespace("sandwich", quietly = TRUE)) {
        stop(
            "Package 'sandwich' is required for vce = \"", vce, "\". ",
            "Install it with: install.packages(\"sandwich\")",
            call. = FALSE
        )
    }

    if (vce == "robust") {
        V <- sandwich::vcovHC(fit, type = "HC1")
    } else if (vce == "cluster") {
        V <- sandwich::vcovCL(fit, cluster = cluster_vec)
    } else {
        stop(sprintf("Unknown vce type: '%s'", vce), call. = FALSE)
    }

    sqrt(V["cut", "cut"])
}


# -----------------------------------------------------------------------------
# .rdtp_vce_prevalidation()
#
# If vce != "ols", runs a quick test regression on one unit's data to
# verify that the VCE specification works. Catches errors early before
# entering the main loop.
# -----------------------------------------------------------------------------

.rdtp_vce_prevalidation <- function(vce, cluster_var, data_one_unit) {

    if (vce == "ols") return(invisible(NULL))

    # Need at least a few observations for a test fit
    if (nrow(data_one_unit) < 5) return(invisible(NULL))

    # Simple test regression
    test_y <- data_one_unit[[1]]  # first column
    test_x <- seq_len(nrow(data_one_unit))

    test_fit <- tryCatch(
        lm(test_y ~ test_x),
        error = function(e) NULL
    )

    if (is.null(test_fit)) return(invisible(NULL))

    # sandwich is needed for both robust and cluster
    if (!requireNamespace("sandwich", quietly = TRUE)) {
        stop(
            "Package 'sandwich' is required for vce = \"", vce, "\". ",
            "Install it with: install.packages(\"sandwich\")",
            call. = FALSE
        )
    }

    # Try the VCE computation
    if (vce == "robust") {
        tryCatch(
            sandwich::vcovHC(test_fit, type = "HC1"),
            error = function(e) {
                stop(sprintf(
                    "vce = \"robust\" failed on test regression: %s", e$message),
                    call. = FALSE)
            }
        )
    } else if (vce == "cluster") {
        clust_vec <- data_one_unit[[cluster_var]]
        tryCatch(
            sandwich::vcovCL(test_fit, cluster = clust_vec),
            error = function(e) {
                stop(sprintf(
                    "vce = \"cluster\" failed on test regression: %s", e$message),
                    call. = FALSE)
            }
        )
    }

    invisible(NULL)
}
