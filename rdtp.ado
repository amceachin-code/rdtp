*! rdtp v1.0.0  23feb2026
*! RD Tipping-Point Search for Discrete Forcing Variables
*! Andrew McEachin
*! Based on McEachin, Domina, & Penner (2020, JPAM)
*! doi:10.1002/pam.22202

// ═══════════════════════════════════════════════════════════════════════════
// rdtp — Main command
// ═══════════════════════════════════════════════════════════════════════════
//
// Searches over candidate cutoff values on a discrete running variable to
// find, for each unit, the cutoff that maximizes R² from a linear-spline
// regression discontinuity (RD) regression.
//
// Model at each candidate cutoff c:
//   depvar = β₀ + β₁·cut + β₂·force + β₃·(cut × force) + ε
//   where cut = 1(forcing ≥ c), force = forcing − c
//
// The cutoff c* that maximizes R² is selected for each unit.
// ═══════════════════════════════════════════════════════════════════════════

program define rdtp, rclass
    version 16.0

    // ─────────────────────────────────────────────────────────────────────
    // 1. SYNTAX PARSING
    // ─────────────────────────────────────────────────────────────────────

    syntax varlist(min=2 max=2 numeric) [if] [in], ///
        BY(varname)                                  ///
        [BANDwidth(real 75)                          ///
         SEARCHRange(numlist min=2 max=2)             ///
         MINobs(integer 10)                          ///
         vce(string)                                 ///
         SAVing(string)                              ///
         REPLACE                                     ///
         NOIsily                                     ///
         Level(cilevel)]

    // Separate the dependent variable from the forcing variable
    gettoken depvar forcing : varlist

    // ─────────────────────────────────────────────────────────────────────
    // 2. INPUT VALIDATION
    // ─────────────────────────────────────────────────────────────────────

    // by() cannot duplicate the outcome or running variable
    if "`by'" == "`depvar'" {
        di as error "by() variable cannot be the same as the dependent variable"
        exit 198
    }
    if "`by'" == "`forcing'" {
        di as error "by() variable cannot be the same as the forcing variable"
        exit 198
    }

    // Bandwidth must be strictly positive
    if `bandwidth' <= 0 {
        di as error "bandwidth() must be a positive number"
        exit 198
    }

    // Minimum observations per side must be at least 1
    if `minobs' < 1 {
        di as error "minobs() must be at least 1"
        exit 198
    }

    // If searchrange() specified, the lower bound must be below the upper
    if "`searchrange'" != "" {
        local sr_min : word 1 of `searchrange'
        local sr_max : word 2 of `searchrange'
        if `sr_min' >= `sr_max' {
            di as error ///
                "searchrange() minimum (`sr_min') must be less than maximum (`sr_max')"
            exit 198
        }
    }

    // If saving() specified without replace, the file must not already exist
    if `"`saving'"' != "" {
        local savefile `"`saving'"'
        // Append .dta extension if not already present
        if substr(`"`savefile'"', -4, 4) != ".dta" {
            local savefile `"`savefile'.dta"'
        }
        if "`replace'" == "" {
            confirm new file `"`savefile'"'
        }
    }

    // Build the estimation-sample marker (drops obs missing depvar,
    // forcing, or by)
    marksample touse
    markout `touse' `by'

    // Ensure non-empty estimation sample
    quietly count if `touse'
    if r(N) == 0 {
        di as error "no observations in estimation sample"
        exit 2000
    }

    // ─────────────────────────────────────────────────────────────────────
    // 3. SETUP
    // ─────────────────────────────────────────────────────────────────────

    // Determine whether the by-variable is string or numeric
    local by_is_string 0
    capture confirm string variable `by'
    if _rc == 0 {
        local by_is_string 1
    }

    // Enumerate unique units
    quietly levelsof `by' if `touse', local(units)
    local n_units : word count `units'
    if `n_units' == 0 {
        di as error "no units found in by() variable"
        exit 2000
    }

    // Build the vce() clause for regress (empty string if default OLS)
    local vce_clause ""
    if `"`vce'"' != "" {
        local vce_clause `"vce(`vce')"'
    }

    // ─────────────────────────────────────────────────────────────────────
    // Postfile setup: one row per unit
    // ─────────────────────────────────────────────────────────────────────

    tempname pf
    tempfile resultsfile

    postfile `pf'                              ///
        str244  unit                           ///
        double  cutoff                         ///
        double  r2                             ///
        double  beta                           ///
        double  se                             ///
        double  tstat                          ///
        long    n_left                         ///
        long    n_right                        ///
        long    n_total                        ///
        double  pred_left                      ///
        double  pred_right                     ///
        using `"`resultsfile'"'

    // ─────────────────────────────────────────────────────────────────────
    // VCE pre-validation
    //   If the user specified vce(), run a test regression on the first
    //   unit's data to catch invalid specifications (e.g., a misspelled
    //   cluster variable) before entering the main loop.  Without this
    //   check, every unit silently fails via `capture regress` and the
    //   user gets n_found = 0 with no useful error message.
    // ─────────────────────────────────────────────────────────────────────

    if `"`vce'"' != "" {
        // Grab the first unit identifier
        local first_unit : word 1 of `units'

        preserve

        // Subset to the first unit's estimation-sample data
        if `by_is_string' {
            quietly keep if `by' == `"`first_unit'"' & `touse'
        }
        else {
            quietly keep if `by' == `first_unit' & `touse'
        }

        // Get the first candidate cutoff for this unit
        if "`searchrange'" != "" {
            quietly levelsof `forcing' ///
                if `forcing' >= `sr_min' & `forcing' <= `sr_max', ///
                local(test_candidates)
        }
        else {
            quietly levelsof `forcing', local(test_candidates)
        }
        local test_c : word 1 of `test_candidates'

        if "`test_c'" != "" {
            // Construct RD variables at this candidate cutoff
            tempvar _tcut _tfrc _tint
            quietly gen byte   `_tcut' = (`forcing' >= `test_c') ///
                if !missing(`forcing')
            quietly gen double `_tfrc' = `forcing' - `test_c'    ///
                if !missing(`forcing')
            quietly gen double `_tint' = `_tcut' * `_tfrc'

            capture regress `depvar' `_tcut' `_tfrc' `_tint' ///
                if abs(`_tfrc') <= `bandwidth', `vce_clause'

            if _rc != 0 {
                local vce_rc = _rc
                restore
                postclose `pf'
                di as error ///
                    "vce(`vce') failed on test regression (error code `vce_rc')."
                di as error ///
                    "Check that all variable names in vce() exist and are spelled correctly."
                exit `vce_rc'
            }
        }
        else {
            // First unit had no candidates in search range — fall back
            // to a simple regress on whatever data we have, just to
            // validate that Stata accepts the vce() specification.
            capture regress `depvar' `forcing', `vce_clause'
            if _rc != 0 {
                local vce_rc = _rc
                restore
                postclose `pf'
                di as error ///
                    "vce(`vce') failed on test regression (error code `vce_rc')."
                di as error ///
                    "Check that all variable names in vce() exist and are spelled correctly."
                exit `vce_rc'
            }
        }

        restore
    }

    // ─────────────────────────────────────────────────────────────────────
    // Display header
    // ─────────────────────────────────────────────────────────────────────

    di _n as text "RD Tipping-Point Search"
    di    as text "{hline 60}"
    di    as text "Dependent variable:  " as result "`depvar'"
    di    as text "Forcing variable:    " as result "`forcing'"
    di    as text "Unit variable:       " as result "`by'"
    di    as text "Bandwidth:           " as result `bandwidth'
    di    as text "Min obs per side:    " as result `minobs'
    if `"`vce'"' != "" {
        di as text "VCE:                 " as result `"`vce'"'
    }
    else {
        di as text "VCE:                 " as result "OLS (default)"
    }
    if "`searchrange'" != "" {
        di as text "Search range:        " as result "[`sr_min', `sr_max']"
    }
    else {
        di as text "Search range:        " as result "(all unique values)"
    }
    di    as text "Units to search:     " as result `n_units'
    di    as text "{hline 60}"

    // ─────────────────────────────────────────────────────────────────────
    // 4. MAIN LOOP — iterate over units, search over cutoffs
    // ─────────────────────────────────────────────────────────────────────

    // Invariant: n_found + n_skipped == n_units (every unit posts exactly one row)
    local n_found   0
    local n_skipped 0
    local uctr      0

    preserve

    foreach u of local units {
        local ++uctr

        // ── Reset to the full dataset ─────────────────────────────────
        restore, preserve

        // ── Subset to this unit within the estimation sample ──────────
        if `by_is_string' {
            quietly keep if `by' == `"`u'"' & `touse'
        }
        else {
            quietly keep if `by' == `u' & `touse'
        }

        quietly count
        if r(N) == 0 {
            // No observations for this unit — post a missing row
            post `pf' (`"`u'"') (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)
            local ++n_skipped
            if "`noisily'" != "" {
                di as text ///
                    "  [`uctr'/`n_units'] `u': no observations — skipped"
            }
            continue
        }

        // ── Identify candidate cutoff values ──────────────────────────
        if "`searchrange'" != "" {
            quietly levelsof `forcing' ///
                if `forcing' >= `sr_min' & `forcing' <= `sr_max', ///
                local(candidates)
        }
        else {
            quietly levelsof `forcing', local(candidates)
        }

        local n_cand : word count `candidates'
        if `n_cand' == 0 {
            post `pf' (`"`u'"') (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)
            local ++n_skipped
            if "`noisily'" != "" {
                di as text ///
                    "  [`uctr'/`n_units'] `u': no candidates in range — skipped"
            }
            continue
        }

        // ── Search over candidate cutoffs ─────────────────────────────
        // Track the best (highest-R²) cutoff for this unit
        local best_r2     = -1
        local best_cutoff = .
        local best_beta   = .
        local best_se     = .
        local best_tstat  = .
        local best_nleft  = .
        local best_nright = .
        local best_ntotal = .
        local best_predL  = .
        local best_predR  = .
        local found_any   0

        // Allocate tempvar names once per unit (reused across candidates)
        tempvar _cut _frc _int

        foreach c of local candidates {

            // ·· Construct the RD variables ··
            quietly gen byte   `_cut' = (`forcing' >= `c') ///
                if !missing(`forcing')
            quietly gen double `_frc' = `forcing' - `c'    ///
                if !missing(`forcing')
            quietly gen double `_int' = `_cut' * `_frc'

            // ·· Count observations on each side within bandwidth ··
            quietly count if `_cut' == 0 & abs(`_frc') <= `bandwidth'
            local nl = r(N)
            quietly count if `_cut' == 1 & abs(`_frc') <= `bandwidth'
            local nr = r(N)

            // Skip if either side has too few observations
            if `nl' < `minobs' | `nr' < `minobs' {
                drop `_cut' `_frc' `_int'
                continue
            }

            // ·· Run the linear-spline RD regression ··
            //    depvar = β₀ + β₁·cut + β₂·force + β₃·(cut×force) + ε
            capture regress `depvar' `_cut' `_frc' `_int' ///
                if abs(`_frc') <= `bandwidth', `vce_clause' level(`level')

            // Check if regression succeeded
            if _rc != 0 {
                drop `_cut' `_frc' `_int'
                continue
            }

            // ·· Skip degenerate regressions ··
            //    - Missing F → model not identified (e.g., too few clusters)
            //    - β₁ = 0   → Stata sets _b[varname] = 0 when a variable
            //      is dropped due to collinearity, which is the only
            //      realistic way this exact-zero condition triggers.
            //      A numerically estimated coefficient is essentially
            //      never exactly 0.0 in floating point.
            if missing(e(F)) | _b[`_cut'] == 0 {
                drop `_cut' `_frc' `_int'
                continue
            }

            // ·· Update best if this R² is the new maximum ··
            //    Tie-breaking: strict > means the first candidate
            //    (lowest forcing value, since levelsof returns sorted)
            //    wins when two cutoffs produce identical R².
            if e(r2) > `best_r2' {
                local best_r2     = e(r2)
                local best_cutoff = `c'
                local best_beta   = _b[`_cut']
                local best_se     = _se[`_cut']
                local best_tstat  = _b[`_cut'] / _se[`_cut']
                local best_nleft  = `nl'
                local best_nright = `nr'
                local best_ntotal = e(N)
                // Predicted outcome at cutoff from the LEFT (cut=0, force=0):
                //   ŷ_left = β₀ = _b[_cons]
                local best_predL  = _b[_cons]
                // Predicted outcome at cutoff from the RIGHT (cut=1, force=0):
                //   ŷ_right = β₀ + β₁ = _b[_cons] + _b[cut]
                local best_predR  = _b[_cons] + _b[`_cut']
                local found_any   1
            }

            drop `_cut' `_frc' `_int'
        }

        // ── Post this unit's results ──────────────────────────────────
        if `found_any' {
            post `pf' (`"`u'"') (`best_cutoff') (`best_r2') ///
                (`best_beta') (`best_se') (`best_tstat')    ///
                (`best_nleft') (`best_nright') (`best_ntotal') ///
                (`best_predL') (`best_predR')
            local ++n_found
            if "`noisily'" != "" {
                di as text "  [`uctr'/`n_units'] `u': " ///
                    as result "cutoff = " %9.0g `best_cutoff' ///
                    as text  ", R{c 178} = " as result %6.4f `best_r2' ///
                    as text  ", {c beta} = " as result %9.3f `best_beta'
            }
        }
        else {
            post `pf' (`"`u'"') (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)
            local ++n_skipped
            if "`noisily'" != "" {
                di as text ///
                    "  [`uctr'/`n_units'] `u': no valid cutoff — skipped"
            }
        }

        // Progress display (every 50 units, unless noisily is on)
        if "`noisily'" == "" & mod(`uctr', 50) == 0 & `uctr' < `n_units' {
            di as text "  ... processed `uctr' of `n_units' units"
        }
    }

    // Done with the unit loop
    restore

    // ─────────────────────────────────────────────────────────────────────
    // 5. POST-PROCESS RESULTS
    // ─────────────────────────────────────────────────────────────────────

    postclose `pf'

    preserve
    quietly use `"`resultsfile'"', clear

    // ── Summary statistics (among units with valid cutoffs) ────────────
    if `n_found' > 0 {
        quietly summarize r2 if !missing(cutoff)
        local mean_r2    = r(mean)
        quietly summarize beta if !missing(cutoff)
        local mean_beta  = r(mean)
        quietly summarize tstat if !missing(cutoff)
        local mean_tstat = r(mean)
    }
    else {
        local mean_r2    = .
        local mean_beta  = .
        local mean_tstat = .
    }

    // ── Build the results matrix (n_units × 10) ──────────────────────
    local nrows = _N
    tempname rmat
    mkmat cutoff r2 beta se tstat n_left n_right n_total ///
        pred_left pred_right, matrix(`rmat')

    // Label rows with (cleaned) unit identifiers
    //   Stata matrix row names: max 32 chars, no spaces
    local rownames ""
    forvalues i = 1/`nrows' {
        local lbl = unit[`i']
        local lbl = subinstr("`lbl'", " ", "_", .)
        local lbl = subinstr("`lbl'", ".", "_", .)
        local lbl = subinstr("`lbl'", ":", "_", .)
        local lbl = subinstr("`lbl'", ",", "_", .)
        local lbl = subinstr("`lbl'", `"""', "_", .)
        local lbl = subinstr("`lbl'", "\", "_", .)
        if strlen("`lbl'") > 32 {
            local lbl = substr("`lbl'", 1, 32)
        }
        if "`lbl'" == "" {
            local lbl "unit_`i'"
        }
        local rownames "`rownames' `lbl'"
    }
    capture matrix rownames `rmat' = `rownames'
    if _rc {
        // Fallback: sequential labels if cleaning produced duplicates
        local rownames ""
        forvalues i = 1/`nrows' {
            local rownames "`rownames' unit_`i'"
        }
        matrix rownames `rmat' = `rownames'
    }

    // ── Display per-unit table (if ≤ 50 units and at least one found) ─
    if `n_units' <= 50 & `n_found' > 0 {
        di _n as text "{hline 77}"
        di as text ///
            %18s "Unit" " {c |}" ///
            %10s "Cutoff"        ///
            %8s  "R{c 178}"     ///
            %11s "Beta"          ///
            %10s "SE"            ///
            %9s  "t"             ///
            %8s  "N"
        di as text "{hline 77}"

        forvalues i = 1/`nrows' {
            local u_lbl = unit[`i']
            // Truncate long names for display
            if strlen(`"`u_lbl'"') > 17 {
                local u_lbl = substr(`"`u_lbl'"', 1, 14) + "..."
            }

            if !missing(cutoff[`i']) {
                di as text %18s `"`u_lbl'"' " {c |}" ///
                    as result                         ///
                    %10.0g cutoff[`i']                ///
                    %8.4f  r2[`i']                    ///
                    %11.3f beta[`i']                  ///
                    %10.3f se[`i']                    ///
                    %9.3f  tstat[`i']                 ///
                    %8.0f  n_total[`i']
            }
            else {
                di as text %18s `"`u_lbl'"' " {c |}" ///
                    as result                         ///
                    %10s "."                           ///
                    %8s  "."                           ///
                    %11s "."                           ///
                    %10s "."                           ///
                    %9s  "."                           ///
                    %8s  "."
            }
        }
        di as text "{hline 77}"
    }

    // ── Save results dataset if requested ─────────────────────────────
    // (savefile local was set during input validation, reuse it here)
    if `"`saving'"' != "" {
        quietly save `"`savefile'"', `replace'
    }

    restore

    // ─────────────────────────────────────────────────────────────────────
    // 6. STORE R-CLASS RESULTS
    // ─────────────────────────────────────────────────────────────────────

    // Scalars
    return scalar n_units    = `n_units'
    return scalar n_found    = `n_found'
    return scalar n_skipped  = `n_skipped'
    return scalar bandwidth  = `bandwidth'
    return scalar minobs     = `minobs'
    return scalar mean_r2    = `mean_r2'
    return scalar mean_beta  = `mean_beta'
    return scalar mean_tstat = `mean_tstat'
    return scalar level      = `level'

    // Macros
    return local cmd         "rdtp"
    return local cmdline     `"rdtp `0'"'
    return local depvar      "`depvar'"
    return local forcing     "`forcing'"
    return local byvar       "`by'"
    if `"`vce'"' != "" {
        return local vce     `"`vce'"'
    }
    if "`searchrange'" != "" {
        return local searchrange "`sr_min' `sr_max'"
    }
    if `"`saving'"' != "" {
        return local saving  `"`savefile'"'
    }

    // Matrix
    return matrix results = `rmat'

    // ─────────────────────────────────────────────────────────────────────
    // 7. DISPLAY SUMMARY
    // ─────────────────────────────────────────────────────────────────────

    di _n as text "{hline 60}"
    di    as text "Results"
    di    as text "{hline 60}"
    di    as text "Units searched:       " as result `n_units'
    di    as text "Cutoffs found:        " as result `n_found'
    di    as text "Units skipped/failed: " as result `n_skipped'

    if `n_found' > 0 {
        di _n as text "Summary (units with valid cutoffs):"
        di    as text "  Mean R{c 178}:            " as result %9.4f `mean_r2'
        di    as text "  Mean discontinuity: " as result %9.4f `mean_beta'
        di    as text "  Mean t-statistic:   " as result %9.4f `mean_tstat'
    }

    if `n_units' > 50 & `n_found' > 0 {
        di _n as text "(Per-unit table suppressed: `n_units' units > 50." ///
            " Use {bf:saving()} to inspect individual results.)"
    }

    if `"`saving'"' != "" {
        di _n as text "Results saved to: " as result `"`savefile'"'
    }

    di ""

end
