*! test_rdtp.do — Automated test suite for rdtp
*! Run from Stata:  do tests/test_rdtp.do
*!
*! Each test block generates synthetic data with known properties, runs
*! rdtp, and asserts that return values match expectations.  If any
*! assertion fails, Stata stops with an error message identifying the
*! failing test.
*!
*! Exit code 0 = all tests passed.

clear all
set more off

// ─────────────────────────────────────────────────────────────────────────
// Ensure rdtp is on the adopath (handles running from the tests/ dir
// or from the project root)
// ─────────────────────────────────────────────────────────────────────────
local test_dir = c(pwd)
capture which rdtp
if _rc {
    // Try adding the parent directory (project root)
    adopath + "../"
    capture which rdtp
    if _rc {
        di as error "rdtp.ado not found. Run this from the project root or tests/ directory."
        exit 111
    }
}

local n_tests  = 0
local n_passed = 0


// ═════════════════════════════════════════════════════════════════════════
// TEST 1: Basic functionality — 5 units with known cutoffs
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 1: Basic functionality (5 units, known cutoffs) ═══"

clear
set seed 12345
set obs 500
gen school = ceil(_n / 100)
gen test_score = 200 + int(runiform() * 301)

// True cutoffs: 300, 350, 400, 325, 375
gen true_cutoff = cond(school==1, 300, cond(school==2, 350, ///
    cond(school==3, 400, cond(school==4, 325, 375))))
gen above = (test_score >= true_cutoff)

// Strong signal (large discontinuity, low noise) so cutoffs are recoverable
gen outcome = 50 + 0.1 * test_score + 20 * above + rnormal(0, 3)

rdtp outcome test_score, by(school) searchrange(250 450)

// Assertions
local ++n_tests
assert r(n_units)   == 5
assert r(n_found)   == 5
assert r(n_skipped) == 0
assert r(bandwidth) == 75
assert r(minobs)    == 10
assert r(mean_r2)    > 0     // R² should be positive
assert r(mean_r2)    < 1     // and less than 1
assert r(mean_beta)  > 0     // discontinuity should be positive (we set it to +20)
local ++n_passed
di as result "  PASSED"


// ═════════════════════════════════════════════════════════════════════════
// TEST 2: Cutoff accuracy — verify detected cutoffs are near true values
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 2: Cutoff accuracy (detected ≈ true) ═══"

// Use results matrix from Test 1 (still in r(results))
matrix R = r(results)

// Each detected cutoff should be within ±25 of the true value
// (with strong signal and 100 obs per unit, this is conservative)
local ++n_tests
assert abs(R[1,1] - 300) <= 25   // school 1: true = 300
assert abs(R[2,1] - 350) <= 25   // school 2: true = 350
assert abs(R[3,1] - 400) <= 25   // school 3: true = 400
assert abs(R[4,1] - 325) <= 25   // school 4: true = 325
assert abs(R[5,1] - 375) <= 25   // school 5: true = 375
local ++n_passed
di as result "  PASSED"


// ═════════════════════════════════════════════════════════════════════════
// TEST 3: String by-variable
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 3: String by-variable ═══"

clear
set seed 99999
set obs 300
gen str10 school_name = cond(_n <= 100, "Alpha", ///
    cond(_n <= 200, "Beta", "Gamma"))
gen score = 1 + int(runiform() * 100)
gen true_c = cond(school_name == "Alpha", 40, ///
    cond(school_name == "Beta", 60, 50))
gen above = (score >= true_c)
gen y = 10 + 5 * above + rnormal(0, 2)

rdtp y score, by(school_name) bandwidth(50)

local ++n_tests
assert r(n_units)  == 3
assert r(n_found)  == 3
assert "`r(byvar)'" == "school_name"
local ++n_passed
di as result "  PASSED"


// ═════════════════════════════════════════════════════════════════════════
// TEST 4: No valid cutoffs (all candidates have too few obs)
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 4: No valid cutoffs (minobs too high) ═══"

clear
set seed 11111
set obs 20
gen unit = 1
gen x = _n
gen y = rnormal()

// With only 20 obs and minobs=50, no candidate can have enough on both sides
rdtp y x, by(unit) minobs(50)

local ++n_tests
assert r(n_units)   == 1
assert r(n_found)   == 0
assert r(n_skipped) == 1
assert r(mean_r2)   == .
local ++n_passed
di as result "  PASSED"


// ═════════════════════════════════════════════════════════════════════════
// TEST 5: Single unit
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 5: Single unit ═══"

clear
set seed 22222
set obs 200
gen unit = 1
gen x = int(runiform() * 100)
gen above = (x >= 50)
gen y = 3 + 8 * above + rnormal(0, 2)

rdtp y x, by(unit) bandwidth(50)

local ++n_tests
assert r(n_units)  == 1
assert r(n_found)  == 1
assert r(mean_beta) > 0
local ++n_passed
di as result "  PASSED"


// ═════════════════════════════════════════════════════════════════════════
// TEST 6: searchrange() restricts candidates
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 6: searchrange() restricts candidates ═══"

clear
set seed 33333
set obs 200
gen unit = 1
gen x = int(runiform() * 100)
gen above = (x >= 50)
gen y = 3 + 8 * above + rnormal(0, 2)

// Search only 60-90: the true cutoff at 50 is excluded, so the best
// fit will be worse and the cutoff will be in [60, 90]
rdtp y x, by(unit) searchrange(60 90) bandwidth(50)

local ++n_tests
assert r(n_found) == 1
matrix R = r(results)
assert R[1,1] >= 60 & R[1,1] <= 90
local ++n_passed
di as result "  PASSED"


// ═════════════════════════════════════════════════════════════════════════
// TEST 7: Non-integer searchrange (the fix we just made)
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 7: Non-integer searchrange ═══"

clear
set seed 44444
set obs 300
gen unit = 1
// GPA-like forcing variable with 0.1 increments
gen gpa = 1.0 + int(runiform() * 30) / 10
gen above = (gpa >= 2.5)
gen y = 10 + 4 * above + rnormal(0, 2)

// This would have failed with the old `integer` constraint
rdtp y gpa, by(unit) searchrange(1.5 3.5) bandwidth(2)

local ++n_tests
assert r(n_found) == 1
local ++n_passed
di as result "  PASSED"


// ═════════════════════════════════════════════════════════════════════════
// TEST 8: vce(robust) works
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 8: vce(robust) ═══"

clear
set seed 55555
set obs 200
gen unit = 1
gen x = int(runiform() * 100)
gen above = (x >= 50)
gen y = 3 + 8 * above + rnormal(0, 2)

rdtp y x, by(unit) vce(robust) bandwidth(50)

local ++n_tests
assert r(n_found) == 1
assert "`r(vce)'" == "robust"
local ++n_passed
di as result "  PASSED"


// ═════════════════════════════════════════════════════════════════════════
// TEST 9: vce(cluster varname) works
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 9: vce(cluster varname) ═══"

clear
set seed 66666
set obs 500
gen school = ceil(_n / 100)
gen classroom = ceil(_n / 25)
gen x = 200 + int(runiform() * 301)
gen true_c = cond(school==1, 300, cond(school==2, 350, ///
    cond(school==3, 400, cond(school==4, 325, 375))))
gen above = (x >= true_c)
gen y = 50 + 0.1 * x + 15 * above + rnormal(0, 5)

rdtp y x, by(school) vce(cluster classroom) searchrange(250 450)

local ++n_tests
assert r(n_found) == 5
assert "`r(vce)'" == "cluster classroom"
local ++n_passed
di as result "  PASSED"


// ═════════════════════════════════════════════════════════════════════════
// TEST 10: VCE pre-validation catches invalid cluster variable
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 10: VCE pre-validation catches bad vce() ═══"

clear
set seed 77777
set obs 200
gen unit = 1
gen x = int(runiform() * 100)
gen y = rnormal()

// This should error immediately (not silently skip)
capture rdtp y x, by(unit) vce(cluster nonexistent_var)

local ++n_tests
assert _rc != 0
local ++n_passed
di as result "  PASSED (rdtp correctly errored with _rc = " _rc ")"


// ═════════════════════════════════════════════════════════════════════════
// TEST 11: saving() produces a valid dataset
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 11: saving() produces valid dataset ═══"

clear
set seed 88888
set obs 200
gen school = ceil(_n / 100)
gen x = int(runiform() * 100)
gen true_c = cond(school == 1, 40, 60)
gen above = (x >= true_c)
gen y = 5 + 6 * above + rnormal(0, 2)

tempfile results_file
rdtp y x, by(school) saving(`results_file') bandwidth(50)

local ++n_tests
// Verify the saved file exists and has the right structure
preserve
use `results_file', clear
assert _N == 2                    // 2 units = 2 rows
confirm variable unit cutoff r2 beta se tstat n_left n_right n_total pred_left pred_right
restore
local ++n_passed
di as result "  PASSED"


// ═════════════════════════════════════════════════════════════════════════
// TEST 12: r(results) matrix dimensions and column names
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 12: r(results) matrix structure ═══"

clear
set seed 10101
set obs 300
gen unit = ceil(_n / 100)
gen x = int(runiform() * 100)
gen above = (x >= 50)
gen y = 3 + 5 * above + rnormal(0, 2)

rdtp y x, by(unit) bandwidth(50)

local ++n_tests
matrix R = r(results)
// Should have 3 rows (3 units) and 10 columns
assert rowsof(R) == 3
assert colsof(R) == 10
local ++n_passed
di as result "  PASSED"


// ═════════════════════════════════════════════════════════════════════════
// TEST 13: Empty estimation sample errors gracefully
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 13: Empty estimation sample ═══"

clear
set seed 12121
set obs 100
gen unit = 1
gen x = int(runiform() * 100)
gen y = rnormal()

// Use an if-condition that selects zero observations
capture rdtp y x if x > 999, by(unit)

local ++n_tests
assert _rc != 0
local ++n_passed
di as result "  PASSED (rdtp correctly errored with _rc = " _rc ")"


// ═════════════════════════════════════════════════════════════════════════
// TEST 14: noisily option runs without error
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 14: noisily option ═══"

clear
set seed 14141
set obs 200
gen unit = ceil(_n / 100)
gen x = int(runiform() * 100)
gen above = (x >= 50)
gen y = 3 + 5 * above + rnormal(0, 2)

capture rdtp y x, by(unit) bandwidth(50) noisily

local ++n_tests
assert _rc == 0
local ++n_passed
di as result "  PASSED"


// ═════════════════════════════════════════════════════════════════════════
// TEST 15: pred_left < pred_right when discontinuity is positive
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 15: pred_left < pred_right for positive jump ═══"

clear
set seed 15151
set obs 500
gen unit = 1
gen x = int(runiform() * 100)
gen above = (x >= 50)
// Large positive jump: outcome increases by 20 at cutoff
gen y = 10 + 20 * above + rnormal(0, 1)

rdtp y x, by(unit) bandwidth(50)

matrix R = r(results)

local ++n_tests
// pred_left (col 9) should be less than pred_right (col 10)
assert R[1,9] < R[1,10]
local ++n_passed
di as result "  PASSED"


// ═════════════════════════════════════════════════════════════════════════
// TEST 16: saving() with replace — first save errors without replace,
//          second save succeeds with replace
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 16: saving() with replace ═══"

clear
set seed 16161
set obs 200
gen unit = 1
gen x = int(runiform() * 100)
gen above = (x >= 50)
gen y = 5 + 6 * above + rnormal(0, 2)

tempfile savefile
// First save — should succeed
rdtp y x, by(unit) saving(`savefile') bandwidth(50)

// Second save without replace — should error
capture rdtp y x, by(unit) saving(`savefile') bandwidth(50)
local save_rc = _rc

// Second save with replace — should succeed
rdtp y x, by(unit) saving(`savefile') replace bandwidth(50)

local ++n_tests
assert `save_rc' != 0    // without replace must fail
assert _rc == 0           // with replace must succeed
local ++n_passed
di as result "  PASSED"


// ═════════════════════════════════════════════════════════════════════════
// TEST 17: if/in qualifiers
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 17: if/in qualifiers ═══"

clear
set seed 17171
set obs 400
gen school = ceil(_n / 100)
gen x = int(runiform() * 100)
gen above = (x >= 50)
gen y = 5 + 6 * above + rnormal(0, 2)

// Using if: restrict to schools 1 and 2 only
rdtp y x if school <= 2, by(school) bandwidth(50)

local ++n_tests
assert r(n_units) == 2    // only 2 of 4 schools
assert r(n_found) == 2
local ++n_passed
di as result "  PASSED (if qualifier)"

// Using in: restrict to first 200 observations (schools 1 and 2)
rdtp y x in 1/200, by(school) bandwidth(50)

local ++n_tests
assert r(n_units) == 2
assert r(n_found) == 2
local ++n_passed
di as result "  PASSED (in qualifier)"


// ═════════════════════════════════════════════════════════════════════════
// TEST 18: level() option
// ═════════════════════════════════════════════════════════════════════════

di _n as text "═══ TEST 18: level() option ═══"

clear
set seed 18181
set obs 200
gen unit = 1
gen x = int(runiform() * 100)
gen above = (x >= 50)
gen y = 5 + 6 * above + rnormal(0, 2)

// Run with 90% confidence level
rdtp y x, by(unit) bandwidth(50) level(90)

local ++n_tests
assert r(n_found) == 1
assert r(level)   == 90
local ++n_passed
di as result "  PASSED"


// ═════════════════════════════════════════════════════════════════════════
// SUMMARY
// ═════════════════════════════════════════════════════════════════════════

di _n as text "{hline 60}"
di as text "Test suite complete: " as result "`n_passed'" as text " of " as result "`n_tests'" as text " tests passed"
di as text "{hline 60}"

if `n_passed' == `n_tests' {
    di _n as result "ALL TESTS PASSED" _n
    exit 0
}
else {
    local n_failed = `n_tests' - `n_passed'
    di _n as error "`n_failed' TEST(S) FAILED" _n
    exit 9
}
