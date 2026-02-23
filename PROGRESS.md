# Progress Log — rdtp Stata Package

## Objective
Rewrite the RD cutoff-search algorithm from McEachin, Domina, Penner (2020, JPAM) as a clean, general-purpose Stata .ado command (`rdtp`) with a help file.

## Plan
1. Read and understand original do-file
2. Create `rdtp.ado` — main command with all options
3. Create `rdtp.sthlp` — help file
4. Code review and testing

## Final Status (v1.0.0)
- **Step 1: COMPLETED** — Read and analyzed original do-file
- **Step 2: COMPLETED** — Created rdtp.ado (535 lines) with all options, proper preserve/restore, postfile, tempvars, r-class results
- **Step 3: COMPLETED** — Created rdtp.sthlp (362 lines) with full SMCL help: syntax, description, algorithm, options, stored results, examples, limitations, references
- **Step 4: COMPLETED** — Code review found 10 issues (3 medium, 4 low, 3 trivial), all applied:
  1. Removed duplicate savefile construction
  2. Moved tempvar outside inner loop
  3. Added invariant comment for n_found + n_skipped == n_units
  4. Added more special character cleaning for matrix row names
  5. Added return scalar level
  6. Left {c 178} as-is (standard Stata practice)
  7. Added {p_end} tags to all {phang} paragraphs
  8. Added Lee2010 citation in description
  9. Removed redundant {synoptset} directives
  10. Added {p2colreset} after each stored results subsection

## Post-Review Fixes (v1.0.1) — /datascience-reviewer findings

**Date:** 2026-02-23
**Review verdict:** Satisfied (all items actionable improvements, not blockers)

### Task
Address three actionable items identified by the `/datascience-reviewer` skill,
then perform a code review pass and fix any issues found.

### Plan
1. Add VCE pre-validation block to `rdtp.ado` (after postfile setup)
2. Add self-contained synthetic-data example to `rdtp.sthlp` (before existing examples)
3. Remove original do-file row from `README.md` Files table
4. Code review pass on all modified files

### Status
- **Step 1: COMPLETED** — Added VCE pre-validation block (lines 150-224 of rdtp.ado).
  Subsets to the first unit's actual data (`preserve` / `keep`), constructs the real
  RD variables (treatment indicator, centered forcing variable, interaction) at the
  first candidate cutoff within the search range, and runs a test `regress` with the
  user's `vce()` spec. This catches invalid VCE specifications (misspelled cluster
  variable, bad vcetype, etc.) before entering the main loop, rather than having every
  unit silently fail via `capture regress`. Includes `if !missing()` guards on all
  `gen` statements for consistency with the main loop. Also includes a fallback branch
  (lines 207-222) for the case where the first unit has no candidates in the search
  range: runs a simple `regress depvar forcing, vce(...)` to still validate the VCE
  spec. The postfile handle is properly closed (`postclose`) before any `exit` call.
- **Step 2: COMPLETED** — Added self-contained synthetic-data example (lines 253-273 of rdtp.sthlp).
  Creates 500 observations across 5 schools with known true cutoffs (300, 350, 400, 325, 375),
  generates a discrete test-score forcing variable (200-500) and an outcome with a 10-unit
  discontinuity, then runs `rdtp` with `noisily` and `searchrange(250 450)` so the user
  can see per-unit output. Placed before the existing examples so it appears first.
- **Step 3: COMPLETED** — Removed the `RD_algorithm_para_school_linear_spline_RR1.do`
  row from the Files table in README.md. The original do-file is still in the repo as
  a reference but no longer advertised in the readme.
- **Step 4: COMPLETED** — Code review found 6 issues; all fixed:
  1. (Low) Ensured Files section at bottom of PROGRESS.md shows current line counts
     (612 for rdtp.ado, 389 for rdtp.sthlp); v1.0.0 section retains historical counts
  2. (Low) Same as above for sthlp — cosmetic consistency fix
  3. (Medium) Added `if !missing()` guards to the VCE block's `gen` statements for
     the treatment indicator and centered forcing variable, matching the main loop's
     pattern and avoiding missing-value propagation edge cases
  4. (Medium) Added fallback branch (lines 207-222) for when the first unit has no
     candidates in the search range — runs a simple `regress depvar forcing, vce(...)`
     so the VCE spec is still validated even in that edge case
  5. (Low) Updated stale description in PROGRESS.md Step 1 (was "lines 150-163",
     now "lines 150-224" reflecting the expanded block with fallback branch)
  6. (Low) Updated README.md Files table line counts to match actual file lengths

### Files Modified
- `rdtp.ado` — Added ~77 lines (VCE pre-validation block, lines 150-224; 612 total)
- `rdtp.sthlp` — Added ~21 lines (synthetic-data example, lines 253-273; 389 total)
- `README.md` — Removed 1 row from Files table; updated line counts
- `PROGRESS.md` — This file (updated with post-review section and code review notes)

## Bugs Fixed from Original
1. -99 sentinels replaced with Stata missing (.)
2. `bob` variable pattern replaced with capture + scalar checks
3. pctL/pctR swap fixed with pred_left/pred_right and correct labels
4. Re-read from disk replaced with restore, preserve in memory
5. Incremental append replaced with postfile
6. Hardcoded paths/magic numbers replaced with configurable options

## Files
- `rdtp.ado` — Main command (CREATED v1.0.0, UPDATED v1.0.1, ~612 lines)
- `rdtp.sthlp` — Help file (CREATED v1.0.0, UPDATED v1.0.1, ~389 lines)
- `RD_algorithm_para_school_linear_spline_RR1.do` — Original (REFERENCE)
- `README.md` — Project readme (UPDATED)
- `PROGRESS.md` — This file (UPDATED)
