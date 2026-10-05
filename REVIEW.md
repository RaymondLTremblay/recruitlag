# REVIEW.md: recruitlag 0.1.0, pre-release review

Review date: 2026-10-05. This review is based on `master` at `295a9df` plus the
uncommitted changes to `.Rbuildignore` and `CLAUDE.md`. Nothing in the package
has been changed yet. Each item has a number so that you can tell me which to
apply ("apply B1-B6, S3, skip S10").

Severity: **B** = bug (wrong result or a crash), **C** = affects CRAN
acceptance, **S** = should fix, **M** = minor.

## 0. Check results (cloud container, R 4.3.3, Ubuntu 24.04)

* `roxygenise` + `R CMD build` + `R CMD check --no-manual`:
  **0 errors, 0 warnings, 1 note**. The note is that `cmdstanr` is suggested
  but not available. 132 tests pass and 3 are skipped (Stan and
  `skip_on_cran`).
* The incoming-feasibility check produced one additional note: the Title is
  not in title case ("Be" should be "be"). See C1.
* **`man/lepanthes.Rd` on the Mac is stale.** It does not match what
  roxygen generates from the current `R/data.R`. Run
  `roxygen2::roxygenise(".")` before building the 0.1.0 release.
* Slowest example: `plot.ceiling_calibration`, 2.3 s, under the 5 s limit.
* Vignette render times (one core): `baby-steps` 22 s, `intervals` 9 s,
  `invented-regimes` 46 s, `lepanthes` 35 s, `mathematics` 2 s; about 2 min
  in total. CRAN accepts this, but it is on the long side (see S12).
* Warnings appear inside the rendered vignettes:
  `invented-regimes` shows a dplyr 1.1 deprecation warning (S11) and
  "Removed 30 rows containing missing values"; `baby-steps` shows "Removed
  5 rows containing missing values".

## 1. Vignette rename: `baby-steps` to `recruitlag`

New file `vignettes/recruitlag.Rmd`, title "Introduction to recruitlag: a
worked analysis, step by step". Because the file has the same name as the
package, a pkgdown site will show it as "Get started". These are the edits:

| # | Where | Current | New |
|---|---|---|---|
| R1 | file | `vignettes/baby-steps.Rmd` | `vignettes/recruitlag.Rmd` (`git mv`) |
| R2 | line 2, `title:` | "Baby steps: running the test one line at a time" | "Introduction to recruitlag: a worked analysis, step by step" |
| R3 | line 5, `VignetteIndexEntry` | same as line 2 | same as R2 |
| R4 | lines 14-15 | "This vignette is for someone who has never used the package, and maybe has not used R much either." | "This vignette is for a reader new to the package." |
| R5 | lines 17-18 | 'read the vignette called "The mathematics behind the test" afterwards' (the real title continues ", and the functions that implement it") | 'read `vignette("mathematics", package = "recruitlag")` afterwards' |
| R6 | `README.md:134` | "The vignettes are `baby-steps` (every step, ..." | "`recruitlag` (the introduction: every step, ..." |
| R7 | `NEWS.md:57` | "Vignettes: `baby-steps`, ..." | "Vignettes: `recruitlag` (introduction), ..." |
| R8 | `CLAUDE.md` §2 | "baby-steps, invented-regimes, lepanthes, mathematics" | "recruitlag, invented-regimes, lepanthes, mathematics, intervals" (`intervals` is also missing from that line now) |

No `vignette("baby-steps")` calls and no roxygen references were found.
Optional tone edits (your wording kept, only the classroom register removed):

* l. 34 "Nothing is printed. That is correct: loading a package is silent." to "Loading a package prints nothing."
* l. 99-100 "Write those two numbers down." to "Note those two numbers."
* l. 115 "Read the printout slowly." to "Read the printout line by line."
* l. 185 "Check yourself with a record..." to "A check with a record..."
* l. 242 "Here is the important part, so read it twice." to "This is the important point."

Also see V4: the paragraph at lines 115-120 is factually wrong.

## 2. Bugs in the code

All verified by running the code, except where noted.

**B1. `expected_recruits()` matches units by position, not by name** (`R/host_test.R:69`).
`check_host_inputs()` reorders `X` to the row names of `R` and returns the
result as `chk$X`, but the function goes on to use the original `X`. If you
pass `X` with its rows in a different order from `R`, every unit's recruits
are paired with another unit's reproduction (mean relative difference 0.63
on `lepanthes_hosts`). The documentation says the order does not matter.
`host_lag_test()` is not affected. Fix: add `X <- chk$X` after line 69.

**B2. `missing = "drop"` computes observed and ceiling over different periods** (`R/ceiling.R:127-144`).
The kernel outputs lose their `NA` periods but `R` keeps all of its periods,
so the two Ginis are measured over different numbers of periods (10 against 7
on `lepanthes_census`, K = 4). Fix: drop the same periods from `R` before
computing `obs`.

**B3. `ceiling_calibration()` on a `missing = "drop"` object** (`R/calibration.R:353`).
The mean of an expected series that contains `NA` is `NA`, so `rnbinom()`
returns `NA`s with the warning "NAs produced". Fix: B2, or refuse
`missing = "drop"` with `rl_abort`.

**B4. A flat reproductive record gives an `Inf` or absurd exceedance** (`R/ceiling.R:142-144`, `R/by_unit.R:153`).
`lag_ceiling(rep(0.37, 12), ...)` returns an exceedance of `Inf`. With
`rep(12.3456, 12)` rounding error gives a Gini of -2.2e-16 and an exceedance
of -2.9e15. The `NA` guard added on 2026-09-13 covers only
`as.data.frame.lag_ceiling_list()`. Fix: clamp `rain_gini()` at 0, test for a
flat record with a tolerance, and return `NA` with a note from
`lag_ceiling()` itself.

**B5. `min_periods` counts rows, not scored periods** (`R/by_unit.R:95`).
The documentation says "scored periods". A unit with 8 rows and 2 scored
periods passes `min_periods = 5`. Fix: compare `scored`.

**B6. `convolve_lag()` default `t_R` runs backwards on short records** (`R/convolve.R:57`).
`convolve_lag(1:3, rep(1, 5))` returns `NA NA 1.6 1.2`: `seq.int` counts
down and produces periods past the end of `X`. Fix: refuse when
`length(X) < lag0 + K + 1`, and validate user-supplied `t_R` and `lag0`
(whole numbers, inside the record).

**B7. `recruit_triage()` silently reads long data with a misnamed column as a wide matrix** (`R/dormancy.R:13`).
If any of `plant`, `period` or `state` is missing (for example `Plant`), the
data frame falls through to `as.matrix()`, and the ID, period and state
columns become three "periods". Fix: if some but not all of the three names
are present, abort with `rl_check_columns()`.

**B8. Bootstrap and posterior probabilities use the wrong rows, and print can crash** (`R/uncertainty.R:283-285, 380, 418`; `R/stan.R:189`).
`ok` is "complete over all six statistics", so a replicate where only the CV
fails is also dropped from the Gini `p_boot`. When recruits are scored at a
single period, every replicate is dropped, `p_boot` is `NaN`, and `print`
fails with "missing value where TRUE/FALSE needed". Fix: compute each
probability on the finite values of its own column, and store that n
for the "< 1/R" label.

**B9. `rain_units()` subtitle can read "NA of the N units"** (`R/plots.R:340`).
A unit with no recruits has an `NA` observed Gini. Fix:
`sum(..., na.rm = TRUE)`, and phrase it as "of the n units with recruits".

**B10. Stan diagnostics skip the quantity that is reported** (`R/stan.R:191-199`).
Rhat and ESS are checked on `a_r, s_ur, s_vr, phi_r` only, not on `ER` (the
expected series whose Gini is the reported exceedance) or on the period
effects. Fix: add `ER` and `v_r` to `pars`. (I read the code but could not run
it: Stan cannot be installed in the container.)

**B11. `host_matrices()` turns missing recruits into zeros without saying so** (`R/tidy.R:34`).
"Not censused at that period" becomes "zero recruits", which raises the
observed concentration. `series_from()` (`R/tidy.R:58`) likewise sums only
the units scored in a period. No effect on the bundled data. Fix:
`rl_warn` with the number of affected cells, and document it.

## 3. CRAN items

**C1. Title case.** Change "Be" to "be". The title is 77 characters long;
CRAN prefers 65 or fewer, but that is not enforced. The wording is yours
to choose.

**C2. Vignettes load Suggests packages unconditionally.** `baby-steps.Rmd:30,80`,
`invented-regimes.Rmd:16-17` and `lepanthes.Rmd:16` load `dplyr` and/or
`tidyr`. With `_R_CHECK_DEPENDS_ONLY_` these vignettes fail. Fix: put
`eval = requireNamespace("dplyr", quietly = TRUE)` (and `tidyr`) in each
setup chunk.

**C3. A `regimes` example uses `dplyr`** (`R/data.R:290-291`). Same issue as
C2. Fix: `regimes[regimes$record == "aseasonal_gated", ]`.

**C4. `lepanthes` and `regimes` have no `@format`** (`R/data.R`). Reviewers
ask about this. Fix: add `@format` blocks.

**C5. Data licence.** README and `R/data.R:37` say the data are CC BY 4.0,
but `License:` says MIT only. Fix: add a short `inst/COPYRIGHTS` stating the
data licence, or mention it in `LICENSE`.

**C6. DESCRIPTION Description.** The text leaves out the three
uncertainty routes and `recruit_triage()`. If 'Stan' or 'cmdstanr' are named,
they must be in single quotes. Once the paper has a DOI, add it in the form
`Tremblay (2026) <doi:...>`. The DOI is left blank here.

**C7. `utils::` and `tools::` are used but not declared.** Here this gave no
note, because both are base packages. Adding them to Imports is still
cleaner. Optional.

**C8. `.Rbuildignore` additions:** `^REVIEW\.md$`, `^cran-comments\.md$`,
`^CRAN-SUBMISSION$`, `^\.zenodo\.json$`, `^CITATION\.cff$`,
`^_pkgdown\.yml$`, `^docs$`.

## 4. Vignettes and README

**V1. `mathematics.Rmd:114`** refers to `lepanthes$lag_profile`, which does
not exist. It should be `lepanthes_profile`.

**V2. `lepanthes.Rmd:61`** says the monthly exceedance is "about 10". The
computed value is 9.00 (Gini). Fix: "about 9", or inline
`` `r round(c_m$exceedance["gini"], 1)` ``.

**V3. Hard-coded numbers in prose.** These can go stale:
`lepanthes.Rmd:60-62, 82-83, 143, 156, 173`; `invented-regimes.Rmd:105`
("about 0.65"); `intervals.Rmd:88-89, 117-119, 219` (2.56, 2.1 to 3.2, 0.999,
0.0005). Fix: compute them inline with `` `r ` ``.

**V4. `baby-steps.Rmd:115-120`.** The prose says "the last line divides one by
the other" (observed over the theorem bound). The exceedance actually divides
by the *searched* ceiling (`R/ceiling.R:143`), and the last printed line is the
note about `ceiling_calibration()`. Fix: "The line 'exceedance' divides the
observed Gini by the searched ceiling."

**V5. Two conventions for lag months.** `lepanthes.Rmd:44, 73-75, 195`
describes the bins as 6, 12, ..., 30 months ("6 to 30 months"). The kernel
labels from `lag_kernels(bin = 6)`, `lepanthes_profile$months` and
`mathematics.Rmd` describe them as 0-6 ... 24-30. The `convolve_lag()` example
(`R/convolve.R:54`) calls `c(0,0,0,0,1)` "a 24 to 30 month delay", and
the `R/kernels.R:12-15` scale section says "6 to 30 months". The labels
ignore `lag0`. **This is a decision for you:** should the labels include
`lag0`?

**V6. `baby-steps.Rmd:89-91`** says the Gini runs "between 0 and 1". The
documented maximum is (n - 1)/n.

**V7. Species not italicised:** the `lepanthes.Rmd` title (use
`"The *Lepanthes eltoroensis* analysis"`; the index entry cannot carry
italics); the forest-plot label "L. eltoroensis, 23 hosts"
(`intervals.Rmd:249`, and the `R/forest.R:71, 137` examples); record
names drawn as plain axis text in `rain_forest()` and `rain_strip()`.
Simplest fix: relabel as "23 host trees". Alternatively, add a
plotmath/italic option.

**V8. `README.md:127`** says "Two data sets", but there are three
(`caladenia_dormancy`). `NEWS.md:54` omits `regimes_truth`.

**V9. `data-raw/` is build-ignored but the docs point to it**
(`invented-regimes.Rmd:30, 202`, `R/data.R:59`). On CRAN those paths will not
exist. Fix: give the GitHub URL instead.

**V10. `data-raw/lepanthes.R:12-16` reads five CSVs that are not in
`data-raw/`.** The data cannot be rebuilt from the repository, which matters
for the Zenodo archive. Fix: commit them or state their archived location.

**V11. `mathematics.Rmd:191` function map** leaves out the uncertainty
functions, `recruit_triage()`, `rain_indices()` and `rain_intervals()`.

**V12. Phase numbers disagree:** README and `mathematics.Rmd` say 9b,
`test-calibration.R:1` says 9a.

## 5. Should fix (code)

**S1. A bare `stop()` outside `checks.R`** (`R/ceiling.R:136`), without
`call. = FALSE`. `rain_indices()` and `print.lag_ceiling_list()` recognise this
error by `grepl` on its English text, so rewording it would silently break
the all-zero diagnosis that §5 of CLAUDE.md protects. Fix: a classed
condition via `rlang::abort(class = "recruitlag_zero_lagged")`, tested with
`inherits()`.

**S2. Validators outside `checks.R`:** `check_counts()` and `need_cmdstan()`
(`R/stan.R`) and `check_level()` (`R/uncertainty.R:388`). `check_counts()`
hard-codes the column name `"recruits"` when called from
`calibration.R:85`, so the message names the wrong column. Fix: move them to
`checks.R` and pass the real column name.

**S3. `lag_ceiling_by()` suppresses per-unit warnings** (`R/by_unit.R:103`).
This hides "NA reproduction read as zero", and line 84 drops interior blank
rows, which turns a missed census into a zero-reproduction gap with no
warning. Fix: collect the warnings and report them once.

**S4. Opaque R errors on degenerate input** (`R/calibration.R:346`,
`R/host_test.R:195`): an all-zero recruit series, or `phi_moment()` with no
recruiting unit, ends in "missing value where TRUE/FALSE needed". Fix:
add `rl_abort` checks, and check that a fixed `phi` is positive.

**S5. Lag weights not validated.** `expected_recruits()` does not check `w`.
In `lag_kernels()`, a negative `rho` gives negative weights, and an `extra`
kernel that sums to 0 gives `NaN`. Fix: reuse the `convolve_lag()` weight
checks and require `rho > 0`.

**S6. `rain_forest()`** takes `d$level[1]` for the subtitle, so mixed levels are
mislabelled (`R/forest.R:173`). Its fixed ticks from 0.25 to 32 leave a blank
axis outside that range (`R/forest.R:159`).

**S7. `need_cmdstan()` runs before data validation** (`R/stan.R:146`). A user
without `cmdstanr` who also has a wrong column name sees only the install
message.

**S8. `recruit_triage(species = )` silently overwrites user `s`/`r`**
(`R/dormancy.R:164-171`). The function also assumes that one period equals
the time step of `s` and `r` (annual for `caladenia_dormancy`), which the
scale section should state.

**S9. `rain_strips()` with no finite key** gives a raw `max()` warning
(`R/plots.R:422`).

**S10. Documentation that disagrees with the code:** `host_test.R:254-257`
(the verdict maximum is per `phi_source`, not over all); `intervals.R:29`
("invisibly", but the result is visible); `by_unit.R:21` ("concentration 1",
but the maximum is (n - 1)/n); `by_unit.R:42` (the columns are
`ceiling_gini`/`observed_gini`); `recruitlag-package.R:199` (`p_boot`
can be 0, displayed as "< 1/R").

**S11. dplyr deprecation in a rendered vignette** (`invented-regimes.Rmd:97`):
`across(where(is.numeric), round, 2)`. Fix: `\(x) round(x, 2)`. The
"Removed rows" ggplot warnings in two vignettes would also be better silenced
or the `NA`s dropped before plotting.

**S12. Vignette run time** (about 2 min). The `intervals` vignette repeats the
`lepanthes` bootstrap runs. Option: smaller `nsim` when
`!identical(Sys.getenv("NOT_CRAN"), "true")`.

## 6. Thresholds and statements without a source on hand

Decision (RLT, 2026-10-05): all five are kept as working conventions of the
package, each documented as such where it is used, with no citation. Nothing
in the list below changes a number the package returns.

* `R/stan.R:199`: ESS warning at 100. Vehtari et al. (2021, *Bayesian
  Analysis* 16: 667-718) is the usual source for ESS and Rhat thresholds;
  please check it before citing.
* `R/intervals.R:39-40, 81, 235`: tolerance 0.2 (text says 0.8 and 1.25;
  code gives 0.8 and 1.2; print says "below about a half"). These are
  inconsistent and have no source. Also, an undated middle census makes
  the median `NA`, yet print still says "All intervals are within 20%".
* `intervals.Rmd:84`: the warning below ten units.
* `R/plots.R:283`: a dotted line at 0.05 in `plot.host_lag_test()`. This
  contradicts "The package applies no thresholds" (`recruitlag-package.R:31`).
* `mathematics.Rmd:66`: the Hardy, Littlewood and Pólya theorem needs a
  reference. Line 22, "distributed-lag or antecedent-effect model", also needs
  one. The vignette has no References section.
* `README.md:111-114`: the 1.5-month persistence of a *Lepanthes* fruit.
* `intervals.Rmd:257-267`: references listed but not cited in the text.

## 7. Tests

* Exported functions with no direct test: `flat_kernel`, `rain_indices`
  (including the re-raise behaviour), `rain_intervals`, `series_from`, all
  `rain_*` figures and the `plot()` methods. A smoke test that builds each
  figure (`expect_s3_class(p, "ggplot")`) would catch a broken rename.
* The success paths of `ceiling_calibration()` and `host_lag_test()` run only
  off CRAN. Add one small-`nsim` test on `regimes`.
* `test-calibration.R:19` draws Dirichlet kernels without a seed.
* `reference/verdict.csv` and `reference/zero_inflation.csv` are not read by
  any test.
* Each bug in §2 should get a regression test when it is fixed.

## 8. Metadata, minor

* URL spelling: `DESCRIPTION` and README use `raymondltremblay`;
  `inst/CITATION` uses `RaymondLTremblay`. Pick one.
* `inst/CITATION`: the version is typed in twice; use `meta$Version`. Add the
  Zenodo `doi =` once it is minted. `note = "submitted"` needs updating at
  the time of submission.
* `Authors@R`: add the `cph` role.
* `data-raw/regimes.R:21`: `T <- 24L` masks `TRUE`. Rename it `n_per`.
