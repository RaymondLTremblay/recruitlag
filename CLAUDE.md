# CLAUDE.md: recruitlag

The fast on-ramp for working on this package. For what the method *is*, read
`README.md` and then `vignettes/mathematics.Rmd`. This file is about the code.

Author: Raymond L. Tremblay (RLT). Last revised 2026-09-13 (uncertainty functions added).

## 1. What the package does, in one paragraph

If recruitment is the delayed consequence of reproduction, expected recruitment
is a convolution of the reproductive record. A convolution is a weighted
average, and a weighted average is never more concentrated in time than the
series it averages. So reproduction's own concentration is a **ceiling** on the
concentration of expected recruitment, attained only by a pure delay. The
package states that bound, confirms it by searching families of lag weights,
calibrates the pooled comparison against negative-binomial counts, and runs the
test at the level at which the hypothesis is actually asserted: expected
recruitment on each replicate unit is that unit's own lagged reproduction.

Two diagnoses come out of it, and they are different things:

1. **The ceiling is exceeded.** Recruitment is more concentrated than any delay
   of that reproductive record can predict. A refutation.
2. **The lag is not identifiable.** Reproduction is nearly flat, so every lag
   profile predicts nearly the same expected series and the data cannot say
   which is right. Not a refutation, an absence of power.

## 2. Layout

```
R/
  ceiling.R        lag_ceiling(): the central object
  by_unit.R        lag_ceiling_by(), rain_indices()
  calibration.R    ceiling_calibration()
  host_test.R      host_lag_test(), host_matrices()
  uncertainty.R    lag_ceiling_boot(), lag_ceiling_bayesboot(), the shared
                   interval arithmetic and the lag_ceiling_draws methods
  stan.R           lag_ceiling_stan(); the model is inst/stan/ceiling_nb.stan
  concentration.R  rain_gini(), rain_cv(), concentration(), phi_moment()
  convolve.R       convolve_lag(), expected_recruits(), series_from()
  kernels.R        lag_kernels(), flat_kernel()
  intervals.R      rain_intervals()
  plots.R          every rain_* figure and the plot() methods
  checks.R         EVERY user-facing error and warning. See §5.
  data.R           documentation for the bundled data
data/              lepanthes_census, lepanthes_monthly, lepanthes_hosts,
                   lepanthes_profile, regimes
vignettes/         baby-steps, invented-regimes, lepanthes, mathematics
tests/testthat/    reproduce the companion paper's archived outputs
```

## 3. The exported functions

**Setting up**

| Function | Does |
|---|---|
| `lag_kernels(K, ...)` | build the family of lag weights to search: pure delays, windows, geometric decay, Dirichlet draws. `bin` and `unit` only label the bins, they do not rescale anything |
| `flat_kernel(K)` | the uniform kernel, on its own |
| `rain_intervals(data, date, period, persistence)` | **check the census schedule before anything else.** Reports each interval against the median, refuses on out-of-order dates, and with `persistence` gives the fraction of a structure that a visit of that length can see |

**The check**

| Function | Does |
|---|---|
| `lag_ceiling()` | the ceiling, the searched ceiling, the observed concentration, the exceedance. The central object; almost everything else takes one |
| `lag_ceiling_by(unit = )` | the same, once per replicate unit. Drops rows with neither reproduction nor recruits as unwatched periods before re-indexing |
| `rain_indices(reproduction = c(...))` | one row per candidate reproductive column. **The sensitivity analysis, not a menu** |
| `ceiling_calibration(ce)` | what counts about a mean obeying the ceiling could produce. The pooled comparison usually cannot refute anything |
| `host_lag_test(unit = )` | the test with the null stated as the hypothesis asserts it |
| `host_matrices()` | the per-unit matrices behind it |

**Uncertainty on the exceedance.** Three routes, all returning a
`lag_ceiling_draws` object, because RLT will not impose one philosophy on
users: each says what it estimates and what its interval means.

| Function | Estimand | Interval | Probability statement |
|---|---|---|---|
| `lag_ceiling_boot(unit = )` | exceedance of the realised counts | BCa (default) or percentile, 90% | `p_boot`, the one-sided bootstrap p-value for exceedance <= 1. **Not a posterior probability**, and the print method says so |
| `lag_ceiling_bayesboot(unit = )` | the same | equal-tailed (default) or HDI | `prob`, P(exceedance > 1) under a flat Dirichlet prior on the distribution of units (Rubin 1981) |
| `lag_ceiling_stan(unit = )` | exceedance of the **expected** series, from a negative-binomial model of both records with unit and period effects | equal-tailed or HDI | `prob`, as above. Usually lower than the other two, because sampling noise is removed. Needs `cmdstanr` (Suggests) and CmdStan; compiles into `tools::R_user_dir("recruitlag", "cache")` on first use; refuses rates because the NB is a distribution for counts |

All three resample or model **units**, so a pooled series without a unit
column is refused with a message pointing to `lag_ceiling()`. The BCa
acceleration comes from the leave-one-unit-out estimates. The HDI is offered
but documented as not invariant to a change of scale (the HDI of a ratio
and of its log differ), which is why the equal-tailed interval is the
default. On the *L. eltoroensis* hosts the cluster bootstrap gives
exceedance 2.56 with 90% BCa [2.1, 3.2] and no replicate at or below 1.

**Quantities**

`rain_gini()`, `rain_cv()`, `concentration()` (both at once), `phi_moment()`
(negative-binomial clumping, read through `1/phi`), `convolve_lag()`,
`expected_recruits()`, `series_from()`.

**Figures**, all prefixed `rain_`, never `plot_`:

`rain_series()`, `rain_strip()` (errors above 12 series), `rain_strips()`
(grouped and chunked), `rain_profiles()`, `rain_expected()`, `rain_units()`
(one point per unit, for many units), plus `plot()` methods for
`ceiling_calibration` and `host_lag_test`.

## 4. Naming rules, both standing

- **Figures are `rain_*`, never `plot_*`.** `plot` is too close to base R and
  the prefix has to be unique across packages. The two `plot()` methods are S3
  methods on the package's own classes, which is a different thing.
- **`gini` and `cv` are `rain_gini` and `rain_cv`.** Both bare names collide.
  When renaming, grep for the function *object* too: `apply(M, 1, gini)` broke
  three vignettes once because the rename missed unquoted references.

## 5. Errors and warnings live in `checks.R`

Every user-facing message is written there, and the standard is high on
purpose: RLT's instruction is that a message must say what is wrong, in which
column, and what to do, never the usual opaque R error.

- All messages use `call. = FALSE`.
- Helpers: `rl_abort`, `rl_warn`, `rl_list`, `rl_near` (did-you-mean, via
  `utils::adist`), `rl_check_columns`, `rl_check_long`, `rl_check_kernels`,
  `rl_check_class`.
- A wrong column name lists the columns the data frame does have and suggests
  the nearest.
- **Never swallow an error into a result cell.** `rain_indices()` catches the
  all-zero-reproduction case and reports it in the `best_kernel` column,
  because that is a property of that candidate column. Anything else it
  re-raises, because a call-level or data-level error is wrong for every
  candidate and would be mislabelled as "lagged reproduction is zero". That
  distinction was added after a period-index error spent an hour disguised as
  a biological result.

## 6. Gotchas, all of them learned the hard way

- **`period` must be an integer from 1 with no gaps.** Rows are not periods.
  A missed census is a row of `NA`, not an absent row, because lag weights
  count periods. `lag_ceiling()` errors if the index does not start at 1.
- **A stock is not a flux, on either side.** Adults change only by recruitment
  minus death, so an adult series is smooth whatever the plants do; its ceiling
  is the lowest available and its exceedance the largest, which makes it the
  index most flattering to the paper's own argument and therefore the one to
  trust least. The same applies to recruits: an untagged seedling count is a
  stock and inflates the observed concentration.
- **Measurement error runs the other way and protects the lag hypothesis.** A
  structure caught only sometimes looks more erratic than it is, which raises
  its concentration, raises the ceiling and lowers the exceedance. So a low
  exceedance from a badly observed column is not evidence that a delay
  survives. A refutation from one is still trustworthy.
- **Persistence over interval is the criterion for "well observed".** A
  *Lepanthes* fruit lasts about 1.5 months: against a monthly census the ratio
  is 1.5 and every fruit is seen; against six-monthly it is 0.25 and three in
  four are missed. `rain_intervals(persistence = )` computes it.
- **`NA` reproduction is read as zero and warns.** A period not measured and a
  period with no reproduction are different claims and the second lowers the
  ceiling. Drop unwatched rows before calling.
- **Key on a value, never a position.** Repeated source of bugs in this work,
  in the package and in the extraction scripts both.
- **Sentinel values.** Field workbooks code missing data as `999`, `-999` or
  `Nsurv`. The strings fail loudly; `999` does not, and read as a count it is
  catastrophic. The package cannot detect it. Recode at extraction.

## 7. Building and checking

The Stan route cannot be tested in the Cowork cloud container: the proxy
blocks r-universe, CRAN and GitHub, so neither `cmdstanr` nor CmdStan can be
installed there. `test-uncertainty.R` skips the Stan test when they are
absent. Run `lag_ceiling_stan()` on the Mac after any change to `stan.R` or
to `inst/stan/ceiling_nb.stan`.

R is not assumed present in every environment used for this project; where it
is, the full cycle is:

```sh
Rscript -e 'roxygen2::roxygenise(".")'     # after any roxygen change
R CMD INSTALL .                             # from a terminal, NOT from an R session
R CMD build . && R CMD check --no-manual recruitlag_0.1.0.tar.gz
```

Installing from inside an R session that has the package loaded produces a
corrupt lazy-load database. Restart R after installing.

The vignettes are the real integration tests: they run every exported function
on real data during `R CMD build`, so a rename that misses a reference fails
there rather than silently.

## 8. What the package has been exercised on

Beyond the bundled *L. eltoroensis* data, the companion project
(`lankesteriana/` in the `Lepanthes_eltoroensis_seed_recruitment` repository)
runs it on *Caladenia valida*, *Encyclia bocourtii*, *Lepanthes rubripetala*
and a ten-year *Lepanthes* metapopulation. Each exposed something:

| Data set | What it forced |
|---|---|
| *Caladenia valida* | units with zero reproduction throughout; unwatched leading periods being read as zero-reproduction periods |
| *Encyclia bocourtii* | a candidate column present in one year only; too many units for one strip, hence `rain_units()` and `rain_strips()` |
| *L. rubripetala* | the recruit-definition problem: 29 of 43 first sightings were adults, so most "recruits" were missed plants |
| *Lepanthes* metapopulation | unequal and undated intervals, hence `rain_intervals()`; `999` as a missing code; a recruit stock rather than a flux |

## 9. Standing conventions

- No em-dashes in any output: prose, documentation, figure captions, code
  comments. Commas, colons, semicolons or parentheses instead.
- Species names italicised everywhere, including figure text.
- Every documented index gets a **Scale and how to read it** section saying
  what the bounds are, what a value means, and where the cutoffs come from.
- No invented references or thresholds. If a number cannot be supported by a
  peer-reviewed source on hand, leave it blank and flag it; RLT will find the
  source.
