# recruitlag 0.1.0

First release. The software behind the companion paper on *Lepanthes eltoroensis*
(Tremblay, submitted to *Oikos*, 2026).

## Behaviour settled in the pre-release review (2026-10-05)

* `expected_recruits()` pairs each unit's recruits with its own reproduction by row
  name, whatever the row order of `X`.
* With `missing = "drop"`, `lag_ceiling()` leaves the scored periods whose lag window
  reaches before the record out of the observed series too, so that the observed and
  the expected concentration are measured on the same periods (`n_dropped` says how
  many).
* A flat reproductive record gives an exceedance of `NA` (with `flat = TRUE`), not
  `Inf`; `rain_gini()` cannot return a negative value from rounding.
* `lag_ceiling_by()` applies `min_periods` to scored periods, and reports the per-unit
  warnings once with the units they concern instead of silencing them.
* `convolve_lag()` refuses a record too short for the window and validates `t_R` and
  `lag0`; lag weights are validated everywhere a profile is applied.
* `recruit_triage()` refuses a long data frame with a misnamed column instead of reading
  it as a wide matrix, and refuses `species =` together with `s`/`r`.
* The bootstrap and posterior probabilities use the finite replicates of their own
  statistic; the print method no longer fails when none is finite.
* `host_matrices()` and `series_from()` warn when missing recruits are read as zero.
* The Stan diagnostics cover `ER`, the expected series that is reported, and the period
  effects. CmdStan is looked for after the data are validated.
* Every user-facing check lives in `checks.R`; the all-zero-reproduction diagnosis is a
  classed condition (`recruitlag_zero_lagged`).
* The introductory vignette is `recruitlag` (formerly `baby-steps`). Lag bins are
  labelled from the first lagged census (0-6, 6-12, ..., 24-30 months) in all prose.

## The check

* `lag_ceiling()`: the concentration ceiling that any delay of a reproductive record
  implies, the searched ceiling over four families of lag weights (pure delays, windows,
  geometric decay, Dirichlet draws; `lag_kernels()`), the observed concentration of the
  recruit series, and their ratio, the exceedance.
* `lag_ceiling_by()`: the same for every replicate unit; `rain_indices()`: the ceiling
  under every candidate reproductive column, reported as one table and never as a menu.
* `ceiling_calibration()`: what negative-binomial counts about a mean that obeys the
  ceiling can produce, which shows what the pooled comparison cannot refute.
* `host_lag_test()` and `host_matrices()`: the test at the level at which the lag
  hypothesis is asserted, expected recruitment on each unit being that unit's own lagged
  reproduction scaled to its own total, with the negative-binomial size moment-matched about
  each kernel's own expectations. Both refuse non-count recruit series.

## Uncertainty on the exceedance

* Three routes returning one object class, `lag_ceiling_draws`: `lag_ceiling_boot()`
  (cluster bootstrap over units, BCa or percentile interval, one-sided bootstrap p),
  `lag_ceiling_bayesboot()` (Bayesian bootstrap, equal-tailed or HDI, P(exceedance > 1)),
  and `lag_ceiling_stan()` (negative-binomial model of the recruits with unit and period
  effects, the exceedance of the expected series; needs `cmdstanr` and CmdStan). Each
  print method says what it estimates and what its interval means.
* `ceiling_draws_table()` and `rain_forest()`: one table and one forest plot for several
  records and routes.

## Is a first sighting a recruit?

* `recruit_triage()` and `rain_sightings()`: the state at first sighting, the runs of
  unseen periods between sightings, and the expected number of first sightings that are
  returns of plants present but unseen, calibrated on the published *Caladenia* survival
  and resighting rates (`caladenia_dormancy`).

## Census conventions

* `rain_intervals()`: the census schedule against its median interval, with the fraction
  of a structure that a visit can see for a given persistence.
* `period` is an integer from 1 with no gaps; a missed census is a row of `NA`, not an
  absent row; `NA` reproduction is read as zero with a warning.

## Quantities and figures

* `rain_gini()`, `rain_cv()`, `concentration()`, `phi_moment()`, `convolve_lag()`,
  `expected_recruits()`, `series_from()`.
* `rain_series()`, `rain_strip()`, `rain_strips()`, `rain_profiles()`, `rain_expected()`,
  `rain_units()`, and `plot()` methods for `ceiling_calibration` and `host_lag_test`.

## Data and vignettes

* Bundled: `lepanthes_census`, `lepanthes_hosts`, `lepanthes_monthly`, `lepanthes_profile`,
  `regimes`, `regimes_truth`, `caladenia_dormancy`.
* Vignettes: `recruitlag` (introduction), `invented-regimes`, `lepanthes`, `mathematics`,
  `intervals`.
* Every user-facing message says what is wrong, in which column, and what to do.
