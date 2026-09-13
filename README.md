# recruitlag

**Can episodic recruitment be a delayed consequence of continuous reproduction?**

If recruitment is the delayed consequence of reproduction, then whatever the
delay, expected recruitment is a *convolution* of the reproductive record:
non-negative weights applied to reproduction at successive past periods. A
convolution is a weighted average, and a weighted average is never more
concentrated in time than the series it averages. Reproduction's own
concentration is therefore a ceiling on the concentration of expected
recruitment, attained only by a pure delay. Where reproduction is
continuous, the ceiling is low, every delay predicts a near-constant
expected series, and an episodic recruitment record is not something any
delay can be tracking.

`recruitlag` states that bound, searches families of lag weights to confirm
it on real data, calibrates the pooled comparison against negative-binomial
counts (which shows what the pooled comparison cannot refute), and runs the
test at the level at which the lag hypothesis is asserted: expected
recruitment on each replicate unit is that unit's own lagged reproduction.
It is the software behind the companion paper on the epiphytic orchid
*Lepanthes eltoroensis*, whose census data are included.

## Installation

```r
# install.packages("remotes")
remotes::install_github("raymondltremblay/recruitlag")
```

## The analysis in six calls

The data are tibbles with one row per period, or per unit and period, and
the functions take them directly.

```r
library(recruitlag)
lepanthes_census        # 12 six-monthly censuses; recruits scored at 3 to 12
lepanthes_hosts         # the same survey, one row per host tree and census

# 0. the census schedule, before anything is fitted: are the intervals equal,
#    and can a visit of that length see the thing being counted at all?
rain_intervals(my_survey_dates, persistence = 1.5)   # 1.5 = a Lepanthes fruit, in months

# 1. the ceiling: reproduction's own concentration, confirmed by a kernel search
k  <- lag_kernels(K = 4, n_dirichlet = 2000, bin = 6, unit = "mo")
ce <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4, kernels = k)
ce
rain_strip(`six-monthly 1999-2004` = ce)

# 1b. an interval on the exceedance, three ways, so that you take the one whose
#     logic you accept: a cluster bootstrap over host trees (a confidence interval
#     and a one-sided bootstrap p-value), a Bayesian bootstrap (a posterior under
#     a flat prior on the distribution of hosts), or a negative-binomial model of
#     the counts (a posterior on the concentration of the EXPECTED series, which
#     is the object the bound applies to; needs cmdstanr)
lag_ceiling_boot(lepanthes_hosts, unit = "host", reproduction = "inflorescences", K = 4, kernels = k)
lag_ceiling_bayesboot(lepanthes_hosts, unit = "host", reproduction = "inflorescences", K = 4, kernels = k)
lag_ceiling_stan(lepanthes_hosts, unit = "host", reproduction = "inflorescences", K = 4, kernels = k)

# 2. what the pooled comparison can refute (nothing, here)
ceiling_calibration(ce)

# 3. the host-level test, with the null stated as the lag hypothesis asserts it
ht <- host_lag_test(lepanthes_hosts, unit = "host", reproduction = "inflorescences", K = 4,
                    kernels = lag_kernels(4, bin = 6, unit = "mo",
                                          extra = list(fitted = lepanthes_profile$weight)))
ht
ht$table            # one row per lag profile
plot(ht)

# 4. the picture: profiles that could not differ more give one prediction
pr <- list(`one interval only (projection-matrix null)` = c(1, 0, 0, 0, 0),
           `concentrated early` = c(.7, .2, .07, .03, 0),
           `evenly spread` = rep(.2, 5),
           `concentrated late` = c(0, .03, .07, .2, .7))
rain_profiles(pr, bin = 6, unit = "months")
rain_expected(ce, pr)
```

Your own record goes in as a tibble with columns `unit`, `period` (an
integer from 1), `reproduction` and `recruits` (`NA` where recruits were not
scored); other column names are passed with `unit = `, `period = `,
`reproduction = ` and `recruits = `.

## Which measure of reproduction, and which measure of recruits

Neither choice is a detail, and both are usually made in silence. On real
records each can move the exceedance several fold on otherwise identical
data, so `rain_indices()` reports every candidate at once and the whole
table belongs in a write-up.

Two forces drive the differences and they pull in opposite directions.

**A stock gives a lower ceiling than a flux.** Adults change only by
recruitment minus death, so an adult series is smooth whatever the plants
are doing; the exceedance computed against it is the largest available. It
is the index most favourable to the argument this package exists to make,
and therefore the one to trust least. The same holds on the recruit side:
an untagged seedling count is a stock, and it inflates the observed
concentration rather than the ceiling.

**Measurement error runs the other way and protects the lag hypothesis.** A
structure caught only sometimes looks more erratic than the truth, which
raises its concentration, raises the ceiling and lowers the exceedance. So a
low exceedance computed from a badly observed column is not evidence that a
delay survives. A refutation from one still is.

The criterion that decides "well observed" is how long the structure
persists relative to the census interval, which is what
`rain_intervals(persistence = )` computes. A *Lepanthes* fruit lasts about
a month and a half: against a monthly census the ratio is above 1 and every
fruit is counted; against a six-monthly census it is 0.25 and three fruits
in four were never seen.

## Census conventions

Lag weights are indexed from the first lagged bin. With `lag0 = 1` (the
default, and the convention of a projection matrix) the first bin is the
period preceding the census at which recruits are scored; with `lag0 = 0`
it is the recruit census itself. Lags that fall before the record began are
read as the first period observed (`missing = "backfill"`); zero-filling is
not offered, because it asserts that no seed fell before the record began
and turns a long pure delay into a prediction of no recruits at the first
censuses.

## Two data sets, five vignettes

`lepanthes_monthly`, `lepanthes_census`, `lepanthes_hosts` and
`lepanthes_profile` hold the census records of the paper; `regimes` holds six
invented records in one long tibble (seasonal, intermediate and aseasonal reproduction, each
with a true delay and with gated recruitment) that show where the test has
power and where it does not. The vignettes are `baby-steps` (every step,
no mathematics), `invented-regimes`, `lepanthes`, `intervals` (what a confidence
interval, a credible interval, BCa, percentile, equal-tailed and highest-density
intervals are, and how to write the sentence for each) and `mathematics` (the
theory and the functions that implement it): `browseVignettes("recruitlag")`.

## The function inventory

| | |
|---|---|
| Schedule | `rain_intervals` |
| Kernels | `lag_kernels`, `flat_kernel` |
| The check | `lag_ceiling`, `lag_ceiling_by`, `rain_indices`, `ceiling_calibration`, `host_lag_test`, `host_matrices` |
| Uncertainty on the exceedance | `lag_ceiling_boot` (cluster bootstrap, BCa or percentile), `lag_ceiling_bayesboot` (Bayesian bootstrap), `lag_ceiling_stan` (negative-binomial model of the recruits, via `cmdstanr`); `ceiling_draws_table` (one table for several records and routes) and `rain_forest` (the forest plot of them) |
| Quantities | `rain_gini`, `rain_cv`, `concentration`, `phi_moment`, `convolve_lag`, `expected_recruits`, `series_from` |
| Figures | `rain_series`, `rain_strip`, `rain_strips`, `rain_profiles`, `rain_expected`, `rain_units`, and `plot()` methods for `ceiling_calibration` and `host_lag_test` |

Every figure function is prefixed `rain_` rather than `plot_`, so nothing in
the package collides with base R or with another package's `plot_*`. Both
concentration indices carry the prefix for the same reason: bare `gini` and
`cv` collide.

Each documented index carries a **Scale and how to read it** section giving
its bounds, what a value means, and where any cutoff comes from.

## Reproducing the paper

The tests in `tests/testthat` reproduce the archived outputs of the
companion paper's pipeline (phases 9a, 9b and 9d): the concentration
indices, the searched ceilings and exceedances exactly, and the Monte Carlo
probabilities within simulation error. Run them with `devtools::test()`.

## Citation

Tremblay, R. L. (2026). recruitlag: Can Episodic Recruitment Be a Delayed
Consequence of Continuous Reproduction? R package version 0.1.0.
