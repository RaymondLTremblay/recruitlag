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

## The analysis in five calls

The data are tibbles with one row per period, or per unit and period, and
the functions take them directly.

```r
library(recruitlag)
lepanthes_census        # 12 six-monthly censuses; recruits scored at 3 to 12
lepanthes_hosts         # the same survey, one row per host tree and census

# 1. the ceiling: reproduction's own concentration, confirmed by a kernel search
k  <- lag_kernels(K = 4, n_dirichlet = 2000, bin = 6, unit = "mo")
ce <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4, kernels = k)
ce
plot_ceiling_strip(`six-monthly 1999-2004` = ce)

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
plot_profiles(pr, bin = 6, unit = "months")
plot_expected(ce, pr)
```

Your own record goes in as a tibble with columns `unit`, `period` (an
integer from 1), `reproduction` and `recruits` (`NA` where recruits were not
scored); other column names are passed with `unit = `, `period = `,
`reproduction = ` and `recruits = `.

## Census conventions

Lag weights are indexed from the first lagged bin. With `lag0 = 1` (the
default, and the convention of a projection matrix) the first bin is the
period preceding the census at which recruits are scored; with `lag0 = 0`
it is the recruit census itself. Lags that fall before the record began are
read as the first period observed (`missing = "backfill"`); zero-filling is
not offered, because it asserts that no seed fell before the record began
and turns a long pure delay into a prediction of no recruits at the first
censuses.

## Two data sets, four vignettes

`lepanthes_monthly`, `lepanthes_census`, `lepanthes_hosts` and
`lepanthes_profile` hold the census records of the paper; `regimes` holds six
invented records in one long tibble (seasonal, intermediate and aseasonal reproduction, each
with a true delay and with gated recruitment) that show where the test has
power and where it does not. The vignettes are `baby-steps` (every step,
no mathematics), `invented-regimes`, `lepanthes` and `mathematics` (the
theory and the functions that implement it): `browseVignettes("recruitlag")`.

## Reproducing the paper

The tests in `tests/testthat` reproduce the archived outputs of the
companion paper's pipeline (phases 9a, 9b and 9d): the concentration
indices, the searched ceilings and exceedances exactly, and the Monte Carlo
probabilities within simulation error. Run them with `devtools::test()`.

## Citation

Tremblay, R. L. (2026). recruitlag: Can Episodic Recruitment Be a Delayed
Consequence of Continuous Reproduction? R package version 0.1.0.
