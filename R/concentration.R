#' Concentration of a series in time
#'
#' Indices of how unevenly a non-negative series is spread across its periods.
#' The Gini coefficient is zero when every period carries the same count and
#' approaches one when the whole total falls in a single period. The
#' coefficient of variation is the standard deviation divided by the mean.
#' Both are Schur-convex, which is what makes them usable with the ceiling
#' argument: a weighted average of a series can never be more concentrated,
#' by either index, than the series itself.
#'
#' @param x A numeric vector of non-negative counts or rates, one per period.
#'   `NA` values are dropped.
#'
#' @return `rain_gini()` and `rain_cv()` return a single number (`NA` when the series
#'   sums to zero). `concentration()` returns a named numeric vector with
#'   `n`, `mean`, `gini`, `cv`, `zeros` (number of periods with a zero count)
#'   and `max_share` (the largest single period's share of the total).
#'
#' @details The Gini coefficient uses the plug-in definition
#'   \eqn{G = 2 \sum_i i x_{(i)} / (n \sum_i x_i) - (n + 1) / n}, with
#'   \eqn{x_{(i)}} the values in increasing order. It is biased downward in
#'   short series (Deltas 2003), but the bias applies equally to the observed
#'   series and to the ceiling, so a short series cannot manufacture an
#'   exceedance.
#'
#' @section Scale and how to read it:
#'
#' **Lower end.** Both indices are 0 when every period carries the same
#' amount, and both rise as the total is moved into fewer periods. Neither
#' can be negative.
#'
#' **Upper end, which depends on the number of periods and not on the data.**
#' The largest value the plug-in Gini can take on `n` periods is
#' \eqn{(n-1)/n}, reached when one period carries the whole total: 0.900 on
#' the 10 scored six-monthly censuses, 0.955 on the 22 scored monthly
#' periods. The ceiling of 1 is a limit as `n` grows and is never attained by
#' a finite series (Damgaard and Weiner 2000). `rain_cv()` uses the `n - 1`
#' divisor of [stats::sd()], so its maximum on `n` periods is \eqn{\sqrt{n}}{sqrt(n)}:
#' 3.16 and 4.69 for the same two series. A Gini of 0.65 therefore means
#' something different on 10 periods than on 100, and values are comparable
#' only between series of the same length. Every comparison this package
#' makes, observed against ceiling and observed against simulated, is between
#' series of the same length for that reason.
#'
#' **Bias correction is deliberately not applied.** Damgaard and Weiner
#' (2000) recommend multiplying the sample Gini by \eqn{n/(n-1)} for an
#' unbiased estimate. `rain_gini()` returns the uncorrected value because every
#' use of it here is a ratio or a rank between series of equal length, in
#' which a common factor cancels. Multiply by \eqn{n/(n-1)} yourself if you
#' want an estimate of the population value.
#'
#' **No threshold.** The package supplies no value of either index above
#' which a series counts as episodic, and takes none from the literature.
#' The level of the index carries no verdict on its own. What carries the
#' inference is the comparison: the index of the observed series against the
#' ceiling ([lag_ceiling()]), and against series simulated under a mean that
#' obeys the ceiling ([ceiling_calibration()], [host_lag_test()]).
#'
#' **Why these two indices and not others.** Both are scale invariant, so
#' the unknown scale of expected recruitment never enters, and both are
#' Schur-convex: moving a quantity from a fuller period to an emptier one
#' lowers them, which is the principle of transfers (Atkinson 1970; Marshall,
#' Olkin and Arnold 2011, ch. 2). Schur-convexity is what makes the ceiling
#' argument go through (Hardy, Littlewood and Polya 1952). They are not
#' interchangeable: two series with the same Gini can have differently shaped
#' Lorenz curves, and Damgaard and Weiner (2000) propose the Lorenz asymmetry
#' coefficient to separate them. `concentration()` therefore reports `zeros`
#' and `max_share` alongside, which say plainly how much of the record is
#' empty and how much of the total falls in its largest period. `max_share`
#' runs from \eqn{1/n} (even) to 1 (everything in one period).
#'
#' @references
#' Atkinson, A. B. (1970) On the measurement of inequality. *Journal of
#' Economic Theory* 2: 244-263. \doi{10.1016/0022-0531(70)90039-6}
#'
#' Damgaard, C. and Weiner, J. (2000) Describing inequality in plant size or
#' fecundity. *Ecology* 81: 1139-1142.
#' \doi{10.1890/0012-9658(2000)081[1139:DIIPSO]2.0.CO;2}
#'
#' Deltas, G. (2003) The small sample bias of the Gini coefficient: results
#' and implications for empirical research. *The Review of Economics and
#' Statistics* 85: 226-234. \doi{10.1162/rest.2003.85.1.226}
#'
#' Hardy, G. H., Littlewood, J. E. and Polya, G. (1952) *Inequalities*, 2nd
#' edn. Cambridge University Press, Cambridge.
#'
#' Marshall, A. W., Olkin, I. and Arnold, B. C. (2011) *Inequalities: Theory
#' of Majorization and Its Applications*. Springer, New York.
#' \doi{10.1007/978-0-387-68276-1}
#'
#' @examples
#' rain_gini(c(5, 5, 5, 5))          # 0: perfectly even
#' rain_gini(c(0, 0, 0, 20))         # 0.75, which is (n - 1) / n: the maximum on 4 periods
#' rain_gini(c(0, 0, 0, 0, 0, 20))   # 0.833: the same series shape, more periods, higher maximum
#' rain_cv(c(0, 0, 0, 20))           # 2: the maximum on 4 periods is sqrt(4)
#'
#' # the two Lepanthes records, each read against its own n
#' concentration(lepanthes_census$recruits)   # n = 10, so Gini can reach 0.900
#' concentration(lepanthes_monthly$recruits)  # n = 22, so Gini can reach 0.955
#' @export
rain_gini <- function(x) {
  x <- sort(as.numeric(x[!is.na(x)]))
  n <- length(x); s <- sum(x)
  if (n == 0 || s == 0) return(NA_real_)
  2 * sum(seq_len(n) * x) / (n * s) - (n + 1) / n
}

#' @rdname rain_gini
#' @export
rain_cv <- function(x) {
  x <- as.numeric(x[!is.na(x)])
  if (length(x) < 2 || mean(x) == 0) return(NA_real_)
  stats::sd(x) / mean(x)
}

#' @rdname rain_gini
#' @export
concentration <- function(x) {
  x <- as.numeric(x[!is.na(x)])
  c(n = length(x), mean = mean(x), gini = rain_gini(x), cv = rain_cv(x),
    zeros = sum(x == 0),
    max_share = if (sum(x) > 0) max(x) / sum(x) else NA_real_)
}
