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
#' @return `gini()` and `cv()` return a single number (`NA` when the series
#'   sums to zero). `concentration()` returns a named numeric vector with
#'   `n`, `mean`, `gini`, `cv`, `zeros` (number of periods with a zero count)
#'   and `max_share` (the largest single period's share of the total).
#'
#' @details The Gini coefficient uses the plug-in definition
#'   \eqn{G = 2 \sum_i i x_{(i)} / (n \sum_i x_i) - (n + 1) / n}, with
#'   \eqn{x_{(i)}} the values in increasing order. It is biased downward in
#'   short series, but the bias applies equally to the observed series and to
#'   the ceiling, so a short series cannot manufacture an exceedance.
#'
#' @examples
#' gini(c(5, 5, 5, 5))          # 0: perfectly even
#' gini(c(0, 0, 0, 20))         # 0.75: everything in one period
#' concentration(lepanthes_census$recruits)
#' @export
gini <- function(x) {
  x <- sort(as.numeric(x[!is.na(x)]))
  n <- length(x); s <- sum(x)
  if (n == 0 || s == 0) return(NA_real_)
  2 * sum(seq_len(n) * x) / (n * s) - (n + 1) / n
}

#' @rdname gini
#' @export
cv <- function(x) {
  x <- as.numeric(x[!is.na(x)])
  if (length(x) < 2 || mean(x) == 0) return(NA_real_)
  stats::sd(x) / mean(x)
}

#' @rdname gini
#' @export
concentration <- function(x) {
  x <- as.numeric(x[!is.na(x)])
  c(n = length(x), mean = mean(x), gini = gini(x), cv = cv(x),
    zeros = sum(x == 0),
    max_share = if (sum(x) > 0) max(x) / sum(x) else NA_real_)
}
