#' Apply a lag profile to a reproductive record
#'
#' Computes the expected recruitment series implied by a set of lag weights:
#' at each period at which recruits are scored, the weighted sum of the
#' reproductive record over the preceding bins.
#'
#' @param X Numeric vector: the reproductive record at periods `1, 2, ...`.
#' @param w Numeric vector of `K + 1` non-negative lag weights (they are
#'   renormalised to sum to one).
#' @param t_R Integer vector: the periods (indices into `X`) at which recruits
#'   are scored. Defaults to the last `length(X) - lag0 - K` periods.
#' @param lag0 Offset of the first bin. With `lag0 = 1` (the default, and the
#'   convention of a projection matrix) weight `k + 1` is applied to
#'   `X[t - 1 - k]`, so the first bin is the period preceding the recruit
#'   census. With `lag0 = 0` the first bin is the recruit census itself.
#' @param missing How to treat lags that fall before the first period of `X`.
#'   `"backfill"` (the default) reads them as `X[1]`, which asserts that
#'   reproduction before the record began was at the level first observed.
#'   `"drop"` returns `NA` for any period whose window is not fully within
#'   the record. Zero-filling is deliberately not offered: it asserts that no
#'   seed fell before the record began, which turns a long pure delay into a
#'   prediction of no recruits at the first censuses and is an artefact.
#'
#' @return A numeric vector of the same length as `t_R`.
#'
#' @section Scale and how to read it:
#'
#' The output is on the scale of `X` and in its units, because the weights
#' are renormalised to sum to 1 and the result is a weighted average: every
#' value lies between the smallest and the largest `X` in its own window.
#' Nothing is multiplied by a rate or a recruitment probability, so the
#' output is not a predicted number of recruits. It is the shape of the
#' expected recruitment series, up to an unknown scale, which is all the
#' ceiling argument needs and all that [lag_ceiling()] uses.
#'
#' Two consequences are worth stating because they are where the intuition
#' usually slips. Averaging cannot sharpen: the output can never be more
#' concentrated in time than `X` itself, whatever the weights, and it equals
#' `X`'s concentration only for a pure delay, which permutes the series
#' rather than mixing it. And averaging a nearly constant `X` gives a nearly
#' constant output for every set of weights, which is why a flat reproductive
#' record cannot identify a lag profile. Slutzky (1937) is the classical
#' statement of what a moving average does to a series.
#'
#' @references
#' Koyck, L. M. (1954) *Distributed Lags and Investment Analysis*.
#' North-Holland, Amsterdam.
#'
#' Slutzky, E. (1937) The summation of random causes as the source of cyclic
#' processes. *Econometrica* 5: 105-146. \doi{10.2307/1907241}
#'
#' @examples
#' X <- lepanthes_census$inflorescences
#' # a pure 24 to 30 month delay, five six-month bins, recruits scored at censuses 3 to 12
#' convolve_lag(X, c(0, 0, 0, 0, 1), t_R = 3:12)
#' @export
convolve_lag <- function(X, w, t_R = NULL, lag0 = 1L, missing = c("backfill", "drop")) {
  missing <- match.arg(missing)
  if (!is.numeric(X) || !length(X))
    rl_abort("X must be a non-empty numeric vector: the reproductive record at periods ",
             "1, 2, ... It is ", class(X)[1], " of length ", length(X), ".")
  if (!is.numeric(w) || !length(w))
    rl_abort("w must be a numeric vector of lag weights, one per bin. It is ",
             class(w)[1], " of length ", length(w), ".")
  if (anyNA(w) || any(w < 0))
    rl_abort("The lag weights w must all be non-negative and present. ",
             "Given: ", rl_list(format(w)), ".")
  if (sum(w) == 0)
    rl_abort("The lag weights w sum to zero, so there is no profile to apply. ",
             "At least one bin must carry weight.")
  X <- as.numeric(X); w <- as.numeric(w); w <- w / sum(w)
  K <- length(w) - 1L
  if (is.null(t_R)) t_R <- seq.int(lag0 + K + 1L, length(X))
  vapply(t_R, function(t) {
    src <- t - lag0 - (0:K)
    if (missing == "backfill") sum(w * X[pmax(src, 1L)])
    else if (all(src >= 1L)) sum(w * X[src]) else NA_real_
  }, numeric(1))
}
