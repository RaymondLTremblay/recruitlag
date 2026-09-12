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
#' @examples
#' X <- lepanthes_census$inflorescences
#' # a pure 24 to 30 month delay, five six-month bins, recruits scored at censuses 3 to 12
#' convolve_lag(X, c(0, 0, 0, 0, 1), t_R = 3:12)
#' @export
convolve_lag <- function(X, w, t_R = NULL, lag0 = 1L, missing = c("backfill", "drop")) {
  missing <- match.arg(missing)
  X <- as.numeric(X); w <- as.numeric(w); w <- w / sum(w)
  K <- length(w) - 1L
  if (is.null(t_R)) t_R <- seq.int(lag0 + K + 1L, length(X))
  vapply(t_R, function(t) {
    src <- t - lag0 - (0:K)
    if (missing == "backfill") sum(w * X[pmax(src, 1L)])
    else if (all(src >= 1L)) sum(w * X[src]) else NA_real_
  }, numeric(1))
}
