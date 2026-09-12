#' The ceiling that any delay places on the concentration of recruitment
#'
#' If recruitment is the delayed consequence of reproduction, expected
#' recruitment is a convolution of the reproductive record, and a convolution
#' is a weighted average. A weighted average is never more concentrated in
#' time than the series it averages, so reproduction's own concentration is a
#' ceiling on the concentration of expected recruitment, attained only by a
#' pure delay. `lag_ceiling()` states that bound, searches a family of lag
#' weights to confirm that none exceeds it on the periods at which recruits
#' are scored, and compares the observed recruitment series with it.
#'
#' @param X Either a numeric vector, the reproductive record at periods
#'   `1, 2, ...`, or a data frame with one row per period (or per unit and
#'   period, which is summed) holding the columns named in `period`,
#'   `reproduction` and `recruits`; recruits are `NA` at unscored periods.
#' @param R Numeric vector: recruits scored at the periods `t_R`. Ignored
#'   when `X` is a data frame.
#' @param t_R Periods of `X` at which `R` was scored; defaults to the last
#'   `length(R)` periods of `X`. Ignored when `X` is a data frame.
#' @param period,reproduction,recruits Column names used when `X` is a data
#'   frame.
#' @param K Lag horizon (see [lag_kernels()]).
#' @param lag0 Offset of the first bin (see [convolve_lag()]).
#' @param kernels A `lag_kernels` object. Defaults to every family with
#'   20,000 Dirichlet draws, as in the companion paper.
#' @param missing Treatment of pre-record lags (see [convolve_lag()]).
#'
#' @return An object of class `lag_ceiling` with elements `obs` (the
#'   concentration of `R`), `theorem` (the concentration of `X`, which is the
#'   bound), `ceiling` (the largest Gini and CV any searched kernel produced
#'   on the periods `t_R`), `exceedance` (observed over ceiling), `best` (the
#'   kernel attaining the ceiling), `search` (a data frame with one row per
#'   kernel) and the inputs.
#'
#' @details The searched ceiling can differ from the theorem bound by a small
#'   amount in either direction, because the kernels are scored on the
#'   recruitment periods only and pre-record lags are backfilled. The
#'   exceedance bounds the expected series, not the counts: whether a given
#'   exceedance refutes the delay depends on the sampling variation of the
#'   counts, which [ceiling_calibration()] evaluates for the pooled series and
#'   [host_lag_test()] evaluates across replicate units.
#'
#' @examples
#' ce <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4,
#'                   kernels = lag_kernels(4, n_dirichlet = 200, bin = 6, unit = "mo"))
#' ce
#' # the same from vectors
#' d <- lepanthes_census
#' lag_ceiling(d$inflorescences, d$recruits[3:12], t_R = 3:12, K = 4,
#'             kernels = lag_kernels(4, bin = 6, unit = "mo"))
#' @export
lag_ceiling <- function(X, R = NULL, t_R = NULL, K, lag0 = 1L,
                        kernels = lag_kernels(K, n_dirichlet = 20000L),
                        missing = c("backfill", "drop"),
                        period = "period", reproduction = "reproduction", recruits = "recruits") {
  missing <- match.arg(missing)
  if (is.data.frame(X)) {
    sf <- series_from(X, period, reproduction, recruits)
    X <- sf$X; R <- sf$R; t_R <- sf$t_R
  }
  if (is.null(R)) stop("R is required when X is not a data frame")
  X <- as.numeric(X); R <- as.numeric(R)
  if (is.null(t_R)) t_R <- seq.int(length(X) - length(R) + 1L, length(X))
  stopifnot(length(t_R) == length(R), all(t_R >= 1), all(t_R <= length(X)),
            inherits(kernels, "lag_kernels"), attr(kernels, "K") == K)
  res <- vapply(kernels, function(k) {
    o <- convolve_lag(X, k$w, t_R, lag0, missing)
    o <- o[!is.na(o)]
    c(gini(o), cv(o))
  }, numeric(2))
  search <- tibble::tibble(label = vapply(kernels, function(k) k$label, character(1)),
                           family = vapply(kernels, function(k) k$family, character(1)),
                           gini = res[1, ], cv = res[2, ])
  ib <- which.max(search$gini)
  structure(list(
    obs = c(gini = gini(R), cv = cv(R)),
    theorem = c(gini = gini(X), cv = cv(X)),
    ceiling = c(gini = max(search$gini, na.rm = TRUE), cv = max(search$cv, na.rm = TRUE)),
    exceedance = c(gini = gini(R) / max(search$gini, na.rm = TRUE),
                   cv = cv(R) / max(search$cv, na.rm = TRUE)),
    best = list(label = search$label[ib], family = search$family[ib], w = kernels[[ib]]$w),
    search = search,
    X = X, R = R, t_R = t_R, K = K, lag0 = lag0, missing = missing
  ), class = "lag_ceiling")
}

#' @export
print.lag_ceiling <- function(x, digits = 3, ...) {
  cat("Ceiling on the concentration of expected recruitment\n")
  cat(sprintf("  reproductive record: %d periods; recruits scored at %d periods; K = %d, lag0 = %d\n",
              length(x$X), length(x$R), x$K, x$lag0))
  f <- function(v) formatC(v, digits = digits, format = "f")
  cat(sprintf("  %-28s Gini %s   CV %s\n", "reproduction (theorem bound)", f(x$theorem["gini"]), f(x$theorem["cv"])))
  cat(sprintf("  %-28s Gini %s   CV %s   (%s, %d kernels)\n", "searched ceiling",
              f(x$ceiling["gini"]), f(x$ceiling["cv"]), x$best$label, nrow(x$search)))
  cat(sprintf("  %-28s Gini %s   CV %s\n", "observed recruitment", f(x$obs["gini"]), f(x$obs["cv"])))
  cat(sprintf("  %-28s Gini %s   CV %s\n", "exceedance (obs / ceiling)",
              formatC(x$exceedance["gini"], digits = 2, format = "f"),
              formatC(x$exceedance["cv"], digits = 2, format = "f")))
  cat("  The exceedance bounds the expected series, not the counts; see ceiling_calibration() and host_lag_test().\n")
  invisible(x)
}
