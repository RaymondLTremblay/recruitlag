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
#' @section Scale and how to read it:
#'
#' **`theorem`, `ceiling` and `obs`** are concentration indices, so they
#' carry the scale of [rain_gini()] and [rain_cv()]: 0 for an even series, rising with
#' concentration, with a maximum of \eqn{(n-1)/n} and \eqn{\sqrt{n}}{sqrt(n)}
#' respectively that depends on the number of scored periods. `theorem` is
#' the bound derived from the reproductive series itself; `ceiling` is the
#' largest value actually attained by any kernel searched, on exactly the
#' periods at which recruits were scored.
#'
#' **`exceedance` is a ratio, observed over searched ceiling.** It is
#' therefore centred on 1, not on 0:
#' \itemize{
#'   \item below 1: the observed recruitment record is no more concentrated
#'     than some delay predicts. Nothing is refuted and nothing needs to be.
#'   \item at 1: the record sits exactly at the bound, which is what a pure
#'     delay of a perfectly observed series would give.
#'   \item above 1: the record is more concentrated than the expectation
#'     under any delay searched. In the *Lepanthes eltoroensis* data the
#'     ratios are 2.56 and 2.51 for the six-monthly record and 9.0 and 15.9
#'     for the monthly one.
#' }
#'
#' **An exceedance above 1 is not by itself a refutation, and no threshold
#' is supplied for it.** The bound applies to the expected series. A count
#' series is more concentrated than its own mean, so some exceedance is
#' expected even when the lag hypothesis is true, and how much depends on the
#' counts rather than on the biology. Read the exceedance as a description of
#' the gap to be explained, then take it to [ceiling_calibration()], which
#' shows what the pooled comparison can refute (in this record, nothing), and
#' to [host_lag_test()], which is where the refutation is made.
#'
#' **A high ceiling means a weak test.** The bound is reproduction's own
#' concentration, so a seasonal reproductive record gives a high ceiling and
#' little to detect. The test bites where reproduction is close to
#' continuous. The closest empirical precedent for comparing the variability
#' of the two series is Wright et al. (2005), who measured both for a
#' neotropical forest but derived no bound from them.
#'
#' @references
#' Hardy, G. H., Littlewood, J. E. and Polya, G. (1952) *Inequalities*, 2nd
#' edn. Cambridge University Press, Cambridge.
#'
#' Marshall, A. W., Olkin, I. and Arnold, B. C. (2011) *Inequalities: Theory
#' of Majorization and Its Applications*. Springer, New York.
#' \doi{10.1007/978-0-387-68276-1}
#'
#' Wright, S. J., Muller-Landau, H. C., Calderon, O. and Hernandez, A. (2005)
#' Annual and spatial variation in seedfall and seedling recruitment in a
#' neotropical forest. *Ecology* 86: 848-860. \doi{10.1890/03-0750}
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
  if (is.null(R))
    rl_abort("X was given as a numeric vector, so the recruit counts have to come from R. ",
             "Either pass R = <the recruit counts>, or pass the whole data frame as X and ",
             "name its columns with period =, reproduction = and recruits =.")
  X <- as.numeric(X); R <- as.numeric(R)
  if (is.null(t_R)) t_R <- seq.int(length(X) - length(R) + 1L, length(X))
  if (length(t_R) != length(R))
    rl_abort("t_R and R have different lengths: ", length(t_R), " periods but ", length(R),
             " recruit counts. t_R says at which periods of X the recruits in R were scored, ",
             "so the two must match one for one.")
  if (any(t_R < 1) || any(t_R > length(X)))
    rl_abort("t_R contains period(s) outside the reproductive record: ",
             rl_list(sort(t_R[t_R < 1 | t_R > length(X)])),
             ". X has ", length(X), " periods, so t_R must lie between 1 and ", length(X), ".")
  rl_check_kernels(kernels, K)
  res <- vapply(kernels, function(k) {
    o <- convolve_lag(X, k$w, t_R, lag0, missing)
    o <- o[!is.na(o)]
    c(rain_gini(o), rain_cv(o))
  }, numeric(2))
  search <- tibble::tibble(label = vapply(kernels, function(k) k$label, character(1)),
                           family = vapply(kernels, function(k) k$family, character(1)),
                           gini = res[1, ], cv = res[2, ])
  if (all(is.na(search$gini)))
    stop("every kernel gives an all-zero expected series: the reproductive record is ",
         "zero at every lagged period used. There is no ceiling to compute for this unit.")
  ib <- which.max(search$gini)
  structure(list(
    obs = c(gini = rain_gini(R), cv = rain_cv(R)),
    theorem = c(gini = rain_gini(X), cv = rain_cv(X)),
    ceiling = c(gini = max(search$gini, na.rm = TRUE), cv = max(search$cv, na.rm = TRUE)),
    exceedance = c(gini = rain_gini(R) / max(search$gini, na.rm = TRUE),
                   cv = rain_cv(R) / max(search$cv, na.rm = TRUE)),
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
