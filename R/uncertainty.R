# Uncertainty on the exceedance. Three routes to the same question, so that
# a user can take the one whose logic they accept: a cluster bootstrap over
# units (a confidence interval and a one-sided bootstrap p-value), a Bayesian
# bootstrap over units (a posterior under a noninformative prior on the
# distribution of units), and, in stan.R, a negative-binomial model of the
# counts (a posterior on the concentration of the EXPECTED series, which is
# a different estimand). The shared arithmetic is here; nothing in it is
# fitted.

# ---- internal: the lag index of the scored periods ---------------------------

# Where each scored period reads its lagged reproduction. Computed once per
# call, because it depends on t_R, K and lag0 and not on the series. Under
# "backfill" every row is kept and pre-record lags read period 1; under
# "drop" rows whose window is not fully within the record are excluded from
# the ceiling, exactly as convolve_lag() and lag_ceiling() do.
lag_index <- function(t_R, K, lag0, missing) {
  src <- outer(as.integer(t_R), 0:K, function(t, k) t - lag0 - k)
  keep <- if (missing == "backfill") rep(TRUE, nrow(src)) else apply(src >= 1L, 1, all)
  list(index = pmax(src, 1L)[keep, , drop = FALSE], keep = keep)
}

# Column-wise Gini and CV of a periods-by-kernels matrix, vectorised so that
# thousands of replicates stay cheap.
gini_cols <- function(M) {
  n <- nrow(M); s <- colSums(M)
  Ms <- apply(M, 2, sort)
  if (!is.matrix(Ms)) Ms <- matrix(Ms, nrow = n)
  g <- 2 * colSums(seq_len(n) * Ms) / (n * s) - (n + 1) / n
  g[s == 0] <- NA_real_
  g
}
cv_cols <- function(M) {
  n <- nrow(M); m <- colMeans(M)
  if (n < 2) return(rep(NA_real_, ncol(M)))
  s <- sqrt(colSums((M - rep(m, each = n))^2) / (n - 1))
  out <- s / m
  out[m == 0] <- NA_real_
  out
}

# The six statistics for one pooled pair of series. `li` is lag_index(),
# `W` the kernel matrix (one row per kernel). Mirrors lag_ceiling(): the
# observed concentration is that of the whole recruit series; the ceiling
# is the maximum over kernels on the periods kept.
ceiling_stats <- function(X, R, li, W) {
  S <- matrix(X[li$index], nrow = nrow(li$index))
  out <- S %*% t(W)
  cg <- suppressWarnings(max(gini_cols(out), na.rm = TRUE))
  cc <- suppressWarnings(max(cv_cols(out), na.rm = TRUE))
  if (!is.finite(cg)) cg <- NA_real_
  if (!is.finite(cc)) cc <- NA_real_
  og <- rain_gini(R); oc <- rain_cv(R)
  c(observed_gini = og, ceiling_gini = cg, exceedance_gini = og / cg,
    observed_cv = oc, ceiling_cv = cc, exceedance_cv = oc / cc)
}

stat_labels <- c(observed_gini = "observed Gini", ceiling_gini = "ceiling Gini",
                 exceedance_gini = "exceedance Gini", observed_cv = "observed CV",
                 ceiling_cv = "ceiling CV", exceedance_cv = "exceedance CV")

# The unit matrices, the kernels and the lag index, shared by the three
# routes. Refuses without a unit column, because every route resamples or
# models units, and a pooled series has no replicate to work with.
uncertainty_setup <- function(data, unit, period, reproduction, recruits, K, lag0,
                              kernels, missing, what) {
  data <- as.data.frame(data)
  if (!unit %in% names(data))
    rl_abort(what, ' needs the replicate unit. unit = "', unit, '" is not a column of ',
             "the data (it has: ", rl_list(names(data), 12), ").", rl_near(unit, names(data)),
             "\nEvery interval here comes from resampling or modelling units (host trees, ",
             "plots, populations), so a series that is already pooled over units cannot be ",
             "given one. For the pooled ceiling alone, use lag_ceiling().")
  hm <- host_matrices(data, unit, period, reproduction, recruits)
  rl_check_kernels(kernels, K)
  if (nrow(hm$R) < 2)
    rl_abort(what, " found only one unit (", rownames(hm$R)[1], "). Resampling one unit ",
             "gives the same series every time, so there is no interval to compute.")
  t_R <- as.integer(colnames(hm$R))
  list(R = hm$R, X = hm$X, t_R = t_R, W = as.matrix(kernels),
       li = lag_index(t_R, K, lag0, missing), n_units = nrow(hm$R))
}

# ---- internal: intervals ----------------------------------------------------

# Bias-corrected and accelerated interval (Efron 1987) for one statistic,
# from the bootstrap replicates and the leave-one-unit-out estimates. Falls
# back to the percentile interval, and says so, when the bias correction is
# undefined (every replicate on one side of the estimate).
bca_interval <- function(theta, boot, jack, level) {
  boot <- boot[is.finite(boot)]
  a <- (1 - level) / 2
  if (!length(boot)) return(list(lower = NA_real_, upper = NA_real_, fallback = "no finite replicates"))
  p0 <- mean(boot < theta)
  if (p0 <= 0 || p0 >= 1)
    return(list(lower = unname(stats::quantile(boot, a)), upper = unname(stats::quantile(boot, 1 - a)),
                fallback = "percentile: every replicate fell on one side of the estimate"))
  z0 <- stats::qnorm(p0)
  jack <- jack[is.finite(jack)]
  jm <- mean(jack)
  den <- 6 * sum((jm - jack)^2)^1.5
  acc <- if (den > 0) sum((jm - jack)^3) / den else 0
  adj <- function(alpha) {
    z <- stats::qnorm(alpha)
    stats::pnorm(z0 + (z0 + z) / (1 - acc * (z0 + z)))
  }
  list(lower = unname(stats::quantile(boot, adj(a))),
       upper = unname(stats::quantile(boot, adj(1 - a))), fallback = NA_character_)
}

percentile_interval <- function(boot, level) {
  boot <- boot[is.finite(boot)]
  a <- (1 - level) / 2
  if (!length(boot)) return(list(lower = NA_real_, upper = NA_real_))
  list(lower = unname(stats::quantile(boot, a)), upper = unname(stats::quantile(boot, 1 - a)))
}

# Shortest interval holding `level` of the draws. Not invariant to a change
# of scale: the HDI of a ratio and the HDI of its log, back-transformed, are
# different intervals. Documented at the two functions that offer it.
hdi_interval <- function(draws, level) {
  d <- sort(draws[is.finite(draws)])
  n <- length(d)
  if (n < 2) return(list(lower = NA_real_, upper = NA_real_))
  k <- max(1L, floor(level * n))
  if (k >= n) return(list(lower = d[1], upper = d[n]))
  width <- d[(k + 1):n] - d[1:(n - k)]
  i <- which.min(width)
  list(lower = d[i], upper = d[i + k])
}

# The interval table from a replicate matrix and the point estimates.
interval_table <- function(point, draws, level, type, jack = NULL) {
  rows <- lapply(names(point), function(s) {
    b <- draws[, s]
    iv <- switch(type,
      bca = bca_interval(point[[s]], b, jack[, s], level),
      percentile = percentile_interval(b, level),
      eti = percentile_interval(b, level),
      hdi = hdi_interval(b, level))
    note <- if (!is.null(iv$fallback) && !is.na(iv$fallback)) iv$fallback else ""
    # An interval that misses its own estimate, or has no width, is the BCa
    # correction failing on too few distinct resamples. Say so rather than
    # print it as if it were an interval.
    if (type == "bca" && all(is.finite(c(iv$lower, iv$upper))) &&
        (iv$lower > point[[s]] || iv$upper < point[[s]] || iv$upper == iv$lower))
      note <- "BCa failed: too few distinct resamples; use type = 'percentile' or lag_ceiling_bayesboot()"
    data.frame(quantity = stat_labels[[s]], estimate = unname(point[[s]]),
               lower = iv$lower, upper = iv$upper, note = note,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ---- the cluster bootstrap --------------------------------------------------

#' A confidence interval for the exceedance, by resampling units
#'
#' The ceiling is arithmetic and carries no uncertainty, but the two
#' concentration indices it is compared with are estimated from a short
#' record, and so is their ratio, the exceedance. `lag_ceiling_boot()` puts
#' a confidence interval on them by resampling the replicate units (host
#' trees, plots, populations) with replacement, rebuilding **both** the
#' reproductive series and the recruit series from the resampled units, and
#' recomputing the ceiling and the exceedance in every replicate. Numerator
#' and denominator therefore move together, as they should, since the same
#' units generate both series.
#'
#' @param data A data frame with one row per unit and period.
#' @param unit,period,reproduction,recruits Column names.
#' @param K Lag horizon (see [lag_kernels()]).
#' @param lag0 Offset of the first bin (see [convolve_lag()]).
#' @param kernels A `lag_kernels` object, used in every replicate. The
#'   default is the named families with no random draws: the pure delays
#'   attain the bound, so random kernels can only add time.
#' @param missing Treatment of pre-record lags (see [convolve_lag()]).
#' @param R Number of bootstrap replicates.
#' @param level Coverage of the interval. The default 0.9 matches the 90%
#'   credible intervals of the companion paper.
#' @param type `"bca"` (the default) for the bias-corrected and accelerated
#'   interval of Efron (1987), with the acceleration from the leave-one-unit-
#'   out estimates; `"percentile"` for the plain percentile interval.
#'
#' @return An object of class `lag_ceiling_draws` with `method = "bootstrap"`:
#'   `point` (the six statistics on the full data, identical to
#'   [lag_ceiling()] on the same kernels), `draws` (one row per replicate),
#'   `table` (estimate and interval for each statistic), `p_boot` (the
#'   one-sided bootstrap p-value for exceedance at most 1, on each index),
#'   `n_units`, `R`, `n_failed` (replicates in which a statistic could not be
#'   computed, usually because no resampled unit recruited) and the inputs.
#'
#' @section Scale and how to read it:
#'
#' **The interval is a confidence interval** and is read in the usual way:
#' the procedure covers the true value in `level` of repeated samples of
#' units from the same population of units. It is not a probability that the
#' exceedance lies in the interval. The Bayesian reading is available from
#' [lag_ceiling_bayesboot()], which is built to give it.
#'
#' **`p_boot` is a one-sided bootstrap p-value**, the fraction of replicates
#' in which the exceedance was at most 1: the achieved significance level
#' against the hypothesis that the record is no more concentrated than a
#' delay predicts. It is not the probability that the exceedance is at most
#' 1, although it is easily mistaken for one, and it cannot be resolved
#' below `1/R`. A value printed as `< 1/R` means no replicate reached 1.
#'
#' **Why units and not periods.** The Gini and the CV are permutation
#' invariant, so resampling periods would not disturb them by reordering,
#' but periods within a unit are not independent and the units are the
#' design's replicates. Resampling units is the cluster bootstrap of Davison
#' and Hinkley (1997, ch. 3). With few units (the *Lepanthes eltoroensis*
#' record has 23) the BCa interval is preferred to the percentile one, since
#' it corrects for the skew a ratio carries and for the bias between the
#' bootstrap distribution and the estimate.
#'
#' **What the interval does not say.** A wide interval does not weaken the
#' ceiling, which is a bound and not an estimate. It says how firmly the
#' record itself sits above the bound. Whether an exceedance above 1 refutes
#' the delay is still the question for [host_lag_test()], because the bound
#' applies to the expected series and not to the counts.
#'
#' @references
#' Davison, A. C. and Hinkley, D. V. (1997) *Bootstrap Methods and their
#' Application*. Cambridge University Press, Cambridge.
#'
#' DiCiccio, T. J. and Efron, B. (1996) Bootstrap confidence intervals.
#' *Statistical Science* 11: 189-228. \doi{10.1214/ss/1032280214}
#'
#' Efron, B. (1987) Better bootstrap confidence intervals. *Journal of the
#' American Statistical Association* 82: 171-185.
#' \doi{10.1080/01621459.1987.10478410}
#'
#' @examples
#' set.seed(1)
#' b <- lag_ceiling_boot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
#'                       K = 4, R = 200, kernels = lag_kernels(4, bin = 6, unit = "mo"))
#' b
#' as.data.frame(b)
#' @export
lag_ceiling_boot <- function(data, unit = "unit", period = "period",
                             reproduction = "reproduction", recruits = "recruits",
                             K, lag0 = 1L, kernels = lag_kernels(K),
                             missing = c("backfill", "drop"),
                             R = 2000L, level = 0.9, type = c("bca", "percentile")) {
  missing <- match.arg(missing); type <- match.arg(type)
  check_level(level)
  if (!is.numeric(R) || length(R) != 1 || R < 2 || R != round(R))
    rl_abort("R must be a single whole number of replicates, at least 2. It was given as ",
             paste(format(R), collapse = ", "), ".")
  R <- as.integer(R)
  su <- uncertainty_setup(data, unit, period, reproduction, recruits, K, lag0, kernels,
                          missing, "lag_ceiling_boot()")
  H <- su$n_units
  if (H < 10 && type == "bca")
    rl_warn("Only ", H, " units to resample. The bootstrap distribution then has few distinct ",
            "values and the BCa correction can collapse or miss the estimate; any such row is ",
            "flagged in the table. Compare with type = \"percentile\" and with ",
            "lag_ceiling_bayesboot(), whose Dirichlet weights are continuous.")
  stat <- function(idx) ceiling_stats(colSums(su$X[idx, , drop = FALSE]),
                                      colSums(su$R[idx, , drop = FALSE]), su$li, su$W)
  point <- stat(seq_len(H))
  draws <- t(vapply(seq_len(R), function(b) stat(sample.int(H, H, replace = TRUE)), numeric(6)))
  jack <- t(vapply(seq_len(H), function(h) stat(setdiff(seq_len(H), h)), numeric(6)))
  n_failed <- sum(!stats::complete.cases(draws))
  if (n_failed > 0.05 * R)
    rl_warn(n_failed, " of ", R, " bootstrap replicates (", round(100 * n_failed / R),
            "%) could not compute every statistic, usually because the resampled units ",
            "held no recruits or no lagged reproduction. The intervals are computed from ",
            "the rest. With this many failures the units carrying the record are few, ",
            "and the interval should be read with that in mind.")
  tab <- interval_table(point, draws, level, type, jack)
  ok <- stats::complete.cases(draws)
  p_boot <- c(gini = mean(draws[ok, "exceedance_gini"] <= 1),
              cv = mean(draws[ok, "exceedance_cv"] <= 1))
  structure(list(method = "bootstrap", estimand = "realised counts",
                 point = point, draws = draws, jack = jack, table = tab, p_boot = p_boot,
                 n_units = H, R = R, n_failed = n_failed, level = level, type = type,
                 K = K, lag0 = lag0, missing = missing, kernels = kernels, t_R = su$t_R),
            class = "lag_ceiling_draws")
}

# ---- the Bayesian bootstrap -------------------------------------------------

#' A posterior for the exceedance, by the Bayesian bootstrap over units
#'
#' The same quantity as [lag_ceiling_boot()], the exceedance of the realised
#' counts over the ceiling, but constructed as a posterior. Rubin's (1981)
#' Bayesian bootstrap places a flat Dirichlet prior on the unknown
#' distribution of units and draws, at each iteration, a set of weights
#' `Dirichlet(1, ..., 1)` over the units; both series are rebuilt as the
#' weighted sums and the ceiling and exceedance are recomputed. The result is
#' a posterior distribution for each statistic, from which a credible
#' interval and `P(exceedance > 1)` are read directly.
#'
#' @inheritParams lag_ceiling_boot
#' @param draws Number of posterior draws.
#' @param interval `"eti"` (the default) for the equal-tailed interval, the
#'   central `level` of the posterior; `"hdi"` for the highest-density
#'   interval, the shortest interval holding `level` of the posterior.
#'
#' @return An object of class `lag_ceiling_draws` with
#'   `method = "bayesian bootstrap"`: `point` (the statistics on the full
#'   data), `draws` (one row per posterior draw), `table` (estimate and
#'   credible interval for each statistic), `prob` (`P(exceedance > 1)` on
#'   each index), `n_units` and the inputs.
#'
#' @section Scale and how to read it:
#'
#' **The interval is a credible interval** and `prob` is a posterior
#' probability: under the stated prior, that is the probability that the
#' exceedance exceeds 1. This is the reading the ordinary bootstrap's
#' `p_boot` is often given and does not have. The prior is the flat
#' Dirichlet on the distribution of units, which is noninformative in
#' Rubin's sense and is the same prior for every data set; nothing else is
#' assumed, and no likelihood for the counts is written down.
#'
#' **Numerically the two bootstraps usually agree closely**, because the
#' Bayesian bootstrap is a smoothed version of the ordinary one. Where they
#' differ is in what the sentence built on them may say.
#'
#' **Equal-tailed or highest density.** The equal-tailed interval is
#' invariant to a change of scale: the 5th percentile of the log exceedance
#' is the log of the 5th percentile of the exceedance. The highest-density
#' interval is not: the shortest interval on the ratio and the shortest
#' interval on its log, back-transformed, are different intervals, and for a
#' skewed ratio the difference is visible. `interval = "hdi"` is offered so
#' that a user who prefers it can have it, and the help page says this so
#' that the choice is made knowingly.
#'
#' **What it estimates.** The exceedance of the *realised* counts, as in
#' [lag_ceiling()] and [lag_ceiling_boot()]. A count series is more
#' concentrated than its own expectation, so the exceedance is inflated by
#' sampling variation in a way that no interval on it removes;
#' [lag_ceiling_stan()] estimates the concentration of the expected series
#' instead, which is the object the bound actually applies to.
#'
#' @references
#' Rubin, D. B. (1981) The Bayesian bootstrap. *The Annals of Statistics* 9:
#' 130-134. \doi{10.1214/aos/1176345338}
#'
#' @examples
#' set.seed(1)
#' bb <- lag_ceiling_bayesboot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
#'                             K = 4, draws = 200, kernels = lag_kernels(4, bin = 6, unit = "mo"))
#' bb
#' @export
lag_ceiling_bayesboot <- function(data, unit = "unit", period = "period",
                                  reproduction = "reproduction", recruits = "recruits",
                                  K, lag0 = 1L, kernels = lag_kernels(K),
                                  missing = c("backfill", "drop"),
                                  draws = 4000L, level = 0.9, interval = c("eti", "hdi")) {
  missing <- match.arg(missing); interval <- match.arg(interval)
  check_level(level)
  if (!is.numeric(draws) || length(draws) != 1 || draws < 2 || draws != round(draws))
    rl_abort("draws must be a single whole number of posterior draws, at least 2. It was ",
             "given as ", paste(format(draws), collapse = ", "), ".")
  draws <- as.integer(draws)
  su <- uncertainty_setup(data, unit, period, reproduction, recruits, K, lag0, kernels,
                          missing, "lag_ceiling_bayesboot()")
  H <- su$n_units
  point <- ceiling_stats(colSums(su$X), colSums(su$R), su$li, su$W)
  D <- t(vapply(seq_len(draws), function(b) {
    w <- stats::rexp(H); w <- w / sum(w)
    ceiling_stats(colSums(w * su$X), colSums(w * su$R), su$li, su$W)
  }, numeric(6)))
  n_failed <- sum(!stats::complete.cases(D))
  tab <- interval_table(point, D, level, interval)
  ok <- stats::complete.cases(D)
  prob <- c(gini = mean(D[ok, "exceedance_gini"] > 1), cv = mean(D[ok, "exceedance_cv"] > 1))
  structure(list(method = "bayesian bootstrap", estimand = "realised counts",
                 point = point, draws = D, table = tab, prob = prob,
                 n_units = H, R = draws, n_failed = n_failed, level = level, type = interval,
                 K = K, lag0 = lag0, missing = missing, kernels = kernels, t_R = su$t_R),
            class = "lag_ceiling_draws")
}

check_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1 || is.na(level) || level <= 0 || level >= 1)
    rl_abort("level must be a single number strictly between 0 and 1, the coverage of the ",
             "interval (for example 0.9). It was given as ", paste(format(level), collapse = ", "), ".")
  invisible(TRUE)
}

# ---- methods shared by the three routes -------------------------------------

#' @export
print.lag_ceiling_draws <- function(x, digits = 3, ...) {
  what <- switch(x$method,
    bootstrap = sprintf("Cluster bootstrap of the ceiling: %d units resampled %d times, %d%% %s intervals",
                        x$n_units, x$R, round(100 * x$level), interval_label(x$type)),
    `bayesian bootstrap` = sprintf("Bayesian bootstrap of the ceiling: %d units, %d posterior draws, %d%% %s credible intervals",
                                   x$n_units, x$R, round(100 * x$level), interval_label(x$type)),
    stan = sprintf("Negative-binomial model of the counts: %d units, %d posterior draws, %d%% %s credible intervals",
                   x$n_units, x$R, round(100 * x$level), interval_label(x$type)))
  cat(what, "\n", sep = "")
  cat(sprintf("  estimand: %s; K = %d, lag0 = %d, %d kernels\n",
              x$estimand, x$K, x$lag0, length(x$kernels)))
  tab <- x$table
  f <- function(v, d = digits) ifelse(is.na(v), "   NA", formatC(v, digits = d, format = "f"))
  for (i in seq_len(nrow(tab))) {
    d <- if (grepl("exceedance", tab$quantity[i])) 2 else digits
    cat(sprintf("  %-17s %s   [%s, %s]%s\n", tab$quantity[i], f(tab$estimate[i], d),
                f(tab$lower[i], d), f(tab$upper[i], d),
                if (nzchar(tab$note[i])) paste0("   (", tab$note[i], ")") else ""))
  }
  if (x$method == "bootstrap") {
    fp <- function(p) if (p == 0) sprintf("< %s", format(1 / x$R, digits = 2)) else format(p, digits = 2)
    cat(sprintf("  one-sided bootstrap p-value for exceedance <= 1: Gini %s, CV %s\n",
                fp(x$p_boot[["gini"]]), fp(x$p_boot[["cv"]])))
    cat("  A confidence interval, read in the usual way; p_boot is not a posterior probability.\n")
  } else {
    fp <- function(p) if (p == 1) sprintf("> %s", format(1 - 1 / x$R, digits = 4)) else format(p, digits = 3)
    cat(sprintf("  P(exceedance > 1): Gini %s, CV %s\n", fp(x$prob[["gini"]]), fp(x$prob[["cv"]])))
  }
  if (x$method == "stan" && !is.null(x$diagnostics)) {
    d <- x$diagnostics
    cat(sprintf("  sampler: %d divergent transitions, max Rhat %.3f, min bulk ESS %.0f\n",
                d$divergences, d$max_rhat, d$min_ess))
  }
  if (x$n_failed > 0)
    cat(sprintf("  %d of %d replicates could not compute every statistic and were left out\n",
                x$n_failed, x$R))
  invisible(x)
}

#' @export
as.data.frame.lag_ceiling_draws <- function(x, ...) {
  tab <- x$table
  tab$method <- x$method
  tab$level <- x$level
  tab$interval <- x$type
  tab[, c("method", "quantity", "estimate", "lower", "upper", "level", "interval", "note")]
}

interval_label <- function(type) {
  c(bca = "BCa", percentile = "percentile", eti = "equal-tailed", hdi = "highest-density")[[type]]
}
