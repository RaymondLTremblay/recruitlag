#' Calibrate the pooled ceiling comparison against realised counts
#'
#' The ceiling bounds the expected recruitment series. Observed recruits are
#' counts about that expectation, and small clumped counts look episodic
#' whatever their expectation. `ceiling_calibration()` simulates
#' negative-binomial counts about two means that both obey the ceiling, a flat
#' mean (N0) and the best-kernel mean rescaled to the observed total (N1), and
#' reports how often the simulated series is at least as concentrated as the
#' observed one. In the companion paper this shows that the pooled comparison
#' has no power on its own; the host-level test ([host_lag_test()]) is where
#' the refutation is made.
#'
#' @param x A `lag_ceiling` object.
#' @param vmr Variance-to-mean ratio of the counts, used to set the
#'   negative-binomial size parameter through `phi = mu / (vmr - 1)`.
#'   Defaults to the variance-to-mean ratio of the observed recruit series.
#' @param nsim Number of simulated series per null.
#'
#' @return An object of class `ceiling_calibration`: a data frame with one
#'   row per null (`N0 flat mean`, `N1 best-kernel mean`) giving `phi`, the
#'   median and 5 to 95% range of the simulated Gini and CV, and the
#'   probabilities `p_gini`, `p_cv` and `p_zeros` that a simulated series is
#'   at least as concentrated, or has at least as many zero periods, as the
#'   observed one. The simulated draws are kept in `attr(, "sims")` for
#'   plotting.
#'
#' @section Scale and how to read it:
#'
#' **`vmr`, the variance-to-mean ratio.** 1 is Poisson, above 1 is
#' overdispersed, below 1 is underdispersed. The negative binomial exists
#' only above 1, so a series with `vmr <= 1` raises an error rather than
#' being silently forced. The observed recruit series is used by default.
#'
#' **`phi`, the reported size parameter.** \eqn{\phi = \mu / (v - 1)}, the
#' same quantity as in [phi_moment()] and read the same way: variance
#' \eqn{\mu + \mu^2/\phi}{mu + mu^2/phi}, so small `phi` is strong clumping and
#' \eqn{\phi \to \infty}{phi -> infinity} is Poisson. It is derived here from the pooled
#' series rather than from unit-level expectations, so it is not directly
#' comparable with the `phi` of [host_lag_test()].
#'
#' **`p_gini`, `p_cv`, `p_zeros`.** Proportions in \[0, 1\]: the fraction of
#' `nsim` simulated series that are at least as concentrated as, or have at
#' least as many zero periods as, the observed one. 0.5 is the middle of the
#' simulated distribution. The Monte Carlo resolution floor is
#' \eqn{1/\mathrm{nsim}}{1/nsim}.
#'
#' **How to read a large value, which is the usual result.** A probability
#' near 0.2 to 0.4, as in the *Lepanthes eltoroensis* record, does not
#' support the lag hypothesis. It says that this comparison cannot separate
#' the two, because counts that are small and clumped look episodic whatever
#' their mean. That is a statement about the power of the pooled test and is
#' the reason [host_lag_test()] exists. Reporting only the pooled comparison
#' would understate what the data can do.
#'
#' **The two nulls.** `N0 flat mean` spreads the expected total evenly;
#' `N1 best-kernel mean` uses the shape of the most concentrated kernel found
#' by [lag_ceiling()], rescaled to the observed total. N1 is the more
#' favourable of the two to a delay, so it is the one to read first. Both
#' obey the ceiling by construction.
#'
#' **No threshold.** No cutoff is applied to these probabilities, and none is
#' taken from the literature (Wasserstein and Lazar 2016).
#'
#' @references
#' Bliss, C. I. and Fisher, R. A. (1953) Fitting the negative binomial
#' distribution to biological data. *Biometrics* 9: 176-200.
#' \doi{10.2307/3001850}
#'
#' Gelman, A., Meng, X.-L. and Stern, H. (1996) Posterior predictive
#' assessment of model fitness via realized discrepancies. *Statistica
#' Sinica* 6: 733-760.
#'
#' Wasserstein, R. L. and Lazar, N. A. (2016) The ASA statement on p-values:
#' context, process, and purpose. *The American Statistician* 70: 129-133.
#' \doi{10.1080/00031305.2016.1154108}
#'
#' @examples
#' ce <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4,
#'                   kernels = lag_kernels(4, bin = 6, unit = "mo"))
#' ceiling_calibration(ce, nsim = 2000)
#' @export
ceiling_calibration <- function(x, vmr = NULL, nsim = 20000L) {
  rl_check_class(x, "lag_ceiling", "x", "lag_ceiling()")
  R <- x$R; n <- length(R); mu <- mean(R)
  check_counts(R, "recruits", "recruit",
               "ceiling_calibration() simulates counts about the expected series, so the recruits must be whole numbers.")
  if (is.null(vmr)) vmr <- stats::var(R) / mu
  if (vmr <= 1)
    rl_abort("The recruit counts have a variance-to-mean ratio of ", signif(vmr, 3),
             ", which is at or below 1, so they are not overdispersed and the ",
             "negative-binomial null does not exist. The counts are Poisson or tighter. ",
             "Pass a value above 1 with vmr = if you want to simulate a clumped null anyway.")
  phi <- mu^2 / (mu * (vmr - 1))
  m1 <- convolve_lag(x$X, x$best$w, x$t_R, x$lag0, x$missing)
  m1 <- m1 * mu / mean(m1)
  sim <- function(mu_vec) {
    M <- matrix(stats::rnbinom(nsim * n, size = phi, mu = rep(mu_vec, each = nsim)), nsim, n)
    data.frame(gini = apply(M, 1, rain_gini), cv = apply(M, 1, rain_cv), zeros = rowSums(M == 0))
  }
  s0 <- sim(rep(mu, n)); s1 <- sim(m1)
  obs <- concentration(R)
  row <- function(s, lab) data.frame(
    null = lab, phi = phi,
    gini_med = stats::median(s$gini, na.rm = TRUE),
    gini_q05 = stats::quantile(s$gini, .05, na.rm = TRUE, names = FALSE),
    gini_q95 = stats::quantile(s$gini, .95, na.rm = TRUE, names = FALSE),
    p_gini = mean(s$gini >= obs["gini"], na.rm = TRUE),
    cv_med = stats::median(s$cv, na.rm = TRUE),
    cv_q05 = stats::quantile(s$cv, .05, na.rm = TRUE, names = FALSE),
    cv_q95 = stats::quantile(s$cv, .95, na.rm = TRUE, names = FALSE),
    p_cv = mean(s$cv >= obs["cv"], na.rm = TRUE),
    p_zeros = mean(s$zeros >= obs["zeros"]),
    stringsAsFactors = FALSE)
  out <- tibble::as_tibble(rbind(row(s0, "N0 flat mean"), row(s1, "N1 best-kernel mean")))
  structure(out, class = c("ceiling_calibration", class(out)),
            sims = list(`N0 flat mean` = s0, `N1 best-kernel mean` = s1),
            obs = obs, ceiling = x$ceiling)
}

#' @export
print.ceiling_calibration <- function(x, digits = 3, ...) {
  obs <- attr(x, "obs")
  cat(sprintf("Observed recruitment: Gini %.3f, CV %.3f, %d zero periods of %d\n",
              obs["gini"], obs["cv"], obs["zeros"], obs["n"]))
  cat("Simulated counts about a mean that obeys the ceiling:\n")
  y <- x; class(y) <- setdiff(class(y), "ceiling_calibration")
  attr(y, "sims") <- NULL; attr(y, "obs") <- NULL; attr(y, "ceiling") <- NULL
  print(y, n = Inf, width = Inf)
  cat("p_* is the probability that a simulated series is at least as concentrated as the observed one.\n")
  invisible(x)
}
