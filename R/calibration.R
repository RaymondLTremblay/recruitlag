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
#' @examples
#' ce <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4,
#'                   kernels = lag_kernels(4, bin = 6, unit = "mo"))
#' ceiling_calibration(ce, nsim = 2000)
#' @export
ceiling_calibration <- function(x, vmr = NULL, nsim = 20000L) {
  stopifnot(inherits(x, "lag_ceiling"))
  R <- x$R; n <- length(R); mu <- mean(R)
  if (is.null(vmr)) vmr <- stats::var(R) / mu
  if (vmr <= 1) stop("vmr must exceed 1 for a negative-binomial null; the counts are not overdispersed")
  phi <- mu^2 / (mu * (vmr - 1))
  m1 <- convolve_lag(x$X, x$best$w, x$t_R, x$lag0, x$missing)
  m1 <- m1 * mu / mean(m1)
  sim <- function(mu_vec) {
    M <- matrix(stats::rnbinom(nsim * n, size = phi, mu = rep(mu_vec, each = nsim)), nsim, n)
    data.frame(gini = apply(M, 1, gini), cv = apply(M, 1, cv), zeros = rowSums(M == 0))
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
