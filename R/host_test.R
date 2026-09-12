#' Expected recruitment per replicate unit under a lag profile
#'
#' The lag hypothesis, stated at the level at which it is asserted: expected
#' recruitment on each unit (host tree, plot, population) at each census is
#' that unit's own lagged reproduction, weighted by the lag profile and scaled
#' so that the unit's expected total equals its observed total. Nothing is
#' fitted; the profile fixes the shape and the observed totals fix the scale.
#'
#' @param R Either a data frame with one row per unit and census holding the
#'   columns named in `unit`, `period`, `reproduction` and `recruits`
#'   (recruits `NA` at unscored censuses; see [host_matrices()]), or a numeric
#'   matrix of recruit counts, one row per unit and one column per census at
#'   which recruits are scored, with `colnames(R)` the census indices and
#'   `rownames(R)` the unit identifiers.
#' @param X Numeric matrix of the reproductive record, one row per unit and
#'   one column per census from the first, with the same `rownames` as `R`.
#'   Ignored when `R` is a data frame.
#' @param unit,period,reproduction,recruits Column names used when `R` is a
#'   data frame.
#' @param w Lag weights of length `K + 1`.
#' @param lag0 Offset of the first bin (see [convolve_lag()]).
#' @param missing Treatment of lags before the first column of `X`:
#'   `"backfill"` reads them as the first census; `"renormalise"` drops them
#'   and renormalises the remaining weights; `"zero"` reads them as zero. The
#'   last two force a long delay to predict nothing at the first censuses and
#'   are offered only as sensitivity cases.
#'
#' @return A matrix of expected counts with the dimensions of `R`. A unit
#'   whose lagged reproduction is zero throughout, but which recruited, is
#'   given a flat expectation.
#'
#' @section Scale and how to read it:
#'
#' The returned values are expected counts on the scale of `R`, and each
#' row sums to that unit's observed recruit total. That is the whole of the
#' scaling: no rate, no offset and no fitted coefficient. Two consequences
#' follow and are worth keeping in mind when reading the output.
#'
#' First, the test built on this expectation is about **timing only**. How
#' many recruits a unit received is taken from the data and is never
#' predicted, so a unit that recruited heavily cannot make the lag hypothesis
#' look good or bad by its total, only by when those recruits arrived.
#'
#' Second, a unit that never recruited has an all-zero row and contributes
#' nothing at all, so the effective replication is the number of units with
#' at least one recruit and not `nrow(R)`. In the *Lepanthes eltoroensis*
#' host matrices that is 17 of the 23 trees, although
#' [host_lag_test()] prints all 23. `sum(rowSums(R) > 0)` gives the number
#' that actually carries the test, and it is the same set [phi_moment()]
#' estimates clumping over.
#'
#' A unit with recruits but no lagged reproduction on record is given a flat
#' row rather than being dropped, because dropping it would quietly remove
#' the units least favourable to the lag hypothesis.
#'
#' @examples
#' mu <- expected_recruits(lepanthes_hosts, unit = "host",
#'                         reproduction = "inflorescences", w = flat_kernel(4))
#' round(mu[1:3, ], 2)
#' # each unit's expected total is its observed total, by construction
#' all.equal(rowSums(mu), rowSums(host_matrices(lepanthes_hosts, unit = "host",
#'                                              reproduction = "inflorescences")$R))
#' @export
expected_recruits <- function(R, X = NULL, w, lag0 = 1L, missing = c("backfill", "renormalise", "zero"),
                              unit = "unit", period = "period", reproduction = "reproduction",
                              recruits = "recruits") {
  missing <- match.arg(missing)
  if (is.data.frame(R)) { hm <- host_matrices(R, unit, period, reproduction, recruits); R <- hm$R; X <- hm$X }
  chk <- check_host_inputs(R, X)
  cens <- chk$cens; K <- length(w) - 1L; w <- as.numeric(w) / sum(w)
  S <- matrix(0, nrow(R), ncol(R), dimnames = dimnames(R))
  for (j in seq_along(cens)) {
    src <- cens[j] - lag0 - (0:K)
    if (missing == "backfill") {
      S[, j] <- X[, pmax(src, 1L), drop = FALSE] %*% w
    } else {
      ok <- src >= 1L; ww <- w * ok
      if (missing == "renormalise" && sum(ww) > 0) ww <- ww / sum(ww)
      S[, j] <- if (any(ok)) X[, src[ok], drop = FALSE] %*% ww[ok] else 0
    }
  }
  hs <- rowSums(R); mu <- S * 0
  for (h in seq_len(nrow(R))) {
    if (hs[h] == 0) next
    mu[h, ] <- if (sum(S[h, ]) > 0) S[h, ] * hs[h] / sum(S[h, ]) else hs[h] / ncol(R)
  }
  mu
}

check_host_inputs <- function(R, X) {
  if (!is.matrix(R) || !is.matrix(X))
    rl_abort("R and X must both be matrices of units by periods, as built by ",
             "host_matrices(). They are ", class(R)[1], " and ", class(X)[1], ". ",
             "If you have a long data frame, pass it as R and name its columns.")
  if (is.null(rownames(R)) || is.null(rownames(X)))
    rl_abort("R and X need rownames: the unit identifiers, so that each unit's recruits ",
             "are matched to its own reproduction.")
  if (is.null(colnames(R)))
    rl_abort("R needs colnames: the period index of each scored census, so the lag ",
             "weights know how far back to reach.")
  cens <- as.integer(colnames(R))
  if (any(is.na(cens)))
    rl_abort("The colnames of R must be the census indices as whole numbers. ",
             "These could not be read as numbers: ",
             rl_list(colnames(R)[is.na(cens)]), ".")
  if (!setequal(rownames(R), rownames(X))) {
    only_r <- setdiff(rownames(R), rownames(X)); only_x <- setdiff(rownames(X), rownames(R))
    rl_abort("R and X describe different sets of units. ",
             if (length(only_r)) paste0(length(only_r), " in R but not in X: ", rl_list(only_r), ". ") else "",
             if (length(only_x)) paste0(length(only_x), " in X but not in R: ", rl_list(only_x), ". ") else "",
             "Both matrices must cover the same units, in any order.")
  }
  X <- X[rownames(R), , drop = FALSE]
  if (max(cens) > ncol(X))
    rl_abort("R scores recruits at census ", max(cens), ", but the reproductive record X ",
             "only reaches census ", ncol(X), ". X must cover every period from the first ",
             "to the last census at which recruits were scored.")
  list(cens = cens, X = X)
}

#' Moment estimate of negative-binomial clumping under a heterogeneous mean
#'
#' With expected counts that vary by unit and census, the clumping parameter
#' is identified by the identity `E[(R - mu)^2] = E[mu] + E[mu^2] / phi`, so
#' `phi = mean(mu^2) / (mean((R - mu)^2) - mean(mu))`, over the units that
#' recruited at least once. Estimating `phi` about a single pooled mean
#' instead attributes the between-unit spread of the means to clumping and
#' roughly halves `phi`.
#'
#' @param R Numeric matrix of recruit counts, units by scored censuses.
#' @param mu Matrix of expected counts, as returned by [expected_recruits()].
#' @param cap Upper limit returned when the counts show no overdispersion,
#'   that is, when the observed mean squared deviation does not exceed the
#'   mean and the moment estimator would be negative or undefined. The
#'   default 1e6 is an effectively Poisson value and is a sentinel, not an
#'   estimate: treat a returned `cap` as "no overdispersion detected".
#' @return A single number on the scale described below.
#'
#' @section Scale and how to read it:
#'
#' `phi` is the size (or clumping) parameter of the negative binomial, `k` in
#' the ecological literature, entering the variance as
#' \eqn{\mathrm{Var}(R) = \mu + \mu^2 / \phi}{Var(R) = mu + mu^2 / phi}. Read \eqn{1/\phi}{1/phi} rather than
#' \eqn{\phi}{phi}: it is the excess variance per unit of squared mean.
#'
#' \itemize{
#'   \item \eqn{\phi \to \infty}{phi -> infinity} (\eqn{1/\phi = 0}): no clumping, the counts
#'     are Poisson about their mean. This is the value `cap` stands in for.
#'   \item \eqn{\phi = 1}: variance \eqn{\mu + \mu^2}, strong clumping.
#'   \item \eqn{\phi < 1}: heavier still. The *Lepanthes eltoroensis* host
#'     matrices give \eqn{\phi = 0.19} under the heterogeneous mean, so at a
#'     mean of 1 recruit the variance is about 6.
#'   \item \eqn{\phi \le 0} is not returned: the estimator is capped instead.
#' }
#'
#' Smaller `phi` makes the null more permissive, because clumped counts are
#' themselves episodic. A test run at a small `phi` is therefore the
#' conservative one, which is why [host_lag_test()] accepts a vector of fixed
#' values as a stress test.
#'
#' **No threshold, and no comparison across data sets.** The package supplies
#' no value of `phi` above which a record counts as aggregated. `phi` is also
#' not a fixed constant of a species or a system: it is density dependent
#' (Taylor, Woiwod and Perry 1979), so two estimates are comparable only at
#' comparable means. Report it with the mean it was estimated at, as
#' [host_lag_test()] does.
#'
#' **Which mean it is estimated about matters.** Estimating `phi` about a
#' single pooled mean, rather than about the unit-by-census means of
#' [expected_recruits()], attributes the between-unit spread of the means to
#' clumping and roughly halves the estimate. In the companion paper the
#' single-mean value was 0.089 against the moment value 0.19 under the
#' heterogeneous mean.
#'
#' @references
#' Anscombe, F. J. (1949) The statistical analysis of insect counts based on
#' the negative binomial distribution. *Biometrics* 5: 165-173.
#' \doi{10.2307/3001918}
#'
#' Bliss, C. I. and Fisher, R. A. (1953) Fitting the negative binomial
#' distribution to biological data. *Biometrics* 9: 176-200.
#' \doi{10.2307/3001850}
#'
#' Taylor, L. R., Woiwod, I. P. and Perry, J. N. (1979) The negative binomial
#' as a dynamic ecological model for aggregation, and the density dependence
#' of k. *Journal of Animal Ecology* 48: 289-304. \doi{10.2307/4114}
#' @examples
#' m <- host_matrices(lepanthes_hosts, unit = "host", reproduction = "inflorescences")
#' mu <- expected_recruits(m$R, m$X, w = flat_kernel(4))
#' phi_moment(m$R, mu)
#' @export
phi_moment <- function(R, mu, cap = 1e6) {
  keep <- rowSums(R) > 0
  r2 <- mean((R[keep, ] - mu[keep, ])^2); mb <- mean(mu[keep, ]); m2 <- mean(mu[keep, ]^2)
  if (r2 <= mb) return(cap)
  min(m2 / (r2 - mb), cap)
}

#' The host-level test of the lag hypothesis
#'
#' Tests whether the recruitment record, scored across replicate units, could
#' have arisen under a delayed consequence of each unit's own reproduction.
#' For each lag profile the expected counts are those of
#' [expected_recruits()], the clumping is [phi_moment()] (or a value you
#' supply), and negative-binomial counts are simulated for every unit and
#' summed over units to give a system-wide census series. The test reports,
#' for each profile, the probability that a simulated series is at least as
#' extreme as the observed one on each statistic. The verdict is the maximum
#' of that probability over every profile, that is, the most favourable case
#' the lag hypothesis can make for itself.
#'
#' @inheritParams expected_recruits
#' @param K Lag horizon (see [lag_kernels()]).
#' @param kernels A `lag_kernels` object; defaults to the named families
#'   (pure delays, windows at every start, geometric decay) with no random
#'   draws. Add `n_dirichlet` draws for a fuller search, at a cost in time.
#' @param flat Whether to add the flat-rate reference, in which each unit's
#'   total is spread evenly over the censuses.
#' @param phi `"moment"` for the moment estimate under each profile, or a
#'   numeric vector of fixed values at which to evaluate every profile as
#'   well (a stress test with more clumping than the counts support).
#' @param nsim Simulated series per profile.
#'
#' @return An object of class `host_lag_test` with `obs` (the observed
#'   statistics of the census totals), `table` (one row per profile and phi
#'   setting: `label`, `family`, `phi_source`, `phi`, `gini_expected` (the Gini
#'   of the expected census totals), `p_silent`, `p_gini`, `p_cv`, `p_max`,
#'   `median_silent`), `verdict` (the maximum of each probability over
#'   profiles, with the profile attaining it), and the inputs.
#'
#' @details The statistics are: `silent`, the number of censuses at which no
#'   unit recruited; `gini` and `cv`, the concentration of the census totals;
#'   `max`, the largest census total. The Gini and CV are the pre-specified
#'   statistics of the companion paper; the silent-census count is the most
#'   interpretable and the largest census the least sensitive.
#'
#' @section Scale and how to read it:
#'
#' **The probabilities.** `p_silent`, `p_gini`, `p_cv` and `p_max` are
#' proportions in \[0, 1\]: the fraction of `nsim` simulated records, drawn
#' under the profile in that row, that are at least as extreme as the
#' observed one on that statistic. They are predictive probabilities in the
#' sense of Gelman, Meng and Stern (1996), computed under a fully specified
#' null with nothing fitted. A value near 0.5 means the observed record is an
#' ordinary draw under that profile. A value near 0 means it is not.
#'
#' **The resolution floor.** A probability cannot be resolved below
#' \eqn{1/\mathrm{nsim}}{1/nsim}. With the default `nsim = 10000` a reported 0.0000
#' means "smaller than 1e-4", not zero, and the difference between 0.0000 and
#' 0.0002 is one simulated record. Raise `nsim` before reading anything into
#' a very small value, and note that [plot.host_lag_test()] floors the axis
#' at \eqn{1/\mathrm{nsim}}{1/nsim} for the same reason.
#'
#' **The verdict is a maximum, not an average.** `verdict` takes, for each
#' statistic, the largest probability over every profile and every `phi`
#' setting, together with the profile attaining it. It is therefore the best
#' case the lag hypothesis can make for itself over the whole family
#' searched. A small verdict says that no profile in the family, not merely
#' the fitted one, reproduces the record.
#'
#' **No threshold.** The package applies no cutoff to these probabilities and
#' reports them as numbers. The 0.05 line drawn by [plot.host_lag_test()] is
#' a conventional reference mark and not a decision rule; a bright-line
#' threshold is exactly what the ASA statement on p-values advises against
#' (Wasserstein and Lazar 2016). What the companion paper reports is the
#' magnitude (a verdict of 0.0007 at the moment `phi`, 0.048 at a `phi` fixed
#' three times lower) and leaves the reader the judgement.
#'
#' **Where the test is weak.** The probabilities rise, correctly, when
#' reproduction is seasonal, because a concentrated input can produce a
#' concentrated output and the ceiling is then high. An ordinary verdict is
#' therefore not evidence for a lag; it is the absence of evidence against
#' one. See the `invented-regimes` vignette, which shows both outcomes on
#' constructed records whose mechanism is known.
#'
#' @references
#' Gelman, A., Meng, X.-L. and Stern, H. (1996) Posterior predictive
#' assessment of model fitness via realized discrepancies. *Statistica
#' Sinica* 6: 733-760.
#'
#' Wasserstein, R. L. and Lazar, N. A. (2016) The ASA statement on p-values:
#' context, process, and purpose. *The American Statistician* 70: 129-133.
#' \doi{10.1080/00031305.2016.1154108}
#'
#' @examples
#' ht <- host_lag_test(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
#'                     K = 4, nsim = 500, kernels = lag_kernels(4, bin = 6, unit = "mo"))
#' ht
#' @export
host_lag_test <- function(R, X = NULL, K, lag0 = 1L,
                          kernels = lag_kernels(K, bin = NULL),
                          missing = c("backfill", "renormalise", "zero"),
                          flat = TRUE, phi = "moment", nsim = 10000L,
                          unit = "unit", period = "period", reproduction = "reproduction",
                          recruits = "recruits") {
  missing <- match.arg(missing)
  if (is.data.frame(R)) { hm <- host_matrices(R, unit, period, reproduction, recruits); R <- hm$R; X <- hm$X }
  chk <- check_host_inputs(R, X); X <- chk$X
  rl_check_kernels(kernels, K)
  Tn <- ncol(R); tot <- colSums(R)
  obs <- c(silent = sum(tot == 0), gini = rain_gini(tot), cv = rain_cv(tot), max = max(tot))
  sim_stats <- function(mu, ph) {
    M <- matrix(0L, nsim, Tn)
    for (h in seq_len(nrow(mu))) if (any(mu[h, ] > 0))
      M <- M + matrix(stats::rnbinom(nsim * Tn, size = ph, mu = rep(mu[h, ], each = nsim)), nsim, Tn)
    z <- rowSums(M == 0)
    c(p_silent = mean(z >= obs["silent"]),
      p_gini = mean(apply(M, 1, rain_gini) >= obs["gini"], na.rm = TRUE),
      p_cv = mean(apply(M, 1, rain_cv) >= obs["cv"], na.rm = TRUE),
      p_max = mean(apply(M, 1, max) >= obs["max"]),
      median_silent = stats::median(z))
  }
  phis <- if (identical(phi, "moment")) list(moment = NA) else {
    if (!is.numeric(phi))
      rl_abort('`phi` must be "moment", or a numeric vector of fixed clumping values to ',
               "evaluate as a stress test. It is ", class(phi)[1], ".")
    c(list(moment = NA), stats::setNames(as.list(phi), paste0("fixed ", phi))) }
  one <- function(label, family, mu) {
    ph_m <- phi_moment(R, mu)
    do.call(rbind, lapply(names(phis), function(src) {
      ph <- if (src == "moment") ph_m else phis[[src]]
      s <- sim_stats(mu, ph)
      data.frame(label = label, family = family, phi_source = src, phi = ph,
                 gini_expected = rain_gini(colSums(mu)),
                 p_silent = s[["p_silent"]], p_gini = s[["p_gini"]], p_cv = s[["p_cv"]],
                 p_max = s[["p_max"]], median_silent = s[["median_silent"]],
                 stringsAsFactors = FALSE)
    }))
  }
  rows <- list()
  if (flat) {
    mu_flat <- matrix(rowSums(R) / Tn, nrow(R), Tn, dimnames = dimnames(R))
    rows[[1]] <- one("flat reference", "flat", mu_flat)
  }
  for (k in kernels)
    rows[[length(rows) + 1]] <- one(k$label, k$family, expected_recruits(R, X, k$w, lag0, missing))
  tab <- tibble::as_tibble(do.call(rbind, rows))
  verdict <- tibble::as_tibble(do.call(rbind, lapply(names(phis), function(src) {
    d <- tab[tab$phi_source == src & tab$family != "flat", ]
    do.call(rbind, lapply(c("p_silent", "p_gini", "p_cv", "p_max"), function(st) {
      i <- which.max(d[[st]])
      data.frame(phi_source = src, statistic = sub("^p_", "", st),
                 observed = unname(obs[sub("^p_", "", st)]),
                 max_p = d[[st]][i], kernel = d$label[i], stringsAsFactors = FALSE)
    }))
  })))
  structure(list(obs = obs, table = tab, verdict = verdict, R = R, X = X, K = K,
                 lag0 = lag0, missing = missing, nsim = nsim, kernels = kernels),
            class = "host_lag_test")
}

#' @export
print.host_lag_test <- function(x, digits = 4, ...) {
  cat(sprintf("Host-level test of the lag hypothesis: %d units, %d censuses, %d profiles, %d simulations each\n",
              nrow(x$R), ncol(x$R), length(x$kernels), x$nsim))
  cat(sprintf("  observed census totals: %d silent, Gini %.3f, CV %.3f, largest %d\n",
              x$obs["silent"], x$obs["gini"], x$obs["cv"], x$obs["max"]))
  fl <- x$table[x$table$family == "flat" & x$table$phi_source == "moment", ]
  if (nrow(fl)) cat(sprintf("  flat reference (phi = %.2f): P(silent) %.4f, P(Gini) %.4f, P(CV) %.4f\n",
                             fl$phi, fl$p_silent, fl$p_gini, fl$p_cv))
  cat("  most favourable profile for the lag hypothesis, by statistic:\n")
  v <- x$verdict
  for (i in seq_len(nrow(v)))
    cat(sprintf("    %-12s %-7s P = %s  (%s)\n", v$phi_source[i], v$statistic[i],
                formatC(v$max_p[i], digits = digits, format = "f"), v$kernel[i]))
  invisible(x)
}
