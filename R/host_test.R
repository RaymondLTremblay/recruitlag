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
  stopifnot(is.matrix(R), is.matrix(X), !is.null(rownames(R)), !is.null(rownames(X)),
            !is.null(colnames(R)))
  cens <- as.integer(colnames(R))
  if (any(is.na(cens))) stop("colnames(R) must be the census indices")
  if (!setequal(rownames(R), rownames(X))) stop("R and X must have the same unit identifiers as rownames")
  X <- X[rownames(R), , drop = FALSE]
  if (max(cens) > ncol(X)) stop("R scores a census beyond the last column of X")
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
#' @param cap Upper limit returned when the counts show no overdispersion.
#' @return A single number.
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
  stopifnot(inherits(kernels, "lag_kernels"), attr(kernels, "K") == K)
  Tn <- ncol(R); tot <- colSums(R)
  obs <- c(silent = sum(tot == 0), gini = gini(tot), cv = cv(tot), max = max(tot))
  sim_stats <- function(mu, ph) {
    M <- matrix(0L, nsim, Tn)
    for (h in seq_len(nrow(mu))) if (any(mu[h, ] > 0))
      M <- M + matrix(stats::rnbinom(nsim * Tn, size = ph, mu = rep(mu[h, ], each = nsim)), nsim, Tn)
    z <- rowSums(M == 0)
    c(p_silent = mean(z >= obs["silent"]),
      p_gini = mean(apply(M, 1, gini) >= obs["gini"], na.rm = TRUE),
      p_cv = mean(apply(M, 1, cv) >= obs["cv"], na.rm = TRUE),
      p_max = mean(apply(M, 1, max) >= obs["max"]),
      median_silent = stats::median(z))
  }
  phis <- if (identical(phi, "moment")) list(moment = NA) else {
    stopifnot(is.numeric(phi)); c(list(moment = NA), stats::setNames(as.list(phi), paste0("fixed ", phi))) }
  one <- function(label, family, mu) {
    ph_m <- phi_moment(R, mu)
    do.call(rbind, lapply(names(phis), function(src) {
      ph <- if (src == "moment") ph_m else phis[[src]]
      s <- sim_stats(mu, ph)
      data.frame(label = label, family = family, phi_source = src, phi = ph,
                 gini_expected = gini(colSums(mu)),
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
