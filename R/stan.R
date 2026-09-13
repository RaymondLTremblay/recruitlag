#' A posterior for the exceedance of the expected series, from a model of the counts
#'
#' The bound applies to the **expected** recruitment series: under any delay,
#' `E[R_t]` is a weighted average of the reproductive record and cannot be
#' more concentrated than it. [lag_ceiling()] compares the bound with the
#' realised counts, which are more concentrated than their own expectation
#' whenever counts are small, and [ceiling_calibration()] shows how much that
#' inflates the exceedance. `lag_ceiling_stan()` estimates the object the
#' bound is about. It fits a negative-binomial model to the unit-by-period
#' counts of both records, with a unit effect and a period effect on the log
#' scale, and hands the posterior of the expected system totals back to the
#' same ceiling arithmetic as everywhere else in the package. The result is
#' a posterior for the concentration of the expected recruitment series, for
#' the ceiling, and for their ratio.
#'
#' @inheritParams lag_ceiling_boot
#' @param ceiling_from `"expected"` (the default) computes the ceiling from
#'   the posterior expected reproductive series, so that both sides of the
#'   ratio are treated alike and sampling noise is removed from both.
#'   `"observed"` keeps the reproductive record as observed and models only
#'   the recruits, which reads the record as the driver itself rather than
#'   as a noisy measurement of it.
#' @param chains,iter_warmup,iter_sampling,seed Passed to CmdStan's sampler.
#' @param interval `"eti"` or `"hdi"`, as in [lag_ceiling_bayesboot()].
#' @param refresh How often CmdStan reports progress; 0 is silent.
#' @param ... Further arguments to the `$sample()` method of a
#'   `CmdStanModel`, such as `adapt_delta` or `parallel_chains`.
#'
#' @return An object of class `lag_ceiling_draws` with `method = "stan"`:
#'   `point` (the posterior median of each statistic), `draws` (one row per
#'   posterior draw), `table`, `prob` (`P(exceedance > 1)`), `diagnostics`
#'   (divergent transitions, largest Rhat and smallest bulk ESS over the
#'   model's parameters), `fit` (the `CmdStanMCMC` object, for anything the
#'   summary does not cover) and the inputs.
#'
#' @section Requirements:
#'
#' The model is compiled and run through `cmdstanr`, which is not on CRAN,
#' and CmdStan itself. Install with
#' `install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))`
#' and then `cmdstanr::install_cmdstan()`. The compiled model is cached in
#' `tools::R_user_dir("recruitlag", "cache")` after the first call. Without
#' them the function stops with that message and nothing else in the
#' package is affected.
#'
#' @section Scale and how to read it:
#'
#' **A different estimand from the two bootstraps.** They put an interval on
#' the exceedance of the realised counts; this puts one on the exceedance of
#' the expected series. The second is usually lower, because the negative-
#' binomial model attributes part of the concentration of the counts to
#' sampling variation about a smoother mean, and the difference between the
#' two is itself informative: it is the share of the raw exceedance that the
#' counts alone could have produced.
#'
#' **What the model assumes.** Each record is a count with a unit effect and
#' a period effect, additive on the log scale, partially pooled, with one
#' clumping parameter per record. No lag is fitted and the two records are
#' modelled independently; the ceiling is imposed afterwards, over the same
#' kernels as [lag_ceiling()]. The period effects are what carry the
#' concentration, so with few units per period they are shrunk toward each
#' other and the estimated concentration is conservative, that is, it errs
#' toward a lower exceedance.
#'
#' **Both records must be counts.** A rate (inflorescences per adult, say)
#' has no negative-binomial likelihood, and the function refuses it with a
#' message. Pass the count the rate was made from.
#'
#' **Read the sampler diagnostics before the interval.** Divergent
#' transitions or an Rhat above about 1.01 mean the posterior was not
#' explored and the interval is not to be trusted; raise `adapt_delta` or
#' the iterations. The printed summary carries the three numbers.
#'
#' **Equal-tailed or highest density.** As in [lag_ceiling_bayesboot()]: the
#' equal-tailed interval is invariant to a change of scale and the
#' highest-density interval is not, which matters for a ratio.
#'
#' @references
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A. and
#' Rubin, D. B. (2013) *Bayesian Data Analysis*, 3rd edn. Chapman and Hall,
#' Boca Raton.
#'
#' Stan Development Team (2024) *Stan Modeling Language Users Guide and
#' Reference Manual*. \url{https://mc-stan.org}
#'
#' @examples
#' \dontrun{
#' st <- lag_ceiling_stan(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
#'                        K = 4, kernels = lag_kernels(4, bin = 6, unit = "mo"))
#' st
#' }
#' @export
lag_ceiling_stan <- function(data, unit = "unit", period = "period",
                             reproduction = "reproduction", recruits = "recruits",
                             K, lag0 = 1L, kernels = lag_kernels(K),
                             missing = c("backfill", "drop"),
                             ceiling_from = c("expected", "observed"),
                             chains = 4, iter_warmup = 1000, iter_sampling = 1000,
                             seed = NULL, level = 0.9, interval = c("eti", "hdi"),
                             refresh = 0, ...) {
  missing <- match.arg(missing); interval <- match.arg(interval)
  ceiling_from <- match.arg(ceiling_from)
  check_level(level)
  need_cmdstan()
  su <- uncertainty_setup(data, unit, period, reproduction, recruits, K, lag0, kernels,
                          missing, "lag_ceiling_stan()")
  check_counts(su$X, reproduction, "reproduction")
  check_counts(su$R, recruits, "recruits")

  stan_data <- list(H = nrow(su$X), T = ncol(su$X), J = ncol(su$R),
                    X = unname(round(su$X)), R = unname(round(su$R)),
                    mx = log(mean(su$X) + 0.5), mr = log(mean(su$R) + 0.5))
  storage.mode(stan_data$X) <- "integer"; storage.mode(stan_data$R) <- "integer"
  model <- stan_model_cached("ceiling_nb")
  fit <- model$sample(data = stan_data, chains = chains, iter_warmup = iter_warmup,
                      iter_sampling = iter_sampling, seed = seed, refresh = refresh, ...)

  EX <- fit$draws("EX", format = "draws_matrix")
  ER <- fit$draws("ER", format = "draws_matrix")
  X_obs <- colSums(su$X)
  D <- t(vapply(seq_len(nrow(ER)), function(i) {
    Xi <- if (ceiling_from == "expected") as.numeric(EX[i, ]) else X_obs
    ceiling_stats(Xi, as.numeric(ER[i, ]), su$li, su$W)
  }, numeric(6)))
  n_failed <- sum(!stats::complete.cases(D))
  ok <- stats::complete.cases(D)
  point <- apply(D[ok, , drop = FALSE], 2, stats::median)
  tab <- interval_table(point, D, level, interval)
  prob <- c(gini = mean(D[ok, "exceedance_gini"] > 1), cv = mean(D[ok, "exceedance_cv"] > 1))

  pars <- c("a_x", "a_r", "s_ux", "s_vx", "s_ur", "s_vr", "phi_x", "phi_r")
  sm <- fit$summary(variables = pars)
  ds <- fit$diagnostic_summary(quiet = TRUE)
  diagnostics <- list(divergences = sum(ds$num_divergent),
                      max_treedepth = sum(ds$num_max_treedepth),
                      max_rhat = max(sm$rhat, na.rm = TRUE),
                      min_ess = min(sm$ess_bulk, na.rm = TRUE),
                      summary = sm)
  if (diagnostics$divergences > 0 || diagnostics$max_rhat > 1.01)
    rl_warn("The sampler reported ", diagnostics$divergences, " divergent transition(s) and a ",
            "largest Rhat of ", formatC(diagnostics$max_rhat, digits = 3, format = "f"),
            ". The posterior may not have been explored; raise adapt_delta (for example ",
            "adapt_delta = 0.95) or the iterations before reading the interval.")

  structure(list(method = "stan", estimand = paste0("expected series (ceiling from ",
                                                    ceiling_from, " reproduction)"),
                 point = point, draws = D, table = tab, prob = prob,
                 diagnostics = diagnostics, fit = fit,
                 n_units = su$n_units, R = nrow(D), n_failed = n_failed,
                 level = level, type = interval, ceiling_from = ceiling_from,
                 K = K, lag0 = lag0, missing = missing, kernels = kernels, t_R = su$t_R),
            class = "lag_ceiling_draws")
}

# cmdstanr and CmdStan, or a message saying how to get them.
need_cmdstan <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE))
    rl_abort("lag_ceiling_stan() needs the cmdstanr package, which is not on CRAN. Install it with\n",
             '  install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))\n',
             "and then CmdStan itself with cmdstanr::install_cmdstan(). The two bootstrap ",
             "functions, lag_ceiling_boot() and lag_ceiling_bayesboot(), need neither.")
  v <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)
  if (is.null(v))
    rl_abort("cmdstanr is installed but CmdStan itself was not found. Install it with ",
             "cmdstanr::install_cmdstan(), or point cmdstanr at an existing installation with ",
             "cmdstanr::set_cmdstan_path().")
  invisible(TRUE)
}

# Compile once into the user's cache directory and reuse. The .stan file is
# copied there first because the installed package directory may be
# read-only.
stan_model_cached <- function(name) {
  src <- system.file("stan", paste0(name, ".stan"), package = "recruitlag")
  if (!nzchar(src))
    rl_abort("The Stan model file ", name, ".stan was not found in the installed package. ",
             "Reinstall recruitlag; the file ships in inst/stan/.")
  cache <- tools::R_user_dir("recruitlag", "cache")
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  dst <- file.path(cache, paste0(name, ".stan"))
  if (!file.exists(dst) || !identical(readLines(src), readLines(dst))) file.copy(src, dst, overwrite = TRUE)
  cmdstanr::cmdstan_model(dst, quiet = TRUE)
}

# The negative binomial is a distribution for counts: whole, non-negative.
check_counts <- function(M, column, role) {
  v <- M[!is.na(M)]
  bad <- v[v != round(v)]
  if (length(bad))
    rl_abort('The ', role, ' column "', column, '" has non-whole values (',
             rl_list(format(utils::head(sort(unique(bad)), 4))), "). lag_ceiling_stan() fits a ",
             "negative-binomial model, which is a distribution for counts, so both records must ",
             "be whole numbers. If this column is a rate such as inflorescences per adult, pass the ",
             "count it was made from instead. The two bootstrap functions accept rates.")
  invisible(TRUE)
}
