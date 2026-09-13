#' A posterior for the exceedance of the expected series, from a model of the counts
#'
#' The bound applies to the **expected** recruitment series: under any delay,
#' `E[R_t]` is a weighted average of the reproductive record and cannot be
#' more concentrated than it. [lag_ceiling()] compares the bound with the
#' realised counts, which are more concentrated than their own expectation
#' whenever counts are small, and [ceiling_calibration()] shows how much that
#' inflates the exceedance. `lag_ceiling_stan()` estimates the object the
#' bound is about. It fits a negative-binomial model to the recruit counts on
#' every unit at every scored period, with a unit effect and a period effect
#' on the log scale, and hands the posterior of the expected system totals to
#' the same ceiling arithmetic as everywhere else in the package. The
#' reproductive record enters as observed: it is the driver the hypothesis
#' names, and the ceiling is computed from it exactly as in [lag_ceiling()].
#' The result is a posterior for the concentration of the expected
#' recruitment series and for its ratio to the ceiling.
#'
#' @inheritParams lag_ceiling_boot
#' @param chains,iter_warmup,iter_sampling,seed,adapt_delta,max_treedepth
#'   Passed to CmdStan's sampler. The defaults are those of the companion
#'   paper's hierarchical recruitment model: a hierarchical model with a
#'   period effect that has to reach very low values at silent periods needs
#'   the higher target acceptance and the deeper trees.
#' @param parallel_chains Chains run at once; defaults to `chains`.
#' @param threads_per_chain Threads within each chain. Above 1 the likelihood
#'   is split across threads with Stan's `reduce_sum`, and the model is
#'   compiled with threading on (once, cached separately). Worth it for
#'   records with many units, such as a thousand patches; on a small record
#'   it only adds overhead. `chains * threads_per_chain` should not exceed the
#'   cores available.
#' @param sigma_scale Scale of the half-Student-t(3) priors on the two effect
#'   standard deviations.
#' @param phi_prior Shape and rate of the gamma prior on the negative-binomial
#'   clumping parameter. The default `c(2, 0.1)` is the working prior of the
#'   companion paper, a placeholder pending elicitation; it keeps `phi` away
#'   from both 0 and infinity, which is what lets the sampler separate
#'   clumping from the period effect.
#' @param interval `"eti"` or `"hdi"`, as in [lag_ceiling_bayesboot()].
#' @param refresh How often CmdStan reports progress; 0 is silent.
#' @param ... Further arguments to the `$sample()` method of a
#'   `CmdStanModel`.
#'
#' @return An object of class `lag_ceiling_draws` with `method = "stan"`:
#'   `point` (the posterior median of each statistic), `draws` (one row per
#'   posterior draw), `table`, `prob` (`P(exceedance > 1)`), `diagnostics`
#'   (divergent transitions, transitions that hit the maximum tree depth,
#'   largest Rhat and smallest bulk ESS over the model's parameters), `fit`
#'   (the `CmdStanMCMC` object, for anything the summary does not cover) and
#'   the inputs.
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
#' **What the model assumes.** The recruit count on a unit at a period is
#' negative binomial with a unit effect and a period effect, additive on the
#' log scale, partially pooled, and one clumping parameter. No lag is fitted;
#' the ceiling is imposed afterwards, from the observed reproductive record,
#' over the same kernels as [lag_ceiling()]. The period effects carry the
#' concentration, so with few units per period they are shrunk toward each
#' other and the estimated concentration is conservative, that is, it errs
#' toward a lower exceedance.
#'
#' **Why the reproductive record is not modelled too.** Modelling it and
#' taking the ceiling from its expected series looks symmetrical, but the
#' ceiling is then a modelled quantity that can shrink to zero when the
#' reproductive period effects do, and the exceedance, a ratio, explodes.
#' The hypothesis names the observed record as the driver, so that is what
#' the ceiling is computed from.
#'
#' **The recruits must be counts.** A rate has no negative-binomial
#' likelihood, and the function refuses it with a message.
#'
#' **Cells that were not censused are left out, not read as zero.** A unit
#' with no row, or an `NA` recruit count, at a scored period contributes
#' nothing to the likelihood for that period, whereas [lag_ceiling()] and
#' the bootstraps read it as zero. The expected series is then the total
#' the full set of units would have produced at every period, so on a
#' record where the units visited varied between periods (the *Lepanthes*
#' metapopulation, for instance) it is an effort-corrected series and can
#' sit above the bootstraps, which count what was seen.
#'
#' **Read the sampler diagnostics before the interval.** Divergent
#' transitions, an Rhat above about 1.01 or a bulk ESS below about 100 mean
#' the posterior was not explored and the interval is not to be trusted. The
#' function warns, and the printed summary carries the numbers. Transitions
#' that hit the maximum tree depth are an efficiency warning, not a validity
#' one, but many of them with a poor Rhat mean the chains sat in different
#' regions.
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
                             chains = 4, parallel_chains = chains, threads_per_chain = 1,
                             iter_warmup = 1500, iter_sampling = 2000,
                             adapt_delta = 0.95, max_treedepth = 12,
                             seed = NULL, sigma_scale = 1, phi_prior = c(2, 0.1),
                             level = 0.9, interval = c("eti", "hdi"),
                             refresh = 0, ...) {
  missing <- match.arg(missing); interval <- match.arg(interval)
  check_level(level)
  if (!is.numeric(phi_prior) || length(phi_prior) != 2 || any(phi_prior <= 0))
    rl_abort("phi_prior must be two positive numbers, the shape and rate of the gamma prior ",
             "on the clumping parameter, for example c(2, 0.1). It was given as ",
             paste(format(phi_prior), collapse = ", "), ".")
  if (!is.numeric(sigma_scale) || length(sigma_scale) != 1 || sigma_scale <= 0)
    rl_abort("sigma_scale must be a single positive number, the scale of the half-t priors ",
             "on the effect standard deviations. It was given as ",
             paste(format(sigma_scale), collapse = ", "), ".")
  need_cmdstan()
  su <- uncertainty_setup(data, unit, period, reproduction, recruits, K, lag0, kernels,
                          missing, "lag_ceiling_stan()")
  check_counts(su$R, recruits, "recruits")

  if (!is.numeric(threads_per_chain) || length(threads_per_chain) != 1 || threads_per_chain < 1)
    rl_abort("threads_per_chain must be a single whole number of 1 or more. It was given as ",
             paste(format(threads_per_chain), collapse = ", "), ".")
  threads_per_chain <- as.integer(threads_per_chain)
  cells <- which(su$observed, arr.ind = TRUE)
  N <- nrow(cells)
  Rf <- as.integer(round(su$R[cells]))
  stan_data <- list(H = nrow(su$R), J = ncol(su$R), N = N, Rf = Rf,
                    hh = as.integer(cells[, 1]), jj = as.integer(cells[, 2]),
                    mr = log(mean(Rf) + 0.5),
                    sigma_scale = sigma_scale, phi_shape = phi_prior[1], phi_rate = phi_prior[2],
                    grainsize = if (threads_per_chain > 1) max(1L, N %/% (4L * threads_per_chain)) else 0L)
  model <- stan_model_cached("ceiling_nb", threads = threads_per_chain > 1)
  args <- list(data = stan_data, chains = chains, parallel_chains = parallel_chains,
               iter_warmup = iter_warmup, iter_sampling = iter_sampling,
               adapt_delta = adapt_delta, max_treedepth = max_treedepth,
               seed = seed, refresh = refresh, ...)
  if (threads_per_chain > 1) args$threads_per_chain <- threads_per_chain
  # With refresh = 0 the run is meant to be silent, so CmdStan's chain messages
  # and the informational exceptions it prints during warmup are switched off
  # too, where this version of cmdstanr allows it. They otherwise land in a
  # rendered document by the dozen and say nothing a reader can use.
  can <- names(formals(model$sample))
  if (refresh == 0) {
    if ("show_messages" %in% can && is.null(args$show_messages)) args$show_messages <- FALSE
    if ("show_exceptions" %in% can && is.null(args$show_exceptions)) args$show_exceptions <- FALSE
  }
  fit <- do.call(model$sample, args)

  ER <- fit$draws("ER", format = "draws_matrix")
  X_obs <- colSums(su$X)
  D <- t(vapply(seq_len(nrow(ER)), function(i) {
    ceiling_stats(X_obs, as.numeric(ER[i, ]), su$li, su$W)
  }, numeric(6)))
  n_failed <- sum(!stats::complete.cases(D))
  ok <- stats::complete.cases(D)
  point <- apply(D[ok, , drop = FALSE], 2, stats::median)
  tab <- interval_table(point, D, level, interval)
  prob <- c(gini = mean(D[ok, "exceedance_gini"] > 1), cv = mean(D[ok, "exceedance_cv"] > 1))

  pars <- c("a_r", "s_ur", "s_vr", "phi_r")
  sm <- fit$summary(variables = pars)
  ds <- fit$diagnostic_summary(quiet = TRUE)
  diagnostics <- list(divergences = sum(ds$num_divergent),
                      max_treedepth = sum(ds$num_max_treedepth),
                      max_rhat = max(sm$rhat, na.rm = TRUE),
                      min_ess = min(sm$ess_bulk, na.rm = TRUE),
                      summary = sm)
  if (diagnostics$divergences > 0 || diagnostics$max_rhat > 1.01 || diagnostics$min_ess < 100)
    rl_warn("The sampler reported ", diagnostics$divergences, " divergent transition(s), a ",
            "largest Rhat of ", formatC(diagnostics$max_rhat, digits = 3, format = "f"),
            " and a smallest bulk ESS of ", round(diagnostics$min_ess),
            ". The posterior may not have been explored; raise iter_warmup and iter_sampling, ",
            "or adapt_delta (up to 0.99), before reading the interval. If Rhat stays high, the ",
            "chains are sitting in different regions: compare fit$summary() across chains.")

  structure(list(method = "stan", estimand = "expected recruitment series over the observed ceiling",
                 point = point, draws = D, table = tab, prob = prob,
                 diagnostics = diagnostics, fit = fit,
                 n_units = su$n_units, R = nrow(D), n_failed = n_failed,
                 level = level, type = interval,
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
stan_model_cached <- function(name, threads = FALSE) {
  src <- system.file("stan", paste0(name, ".stan"), package = "recruitlag")
  if (!nzchar(src))
    rl_abort("The Stan model file ", name, ".stan was not found in the installed package. ",
             "Reinstall recruitlag; the file ships in inst/stan/.")
  cache <- tools::R_user_dir("recruitlag", "cache")
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  # The threaded build is a different executable, so it gets its own copy of
  # the source and its own name in the cache.
  dst <- file.path(cache, paste0(name, if (threads) "_threads" else "", ".stan"))
  if (!file.exists(dst) || !identical(readLines(src), readLines(dst))) file.copy(src, dst, overwrite = TRUE)
  if (threads) cmdstanr::cmdstan_model(dst, quiet = TRUE, cpp_options = list(stan_threads = TRUE))
  else cmdstanr::cmdstan_model(dst, quiet = TRUE)
}

# The negative binomial is a distribution for counts: whole, non-negative. Used by every
# function that simulates or models counts; `why` names the caller's reason.
check_counts <- function(M, column, role,
                         why = paste0("lag_ceiling_stan() fits a negative-binomial model, which is a ",
                                      "distribution for counts, so the recruits must be whole numbers.")) {
  v <- M[!is.na(M)]
  bad <- v[v != round(v)]
  if (length(bad))
    rl_abort('The ', role, ' column "', column, '" has non-whole values (',
             rl_list(format(utils::head(sort(unique(bad)), 4))), "). ", why,
             " If this column is a rate, a biomass or a model estimate, pass the count it was ",
             "made from instead; lag_ceiling() and the two bootstrap functions accept any ",
             "non-negative series, because the concentration indices are scale-free.")
  invisible(TRUE)
}
