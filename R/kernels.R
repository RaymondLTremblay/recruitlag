#' Families of lag weights
#'
#' Builds the set of lag profiles (kernels) over which the ceiling is searched
#' and the host-level null is evaluated. Each kernel is a vector of `K + 1`
#' non-negative weights summing to one; weight `k + 1` applies to the
#' reproductive record `k` bins before the first lagged bin (see
#' [convolve_lag()] for the census convention).
#'
#' @param K Lag horizon: the number of bins beyond the first, so that a kernel
#'   has `K + 1` weights. The paper uses `K = 6` for monthly data (lags 0 to 6
#'   months) and `K = 4` for six-monthly data (five bins, 6 to 30 months
#'   before the recruit census).
#' @param families Which families to include. `"delay"`: pure delays, all
#'   weight on one bin (these attain the ceiling). `"window"`: uniform windows
#'   of every width from 2 to `K + 1`, at every start (`windows = "all"`) or
#'   anchored at the first bin only (`windows = "anchored"`). `"geometric"`:
#'   geometric decay at the rates in `rho`. `"dirichlet"`: `n_dirichlet`
#'   random draws from a symmetric Dirichlet distribution with shape `alpha`.
#' @param windows `"all"` or `"anchored"`; see `families`.
#' @param rho Decay rates for the geometric family.
#' @param n_dirichlet Number of random Dirichlet kernels.
#' @param alpha Shape of the Dirichlet distribution (0.3 gives sparse
#'   kernels; 1 is uniform on the simplex).
#' @param bin Width of one lag bin in the time unit of `unit`, used only to
#'   label kernels (for example `bin = 6, unit = "mo"` labels the six-monthly
#'   bins "0-6 mo", "6-12 mo", ...). `NULL` labels bins by index.
#' @param unit Time unit for labels.
#' @param extra A named list of additional weight vectors to append, for
#'   example a posterior-mean lag profile from a fitted model.
#'
#' @return An object of class `lag_kernels`: a list with one element per
#'   kernel, each a list with `label`, `family` and `w`. `as.matrix()` returns
#'   the weights as a matrix with one row per kernel.
#'
#' @section Scale and how to read it:
#'
#' **The weights.** Every kernel is a probability vector: `K + 1`
#' non-negative weights summing to 1. The weights say how memory is
#' distributed over the horizon, never how much recruitment there is; the
#' scale is carried elsewhere and never enters the ceiling.
#'
#' **`K`, the horizon, in bins and in time.** A kernel covers `K + 1` bins,
#' so with `bin` months per bin the horizon reaches `bin * (K + 1)` months
#' back from the first lagged bin. The paper's `K = 4` with `bin = 6` and
#' `lag0 = 1` searches 6 to 30 months before the recruit census. Nothing
#' outside the horizon is tested, so a delay longer than `bin * (K + 1)` is
#' not refuted by a small verdict: state the horizon whenever the result is
#' reported.
#'
#' **`alpha`, the Dirichlet shape.** This is the only parameter here with a
#' scale that is easy to misread. It is the shape of a symmetric Dirichlet
#' on the simplex of weights:
#' \itemize{
#'   \item `alpha < 1` (the default 0.3) pushes draws toward the vertices, so
#'     most random kernels are sparse and close to pure delays. This is the
#'     useful setting, because concentrated kernels are the ones that
#'     approach the ceiling.
#'   \item `alpha = 1` is uniform on the simplex.
#'   \item `alpha > 1` pushes draws toward the centre, giving flat profiles,
#'     which searches the least concentrated corner of the family and will
#'     understate the ceiling.
#' }
#' The construction of a Dirichlet simplex over lag weights follows Ogle et
#' al. (2015), where it is a prior; here the same object is a search grid and
#' nothing is fitted.
#'
#' **`rho`, the geometric decay rates.** Weight `k` is proportional to
#' \eqn{\rho^k}{rho^k}, so `rho` near 0 concentrates on the first bin and `rho` near
#' 1 approaches equal weights. This is the Koyck (1954) form, the shape
#' assumed whenever memory is taken to fade at a constant rate per period.
#'
#' **How many kernels you get.** The delay family contributes `K + 1`, the
#' window family every contiguous width and start (`windows = "all"`) or only
#' those anchored at the first bin, the geometric family one per `rho`, and
#' the Dirichlet family `n_dirichlet`. Print the object to see the total.
#' More kernels can only raise the searched ceiling and so can only make the
#' test more conservative.
#'
#' @references
#' Gasparrini, A. (2014) Modeling exposure-lag-response associations with
#' distributed lag non-linear models. *Statistics in Medicine* 33: 881-899.
#' \doi{10.1002/sim.5963}
#'
#' Koyck, L. M. (1954) *Distributed Lags and Investment Analysis*.
#' North-Holland, Amsterdam.
#'
#' Ogle, K., Barber, J. J., Barron-Gafford, G. A., Bentley, L. P., Young,
#' J. M., Huxman, T. E., Loik, M. E. and Tissue, D. T. (2015) Quantifying
#' ecological memory in plant and ecosystem processes. *Ecology Letters* 18:
#' 221-235. \doi{10.1111/ele.12399}
#'
#' van de Pol, M. and Cockburn, A. (2011) Identifying the critical climatic
#' time window that affects trait expression. *The American Naturalist* 177:
#' 698-707. \doi{10.1086/659101}
#'
#' @examples
#' k <- lag_kernels(4, bin = 6, unit = "mo")
#' k
#' as.matrix(k)[1:3, ]
#' @export
lag_kernels <- function(K, families = c("delay", "window", "geometric", "dirichlet"),
                        windows = c("all", "anchored"),
                        rho = c(0.2, 0.4, 0.6, 0.8, 0.95),
                        n_dirichlet = 0L, alpha = 0.3,
                        bin = NULL, unit = "", extra = NULL) {
  if (!is.numeric(K) || length(K) != 1 || is.na(K) || K < 0 || K != round(K))
    rl_abort("K must be a single whole number of 0 or more: the number of lag bins beyond ",
             "the first, so a kernel has K + 1 weights. It was given as ",
             paste(format(K), collapse = ", "), ".")
  K <- as.integer(K)
  families <- match.arg(families, several.ok = TRUE)
  windows <- match.arg(windows)
  lab <- function(from, to) {
    if (is.null(bin)) sprintf("bins %d-%d", from, to - 1L)
    else sprintf("%g-%g %s", bin * from, bin * to, unit)
  }
  out <- list()
  add <- function(label, family, w) {
    out[[length(out) + 1L]] <<- list(label = label, family = family,
                                     w = as.numeric(w) / sum(w))
  }
  if ("delay" %in% families)
    for (k in 0:K) { w <- numeric(K + 1); w[k + 1] <- 1
      add(paste("pure delay", lab(k, k + 1L)), "pure delay", w) }
  if ("window" %in% families && K >= 1) {
    starts <- if (windows == "all") 0:(K - 1) else 0L
    for (s in starts) for (width in seq(2L, K + 1L - s)) {
      w <- numeric(K + 1); w[(s + 1):(s + width)] <- 1 / width
      add(paste("uniform", lab(s, s + width)), "uniform window", w)
    }
  }
  if ("geometric" %in% families)
    for (r in rho) add(sprintf("geometric rho=%.2f", r), "geometric decay", r^(0:K))
  if ("dirichlet" %in% families && n_dirichlet > 0) {
    W <- matrix(stats::rgamma(n_dirichlet * (K + 1), shape = alpha, rate = 1), ncol = K + 1)
    for (i in seq_len(n_dirichlet))
      add(sprintf("Dirichlet(%g) draw %d", alpha, i), sprintf("Dirichlet(%g)", alpha), W[i, ])
  }
  if (!is.null(extra)) {
    if (!is.list(extra) || is.null(names(extra)) || any(names(extra) == ""))
      rl_abort("`extra` must be a NAMED list of weight vectors, for example ",
               "extra = list(`fitted profile` = c(0.4, 0.3, 0.2, 0.1)). ",
               "The names become the kernel labels.")
    for (nm in names(extra)) {
      if (length(extra[[nm]]) != K + 1)
        rl_abort('The extra kernel "', nm, '" has ', length(extra[[nm]]),
                 " weights, but K = ", K, " needs exactly ", K + 1, ".")
      if (any(extra[[nm]] < 0, na.rm = TRUE) || anyNA(extra[[nm]]))
        rl_abort('The extra kernel "', nm, '" has negative or missing weights. ',
                 "Lag weights are non-negative and are renormalised to sum to 1.")
      add(nm, "supplied", extra[[nm]])
    }
  }
  structure(out, class = "lag_kernels", K = K)
}

#' @export
print.lag_kernels <- function(x, ...) {
  fam <- table(vapply(x, function(k) k$family, character(1)))
  cat(sprintf("%d lag kernels over %d bins\n", length(x), attr(x, "K") + 1L))
  for (f in names(fam)) cat(sprintf("  %-16s %d\n", f, fam[[f]]))
  invisible(x)
}

#' @export
as.matrix.lag_kernels <- function(x, ...) {
  m <- do.call(rbind, lapply(x, function(k) k$w))
  rownames(m) <- vapply(x, function(k) k$label, character(1))
  m
}

#' The flat reference: every bin weighted equally
#'
#' A convenience for the limiting case in which reproduction contributes
#' equally at every lag, which for a near-constant reproductive record gives
#' a near-constant expected recruitment series.
#' @param K Lag horizon, as in [lag_kernels()].
#' @return A numeric vector of `K + 1` equal weights.
#'
#' @section Scale and how to read it:
#'
#' Each weight is \eqn{1/(K+1)}, so `flat_kernel(4)` is five weights of 0.2
#' and the horizon is covered evenly. It is one end of the range that
#' [lag_kernels()] searches: the flat kernel gives the least concentrated
#' expected series available, a pure delay the most concentrated, and every
#' other profile lies between.
#'
#' **It is not a null of "no lag".** The flat kernel is the lag hypothesis
#' with memory spread evenly over the whole horizon, which is a strong
#' biological claim in its own right (seed from every one of the last `K + 1`
#' bins contributing equally). The projection-matrix null, recruitment set by
#' the immediately preceding interval alone, is the opposite extreme,
#' `c(1, 0, 0, 0, 0)`.
#'
#' **Do not confuse it with the flat reference of [host_lag_test()].** The
#' `flat = TRUE` row of that function spreads each unit's recruit total
#' evenly over the censuses, which is a statement about recruitment. The flat
#' kernel spreads the weights evenly over lags, which is a statement about
#' reproduction's memory. Where reproduction is close to constant the two
#' give nearly the same expected series, and that near-coincidence is itself
#' the second diagnosis of the companion paper.
#'
#' @examples
#' flat_kernel(4)
#' sum(flat_kernel(4))
#' # the two extremes of the searched family, on the same record
#' X <- lepanthes_census$inflorescences
#' convolve_lag(X, flat_kernel(4), t_R = 3:12)      # least concentrated
#' convolve_lag(X, c(0, 0, 0, 0, 1), t_R = 3:12)    # a pure delay, most concentrated
#' @export
flat_kernel <- function(K) rep(1 / (K + 1), K + 1)
