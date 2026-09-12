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
  K <- as.integer(K); stopifnot(K >= 0)
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
    stopifnot(is.list(extra), !is.null(names(extra)))
    for (nm in names(extra)) {
      stopifnot(length(extra[[nm]]) == K + 1, all(extra[[nm]] >= 0))
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
#' @examples
#' flat_kernel(4)
#' sum(flat_kernel(4))
#' @export
flat_kernel <- function(K) rep(1 / (K + 1), K + 1)
