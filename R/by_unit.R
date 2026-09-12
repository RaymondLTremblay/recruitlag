#' The ceiling for every unit in a long data frame
#'
#' Runs [lag_ceiling()] once per replicate unit (population, plot, host tree)
#' and returns the results together, so that a figure of one strip per unit,
#' or a table of one row per unit, is a single call. This is the comparison
#' the ceiling is usually wanted for: not whether one record exceeds its
#' bound, but how far above the bound each unit sits and how much that varies
#' within a study.
#'
#' @param data A data frame with one row per unit and period.
#' @param unit,period,reproduction,recruits Column names.
#' @param K Lag horizon (see [lag_kernels()]).
#' @param lag0 Offset of the first bin (see [convolve_lag()]).
#' @param kernels A `lag_kernels` object, used for every unit so that the
#'   units are compared on the same family.
#' @param missing Treatment of pre-record lags (see [convolve_lag()]).
#' @param min_periods Units with fewer scored periods than this are skipped.
#'   The default `K + 3` is the shortest record on which a horizon of `K + 1`
#'   bins leaves anything to compare.
#' @param min_recruits Units with fewer recruits than this are skipped, since
#'   a unit with one recruit has concentration 1 by arithmetic rather than by
#'   biology.
#'
#' @details Each unit's periods are re-indexed to start at 1 before its
#'   ceiling is computed, so that a unit watched from 2015 to 2021 is treated
#'   as a seven-period record rather than as a 22-period record with fifteen
#'   silent years in front of it. Units are otherwise independent: nothing is
#'   pooled, and the kernels are the same for all of them.
#'
#' @return A named list of `lag_ceiling` objects, of class
#'   `lag_ceiling_list`, with the skipped units and the reason in
#'   `attr(, "skipped")`. It can be passed straight to
#'   [rain_strip()], and `as.data.frame()` gives one row per unit.
#'
#' @section Scale and how to read it:
#'
#' The columns of `as.data.frame()` carry the scales of [lag_ceiling()]:
#' `ceiling` and `observed` are concentration indices bounded above by
#' \eqn{(n-1)/n} for the Gini, where `n` is that unit's own number of scored
#' periods, and `exceedance` is their ratio, centred on 1.
#'
#' Because `n` differs between units, two cautions apply to reading the table
#' or the figure down a column. A unit with few periods has a lower attainable
#' maximum, so its Gini is not on the same footing as a longer unit's, and
#' `n` is reported for that reason. And the exceedance is a ratio of two
#' numbers computed on the same `n`, so it is the safer quantity to compare
#' across units of different lengths.
#'
#' What varies between units is usually the ceiling rather than the
#' observation: the ceiling is a property of that unit's reproductive record,
#' so a unit with seasonal or patchy reproduction sits higher and has less to
#' detect. A column of exceedances close to 1 is not a weak result, it is a
#' statement that those units' reproduction is itself episodic.
#'
#' @examples
#' ces <- lag_ceiling_by(lepanthes_hosts, unit = "host",
#'                       reproduction = "inflorescences", K = 4,
#'                       kernels = lag_kernels(4, bin = 6, unit = "mo"))
#' ces
#' as.data.frame(ces)
#' @export
lag_ceiling_by <- function(data, unit = "unit", period = "period",
                           reproduction = "reproduction", recruits = "recruits",
                           K, lag0 = 1L, kernels = lag_kernels(K),
                           missing = c("backfill", "drop"),
                           min_periods = K + 3L, min_recruits = 2) {
  missing <- match.arg(missing)
  data <- as.data.frame(data)
  need <- c(unit, period, reproduction, recruits)
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("column(s) not found in data: ", paste(miss, collapse = ", "))

  parts <- split(data, as.character(data[[unit]]))
  out <- list(); skipped <- character(0)
  for (nm in names(parts)) {
    s <- parts[[nm]]
    s <- s[!is.na(s[[period]]), , drop = FALSE]
    scored <- sum(!is.na(s[[recruits]]))
    nrec <- sum(s[[recruits]], na.rm = TRUE)
    if (nrow(s) < min_periods) {
      skipped[nm] <- sprintf("%d periods, fewer than min_periods = %d", nrow(s), min_periods); next
    }
    if (scored == 0) { skipped[nm] <- "no period with a recruit count"; next }
    if (nrec < min_recruits) {
      skipped[nm] <- sprintf("%g recruits, fewer than min_recruits = %g", nrec, min_recruits); next
    }
    s[[period]] <- as.integer(s[[period]]) - min(as.integer(s[[period]])) + 1L
    out[[nm]] <- lag_ceiling(s, K = K, lag0 = lag0, kernels = kernels, missing = missing,
                             period = period, reproduction = reproduction, recruits = recruits)
  }
  if (!length(out)) stop("no unit met min_periods and min_recruits")
  structure(out, class = "lag_ceiling_list", skipped = skipped)
}

#' @export
print.lag_ceiling_list <- function(x, digits = 3, ...) {
  cat(sprintf("Ceiling for %d units, K = %d, lag0 = %d\n",
              length(x), x[[1]]$K, x[[1]]$lag0))
  print(as.data.frame(x), row.names = FALSE, digits = digits)
  sk <- attr(x, "skipped")
  if (length(sk)) {
    cat(sprintf("\n%d unit(s) skipped:\n", length(sk)))
    for (nm in names(sk)) cat(sprintf("  %-12s %s\n", nm, sk[[nm]]))
  }
  invisible(x)
}

#' @export
as.data.frame.lag_ceiling_list <- function(x, ...) {
  do.call(rbind, lapply(names(x), function(nm) {
    o <- x[[nm]]
    data.frame(unit = nm,
               n = length(o$R),
               mean_recruits = mean(o$R),
               ceiling_gini = unname(o$ceiling["gini"]),
               observed_gini = unname(o$obs["gini"]),
               exceedance_gini = unname(o$exceedance["gini"]),
               exceedance_cv = unname(o$exceedance["cv"]),
               best_kernel = o$best$label,
               stringsAsFactors = FALSE)
  }))
}
