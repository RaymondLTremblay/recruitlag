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
#' @details Rows in which both reproduction and recruits are missing are
#'   taken as periods in which that unit was not watched, and are dropped
#'   before its periods are re-indexed to start at 1. A unit followed from
#'   2015 to 2021 is therefore a seven-period record, not a 22-period record
#'   with fifteen silent years in front of it. Keeping them would read those
#'   years as periods with no reproduction, which raises the unit's ceiling
#'   and makes the test look weaker than it is. Units are otherwise
#'   independent: nothing is pooled, and the kernels are the same for all of
#'   them.
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
  rl_check_long(data, period, reproduction, recruits, unit = unit, what = "lag_ceiling_by()")
  rl_check_kernels(kernels, K)

  # A row carrying neither reproduction nor recruits is a period in which this
  # unit was not watched. Such rows are dropped before the periods are
  # re-indexed, so a unit followed from 2015 to 2021 is a seven-period record
  # and not a 22-period record with fifteen silent years in front of it. They
  # would otherwise be read as periods with no reproduction, which raises that
  # unit's ceiling and makes the test look weaker than it is.
  blank <- is.na(data[[reproduction]]) & is.na(data[[recruits]])
  dropped_rows <- sum(blank)
  data <- data[!blank, , drop = FALSE]
  if (!nrow(data))
    rl_abort("Every row has neither reproduction nor recruits, so there is nothing to do.")

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
    ce <- tryCatch(suppressWarnings(
            lag_ceiling(s, K = K, lag0 = lag0, kernels = kernels, missing = missing,
                        period = period, reproduction = reproduction, recruits = recruits)),
                   error = function(e) e)
    if (inherits(ce, "error")) { skipped[nm] <- conditionMessage(ce); next }
    out[[nm]] <- ce
  }
  if (!length(out))
    rl_abort("No unit could be given a ceiling. All ", length(parts),
             " were skipped, for these reasons:\n  ",
             paste(sprintf("%s: %s", names(skipped), skipped), collapse = "\n  "),
             "\nLower min_periods or min_recruits if the units really are this short.")
  structure(out, class = "lag_ceiling_list", skipped = skipped, dropped_rows = dropped_rows)
}

#' @export
print.lag_ceiling_list <- function(x, digits = 3, ...) {
  cat(sprintf("Ceiling for %d units, K = %d, lag0 = %d\n",
              length(x), x[[1]]$K, x[[1]]$lag0))
  print(as.data.frame(x), row.names = FALSE, digits = digits)
  dr <- attr(x, "dropped_rows")
  if (!is.null(dr) && dr > 0)
    cat(sprintf("  %d row(s) dropped as unwatched periods (no reproduction and no recruits)\n", dr))
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

#' Compare the ceiling across candidate reproductive measures
#'
#' A census usually offers several measures of reproduction: adults,
#' reproductive plants, inflorescences, flowers, fruits, and rates derived
#' from them. They do not give the same ceiling, and the choice can change
#' the exceedance several fold, so it is a decision to report rather than to
#' make quietly. `rain_indices()` runs [lag_ceiling()] once per candidate
#' column and tabulates the result.
#'
#' @param data A data frame with one row per period, or per unit and period,
#'   which is summed.
#' @param reproduction Character vector of candidate reproductive columns.
#' @param period,recruits Column names.
#' @param K,lag0,kernels,missing As in [lag_ceiling()]. The same kernels are
#'   used for every candidate, so the comparison is like for like.
#'
#' @return A data frame with one row per candidate: the concentration of that
#'   reproductive series, the searched ceiling, the exceedance on both
#'   indices, the kernel that attained the ceiling, and how many values of
#'   that column were missing (and so read as zero).
#'
#' @section Scale and how to read it:
#'
#' Read the table down the `exceedance` column, and read it as a sensitivity
#' analysis rather than as a menu. Two things drive the differences.
#'
#' **A stock gives a lower ceiling than a flux.** Adults change only by
#' recruitment minus death, so the adult series is smooth whatever the plants
#' are doing, and the exceedance computed from it is the largest available.
#' It is the most favourable index to the paper's own argument and therefore
#' the one to trust least.
#'
#' **Measurement error runs the other way, and protects the lag
#' hypothesis.** A reproductive measure that is caught only sometimes, such
#' as a fruit that persists a month against a six-month census, looks more
#' erratic than the truth. That raises its concentration, raises the ceiling,
#' and lowers the exceedance. So a badly observed index makes the test
#' conservative, and a low exceedance from such a column is not evidence that
#' the lag hypothesis survives.
#'
#' The consequence for reporting: give the whole table, say which columns are
#' well observed in that particular study, and let the reader see the range.
#' No index is the right one in general.
#'
#' @examples
#' rain_indices(lepanthes_census, c("adults", "inflorescences", "flowers", "fruits"),
#'              K = 4, kernels = lag_kernels(4, bin = 6, unit = "mo"))
#' @export
rain_indices <- function(data, reproduction, period = "period", recruits = "recruits",
                         K, lag0 = 1L, kernels = lag_kernels(K),
                         missing = c("backfill", "drop")) {
  missing <- match.arg(missing)
  data <- as.data.frame(data)
  if (!is.character(reproduction) || !length(reproduction))
    rl_abort("`reproduction` must be a character vector naming the candidate columns, ",
             'for example c("n_adults", "n_inflorescences", "n_fruits").')
  rl_check_columns(data, stats::setNames(c(period, recruits), c("period", "recruits")))
  miss <- setdiff(reproduction, names(data))
  if (length(miss))
    rl_abort("These candidate reproductive columns are not in the data: ",
             rl_list(miss), ".\nThe data frame has: ", rl_list(names(data), 12), ".")
  rl_check_kernels(kernels, K)

  do.call(rbind, lapply(reproduction, function(v) {
    ce <- tryCatch(suppressWarnings(
      lag_ceiling(data, period = period, reproduction = v, recruits = recruits,
                  K = K, lag0 = lag0, kernels = kernels, missing = missing)),
      error = function(e) e)
    if (inherits(ce, "error"))
      return(data.frame(index = v, gini = NA_real_, ceiling_gini = NA_real_,
                        exceedance_gini = NA_real_, exceedance_cv = NA_real_,
                        best_kernel = "no ceiling: lagged reproduction is zero throughout",
                        n_missing = sum(is.na(data[[v]])),
                        stringsAsFactors = FALSE))
    data.frame(index = v,
               gini = unname(ce$theorem["gini"]),
               ceiling_gini = unname(ce$ceiling["gini"]),
               exceedance_gini = unname(ce$exceedance["gini"]),
               exceedance_cv = unname(ce$exceedance["cv"]),
               best_kernel = ce$best$label,
               n_missing = sum(is.na(data[[v]])),
               stringsAsFactors = FALSE)
  }))
}
