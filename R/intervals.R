#' Check that the census intervals are what the analysis assumes
#'
#' Everything in this package is indexed by period, and a period is assumed to
#' be a fixed length of time: bin 6 months, and lag bin 3 means 18 to 24 months
#' back. Nothing in a table of counts reveals whether that assumption holds,
#' because a data frame of periods 1 to 9 looks the same whether the visits
#' were six months apart or scattered between two months and eighteen. If the
#' study recorded the dates of its visits, this function reads them and says.
#'
#' Run it before [lag_ceiling()] on any record whose visits were not on a
#' schedule, and always on a record you did not collect yourself.
#'
#' @param data A data frame, or a vector of dates.
#' @param date Name of the date column, when `data` is a data frame. The column
#'   may be `Date`, `POSIXct`, or a character vector in ISO order
#'   (`"2001-06-01"`).
#' @param period Name of the period column. If given, one date is expected per
#'   period and the function complains if a period carries more than one.
#' @param persistence How long the thing being counted stays countable, in the
#'   same units as `unit`: about 1.5 months for a *Lepanthes* fruit, 8 for an
#'   *Encyclia* fruit, indefinitely for an adult plant. Optional. When given,
#'   the ratio of persistence to interval is reported for each interval, which
#'   is what decides whether that column was well observed.
#' @param unit Time unit for the report, `"mo"` (default), `"yr"` or `"day"`.
#' @param tolerance An interval is called irregular when it differs from the
#'   median interval by more than this fraction of the median. The default 0.2
#'   allows a visit a fortnight early in a six-monthly schedule.
#'
#' @return An object of class `rain_intervals`, invisibly a data frame with one
#'   row per interval: the two periods, the two dates, the gap, the ratio to
#'   the median gap, and, when `persistence` is given, the persistence to
#'   interval ratio. `print()` gives the report.
#'
#' @section Scale and how to read it:
#'
#' Three numbers matter, and each answers a different question.
#'
#' **The ratio to the median interval** says whether the period index is
#' honest. At 1 the visit fell where the schedule says. Intervals below about
#' 0.8 or above about 1.25 mean the period index is not a time axis, and both
#' the reproduction and the recruit series are then partly a record of how long
#' the crew waited rather than of what the plants did. A short interval catches
#' fewer recruits and fewer fruits than a long one for reasons that have
#' nothing to do with biology, and no kernel can undo that: the convolution is
#' defined on equally spaced bins. There is no correction inside this package
#' for unequal intervals, and dividing the counts by the interval length is not
#' one, because a count divided by time is a rate and the ceiling argument is
#' about counts.
#'
#' **A negative or zero gap** is a data error, not a short interval: the visits
#' are out of order, or a year has been mistyped. It is reported as an error
#' rather than a warning because every downstream number would be meaningless.
#'
#' **The persistence to interval ratio** says whether the column was well
#' observed, and it is the criterion for choosing between candidate
#' reproductive measures in [rain_indices()]. At 1 or above, nothing appears
#' and disappears between visits and the count is complete. At 0.25, a fruit
#' that persists six weeks against a six-monthly visit, roughly three in four
#' are never seen, and the series is erratic through being badly seen rather
#' than through being episodic. That inflates its concentration, inflates the
#' ceiling, and lowers the exceedance, so a badly observed index makes the
#' check conservative: a low exceedance computed from one is not evidence that
#' the lag hypothesis survives. Note the direction. Poor observation protects
#' the hypothesis, so it is not a reason to distrust a refutation.
#'
#' A record whose intervals are unequal is not useless. It can still support
#' the reproductive comparison of [rain_indices()] within a single interval,
#' and it can support anything asked of the whole record at once. What it
#' cannot support is a lag, because a lag is a statement about a number of
#' periods and the periods are not the same length.
#'
#' @examples
#' # A six-monthly record, as the census intended
#' rain_intervals(data.frame(period = 1:5,
#'                           date = as.Date(c("2000-01-01", "2000-07-01", "2001-01-01",
#'                                            "2001-07-01", "2002-01-01"))),
#'                persistence = 1.5)
#' @export
rain_intervals <- function(data, date = "date", period = NULL,
                           persistence = NULL, unit = c("mo", "yr", "day"),
                           tolerance = 0.2) {
  unit <- match.arg(unit)
  per_day <- c(mo = 30.4375, yr = 365.25, day = 1)[[unit]]
  unit_long <- c(mo = "months", yr = "years", day = "days")[[unit]]

  if (inherits(data, c("Date", "POSIXct", "POSIXt"))) {
    d <- as.Date(data); p <- seq_along(d)
  } else {
    data <- as.data.frame(data)
    if (!date %in% names(data))
      rl_abort('There is no column named "', date, '" in the data.',
               "\nThe data frame has: ", rl_list(names(data), 12),
               '.\nIf the survey dates are in another column, name it: date = "<column>".',
               "\nIf the dates were never recorded, this check cannot be run, and ",
               "whether the intervals are equal has to be established from the field notes.")
    d <- data[[date]]
    # With no period column and a long data frame, every unit contributes a row
    # per census, so the dates repeat. Taking the distinct dates in order is
    # then right and taking the row index is not.
    p <- if (is.null(period)) match(d, sort(unique(d))) else {
      if (!period %in% names(data))
        rl_abort('There is no column named "', period, '" in the data.',
                 "\nThe data frame has: ", rl_list(names(data), 12), ".")
      data[[period]]
    }
  }

  if (is.character(d) || is.factor(d)) {
    ch <- as.character(d)
    ch[!is.na(ch) & trimws(ch) == ""] <- NA   # an unrecorded date read from csv
    d2 <- suppressWarnings(as.Date(ch))
    if (any(is.na(d2) & !is.na(ch)))
      rl_abort('The date column "', date, '" could not be read as dates. ',
               "The first value that failed is ",
               rl_list(ch[is.na(d2) & !is.na(ch)][1]),
               '.\nDates must be in year-month-day order, for example "2001-06-01". ',
               'Convert them first, for example with as.Date(x, format = "%d/%m/%Y").')
    # as.Date() will happily read "01/02/2000" as year 1, month 2, day 20,
    # which is worse than failing because the analysis then proceeds on
    # nonsense. Anything outside a plausible range of field-study years is
    # treated as a format that was misread rather than as a date.
    yr <- as.integer(format(d2, "%Y"))
    if (any(yr < 1900 | yr > 2100, na.rm = TRUE)) {
      i <- which(yr < 1900 | yr > 2100)[1]
      rl_abort('The date column "', date, '" was read as dates outside any plausible ',
               "range: ", rl_list(ch[i]), " came out as ", as.character(d2[i]), ".",
               '\nThis happens when dates are written day-first or month-first. ',
               'as.Date() assumes year-first, so "01/02/2000" becomes the year 1.',
               '\nGive the format explicitly before calling this, for example ',
               'data$', date, ' <- as.Date(data$', date, ', format = "%d/%m/%Y").')
    }
    d <- d2
  }
  d <- as.Date(d)
  if (!length(d) || all(is.na(d)))
    rl_abort('The date column "', date, '" is empty or all NA, so the intervals ',
             "cannot be checked.")

  # One date per period. A long data frame has one row per unit and period, so
  # the same date repeats legitimately; two DIFFERENT dates in one period does
  # not, and is worth saying plainly because it is usually a merge gone wrong.
  # A period with no recorded date is kept, not dropped. Dropping it would
  # silently merge the two intervals on either side into one and report that
  # merged gap as if it were a single census interval.
  keep <- !is.na(p)
  d <- d[keep]; p <- p[keep]
  parts <- lapply(split(d, p), function(x) x[!is.na(x)])
  undated <- names(parts)[lengths(parts) == 0]
  n_multi <- vapply(parts, function(x) length(unique(x)), integer(1))
  if (any(n_multi > 1)) {
    bad <- names(n_multi)[n_multi > 1]
    rl_abort("These periods carry more than one date: ", rl_list(bad), ".",
             "\nA period is one visit, so it has one date. Period ", bad[1], " has ",
             n_multi[[bad[1]]], ": ", rl_list(sort(unique(parts[[bad[1]]]))), ".",
             "\nIf the visits within a period really were on different days, use the ",
             "first or the median of them, and say which in the write-up.")
  }
  per <- suppressWarnings(as.numeric(names(parts)))
  o <- order(if (all(!is.na(per))) per else seq_along(parts))
  per <- per[o]
  dts <- as.Date(vapply(parts, function(x) if (length(x)) as.character(x[1]) else NA_character_,
                        character(1))[o])
  if (sum(!is.na(dts)) < 2)
    rl_abort("There are fewer than two census dates, so there is no interval to check.")

  gap <- as.numeric(diff(dts)) / per_day
  med <- stats::median(gap[!is.na(gap) & gap > 0])
  out <- data.frame(from = per[-length(per)], to = per[-1],
                    date_from = dts[-length(dts)], date_to = dts[-1],
                    interval = round(gap, 2),
                    ratio_to_median = round(gap / med, 2),
                    stringsAsFactors = FALSE)
  if (!is.null(persistence)) {
    if (!is.numeric(persistence) || length(persistence) != 1 || persistence <= 0)
      rl_abort("`persistence` must be a single positive number, in ", unit_long,
               ": how long the thing being counted stays countable. ",
               "About 1.5 for a Lepanthes fruit at unit = \"mo\".")
    out$observed_fraction <- round(pmin(1, persistence / gap), 2)
  }

  bad <- which(!is.na(gap) & gap <= 0)
  if (length(bad))
    rl_abort("The censuses are not in chronological order. ",
             if (length(bad) == 1) "This interval goes " else "These intervals go ",
             "backwards or nowhere:\n  ",
             paste(sprintf("period %s (%s) to period %s (%s): %+.1f %s",
                           out$from[bad], out$date_from[bad], out$to[bad],
                           out$date_to[bad], gap[bad], unit_long), collapse = "\n  "),
             "\nA later period must have a later date. This is usually a mistyped year ",
             "in one row, or two surveys numbered the wrong way round. ",
             "Fix the dates before anything else: every lag in the analysis is counted ",
             "in periods, so an out-of-order period makes all of them wrong.")

  irr <- !is.na(gap) & abs(gap / med - 1) > tolerance
  structure(out, class = c("rain_intervals", "data.frame"),
            median_interval = med, unit = unit_long, tolerance = tolerance,
            irregular = irr, persistence = persistence, undated = undated)
}

#' @export
print.rain_intervals <- function(x, ...) {
  u   <- attr(x, "unit"); med <- attr(x, "median_interval")
  irr <- attr(x, "irregular"); tol <- attr(x, "tolerance")
  cat(sprintf("Census intervals: %d intervals over %d censuses, median %.2g %s\n",
              nrow(x), nrow(x) + 1, med, u))
  print(as.data.frame(x), row.names = FALSE)
  und <- attr(x, "undated")
  if (length(und))
    cat(sprintf(paste("\n%d census(es) have no recorded date: %s. The intervals on either",
                      "side\nare unknown, not zero, and are shown as NA. Every lag counted",
                      "across such a\ncensus is unverifiable.\n"),
                length(und), rl_list(und)))
  if (any(irr)) {
    cat(sprintf("\n%d of %d intervals differ from the median by more than %.0f%%: %s.\n",
                sum(irr), nrow(x), 100 * tol,
                paste(sprintf("%s to %s (%.2g %s)", x$from[irr], x$to[irr],
                              x$interval[irr], u), collapse = ", ")))
    cat("The period index is therefore not a time axis, and a lag counted in periods\n",
        "does not correspond to a fixed delay. Counts are not comparable between\n",
        "intervals either, since a longer interval accumulates more of everything.\n",
        "There is no correction for this in the package: report the whole record, or\n",
        "analyse a run of intervals that really are equal.\n", sep = "")
  } else {
    cat(sprintf("\nAll intervals are within %.0f%% of the median, so the period index is a time axis.\n",
                100 * tol))
  }
  p <- attr(x, "persistence")
  if (!is.null(p)) {
    f <- x$observed_fraction[!is.na(x$observed_fraction)]
    cat(sprintf("\nA structure persisting %.2g %s is fully counted at %d of %d dated intervals",
                p, u, sum(f >= 1), length(f)))
    if (any(f < 1))
      cat(sprintf(", and at the worst\nof them about %.0f%% of it is seen.\n", 100 * min(f)))
    else cat(".\n")
    if (any(f < 0.5))
      cat("Below about a half, the column is erratic through being badly seen rather than\n",
          "through being episodic. That raises its ceiling and lowers the exceedance, so it\n",
          "makes the check conservative: a low exceedance from such a column is not evidence\n",
          "that the lag hypothesis survives.\n", sep = "")
  }
  invisible(x)
}
