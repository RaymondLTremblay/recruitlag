#' Build the unit-by-period matrices from a long data frame
#'
#' The host-level functions work on two matrices, recruits (units by scored
#' periods) and reproduction (units by all periods). Most census data live in
#' a long data frame with one row per unit and period; this function makes
#' the matrices from it. Periods at which `recruits` is `NA` for every unit
#' are taken as unscored.
#'
#' @param data A data frame with one row per unit and period.
#' @param unit,period,reproduction,recruits Names of the columns holding the
#'   unit identifier, the period index (integer, from 1), the reproductive
#'   count and the recruit count.
#' @return A list with `R` (recruits, units by scored periods, `colnames` the
#'   period indices) and `X` (reproduction, units by periods 1 to the last
#'   period in `data`). Missing reproduction is treated as zero.
#' @examples
#' m <- host_matrices(lepanthes_hosts, unit = "host", reproduction = "inflorescences")
#' dim(m$R); dim(m$X)
#' @export
host_matrices <- function(data, unit = "unit", period = "period",
                          reproduction = "reproduction", recruits = "recruits") {
  data <- as.data.frame(data)
  rl_check_long(data, period, reproduction, recruits, unit = unit, what = "host_matrices()")
  u <- as.character(data[[unit]]); p <- as.integer(data[[period]])
  units <- unique(u); periods <- seq_len(max(p))
  X <- matrix(0, length(units), length(periods), dimnames = list(units, periods))
  X[cbind(match(u, units), p)] <- ifelse(is.na(data[[reproduction]]), 0, data[[reproduction]])
  Rfull <- matrix(NA_real_, length(units), length(periods), dimnames = list(units, periods))
  Rfull[cbind(match(u, units), p)] <- data[[recruits]]
  scored <- which(colSums(!is.na(Rfull)) > 0)
  if (!length(scored))
    rl_abort('No period has a recruit count: "', recruits, '" is NA in every row.')
  R <- Rfull[, scored, drop = FALSE]
  R[is.na(R)] <- 0
  list(R = R, X = X)
}

#' Pull the two system-wide series from a data frame
#'
#' Sums (or takes, if there is one row per period) reproduction and recruits
#' by period. Periods at which `recruits` is `NA` are unscored.
#' @param data A data frame with one row per period, or per unit and period.
#' @param period,reproduction,recruits Names of the columns holding the
#'   period index (integer, from 1), the reproductive count and the recruit
#'   count.
#' @return A list with `X` (reproduction at periods 1 to the last), `R`
#'   (recruits at the scored periods) and `t_R` (those periods).
#' @examples
#' series_from(lepanthes_census, reproduction = "inflorescences")
#' @export
series_from <- function(data, period = "period", reproduction = "reproduction",
                        recruits = "recruits") {
  data <- as.data.frame(data)
  rl_check_long(data, period, reproduction, recruits, what = "series_from()")
  p <- as.integer(data[[period]]); periods <- seq_len(max(p))
  X <- tapply(ifelse(is.na(data[[reproduction]]), 0, data[[reproduction]]), factor(p, periods), sum)
  X[is.na(X)] <- 0
  Rs <- tapply(data[[recruits]], factor(p, periods), function(v) if (all(is.na(v))) NA_real_ else sum(v, na.rm = TRUE))
  t_R <- which(!is.na(Rs))
  list(X = as.numeric(X), R = as.numeric(Rs[t_R]), t_R = as.integer(t_R))
}
