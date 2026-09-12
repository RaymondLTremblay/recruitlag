#' Census records of the epiphytic orchid *Lepanthes eltoroensis*
#'
#' Four tibbles from the two surveys analysed in the companion paper, in
#' which reproduction is continuous in time and recruitment is episodic.
#' Luquillo Mountains, Puerto Rico.
#'
#' @details
#' `lepanthes_monthly` (23 rows): the monthly survey of four populations on
#' the El Toro trail, September 1994 to November 1996. `period` is the census
#' index, `reproductive_adults` the number of adults bearing an inflorescence
#' (summed over populations) and `recruits` the new seedlings first seen at
#' that census (`NA` at census 1, where first appearance is undefined).
#'
#' `lepanthes_census` (12 rows): system totals of the six-monthly survey of 23
#' host trees begun after Hurricane Georges, March 1999 to September 2004:
#' `period`, `census` (label), `years_since_georges`, stage counts (`adults`,
#' `juveniles`, `seedlings`), reproductive totals (`inflorescences`,
#' `flowers`, `fruits`) and `recruits`, scored at censuses 3 to 12 (`NA` at
#' censuses 1 and 2, which are the initial enumeration and relocation
#' additions).
#'
#' `lepanthes_hosts` (276 rows): the same survey with one row per host tree
#' and census: `host` (identifier), `period`, `adults`, `inflorescences` and
#' `recruits` (`NA` at censuses 1 and 2). This is the input of
#' [host_lag_test()], through `unit = "host"` and
#' `reproduction = "inflorescences"`.
#'
#' `lepanthes_profile` (5 rows): the posterior-mean lag profile of the
#' distributed-lag model fitted in the companion paper, one row per six-month
#' bin (`bin`, `months`, `weight`), weights summing to one. An example of a
#' fitted profile to pass to [lag_kernels()] through `extra`.
#'
#' @source Tremblay, R. L. Sex every month, seedlings once in a while: can
#'   episodic recruitment be a delayed consequence of continuous
#'   reproduction? (companion paper); Tremblay & Ackerman (2001) Biological
#'   Journal of the Linnean Society 72: 47-62 for the monthly survey.
#' @examples
#' lepanthes_census
#' head(lepanthes_hosts)
#' @name lepanthes
#' @aliases lepanthes_monthly lepanthes_census lepanthes_hosts lepanthes_profile
NULL

#' Six invented census records spanning the cases the test separates
#'
#' Constructed records, not data, in one long tibble. Three reproductive
#' regimes (seasonal, intermediate, aseasonal) crossed with two recruitment
#' mechanisms. Under `"delayed"`, recruits on each unit are a true delayed
#' consequence of that unit's own reproduction: a fixed lag profile over five
#' bins (weights 0.10, 0.25, 0.35, 0.20, 0.10 at one to five periods back)
#' with negative-binomial counts. Under `"gated"`, recruits arrive only at
#' three of the 24 periods, when an external gate is open, at a rate
#' proportional to each unit's average reproduction, so the record is
#' episodic and does not track reproduction. Each record has 20 units and 24
#' periods; recruits are scored at periods 6 to 24 (`K = 4`, `lag0 = 1`) and
#' are `NA` elsewhere. The generating script is `data-raw/regimes.R`.
#'
#' @format `regimes` has one row per record, unit and period (2,880 rows):
#'   `record` (`seasonal_delayed`, `seasonal_gated`, `intermediate_delayed`,
#'   `intermediate_gated`, `aseasonal_delayed`, `aseasonal_gated`), `label`,
#'   `seasonality` (0, 0.5 or 1), `mechanism`, `unit`, `period`,
#'   `reproduction` and `recruits`. `regimes_truth` holds the generating
#'   rules: the five weights of the true lag profile and the 24 gate values.
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' r <- filter(regimes, record == "aseasonal_gated")
#' ce <- lag_ceiling(r, K = 4, kernels = lag_kernels(4))
#' ce
#' @name regimes
#' @aliases regimes_truth
NULL
