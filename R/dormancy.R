# Is a first sighting a recruit? The triage a tagged-plant census needs before
# its first sightings are called recruits: what state plants were in when
# first seen, how much dormancy the record shows, how long each plant had
# been watched for, and how many first sightings the record's own dormancy
# would produce from plants that were there all along.

# ---- the input: a plant-by-period state matrix ------------------------------

# Accepts a matrix (plants by periods) or a long data frame (plant, period,
# state) and returns a character matrix of states with plants as rows and
# periods 1..T as columns, plus the period labels.
state_matrix <- function(x, plant = "plant", period = "period", state = "state") {
  if (is.data.frame(x) && all(c(plant, period, state) %in% names(x))) {
    if (!nrow(x)) rl_abort("The data frame has no rows.")
    p <- x[[period]]
    if (!is.numeric(p))
      rl_abort('The period column "', period, '" must be numeric (a year or a census index). It is ',
               class(p)[1], ".")
    periods <- sort(unique(p[!is.na(p)]))
    if (any(diff(periods) != diff(periods)[1]))
      rl_warn("The periods are not evenly spaced (", rl_list(periods, 8), "). A gap is read as a ",
              "period in which nothing was seen, so dormancy runs will be overstated across it. ",
              "If the census skipped a year, add a row of NA for every plant at that period.")
    u <- as.character(x[[plant]]); plants <- unique(u[!is.na(u)])
    M <- matrix(NA_character_, length(plants), length(periods), dimnames = list(plants, periods))
    ok <- !is.na(u) & !is.na(p)
    M[cbind(match(u[ok], plants), match(p[ok], periods))] <- as.character(x[[state]][ok])
    return(list(M = M, periods = periods))
  }
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x))
    rl_abort("x must be a matrix with one row per plant and one column per period, or a long ",
             'data frame with columns "', plant, '", "', period, '" and "', state, '". It is ',
             class(x)[1], ".")
  M <- matrix(as.character(x), nrow(x), ncol(x), dimnames = dimnames(x))
  if (is.null(rownames(M))) rownames(M) <- paste0("plant", seq_len(nrow(M)))
  periods <- suppressWarnings(as.numeric(gsub("[^0-9]", "", colnames(M))))
  if (is.null(colnames(M)) || anyNA(periods)) periods <- seq_len(ncol(M))
  list(M = M, periods = periods)
}

# ---- the triage -------------------------------------------------------------

#' Are the first sightings recruits? Dormancy, detection and the discovery series
#'
#' A tagged-plant census records, for each plant, the state it was in at each
#' period: flowering, vegetative, not seen. The recruit series that such a
#' record offers is the number of plants seen for the first time at each
#' period, and that series is a recruitment series only if a first sighting
#' is a new plant. Two things break the equivalence. Dormancy: a plant that
#' was there all along can stay underground, or be missed, for a year or
#' more and then be "discovered". Search: a plot that is enlarged, or
#' searched harder, discovers established plants in bulk. This function puts
#' the numbers a reader needs to judge a first-sighting series beside it.
#'
#' @param x A matrix with one row per plant and one column per period, or a
#'   long data frame with one row per plant and period.
#' @param plant,period,state Column names when `x` is a long data frame.
#' @param seen The state codes that mean the plant was above ground. The
#'   default covers `F`/`V` and `1`/`2` codings. Anything else (`NA`, `"NA"`,
#'   `"D"`, `0`, blank) is read as not seen.
#' @param flowering The codes, within `seen`, that mean the plant was
#'   flowering.
#' @param s,r Annual survival and the probability of being seen given
#'   survival, used to compute how many first sightings the record's own
#'   dormancy would produce from plants present at the first period. Supply
#'   them from a capture-recapture study of the species (for nine
#'   *Caladenia* species, [caladenia_dormancy] carries the published values)
#'   or leave `r` as `NULL` to use the record's own estimate (below) with
#'   the `s` you supply; leave both `NULL` and that table is not computed.
#' @param species A row name of [caladenia_dormancy], as a shortcut for `s`
#'   and `r`.
#'
#' @return An object of class `recruit_triage` with: `first` (one row per
#'   period: plants first seen, split by state at first sighting, the number
#'   of periods each had been watched for, the expected number of first
#'   sightings that are returns of plants present but unseen at the first
#'   period, and the excess over that), `dormancy` (the run-length table of
#'   unseen periods between sightings, as in Coates et al. 2006 Table 2, the
#'   share of plants with at least one such run, and the longest), `detection`
#'   (the record's own estimate of `r`: the share of periods between a
#'   plant's first and last sighting at which it was seen), `watched` (how
#'   many periods the plot had been watched before each first sighting) and
#'   the inputs.
#'
#' @section Scale and how to read it:
#'
#' **State at first sighting.** A plant that is flowering when first seen
#' has, in most orchids, been established for years. A first-sighting series
#' in which most plants arrive flowering is a discovery series, not a
#' recruitment series, and the share flowering is the first number to read.
#' In *Lepanthes rubripetala* it was 29 of 43; in the Victorian terrestrials
#' it ran from 57% to 96%.
#'
#' **The dormancy runs.** The table counts, for every plant, every run of
#' consecutive unseen periods that falls between two sightings, so each run
#' is a plant known to have been alive and present while unseen. The share
#' of plants with at least one run, and the longest run, say how much a
#' first sighting after a few years of watching can be trusted. Runs at the
#' end of a plant's record are not counted: a plant not seen again may be
#' dead, and the record cannot tell.
#'
#' **The record's own detection estimate.** `detection` is the share of
#' periods, between each plant's first and last sighting, at which it was
#' seen: an estimate of `r`, the probability of being seen given survival,
#' by the conventional method (Kery et al. 2005; Tremblay et al. 2009). It
#' is biased upward, because a plant is only known to be alive while it is
#' being resighted, and it says nothing about survival.
#'
#' **Expected returns.** With `s` and `r` given, a plant alive at period 1
#' was unseen there with probability `1 - r`, survives and stays unseen for a
#' further period with probability `q = s (1 - r)`, and is seen at period `k`
#' with probability `s r` after `k - 2` such periods. If `N1` plants were
#' seen at period 1, then `N1 / r` were alive, `N1 (1 - r) / r` of them were
#' unseen, and the expected number of those discovered at period `k` is
#' `N1 (1 - r) s q^(k - 2)`. That is the number of first sightings that the
#' record's own dormancy accounts for, and it falls fast: for *Caladenia
#' valida* (survival 0.82, resighting 0.73) it is 22% of `N1` at the second period,
#' 5% at the third and 1% at the fourth. The `excess` column is what is left
#' after subtracting it, and is the most that recruitment plus search effort
#' can have contributed. A large excess in a period in which the plot was
#' enlarged is search; a large excess after several settled years of
#' watching is the nearest thing to a recruit count this kind of record can
#' give. Returns of plants that arrived after period 1 are not modelled:
#' those are new plants whenever they are first seen, only late.
#'
#' **What this does not do.** It does not make a first-sighting series into
#' a recruitment series. It says how far the two are apart, and where. A
#' record with a seedling or protocorm class does not need it.
#'
#' @references
#' Coates, F., Lunt, I. D. and Tremblay, R. L. (2006) Effects of disturbance
#' on population dynamics of the threatened orchid *Prasophyllum correctum*
#' D.L. Jones and implications for grassland management in south-eastern
#' Australia. *Biological Conservation* 129: 59-69.
#' \doi{10.1016/j.biocon.2005.06.037}
#'
#' Kery, M., Gregg, K. B. and Schaub, M. (2005) Demographic estimation
#' methods for plants with unobservable life-states. *Oikos* 108: 307-320.
#' \doi{10.1111/j.0030-1299.2005.13589.x}
#'
#' Tremblay, R. L., Perez, M.-E., Larcombe, M., Brown, A., Quarmby, J.,
#' Bickerton, D., French, G. and Bould, A. (2009) Dormancy in *Caladenia*: a
#' Bayesian approach to evaluating latency. *Australian Journal of Botany*
#' 57: 340-350. \doi{10.1071/BT08163}
#'
#' @examples
#' # a small invented record: 30 plants, 8 years, F/V/NA
#' set.seed(7)
#' M <- matrix(NA_character_, 30, 8, dimnames = list(paste0("p", 1:30), 2000:2007))
#' for (i in 1:30) { start <- sample(1:5, 1)
#'   for (j in start:8) M[i, j] <- sample(c("F", "V", NA), 1, prob = c(.4, .35, .25)) }
#' tr <- recruit_triage(M, species = "C. valida")
#' tr
#' rain_sightings(tr)
#' @export
recruit_triage <- function(x, plant = "plant", period = "period", state = "state",
                           seen = c("F", "V", "1", "2"), flowering = c("F", "2"),
                           s = NULL, r = NULL, species = NULL) {
  sm <- state_matrix(x, plant, period, state)
  M <- sm$M; periods <- sm$periods
  H <- nrow(M); Tn <- ncol(M)
  if (Tn < 2) rl_abort("The record has ", Tn, " period(s); at least two are needed to see a first sighting.")
  if (!is.null(species)) {
    cd <- NULL; utils::data("caladenia_dormancy", package = "recruitlag", envir = environment())
    cd <- get("caladenia_dormancy", envir = environment())
    if (!species %in% cd$species)
      rl_abort('species = "', species, '" is not in caladenia_dormancy. It has: ',
               rl_list(cd$species, 9), ".")
    s <- cd$p_survival[cd$species == species]; r <- cd$p_resight[cd$species == species]
  }
  S <- matrix(!is.na(M) & M %in% seen, H, Tn)
  Fl <- matrix(!is.na(M) & M %in% flowering, H, Tn)
  ever <- rowSums(S) > 0
  if (!any(ever))
    rl_abort("No plant was ever seen: none of the states matches `seen` (", rl_list(seen),
             "). The record's codes are: ", rl_list(sort(unique(as.character(M[!is.na(M)])))),
             ". Set `seen` to the codes that mean above ground.")
  first <- apply(S, 1, function(v) if (any(v)) which(v)[1] else NA_integer_)
  last  <- apply(S, 1, function(v) if (any(v)) max(which(v)) else NA_integer_)

  # dormancy runs: unseen periods strictly between first and last sighting
  runs <- integer(0); ever_dormant <- logical(H)
  for (i in which(ever)) {
    if (last[i] - first[i] < 2) next
    v <- S[i, first[i]:last[i]]
    rl <- rle(v); z <- rl$lengths[!rl$values]
    if (length(z)) { runs <- c(runs, z); ever_dormant[i] <- TRUE }
  }
  run_tab <- if (length(runs)) table(factor(pmin(runs, 6L), levels = 1:6,
                                            labels = c("1", "2", "3", "4", "5", ">5"))) else
    table(factor(character(0), levels = c("1", "2", "3", "4", "5", ">5")))
  # detection by the conventional method: seen / periods between first and last sighting
  bracket <- sum(last[ever] - first[ever] + 1L); seen_in <- sum(S[ever, ])
  r_data <- if (bracket > 0) seen_in / bracket else NA_real_

  # first sightings by period, with the state they were in and how long watched
  fs <- data.frame(period = periods, first_seen = 0L, first_seen_flowering = 0L,
                   first_seen_vegetative = 0L, stringsAsFactors = FALSE)
  for (i in which(ever)) {
    j <- first[i]; fs$first_seen[j] <- fs$first_seen[j] + 1L
    if (Fl[i, j]) fs$first_seen_flowering[j] <- fs$first_seen_flowering[j] + 1L
    else fs$first_seen_vegetative[j] <- fs$first_seen_vegetative[j] + 1L
  }
  fs$periods_watched_before <- seq_len(Tn) - 1L
  N1 <- fs$first_seen[1]
  if (is.null(r) && !is.null(s)) r <- r_data
  if (!is.null(s) && !is.null(r)) {
    if (any(c(s, r) <= 0) || any(c(s, r) > 1))
      rl_abort("s and r are probabilities and must lie in (0, 1]. Given s = ", s, ", r = ", r, ".")
    q <- s * (1 - r)
    fs$expected_returns <- c(NA_real_, N1 * (1 - r) * s * q^(seq_len(Tn - 1) - 1))
    fs$excess <- fs$first_seen - fs$expected_returns
    fs$excess[1] <- NA_real_
    # A run of k unseen periods between two sightings has probability q^k for a
    # living plant. Runs so long that q^k < 0.01 are more likely a tag reused or
    # a plant misidentified than a dormancy, and are counted.
    k_star <- if (q > 0 && q < 1) ceiling(log(0.01) / log(q)) else NA_integer_
    improbable <- if (!is.na(k_star)) sum(runs >= k_star) else NA_integer_
  } else { q <- NA_real_; k_star <- NA_integer_; improbable <- NA_integer_ }
  structure(list(first = fs,
                 dormancy = list(runs = run_tab, n_runs = length(runs),
                                 share_plants_ever_dormant = mean(ever_dormant[ever]),
                                 longest_run = if (length(runs)) max(runs) else 0L,
                                 q = q, improbable_from = k_star, improbable_runs = improbable),
                 detection = r_data, s = s, r = r, species = species,
                 n_plants = sum(ever), n_periods = Tn, periods = periods,
                 first_state = c(flowering = sum(fs$first_seen_flowering[-1]),
                                 vegetative = sum(fs$first_seen_vegetative[-1])),
                 M = M, seen = seen, flowering = flowering),
            class = "recruit_triage")
}

#' @export
print.recruit_triage <- function(x, ...) {
  cat(sprintf("Recruit triage: %d plants over %d periods (%s to %s)\n", x$n_plants, x$n_periods,
              format(x$periods[1]), format(x$periods[x$n_periods])))
  fsn <- sum(x$first_state)
  cat(sprintf("  first sightings after period 1: %d, of which %d (%.0f%%) were flowering when first seen\n",
              fsn, x$first_state[["flowering"]], if (fsn) 100 * x$first_state[["flowering"]] / fsn else 0))
  cat(sprintf("  dormancy: %.0f%% of plants have at least one unseen run between sightings; longest run %d period(s)\n",
              100 * x$dormancy$share_plants_ever_dormant, x$dormancy$longest_run))
  rt <- x$dormancy$runs
  cat("  runs of unseen periods (1, 2, 3, 4, 5, >5): ", paste(as.integer(rt), collapse = ", "), "\n", sep = "")
  cat(sprintf("  detection, conventional estimate of r from the record: %.3f\n", x$detection))
  if (!is.null(x$s) && !is.null(x$r)) {
    cat(sprintf("  expected returns computed with s = %.3f, r = %.3f%s\n", x$s, x$r,
                if (!is.null(x$species)) paste0(" (", x$species, ", caladenia_dormancy)") else ""))
    if (!is.na(x$dormancy$improbable_runs))
      cat(sprintf("  runs of %d or more unseen periods have probability below 1%% for a living plant at these rates: %d such run(s), more likely a tag reused or a plant misidentified\n",
                  x$dormancy$improbable_from, x$dormancy$improbable_runs))
  }
  fs <- x$first
  if (!is.null(fs$expected_returns)) {
    fs$expected_returns <- round(fs$expected_returns, 2); fs$excess <- round(fs$excess, 1)
  }
  print(fs, row.names = FALSE)
  cat("  first_seen: plants seen for the first time at that period. expected_returns: how many of them the\n",
      " record's own dormancy accounts for, from plants present but unseen at period 1. excess: the rest,\n",
      " which is recruitment plus search effort, and the record cannot tell which.\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.recruit_triage <- function(x, ...) x$first

#' First sightings by period, against what dormancy alone would produce
#'
#' Bars of the plants first seen at each period, stacked by their state at
#' first sighting, with the expected number of returns of plants present but
#' unseen at the first period drawn over them when `s` and `r` were given.
#'
#' @param x A `recruit_triage` object.
#' @param unit Label for the period axis.
#' @return A ggplot.
#'
#' @section Scale and how to read it:
#'
#' The vertical axis is a count of plants. Bars far above the line, in the
#' first periods, are the plot being searched or enlarged; bars far above it
#' after several settled periods are the closest thing to recruitment this
#' record offers; bars at or below it are dormancy returns and nothing more.
#' The colour split says whether the plants arrived flowering, which in most
#' orchids means they were established already.
#'
#' @examples
#' set.seed(7)
#' M <- matrix(NA_character_, 30, 8, dimnames = list(paste0("p", 1:30), 2000:2007))
#' for (i in 1:30) { start <- sample(1:5, 1)
#'   for (j in start:8) M[i, j] <- sample(c("F", "V", NA), 1, prob = c(.4, .35, .25)) }
#' rain_sightings(recruit_triage(M, species = "C. valida"))
#' @export
rain_sightings <- function(x, unit = "period") {
  rl_check_class(x, "recruit_triage", "x", "recruit_triage()")
  d <- x$first[-1, ]
  long <- rbind(data.frame(period = d$period, state = "flowering when first seen", n = d$first_seen_flowering),
                data.frame(period = d$period, state = "vegetative when first seen", n = d$first_seen_vegetative))
  long$state <- factor(long$state, levels = c("vegetative when first seen", "flowering when first seen"))
  p <- ggplot2::ggplot(long, ggplot2::aes(.data$period, .data$n, fill = .data$state)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_fill_manual(values = c(`vegetative when first seen` = "#56B4E9",
                                          `flowering when first seen` = rl_verm), name = NULL) +
    ggplot2::labs(x = unit, y = "plants seen for the first time",
                  subtitle = "First sightings by period. Plants first seen at the first period are omitted.") +
    rl_theme() + ggplot2::theme(legend.position = "top")
  if (!is.null(d$expected_returns))
    p <- p + ggplot2::geom_line(data = d, ggplot2::aes(.data$period, .data$expected_returns),
                                inherit.aes = FALSE, colour = rl_grey, linewidth = 0.9) +
      ggplot2::geom_point(data = d, ggplot2::aes(.data$period, .data$expected_returns),
                          inherit.aes = FALSE, colour = rl_grey, size = 1.8) +
      ggplot2::labs(subtitle = sprintf("Bars: first sightings by state at first sighting.\nLine: returns expected from plants present but unseen at the first period (survival %.2f, resighting %.2f).",
                                       x$s, x$r))
  p
}

#' Published dormancy, survival and detection for nine *Caladenia* species
#'
#' Posterior means and 95% credible intervals from the multi-state
#' capture-recapture analysis of Tremblay et al. (2009), Table 2: the
#' probability of surviving from one year to the next, the probability of
#' being seen above ground given survival, and the probability of dormancy,
#' which is `p_survival * (1 - p_resight)`. These are the `s` and `r` that
#' [recruit_triage()] uses to compute how many first sightings a record's own
#' dormancy would produce.
#'
#' @format A data frame with 9 rows: `species`, `p_resight`, `p_resight_lo`,
#'   `p_resight_hi`, `p_survival`, `p_survival_lo`, `p_survival_hi`,
#'   `p_dormancy`, `p_dormancy_lo`, `p_dormancy_hi`, `n_plants`, `years`.
#'
#' @section Scale and how to read it:
#'
#' All are annual probabilities. `p_resight` near 1 means a living plant is
#' almost always above ground; near 0.6 means a third of the living plants
#' are unseen in any year. Dormancy of two or more consecutive years was rare
#' in every species (below 3% except *C. macroclavia*, 12%), which is why the
#' expected-returns series in [recruit_triage()] falls off within three
#' periods. The values come from one or a few populations per species and
#' are not the species' constants; use them as a calibration, and prefer a
#' capture-recapture estimate from the record itself when it is long enough
#' to give one.
#'
#' @source Tremblay, R. L., Perez, M.-E., Larcombe, M., Brown, A., Quarmby,
#'   J., Bickerton, D., French, G. and Bould, A. (2009) Dormancy in
#'   *Caladenia*: a Bayesian approach to evaluating latency. *Australian
#'   Journal of Botany* 57: 340-350. \doi{10.1071/BT08163}
#' @examples
#' caladenia_dormancy
"caladenia_dormancy"
