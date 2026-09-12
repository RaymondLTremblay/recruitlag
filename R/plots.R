# Colours are semantic and used for nothing else: blue is reproduction and
# anything derived from it (the ceiling is reproduction's own concentration),
# vermillion is recruitment as observed, grey marks regions of an axis.
rl_blue <- "#0072B2"; rl_verm <- "#D55E00"; rl_grey <- "#4d4d4d"
rl_profile_cols <- c("#000000", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#F0E442")
rl_theme <- function(base_size = 11)
  ggplot2::theme_minimal(base_size = base_size) +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))

#' Plot the two series on a common time axis
#'
#' Reproduction (blue) and recruitment (vermillion) as bars, one panel each,
#' with the concentration of each printed in the panel title.
#' @param x A `lag_ceiling` object.
#' @param unit Label for the time axis.
#' @return A ggplot.
#' @examples
#' ce <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4,
#'                   kernels = lag_kernels(4, bin = 6, unit = "mo"))
#' plot_series(ce, unit = "census")
#' @export
plot_series <- function(x, unit = "period") {
  stopifnot(inherits(x, "lag_ceiling"))
  d <- rbind(
    data.frame(series = sprintf("reproduction, Gini %.2f (the ceiling)", x$theorem["gini"]),
               t = seq_along(x$X), y = x$X, col = rl_blue),
    data.frame(series = sprintf("recruitment, Gini %.2f", x$obs["gini"]),
               t = x$t_R, y = x$R, col = rl_verm))
  d$series <- factor(d$series, levels = unique(d$series))
  ggplot2::ggplot(d, ggplot2::aes(.data$t, .data$y, fill = .data$series)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_wrap(~ series, ncol = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(rl_blue, rl_verm), guide = "none") +
    ggplot2::labs(x = unit, y = "count") +
    rl_theme() +
    ggplot2::theme(strip.text = ggplot2::element_text(hjust = 0))
}

#' The ceiling strip: observed recruitment against the bound
#'
#' Each strip is the concentration axis from 0 to 1. The blue tick is the
#' ceiling (reproduction's own concentration); the region beyond it is
#' unreachable by any delay; the vermillion point is observed recruitment.
#' @param ... One or more `lag_ceiling` objects, named to label the strips.
#' @param index `"gini"` or `"cv"` (for `"cv"` the axis runs to the largest
#'   value shown).
#' @return A ggplot.
#' @examples
#' k6 <- lag_kernels(4, bin = 6, unit = "mo"); k1 <- lag_kernels(6, bin = 1, unit = "mo")
#' c6 <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4, kernels = k6)
#' c1 <- lag_ceiling(lepanthes_monthly, reproduction = "reproductive_adults", K = 6, lag0 = 0,
#'                   kernels = k1)
#' plot_ceiling_strip(`monthly 1994` = c1, `six-monthly 1999-2004` = c6)
#' @export
plot_ceiling_strip <- function(..., index = c("gini", "cv")) {
  index <- match.arg(index)
  objs <- list(...)
  if (length(objs) == 1 && is.list(objs[[1]]) && !inherits(objs[[1]], "lag_ceiling")) objs <- objs[[1]]
  stopifnot(all(vapply(objs, inherits, logical(1), "lag_ceiling")))
  if (is.null(names(objs))) names(objs) <- paste("series", seq_along(objs))
  d <- do.call(rbind, lapply(names(objs), function(nm) {
    o <- objs[[nm]]
    data.frame(strip = nm, ceiling = unname(o$theorem[index]), observed = unname(o$obs[index]),
               exceedance = unname(o$obs[index] / o$theorem[index]))
  }))
  d$strip <- factor(d$strip, levels = rev(names(objs)))
  xmax <- if (index == "gini") 1 else max(1, d$observed, d$ceiling) * 1.05
  ggplot2::ggplot(d) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = .data$ceiling,
                                    ymin = as.numeric(.data$strip) - 0.3, ymax = as.numeric(.data$strip) + 0.3),
                       fill = "#efefef") +
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$ceiling, xmax = xmax,
                                    ymin = as.numeric(.data$strip) - 0.3, ymax = as.numeric(.data$strip) + 0.3),
                       fill = "#d7d7d7") +
    ggplot2::geom_segment(ggplot2::aes(x = .data$ceiling, xend = .data$ceiling,
                                       y = as.numeric(.data$strip) - 0.42, yend = as.numeric(.data$strip) + 0.42),
                          colour = rl_blue, linewidth = 1.6) +
    ggplot2::geom_text(ggplot2::aes(x = .data$ceiling, y = as.numeric(.data$strip) + 0.40,
                                    label = sprintf("ceiling %.2f", .data$ceiling),
                                    hjust = ifelse(.data$ceiling < 0.15 * xmax, 0, 0.5)), colour = rl_blue, size = 3.2) +
    ggplot2::geom_point(ggplot2::aes(x = .data$observed, y = as.numeric(.data$strip)),
                        colour = rl_verm, size = 4) +
    # The observed label sits beside its point rather than below it. Below, it
    # lands on the ceiling label of the next strip as soon as there are more
    # than two strips, and a figure of one strip per population or per species
    # is the main use.
    ggplot2::geom_text(ggplot2::aes(x = .data$observed + ifelse(.data$observed > 0.8 * xmax, -1, 1) * 0.02 * xmax,
                                    y = as.numeric(.data$strip),
                                    label = sprintf("%.2f (%.1fx)", .data$observed, .data$exceedance),
                                    hjust = ifelse(.data$observed > 0.8 * xmax, 1, 0)),
                       colour = rl_verm, size = 3.2) +
    ggplot2::scale_y_continuous(breaks = seq_along(levels(d$strip)), labels = levels(d$strip),
                                limits = c(0.4, length(objs) + 0.7)) +
    ggplot2::scale_x_continuous(limits = c(0, xmax), expand = ggplot2::expansion(0)) +
    ggplot2::labs(x = if (index == "gini") "concentration in time (Gini coefficient)\n0 = the same count every period, 1 = the whole total in one period"
                  else "concentration in time (coefficient of variation)",
                  y = NULL,
                  subtitle = paste("Each strip is a number line. The light region is reachable by a delay; the dark region is not.",
                                   "Beside each point, its concentration and how many times the ceiling that is.", sep = "\n")) +
    rl_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 10))
}

#' Lag profiles, and what they predict on the real record
#'
#' `plot_profiles()` draws a set of lag profiles as lollipops.
#' `plot_expected()` applies them to the reproductive record of a
#' `lag_ceiling` object, scales each expected series to the observed total,
#' and draws them over the observed recruits. The point of the pair is that
#' profiles which could not differ more give expected series that all but
#' coincide when reproduction is near-constant, and none resembles an
#' episodic record.
#'
#' @param profiles A named list of weight vectors, or a `lag_kernels` object.
#' @param bin,unit Bin width and time unit for the lag axis labels; with
#'   `bin = NULL` bins are labelled by index from the first lagged bin.
#' @param lag0 Offset of the first bin, used for the axis labels.
#' @return A ggplot.
#' @examples
#' pr <- list(`one interval only (projection-matrix null)` = c(1, 0, 0, 0, 0),
#'            `concentrated early` = c(.7, .2, .07, .03, 0),
#'            `evenly spread` = rep(.2, 5),
#'            `concentrated late` = c(0, .03, .07, .2, .7))
#' plot_profiles(pr, bin = 6, unit = "months")
#' ce <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4,
#'                   kernels = lag_kernels(4, bin = 6, unit = "mo"))
#' plot_expected(ce, pr)
#' @export
plot_profiles <- function(profiles, bin = NULL, unit = "", lag0 = 1L) {
  pl <- as_profile_list(profiles)
  K <- length(pl[[1]]) - 1L
  d <- do.call(rbind, lapply(names(pl), function(nm)
    data.frame(profile = nm, k = 0:K, w = pl[[nm]] / sum(pl[[nm]]))))
  d$profile <- factor(d$profile, levels = names(pl))
  npr <- length(pl); off <- (seq_len(npr) - (npr + 1) / 2) * (0.7 / npr)
  d$x <- d$k + off[as.integer(d$profile)]
  labs_k <- if (is.null(bin)) as.character(0:K + lag0) else sprintf("%g", bin * (0:K + lag0))
  ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$w, colour = .data$profile)) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$x, yend = 0), linewidth = 1.2) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_x_continuous(breaks = 0:K, labels = labs_k) +
    ggplot2::scale_colour_manual(values = rl_profile_cols[seq_len(npr)], name = NULL) +
    ggplot2::labs(x = if (is.null(bin)) "lag (bins before the recruit census)" else sprintf("%s before the recruit census (lag bin)", unit),
                  y = "weight on past reproduction") +
    rl_theme() +
    ggplot2::theme(legend.position = "top", legend.direction = "vertical",
                   panel.grid.minor = ggplot2::element_blank())
}

#' @rdname plot_profiles
#' @param x A `lag_ceiling` object supplying the reproductive record, the
#'   recruit series and the lag convention.
#' @param period Label for the time axis.
#' @export
plot_expected <- function(x, profiles, period = "census") {
  stopifnot(inherits(x, "lag_ceiling"))
  pl <- as_profile_list(profiles)
  stopifnot(all(lengths(pl) == x$K + 1))
  tot <- sum(x$R)
  d <- do.call(rbind, lapply(names(pl), function(nm) {
    e <- convolve_lag(x$X, pl[[nm]], x$t_R, x$lag0, x$missing)
    data.frame(profile = nm, t = x$t_R, expected = e * tot / sum(e))
  }))
  d$profile <- factor(d$profile, levels = names(pl))
  obs <- data.frame(t = x$t_R, R = x$R)
  ggplot2::ggplot() +
    ggplot2::geom_col(data = obs, ggplot2::aes(.data$t, .data$R), fill = rl_verm, alpha = 0.85, width = 0.6) +
    ggplot2::geom_line(data = d, ggplot2::aes(.data$t, .data$expected, colour = .data$profile), linewidth = 0.9) +
    ggplot2::geom_point(data = d, ggplot2::aes(.data$t, .data$expected, colour = .data$profile), size = 1.6) +
    ggplot2::scale_colour_manual(values = rl_profile_cols[seq_along(pl)], name = NULL) +
    ggplot2::scale_x_continuous(breaks = x$t_R) +
    ggplot2::labs(x = period, y = "recruits per period",
                  subtitle = "Bars: observed recruits. Lines: expected recruits under each lag profile,\nscaled to the observed total.") +
    rl_theme() +
    ggplot2::theme(legend.position = "top", legend.direction = "vertical",
                   panel.grid.minor = ggplot2::element_blank())
}

as_profile_list <- function(profiles) {
  if (inherits(profiles, "lag_kernels"))
    profiles <- stats::setNames(lapply(profiles, function(k) k$w), vapply(profiles, function(k) k$label, character(1)))
  stopifnot(is.list(profiles), !is.null(names(profiles)), length(unique(lengths(profiles))) == 1)
  profiles
}

#' Plot a calibration or a host-level test
#'
#' For a `ceiling_calibration`, the distribution of the simulated Gini under
#' each null with the observed value and the ceiling marked. For a
#' `host_lag_test`, the probability of a record at least as extreme as the
#' observed one, by kernel family, on a log scale.
#' @param x The object.
#' @param statistic For `host_lag_test`: which probability to plot.
#' @param ... Ignored.
#' @return A ggplot.
#'
#' @section Scale and how to read it:
#'
#' For a `host_lag_test` the vertical axis is a probability on a log scale,
#' floored at \eqn{1/\mathrm{nsim}}{1/nsim} because no simulation can resolve below
#' that. Points sitting on the floor mean "smaller than
#' \eqn{1/\mathrm{nsim}}{1/nsim}", not zero. Two horizontal lines are drawn: the
#' dashed line is the flat reference, which is the comparison that matters,
#' and the dotted line at 0.05 is a conventional reference mark and not a
#' decision rule (Wasserstein and Lazar 2016). Read the plot as the spread of
#' the whole kernel family, since the verdict is the highest point on it.
#'
#' For a `ceiling_calibration` the axis is the Gini of simulated count
#' series, with the observed value and the ceiling on the mean marked. The
#' distance between those two marks is the question the plot answers: how
#' much of the gap between the ceiling and the observation is explained by
#' counting alone.
#'
#' @references
#' Wasserstein, R. L. and Lazar, N. A. (2016) The ASA statement on p-values:
#' context, process, and purpose. *The American Statistician* 70: 129-133.
#' \doi{10.1080/00031305.2016.1154108}
#' @examples
#' ce <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4,
#'                   kernels = lag_kernels(4, bin = 6, unit = "mo"))
#' plot(ceiling_calibration(ce, nsim = 2000))
#'
#' ht <- host_lag_test(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
#'                     K = 4, nsim = 500, kernels = lag_kernels(4, bin = 6, unit = "mo"))
#' plot(ht, statistic = "silent")
#' @export
plot.ceiling_calibration <- function(x, ...) {
  sims <- attr(x, "sims"); obs <- attr(x, "obs"); ce <- attr(x, "ceiling")
  d <- do.call(rbind, lapply(names(sims), function(nm) data.frame(null = nm, gini = sims[[nm]]$gini)))
  v <- data.frame(what = c("observed recruitment", "ceiling on the mean"),
                  x = c(unname(obs["gini"]), unname(ce["gini"])))
  ggplot2::ggplot(d, ggplot2::aes(.data$gini, fill = .data$null)) +
    ggplot2::geom_histogram(bins = 60, alpha = 0.55, position = "identity") +
    ggplot2::geom_vline(data = v, ggplot2::aes(xintercept = .data$x, colour = .data$what), linewidth = 0.8) +
    ggplot2::scale_fill_manual(values = c("#56B4E9", "#009E73"), name = NULL) +
    ggplot2::scale_colour_manual(values = c(rl_blue, rl_verm), name = NULL) +
    ggplot2::labs(x = "Gini coefficient of the simulated count series", y = NULL,
                  subtitle = "What a mean that obeys the ceiling looks like once it is sampled as counts") +
    rl_theme() + ggplot2::theme(legend.position = "top")
}

#' @rdname plot.ceiling_calibration
#' @export
plot.host_lag_test <- function(x, statistic = c("silent", "gini", "cv", "max"), ...) {
  statistic <- match.arg(statistic)
  col <- paste0("p_", statistic)
  d <- x$table[x$table$phi_source == "moment", ]
  floor <- 1 / x$nsim
  d$p <- pmax(d[[col]], floor)
  d$family <- factor(d$family, levels = unique(d$family))
  flat <- d$p[d$family == "flat"]
  ggplot2::ggplot(d, ggplot2::aes(.data$family, .data$p, colour = .data$family)) +
    { if (length(flat)) ggplot2::geom_hline(yintercept = flat, linetype = 2, colour = "grey45") } +
    ggplot2::geom_hline(yintercept = 0.05, linetype = 3, colour = "grey60") +
    ggplot2::geom_jitter(width = 0.18, height = 0, alpha = 0.7, size = 1.8) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = rep(rl_profile_cols, length.out = nlevels(d$family)), guide = "none") +
    ggplot2::labs(x = NULL,
                  y = sprintf("P(record at least as extreme as observed: %s), log scale", statistic),
                  subtitle = "Dashed: flat reference. Dotted: 0.05. Values at the floor were never reached in simulation.") +
    rl_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
}
