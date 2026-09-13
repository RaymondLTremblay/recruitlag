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
#' rain_series(ce, unit = "census")
#' @export
rain_series <- function(x, unit = "period") {
  rl_check_class(x, "lag_ceiling", "x", "lag_ceiling()")
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
#' rain_strip(`monthly 1994` = c1, `six-monthly 1999-2004` = c6)
#' @export
rain_strip <- function(..., index = c("gini", "cv")) {
  index <- match.arg(index)
  objs <- list(...)
  if (length(objs) == 1 && is.list(objs[[1]]) && !inherits(objs[[1]], "lag_ceiling")) objs <- objs[[1]]
  bad <- !vapply(objs, inherits, logical(1), "lag_ceiling")
  if (any(bad))
    rl_abort(sum(bad), " of the ", length(objs), " things passed to rain_strip() ",
             if (sum(bad) == 1) "is" else "are", " not a lag_ceiling object (",
             rl_list(unique(vapply(objs[bad], function(z) class(z)[1], character(1)))), "). ",
             "Pass lag_ceiling() results, named, or the whole list from lag_ceiling_by().")
  if (is.null(names(objs))) names(objs) <- paste("series", seq_along(objs))
  if (length(objs) > 12)
    rl_abort("rain_strip() was given ", length(objs), " series. One strip per series is ",
             "readable up to about a dozen; beyond that the strips are too thin to read ",
             "and the labels collide. For many units use rain_units(), which puts every ",
             "unit's ceiling against its observed concentration on one pair of axes. ",
             "To keep the strips, pass a subset, for example rain_strip(ces[1:8]).")
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
                                    label = ifelse(is.finite(.data$exceedance),
                                                   sprintf("%.2f (%.1fx)", .data$observed, .data$exceedance),
                                                   sprintf("%.2f (ceiling 0)", .data$observed)),
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
#' `rain_profiles()` draws a set of lag profiles as lollipops.
#' `rain_expected()` applies them to the reproductive record of a
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
#' rain_profiles(pr, bin = 6, unit = "months")
#' ce <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4,
#'                   kernels = lag_kernels(4, bin = 6, unit = "mo"))
#' rain_expected(ce, pr)
#' @export
rain_profiles <- function(profiles, bin = NULL, unit = "", lag0 = 1L) {
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

#' @rdname rain_profiles
#' @param x A `lag_ceiling` object supplying the reproductive record, the
#'   recruit series and the lag convention.
#' @param period Label for the time axis.
#' @export
rain_expected <- function(x, profiles, period = "census") {
  rl_check_class(x, "lag_ceiling", "x", "lag_ceiling()")
  pl <- as_profile_list(profiles)
  if (!all(lengths(pl) == x$K + 1))
    rl_abort("The profiles have ", rl_list(unique(lengths(pl))), " weights, but this ",
             "lag_ceiling was built with K = ", x$K, ", which needs ", x$K + 1,
             " weights in every profile.")
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
  if (!is.list(profiles) || is.null(names(profiles)) || any(names(profiles) == ""))
    rl_abort("`profiles` must be a NAMED list of weight vectors, for example ",
             "list(`evenly spread` = rep(0.25, 4)), or a lag_kernels object. ",
             "The names are what the legend shows.")
  if (length(unique(lengths(profiles))) != 1)
    rl_abort("The profiles have different lengths (", rl_list(unique(lengths(profiles))),
             "). Every profile must cover the same horizon.")
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
  d <- d[is.finite(d$gini), , drop = FALSE]      # an all-zero simulated series has no Gini
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


#' Every unit at once: its ceiling against what it recorded
#'
#' The companion figure to [rain_strip()], for when there are more units than
#' strips can show. Each point is one unit: its own ceiling on the horizontal
#' axis, the concentration of its observed recruitment on the vertical. The
#' diagonal is equality, so a unit above the line recorded recruitment more
#' concentrated than any delay of its own reproduction can predict, and the
#' vertical distance from the line is the size of that gap.
#'
#' @param x A `lag_ceiling_list`, as returned by [lag_ceiling_by()], or a
#'   named list of `lag_ceiling` objects.
#' @param index `"gini"` or `"cv"`.
#' @param label Units to label, as a character vector of names, or `NULL` for
#'   none. Labelling more than a handful defeats the purpose of this figure.
#' @return A ggplot.
#'
#' @section Scale and how to read it:
#'
#' Both axes are concentration indices on the scale of [rain_gini()], so both
#' run from 0 to a maximum set by the number of periods, and the diagonal is
#' the only line that matters: it is the ceiling. The spread along the
#' horizontal axis is the part people do not expect. It shows how much the
#' bound itself varies between units, because the ceiling is a property of
#' each unit's own reproductive record, and a unit whose reproduction is
#' patchy sits far to the right and can never be refuted.
#'
#' A unit at ceiling 0, whose reproduction is perfectly flat, sits on the
#' vertical axis. Any recruitment concentration at all is infinitely above
#' its ceiling, which is correct and is why the exceedance ratio is not used
#' as the vertical axis here.
#'
#' @examples
#' ces <- lag_ceiling_by(lepanthes_hosts, unit = "host",
#'                       reproduction = "inflorescences", K = 4,
#'                       kernels = lag_kernels(4, bin = 6, unit = "mo"))
#' rain_units(ces)
#' @export
rain_units <- function(x, index = c("gini", "cv"), label = NULL) {
  index <- match.arg(index)
  if (!is.list(x) || !length(x)) rl_abort("`x` must be a non-empty list of lag_ceiling objects.")
  bad <- !vapply(x, inherits, logical(1), "lag_ceiling")
  if (any(bad))
    rl_abort(sum(bad), " of the ", length(x), " elements are not lag_ceiling objects. ",
             "Pass the result of lag_ceiling_by().")
  if (is.null(names(x))) names(x) <- paste("unit", seq_along(x))
  d <- do.call(rbind, lapply(names(x), function(nm) {
    o <- x[[nm]]
    data.frame(unit = nm, ceiling = unname(o$theorem[index]),
               observed = unname(o$obs[index]), n = length(o$R),
               stringsAsFactors = FALSE)
  }))
  above <- sum(d$observed > d$ceiling)
  lim <- c(0, max(d$ceiling, d$observed, na.rm = TRUE) * 1.05)
  d$lab <- ifelse(d$unit %in% label, d$unit, NA_character_)
  ggplot2::ggplot(d, ggplot2::aes(.data$ceiling, .data$observed)) +
    ggplot2::geom_polygon(data = data.frame(x = c(lim[1], lim[2], lim[1]),
                                            y = c(lim[1], lim[2], lim[2])),
                          ggplot2::aes(.data$x, .data$y), fill = "#d7d7d7", alpha = 0.55) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = rl_blue, linewidth = 0.8) +
    ggplot2::geom_point(colour = rl_verm, size = 2.4, alpha = 0.75) +
    { if (any(!is.na(d$lab)))
        ggplot2::geom_text(ggplot2::aes(label = .data$lab), na.rm = TRUE,
                           hjust = -0.15, size = 3, colour = rl_grey) } +
    ggplot2::coord_equal(xlim = lim, ylim = lim, expand = FALSE) +
    ggplot2::labs(
      x = sprintf("this unit's ceiling (%s of its own reproduction)",
                  if (index == "gini") "Gini" else "CV"),
      y = sprintf("observed recruitment (%s)", if (index == "gini") "Gini" else "CV"),
      subtitle = sprintf("One point per unit. The shaded region is unreachable by any delay.\n%d of the %d units are in it.",
                         above, nrow(d))) +
    rl_theme()
}


#' Several strip figures instead of one unreadable one
#'
#' [rain_strip()] is readable up to about a dozen series. With more units than
#' that, `rain_strips()` splits them and returns one figure per group, either
#' by a grouping you supply (tree species, transect, site) or in chunks of
#' `size`. A group that is still too large is split further, and the figures
#' are named so the caption can say which is which.
#'
#' @param x A `lag_ceiling_list` from [lag_ceiling_by()], or a named list of
#'   `lag_ceiling` objects.
#' @param group Optional grouping for the units: either a vector as long as
#'   `x`, or a named vector or one-to-one lookup whose names are unit names.
#'   Units whose group is missing are collected as "ungrouped".
#' @param size Largest number of strips in one figure. Groups larger than this
#'   are split into "... (1 of 3)" and so on.
#' @param index `"gini"` or `"cv"`, passed to [rain_strip()].
#' @param order_by `"exceedance"` (default), `"ceiling"`, `"observed"` or
#'   `"name"`: how units are ordered within each figure.
#' @return A named list of ggplots. Print them in a loop, one per figure.
#'
#' @examples
#' ces <- lag_ceiling_by(lepanthes_hosts, unit = "host",
#'                       reproduction = "inflorescences", K = 4,
#'                       kernels = lag_kernels(4, bin = 6, unit = "mo"))
#' figs <- rain_strips(ces, size = 8)
#' names(figs)
#' figs[[1]]
#' @export
rain_strips <- function(x, group = NULL, size = 10,
                        index = c("gini", "cv"),
                        order_by = c("exceedance", "ceiling", "observed", "name")) {
  index <- match.arg(index); order_by <- match.arg(order_by)
  if (!is.list(x) || !length(x)) rl_abort("`x` must be a non-empty list of lag_ceiling objects.")
  bad <- !vapply(x, inherits, logical(1), "lag_ceiling")
  if (any(bad))
    rl_abort(sum(bad), " of the ", length(x), " elements are not lag_ceiling objects. ",
             "Pass the result of lag_ceiling_by().")
  if (is.null(names(x))) names(x) <- paste("unit", seq_along(x))
  if (size < 1) rl_abort("`size` must be at least 1; it is ", size, ".")

  g <- rep("all units", length(x))
  if (!is.null(group)) {
    if (length(group) == length(x) && is.null(names(group))) {
      g <- as.character(group)
    } else if (!is.null(names(group))) {
      g <- unname(as.character(group)[match(names(x), names(group))])
    } else {
      rl_abort("`group` must be either a vector as long as x (", length(x),
               " units), or a named vector whose names are the unit names. ",
               "It has length ", length(group), " and no names.")
    }
    g[is.na(g) | g == ""] <- "ungrouped"
  }

  key <- switch(order_by,
    exceedance = vapply(x, function(o) unname(o$obs[index] / o$theorem[index]), numeric(1)),
    ceiling    = vapply(x, function(o) unname(o$theorem[index]), numeric(1)),
    observed   = vapply(x, function(o) unname(o$obs[index]), numeric(1)),
    name       = seq_along(x))
  if (order_by != "name") key[!is.finite(key)] <- max(key[is.finite(key)], na.rm = TRUE) + 1

  out <- list()
  for (grp in unique(g)) {
    idx <- which(g == grp)
    idx <- idx[order(key[idx], decreasing = order_by != "name")]
    chunks <- split(idx, ceiling(seq_along(idx) / size))
    for (i in seq_along(chunks)) {
      nm <- if (length(chunks) == 1) grp else sprintf("%s (%d of %d)", grp, i, length(chunks))
      out[[nm]] <- rain_strip(x[chunks[[i]]], index = index) +
        ggplot2::labs(title = nm)
    }
  }
  out
}
