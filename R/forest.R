# One table and one figure for any number of lag_ceiling_draws objects, so
# that records and routes can be read side by side: a tidy table for a
# document, and a forest plot of the exceedance intervals against the line
# at 1.

collect_draws <- function(dots, .list, what) {
  x <- c(dots, .list)
  x <- x[!vapply(x, is.null, logical(1))]        # a route that was not run may be passed as NULL
  if (!length(x))
    rl_abort(what, " needs at least one lag_ceiling_draws object, as returned by ",
             "lag_ceiling_boot(), lag_ceiling_bayesboot() or lag_ceiling_stan(). ",
             "Pass them named, for example ", what, "(`L. eltoroensis` = b, `Caladenia` = bb).")
  bad <- !vapply(x, inherits, logical(1), "lag_ceiling_draws")
  if (any(bad))
    rl_abort(sum(bad), " of the ", length(x), " objects are not lag_ceiling_draws objects: ",
             rl_list(vapply(x[bad], function(o) class(o)[1], character(1))), ". ",
             "Every argument must come from lag_ceiling_boot(), lag_ceiling_bayesboot() or ",
             "lag_ceiling_stan().")
  if (is.null(names(x))) names(x) <- rep("", length(x))
  blank <- names(x) == ""
  if (any(blank)) names(x)[blank] <- paste("record", which(blank))
  x
}

method_label <- function(m) c(bootstrap = "cluster bootstrap", `bayesian bootstrap` = "Bayesian bootstrap",
                              stan = "negative-binomial model")[[m]]

#' One table for several records and routes
#'
#' Stacks the interval tables of any number of `lag_ceiling_draws` objects
#' (from [lag_ceiling_boot()], [lag_ceiling_bayesboot()] or
#' [lag_ceiling_stan()]) into one tidy data frame with a `record` column, so
#' that a document can print one table and a reader can compare records and
#' routes down a column. The probability statement each route makes is
#' carried in its own column with its name beside it, because the bootstrap's
#' one-sided p-value and the Bayesian routes' posterior probability are
#' different things and must not share a heading.
#'
#' @param ... `lag_ceiling_draws` objects, named. The names become the
#'   `record` column, so name them by the record or the route as the table
#'   needs: `` `L. eltoroensis` = b `` or `` `L. eltoroensis, bootstrap` = b ``.
#' @param .list Alternatively, a named list of such objects. `NULL` entries,
#'   for a route that was not run, are dropped.
#' @param quantity Which rows to keep: `"exceedance"` (the default), `"all"`,
#'   or any of the six quantity labels (`"observed Gini"`, `"ceiling Gini"`,
#'   `"exceedance Gini"`, `"observed CV"`, `"ceiling CV"`, `"exceedance CV"`).
#'
#' @return A data frame with columns `record`, `method`, `estimand`,
#'   `quantity`, `estimate`, `lower`, `upper`, `level`, `interval`,
#'   `probability`, `probability_is` and `note`. `probability` is filled on
#'   the exceedance rows only.
#'
#' @section Scale and how to read it:
#'
#' `estimate`, `lower` and `upper` are on the scale of the quantity: a
#' concentration index for the observed and ceiling rows, a ratio centred on
#' 1 for the exceedance rows. `interval` says which interval it is (BCa,
#' percentile, equal-tailed, highest-density) and `method` which route made
#' it, and the two together fix how the interval may be read: a confidence
#' interval from the cluster bootstrap, a credible interval from the other
#' two. `probability_is` names the statement in `probability`. Print the
#' table with `knitr::kable(digits = 2)` in a document.
#'
#' @examples
#' set.seed(1)
#' k <- lag_kernels(4, bin = 6, unit = "mo")
#' b  <- lag_ceiling_boot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
#'                        K = 4, kernels = k, R = 200)
#' bb <- lag_ceiling_bayesboot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
#'                             K = 4, kernels = k, draws = 200)
#' ceiling_draws_table(`L. eltoroensis` = b, `L. eltoroensis` = bb)
#' @export
ceiling_draws_table <- function(..., .list = NULL, quantity = "exceedance") {
  x <- collect_draws(list(...), .list, "ceiling_draws_table()")
  labs <- unname(stat_labels)
  keep <- switch(quantity, exceedance = labs[grepl("exceedance", labs)], all = labs,
                 if (quantity %in% labs) quantity else
                   rl_abort('quantity must be "exceedance", "all" or one of: ', rl_list(labs), '. It was "', quantity, '".'))
  out <- do.call(rbind, lapply(seq_along(x), function(i) {
    nm <- names(x)[i]; o <- x[[i]]; tab <- o$table
    is_ex <- grepl("exceedance", tab$quantity)
    idx <- ifelse(grepl("Gini", tab$quantity), "gini", "cv")
    prob <- rep(NA_real_, nrow(tab)); pis <- rep("", nrow(tab))
    if (o$method == "bootstrap") {
      prob[is_ex] <- o$p_boot[idx[is_ex]]; pis[is_ex] <- "one-sided bootstrap p for exceedance <= 1"
    } else {
      prob[is_ex] <- o$prob[idx[is_ex]]; pis[is_ex] <- "P(exceedance > 1)"
    }
    data.frame(record = nm, method = method_label(o$method), estimand = o$estimand,
               quantity = tab$quantity, estimate = tab$estimate, lower = tab$lower,
               upper = tab$upper, level = o$level, interval = interval_label(o$type),
               probability = prob, probability_is = pis, note = tab$note,
               stringsAsFactors = FALSE)
  }))
  out <- out[out$quantity %in% keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' The exceedance and its interval, for several records and routes on one axis
#'
#' A forest plot: one row per record, one point and interval per route, and
#' a vertical line at 1, which is the bound. Because the exceedance is a
#' ratio the axis is logarithmic, so that a record twice above the bound and
#' one twice below it sit at equal distances from the line.
#'
#' @inheritParams ceiling_draws_table
#' @param index `"gini"` (the default), `"cv"` or `"both"` (two panels).
#' @param order `"estimate"` (the default, largest exceedance at the top) or
#'   `"name"` (the order the records were given).
#' @return A ggplot.
#'
#' @section Scale and how to read it:
#'
#' The horizontal axis is the exceedance, observed concentration over the
#' ceiling, on a log scale with 1 marked. An interval wholly to the right
#' of 1 is a record more concentrated than any delay of its reproductive
#' series predicts, under that route's reading of "interval"; one that
#' crosses 1 is a record on which the check cannot decide. Routes are
#' distinguished by shape and colour and are read as their help pages say:
#' a confidence interval for the cluster bootstrap, a credible interval for
#' the Bayesian bootstrap and the negative-binomial model, and the last of
#' these estimates the expected series rather than the counts, so it usually
#' sits to the left of the other two. A row whose routes disagree is worth a
#' sentence; a row whose routes agree and cross 1 is a record on which the
#' reproductive series is itself episodic, or the recruits too few, and the
#' figure cannot tell which. The rows are ordered by the first route's
#' estimate unless `order = "name"`.
#'
#' @examples
#' set.seed(1)
#' k <- lag_kernels(4, bin = 6, unit = "mo")
#' b  <- lag_ceiling_boot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
#'                        K = 4, kernels = k, R = 200)
#' bb <- lag_ceiling_bayesboot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
#'                             K = 4, kernels = k, draws = 200)
#' rain_forest(`L. eltoroensis` = b, `L. eltoroensis` = bb, index = "both")
#' @export
rain_forest <- function(..., .list = NULL, index = c("gini", "cv", "both"),
                        order = c("estimate", "name")) {
  index <- match.arg(index); order <- match.arg(order)
  d <- ceiling_draws_table(..., .list = .list, quantity = "exceedance")
  d$index <- ifelse(grepl("Gini", d$quantity), "Gini", "CV")
  if (index != "both") d <- d[d$index == (if (index == "gini") "Gini" else "CV"), , drop = FALSE]
  d <- d[is.finite(d$estimate), , drop = FALSE]
  if (!nrow(d)) rl_abort("No finite exceedance to draw.")
  d$method <- factor(d$method, levels = c("cluster bootstrap", "Bayesian bootstrap", "negative-binomial model"))
  recs <- unique(d$record)
  if (order == "estimate") {
    first <- vapply(recs, function(r) {
      dd <- d[d$record == r & d$index == d$index[1], ]
      dd$estimate[order(as.integer(dd$method))][1]
    }, numeric(1))
    recs <- recs[order(first)]
  } else recs <- rev(recs)
  d$record <- factor(d$record, levels = recs)
  d$index <- factor(d$index, levels = c("Gini", "CV"))
  lo <- min(c(d$lower, d$estimate, 1), na.rm = TRUE); hi <- max(c(d$upper, d$estimate, 1), na.rm = TRUE)
  brk <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32); brk <- brk[brk >= lo / 2 & brk <= hi * 2]
  ggplot2::ggplot(d, ggplot2::aes(.data$estimate, .data$record, colour = .data$method, shape = .data$method)) +
    ggplot2::geom_vline(xintercept = 1, colour = rl_grey, linewidth = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
                           width = 0.25, position = ggplot2::position_dodge(width = 0.6),
                           linewidth = 0.6, na.rm = TRUE, orientation = "y") +
    ggplot2::geom_point(size = 2.6, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_x_log10(breaks = brk) +
    ggplot2::scale_colour_manual(values = c(`cluster bootstrap` = rl_blue, `Bayesian bootstrap` = rl_verm,
                                            `negative-binomial model` = "#009E73"), name = NULL, drop = TRUE) +
    ggplot2::scale_shape_manual(values = c(`cluster bootstrap` = 16, `Bayesian bootstrap` = 17,
                                           `negative-binomial model` = 15), name = NULL, drop = TRUE) +
    { if (index == "both") ggplot2::facet_wrap(~ index, scales = "free_x") } +
    ggplot2::labs(x = "exceedance: observed concentration / ceiling (log scale)", y = NULL,
                  subtitle = sprintf("%d%% intervals. Right of 1: more concentrated than any delay predicts.",
                                     round(100 * d$level[1]))) +
    rl_theme() + ggplot2::theme(legend.position = "top")
}
