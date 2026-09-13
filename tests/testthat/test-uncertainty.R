# The three routes to an interval on the exceedance. The point estimate of
# every route must equal lag_ceiling() on the same kernels, the intervals
# must behave as intervals, and on the Lepanthes eltoroensis hosts the
# exceedance must stay above 1 in every replicate, which is what the
# companion paper reports.
k <- lag_kernels(4, bin = 6, unit = "mo")
ce <- lag_ceiling(lepanthes_hosts, reproduction = "inflorescences", K = 4, kernels = k)

test_that("the bootstrap point estimate is lag_ceiling() and the interval excludes 1", {
  set.seed(1)
  b <- lag_ceiling_boot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                        K = 4, kernels = k, R = 400)
  expect_s3_class(b, "lag_ceiling_draws")
  expect_equal(unname(b$point["exceedance_gini"]), unname(ce$exceedance["gini"]))
  expect_equal(unname(b$point["exceedance_cv"]),   unname(ce$exceedance["cv"]))
  expect_equal(unname(b$point["ceiling_gini"]),    unname(ce$ceiling["gini"]))
  expect_equal(b$n_units, 23); expect_equal(nrow(b$draws), 400); expect_equal(nrow(b$jack), 23)
  tab <- b$table
  expect_true(all(tab$lower <= tab$estimate & tab$estimate <= tab$upper))
  ex <- tab[tab$quantity == "exceedance Gini", ]
  expect_gt(ex$lower, 1)
  expect_equal(unname(b$p_boot["gini"]), 0)
  expect_output(print(b), "BCa")
  expect_output(print(b), "not a posterior probability")
  df <- as.data.frame(b)
  expect_equal(nrow(df), 6); expect_equal(unique(df$method), "bootstrap")
})

test_that("percentile and BCa intervals both bracket the estimate, and BCa differs", {
  set.seed(2)
  bp <- lag_ceiling_boot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                         K = 4, kernels = k, R = 400, type = "percentile")
  set.seed(2)
  bc <- lag_ceiling_boot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                         K = 4, kernels = k, R = 400, type = "bca")
  expect_equal(bp$draws, bc$draws)                      # same replicates, same seed
  expect_false(isTRUE(all.equal(bp$table$lower, bc$table$lower)))
  expect_true(all(bp$table$lower <= bp$table$estimate & bp$table$estimate <= bp$table$upper))
})

test_that("the Bayesian bootstrap agrees with the ordinary one and gives a probability", {
  set.seed(3)
  bb <- lag_ceiling_bayesboot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                              K = 4, kernels = k, draws = 400)
  expect_equal(unname(bb$point["exceedance_gini"]), unname(ce$exceedance["gini"]))
  expect_equal(unname(bb$prob["gini"]), 1)
  ex <- bb$table[bb$table$quantity == "exceedance Gini", ]
  expect_gt(ex$lower, 1.5); expect_lt(ex$upper, 4)
  expect_output(print(bb), "P\\(exceedance > 1\\)")
  set.seed(3)
  bh <- lag_ceiling_bayesboot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                              K = 4, kernels = k, draws = 400, interval = "hdi")
  expect_equal(bh$draws, bb$draws)
  eh <- bh$table[bh$table$quantity == "exceedance Gini", ]
  expect_lte(eh$upper - eh$lower, ex$upper - ex$lower)   # the HDI is the shortest interval
})

test_that("the interval helpers do what they say", {
  set.seed(4)
  d <- rexp(5000)
  h <- hdi_interval(d, 0.9); e <- percentile_interval(d, 0.9)
  expect_lt(h$upper - h$lower, e$upper - e$lower)
  expect_equal(h$lower, min(d), tolerance = 0.01)         # an exponential's HDI starts at 0
  # BCa on a symmetric, unbiased bootstrap distribution is the percentile interval
  b <- rnorm(20000); j <- rnorm(30)
  bc <- bca_interval(0, b, j, 0.9); pc <- percentile_interval(b, 0.9)
  expect_equal(bc$lower, pc$lower, tolerance = 0.05); expect_equal(bc$upper, pc$upper, tolerance = 0.05)
  # the column-wise indices match the scalar ones
  M <- matrix(rpois(40, 3), 10, 4)
  expect_equal(gini_cols(M), apply(M, 2, rain_gini))
  expect_equal(cv_cols(M), apply(M, 2, rain_cv))
})

test_that("the routes refuse what they cannot do, and say why", {
  expect_error(lag_ceiling_boot(lepanthes_census, reproduction = "inflorescences", K = 4, kernels = k),
               "needs the replicate unit")
  expect_error(lag_ceiling_boot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                                K = 4, kernels = k, level = 90), "strictly between 0 and 1")
  expect_error(lag_ceiling_boot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                                K = 4, kernels = k, R = 1), "at least 2")
  one <- lepanthes_hosts[lepanthes_hosts$host == lepanthes_hosts$host[1], ]
  expect_error(lag_ceiling_bayesboot(one, unit = "host", reproduction = "inflorescences",
                                     K = 4, kernels = k), "only one unit")
  expect_error(lag_ceiling_boot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                                K = 3, kernels = k), "K = 3")
})

test_that("lag_ceiling_stan() needs cmdstanr and refuses rates", {
  if (!requireNamespace("cmdstanr", quietly = TRUE) ||
      is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL))) {
    expect_error(lag_ceiling_stan(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                                  K = 4, kernels = k), "cmdstanr|CmdStan")
    skip("cmdstanr or CmdStan not available")
  }
  d <- lepanthes_hosts; d$rate <- d$recruits / pmax(d$adults, 1)
  expect_error(lag_ceiling_stan(d, unit = "host", reproduction = "inflorescences", recruits = "rate",
                                K = 4, kernels = k), "non-whole values")
  expect_error(lag_ceiling_stan(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                                K = 4, kernels = k, phi_prior = 2), "shape and rate")
  skip_on_cran()
  st <- lag_ceiling_stan(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                         K = 4, kernels = k, chains = 2, iter_warmup = 500, iter_sampling = 300,
                         seed = 1, refresh = 0)
  expect_s3_class(st, "lag_ceiling_draws")
  expect_equal(st$method, "stan")
  expect_equal(nrow(st$draws), 600)
  ex <- st$table[st$table$quantity == "exceedance Gini", ]
  expect_true(ex$lower <= ex$estimate && ex$estimate <= ex$upper)
  expect_true(is.finite(st$diagnostics$max_rhat))
  expect_output(print(st), "Negative-binomial")
})

test_that("few units warn, and a collapsed BCa row is flagged rather than printed as an interval", {
  five <- lepanthes_hosts[lepanthes_hosts$host %in% unique(lepanthes_hosts$host)[1:5], ]
  set.seed(5)
  expect_warning(b <- lag_ceiling_boot(five, unit = "host", reproduction = "inflorescences",
                                       K = 4, kernels = k, R = 300), "Only 5 units")
  tab <- b$table
  bad <- tab$lower > tab$estimate | tab$upper < tab$estimate | tab$upper == tab$lower
  expect_true(all(nzchar(tab$note[bad])))
  set.seed(5)
  expect_silent(lag_ceiling_boot(five, unit = "host", reproduction = "inflorescences",
                                 K = 4, kernels = k, R = 300, type = "percentile"))
})

test_that("ceiling_draws_table() stacks routes and names each probability, rain_forest() draws", {
  set.seed(6)
  b  <- lag_ceiling_boot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                         K = 4, kernels = k, R = 200)
  bb <- lag_ceiling_bayesboot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                              K = 4, kernels = k, draws = 200)
  tab <- ceiling_draws_table(`A` = b, `A` = bb)
  expect_equal(nrow(tab), 4)
  expect_setequal(unique(tab$method), c("cluster bootstrap", "Bayesian bootstrap"))
  expect_true(all(tab$probability_is[tab$method == "cluster bootstrap"] == "one-sided bootstrap p for exceedance <= 1"))
  expect_true(all(tab$probability_is[tab$method == "Bayesian bootstrap"] == "P(exceedance > 1)"))
  expect_equal(nrow(ceiling_draws_table(.list = list(x = b), quantity = "all")), 6)
  expect_error(ceiling_draws_table(), "at least one")
  expect_error(ceiling_draws_table(a = 1), "not lag_ceiling_draws")
  expect_error(ceiling_draws_table(a = b, quantity = "nope"), "quantity must be")
  p <- rain_forest(`A` = b, `A` = bb, index = "both")
  expect_s3_class(p, "ggplot")
  expect_s3_class(rain_forest(.list = list(A = b), order = "name"), "ggplot")
})

test_that("a unit with a flat reproductive record gets NA, not Inf, and many skipped units are summarised", {
  d <- lepanthes_hosts
  h1 <- unique(d$host)[1]
  d$inflorescences[d$host == h1] <- 5                     # perfectly flat: ceiling 0
  ces <- suppressWarnings(lag_ceiling_by(d, unit = "host", reproduction = "inflorescences", K = 4, kernels = k))
  df <- as.data.frame(ces)
  expect_true(is.na(df$exceedance_gini[df$unit == as.character(h1)]))
  expect_match(df$best_kernel[df$unit == as.character(h1)], "undefined")
  expect_false(any(is.infinite(df$exceedance_gini)))
  ces2 <- suppressWarnings(lag_ceiling_by(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                                          K = 4, kernels = k, min_recruits = 8))
  expect_gt(length(attr(ces2, "skipped")), 10)
  expect_output(print(ces2), "of 23 units skipped")
})
