# Regression tests for the bugs found in the 2026-10-05 pre-release review
# (REVIEW.md, items B1 to B11). Each test reproduces the failure that was
# observed and checks the behaviour that replaced it.

k4 <- lag_kernels(4, bin = 6, unit = "mo")

test_that("B1: expected_recruits() pairs units by name, whatever the row order of X", {
  m <- host_matrices(lepanthes_hosts, unit = "host", reproduction = "inflorescences")
  a <- expected_recruits(m$R, m$X, w = flat_kernel(4))
  b <- expected_recruits(m$R, m$X[rev(rownames(m$X)), ], w = flat_kernel(4))
  expect_equal(a, b)
})

test_that("B2 and B3: with missing = 'drop' the observed series is cut to the same periods as the ceiling", {
  ce <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4,
                    kernels = k4, missing = "drop")
  expect_equal(ce$t_R, 6:12)                 # periods 3 to 5 have windows before the record
  expect_equal(length(ce$R), 7)
  expect_equal(ce$n_dropped, 3L)
  expect_equal(unname(ce$obs["gini"]), rain_gini(lepanthes_census$recruits[6:12]))
  expect_output(print(ce), "3 scored period\\(s\\) left out")
  cal <- ceiling_calibration(ce, nsim = 200)
  expect_false(anyNA(cal$p_gini))
  expect_error(lag_ceiling(1:6, R = c(2, 0, 3), t_R = 3:5, K = 4, kernels = k4, missing = "drop"),
               "no scored period has a window")
})

test_that("B4: a flat reproductive record gives NA exceedance, never Inf or a negative Gini", {
  expect_equal(rain_gini(rep(12.3456, 12)), 0)
  expect_gte(rain_gini(rep(0.37, 12)), 0)
  ce <- lag_ceiling(rep(0.37, 12), R = c(0, 0, 5, 0, 0, 0, 9, 0, 0, 0), t_R = 3:12,
                    K = 2, kernels = lag_kernels(2))
  expect_true(ce$flat)
  expect_true(all(is.na(ce$exceedance)))
  expect_output(print(ce), "undefined")
  ce2 <- lag_ceiling(rep(12.3456, 12), R = c(0, 0, 5, 0, 0, 0, 9, 0, 0, 0), t_R = 3:12,
                     K = 2, kernels = lag_kernels(2))
  expect_true(all(is.na(ce2$exceedance)))
  # the per-unit table and rain_indices() carry it through as NA
  d <- lepanthes_hosts; h1 <- unique(d$host)[1]
  d$inflorescences[d$host == h1] <- 5
  df <- as.data.frame(suppressWarnings(lag_ceiling_by(d, unit = "host", reproduction = "inflorescences",
                                                      K = 4, kernels = k4)))
  expect_true(is.na(df$exceedance_gini[df$unit == as.character(h1)]))
})

test_that("B5: min_periods counts scored periods, not rows", {
  d <- data.frame(unit = "a", period = 1:8, reproduction = c(3, 4, 5, 4, 3, 5, 4, 3),
                  recruits = c(NA, NA, NA, NA, NA, NA, 2, 3))
  expect_error(lag_ceiling_by(d, K = 2, kernels = lag_kernels(2), min_periods = 5),
               "2 scored periods, fewer than min_periods = 5")
  ces <- lag_ceiling_by(d, K = 2, kernels = lag_kernels(2), min_periods = 2)
  expect_equal(length(ces), 1)
})

test_that("B6: convolve_lag() refuses a record too short for the window, and validates t_R and lag0", {
  expect_error(convolve_lag(1:3, rep(1, 5)), "no period to compute")
  expect_error(convolve_lag(1:10, rep(1, 3), t_R = c(2, 11)), "outside the reproductive record")
  expect_error(convolve_lag(1:10, rep(1, 3), t_R = 4.5), "whole numbers")
  expect_error(convolve_lag(1:10, rep(1, 3), lag0 = -1), "lag0")
  expect_equal(convolve_lag(1:10, c(0, 0, 1)), 1:7)   # lag0 = 1, K = 2: periods 4 to 10 read X[t - 3]
})

test_that("B7: long data with a misnamed column is refused, not read as a wide matrix", {
  d <- data.frame(Plant = rep(c("a", "b", "c"), each = 3), period = rep(2001:2003, 3),
                  state = c("V", "F", "V", "V", "V", "F", "F", "V", "V"))
  expect_error(recruit_triage(d), 'plant = "plant"  was not found')
  tr <- recruit_triage(d, plant = "Plant")
  expect_equal(tr$n_plants, 3)
})

test_that("B8: each probability uses the finite replicates of its own statistic, and printing never crashes", {
  # recruits scored at a single period: no replicate can compute either index,
  # so both probabilities are NA and the print method says so instead of failing
  d <- lepanthes_hosts
  d$recruits[d$period != 5] <- NA
  set.seed(3)
  b <- suppressWarnings(lag_ceiling_boot(d, unit = "host", reproduction = "inflorescences",
                                         K = 2, kernels = lag_kernels(2), R = 50))
  expect_true(is.na(b$p_boot[["cv"]]))
  expect_true(is.na(b$p_boot[["gini"]]))
  expect_output(print(b), "no replicate could compute it")
  expect_equal(unname(b$n_used["cv"]), 0L)
  # the rule itself: a replicate with a finite Gini but NA CV still counts for the Gini
  expect_equal(recruitlag:::tail_prob(c(2, 3, NA, Inf, 0.5), upper = FALSE), 1 / 3)
  expect_equal(recruitlag:::tail_prob(c(NA, NaN), upper = TRUE), NA_real_)
  bb <- suppressWarnings(lag_ceiling_bayesboot(d, unit = "host", reproduction = "inflorescences",
                                               K = 2, kernels = lag_kernels(2), draws = 50))
  expect_output(print(bb), "no draw could compute it")
})

test_that("B9: rain_units() counts only the units that can be placed", {
  d <- lepanthes_hosts
  ces <- suppressWarnings(lag_ceiling_by(d, unit = "host", reproduction = "inflorescences",
                                         K = 4, kernels = k4, min_recruits = 0))
  ces[[1]]$obs[] <- NA_real_                   # a unit with no recruits has NA observed indices
  p <- rain_units(ces)
  expect_s3_class(p, "ggplot")
  expect_false(grepl("NA of", p$labels$subtitle))
  expect_match(p$labels$subtitle, "units with recruits")
})

test_that("B11: host_matrices() and series_from() say when missing recruits are read as zero", {
  d <- lepanthes_hosts
  d$recruits[d$host == unique(d$host)[1] & d$period == 5] <- NA
  expect_warning(host_matrices(d, unit = "host", reproduction = "inflorescences"),
                 "read as zero recruits")
  expect_warning(series_from(d, reproduction = "inflorescences"), "some rows have NA")
  expect_silent(host_matrices(lepanthes_hosts, unit = "host", reproduction = "inflorescences"))
})

test_that("S1: the all-zero diagnosis is a classed condition that rain_indices() recognises", {
  d <- lepanthes_census; d$zero <- 0
  expect_error(lag_ceiling(d, reproduction = "zero", K = 4, kernels = k4),
               class = "recruitlag_zero_lagged")
  tab <- suppressWarnings(rain_indices(d, c("inflorescences", "zero"), K = 4, kernels = k4))
  expect_match(tab$best_kernel[tab$index == "zero"], "zero throughout")
  expect_false(is.na(tab$exceedance_gini[tab$index == "inflorescences"]))
  # a data-level error is re-raised, not swallowed into the table
  expect_error(rain_indices(d, "inflorescences", period = "nope", K = 4, kernels = k4), "not in the data")
})

test_that("S3: lag_ceiling_by() reports the per-unit warnings once, with the units", {
  d <- lepanthes_hosts
  d$inflorescences[d$host == unique(d$host)[2] & d$period == 4] <- NA
  w <- capture_warnings(lag_ceiling_by(d, unit = "host", reproduction = "inflorescences", K = 4, kernels = k4))
  expect_true(any(grepl("For unit\\(s\\) .*missing values, and they are read as zero", w)))
  expect_equal(sum(grepl("For unit", w)), 1)
})

test_that("S4 and S5: degenerate inputs get a sentence, not an R error", {
  ce <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4, kernels = k4)
  ce0 <- ce; ce0$R[] <- 0
  expect_error(ceiling_calibration(ce0, nsim = 10), "Every recruit count is zero")
  m <- host_matrices(lepanthes_hosts, unit = "host", reproduction = "inflorescences")
  expect_error(phi_moment(m$R * 0, m$R), "No unit has any recruits")
  expect_error(host_lag_test(m$R, m$X, K = 4, kernels = k4, nsim = 10, phi = -1), "must be positive")
  expect_error(expected_recruits(m$R, m$X, w = c(1, -1, 0, 0, 0)), "non-negative")
  expect_error(lag_kernels(4, rho = -0.5), "rho must be positive")
  expect_error(lag_kernels(4, extra = list(z = c(0, 0, 0, 0, 0))), "sum to zero")
})

test_that("S6: rain_forest() refuses mixed levels and keeps axis ticks for extreme exceedances", {
  set.seed(1)
  b1 <- lag_ceiling_boot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                         K = 4, kernels = k4, R = 30, level = 0.9)
  b2 <- lag_ceiling_boot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                         K = 4, kernels = k4, R = 30, level = 0.5)
  expect_error(rain_forest(a = b1, b = b2), "different interval levels")
  b3 <- b1; b3$table$estimate <- b3$table$estimate * 1000
  b3$table$lower <- b3$table$lower * 1000; b3$table$upper <- b3$table$upper * 1000
  expect_s3_class(rain_forest(a = b3), "ggplot")
})

test_that("S8: species = and s/r together are refused; rain_intervals() does not vouch for undated gaps", {
  d <- data.frame(plant = rep(c("a", "b"), each = 3), period = rep(2001:2003, 2),
                  state = c("V", "F", "V", "V", "V", "F"))
  expect_error(recruit_triage(d, species = "C. valida", s = 0.9), "not both")
  cens <- data.frame(period = 1:4, date = as.Date(c("2001-01-01", "2001-07-01", NA, "2002-07-01")))
  ri <- rain_intervals(cens, date = "date", period = "period")
  expect_output(print(ri), "cannot be checked")
})

test_that("figures and remaining exports run (smoke tests)", {
  ce <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4, kernels = k4)
  expect_s3_class(rain_series(ce), "ggplot")
  expect_s3_class(rain_strip(a = ce), "ggplot")
  expect_s3_class(rain_profiles(list(flat = flat_kernel(4)), bin = 6, unit = "mo"), "ggplot")
  expect_s3_class(rain_expected(ce, list(flat = flat_kernel(4))), "ggplot")
  ces <- suppressWarnings(lag_ceiling_by(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                                         K = 4, kernels = k4))
  expect_s3_class(rain_strips(ces)[[1]], "ggplot")
  expect_s3_class(rain_units(ces), "ggplot")
  expect_equal(sum(flat_kernel(4)), 1); expect_length(flat_kernel(4), 5)
  sf <- series_from(lepanthes_census, reproduction = "inflorescences")
  expect_equal(sf$t_R, 3:12)
  cal <- ceiling_calibration(ce, nsim = 100)
  expect_s3_class(plot(cal), "ggplot")
  set.seed(2)
  ht <- host_lag_test(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                      K = 4, nsim = 50, kernels = lag_kernels(4, families = "delay", bin = 6, unit = "mo"))
  expect_s3_class(plot(ht), "ggplot")
  expect_true(all(ht$verdict$max_p >= 0 & ht$verdict$max_p <= 1))
})
