# Reproduces output/phase9a of the companion paper (concentration.csv,
# ceiling_verdict.csv). The ceiling is attained by a pure delay, which is in
# the named families, so no random draws are needed to reproduce it exactly.
ref_conc <- read.csv(test_path("reference", "concentration.csv"), check.names = FALSE)
ref_verd <- read.csv(test_path("reference", "ceiling_verdict.csv"))

c_m <- lag_ceiling(lepanthes_monthly, reproduction = "reproductive_adults", K = 6, lag0 = 0,
                   kernels = lag_kernels(6, bin = 1, unit = "mo"))
c_6 <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4, lag0 = 1,
                   kernels = lag_kernels(4, bin = 6, unit = "mo"))

test_that("the data-frame and vector interfaces agree", {
  d <- lepanthes_census
  c_v <- lag_ceiling(d$inflorescences, d$recruits[3:12], t_R = 3:12, K = 4,
                     kernels = lag_kernels(4, bin = 6, unit = "mo"))
  expect_equal(c_v$obs, c_6$obs); expect_equal(c_v$ceiling, c_6$ceiling)
  expect_equal(c_6$t_R, 3:12); expect_equal(c_m$t_R, 2:23)
  # a long unit-by-period frame is summed to the system series
  c_h <- lag_ceiling(lepanthes_hosts, reproduction = "inflorescences", K = 4,
                     kernels = lag_kernels(4, bin = 6, unit = "mo"))
  expect_equal(c_h$obs, c_6$obs); expect_equal(c_h$theorem, c_6$theorem)
})

test_that("concentration indices match the paper", {
  g <- function(ds, ser, col) ref_conc[[col]][ref_conc$dataset == ds & ref_conc$series == ser]
  expect_equal(unname(c_m$theorem["gini"]), g("monthly (1994)", "reproduction", "Gini"))
  expect_equal(unname(c_m$obs["gini"]),     g("monthly (1994)", "recruitment", "Gini"))
  expect_equal(unname(c_m$obs["cv"]),       g("monthly (1994)", "recruitment", "CV"))
  expect_equal(unname(c_6$theorem["gini"]), g("six-monthly (1999-2004)", "reproduction", "Gini"))
  expect_equal(unname(c_6$obs["gini"]),     g("six-monthly (1999-2004)", "recruitment", "Gini"))
  expect_equal(unname(c_6$theorem["cv"]),   g("six-monthly (1999-2004)", "reproduction", "CV"))
})

test_that("the searched ceiling and exceedance match the paper", {
  v <- function(ds, col) ref_verd[[col]][ref_verd$dataset == ds]
  expect_equal(unname(c_m$ceiling["gini"]), v("monthly (1994)", "ceiling_gini"))
  expect_equal(unname(c_m$ceiling["cv"]),   v("monthly (1994)", "ceiling_cv"))
  expect_equal(unname(c_6$ceiling["gini"]), v("six-monthly (1999-2004)", "ceiling_gini"))
  expect_equal(unname(c_6$ceiling["cv"]),   v("six-monthly (1999-2004)", "ceiling_cv"))
  expect_equal(unname(c_m$exceedance["gini"]), v("monthly (1994)", "gini_exceedance"))
  expect_equal(unname(c_6$exceedance["cv"]),   v("six-monthly (1999-2004)", "cv_exceedance"))
  expect_equal(c_m$best$family, "pure delay")
  expect_equal(c_6$best$family, "pure delay")
})

test_that("convolve_lag follows the census convention", {
  X <- c(5, 9, 2, 7, 4, 8)
  # pure delay at bin k reads X[t - lag0 - k]; pre-record lags backfilled with X[1]
  expect_equal(convolve_lag(X, c(0, 1, 0), t_R = 3:6, lag0 = 1), X[c(1, 2, 3, 4)])
  expect_equal(convolve_lag(X, c(0, 1, 0), t_R = 3:6, lag0 = 0), X[c(2, 3, 4, 5)])
  expect_equal(convolve_lag(X, c(0, 0, 1), t_R = 3:6, lag0 = 1, missing = "drop"),
               c(NA, 5, 9, 2))
  expect_equal(convolve_lag(X, c(0, 0, 1), t_R = 2:6, lag0 = 1), X[c(1, 1, 1, 2, 3)])
  # a flat kernel on a constant record is constant, whatever the profile
  expect_equal(convolve_lag(rep(3, 10), c(.2, .5, .3), t_R = 4:10), rep(3, 7))
  # the searched maximum over kernels is attained by a pure delay
  set.seed(2)
  Y <- rpois(30, 15); k <- lag_kernels(4, n_dirichlet = 300)
  ce <- lag_ceiling(Y, Y[6:30], t_R = 6:30, K = 4, kernels = k)
  expect_equal(ce$best$family, "pure delay")
})

test_that("gini behaves at the extremes", {
  expect_equal(rain_gini(rep(3, 8)), 0)
  expect_equal(rain_gini(c(0, 0, 0, 20)), 0.75)
  expect_true(is.na(rain_gini(c(0, 0))))
})
