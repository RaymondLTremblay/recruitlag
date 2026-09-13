# Reproduces output/phase9a/ceiling_calibration.csv within simulation error.
ref <- read.csv(test_path("reference", "ceiling_calibration.csv"))

test_that("the pooled calibration reproduces the paper's phi and p-values", {
  skip_on_cran()
  ce <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4,
                    kernels = lag_kernels(4, bin = 6, unit = "mo"))
  set.seed(20260903)
  cal <- ceiling_calibration(ce, nsim = 5000)
  r6 <- ref[ref$dataset == "six-monthly (1999-2004)", ]
  expect_equal(cal$phi[1], r6$phi[1], tolerance = 1e-8)
  expect_lt(abs(cal$p_gini[cal$null == "N0 flat mean"] - r6$p_gini[r6$null == "N0 flat mean"]), 0.03)
  expect_lt(abs(cal$p_gini[cal$null == "N1 best-kernel mean"] - r6$p_gini[r6$null == "N1 best-kernel mean"]), 0.03)
  # the pooled comparison has no power: neither null is refuted
  expect_gt(min(cal$p_gini), 0.1)
})

test_that("ceiling_calibration() and host_lag_test() refuse recruit series that are not counts", {
  k <- lag_kernels(2, n_dirichlet = 20, bin = 1, unit = "yr")
  ce <- lag_ceiling(c(10, 12, 9, 11, 10, 13, 9, 10), R = c(1.5, 30.2, 2.1, 0.4, 8.8, 1.2), t_R = 3:8,
                    K = 2, lag0 = 0, kernels = k)
  expect_equal(unname(ce$exceedance["gini"]) > 1, TRUE)          # the ceiling itself is scale-free
  expect_error(ceiling_calibration(ce, nsim = 10), "non-whole values")
  d <- data.frame(unit = rep(c("a", "b"), each = 6), period = rep(1:6, 2), reproduction = 10,
                  recruits = c(1.5, 2, 0, 3, 1, 0, 2, 2, 1.2, 0, 4, 1))
  expect_error(host_lag_test(d, K = 1, lag0 = 0, kernels = lag_kernels(1, bin = 1, unit = "yr"), nsim = 10),
               "non-whole values")
})
