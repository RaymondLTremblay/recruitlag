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
