# Reproduces output/phase9d of the companion paper (host_level_lag_null_named.csv):
# the deterministic quantities exactly, the Monte Carlo probabilities within
# their simulation error.
ref <- read.csv(test_path("reference", "host_level_lag_null_named.csv"))
ref <- ref[ref$input == "inflorescences" & ref$missing == "backfill", ]
h <- host_matrices(lepanthes_hosts, unit = "host", reproduction = "inflorescences")

test_that("host_matrices rebuilds the paper's matrices", {
  expect_equal(dim(h$R), c(23, 10)); expect_equal(as.integer(colnames(h$R)), 3:12)
  expect_equal(dim(h$X), c(23, 12)); expect_equal(sum(h$R), 149)
  expect_equal(sum(h$X[, 1]), 509)
})

test_that("expected recruits and the moment phi match the paper for named kernels", {
  for (lab in c("pure delay 0-6 mo", "pure delay 24-30 mo", "uniform 0-30 mo", "geometric rho=0.60")) {
    k <- lag_kernels(4, bin = 6, unit = "mo")
    w <- k[[which(vapply(k, function(z) z$label, character(1)) == lab)]]$w
    mu <- expected_recruits(h$R, h$X, w)
    expect_equal(rowSums(mu), rowSums(h$R))        # scaled to each host's total
    expect_equal(phi_moment(h$R, mu), ref$phi_moment[ref$label == lab], tolerance = 1e-10)
    expect_equal(rain_gini(colSums(mu)), ref$gini_expected_totals[ref$label == lab], tolerance = 1e-10)
  }
})

test_that("the flat reference phi matches", {
  mu <- matrix(rowSums(h$R) / 10, 23, 10, dimnames = dimnames(h$R))
  expect_equal(phi_moment(h$R, mu), ref$phi_moment[ref$label == "flat reference"], tolerance = 1e-10)
})

test_that("host-level probabilities agree with the paper within simulation error", {
  skip_on_cran()
  set.seed(20260910)
  ht <- host_lag_test(lepanthes_hosts, unit = "host", reproduction = "inflorescences", K = 4, nsim = 4000,
                      kernels = lag_kernels(4, bin = 6, unit = "mo",
                                            extra = list(`4d posterior-mean profile` = lepanthes_profile$weight)),
                      phi = 0.08903799184256246)
  tab <- ht$table[ht$table$phi_source == "moment", ]
  for (lab in c("flat reference", "pure delay 24-30 mo", "uniform 0-30 mo", "4d posterior-mean profile")) {
    r <- ref[ref$label == lab, ]; t1 <- tab[tab$label == lab, ]
    expect_lt(abs(t1$p_silent - r$p_zero_moment), 0.02)
    expect_lt(abs(t1$p_gini - r$p_gini_moment), 0.02)
    expect_lt(abs(t1$p_cv - r$p_cv_moment), 0.02)
  }
  # the verdict: no profile makes three silent censuses ordinary
  expect_lt(ht$verdict$max_p[ht$verdict$phi_source == "moment" & ht$verdict$statistic == "silent"], 0.01)
  expect_lt(ht$verdict$max_p[ht$verdict$phi_source == "moment" & ht$verdict$statistic == "gini"], 0.03)
})
