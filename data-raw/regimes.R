# Builds data/regimes.rda: six invented census records that span the cases the
# test is meant to separate, as ONE long tibble (`regimes`, one row per record,
# unit and period) plus the generating rules (`regimes_truth`).
# Three reproductive regimes (seasonal, intermediate, aseasonal) crossed with
# two recruitment mechanisms:
#   "delayed": recruits are a true delayed consequence of each unit's own
#              reproduction (a fixed lag profile, negative-binomial counts);
#   "gated":   recruits arrive only when an external gate is open, at a rate
#              proportional to each unit's average reproduction, so they are
#              episodic and do not track the reproductive record.
# Every record has H units censused at T periods, reproduction at all T
# periods and recruits scored at periods K + 2 to T (lag0 = 1); recruits are
# NA at unscored periods. Nothing here is data.
# Run from the package root: source("data-raw/regimes.R")
library(tibble)

# The seed was chosen so that the true-delay records are typical draws (their
# host-level P lies in the middle of its null distribution); the invented-data
# vignette shows that distribution over replicate draws.
set.seed(2)
H <- 20L; T <- 24L; K <- 4L; lag0 <- 1L
t_R <- seq.int(K + 2L, T)
w_true <- c(0.10, 0.25, 0.35, 0.20, 0.10)        # the true delay, 1 to 5 periods back
gate <- rep(0.03, T); gate[c(11, 12, 20)] <- 1    # open at three of 24 periods
base <- round(exp(rnorm(H, log(12), 0.6)))         # units differ in size

season <- function(a) {                            # a = 0 aseasonal, 1 fully seasonal
  s <- exp(-((seq_len(T) %% 12 - 4)^2) / 3)        # one peak per 12 periods
  (1 - a) + a * s / mean(s)
}
make <- function(record, a, mechanism, label) {
  s <- season(a)
  X <- t(sapply(base, function(b) rpois(T, b * s)))           # reproduction, unit x period
  mu <- switch(mechanism,
    delayed = t(sapply(seq_len(H), function(h)
      0.08 * recruitlag::convolve_lag(X[h, ], w_true, t_R, lag0))),
    gated   = t(sapply(seq_len(H), function(h) 0.35 * mean(X[h, ]) * gate[t_R])))
  R <- matrix(rnbinom(H * length(t_R), size = 2, mu = mu), H, length(t_R))
  Rfull <- matrix(NA_real_, H, T); Rfull[, t_R] <- R
  tibble(record = record, label = label, seasonality = a, mechanism = mechanism,
         unit = rep(sprintf("u%02d", seq_len(H)), each = T),
         period = rep(seq_len(T), times = H),
         reproduction = as.vector(t(X)), recruits = as.vector(t(Rfull)))
}
regimes <- rbind(
  make("seasonal_delayed",     1.0, "delayed", "seasonal reproduction, true delay"),
  make("seasonal_gated",       1.0, "gated",   "seasonal reproduction, gated recruitment"),
  make("intermediate_delayed", 0.5, "delayed", "intermediate seasonality, true delay"),
  make("intermediate_gated",   0.5, "gated",   "intermediate seasonality, gated recruitment"),
  make("aseasonal_delayed",    0.0, "delayed", "aseasonal reproduction, true delay"),
  make("aseasonal_gated",      0.0, "gated",   "aseasonal reproduction, gated recruitment"))
regimes_truth <- tibble(
  mechanism = c(rep("delayed", K + 1L), rep("gated", T)),
  index = c(seq_len(K + 1L), seq_len(T)),
  what = c(rep("lag bin (periods before the recruit census)", K + 1L), rep("period", T)),
  value = c(w_true, gate))
save(regimes, regimes_truth, file = "data/regimes.rda", compress = "xz", version = 2)
