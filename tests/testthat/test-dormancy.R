# recruit_triage(): a first-sighting series against the record's own dormancy.
mk <- function() {
  M <- matrix(NA_character_, 6, 6, dimnames = list(paste0("p", 1:6), 2001:2006))
  M[1, ] <- c("F", "F", NA, "F", "V", "V")     # one unseen run of 1 between sightings
  M[2, ] <- c("V", NA, NA, NA, "F", NA)        # run of 3, then a trailing gap (not counted)
  M[3, ] <- c(NA, "F", "F", "F", "F", "F")     # first seen flowering at period 2
  M[4, ] <- c(NA, NA, "V", "V", "F", "F")      # first seen vegetative at period 3
  M[5, ] <- c(NA, NA, NA, NA, NA, "F")         # first seen flowering at period 6
  M[6, ] <- c("F", "F", "F", "F", "F", "F")
  M
}

test_that("first sightings, states and dormancy runs are counted as documented", {
  tr <- recruit_triage(mk())
  expect_s3_class(tr, "recruit_triage")
  expect_equal(tr$n_plants, 6)
  expect_equal(tr$first$first_seen, c(3, 1, 1, 0, 0, 1))
  expect_equal(tr$first$first_seen_flowering, c(2, 1, 0, 0, 0, 1))
  expect_equal(unname(tr$first_state), c(2, 1))
  expect_equal(as.integer(tr$dormancy$runs), c(1, 0, 1, 0, 0, 0))   # runs of 1 and 3, trailing gap ignored
  expect_equal(tr$dormancy$longest_run, 3)
  expect_equal(tr$dormancy$share_plants_ever_dormant, 2 / 6)
  # detection: seen / bracketed periods = (5 + 2 + 5 + 4 + 1 + 6) / (6 + 5 + 5 + 4 + 1 + 6)
  expect_equal(tr$detection, 23 / 27)
  expect_null(tr$first$expected_returns)
  expect_output(print(tr), "first sightings after period 1: 3")
})

test_that("expected returns follow N1 (1 - r) s q^(k - 2), and the species shortcut works", {
  tr <- recruit_triage(mk(), s = 0.8, r = 0.5)
  N1 <- 3; q <- 0.8 * 0.5
  expect_equal(tr$first$expected_returns[2:6], N1 * 0.5 * 0.8 * q^(0:4))
  expect_equal(tr$first$excess[2], 1 - N1 * 0.5 * 0.8)
  expect_equal(tr$dormancy$improbable_from, ceiling(log(0.01) / log(q)))
  tv <- recruit_triage(mk(), species = "C. valida")
  cd <- caladenia_dormancy
  expect_equal(tv$s, cd$p_survival[cd$species == "C. valida"])
  expect_equal(tv$r, cd$p_resight[cd$species == "C. valida"])
  expect_error(recruit_triage(mk(), species = "C. nope"), "not in caladenia_dormancy")
  expect_error(recruit_triage(mk(), s = 1.2, r = 0.5), "probabilities")
  # r from the record when only s is given
  ts <- recruit_triage(mk(), s = 0.9)
  expect_equal(ts$r, ts$detection)
  expect_s3_class(rain_sightings(tv), "ggplot")
})

test_that("the long-data interface agrees with the matrix, and the codes are checked", {
  M <- mk()
  long <- data.frame(plant = rep(rownames(M), each = ncol(M)), year = rep(2001:2006, times = nrow(M)),
                     st = as.vector(t(M)), stringsAsFactors = FALSE)
  a <- recruit_triage(M); b <- recruit_triage(long, plant = "plant", period = "year", state = "st")
  expect_equal(a$first, b$first); expect_equal(a$detection, b$detection)
  expect_error(recruit_triage(M, seen = c("x", "y")), "none of the states matches")
  expect_error(recruit_triage(M[, 1, drop = FALSE]), "at least two")
  expect_error(recruit_triage(long[0, ], plant = "plant", period = "year", state = "st"), "no rows")
})

test_that("caladenia_dormancy is internally consistent", {
  cd <- caladenia_dormancy
  expect_equal(nrow(cd), 9)
  expect_true(all(abs(cd$p_dormancy - cd$p_survival * (1 - cd$p_resight)) < 0.006))
  expect_true(all(cd$p_resight_lo <= cd$p_resight & cd$p_resight <= cd$p_resight_hi))
})
