# Builds the Lepanthes eltoroensis data sets (tibbles) from the archived outputs
# of the companion paper's pipeline (project Lepanthes_eltoroensis_seed_recruitment).
# The CSVs in this folder are copies of:
#   output/phase8a/monthly_series.csv        monthly survey 1994-1996 (system totals)
#   output/flowering_per_census.csv          six-monthly survey 1999-2004 (system totals)
#   output/phase4c/decomposition_infl.csv    recruits per host tree per census
#   output/host_by_census_inputs.csv         inflorescences and adults per host per census
#   output/phase4d/lag_profile.csv           posterior-mean lag profile of the fitted model
# Run from the package root: source("data-raw/lepanthes.R")
library(tibble)

m  <- read.csv("data-raw/monthly_series.csv")
fl <- read.csv("data-raw/flowering_per_census.csv")
d6 <- read.csv("data-raw/decomposition_infl.csv")
hb <- read.csv("data-raw/host_by_census_inputs.csv")
lp <- read.csv("data-raw/lag_profile.csv")

# ---- monthly survey, system totals ------------------------------------------
lepanthes_monthly <- tibble(period = m$time, reproductive_adults = m$repro_ad,
                            recruits = m$recruits)

# ---- six-monthly survey, system totals --------------------------------------
rec6 <- tapply(d6$R, d6$t_plus_1, sum)
lepanthes_census <- tibble(
  period = fl$time, census = fl$census_label, years_since_georges = fl$years_since_georges,
  adults = fl$n_adults_total, juveniles = fl$n_juveniles_total, seedlings = fl$n_seedlings_total,
  inflorescences = fl$total_inflorescences, flowers = fl$total_flowers, fruits = fl$total_fruits,
  recruits = NA_real_)
lepanthes_census$recruits[match(as.integer(names(rec6)), lepanthes_census$period)] <- as.numeric(rec6)

# ---- six-monthly survey, one row per host tree and census -------------------
hosts_id <- sort(unique(d6$pop_num))
hb <- hb[hb$pop_num %in% hosts_id, ]
lepanthes_hosts <- tibble(host = as.character(hb$pop_num), period = hb$time,
                          adults = hb$adults, inflorescences = hb$inflorescences,
                          recruits = NA_real_)
key <- paste(d6$pop_num, d6$t_plus_1)
i <- match(paste(lepanthes_hosts$host, lepanthes_hosts$period), key)
lepanthes_hosts$recruits[!is.na(i)] <- d6$R[i[!is.na(i)]]
lepanthes_hosts <- lepanthes_hosts[order(as.integer(lepanthes_hosts$host), lepanthes_hosts$period), ]

# ---- posterior-mean lag profile of the fitted model -------------------------
lepanthes_profile <- tibble(bin = lp$lag[order(lp$lag)],
                            months = sprintf("%d-%d", 6 * lp$lag[order(lp$lag)], 6 * lp$lag[order(lp$lag)] + 6),
                            weight = lp$mean[order(lp$lag)] / sum(lp$mean))

stopifnot(sum(lepanthes_hosts$recruits, na.rm = TRUE) == 149,
          length(unique(lepanthes_hosts$host)) == 23,
          sum(lepanthes_hosts$inflorescences[lepanthes_hosts$period == 1]) == 509,
          sum(lepanthes_hosts$adults[lepanthes_hosts$period == 1]) == 264,
          sum(lepanthes_monthly$recruits, na.rm = TRUE) == 81,
          all(is.na(lepanthes_hosts$recruits[lepanthes_hosts$period <= 2])))

save(lepanthes_monthly, lepanthes_census, lepanthes_hosts, lepanthes_profile,
     file = "data/lepanthes.rda", compress = "xz", version = 2)
