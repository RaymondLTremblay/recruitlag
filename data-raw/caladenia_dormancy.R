# caladenia_dormancy: Table 2 of Tremblay, Perez, Larcombe, Brown, Quarmby,
# Bickerton, French and Bould (2009) Dormancy in Caladenia: a Bayesian approach
# to evaluating latency. Australian Journal of Botany 57: 340-350.
# Posterior means and 95% credible intervals; n and years from Table 1.
caladenia_dormancy <- data.frame(
  species = c("C. amoena", "C. argocalla", "C. clavigera", "C. elegans", "C. graniticola",
              "C. macroclavia", "C. oenochila", "C. rosella", "C. valida"),
  p_resight    = c(0.946, 0.670, 0.895, 0.830, 0.793, 0.630, 0.843, 0.880, 0.732),
  p_resight_lo = c(0.908, 0.630, 0.786, 0.686, 0.641, 0.554, 0.750, 0.743, 0.686),
  p_resight_hi = c(0.976, 0.709, 0.970, 0.935, 0.912, 0.703, 0.918, 0.974, 0.775),
  p_survival    = c(0.901, 0.875, 0.916, 0.822, 0.892, 0.932, 0.881, 0.891, 0.823),
  p_survival_lo = c(0.864, 0.849, 0.840, 0.715, 0.804, 0.884, 0.816, 0.797, 0.790),
  p_survival_hi = c(0.936, 0.900, 0.977, 0.902, 0.964, 0.976, 0.934, 0.967, 0.853),
  p_dormancy    = c(0.049, 0.289, 0.097, 0.141, 0.185, 0.346, 0.138, 0.108, 0.221),
  p_dormancy_lo = c(0.022, 0.250, 0.028, 0.051, 0.076, 0.269, 0.071, 0.023, 0.183),
  p_dormancy_hi = c(0.084, 0.329, 0.199, 0.271, 0.331, 0.427, 0.224, 0.236, 0.262),
  n_plants = c(80L, 429L, 6L, 22L, 18L, 98L, 22L, 17L, 188L),
  years    = c(12L, 5L, 11L, 5L, 4L, 7L, 11L, 5L, 8L),
  stringsAsFactors = FALSE)
stopifnot(all(abs(caladenia_dormancy$p_dormancy -
                  caladenia_dormancy$p_survival * (1 - caladenia_dormancy$p_resight)) < 0.006))
save(caladenia_dormancy, file = "data/caladenia_dormancy.rda", compress = "bzip2")
