## Test environments

* local: macOS 26 (Apple silicon), R 4.6.1, CmdStan 2.39.0 (the Stan route exercised);
  R CMD check --as-cran, 2026-10-05
* Ubuntu 24.04, R 4.3.3, cmdstanr not installed (R CMD check --as-cran, 2026-10-05)
* CRAN incoming pretest (2026-10-05): Windows R-devel (r90605 ucrt), 1 NOTE;
  Debian R-devel, 1 NOTE. Both notes are the incoming-feasibility note below.

## R CMD check results

0 errors | 0 warnings | 2 notes

* checking CRAN incoming feasibility ... NOTE
  New submission
  Suggests or Enhances not in mainstream repositories: cmdstanr
  Availability using Additional_repositories specification: cmdstanr yes

  This is a new submission. cmdstanr is not on CRAN. It is declared in Suggests with
  Additional_repositories: https://stan-dev.r-universe.dev, and every use is
  conditional: lag_ceiling_stan() tests requireNamespace() and gives an
  installation message, the test skips, and the vignette chunk that calls it
  has eval = FALSE. The package installs, checks and runs its vignettes
  without cmdstanr; the two bootstrap routes to the same quantity need
  nothing beyond Imports.

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.

  A property of the local macOS check machine (an old HTML Tidy), not of the
  package. The Ubuntu check did not produce it.

The incoming-feasibility check flags as possibly misspelled: Lepanthes and
eltoroensis (the species), epiphytic (an ordinary word) and reproduction's.

## Downstream dependencies

This is a new submission; there are no reverse dependencies.

## Submission record

* 2026-10-05: 0.1.0 submitted (devtools::submit_cran, confirmed). Auto-processed
  the same day, pretests clean (1 NOTE each on Windows and Debian), pending
  manual inspection by the CRAN team (up to 10 working days).
