# RELEASE_PLAN.md: recruitlag to Zenodo and CRAN

Handoff written 2026-10-05, when `recruitlag` moved out of the Claude project
"Lepanthes_eltoroensis_seed_recruitment" into a project of its own. A new
Claude project cannot see the old project's memory, so everything a session
needs to continue the release is here or in `CLAUDE.md`.

## Where things stand (checked 2026-10-05)

* `master` is at `295a9df` and equals `origin/master` on GitHub
  (https://github.com/RaymondLTremblay/recruitlag, public since 2026-09-21).
* R CMD check on 2026-09-21: 0 errors, 0 warnings, 1 note (cmdstanr suggested
  but not available in the checking environment). 132 tests pass.
* **Tags are inconsistent.** Locally there are two tags on `295a9df`: `v0.1.0`
  (annotated, "the version cited in the Oikos paper") and `v.0.1.0` (stray,
  with a dot). **GitHub has only `v.0.1.0`.** The annotated `v0.1.0` was never
  pushed. Decide before anything else: delete `v.0.1.0` locally and on
  GitHub, and push `v0.1.0`. Zenodo takes its version label from the GitHub
  release, so the stray name would end up in the DOI record.
* No GitHub release and no Zenodo DOI yet.
* The Oikos paper (manuscript.qmd in the Lepanthes repository) cites
  version 0.1.0 through the `@recruitlag` bib entry, whose DOI is pending
  this release. `09e_ceiling_intervals.qmd` loads the package. The paper is
  NOT submitted yet: the package must be public with a DOI first.

## Order of work

### A. Zenodo (needed first: the paper waits on this DOI)

1. Fix the tags (above).
2. Log in to Zenodo with GitHub, go to the GitHub settings page in Zenodo and
   switch the `recruitlag` repository ON. This must happen BEFORE the release
   is made; Zenodo archives only releases made after the switch.
3. Optional but recommended: add a `.zenodo.json` (or `CITATION.cff`) so the
   Zenodo record carries the ORCID 0000-0002-8588-4372, the licence and
   keywords instead of whatever Zenodo guesses. Add `^\.zenodo\.json$` or
   `^CITATION\.cff$` to `.Rbuildignore`.
4. On GitHub, create a release from tag `v0.1.0`. Zenodo then mints two DOIs:
   one for 0.1.0 and a concept DOI for all versions. The paper cites the
   0.1.0 DOI; README and `inst/CITATION` can carry the concept DOI.
5. Put the DOI into: `inst/CITATION`, README (badge), and the `@recruitlag`
   entry of `references.bib` in the Lepanthes repository. Then the Lepanthes
   Zenodo data deposit can be finished (its own checklist:
   `zenodo_upload_checklist.md` in that repository).

### B. CRAN (after Zenodo; may well become 0.1.1)

CRAN reviewers often ask for changes, and any change means a new version.
That is fine: the paper cites 0.1.0 on Zenodo, and CRAN can receive 0.1.1.

Before submitting, items to check (none verified yet):

* Title case: CRAN's checker wants "be" lowercase in the Title field
  (seen in an earlier check). Title must also not start with the package name.
* `inst/CITATION` and `NEWS.md` say the paper is "submitted"; update to match
  the real state at the time of submission.
* `cmdstanr` is in Suggests with `Additional_repositories:
  https://stan-dev.r-universe.dev`. Every use must stay conditional (tests
  skip, `lag_ceiling_stan()` errors cleanly, the vignette chunk is
  `eval = FALSE`). Confirm on a machine WITHOUT cmdstanr.
* `URL` and `BugReports` in DESCRIPTION use lowercase
  `raymondltremblay`; the remote is `RaymondLTremblay`. GitHub is case
  insensitive, but use one spelling everywhere.
* `README.md` is excluded by `.Rbuildignore`; that is allowed.
* Run checks off the Mac too: `devtools::check_win_devel()`,
  `rhub::rhub_check()`, and `R CMD check --as-cran` on the built tarball.
* Write `cran-comments.md` (test environments, the check results, "This is a
  new submission"); add `^cran-comments\.md$` and `^CRAN-SUBMISSION$` to
  `.Rbuildignore`.
* Submit with `devtools::submit_cran()`, confirm the email CRAN sends.

## Working environment notes (from the old project's memory)

* R is not preinstalled in the Cowork cloud container but can be installed:
  `apt-get update` first (it warns about a blocked docker repo but works),
  then `apt-get install -y --fix-missing r-base-core r-base-dev
  r-cran-ggplot2 r-cran-tibble r-cran-rlang r-cran-testthat r-cran-knitr
  r-cran-rmarkdown r-cran-dplyr r-cran-tidyr r-cran-roxygen2 r-cran-devtools
  r-cran-posterior`. Stage the package files (skip `man/` and `data-raw/`;
  roxygen regenerates `man/`), copy out of the read-only uploads folder,
  `chmod -R u+w`, then roxygenise, build, and
  `_R_CHECK_FORCE_SUGGESTS_=false _R_CHECK_CRAN_INCOMING_=false R CMD check --no-manual`.
* The container's proxy blocks CRAN, r-universe and GitHub, so cmdstanr,
  CmdStan, win-builder and rhub are all run from the Mac.
* git on the mounted folder from the Cowork device shell needs delete
  permission for the folder (otherwise `git add` fails on
  `.git/objects/tmp_obj_*`), and a stale `.git/index.lock` must then be
  removed. Commit with `-c user.name=... -c user.email=...`. RLT pushes with
  GitHub Desktop.
* The repository lives in Dropbox. Watch for `*conflicted copy*` files
  (already in `.gitignore`).
