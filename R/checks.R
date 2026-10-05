# Every message the user can hit is written here, in one place, to the same
# rule: say what is wrong, say what was found, and say what to do about it.
# All of them use call. = FALSE, so the reader gets a sentence rather than
# "Error in series_from(X, period, reproduction, recruits) :".

# `class` adds a condition class, so that a caller can recognise one error by
# inherits() instead of by matching its wording.
rl_abort <- function(..., class = NULL) {
  stop(structure(class = c(class, "error", "condition"),
                 list(message = paste0(...), call = NULL)))
}
rl_warn  <- function(...) warning(paste0(...), call. = FALSE)

# A compact list for a message: "3, 7, 11 and 4 more".
rl_list <- function(x, n = 6) {
  x <- as.character(x)
  if (length(x) <= n) return(paste(x, collapse = ", "))
  paste0(paste(x[seq_len(n)], collapse = ", "), " and ", length(x) - n, " more")
}

# "Did you mean ...?" from the column names actually present.
rl_near <- function(name, candidates) {
  if (!length(candidates)) return("")
  d <- utils::adist(tolower(name), tolower(candidates))[1, ]
  best <- candidates[d <= max(2, ceiling(nchar(name) / 3))]
  if (!length(best)) return("")
  paste0(" Did you mean ", paste(sprintf('"%s"', utils::head(best, 3)), collapse = " or "), "?")
}

# roles is a named character vector: argument name -> column name given.
rl_check_columns <- function(data, roles) {
  if (!is.data.frame(data))
    rl_abort("The data must be a data frame with one row per period, ",
             "or per unit and period. It is ", class(data)[1], ".")
  if (!nrow(data)) rl_abort("The data frame has no rows.")
  have <- names(data)
  miss <- roles[!roles %in% have]
  if (length(miss)) {
    one <- function(arg, col) paste0(
      '  ', arg, ' = "', col, '"  was not found.', rl_near(col, have))
    rl_abort("These column names are not in the data:\n",
             paste(mapply(one, names(miss), miss), collapse = "\n"),
             "\nThe data frame has: ", rl_list(have, 12), "\n",
             "Name your own columns with ",
             paste(sprintf('%s = "..."', names(roles)), collapse = ", "), ".")
  }
  invisible(TRUE)
}

# The checks that apply to any long data frame the package accepts.
# Called once at each entry point. `what` names the calling function so the
# message says where the problem was noticed.
rl_check_long <- function(data, period, reproduction, recruits, unit = NULL, what = "") {
  roles <- c(period = period, reproduction = reproduction, recruits = recruits)
  if (!is.null(unit)) roles <- c(unit = unit, roles)
  rl_check_columns(data, roles)

  p <- data[[period]]
  if (!is.numeric(p))
    rl_abort('The period column "', period, '" must be numeric: the census index, ',
             "an integer counting from 1. It is ", class(p)[1], ". ",
             "If it holds calendar years, convert with period = year - min(year) + 1.")
  if (anyNA(p))
    rl_abort('The period column "', period, '" has ', sum(is.na(p)), " missing values. ",
             "Every row needs a census index; a census with no data is a row of blanks, ",
             "not a row with a blank period.")
  if (any(p != round(p)))
    rl_abort('The period column "', period, '" has non-whole values (',
             rl_list(unique(p[p != round(p)])), "). Periods are census indices, ",
             "not dates or decimal years.")
  if (any(p < 1))
    rl_abort('The period column "', period, '" has values below 1 (',
             rl_list(sort(unique(p[p < 1]))), "). Periods count from 1.")
  if (min(p) != 1)
    rl_abort('The period column "', period, '" starts at ', min(p), ", not at 1, so the ",
             "record would be read as ", max(p), " periods of which the first ", min(p) - 1,
             " are empty. Periods are census indices counting from 1. ",
             if (min(p) > 1800 && min(p) < 2200)
               paste0("These look like calendar years: convert with period = ", period,
                      " - ", min(p), " + 1.")
             else "Re-index with period = period - min(period) + 1.")

  if (!is.null(unit)) {
    u <- data[[unit]]
    if (anyNA(u) || any(trimws(as.character(u)) == ""))
      rl_abort('The unit column "', unit, '" has ',
               sum(is.na(u) | trimws(as.character(u)) == ""), " blank values. ",
               "Every row must say which unit it belongs to. If your file names the ",
               "unit only on its first row, carry the name down before calling this.")
    key <- paste(as.character(u), p, sep = "\r")
    if (anyDuplicated(key)) {
      dup <- unique(key[duplicated(key)])
      rl_abort(length(dup), " unit and period combination(s) appear more than once: ",
               rl_list(gsub("\r", " at period ", dup)), ". ",
               "Each unit is expected once per period. Collapse the duplicates first, ",
               "and decide deliberately how to combine them.")
    }
  }

  # Interior gaps. A period that is absent is read as a period with no
  # reproduction, which is a claim about the plants rather than about the
  # fieldwork, so it is worth saying out loud.
  gaps <- setdiff(seq(min(p), max(p)), unique(p))
  if (length(gaps))
    rl_warn("Period(s) ", rl_list(gaps), " are absent from the data but lie inside ",
            "the record (", min(p), " to ", max(p), "). They will be read as periods ",
            "with no reproduction and no recruits, which is not the same as periods ",
            "that were never censused. Add them as rows of blanks if they were not censused.")

  x <- data[[reproduction]]
  if (!is.numeric(x))
    rl_abort('The reproduction column "', reproduction, '" must be numeric. It is ',
             class(x)[1], ".")
  if (any(x < 0, na.rm = TRUE))
    rl_abort('The reproduction column "', reproduction, '" has negative values (',
             rl_list(sort(unique(x[!is.na(x) & x < 0]))), "). ",
             "Reproduction is a count or a rate and cannot be negative.")
  if (anyNA(x))
    rl_warn("The reproduction column \"", reproduction, "\" has ", sum(is.na(x)),
            " missing values, and they are read as zero. A year that was not measured ",
            "and a year with no reproduction are different claims, and the second one ",
            "lowers the ceiling. If those rows are unmeasured rather than empty, ",
            "drop them before calling ", if (nzchar(what)) what else "this function", ".")

  r <- data[[recruits]]
  if (!is.numeric(r))
    rl_abort('The recruits column "', recruits, '" must be numeric, with NA where ',
             "recruits were not scored. It is ", class(r)[1], ".")
  if (any(r < 0, na.rm = TRUE))
    rl_abort('The recruits column "', recruits, '" has negative values (',
             rl_list(sort(unique(r[!is.na(r) & r < 0]))), ").")
  if (all(is.na(r)))
    rl_abort('No period has a recruit count: "', recruits, '" is NA in all ',
             length(r), " rows. There is nothing to compare against the ceiling.")
  invisible(TRUE)
}

# Kernels, profiles and objects: the arguments users get wrong second most often.
rl_check_kernels <- function(kernels, K, arg = "kernels") {
  if (!inherits(kernels, "lag_kernels"))
    rl_abort("`", arg, "` must be a lag_kernels object, as built by lag_kernels(). It is ",
             class(kernels)[1], ". Build one with ", arg, " = lag_kernels(", K, ").")
  if (attr(kernels, "K") != K)
    rl_abort("The kernels were built for K = ", attr(kernels, "K"), " but K = ", K,
             " was passed, so the weights and the horizon disagree. ",
             "Build the kernels with the same K: ", arg, " = lag_kernels(", K, ").")
  invisible(TRUE)
}

rl_check_class <- function(x, cls, arg, made_by) {
  if (!inherits(x, cls))
    rl_abort("`", arg, "` must be ", if (grepl("^[aeiou]", cls)) "an " else "a ", cls,
             " object, as returned by ", made_by, ". It is ", class(x)[1], ".")
  invisible(TRUE)
}

# Lag weights: one per bin, non-negative, present, not all zero. Used by every
# function that applies a profile, so that they all refuse the same inputs.
rl_check_weights <- function(w, arg = "w") {
  if (!is.numeric(w) || !length(w))
    rl_abort("`", arg, "` must be a numeric vector of lag weights, one per bin. It is ",
             class(w)[1], " of length ", length(w), ".")
  if (anyNA(w) || any(w < 0))
    rl_abort("The lag weights `", arg, "` must all be non-negative and present. ",
             "Given: ", rl_list(format(w)), ".")
  if (sum(w) == 0)
    rl_abort("The lag weights `", arg, "` sum to zero, so there is no profile to apply. ",
             "At least one bin must carry weight.")
  invisible(TRUE)
}

# A whole-number argument (a period offset, a horizon, a number of replicates).
rl_check_whole <- function(x, arg, min = 0) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x != round(x) || x < min)
    rl_abort("`", arg, "` must be a single whole number of ", min, " or more. It was given as ",
             paste(format(x), collapse = ", "), ".")
  invisible(TRUE)
}

# The negative binomial is a distribution for counts: whole, non-negative. Used by every
# function that simulates or models counts; `why` names the caller's reason.
rl_check_counts <- function(M, column, role,
                            why = paste0("lag_ceiling_stan() fits a negative-binomial model, which is a ",
                                         "distribution for counts, so the recruits must be whole numbers.")) {
  v <- M[!is.na(M)]
  bad <- v[v != round(v)]
  what <- if (is.null(column)) paste0("The ", role, " series") else paste0('The ', role, ' column "', column, '"')
  if (length(bad))
    rl_abort(what, " has non-whole values (",
             rl_list(format(utils::head(sort(unique(bad)), 4))), "). ", why,
             " If this column is a rate, a biomass or a model estimate, pass the count it was ",
             "made from instead; lag_ceiling() and the two bootstrap functions accept any ",
             "non-negative series, because the concentration indices are scale-free.")
  invisible(TRUE)
}

rl_check_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1 || is.na(level) || level <= 0 || level >= 1)
    rl_abort("level must be a single number strictly between 0 and 1, the coverage of the ",
             "interval (for example 0.9). It was given as ", paste(format(level), collapse = ", "), ".")
  invisible(TRUE)
}

# cmdstanr and CmdStan, or a message saying how to get them.
rl_need_cmdstan <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE))
    rl_abort("lag_ceiling_stan() needs the cmdstanr package, which is not on CRAN. Install it with\n",
             '  install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))\n',
             "and then CmdStan itself with cmdstanr::install_cmdstan(). The two bootstrap ",
             "functions, lag_ceiling_boot() and lag_ceiling_bayesboot(), need neither.")
  v <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)
  if (is.null(v))
    rl_abort("cmdstanr is installed but CmdStan itself was not found. Install it with ",
             "cmdstanr::install_cmdstan(), or point cmdstanr at an existing installation with ",
             "cmdstanr::set_cmdstan_path().")
  invisible(TRUE)
}
