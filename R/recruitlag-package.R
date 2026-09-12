#' @keywords internal
#'
#' @section Reading the numbers:
#'
#' Every quantity the package reports is listed here with its range and the
#' direction that means "more episodic". Each function's help page carries
#' the full account under "Scale and how to read it".
#'
#' \tabular{lll}{
#'   **Quantity (function)** \tab **Range** \tab **Reading** \cr
#'   `gini` ([gini()]) \tab 0 to \eqn{(n-1)/n} \tab 0 is even; the maximum is reached when one period holds the whole total. The maximum depends on the number of periods, so compare only series of equal length. \cr
#'   `cv` ([cv()]) \tab 0 to \eqn{\sqrt{n}}{sqrt(n)} \tab As for `gini`, computed with the `n - 1` divisor of [stats::sd()]. \cr
#'   `zeros` ([concentration()]) \tab 0 to `n` \tab Periods with no count at all. \cr
#'   `max_share` ([concentration()]) \tab \eqn{1/n} to 1 \tab Share of the total falling in the largest single period. \cr
#'   weights `w` ([lag_kernels()]) \tab sum to 1 \tab How memory is spread over the horizon, never how much recruitment there is. \cr
#'   `alpha` ([lag_kernels()]) \tab above 0 \tab Below 1 gives sparse kernels close to pure delays; 1 is uniform on the simplex; above 1 gives flat kernels. \cr
#'   `rho` ([lag_kernels()]) \tab 0 to 1 \tab Geometric decay. Near 0 puts all weight on the first bin, near 1 approaches equal weights. \cr
#'   `theorem`, `ceiling` ([lag_ceiling()]) \tab as `gini` and `cv` \tab Reproduction's own concentration, which no delay can exceed. \cr
#'   `exceedance` ([lag_ceiling()]) \tab above 0, centred on 1 \tab 1 is exactly at the bound; above 1 is more concentrated than any delay predicts. A description of the gap, not a test. \cr
#'   `vmr` ([ceiling_calibration()]) \tab above 1 \tab 1 is Poisson and larger is more clumped. Below 1 the negative binomial does not exist. \cr
#'   `phi` ([phi_moment()]) \tab above 0 \tab Variance is \eqn{\mu + \mu^2/\phi}{mu + mu^2/phi}. Small is strongly clumped and very large is Poisson. Density dependent, so not comparable across means. \cr
#'   `p_gini`, `p_cv`, `p_zeros` ([ceiling_calibration()]) \tab 0 to 1 \tab Fraction of simulated series at least as extreme as the observed one. Near 0.5 is ordinary. \cr
#'   `p_silent`, `p_gini`, `p_cv`, `p_max`, `verdict` ([host_lag_test()]) \tab 1/nsim to 1 \tab As above, with the verdict the maximum over profiles, that is, the best case the lag hypothesis can make for itself. Nothing below 1/nsim can be resolved. \cr
#' }
#'
#' **The package applies no thresholds.** There is no Gini above which a
#' record is declared episodic, no `phi` above which counts are declared
#' clumped, and no probability below which a delay is declared refuted. The
#' functions return numbers. A bright-line cutoff on the probabilities is
#' what the ASA statement on p-values advises against (Wasserstein and Lazar
#' 2016), and for the concentration indices no cutoff is taken from the
#' literature, because the level of an index carries no verdict on its own:
#' the inference lives in the comparison between series of equal length.
#'
#' @references
#' Wasserstein, R. L. and Lazar, N. A. (2016) The ASA statement on p-values:
#' context, process, and purpose. *The American Statistician* 70: 129-133.
#' \doi{10.1080/00031305.2016.1154108}
#'
#' @importFrom rlang .data
"_PACKAGE"
