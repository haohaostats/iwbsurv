\name{iwb_fit}
\alias{iwb_fit}
\title{Fit information-weighted survival extrapolation}
\description{Combine an observation-side Royston--Parmar model with an uncertain, conditionally anchored Gompertz tail.}
\usage{
iwb_fit(time, status, t_obs, t_max, t_a = t_max, Sa_mean, ess = 50,
  step = (t_max - t_obs) / 200, nsim = 2000, k_rp = 2, k_w = 10,
  loggamma_mean = log(0.08), loggamma_sd = 0.5,
  level = 0.95, seed = 1L, obs = NULL, keep_draws = FALSE)
}
\arguments{
\item{time}{Numeric, non-negative follow-up times. Zero times are replaced by \eqn{10^{-6}} during fitting. Missing values are not accepted.}
\item{status}{Event indicator: 1 for an event, 0 for right censoring.}
\item{t_obs}{Positive administrative observation cutoff.}
\item{t_max}{Extrapolation horizon, greater than the cutoff.}
\item{t_a}{Anchor time in \eqn{(t_{\rm obs},t_{\max}]}.}
\item{Sa_mean}{Mean unconditional survival anchor, strictly between zero and one.}
\item{ess}{Positive Beta concentration parameter, also denoted \eqn{\nu}.}
\item{step}{Grid spacing. The cutoff, anchor time and horizon are included explicitly. At least seven points at or beyond the cutoff are required.}
\item{nsim}{Initial coefficient draw count, an integer of at least 100. Ignored when \code{obs} is provided.}
\item{k_rp}{Number of internal RP spline knots. Ignored when \code{obs} is provided.}
\item{k_w}{Monotone weight spline basis dimension, at least five; reduced when required by the grid size.}
\item{loggamma_mean,loggamma_sd}{Mean and standard deviation of the normal distribution of log Gompertz shape. The shape has inverse-time units; choose parameters consistent with \code{time}.}
\item{level}{Pointwise interval probability, strictly between zero and one.}
\item{seed}{Non-negative integer seed no greater than \code{.Machine$integer.max - 7919}. The caller's random-number state is restored.}
\item{obs}{Optional result of \code{\link{obs_fit_rp}} for the same administratively censored data and cutoff.}
\item{keep_draws}{Retain cumulative-hazard matrices for the compatible draws.}
}
\details{
The model is intercept-only. Observations after the cutoff are administratively censored. Coefficient draws use the asymptotic multivariate normal distribution from the RP fit and are screened for positive finite hazards on the positive-time grid.

For each retained observation draw, an anchor is drawn from
\eqn{\mathrm{Beta}(\nu\bar S_a,\nu(1-\bar S_a))} and a positive Gompertz shape from a log-normal distribution. Pairs with \eqn{0<S_a/S_c<1} are retained. Subsequent calculations use the same paired draws. At least 100 compatible pairs are required.

The variance ratio of log cumulative-hazard increments is smoothed by a monotone increasing P-spline on the logit scale. At the cutoff, log-hazard variances supply the limiting ratio. The transformed weight starts at zero. A common factor \eqn{\kappa^\star} constrains the evaluated blended hazards to be non-negative. Both the estimated weight and this factor are fixed across draws for interval calculation.

All four reported components use the retained paired draws. Their cumulative hazards coincide within the observation window. The observation component therefore summarizes this retained ensemble, which can differ from the original unfiltered RP ensemble. Intervals are pointwise empirical quantiles (R quantile type 7), not simultaneous bands.

The Beta(3,3) comparator blends the same components without the common hazard scaling. Hazard validity is evaluated on the computational grid using numerical weight derivatives. The package does not fit covariate effects. Cutoff, anchor and horizon values must be representable to ten decimal places.
}
\value{
An \code{iwb_fit} object. \code{time} contains the grid; \code{iwb}, \code{obs}, \code{ext} and \code{beta} are matrices with columns \code{mean}, \code{lo}, \code{hi}.
\code{weight$pi} is the unscaled weight, with raw and smoothed ratios on \code{weight$t}; \code{kappa * weight$pi} is effective borrowing.
\code{diag} contains the compatible draw count, observation hazard-screen rejection proportion, minimum evaluated IWB and Beta hazards, and terminal unscaled weight.
\code{p_inc} is the incompatible-anchor proportion among observation-valid draws.
Other entries include the training data, fitted RP model, anchor specification, interval level, call and computational settings.
If requested, \code{draws} contains \code{H_obs}, \code{H_ext}, \code{H_iwb} and \code{H_beta}, with rows corresponding to retained pairs.
}
\examples{
d <- rotterdam_example()
fit <- iwb_fit(d$time, d$status,
               t_obs = 5, t_max = 15, Sa_mean = 0.39,
               nsim = 300, step = 0.1, seed = 2026)
predict(fit, times = c(5, 10, 15))
restricted_survival(fit)
}
\seealso{\code{\link{obs_fit_rp}}, \code{\link{predict.iwb_fit}}, \code{\link{plot.iwb_fit}}}
