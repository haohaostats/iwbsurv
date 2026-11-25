\name{obs_fit_rp}
\alias{obs_fit_rp}
\title{Fit and sample the observation-side survival model}
\description{Fit an intercept-only RP log cumulative-hazard spline after administrative censoring. The returned draws can be reused across anchor scenarios.}
\usage{obs_fit_rp(time, status, t_obs, k = 2, nsim = 2000, seed = 1L)}
\arguments{
\item{time,status}{Follow-up times and binary event indicators; see \code{\link{iwb_fit}}.}
\item{t_obs}{Positive observation cutoff.}
\item{k}{Non-negative integer number of internal spline knots.}
\item{nsim}{Number of coefficient draws, at least 100.}
\item{seed}{Non-negative integer random seed; the caller's random state is restored.}
}
\value{An \code{iwb_observation} list containing \code{fit}, \code{knots}, \code{theta}, \code{theta_hat}, \code{train}, \code{t_obs}, \code{k} and \code{nsim}. Draws are screened on the extrapolation grid by \code{iwb_fit}, not by this function.}
\seealso{\code{\link{iwb_fit}}}
