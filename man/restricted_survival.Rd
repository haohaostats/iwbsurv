\name{restricted_survival}
\alias{restricted_survival}
\title{Integrate a fitted survival curve over a time window}
\description{Calculate the restricted survival contribution by trapezoidal integration of the interpolated mean survival curve.}
\usage{
restricted_survival(object, from = object$t_obs, to = object$t_max,
  component = c("iwb", "obs", "ext", "beta"), n = 401L)
}
\arguments{
\item{object}{An \code{iwb_fit} object.}
\item{from,to}{Window endpoints, with \eqn{0\leq from<to\leq t_{\max}}.}
\item{component}{Survival component to integrate.}
\item{n}{Number of equally spaced integration points, at least two.}
}
\value{A scalar estimated restricted survival contribution, in the input time units. This is a point estimate; no interval is calculated by this function.}
