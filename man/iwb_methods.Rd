\name{iwb_methods}
\alias{print.iwb_fit}
\alias{summary.iwb_fit}
\alias{predict.iwb_fit}
\alias{as.data.frame.iwb_fit}
\title{Summarize and extract IWB results}
\description{Summarize terminal survival, interpolate grid summaries or extract all component curves.}
\usage{
\method{print}{iwb_fit}(x, ...)
\method{summary}{iwb_fit}(object, ...)
\method{predict}{iwb_fit}(object, times = object$time,
  component = c("iwb", "obs", "ext", "beta"),
  type = c("survival", "weight"), ...)
\method{as.data.frame}{iwb_fit}(x, row.names = NULL, optional = FALSE, ...)
}
\arguments{
\item{x,object}{An \code{iwb_fit} object.}
\item{times}{Prediction times between zero and the fitted horizon.}
\item{component}{Curve to extract: IWB, RP continuation, attached external tail, or Beta blending.}
\item{type}{Return survival summaries or borrowing weights.}
\item{row.names}{Optional row names for the combined data frame.}
\item{optional}{Standard data-frame method argument; no effect.}
\item{...}{Reserved additional arguments.}
}
\details{Predictions linearly interpolate stored summaries; they do not refit or extend the fitted horizon. For weights, \code{pi} is unscaled and \code{effective} equals \eqn{\kappa^\star\pi}. Component selection has no effect for \code{type = "weight"}.}
\value{\code{print} invisibly returns the input. \code{summary} returns terminal survival means and interval limits for all four components. \code{predict} returns a data frame with \code{time}, \code{mean}, \code{lower}, \code{upper}, or with \code{time}, \code{pi}, \code{effective} for weights. \code{as.data.frame} returns all survival summaries with an additional \code{component} column.}
