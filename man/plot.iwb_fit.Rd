\name{plot.iwb_fit}
\alias{plot.iwb_fit}
\title{Plot survival curves or borrowing weights}
\description{Plot fitted curves with an optional pointwise IWB interval, or plot unscaled and effective borrowing weights. Legends appear below the plotting region; no background grid is drawn.}
\usage{
\method{plot}{iwb_fit}(x, type = c("survival", "weight"),
  components = c("iwb", "obs", "ext", "beta"), interval = TRUE,
  xlab = "Time", ylab = NULL, ...)
}
\arguments{
\item{x}{An \code{iwb_fit} object.}
\item{type}{Survival or borrowing-weight plot.}
\item{components}{One or more survival components to display.}
\item{interval}{Display the pointwise IWB interval when IWB is selected.}
\item{xlab,ylab}{Axis labels. A default y label is chosen when \code{NULL}.}
\item{...}{Additional arguments to the initial \code{plot} call, excluding arguments already set by this method (including \code{xlim}, \code{ylim} and \code{type}).}
}
\value{The input object, invisibly.}
\details{Use a standard R vector graphics device such as \code{pdf} or \code{svg} to export publication-quality plots. A device of at least 7 by 5 inches is recommended for the bottom legend.}
