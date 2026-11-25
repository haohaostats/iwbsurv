\name{rotterdam_example}
\alias{rotterdam_example}
\title{Rotterdam breast cancer example data}
\description{Access the public Rotterdam breast cancer cohort from the installed \pkg{survival} package, prepared for the IWB example. No separate download is required.}
\usage{rotterdam_example()}
\details{Death is the event; recurrence is not used. Follow-up is converted from days to years by dividing by 365.25, with a minimum of one day, matching the application preprocessing. The example fits after administrative censoring at five years and extrapolates to fifteen years. Original follow-up remains available for retrospective benchmarking.}
\value{A data frame with 2982 rows and two columns: \code{time} (years to death or censoring) and \code{status} (1 for death, 0 for censoring).}
\source{\code{survival::rotterdam}; see \code{\link[survival]{rotterdam}} for the original data documentation and references.}
\examples{
d <- rotterdam_example()
head(d)
}
