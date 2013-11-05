\name{lm_stats}
\alias{lm_stats}
\title{Returns statistics for linear regressions (e.g. slope, intercept, error, etc)}
\usage{
  lm_stats(d, xvar, yvar)
}
\arguments{
  \item{d}{Data frame containing at least two numeric
  columns}

  \item{xvar}{Numeric column to be used as the independant
  variable for the regression}

  \item{yvar}{Numeric column to be used as the dependant
  variable for the regression}
}
\value{
  Returns a one-row data frame containing (at present):
  slope, intercept, slope standard error, intercept
  standard error, p value, r-squared, and number of points.
}
\description{
  Given data frame d, returns a data frame of the various
  regression statistics
}
\details{
  include some details here
}
\examples{
require(ggplot2)
lm_stats(mpg, xvar="cty", yvar="hwy")
}

