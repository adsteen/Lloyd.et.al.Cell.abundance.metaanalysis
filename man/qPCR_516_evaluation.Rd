\name{qPCR_516_evaluation}
\alias{qPCR_516_evaluation}
\title{Makes plot and calculates summary stats for qPCR using or not using 516 as a primer}
\usage{
  qPCR_516_evaluation(corrected_seds)
}
\arguments{
  \item{corrected_seds}{data frame of corrected seds}
}
\value{
  Returns a list with the following components:

  * `p_516`, the plot object * med_and_IQR, the data frame
  of median yield and interquartile range for each depth
  and primer set (uses or doesn't use 516) * t_test_516 An
  object of class `htest`, containing the results of the
  student's t test of yield for samples using vs not using
  516 (all depths) * n_points The number of data points for
  each combination of depth and uses/doesn't use 516
}
\description{
  Makes plot and calculates summary stats for qPCR using or
  not using 516 as a primer
}

