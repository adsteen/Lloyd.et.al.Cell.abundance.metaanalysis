\name{make_sw_yield_boxplots}
\alias{make_sw_yield_boxplots}
\title{Create boxplots of yield in seawater, split by various variables}
\usage{
  make_sw_yield_boxplots(corrected_sw, n_perm = 1000)
}
\arguments{
  \item{corrected_sw}{The corrected seawater boxplot}

  \item{n_perm}{Number of permutations for permutation
  test, passed through to `aov_perm_test`}
}
\description{
  Create boxplots of yield in seawater, split by various
  variables
}
\details{
  Note that this throws some warnings about points being
  omitted. This is simply because the outliers are turned
  off on the boxplots, which is fine because the outliers
  are nevertheless represented in the backgound points
}

