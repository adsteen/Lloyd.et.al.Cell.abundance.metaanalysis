\name{aov_perm_test}
\alias{aov_perm_test}
\title{AOV permutation test}
\usage{
  aov_perm_test(n, d, cat_var, data_var,
    print_hist = FALSE, save_hist = TRUE, title = NULL,
    fn = NA)
}
\arguments{
  \item{n}{Number of permutations}

  \item{d}{Data frame to permute}

  \item{cat_var}{Name of categorical (independant)
  variable, as a string. Should work as a column number as
  well.}

  \item{data_var}{Name of random (dependant) variable, as a
  string. Should work as a column number as well.}

  \item{print_hist}{If `true`, plot the histogram of F
  values.}

  \item{save_hist}{If `true`, save the histogram of F
  values}

  \item{title}{Title of the histogram}

  \item{fn}{Filename, as a string, to save the histogram}
}
\description{
  Performs a permutation test to assess variability among
  groups in data with non-normally-distributed variance
}

