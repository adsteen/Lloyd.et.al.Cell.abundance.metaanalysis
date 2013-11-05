\name{reproduce_research}
\alias{reproduce_research}
\title{Reproduces analysis from Lloyd et al (submitted)}
\usage{
  reproduce_research(print_plots = TRUE,
    save_plots = FALSE, reload_data = FALSE,
    fast_calc = TRUE, n_reps = 1000)
}
\arguments{
  \item{print_plots}{If true, print the plots to screen}

  \item{save_plots}{If true, save plots. There is currently
  no option to change the filenames of the saved plots.}

  \item{reload_data}{If true, reload the data from .xlsx
  files. Otherwise read data from .Rdata files. This is
  passed along to `read_data`.}

  \item{fast_calc}{If true, skip the bootstrapped parts of
  the analysis, which take the lion's share of
  computational time.}

  \item{n_reps}{Number of reps for bootstrapped analyses.}
}
\value{
  Returns the three central data frames used in the
  analysis: `all_data`, `corrected_seds`, and
  `corrected_sw`, as named elements in a list. Note that
  `corrected_seds` and `corrected_sw` are, collectively,
  redundant with `all_data`
}
\description{
  Reproduces entire analysis from Lloyd et al, 2013,
  Applied and Environmental Microbiology.
}

