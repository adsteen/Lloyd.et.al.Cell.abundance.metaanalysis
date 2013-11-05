\name{boxplots_by_perm}
\alias{boxplots_by_perm}
\title{Creates boxplots of yield and %Archaea by permeabilization method
Also runs statistics on those plots}
\usage{
  boxplots_by_perm(corrected_seds, sw_quant, yaxsize = 7,
    colors = brewer.pal(n = 5, name = "RdBu"),
    print_plot = FALSE, save_plot = FALSE, height = 4,
    width = clm, res = 900)
}
\arguments{
  \item{corrected_seds}{Data frame of sediments data}

  \item{sw_quant}{Vector of 25% percentile, median, and 75%
  percentile of all sw yields}

  \item{yaxsize}{Font size of y axis text}

  \item{colors}{Vector of colors for box fills, should be
  length 5}

  \item{print_plot}{If true, print the plot}

  \item{save_plot}{If true, save the plot}

  \item{height}{Height of the plot, in inches, if plot is
  saved}

  \item{width}{Width of the plot, in inches, if plot is
  saved}

  \item{res}{Resolution of the plot, if saved, in dpi}
}
\description{
  Creates plots in Fig 3, and calculates relevant
  statistics
}

