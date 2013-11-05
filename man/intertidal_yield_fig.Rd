\name{intertidal_yield_fig}
\alias{intertidal_yield_fig}
\title{Makes a supplemental figure of FISH yield in intertidal sediments}
\usage{
  intertidal_yield_fig(corrected_seds, corrected_sw,
    yaxsize = 7, colors = brewer.pal(n = 5, name = "RdBu"))
}
\arguments{
  \item{corrected_seds}{Data frame of sediments data}

  \item{corrected_sw}{Data frame of seawater data}

  \item{yaxsize}{Font size for yaxis font}

  \item{colors}{Vector of colors for box fill, should have
  length 5}
}
\description{
  Makes a figure of yields of FISH/CARD-FISH for intertidal
  sediments, analogous to Fig 2
}

