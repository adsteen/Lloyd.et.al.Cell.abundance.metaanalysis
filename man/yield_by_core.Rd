\name{yield_by_core}
\alias{yield_by_core}
\title{Create a boxplot of yield (total cells by *-FISH relative to total cells by direct count) for each core in the database}
\usage{
  yield_by_core(corrected_seds, sw_quant,
    yield_label = "yield", yaxsize = 7,
    colors = brewer.pal(n = 5, name = "RdBu"))
}
\arguments{
  \item{corected_seds}{The corrected_seds dataframe
  containing all the sediments data}

  \item{sw_quant}{Vector of median and interquartile range
  if seawater yield}
}
\value{
  The plot (a ggplot2 object)
}
\description{
  Create a boxplot of yield (total cells by *-FISH relative
  to total cells by direct count) for each core in the
  database
}

