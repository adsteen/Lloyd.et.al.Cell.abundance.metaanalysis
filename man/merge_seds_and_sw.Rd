\name{merge_seds_and_sw}
\alias{merge_seds_and_sw}
\title{Merges the sediments and sw data frames into a single allData data frame}
\usage{
  merge_seds_and_sw(corrected_seds, corrected_sw)
}
\arguments{
  \item{corrected_seds}{The 'corrected' sedimetents data
  frame}

  \item{corrected_sw}{The 'corrected' seawater data frame}
}
\value{
  Returns a merged data frame with the name of the
  environment (sediments or seawater) in the `environment`
  column
}
\description{
  Merges the sediments and sw data frames into a single
  allData data frame
}
\details{
  See `obsolete_raw_data_processing_script.R` to see how we
  went from our original .xlsx files to the 'corrected'
  frames. This process involved very minimal processing,
  like setting depth and cell count values of zero to a
  very low number (e.g. 1 mm or 1 cell/ml) to avoid
  problems with log scales
}

