\name{calc_cn_ratio}
\alias{calc_cn_ratio}
\title{Calculates the copy number ratio}
\usage{
  calc_cn_ratio(corrected_seds)
}
\arguments{
  \item{corrected_seds}{corrected sediments data frame}
}
\description{
  Calculates the copy number ratio
}
\details{
  Uses the equation cn(a) / cn(b) = (Bq / Ba) * ((cell ct *
  yield)/B(estimated) - 1)
}

