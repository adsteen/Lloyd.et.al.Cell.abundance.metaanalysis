\name{AIC_lik}
\alias{AIC_lik}
\title{Returns liklihood of models relative to the best (based on AIC)}
\usage{
  AIC_lik(AIC_df)
}
\arguments{
  \item{AIC_df}{data frame of AIC results, as provided by
  the AIC function}
}
\value{
  The input data frame, with a `liklihood` column
}
\description{
  Calculates relative liklihood of models based on Akaike
  Information Criterion
}

