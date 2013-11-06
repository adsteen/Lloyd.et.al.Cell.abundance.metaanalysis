\name{Lloyd.et.al.Cell.abundance.metaanalysis-package}
\alias{Lloyd.et.al.Cell.abundance.metaanalysis-package}
\alias{Lloyd.et.al.Cell.abundance.metaanalysis}
\docType{package}
\title{
Lloyd.et.al.Cell.abundance.metaanalysis
}
\description{
Contains the database underlying the paper, and includes functions to reproduce all the analysis from the paper.
}
\details{
\tabular{ll}{
Package: \tab Lloyd.et.al.Cell.abundance.metaanalysis\cr
Type: \tab Package\cr
Version: \tab 0.9\cr
Date: \tab 2013-11-05\cr
License: \tab MIT License\cr
}
This package contains the date contained in Lloyd et al 
}
\author{
Andrew D. Steen
Maintainer: Andrew D. Steen <andrew.decker.steen@gmail.com>
}
\references{
Lloyd, Karen G., May, Megan T., Kevorkian, Richard T., and Steen, Andrew D. Metaanalysis of quantification methods shows archaea and bacteria to be similarly abundant in the subseafloor. Applied and Environmental Microbiology.
}

\keyword{ package }
\seealso{
}
\examples{
data_list <- read_data()
corrected_seds <- data_list$corrected_seds
corrected_sw <- data_list$corrected_sw
all_data <- data_list$all_data
ggplot(all_data, aes(x=))
}
