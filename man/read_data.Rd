\name{read_data}
\alias{read_data}
\title{Reads and formats cell abundance database}
\usage{
  read_data(reload.from.xlsx = FALSE, seds_fn,
    seds_sheet_name, sw_fn, sw_sheet_name,
    save_seds_ref_table = TRUE, seds_paper_table_fn = NA,
    all_data_fn = "data/all_data.RData",
    corrected_seds_fn = "data/corrected_seds.RData",
    corrected_sw_fn = "data/corrected_sw.RData")
}
\arguments{
  \item{reload.from.xlsx}{If TRUE, reload & reprocess the
  data from the .csv files. Otherwise, rely on .Rdata file}

  \item{seds_fn}{Filename for sediments data, as .xlsx.
  SHOULD INCLUDE .xlsx extension}

  \item{seds_sheet_name}{Name of worksheet for data in
  sediments file}

  \item{sw_fn}{Filename for seawater data, as .xlsx SHOULD
  INCLUDE .xlsx extension}

  \item{sw_sheet_name}{Name of worksheet for data in
  seawater file}

  \item{save_data}{Whether to save the compiled data for
  future use (as *both* .RData and .csv files)}

  \item{all_data_fn}{Filename for saved data - SHOULD HAVE
  NO EXTENSION as two files will be saved (one .Rdata file,
  one .csv file)}
}
\description{
  Loads the data & applies some basic corrections. If
  reload.from.xlsx is TRUE, will read data from an
  appropriately-formatted xlsx file, which is
  time-consuming. Otherwise, will reload data from an R
  file, which may be the original, package-supplied
  database, or a user-modifed database (much faster). We
  recommend that users who wish to extend the database save
  the database as .RData files, rather than repeatedly
  reloading it from Excel files.
}
\details{
  In order to reproduce the analysis from our original
  database, simply use `invisible(read_data())` - you don't
  have to assign the results to a variable, because
  `read_data()` creates global variables for the
}
\examples{
data_list <- read_data()
corrected_seds <- data_list$corrected_seds
}

