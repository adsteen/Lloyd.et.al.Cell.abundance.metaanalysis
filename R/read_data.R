##' Reads and formats cell abundance database
##' 
##' @description Loads the data & applies some basic corrections. If reload.from.xlsx is TRUE, will read data from an appropriately-formatted xlsx file, which is time-consuming. Otherwise, will reload data from an R file, which may be the original, package-supplied database, or a user-modifed database (much faster). We recommend that users who wish to extend the database save the database as .RData files, rather than repeatedly reloading it from Excel files. 
##' @param reload.from.xlsx If TRUE, reload & reprocess the data from the .csv files. Otherwise, rely on .Rdata file
##' @param seds_fn Filename for sediments data, as .xlsx. SHOULD INCLUDE .xlsx extension
##' @param seds_sheet_name Name of worksheet for data in sediments file
##' @param sw_fn Filename for seawater data, as .xlsx SHOULD INCLUDE .xlsx extension
##' @param sw_sheet_name Name of worksheet for data in seawater file
##' @param save_data Whether to save the compiled data for future use (as *both* .RData and .csv files)
##' @param all_data_fn Filename for saved data - SHOULD HAVE NO EXTENSION as two files will be saved (one .Rdata file, one .csv file)
##' @details In order to reproduce the analysis from our original database, simply use `invisible(read_data())` - you don't have to assign the results to a variable, because `read_data()` creates global variables for the 
##' 
##' @export
##' @examples
##' data_list <- read_data()
##' corrected_seds <- data_list$corrected_seds

read_data <- function(reload.from.xlsx=FALSE, seds_fn, seds_sheet_name, 
                      sw_fn, sw_sheet_name, 
                      save_seds_ref_table=TRUE, seds_paper_table_fn=NA,
                      all_data_fn="data/all_data.RData", 
                      corrected_seds_fn="data/corrected_seds.RData", 
                      corrected_sw_fn="data/corrected_sw.RData") {
  

  if (reload.from.xlsx) {
    print("Reading from xlsx files. This could take several minutes.")
    
    #======================
    # Read & munge sediment data
    #=====================
    
    raw_seds <- read.xlsx(seds_fn, sheetName=seds_sheet_name)
    
    #=================
    # Fix column names
    #=================
    names(raw_seds)[names(raw_seds)=="Depth..m."] <- "depth"
    attr(raw_seds$depth, "units") <- "meters"
    
    names(raw_seds)[names(raw_seds)=="Cells.per.cc"] <- "totalcells"
    attr(raw_seds$totalcells, "units") <- "cells per cc"
    
    names(raw_seds)[names(raw_seds)=="qPCR.Total.per.cc"] <- "qPCRtotal"
    attr(raw_seds$qPCRtotal, "units") <- "qPCR total copies per cc"
    
    names(raw_seds)[names(raw_seds)=="Fraction.Arc.qPCR"] <- "percentqPCR"
    attr(raw_seds$percentqPCR, "units") <- "fraction of 1"
    
    names(raw_seds)[names(raw_seds)=="Sulfate..mM."] <- "sulfate"
    attr(raw_seds$sulfate, "units") <- "mM"
    
    names(raw_seds)[names(raw_seds)=="qPCR.with.universal.primers.per.cc"] <- "qPCRuniversal"
    attr(raw_seds$qPCRuniversal, "units") <- "qPCR with universal primers, per cc"
    
    names(raw_seds)[names(raw_seds)=="CARDFISH.Arc.below.detection.limit.or.not.measured."] <- "CARDFISH.Arc.nd"
    
    
    #============
    # Convert columns into the correct data type
    #============
    
    raw_seds$Mud.volcano.or.seep.[raw_seds$Mud.volcano.or.seep.==FALSE | is.na(raw_seds$Mud.volcano.or.seep.)] <- FALSE
    raw_seds$Environment.Type <- as.character(raw_seds$Environment.Type)
    raw_seds$Environment.Type[is.na(raw_seds$Environment.Type)] <- ""
    
    corrected_seds <- raw_seds
    
    names(corrected_seds)[names(corrected_seds)=="Fraction.Arc.qPCR.out.of.bac.plus.arc"] <- "percentqPCR"
    
    #==========
    # Create some useful columns
    #==========
    
    # Convert all zero depths to 1 cm & log-transform the depth
    corrected_seds$depth[(corrected_seds$depth<=0.01)&(is.na(corrected_seds$depth)==FALSE)] <- 0.01
    attr(corrected_seds$depth, "units") <- "m"
    corrected_seds$depth_log10 <- log10(corrected_seds$depth)
    
    # Append a number corresponding to each paper, to indicate the paper the core came from
    corrected_seds$paperNumber <- as.numeric(corrected_seds$paper)
    corrected_seds$core <- paste(as.character(corrected_seds$core), " [", corrected_seds$paperNumber, "]", sep="")
    corrected_seds$Uses.516.for.Arc[is.na(corrected_seds$Uses.516.for.Arc)] <- FALSE
    
    #=============
    # Read & munge seawater data
    #=============
    
    # div/0! and blanks to be read as NA
    na.strings <- c("#DIV/0!", "")
    water_raw_data <- read.xlsx(sw_fn, sheetName=sw_sheet_name)
    corrected_sw <- water_raw_data
    
    names(corrected_sw)[names(corrected_sw)=="Depth..m."] <- "depth"
    
    names(corrected_sw)[names(corrected_sw)=="qPCR.Bacteria..copies.mL.water."] <- "qPCRbac"
    attr(corrected_sw$qPCRbac, "units") <- "copiesPerMl"
    
    names(corrected_sw)[names(corrected_sw)=="qPCR.Archaea..copies.mL.water."] <- "qPCRarc"
    attr(corrected_sw$qPCRarc, "units") <- "copiesPerMl"
    
    names(corrected_sw)[names(corrected_sw)=="Total.qPCR.copies.mL.water."] <- "qPCRtotal"
    attr(corrected_sw$qPCRarc, "units") <- "copiesPerMl"
    
    names(corrected_sw)[names(corrected_sw)=="Cells.per.cc"] <- "totalcells"
    attr(corrected_sw$qPCRarc, "units") <- "cellsPerMl"
    
    #names(corrected_sw)[names(corrected_sw)=="Fraction.Arc.CARDFISH"] <- "percentArcCARD"
    #attr(corrected_sw$qPCRarc, "units") <- "percent"
    
    # Convert all seawater "0" depths less than 1 m (mostly zeros) to 1 m to prevent log problems
    corrected_sw$depth[(corrected_sw$depth<=0)&(is.na(corrected_sw$depth)==FALSE)] <- 1
    
    # Convert all seawater Archaea "0" counts to 1 to prevent log problems
    corrected_sw$CARDFISH.Arc.per.cc[corrected_sw$CARDFISH.Arc.per.cc<=0 & !is.na(corrected_sw$CARDFISH.Arc.per.cc)] <- 1
    
    # Log-transform the depths
    corrected_sw$depth_log10 <- log10(corrected_sw$depth)
    
    # Convert all of the Schattenhoffer permeabilization data to "proteinase K"
    corrected_sw$Arc.permeabilization[corrected_sw$paper=="Schattenhofer 2009"] <- "proteinase K"
  
    # Merge sediments and seawater data (also drop and rename some columns from each data set)
    all_data <- merge_seds_and_sw(corrected_seds=corrected_seds, corrected_sw=corrected_sw)
    
    # Potentially save the data for later, quick loading
    save(all_data, file="data/all_data.RData")
    save(corrected_seds, file="data/corrected_seds.RData")
    save(corrected_sw, file="data/corrected_sw.RData")
    write.csv(all_data, file="data/all_data.csv") 
  }
  
  else {
    ######################
    # Rely on pre-saved versions of the RData file
    ######################
    print("Loading pre-saved data.")
    load(all_data_fn)
    load(corrected_seds_fn)
    load(corrected_sw_fn)
  }
  
  # Strip out bad rows (with NA for paper; these are mostly empty rows that Bill Gates helpfully inserts at the end of the Excel file)
  all_data <- all_data[!is.na(all_data$paper), ]
  corrected_sw <- corrected_sw[!is.na(corrected_sw$paper), ]
  corrected_seds <- corrected_seds[!is.na(corrected_seds$paper), ]
  
  if (save_seds_ref_table) {
    # Save the table of papers & the reference numbers we assign them as a pdf
    seds_paper_table <- ddply(corrected_seds, c("paper", "paperNumber"), function(x) "")
    
    if(is.na(seds_paper_table_fn)) {
      # Automatically create filename if not passed as an argument
      seds_paper_table_fn <- "paper_index_numbers.pdf"
    }
    # Create the pdf
    pdf(seds_paper_table_fn) 
    grid.table(seds_paper_table)
    dev.off()
    
  }
  
  # Return the master data frame
  # Note that this is kind of redundant, since 
  list(all_data=all_data, 
       corrected_seds=corrected_seds, 
       corrected_sw=corrected_sw,
       seds_paper_table=seds_paper_table)
}