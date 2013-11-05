##' Make a plot of totalcells vs CARDFISH (or other FISH) counts in seawater & sediments
##' 
##' @param allData The data frame of all data from merge_seds_and_sw.R
##' @param point_size The size of points on the printed plot
##' @return Returns the ggplot object 
##' @export
##' 
plot_cell_vs_fish <- function(allData, point_size=0.5) {
  # Fig 1 excludes polyribonucleotide-FISH
  fig1data <- allData[!is.na(allData$totalcells) & 
                        !is.na(allData$paper) &
                        !is.na(allData$CARDFISH.Total.per.cc) &
                        !is.na(allData$Fish.or.cardFish) & 
                        !is.na(allData$Arc.permeabilization), ]
  
  # Create vector of permeabilization method and/or FISH 
  fig1data$permForFig1 <- NA
  fig1data$permForFig1[fig1data$Fish.or.cardFish == "FISH"] <- "FISH"
  fig1data$permForFig1[fig1data$Fish.or.cardFish == "Polyribonucleotide FISH"] <- "Polyribonucleotide_FISH"
  fig1data$permForFig1[fig1data$Fish.or.cardFish == "CARDFISH" & fig1data$Arc.permeabilization == "proteinase K"] <- "CARD_w_protK"
  fig1data$permForFig1[fig1data$Fish.or.cardFish == "CARDFISH" & fig1data$Arc.permeabilization != "proteinase K"] <- "CARD_other_perm"
  fig1data <- fig1data[fig1data$totalcells!=0 &
                         fig1data$CARDFISH.Total.per.cc != 0, ]
  
  fig1data$permForFig1 <- factor(fig1data$permForFig1, 
                                 levels=c("FISH", "CARD_w_protK", "CARD_other_perm", "Polyribonucleotide_FISH"), ordered=TRUE)
  
  # Create a 1:1 line
  # oneOne <- data.frame(x=c(1e4, 1e11), y=c(1e4, 1e11))
  
  # How many data points from each method and environment
  nPoints <- ddply(fig1data, c("environment", "permForFig1"), summarise, nPoints=length(paper))
  nPoints$x <- 1e4
  #nPoints$y <- rep(c(1e10, 1e9, 1e8), 2)
  nPoints$y <- c(1e10, 1e9, 1e8, 1e7, 1e10, 1e9, 1e8)
  nPoints$label <- paste("n=", nPoints$nPoints, sep="")
  
#   fig1col <- brewer.pal(n=3, name="Paired")
#   fig1col[3] <- brewer.pal(n=3, name="OrRd")[3]
#   fig1col <- fig1col[c(1, 3, 2)]
  fig1col <- brewer.pal(n=4, name="Paired")
  fig1col[3] <- brewer.pal(n=3, name="OrRd")[3]
  fig1col[1:3] <- fig1col[c(1, 3, 2)] 
  
  # Build the plot
  p_cell_vs_fish <- ggplot(fig1data) + 
    geom_point(aes(x=totalcells, y=CARDFISH.Total.per.cc, colour=permForFig1), size=point_size, alpha=1, guide=FALSE) +
    geom_abline(slope=1, intercept=0, colour="black") +
    geom_text(data=nPoints, aes(x=x, y=y, label=label, colour=permForFig1), size=2, hjust=0, guide=FALSE) +
    scale_x_log10(expression(paste("Total Cells (cells ", ml^{-1}, ")")), breaks=10^(4:11)) + 
    scale_y_log10(expression(paste("FISH or CARD-FISH\nSum of Bac + Arc (cells ", ml^{-1}, ")")), breaks=10^(4:11)) +
    scale_colour_manual(values=fig1col, 
                        name="Method",
                        breaks=c("FISH", "CARD_w_protK", "CARD_other_perm", "Polyribonucleotide_FISH"),
                        labels=c("FISH", "CARD-FISH,\nproteinase K", "CARD-FISH,\nother permeabilization", "Polyribonucleotide\nFISH")
    ) +
    coord_fixed() +
    facet_wrap(~environment) +
    guides(colour = guide_legend(override.aes = list(size=4, alpha=1), nrow=2)) +
    theme(legend.position="top",
          axis.title.y=element_text(hjust=0.5),
          axis.text.x=element_text(angle=-45, hjust=0))
  
  p_cell_vs_fish
}