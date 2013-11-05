##' Creates boxplots of yield and %Archaea by permeabilization method
##' Also runs statistics on those plots
##'
##'@description Creates plots in Fig 3, and calculates relevant statistics
##'@param corrected_seds Data frame of sediments data
##'@param sw_quant Vector of 25% percentile, median, and 75% percentile of all sw yields
##'@param yaxsize Font size of y axis text
##'@param colors Vector of colors for box fills, should be length 5
##'@param print_plot If true, print the plot
##'@param save_plot If true, save the plot
##'@param height Height of the plot, in inches, if plot is saved
##'@param width Width of the plot, in inches, if plot is saved
##'@param res Resolution of the plot, if saved, in dpi
##'@export

boxplots_by_perm <- function(corrected_seds, 
                             sw_quant, 
                             yaxsize=7,
                             colors=brewer.pal(n=5, name="RdBu"),
                             print_plot=FALSE, save_plot=FALSE, 
                             height=4, width=clm, res=900) {
  
  # Colors to indicate permeabilization method
  box_col <- colors[c(1, 4, 3, 5)]
  
  # Y-axis labels
  percent_arc_label <- expression(frac("Archaea (CARD-FISH)", "Bacteria + Archaea (CARD-FISH)"))
  yield_label <- expression(frac("Bacteria + Archaea (CARD-FISH)", "Total Cell Count"))
  
  
  # Want to replace with boxplot of all cardfish % Archaea + Fish, colored by permeabilization method, for sediments only
  p_percent_arc_data <- corrected_seds[!is.na(corrected_seds$Arc.permeabilization) &
                                 !is.na(corrected_seds$Fish.or.cardFish) &
                                 !is.na(corrected_seds$Fraction.Arc.CARDFISH) &
                                 !corrected_seds$Arc.permeabilization == "unknown" &
                                 !corrected_seds$Arc.permeabilization == "not measured"  &
                                 !corrected_seds$Fish.or.cardFish == "FISH" &
                                 corrected_seds$Environment.Type != "Intertidal" &
                                 corrected_seds$Environment.Type != "Salt marsh", 
                               c("depth", "Arc.permeabilization", "Fraction.Arc.CARDFISH", "paper")]
  
  # Re-order the factor for permeabilization method
  p_percent_arc_data$perm <- as.character(p_percent_arc_data$Arc.permeabilization)
  p_percent_arc_data$perm <- factor(p_percent_arc_data$perm, 
                            levels=c("proteinase K", "detergent", 
                                     "lysozyme/achromopeptidase", "lysozyme"),
                            labels=c("proteinase K", "detergent", 
                                     "lysozyme/\nachromopeptidase", "lysozyme"),
                            ordered=TRUE)
  
  p_percent_arc_medians <- ddply(p_percent_arc_data, .(perm), summarise, med=median(Fraction.Arc.CARDFISH))
  p_percent_arc_data$perm <- factor(p_percent_arc_data$perm, levels=p_percent_arc_medians$perm[order(p_percent_arc_medians$med, decreasing=TRUE)])
  
  pointsAndStudies <- ddply(p_percent_arc_data, .(perm),
                            function(x) data.frame(nPoints=nrow(x), 
                                                   nStudies=length(unique(x$paper))))
  pointsAndStudies$label <- paste(pointsAndStudies$nPoints, " Points\n", pointsAndStudies$nStudies, " Studies", sep="")
  
  # Add the letters to a data frame in order to display them
  pointsAndStudies$letter <- NA
  pointsAndStudies[pointsAndStudies$perm=="proteinase K", "letter"] <- "a"
  pointsAndStudies[pointsAndStudies$perm=="lysozyme/\nachromopeptidase", "letter"] <- "b"
  pointsAndStudies[pointsAndStudies$perm=="detergent", "letter"] <- "b"
  pointsAndStudies[pointsAndStudies$perm=="lysozyme", "letter"] <- "c"
  
  
  # Plot of % Archaea by permeabilization method, all depths
  p_percent_arc_all_depths <- ggplot(p_percent_arc_data) + 
    geom_boxplot(data=p_percent_arc_data, aes(x=perm, y=Fraction.Arc.CARDFISH, fill=perm), outlier.size=1) +
    geom_text(data=pointsAndStudies, aes(x=perm, y=1.2, label=label), vjust=1, size=2) +
    geom_text(data=pointsAndStudies, aes(x=perm, y=-0.03, label=letter), vjust=1, size=2.5) +
    scale_y_continuous(percent_arc_label, limits=c(-0.1, 1.25), breaks=seq(from=0, to=1, by=0.25)) +
    scale_x_discrete("All Data") +
    scale_fill_manual(values=box_col, name="permeabilization\nmethod") +
    theme(legend.position="none",
          text=element_text(size=8), 
          axis.text.x=element_text(angle=0, hjust=0.5),
          axis.title.y=element_text(size=yaxsize))#,
          #plot.margin=unit(fig2bc_spacing, "in")) 
  
  # Calculate % Archaea for all data
  percent_arc_all_depths <- ddply(p_percent_arc_data, .(perm), summarise,
                                  median=quantile(Fraction.Arc.CARDFISH, probs=0.50, na.rm=TRUE), 
                                  twenty.fifth.percentile=quantile(Fraction.Arc.CARDFISH, probs=0.25, na.rm=TRUE),
                                  seventy.fifth.percentile=quantile(Fraction.Arc.CARDFISH, probs=0.75, na.rm=TRUE))
  # Create a data frame to hold all data and < 1 m depth
  percent_arc_all_depths$quantity <- "fraction archaea"
  percent_arc_all_depths$dataset <- "all depths"
  
  #=========
  # Fig 2c: boxplots of percent ARC for only seds < 1 m deep
  #=========
  
  p_percent_arc_shallow_data <- p_percent_arc_data[p_percent_arc_data$depth <= 1, ]
  
  
  # Calculate the number of data points and studies
  pointsAndStudies1m <- ddply(p_percent_arc_shallow_data, .(perm),
                              function(x) data.frame(nPoints=nrow(x), 
                                                     nStudies=length(unique(x$paper))))
  pointsAndStudies1m$label <- paste(pointsAndStudies1m$nPoints, " Points\n", pointsAndStudies1m$nStudies, " Studies", sep="")
  
  # Add letter labels for signifciantly different permeabilization methods
  pointsAndStudies1m$letter <- NA
  pointsAndStudies1m[pointsAndStudies1m$perm=="proteinase K", "letter"] <- "a"
  pointsAndStudies1m[pointsAndStudies1m$perm=="lysozyme/\nachromopeptidase", "letter"] <- "b"
  pointsAndStudies1m[pointsAndStudies1m$perm=="detergent", "letter"] <- "b"
  pointsAndStudies1m[pointsAndStudies1m$perm=="lysozyme", "letter"] <- "b"
  
  
  
  
  p_percent_arc_shallow <- ggplot(p_percent_arc_shallow_data) + 
    geom_boxplot(data=p_percent_arc_shallow_data, aes(x=perm, y=Fraction.Arc.CARDFISH, fill=perm), outlier.size=1) +
    geom_text(data=pointsAndStudies1m, aes(x=perm, y=1.2, label=label), vjust=1, size=2) +
    geom_text(data=pointsAndStudies1m, aes(x=perm, y=-0.03, label=letter), vjust=1, size=2.5) +
    scale_y_continuous(percent_arc_label,limits=c(-0.1, 1.25), breaks=seq(from=0, to=1, by=0.25)) +
    scale_x_discrete("Shallower than 1 m") +
    scale_fill_manual(values=box_col) +
    theme(legend.position="none",
          text=element_text(size=8), 
          axis.text.x=element_text(angle=0, hjust=0.5),
          axis.text.x=element_blank(),
          axis.title.y=element_text(size=yaxsize))#,
          #plot.margin=unit(fig2bc_spacing, "in")) 
  
  # Calculate % Archaea for data < 1m deep
  percent_arc_shallow <- ddply(p_percent_arc_shallow_data, .(perm), summarise,
                                  median=quantile(Fraction.Arc.CARDFISH, probs=0.50, na.rm=TRUE), 
                                  twenty.fifth.percentile=quantile(Fraction.Arc.CARDFISH, probs=0.25, na.rm=TRUE),
                                  seventy.fifth.percentile=quantile(Fraction.Arc.CARDFISH, probs=0.75, na.rm=TRUE))
  percent_arc_shallow$quantity <- "fraction archaea"
  percent_arc_shallow$dataset <- "< 1 m"
  
  #====================
  # Fig 2d and 2e: same as Fig 2b and c, but with yield
  #====================
  
  p_yield_data <- corrected_seds[!is.na(corrected_seds$Arc.permeabilization) &
                                               !is.na(corrected_seds$Fish.or.cardFish) &
                                               !is.na(corrected_seds$FISH.yield) &
                                               !corrected_seds$Arc.permeabilization=="unknown" &
                                               !corrected_seds$Arc.permeabilization=="not measured"  &
                                               !corrected_seds$Fish.or.cardFish=="FISH" &
                                               corrected_seds$Environment.Type != "Intertidal" &
                                               corrected_seds$Environment.Type != "Salt marsh", ]
  
  #p_yield_data <- fig2b_data
  p_yield_data$perm <- as.character(p_yield_data$Arc.permeabilization)
  #fig2b_data$perm[fig2b_data$Fish.or.cardFish=="FISH"] <- "FISH"
  p_yield_data$perm <- factor(p_yield_data$perm, 
                            levels=c("proteinase K", "lysozyme/achromopeptidase", "detergent", "lysozyme"), 
                            labels=c("proteinase K", "lysozyme/\nachromopeptidase", "detergent", "lysozyme"),
                            ordered=TRUE)
  
  p_yield_shallow_data <- p_yield_data[!is.na(p_yield_data$depth) &
                             p_yield_data$depth <= 1, ]
  
  # proteinase K gets A
  # lysozyme-achromopeptidase gets B
  # detergent gets B
  # lysozyme gets C
  pointsAndStudiesYield <- ddply(p_yield_data, .(perm),
                                 function(x) data.frame(nPoints=nrow(x), 
                                                        nStudies=length(unique(x$paper))))
  
  pointsAndStudiesYield$label <- paste(pointsAndStudies$nPoints, " Points\n", pointsAndStudies$nStudies, " Studies", sep="")
  pointsAndStudiesYield$letter <- NA
  pointsAndStudiesYield[pointsAndStudiesYield$perm=="proteinase K", "letter"] <- "a"
  pointsAndStudiesYield[pointsAndStudiesYield$perm=="lysozyme/\nachromopeptidase", "letter"] <- "b"
  pointsAndStudiesYield[pointsAndStudiesYield$perm=="detergent", "letter"] <- "b"
  pointsAndStudiesYield[pointsAndStudiesYield$perm=="lysozyme", "letter"] <- "c"
  
  swDF <- data.frame(perm=c(0, 1, 2, 5), #perm = levels(p_yield_data$perm), #
                     ymin=rep(sw_quant[1], nlevels(p_yield_data$perm)), 
                     ymed=rep(sw_quant[2], nlevels(p_yield_data$perm)),
                     ymax=rep(sw_quant[3], nlevels(p_yield_data$perm)))
  
  p_yield_all <- ggplot(p_yield_data) + 
    geom_blank(data=p_yield_data, aes(x=perm, y=FISH.yield, fill=perm)) +
    geom_ribbon(data=swDF, aes(x=perm, ymin=ymin, ymax=ymax), fill="skyblue") +
    geom_line(data=swDF, aes(x=perm, y=ymed)) +
    geom_boxplot(data=p_yield_data, aes(x=perm, y=FISH.yield, fill=perm), outlier.size=1) +
    geom_text(data=pointsAndStudiesYield, aes(x=perm, y=1.4, label=label), vjust=1, size=2) +
    geom_text(data=pointsAndStudiesYield, aes(x=perm, y=0, label=letter), vjust=1, size=2.5) +
    #scale_y_continuous("yield", limits=c(-0.1, 1.25)) +
    coord_cartesian(xlim=c(0.25, 4.75), ylim=c(-0.1, 1.5)) +
    xlab("All data") +
    ylab(yield_label) +
    scale_fill_manual(values=box_col, name="permeabilization\nmethod") +
    theme(legend.position="none",
          text=element_text(size=8), 
          axis.text.x=element_text(angle=0, hjust=0.5),
          axis.text.x=element_blank(),
          axis.title.y=element_text(size=yaxsize),
          #plot.margin=unit(fig2bc_spacing, "in"),
          axis.title.x=element_blank())  
  
  # Calculate median and IQR of yield for all data
  yield_all_data <- ddply(p_yield_data, .(perm), summarise,
                               median=quantile(FISH.yield, probs=0.50, na.rm=TRUE), 
                               twenty.fifth.percentile=quantile(FISH.yield, probs=0.25, na.rm=TRUE),
                               seventy.fifth.percentile=quantile(FISH.yield, probs=0.75, na.rm=TRUE))
  yield_all_data$quantity <- "yield"
  yield_all_data$dataset <- "all data"
  
  
  
  pointsAndStudiesYield1m <- ddply(p_yield_shallow_data, .(perm),
                                   function(x) data.frame(nPoints=nrow(x), 
                                                          nStudies=length(unique(x$paper))))
  
  pointsAndStudiesYield1m$label <- paste(pointsAndStudies$nPoints, " Points\n", pointsAndStudies$nStudies, " Studies", sep="")
  pointsAndStudiesYield1m$letter <- NA
  pointsAndStudiesYield1m[pointsAndStudiesYield1m$perm=="proteinase K", "letter"] <- "a"
  pointsAndStudiesYield1m[pointsAndStudiesYield1m$perm=="lysozyme/\nachromopeptidase", "letter"] <- "b"
  pointsAndStudiesYield1m[pointsAndStudiesYield1m$perm=="detergent", "letter"] <- "b"
  pointsAndStudiesYield1m[pointsAndStudiesYield1m$perm=="lysozyme", "letter"] <- "c"
  
  
  p_yield_shallow <- ggplot(p_yield_shallow_data) + 
    geom_blank(data=p_yield_shallow_data, aes(x=perm, y=FISH.yield, fill=perm)) +
    geom_ribbon(data=swDF, aes(x=perm, ymin=ymin, ymax=ymax), fill="skyblue") +
    geom_line(data=swDF, aes(x=perm, y=ymed)) +
    geom_boxplot(data=p_yield_shallow_data, aes(x=perm, y=FISH.yield, fill=perm), outlier.size=1) +
    geom_text(data=pointsAndStudies1m, aes(x=perm, y=1.4, label=label), vjust=1, size=2) +
    geom_text(data=pointsAndStudies1m, aes(x=perm, y=0, label=letter), vjust=1, size=2.5) +
    #scale_y_continuous("yield", limits=c(-0.1, 1.25)) +
    coord_cartesian(xlim=c(0.25, 4.75), ylim=c(-0.1, 1.5)) +
    #xlab("Shallower than 1 m") +
    ylab(yield_label) +
    scale_fill_manual(values=box_col, name="permeabilization\nmethod") +
    theme(legend.position="none",
          text=element_text(size=8), 
          axis.text.x=element_text(angle=0, hjust=0.5),
          axis.text.x=element_blank(),
          axis.title.y=element_text(size=yaxsize),
          #plot.margin=unit(fig2bc_spacing, "in"),
          #plot.margin=unit(c(0, 0, 0, 0), "in"),
          axis.title.x=element_blank())  
  
  # Calculate median and IQR opf yield for <1 m data
  yield_shallow <- ddply(p_yield_shallow_data, .(perm), summarise,
                     median=quantile(FISH.yield, probs=0.50, na.rm=TRUE), 
                     twenty.fifth.percentile=quantile(FISH.yield, probs=0.25, na.rm=TRUE),
                     seventy.fifth.percentile=quantile(FISH.yield, probs=0.75, na.rm=TRUE))
  yield_shallow$quantity <- "yield"
  yield_shallow$dataset <- "< 1 m"
  
  # Bind the summary data frames into a single frame
  permeabilization_methods_summary <- rbind(percent_arc_all_depths,
                                            percent_arc_shallow)#,
  permeabilization_methods_summary <- rbind(permeabilization_methods_summary,
                                            yield_all_data)
  permeabilization_methods_summary <- rbind(permeabilization_methods_summary,
                                            yield_shallow)
  permeabilization_methods_summary$perm <- as.character(permeabilization_methods_summary$perm)
  permeabilization_methods_summary$perm[permeabilization_methods_summary$perm=="lysozyme/\nachromopeptidase"] <- "lys/achrom"
                                            
  # Return values
  list(p_percent_arc_all_depths=p_percent_arc_all_depths, 
       p_percent_arc_shallow=p_percent_arc_shallow, 
       p_yield_all=p_yield_all, 
       p_yield_shallow=p_yield_shallow,
       permeabilization_methods_summary=permeabilization_methods_summary)
}