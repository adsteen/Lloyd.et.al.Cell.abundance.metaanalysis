##' Create a boxplot of yield (total cells by *-FISH relative to total cells by direct count) for each core in the database
##' 
##' @param corected_seds The corrected_seds dataframe containing all the sediments data
##' @param sw_quant Vector of median and interquartile range if seawater yield
##' @return The plot (a ggplot2 object)
##' @export


#yield_by_core <- function(corrected_seds, height=5, width=7.08, res=900, print_plot=FALSE, save_plot=FALSE) {
yield_by_core <- function(corrected_seds, 
                          sw_quant,  
                          yield_label="yield", 
                          yaxsize=7,
                          colors=brewer.pal(n=5, name="RdBu")) {
  
  # Create a smaller data frame, excluding irrelevant data (e.g. missing yield information)
  yield_by_core_data <- corrected_seds[!is.na(corrected_seds$Arc.permeabilization) &
                                 !is.na(corrected_seds$FISH.yield) &
                                 !is.na(corrected_seds$Fish.or.cardFish) &
                                 !corrected_seds$Arc.permeabilization == "CARD-FISH, archaea not measured" &
                                 !corrected_seds$Arc.permeabilization == "FISH, archaea not measured" &
                                 (corrected_seds$Fish.or.cardFish == "CARDFISH" | corrected_seds$Fish.or.cardFish=="FISH") &
                                 corrected_seds$Environment.Type != "Intertidal" &
                                 corrected_seds$Environment.Type != "Salt marsh", ]
  
  # Make a permeabilization column, to reflect that FISH doesn't need a permeabilization method
  yield_by_core_data$perm <- as.character(yield_by_core_data$Arc.permeabilization)
  yield_by_core_data$perm[yield_by_core_data$Fish.or.cardFish=="FISH"] <- "FISH"
  
  # Re-level the permeabilization methods
  yield_by_core_data$perm <- factor(yield_by_core_data$perm, 
                            levels=c("proteinase K", "FISH", 
                                     "detergent", "unknown",  
                                     "lysozyme/achromopeptidase", "lysozyme"), 
                            labels=c("proteinase K", "none (FISH)", "detergent", 
                                     "unknown", "lysozyme/\nachromopeptidase", "lysozyme"),
                            ordered=TRUE)
  
  ### FIg 2a reduced: cut out all the cores with only one data point
  # Count the # points in the core, remove all < 1
  yield_by_core_data2 <- ddply(yield_by_core_data, .(core), transform, count=length(FISH.yield))
  yield_by_core_data2 <- yield_by_core_data2[yield_by_core_data2$count > 1, ]
  yield_by_core_data2$core <- as.factor(as.character(yield_by_core_data2$core))
  
  # Order the papers in order of decreasing median yield
  yield_by_core_dataMedianYields <- ddply(yield_by_core_data2, .(core), summarise, medianYield=median(FISH.yield, na.rm=TRUE))
  yield_by_core_dataMedianYields <- yield_by_core_dataMedianYields[order(yield_by_core_dataMedianYields$medianYield, decreasing=TRUE), ]
  yield_by_core_data2$core <- factor(yield_by_core_data2$core, levels=yield_by_core_dataMedianYields$core, ordered=TRUE)
  
  
  # Determine the quartiles of yield from SW data; overplot that on sediments plot
  swYieldDF <- data.frame(ymin=sw_quant[1], ymed=sw_quant[2], ymax=sw_quant[3])
  rownames(swYieldDF) <- NULL
  swYieldDF <- data.frame(x=c(0, length(unique(yield_by_core_data2$core))+0.5), rbind(swYieldDF, swYieldDF))
  
  ## Count number of points in each core
  nPoints <- ddply(yield_by_core_data2, .(core), function(x) data.frame(nrow=nrow(x), perm=x$perm[1]))
  nPoints$core <- factor(nPoints$core, levels=levels(yield_by_core_data2$core), ordered=TRUE)
  nPoints$perm <- factor(nPoints$perm, levels=levels(yield_by_core_data2$perm), ordered=TRUE)
  
  envDF <- ddply(yield_by_core_data2, .(core), function(x) as.character(x$Environment.Type[1]))
  
  p_yield_by_core <- ggplot() +
    geom_blank(data=yield_by_core_data2, aes(x=core, y=FISH.yield, fill=perm)) +
    geom_ribbon(data=swYieldDF, aes(x=x, ymin=ymin, ymax=ymax), fill="skyblue", alpha=1) +
    geom_line(data=swYieldDF, aes(x=x, y=ymed), colour="black") +
    geom_boxplot(data=yield_by_core_data2, aes(x=core, y=FISH.yield, fill=perm), colour="black", outlier.size=1) +
    geom_rect(data=yield_by_core_data2, aes(xmin=as.numeric(core)-0.5, xmax=as.numeric(core)+0.5, ymin=1.5-0.05, ymax=1.5+0.05, fill=perm)) +
    geom_text(data=nPoints, aes(x=core, y=1.5, label=nrow), size=2) +
    geom_text(data=envDF, aes(x=core, y=1.5-0.075, label=V1), size=2, angle=-90, hjust=0) +
    ylab(yield_label) +
    coord_cartesian(ylim=c(0,1.6)) + 
    scale_fill_manual(values=colors, name="Archaeal\nPermeabilization Method") +
    theme(axis.text.x=element_text(angle=-45, hjust=0, size=6),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=yaxsize),
          legend.position="top",
          legend.title.align=1,
          plot.margin=unit(c(0, 0.9, 0, 0), "in")) #default margin is c(1, 1, 0.5, 0.5) 
#   
#   
#   if (print_plot) {
#     print(fig2a)
#   }
#   if (save_plot) {
#     if (is.na(fn)) {
#       fn <- "plot_of_yield_by_core.tiff"
#     }
#     tiff(fn, height=height, width=width, units="in", res=myDPI, compression="lzw", type="cairo")
#     print(fig2a)
#     dev.off()
#   }
  p_yield_by_core
  
}