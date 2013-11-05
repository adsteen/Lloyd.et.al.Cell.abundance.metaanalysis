##' Makes a supplemental figure of FISH yield in intertidal sediments
##' 
##' @description Makes a figure of yields of FISH/CARD-FISH for intertidal sediments, analogous to Fig 2
##' @param corrected_seds Data frame of sediments data
##' @param corrected_sw Data frame of seawater data
##' @param yaxsize Font size for yaxis font
##' @param colors Vector of colors for box fill, should have length 5
##' @export

intertidal_yield_fig <- function(corrected_seds, 
                                 corrected_sw, 
                                 yaxsize=7, 
                                 colors=brewer.pal(n=5, name="RdBu")) {
  
  #browser()
  # Set up a data frame containing only samples from the intertidal
  intertidal_yields <- corrected_seds[!is.na(corrected_seds$Arc.permeabilization) &
                                  !is.na(corrected_seds$FISH.yield) &
                                  !is.na(corrected_seds$Fish.or.cardFish) &
                                  #!corrected_seds$Arc.permeabilization == "CARD-FISH, archaea not measured" &
                                  #!corrected_seds$Arc.permeabilization == "FISH, archaea not measured" &
                                  (corrected_seds$Fish.or.cardFish == "CARDFISH" | corrected_seds$Fish.or.cardFish=="FISH") &
                                  corrected_seds$Environment.Type == "Intertidal", ]
  
  
  # Make a permeabilization column, to reflect that FISH doesn't need a permeabilization method
  intertidal_yields$perm <- as.character(intertidal_yields$Arc.permeabilization)
  intertidal_yields$perm[intertidal_yields$perm=="none"] <- "FISH"
  
  # Re-level the permeabilization methods
  intertidal_yields$perm <- factor(intertidal_yields$perm, 
                             levels=c("proteinase K", 
                                      "FISH", 
                                      "FISH, archaea not measured",
                                      "CARD-FISH, archaea not measured",
                                      #"detergent", 
                                      #"unknown",  
                                      #"lysozyme",
                                      "lysozyme/achromopeptidase"), 
                             #labels=c("proteinase K", "none (FISH)", "detergent", 
                             #         "unknown", "lysozyme/achromopeptidase", "lysozyme"),
                             ordered=TRUE)
  
  ### FIg 2a reduced: cut out all the cores with only one data point
  # Count the # points in the core, remove all < 1
  intertidal_yields_no_sing <- ddply(intertidal_yields, .(core), transform, count=length(FISH.yield))
  intertidal_yields_no_sing <- intertidal_yields_no_sing[intertidal_yields_no_sing$count > 1, ]
  intertidal_yields_no_sing$core <- as.factor(as.character(intertidal_yields_no_sing$core))
  
  # Order the papers in order of decreasing median yield
  intertidal_yieldsMedianYields <- ddply(intertidal_yields_no_sing, .(core), summarise, medianYield=median(FISH.yield, na.rm=TRUE))
  intertidal_yieldsMedianYields <- intertidal_yieldsMedianYields[order(intertidal_yieldsMedianYields$medianYield, decreasing=TRUE), ]
  intertidal_yields_no_sing$core <- factor(intertidal_yields_no_sing$core, levels=intertidal_yieldsMedianYields$core, ordered=TRUE)
  
  
  # Determine the quartiles of yield from SW data; overplot that on sediments plot
  #swQuant <- quantile(all_data$FISH.yield[all_data$environment=="seawater"], probs=c(0.5-0.341, 0.5, 0.5+0.341), na.rm=TRUE)
  swQuant <- quantile(corrected_sw$FISH.yield, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
  swYieldDF <- data.frame(ymin=swQuant[1], ymed=swQuant[2], ymax=swQuant[3])
  rownames(swYieldDF) <- NULL
  swYieldDF <- data.frame(x=c(0, length(unique(intertidal_yields_no_sing$core))+0.5), rbind(swYieldDF, swYieldDF))
  
  ## Count number of points in each core
  nPoints <- ddply(intertidal_yields_no_sing, .(core), function(x) data.frame(nrow=nrow(x), perm=x$perm[1]))
  nPoints$core <- factor(nPoints$core, levels=levels(intertidal_yields_no_sing$core), ordered=TRUE)
  nPoints$perm <- factor(nPoints$perm, levels=levels(intertidal_yields_no_sing$perm), ordered=TRUE)
  
  envDF <- ddply(intertidal_yields_no_sing, .(core), function(x) as.character(x$Environment.Type[1]))
  
  # Set the text position
  colors <- brewer.pal(n=6, name="RdBu")
  colors[4] <- "#999999"
  textPos <- 1.4
  #levels(intertidal_yields_no_sing$perm)
  
  #browser()
  fig2aS_2 <- ggplot() +
    geom_blank(data=intertidal_yields_no_sing, aes(x=core, y=FISH.yield, fill=perm)) +
    geom_ribbon(data=swYieldDF, aes(x=x, ymin=ymin, ymax=ymax), fill="skyblue", alpha=1) +
    geom_line(data=swYieldDF, aes(x=x, y=ymed), colour="black") +
    geom_boxplot(data=intertidal_yields_no_sing, aes(x=core, y=FISH.yield, fill=perm), colour="black", outlier.size=1) +
    #geom_rect(data=intertidal_yields_no_sing, aes(xmin=as.numeric(core)-0.5, xmax=as.numeric(core)+0.5, ymin=textPos-0.05, ymax=textPos+0.05, fill=perm)) +
    geom_rect(data=intertidal_yields_no_sing, aes(xmin=as.numeric(core)-0.5, xmax=as.numeric(core)+0.5, ymin=1.4-0.05, ymax=1.4+0.05, fill=perm)) +
    geom_text(data=nPoints, aes(x=core, y=1.4, label=nrow), size=2) +
    geom_text(data=envDF, aes(x=core, y=1.4-0.075, label=V1), size=2, angle=-90, hjust=0) +
    #geom_text(data=nPoints, aes(x=core, y=textPos, label=nrow), size=2) +
    #geom_text(data=envDF, aes(x=core, y=textPos-0.075, label=V1), size=2, angle=-90, hjust=0) +
    ylab(expression(frac("Bacteria + Archaea (FISH or CARD-FISH)", "Total Cell Count"))) +
    coord_cartesian(ylim=c(0,1.6)) + 
    scale_fill_manual(values=colors, name="Archaeal\nPermeabilization\nMethod") +
    theme(axis.text.x=element_text(angle=-45, hjust=0, size=6),
          axis.title.x=element_blank(),
          legend.position="top",
          legend.title.align=1,
          plot.margin=unit(c(0, 0.9, 0, 0), "in")) #default margin is c(1, 1, 0.5, 0.5) 
  
  #ggsave("plots/fig2aS_intertidal_only_yield_tosubmit.png", height=5, width=clm2, dpi=300)
  fig2aS_2
}