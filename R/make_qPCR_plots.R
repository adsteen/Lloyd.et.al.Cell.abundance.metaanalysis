##' Creates Fig 4, the qPCR plots
##' 
##' @param all_data The `all_data` data frame
##' @param corrected_seds The `corrected_seds` data frame. Yes, these are redundant
##' @param ps The point size for the plots b and c
##' @export
##' 
make_qPCR_plots <- function(all_data, corrected_seds, ps=0.75) {
  
  ##############
  # Plot of sum (Bacterial + Archaeal gene copies) vs gene copies by universal primers
  #############
  
  fig3a_df <- corrected_seds[!is.na(corrected_seds$qPCR.Bac.per.cc) &
                               !is.na(corrected_seds$qPCR.Arc.per.cc) &
                               !is.na(corrected_seds$qPCRuniversal) ,]
  
  oneOne <- data.frame(x=10^(3:9), y=10^(3:9))
  
  # Not colored
  p_qPCR_sum_bac_arc_v_universal <- ggplot(fig3a_df, aes(x=qPCRuniversal, y=(qPCR.Bac.per.cc + qPCR.Arc.per.cc))) + 
    geom_point() +
    geom_smooth(method="lm", se=TRUE, colour="black", linetype=2) +
    geom_abline(slope=1, intercept=0) +
    scale_x_log10(breaks=10^c(3:9)) +
    scale_y_log10(breaks=10^c(3:9)) +
    ylab(expression(atop(paste("Gene Copies ", ml^{-1}), " Bacterial + Archaeal primers", sep=""))) +
    xlab(expression(paste("Gene Copies ", ml^{-1}, ", Universal Primers"))) +
    coord_fixed() +
    theme(axis.text.x=element_text(angle=-45, hjust=0))# +
    #annotation_custom(grob=textGrob("B", x=0.05, y=0.95))
  #print(p_qPCR_sum_bac_arc_v_universal)
  
  ##############
  # Plot of sum (bacterial + archaeal copy number) vs direct-count total cells
  ##############
  
  oneOneLim <- 10^seq(from=3, to=11, by=0.1)
  reasonableCopyNumber <- 3.04 # average of bacterial copy numbers in sequenced genomes
  oneOne <- data.frame(x=oneOneLim, y=oneOneLim, ymax=24*oneOneLim, yReasonable=reasonableCopyNumber*oneOneLim)
  fig3b_data <- all_data[!is.na(all_data$qPCRtotal) &
                           !is.na(all_data$totalcells) &
                           !is.na(all_data$core) &
                           all_data$environment=="sediments" &
                           all_data$qPCRtotal!=0 & all_data$totalcells != 0, ]
  
  #For some reason the papers are out of alphabetical order.
  fig3b_data$core <- factor(fig3b_data$core, levels=sort(unique(fig3b_data$core)), ordered=TRUE)
  
  fig3bm <- lm(log10(qPCRtotal) ~ log10(totalcells), data=fig3b_data)
  rng <- seq(from=min(log10(fig3b_data$totalcells)), to=max(log10(fig3b_data$totalcells)), length.out=100)
  rngDF <- data.frame(totalcells=rng)
  fig3b_pred <- predict(fig3bm,
                        #newdata=rngDF,
                        interval="prediction")
  fig3b_predDF <- as.data.frame(fig3b_pred)
  fig3b_predDF$x <- fig3b_data$totalcells
  fig3b_predDF$span <- fig3b_predDF$upr - fig3b_predDF$lwr
  #ggplot(fig3b_predDF, aes(x=x, y=span)) + geom_line() + scale_y_continuous(limits=c(0, 5))
  mean(fig3b_predDF$span, na.omit=TRUE)
  
  
  fig3b_predDF$lwr <- 10^fig3b_predDF$lwr
  fig3b_predDF$fit <- 10^fig3b_predDF$fit
  fig3b_predDF$upr <- 10^fig3b_predDF$upr
  #fig3b_predDF$lwr[fig3b_predDF$lwr < 1000] <- 1000
  
  #fig3b_predDF
  fig3b_predDF$span <- fig3b_predDF$upr / fig3b_predDF$lwr
  
  # Make the various data frame have max values at 1e11
  oneOne[oneOne$ymax > 1e11, "ymax"] <- 1e11
  oneOne[oneOne$yReasonable > 1e11, "yReasonable"] <- 1e11
  
  fig3b_predDF[fig3b_predDF$lwr < 1e4, "lwr"] <- 1e4
  fig3b_predDF[fig3b_predDF$upr > 1e11, "upr"] <- 1e11
  
  
  p_qPCRtotal_v_totalcells <- ggplot(fig3b_data) + 
    geom_blank(data=fig3b_data, aes(x=totalcells, y=qPCRtotal, colour=core, shape=core)) +
    #geom_blank() +
    #geom_line(data=oneOne, aes(x=x, y=y)) +
    geom_abline(slope=1, intercept=0) +
    geom_ribbon(data=fig3b_predDF, aes(x=x, ymin=lwr, ymax=upr), fill="black", alpha=0.1) +
    geom_line(data=fig3b_predDF, aes(x=x, y=fit), linetype=2) +
    geom_ribbon(data=oneOne, aes(x=x, ymin=y, ymax=ymax), fill="green", alpha=0.2) +
    geom_ribbon(data=oneOne, aes(x=x, ymin=y, ymax=yReasonable), fill="darkgreen", alpha=0.4) +
    geom_point(data=fig3b_data, aes(x=totalcells, y=qPCRtotal, colour=core, shape=core), size=1) +
    geom_smooth(data=fig3b_data, aes(x=totalcells, y=qPCRtotal, colour=core), method="lm", se=FALSE, size=0.2) +
    scale_x_log10(limits = c(1e5, 1e10), breaks=10^(5:10)) +
    scale_y_log10(limits=c(1e4, 1e11), breaks=10^(4:11)) +
    ylab(expression(paste("Gene Copies ", ml^{-1}, ", Bacterial + Archaeal primers"))) +
    scale_shape_manual(values=rep(15:18, nlevels(fig3b_data$core))) +
    coord_fixed() +
    xlab(expression(paste("Total Cells (cells ", ml^{-1}, ")"))) +
    ylab(expression(paste("Gene copies ", ml^{-1}, ", Bacterial + Archaeal primers"))) +
    guides(colour = guide_legend(ncol = 2), shape=guide_legend(ncol=2)) +
    theme(axis.text.x=element_text(angle=-45, hjust=0),
          legend.position="none") #+
   # annotation_custom(grob=textGrob("A", x=0.05, y=0.95))
  #print(p_qPCRtotal_v_totalcells)
 
  
  ####
  #Plot of variance in qPCR and CARDFISH vs variance in total cells
  #####
  
  # Way 1 to pull out desired data. Problem is, the qPCR data set includes CARDFISH data with bad methods
  varVvar_data_qPCR <- corrected_seds[!is.na(corrected_seds$totalcells) &
                                        !is.na(corrected_seds$qPCRtotal), ]
  varVvar_data_CARDFISH <- corrected_seds[!is.na(corrected_seds$totalcells) &
                                            !is.na(corrected_seds$Fish.or.cardFish) &
                                            !is.na(corrected_seds$CARDFISH.Total.per.cc) &
                                            corrected_seds$Fish.or.cardFish == "CARDFISH" &
                                            corrected_seds$Arc.permeabilization == "proteinase K", ]
  varVvar_data <- rbind(varVvar_data_qPCR, varVvar_data_CARDFISH)
  
  # Way 2 to pull out desired data, which is any row with EITHER valid qPCR and totalcells data, 
  #    OR any row with valid CARDFISH data AND the right methods (Arc.permeabilization=="proteinase K")
  varVvar_data <- corrected_seds[ #get rid of bad qPCR data
    (!is.na(corrected_seds$totalcells) & !is.na(corrected_seds$qPCRtotal) |
       # get rid of bad CARDFISH data
       !is.na(corrected_seds$totalcells) &
       !is.na(corrected_seds$CARDFISH.Total.per.cc) &
       !is.na(corrected_seds$Arc.permeabilization) &
       corrected_seds$Arc.permeabilization != "proteinase K"), ]
  
  varVvar_data_CARDFISH <- corrected_seds[!is.na(corrected_seds$totalcells) &
                                            !is.na(corrected_seds$CARDFISH.Total.per.cc) &
                                            !is.na(corrected_seds$Arc.permeabilization) &
                                            corrected_seds$Arc.permeabilization == "proteinase K", ]
  vVv_CF <- ddply(varVvar_data_CARDFISH, .(core), summarise, # there are 4 NAs in here b/c 4 cores only have 1 data point
                  cellsSD=sd(log10(totalcells), na.rm=TRUE), 
                  CARDsd=sd(log10(CARDFISH.Total.per.cc), na.rm=TRUE),
                  nSamples = length(CARDFISH.Total.per.cc))
  
  nrow(vVv_CF[!is.na(vVv_CF$cellsSD) & !is.na(vVv_CF$CARDsd), ])
#   cfm <- lm(CARDsd ~ cellsSD, data=vVv_CF)
#   summary(cfm)
  
  # Variance in total cells vs variance in qPCR
  varVvar_data_qPCR <- corrected_seds[!is.na(corrected_seds$totalcells) &
                                        !is.na(corrected_seds$qPCRtotal), ]
  vVv_q <- ddply(varVvar_data_qPCR, .(core), summarise,
                 cellsSD=sd(log10(totalcells), na.rm=TRUE), 
                 qPCRsd=sd(log10(qPCRtotal), na.rm=TRUE),
                 nSamples = length(qPCRtotal))
  
  nrow(vVv_q[!is.na(vVv_q$cellsSD) & !is.na(vVv_q$qPCRsd), ])
#   qm <- lm(qPCRsd ~ cellsSD, data=vVv_q)
#   summary(qm)
  
  vVv_CF$qPCRsd <- NA
  vVv_q$CARDsd <- NA
  
  vVv_df <- rbind(vVv_CF, vVv_q)
  vVv_m <- melt(vVv_df, measure.vars=c("CARDsd", "qPCRsd"), variable.name="enumeration_method", value.name="sd")
  
  p_var_v_var <- ggplot(vVv_m, aes(x=cellsSD, y=sd, color=enumeration_method, shape=enumeration_method, fill=enumeration_method)) + 
  #p_var_v_var <- ggplot(vVv_m, aes(x=cellsSD, y=sd, color=enumeration_method, shape=enumeration_method)) + 
    geom_point(pointsize=ps) +
    #geom_smooth(method="lm", se=FALSE, linetype=2, colour="black") +
    geom_smooth(method="lm", se=TRUE, linetype=2, alpha=0.25) +
    geom_abline(slope=1, intercept=0) +
    scale_colour_manual(values=c("red", "blue"),
                        labels=c("CARD-FISH\nWith Proteinase K", "qPCR"),
                        name="Enumeration Method") +
    scale_fill_manual(values=c("red", "blue"),
                        labels=c("CARD-FISH\nWith Proteinase K", "qPCR"),
                        name="Enumeration Method") +
    scale_shape_manual(values=c(19, 1),
                       labels=c("CARD-FISH\nWith Proteinase K", "qPCR"),
                       name="Enumeration Method") +
    xlab("Std. dev. of log of total cells") +
    ylab("Std. dev. of log of qPCR or CARD-FISH") +
    coord_cartesian(xlim=seq(0, 1.25, by=0.25), ylim=seq(0, 1.25, by=0.25)) +
    theme(legend.position="top") #+
    #annotation_custom(grob=textGrob("B", x=0.05, y=0.95))
  #print(fig3c)
  ggsave("~/Dropbox/Metadata Analysis/Revision for AEM/plots/AAA_experimentalFig4c.png", height=3.5, width=3.04, units="in", dpi=900, type="cairo")
  
  # Used later for statistical analysis
  save(fig3a_df, fig3b_data, vVv_q, vVv_CF, file="data/qPCR_stats_data.RData")
  
  # qPCR total vs qPCR Universal
  m_tot_v_univ <- lm(log10(qPCR.Bac.per.cc + qPCR.Arc.per.cc) ~ depth_log10, data=fig3a_df)
  summary(m_tot_v_univ)
  
  list(p_qPCR_sum_bac_arc_v_universal=p_qPCR_sum_bac_arc_v_universal, 
       p_qPCRtotal_v_totalcells=p_qPCRtotal_v_totalcells, 
       p_var_v_var=p_var_v_var,
       qPCR_total_v_totalcell_df=fig3b_data,
       qPCR_total_v_qPCR_universal_df=fig3a_df,
       qPCR_variance_df=vVv_m
       )
}