##' Plots depth profiles of various things in seawater
##' 
##' @description Creates depth profiles in figures 5 a-c: Bacteria, Archaea, and %Archaea as a function of depth in seawater. Also performs breakpoint analysis
##' @param corrected_sw Data frame of seawater data
##' @param point_size point size in plot
##' @export

sw_depth_profiles <- function(corrected_sw, point_size=0.75) {
  #browser()
  # Point size for all figures
  ps <- point_size
  
  # Set up seawater data frames:
  # Bacteria in seawater, by CARDFISH
  s3bac_df <- corrected_sw[!is.na(corrected_sw$depth_log10) &
                             !is.na(corrected_sw$CARDFISH.Bac.per.cc), 
                           c("paper", "depth_log10", "CARDFISH.Bac.per.cc")]
  s3bac_df$CARDFISH.Bac.per.cc[s3bac_df$CARDFISH.Bac.per.cc == 0] <- 1 # set zeros to one, to avoid problems with log fits
  s3bac_df$logBac <- log10(s3bac_df$CARDFISH.Bac.per.cc)
  
  # Archaea in seawater, by CARDFISH
  s3arc_df <- corrected_sw[!is.na(corrected_sw$depth_log10) &
                             !is.na(corrected_sw$CARDFISH.Arc.per.cc), 
                           c("paper", "depth_log10", "CARDFISH.Arc.per.cc")]
  s3arc_df$CARDFISH.Bac.per.cc[s3arc_df$CARDFISH.Arc.per.cc == 0] <- 1 # set zeros to one, to avoid problems with log fits
  s3arc_df$logArc <- log10(s3arc_df$CARDFISH.Arc.per.cc)
  
  # % Archaea (not exactly the same data points as above, since some studies report % Archaea but not actual numbers)
  s3percent_df <- corrected_sw[!is.na(corrected_sw$depth_log10) &
                                 !is.na(corrected_sw$Fraction.Arc.CARDFISH), 
                               c("paper", "depth_log10", "Fraction.Arc.CARDFISH")]
  
  # Breakpoint analysis: Bacteria
  bac_depthM <- lm(logBac ~ depth_log10, data=s3bac_df)
  seg_bac_depthM <- segmented(bac_depthM, seg.Z = ~ depth_log10, psi=2)
  summary(seg_bac_depthM)
  bacBreak <- summary(seg_bac_depthM)$psi[2] #1.81
  AIC_lik(AIC(bac_depthM, seg_bac_depthM))
  
  # Fakey way to give linear models (really I should get them from summary(seg_bac_depthM), but I don't totally understand the output)
  shallowBac <- lm(logBac ~ depth_log10, data=subset(s3bac_df, depth_log10 <= bacBreak))
  summary(shallowBac)
  deepBac <- lm(logBac ~ depth_log10, data=subset(s3bac_df, depth_log10 > bacBreak))
  summary(deepBac)
  
  
  # Breakpoint analysis: Archaea
  arc_depthM <- lm(logArc ~ depth_log10, data=s3arc_df)
  seg_arc_depthM <- segmented(arc_depthM, seg.Z = ~ depth_log10, psi=2)
  summary(seg_arc_depthM)
  arcBreak <- summary(seg_arc_depthM)$psi[2] #2.57
  
  # Fakey way to give linear models (really I should get them from summary(seg_bac_depthM), but I don't totally understand the output)
  shallowArc <- lm(logArc ~ depth_log10, data=subset(s3arc_df, depth_log10 <= arcBreak))
  summary(shallowArc)
  deepArc <- lm(logArc ~ depth_log10, data=subset(s3arc_df, depth_log10 > arcBreak))
  summary(deepArc)
  
  # Breakpoint analysis: Percent Archaea
  percent_depthM <- lm(Fraction.Arc.CARDFISH ~ depth_log10, s3percent_df)
  seg_percent_depthM <- segmented(percent_depthM, seg.Z= ~ depth_log10, psi=2)
  summary(seg_percent_depthM)
  percentBreak <- summary(seg_percent_depthM)$psi[2] #1.71
  
  # Fakey way to give linear models (really I should get them from summary(seg_bac_depthM), but I don't totally understand the output)
  shallowPer <- lm(Fraction.Arc.CARDFISH ~ depth_log10, data=subset(s3percent_df, depth_log10 <= percentBreak))
  summary(shallowPer)
  deepPer <- lm(Fraction.Arc.CARDFISH ~ depth_log10, data=subset(s3percent_df, depth_log10 > percentBreak))
  summary(deepPer)
  
  # Seawater figure limits and breaks
  sw_cell_limits <- 10^c(2.8, 6.3)
  sw_cell_breaks <- 10^(3:6)
  sw_depth_limits <- c(4, 0)
  
  sw_depth_breaks <- sw_depth_limits[1]:sw_depth_limits[2]
  fig5sw_bac <- ggplot() + 
    geom_point(data=s3bac_df, aes(x=depth_log10, y=CARDFISH.Bac.per.cc), size=ps) +
    geom_smooth(data=subset(s3bac_df, depth_log10 <=bacBreak), 
                aes(x=depth_log10, y=CARDFISH.Bac.per.cc), method="lm", colour="black", linetype=2) +
    geom_smooth(data=subset(s3bac_df, depth_log10 > bacBreak), 
                aes(x=depth_log10, y=CARDFISH.Bac.per.cc), method="lm",  colour="black", linetype=2) +
    scale_y_log10(expression(paste("Bacteria, cells ", ml^{-1})), limits=sw_cell_limits, breaks=sw_cell_breaks) +
    scale_x_reverse(expression(paste(log[10], " depth, m")), limits=sw_depth_limits) +
    coord_flip()
  #print(fig5sw_bac)
  
  fig5sw_arc <- ggplot(s3arc_df, aes(x=depth_log10, y=CARDFISH.Arc.per.cc)) + geom_point(size=ps) +
    geom_smooth(data=subset(s3arc_df, depth_log10 <=arcBreak), 
                aes(x=depth_log10, y=CARDFISH.Arc.per.cc), method="lm", colour="black", linetype=2) +
    geom_smooth(data=subset(s3arc_df, depth_log10 > arcBreak), 
                aes(x=depth_log10, y=CARDFISH.Arc.per.cc), method="lm",  colour="black", linetype=2) +
    scale_y_log10(expression(paste("Archaea, cells ", ml^{-1})), limits=sw_cell_limits, breaks=sw_cell_breaks) +
    scale_x_reverse(expression(paste(log[10], " depth, m")), limits=sw_depth_limits) +
    coord_flip()
  #print(fig5sw_arc)
  
  fig5sw_percent <- ggplot(s3percent_df, aes(x=depth_log10, y=Fraction.Arc.CARDFISH)) + geom_point(size=ps) +
    geom_smooth(data=subset(s3percent_df, depth_log10 <=percentBreak), 
                aes(x=depth_log10, y=Fraction.Arc.CARDFISH), method="lm", colour="black", linetype=2) +
    geom_smooth(data=subset(s3percent_df, depth_log10 > percentBreak), 
                aes(x=depth_log10, y=Fraction.Arc.CARDFISH), method="lm",  colour="black", linetype=2) +
    scale_y_continuous(expression(Archaea / (Archaea + Bacteria)), 
                       limits=c(-0.02, 1),
                       labels=percent) +
    scale_x_reverse(expression(paste(log[10], " depth, m")), limits=sw_depth_limits) +
    coord_flip()
  #print(fig5sw_percent)
  #browser()
  
  sw_bac_shallow_stats <- lm_stats(subset(s3bac_df, depth_log10 <= bacBreak), xvar="depth_log10", yvar="logBac")
  sw_bac_deep_stats <- lm_stats(subset(s3bac_df, depth_log10 > bacBreak), xvar="depth_log10", yvar="logBac")
  sw_arc_shallow_stats <- lm_stats(subset(s3arc_df, depth_log10 <= arcBreak), xvar="depth_log10", yvar="logArc")
  sw_arc_deep_stats <- lm_stats(subset(s3arc_df, depth_log10 > arcBreak), xvar="depth_log10", yvar="logArc")
  sw_frac_shallow_stats <- lm_stats(subset(s3percent_df, depth_log10 <= percentBreak), xvar="depth_log10", yvar="Fraction.Arc.CARDFISH")
  sw_frac_deep_stats <- lm_stats(subset(s3percent_df, depth_log10 > percentBreak), xvar="depth_log10", yvar="Fraction.Arc.CARDFISH")
  
  #AIC data frames
  AIC_sw_bac <- AIC_lik(AIC(bac_depthM, seg_bac_depthM))
  AIC_sw_arc <- AIC_lik(AIC(arc_depthM, seg_arc_depthM))
  AIC_sw_frac <- AIC_lik(AIC(percent_depthM, seg_percent_depthM))
  AIC_sw <- rbind(AIC_sw_bac, AIC_sw_arc)
  AIC_sw <- rbind(AIC_sw, AIC_sw_frac)
  
  sw_lm_list <- list(sw_bac_shallow_stats=sw_bac_shallow_stats,
                     sw_bac_deep_stats=sw_bac_deep_stats,
                     sw_arc_shallow_stats=sw_arc_shallow_stats,
                     sw_arc_deep_stats=sw_arc_deep_stats,
                     sw_frac_shallow_stats=sw_frac_shallow_stats,
                     sw_frac_deep_stats=sw_frac_deep_stats) 
  sw_lm_stats <- do.call("rbind", sw_lm_list)
  #browser()
  list(p_bac_profile=fig5sw_bac, 
       p_arc_profile=fig5sw_arc, 
       p_percent_profile=fig5sw_percent, 
       sw_lm_stats=sw_lm_stats, 
       bac_break=bacBreak,
       arc_break=arcBreak,
       frac_break=percentBreak,
       AIC_sw=AIC_sw)
}