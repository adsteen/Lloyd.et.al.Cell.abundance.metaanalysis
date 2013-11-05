##' Make depth profiles and calculate linear models for sediment percent Archaea vs depth
##' 
##' @description `sed_percent_arc_v_depth`
##' @param corrected_seds Data frame of sediments data
##' @param point_size point size for the saved plots
##' @return Returns a list with the following elements
##' 
##' * p_sed_frac_arc_v_depth, a plot of the fraction of archaea vs depth in sediments
##' * m_sed_frac_arc_depth_qPCR, linear model of fraction of archaea vs depth by qPCR
##' * sed_frac_CARD_all_stats, correlation coefficients & other statistics from linear models of CARD-FISH counts vs depth in sediments
##' * sed_frac_est_all_stats, correlation coefficients & other statistics from linear models of fraction of Archaea vs depth in sediments
##' * frac_est_AIC, data frame of Akaike Inforamtion Criterion results from segmented and unsegmented linear models
##' * frac_est_break, breakpoint (depth, m) in fraction of Archaea by qPCR
##' @export
##' 
sed_percent_arc_v_depth <- function(corrected_seds, point_size=0.75) {
  ##########
  # Percent Archaea by depth, based on best-practices qPCR and CARDFISH data
  ##########
  
  # Cut out all the points that have either percentqPCR or Fraction.Arc.CARDFISH, keep only the relevant columns
  sed_frac_arc <- corrected_seds[
    # Omit intertidal and salt marsh sediments
    corrected_seds$Environment.Type!="Intertidal" & corrected_seds$Environment.Type!="Salt marsh" & 
      # Omit if depth is missing
      !is.na(corrected_seds$depth_log10) ,
    c("depth_log10", "percentqPCR", "Fraction.Arc.CARDFISH", "Arc.permeabilization", "Fish.or.cardFish", "Uses.516.for.Arc")]
  
  # Melt for plotting
  sed_frac_arc_m <- melt(sed_frac_arc, measure.vars=c("percentqPCR", "Fraction.Arc.CARDFISH"), variable.name="method", value.name="fraction.Arc")
  
  
  # KEEP ONLY BEST PRACTICE POINTS:
  #    qPCR points that DON'T use 516 for Arc
  #    CARDFISH points that use proteinase K for permeabilization
  sed_frac_arc_m_edit <- sed_frac_arc_m[!is.na(sed_frac_arc_m$fraction.Arc) &
                                          (sed_frac_arc_m$method=="percentqPCR" & 
                                             sed_frac_arc_m$Uses.516.for.Arc == FALSE) |
                                          (sed_frac_arc_m$method=="Fraction.Arc.CARDFISH" & 
                                             sed_frac_arc_m$Arc.permeabilization=="proteinase K"), ]
  
  # Set up blue for qPCR, red for CARD_FISH
  qPCR_color <- brewer.pal(n=3, name="Paired")[2]
  CARD_color <- brewer.pal(n=3, name="OrRd")[3]
  
  # Change the order of the methods, to match the direct count panels
  sed_frac_arc_m_edit$method <- factor(sed_frac_arc_m_edit$method, 
                                       levels=levels(sed_frac_arc_m_edit$method)[2:1],
                                       ordered=TRUE)
  
 
  #######
  # Make linear models for % Archaea vs depth (segmented and unsegmented)
  #######
  
  ### Unsegmented
  # By qPCR
  m_sed_frac_arc_depth_qPCR <- lm(fraction.Arc ~ depth_log10, 
                                  data=subset(sed_frac_arc_m_edit, method=="percentqPCR"))
  #print("Model of % Archaea by qPCR in sediments")
  #print(summary(m_sed_frac_arc_depth_qPCR))
  
  # By CARDFISH
  m_sed_frac_arc_depth_CARD <- lm(fraction.Arc ~ depth_log10,
                                  data=subset(sed_frac_arc_m_edit, method=="Fraction.Arc.CARDFISH"))
  #print("Model of % Archaea by CARDFISH in sediments")
  #print(summary(m_sed_frac_arc_depth_CARD))
  
  #####
  ## Make segmented models (although I don't think I will use them)
  #####
  # Segmented model of % Archaea by CARDFISH
  #m_sed_frac_arc_CARD_seg <- segmented(m_sed_frac_arc_depth_CARD, seg.Z = ~ depth_log10, psi=-1) #This throws an error b.c breakpoint is at edge of data
  #sed_frac_arc_breakpoint <- summary(m_sed_frac_arc_CARD_seg)$psi[2]
  #davies.test(m_sed_frac_arc_depth_CARD, seg.Z=~depth_log10)
  
  #frac_meas_AIC <- AIC_lik(AIC(m_sed_frac_arc_CARD_seg, m_sed_frac_arc_depth_CARD))
  #BIC(m_sed_frac_arc_CARD_seg, m_sed_frac_arc_depth_CARD)
  ## The more _negative_ AIC is the best model: so this function seems to be working 
  
  # Segmented model of % Archaea by qPCR
  m_sed_frac_arc_qPCR_seg <- segmented(m_sed_frac_arc_depth_qPCR, seg.Z = ~ depth_log10, psi=0)
  frac_qPCR_break <- summary(m_sed_frac_arc_qPCR_seg)$psi[2] # 0.244
  frac_est_AIC <- AIC_lik(AIC(m_sed_frac_arc_qPCR_seg, m_sed_frac_arc_depth_qPCR))
  
  
  #########
  # Build the plot
  ########
  sed_frac_arc_m_edit <- sed_frac_arc_m_edit[!is.na(sed_frac_arc_m_edit$method), ]
  
  p_sed_frac_arc_v_depth <- ggplot(sed_frac_arc_m_edit, aes(x=depth_log10, y=fraction.Arc, colour=method)) + 
    geom_point(size=point_size) + 
    geom_smooth(data=sed_frac_arc_m_edit[sed_frac_arc_m_edit$depth_log10<=1, ], aes(x=depth_log10, y=fraction.Arc, colour=method) ,
                method="lm", linetype=2) +
    scale_x_reverse(expression(paste(log[10], " of depth, m"))) +
    scale_y_continuous(expression(Archaea / (Archaea + Bacteria)),
                       labels=percent) +
    scale_colour_manual(values=c(CARD_color, qPCR_color),
                        breaks=c("Fraction.Arc.CARDFISH", "percentqPCR"),
                        labels=c("measured by\nCARD-FISH", "estimated\nfrom qPCR")) +
    coord_flip() +
    theme(legend.position="top")
  
  ##########
  # Return the plot and the linear model and lm_stats
  #########
  # Note: The density of data falls off a cliff at about 12 m (log10(12)~=1.08) so I've cut these out
  sed_frac_CARD_shallow_stats <- lm_stats(subset(sed_frac_arc_m_edit, depth_log10<=1.08 & method=="Fraction.Arc.CARDFISH"), 
                                          xvar="depth_log10", yvar="fraction.Arc")
  sed_frac_est_shallow_stats <- lm_stats(subset(sed_frac_arc_m_edit, depth_log10<=1.08 & method=="percentqPCR"),
                                         xvar="depth_log10", yvar="fraction.Arc")
  
  list(p_sed_frac_arc_v_depth=p_sed_frac_arc_v_depth,
       m_sed_frac_arc_depth_qPCR=m_sed_frac_arc_depth_qPCR,
       sed_frac_CARD_all_stats=sed_frac_CARD_shallow_stats,
       sed_frac_est_all_stats=sed_frac_est_shallow_stats,
       frac_est_AIC=frac_est_AIC,
       frac_est_break=frac_qPCR_break)
}