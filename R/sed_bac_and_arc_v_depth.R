##' Makes plots and models of sediments bac and arc concentrations vs depth
##' 
##' @param corrected_sed Data frame of sediments data
##' @param point_size Point size for the plot
##' @param reproduce_exactly If TRUE, set the state of the random number generator (set.seed(2112))
##' @export

sed_bac_and_arc_v_depth <- function(corrected_sed, point_size=0.75, reproduce_exactly=FALSE) {
  
  # If you want to reproduce results from the manuscript precisely, set the seed
  if (reproduce_exactly) {
    set.seed(2112)
  }
  
  # For 'estimate' absolute quantity arc and bac, calculate median best-practice FISH yield
  protK_yield <- median(corrected_seds$FISH.yield[corrected_seds$Arc.permeabilization=="proteinase K"], na.rm=TRUE)
  
  ######
  # Measured (by best-practices) and estimated absolute quantites of Archaea in sediments
  ######
  
  sed_arc_CARD_counts <- corrected_seds[!is.na(corrected_seds$Environment.Type) &
                                          !is.na(corrected_seds$depth_log10) &
                                          #!is.na(corrected_seds$totalcells) &
                                          !is.na(corrected_seds$CARDFISH.Arc.per.cc) &
                                          !is.na(corrected_seds$Arc.permeabilization) &
                                          corrected_seds$Environment.Type != "Intertidal" &
                                          corrected_seds$Environment.Type != "Salt marsh" &
                                          corrected_seds$Arc.permeabilization=="proteinase K",
                                        c("paper", "depth_log10", "CARDFISH.Arc.per.cc")]
  sed_arc_CARD_counts$method <- "measured"
  
  # Rename the 'value' column so it will correspond to the 'estimtaed' data frame
  sed_arc_CARD_counts <- rename(sed_arc_CARD_counts, c("CARDFISH.Arc.per.cc" = "Arc.per.cc"))
  
  # Create a data frame to use to calculate _estimated_ Archaea
  sed_arc_est <- corrected_seds[!is.na(corrected_seds$Environment.Type) &
                                        !is.na(corrected_seds$depth_log10) &
                                       # !is.na(corrected_seds$totalcells) &
                                        !is.na(corrected_seds$percentqPCR) &
                                        corrected_seds$Uses.516.for.Arc==FALSE &
                                        corrected_seds$Environment.Type != "Intertidal" &
                                        corrected_seds$Environment.Type != "Salt marsh" ,
                                      c("paper", "depth_log10", "totalcells", "percentqPCR")]
  sed_arc_est$method <- "estimated"
  sed_arc_est$Arc.per.cc <- sed_arc_est$totalcells * sed_arc_est$percentqPCR * protK_yield
  sed_arc_est <- sed_arc_est[ , -which(names(sed_arc_est) %in% c("totalcells", "percentqPCR"))]
  
  # Set up blue for qPCR (estimate), red for CARD_FISH (measured)
  CARD_color <- brewer.pal(n=3, name="Paired")[2]
  qPCR_color <- brewer.pal(n=3, name="OrRd")[3]
  
  
  # Bind the estimated and CARDFISH counted data frames
  sed_arc <- rbind(sed_arc_CARD_counts, sed_arc_est)
  
  # Sample to avoid systematic overplotting
  sed_arc <- sed_arc[sample(1:nrow(sed_arc)), ]
  
  
  #print(p_sed_arc)
  #ggsave("~/Dropbox/Metadata Analysis/Revision for AEM/plots/draft_fig_arc_by_depth_seds.png", p_sed_arc, height=3, width=clm2/3, units="in", dpi=myDPI, type="cairo")
  
  ### Calculate the linear model for each method
  # Measured by CARD
  m_arc_meas <- lm(log10(Arc.per.cc) ~ depth_log10, data=sed_arc_CARD_counts)
  #print("Linear model of absolute Archaea (CARDFISH counts) vs. depth")
  #print(summary(m_arc_meas))
  
  # Estimated from qPCR and total cells
  m_arc_est <- lm(log10(Arc.per.cc) ~ depth_log10, data=sed_arc_est)
  #print("Linear model of estimated Archaea from qPCR-% Archaea and total cells")
  #print(summary(m_arc_est))
 
  
  sed_arc_est$log10_arc <- log10(sed_arc_est$Arc.per.cc)
  m_sed_arc_est <- lm(log10_arc ~ depth_log10, data=sed_arc_est)
  m_sed_arc_est_seg <- segmented(m_sed_arc_est, seg.Z=~depth_log10, psi=-1)
  break_sed_arc_est <- summary(m_sed_arc_est_seg)$psi[2]
  
  arc_est_AIC <- AIC_lik(AIC(m_sed_arc_est, m_sed_arc_est_seg))
  
  ### Segmented models
  #m_arc_meas_seg <- segmented(m_arc_meas, seg.Z=~depth_log10, psi=-1) This does not work for some reason
  sed_arc_CARD_counts$log10_arc <- log10(sed_arc_CARD_counts$Arc.per.cc)
  m_arc_meas2 <- lm(log10_arc ~ depth_log10, data=subset(sed_arc_CARD_counts, method="percentqPCR"))
  m_arc_meas_seg <- segmented(m_arc_meas2, seg.Z=~depth_log10, psi=-1)
  sed_arc_break_meas <- summary(m_arc_meas_seg)$psi[2] #12 cm
  arc_meas_AIC <- AIC_lik(AIC(m_arc_meas2, m_arc_meas_seg))
  
  ######
  # Measured (by best-practices) and estimated absolute quantites of Bacteria in sediments
  ######
  
  # Counts of bacteria based on CARD_FISH best practices
  sed_bac_CARD_counts <- corrected_seds[!is.na(corrected_seds$Environment.Type) &
                                                 !is.na(corrected_seds$depth_log10) &
                                                 !is.na(corrected_seds$totalcells) &
                                                 !is.na(corrected_seds$CARDFISH.Bac.per.cc) &
                                                 corrected_seds$CARDFISH.Bac.per.cc != 0 &
                                                 corrected_seds$Environment.Type != "Intertidal" &
                                                 corrected_seds$Environment.Type != "Salt marsh",
                                               c("paper", "depth_log10", "CARDFISH.Bac.per.cc")]
  sed_bac_CARD_counts$method <- "measured"
  sed_bac_CARD_counts <- rename(sed_bac_CARD_counts, c("CARDFISH.Bac.per.cc" = "Bac.per.cc"))
  
  
  # Estimates of bacteria based on total cell counts and % Bacteria from qPCR (not using 516 for Archaea)
  sed_bac_est <- corrected_seds[!is.na(corrected_seds$Environment.Type) &
                                        !is.na(corrected_seds$depth_log10) &
                                        !is.na(corrected_seds$totalcells) &
                                        !is.na(corrected_seds$percentqPCR) &
                                        corrected_seds$Uses.516.for.Arc==FALSE &
                                        corrected_seds$Environment.Type != "Intertidal" &
                                        corrected_seds$Environment.Type != "Salt marsh" ,
                                      c("paper", "depth_log10", "totalcells", "percentqPCR")]
  sed_bac_est$method <- "estimated"
  sed_bac_est$Bac.per.cc <- sed_bac_est$totalcells * (1-sed_bac_est$percentqPCR) * protK_yield
  sed_bac_est <- sed_bac_est[ , -which(names(sed_bac_est) %in% c("totalcells", "percentqPCR"))]
  
  # Bind the resulting data frames
  sed_bac <- rbind(sed_bac_CARD_counts, sed_bac_est)
  
  # Sample to avoid systematic overplotting
  sed_bac <- sed_bac[sample(1:nrow(sed_bac)), ]
  
  # Some linear models
  m_bac_meas <- lm(log10(Bac.per.cc) ~ depth_log10, data=sed_bac_CARD_counts)
  #print("Linear model of MEASURED log10 of Bac per cc as a function of log10 of depth")
  summary(m_bac_meas)
  
  m_bac_est <- lm(log10(Bac.per.cc) ~ depth_log10, data=sed_bac_est)
  #print("Linear model of ESTIMATED log10 of Bac per cc as a function of log10 of depth")
  #print(summary(m_bac_est))
  
  sed_bac_CARD_counts$log10_bac <- log10(sed_bac_CARD_counts$Bac.per.cc)
  m_bac_meas2 <- lm(log10_bac ~ depth_log10, data=sed_bac_CARD_counts)
  m_bac_meas_seg <- segmented(m_bac_meas2, seg.Z = ~depth_log10, psi=0)
  bac_meas_break <- summary(m_bac_meas_seg)$psi[2]
  bac_meas_AIC <- AIC_lik(AIC(m_bac_meas2, m_bac_meas_seg))
  
  sed_bac_est$log10_bac <- log10(sed_bac_est$Bac.per.cc)
  m_bac_est2 <- lm(log10_bac ~ depth_log10, data=sed_bac_est)
  m_bac_est_seg <- segmented(m_bac_est2, seg.Z=~depth_log10, psi=1)
  bac_est_break <- summary(m_bac_est_seg)$psi[2]
  bac_est_AIC <- AIC_lik(AIC(m_bac_est2, m_bac_est_seg))
  
  bac_depth_meas_shallow <- lm(log10_bac ~ depth_log10, data=subset(sed_bac_CARD_counts, depth_log10 <= bac_meas_break))
  bac_depth_meas_deep <- lm(log10_bac ~ depth_log10, data=subset(sed_bac_CARD_counts, depth_log10 >= bac_meas_break))
  bac_depth_est_shallow <- lm(log10_bac ~ depth_log10, data=subset(sed_bac_est, depth_log10 <= bac_est_break))
  
  sed_bac$log10_bac <- log10(sed_bac$Bac.per.cc)
  
  # Plot of Archaea vs depth (no breakpoints)
  p_sed_arc <- ggplot(sed_arc, aes(x=depth_log10, y=Arc.per.cc, colour=method)) + 
    geom_point(size=point_size) +
    geom_smooth(data=sed_arc[sed_arc$depth_log10<=bac_meas_break, ],
                             method="lm", linetype=2) +
    scale_colour_manual(values=c(CARD_color, qPCR_color),
                        breaks=c("measured", "estimated"),
                        labels=c("measured by\nCARD-FISH", "estimated\nfrom qPCR")) +
    scale_x_reverse(limits=c(2.7, -2)) + 
    scale_y_log10(limits=c(1e4 , 1e11), breaks=10^(5:10)) +
    xlab(expression(paste(log[10], " depth, m"))) +
    ylab(expression(paste("Archaea, cells ", ml^{-1}))) +
    coord_flip() +
    theme(legend.position="top",
          axis.text.x=element_text(angle=-45, hjust=0))
  
  # Segmented plot of bacterial numbers with depth
  p_sed_bac <- ggplot(sed_bac, aes(x=depth_log10, y=Bac.per.cc, colour=method)) + 
    geom_point(size=point_size) + 
    geom_smooth(data=sed_bac_CARD_counts[sed_bac_CARD_counts$depth_log10 <= bac_meas_break, ],
                aes(x=depth_log10, y=Bac.per.cc), method="lm", linetype=2) +
    geom_smooth(data=sed_bac_CARD_counts[sed_bac_CARD_counts$depth_log10 >= bac_meas_break, ],
                aes(x=depth_log10, y=Bac.per.cc), method="lm", linetype=2) +
    geom_smooth(data=sed_bac_est[sed_bac_est$depth_log10 <= bac_est_break, ], 
                method="lm", linetype=2) +
    scale_colour_manual(values=c(CARD_color, qPCR_color),
                        breaks=c("measured", "estimated"),
                        labels=c("measured by\nCARD-FISH", "estimated\nfrom qPCR")) +
    xlab(expression(paste(log[10], " depth, m"))) +
    ylab(expression(paste("Bacteria, cells ", ml^{-1}))) +
    scale_x_reverse(limits=c(2.7, -2)) + 
    scale_y_log10(limits=c(1e4 , 1e11), breaks=10^(5:10)) +
#     scale_linetype_manual(values=c(1, 2)),
#                            breaks=c("Fraction.Arc.CARDFISH", "percentqPCR"),
#                           labels=c("CARD-FISH\ncounts", "estimated\nfrom qPCR")) +
#     scale_shape_manual(values=c(19, 1), ,
#                        breaks=c("Fraction.Arc.CARDFISH", "percentqPCR"),
#                        labels=c("CARD-FISH\ncounts", "estimated\nfrom qPCR")) +
    coord_flip() +
    theme(legend.position="top",
          axis.text.x=element_text(angle=-45, hjust=0))
  
  #####
  ## Create data frame of slopes and standard errors
  ####
  sed_bac_CARD_counts$log10.cells.per.cc <- log10(sed_bac_CARD_counts$Bac.per.cc)
  sed_bac_est$log10.cells.per.cc <- log10(sed_bac_est$Bac.per.cc)
  sed_arc_est <- rename(sed_arc_est, c("log10_arc" = "log10.cells.per.cc"))
  sed_arc_meas <- rename(sed_arc_CARD_counts, c("log10_arc" = "log10.cells.per.cc"))
  #sed_arc_est_stats <- lm_stats(sed_arc_est, xvar="depth_log10", yvar="log10_arc")
  lms_to_make <- list(bac_meas_shallow = subset(sed_bac_CARD_counts, depth_log10 <= bac_meas_break),
                      bac_meas_deep = subset(sed_bac_CARD_counts, depth_log10 > bac_meas_break),
                      bac_est_shallow = subset(sed_bac_est, depth_log10 <= bac_est_break),
                      #arc_est_all = sed_arc_est,
                      arc_est_all = subset(sed_arc_est, depth_log10 <= log10(12)),  # Density decreases badly below 12 m; don't want to include those points in the lm
                      arc_meas_all = sed_arc_meas) # No reason to subset the CARDFISH deptsh, because they only go down to ~= 0.75 m
  sed_lm_summary <- ldply(lms_to_make, lm_stats, xvar="depth_log10", yvar="log10.cells.per.cc")
  
  list(p_sed_arc=p_sed_arc,
       p_sed_bac=p_sed_bac,
       #m_arc_meas=m_arc_meas,
       #m_arc_est=m_arc_est,
       #m_bac_meas=m_bac_meas,
       #m_bac_est=m_bac_est,
       bac_depth_meas_shallow=bac_depth_meas_shallow,
       bac_depth_meas_deep=bac_depth_meas_deep,
       bac_depth_est_shallow=bac_depth_est_shallow,
       sed_lm_summary=sed_lm_summary,
       arc_meas_AIC=arc_meas_AIC, 
       arc_meas_break=sed_arc_break_meas,
       arc_est_AIC=arc_est_AIC,
       arc_est_break=break_sed_arc_est,
       bac_meas_AIC=bac_meas_AIC,
       bac_meas_break=bac_meas_break,
       bac_est_AIC=bac_est_AIC,
       bac_est_break=bac_est_break)  
  
}