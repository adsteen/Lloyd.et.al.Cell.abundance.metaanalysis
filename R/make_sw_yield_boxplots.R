##' Create boxplots of yield in seawater, split by various variables
##' 
##' @param corrected_sw The corrected seawater boxplot
##' @param n_perm Number of permutations for permutation test, passed through to `aov_perm_test`
##' @details Note that this throws some warnings about points being omitted. This is simply because the outliers are turned off on the boxplots, which is fine because the outliers are nevertheless represented in the backgound points
##' @export

make_sw_yield_boxplots <- function(corrected_sw, n_perm=1000) {
  
  #########
  # Yield by hybridization method
  ##########
  sw_yield_hyb_method <- corrected_sw[!is.na(corrected_sw$Fish.or.cardFish) &
                              !is.na(corrected_sw$FISH.yield), ]
  sw_yield_hyb_method$Fish.or.cardFish <- factor(sw_yield_hyb_method$Fish.or.cardFish, 
                                       levels=c("FISH", "CARDFISH", "Polyribonucleotide FISH"), ordered=TRUE)
  # AOV of yield by hybridization method
  m_hyb_method <- aov(FISH.yield ~ Fish.or.cardFish, data=sw_yield_hyb_method)
  summary(m_hyb_method)
  #print("Tukey HSD test of seawater yield as a function of hybdridization method")
  #print(TukeyHSD(m_hyb_method))
  tukey_data <- data.frame(x=levels(sw_yield_hyb_method$Fish.or.cardFish), 
                           y=rep(0.1, 3), 
                           label=c("a", "a", "b"))
  perm_p <- aov_perm_test(n_perm, sw_yield_hyb_method, data_var="FISH.yield", cat_var="Fish.or.cardFish", 
                          title="seawater: hybridization method",
                          fn="plots/sw_hybridization_method.png")
  #sig_text <- significance_labeller(summary(aov_FISH_yield)[[1]]$"Pr(>F)"[1])
  sig_text <- significance_labeller(perm_p)
  
  p_FISH_method <- single_sw_yield_boxplot(sw_yield_hyb_method, "Fish.or.cardFish", "FISH.yield", sig_text=sig_text) + 
    xlab("hybridization method") +
    ylab("yield") +
    theme(axis.text.x=element_text(angle=-22.5)) +
    geom_text(data=tukey_data, aes(x=x, y=y, label=label))
    
  
  
  ###########
  # Yield by stain
  ##########
  sw_yield_stain <- corrected_sw[!is.na(corrected_sw$Cell.stain) &
                              !is.na(corrected_sw$FISH.yield), ]
  sw_yield_stain$stain_reduced <- NA
  sw_yield_stain$stain_reduced[sw_yield_stain$Cell.stain=="DAPI"] <- "DAPI"
  sw_yield_stain$stain_reduced[sw_yield_stain$Cell.stain=="AO" | sw_yield_stain$Cell.stain=="AO/isopropanol"] <- "AO"
  sw_yield_stain$stain_reduced <- factor(sw_yield_stain$stain_reduced)
  
  m_stain <- aov(FISH.yield ~ stain_reduced, data=sw_yield_stain)
  #print("AOV of seawater yield as a function of stain")
  #print(summary(m_stain))
  
  
  perm_p <- aov_perm_test(n_perm, sw_yield_stain, data_var="FISH.yield", cat_var="stain_reduced", 
                          title="seawater: stain", fn="plots/sw_stain.png")
  sig_text <- significance_labeller(perm_p)
  
  # Note: The warning here is just becasue boxplot outliers are turned off
  p_stain <- single_sw_yield_boxplot(sw_yield_stain, x="stain_reduced", y="FISH.yield", sig_text=sig_text) +
    xlab("cell stain") +
    ylab("yield") +
    theme(axis.text.x=element_text(angle=-22.5))
  
  
  #############
  # Automated vs manual counting
  ##############
  
  sw_yield_automated <- corrected_sw[!is.na(corrected_sw$Counting.method) & 
                              !is.na(corrected_sw$FISH.yield), ]
  
  m_automated <- aov(FISH.yield ~ Counting.method, data=sw_yield_automated)
  
  perm_p_automated <- aov_perm_test(n_perm, sw_yield_automated, data_var="FISH.yield", cat_var="Counting.method", 
                                    title="seawater: Automated/manual count", fn="plots/sw_automated.png")
  sig_text <- significance_labeller(perm_p_automated)
  
  p_automated <- single_sw_yield_boxplot(sw_yield_automated, "Counting.method", "FISH.yield", 
                                         sig_text=sig_text) +
    xlab("counting method") +
    ylab("yield") +
    theme(axis.text.x=element_text(angle=-22.5))
  #print(p_automated)
  
  
  ###############
  # Archaeal probe
  ###############
  sw_yield_arc_probe <- corrected_sw[!is.na(corrected_sw$Arc.probe) & 
                              !is.na(corrected_sw$FISH.yield) &
                              corrected_sw$Arc.probe != "not tried" &
                              corrected_sw$Arc.permeabilization != "lysozyme", ]
  
  # Shorten the name for the polyribonucleotide FISH method
  sw_yield_arc_probe$corr_probe <- as.character(sw_yield_arc_probe$Arc.probe)
  sw_yield_arc_probe$corr_probe[sw_yield_arc_probe$corr_probe=="Collection of amplicons from SSU and LSU from fosmids of Eurys and Crens"] <-
    "mixed Eury and Cren amplicons"
  
#   # How many samples for each of the Arc.permeabilization methods?
#   count(sw_yield_arc_probe, vars="Arc.permeabilization")
  
  # Order the factor: calculate the median yield for each method
  arc_probe_med <- ddply(sw_yield_arc_probe, .(corr_probe), summarise, medYield = median(FISH.yield))
  sw_yield_arc_probe$corr_probe <- factor(sw_yield_arc_probe$corr_probe, 
                                 levels=arc_probe_med$corr_probe[order(arc_probe_med$medYield, decreasing=TRUE)], 
                                 ordered=TRUE)
  
  # Calculate aov of FISH yield by stain
  aov_probe <- aov(FISH.yield~corr_probe, data=sw_yield_arc_probe)
  #print("Tukey HSD test of FISH yield as a function of probe")
  #print(TukeyHSD(aov_probe)) # Each difference is signficant
  tukey_data <- data.frame(x=levels(sw_yield_arc_probe$corr_probe), 
                           y=rep(0.1, nlevels(sw_yield_arc_probe$corr_probe)), 
                           label=letters[1:nlevels(sw_yield_arc_probe$corr_probe)])
  
  # Calculate p value by permutation test
  perm_p_arc_probe <- aov_perm_test(n_perm, sw_yield_arc_probe, data_var="FISH.yield", cat_var="corr_probe", 
                                title="seawater: Archaeal probe", fn="plots/sw_Archaeal_probe.png")
  sig_text <- significance_labeller(perm_p_arc_probe)
  
  p_arc_probe <- single_sw_yield_boxplot(sw_yield_arc_probe, "corr_probe", "FISH.yield", sig_text=sig_text) +
    xlab("probe") + 
    ylab("yield") +
    theme(axis.text.x=element_text(angle=-22.5)) +
    geom_text(data=tukey_data, aes(x=x, y=y, label=label))
  
  
  ##########
  # Bacterial probe
  #########
  sw_yield_bac_probe <- corrected_sw[!is.na(corrected_sw$Bac.probe) &
                                       !is.na(corrected_sw$FISH.yield), ]
  
  # Replace the long label on one of the probes
  sw_yield_bac_probe$Bac.probe <- as.character(sw_yield_bac_probe$Bac.probe)
  sw_yield_bac_probe$Bac.probe[sw_yield_bac_probe$Bac.probe=="Collection of amplicons from SSU and LSU from natural planktonic samples"] <- 
    "mixed natural amplicons"
  sw_yield_bac_probe$Bac.probe <- as.factor(sw_yield_bac_probe$Bac.probe)
  
  # ORder the factors in decreasing order
  bac_probe_med <- ddply(sw_yield_bac_probe, .(Bac.probe), summarise, medYield = median(FISH.yield))
  sw_yield_bac_probe$Bac.probe <- factor(sw_yield_bac_probe$Bac.probe, 
                                         levels=bac_probe_med$Bac.probe[order(bac_probe_med$medYield, decreasing=TRUE)], 
                                         ordered=TRUE)
  
  m_bac_probe <- aov(FISH.yield ~ Bac.probe, data=sw_yield_bac_probe)
  #print("Tukey HSD testing of FISH yield as a function of bacterial probe")
  #print(TukeyHSD(m_bac_probe))
  tukey_data <- data.frame(x=levels(sw_yield_bac_probe$Bac.probe), 
                           y=rep(0.1, nlevels(sw_yield_bac_probe$Bac.probe)), 
                           label=c("a", "b", "b"))
  
  
  
  perm_p_bac_probe <- aov_perm_test(n_perm, sw_yield_bac_probe, data_var="FISH.yield", cat_var="Bac.probe", 
                                    title="seawater: Bacterial probe",
                                    fn="plots/sw_Bacterial_probe.png")
  sig_text <- significance_labeller(perm_p_arc_probe)
  
  p_bac_probe <- single_sw_yield_boxplot(sw_yield_bac_probe, "Bac.probe", "FISH.yield", sig_text=sig_text) +
    xlab("Bacterial probe") +
    ylab("yield") +
    theme(axis.text.x=element_text(angle=-22.5)) +
    geom_text(data=tukey_data, aes(x=x, y=y, label=label))
  
  #browser()
  
  ##########
  # Set up data frames of median values
  ##########
  
  hyb_method_summary <- ddply(sw_yield_hyb_method, .(Fish.or.cardFish), summarise, medianYield=median(FISH.yield), interquartile.range=IQR(FISH.yield))
  stain_summary <- ddply(sw_yield_stain, .(stain_reduced), summarise, medianYield=median(FISH.yield), interquartile.range=IQR(FISH.yield))
  automated_summary <- ddply(sw_yield_automated, .(Counting.method), summarise, medianYield=median(FISH.yield), interquartile.range=IQR(FISH.yield))
  arc_probe_summary <- ddply(sw_yield_arc_probe, .(corr_probe), summarise, medianYield=median(FISH.yield), interquartile.range=IQR(FISH.yield))                       
  bac_probe_summary <- ddply(sw_yield_bac_probe, .(Bac.probe), summarise, medianYield=median(FISH.yield), interquartile.range=IQR(FISH.yield))
  
  save(sw_yield_hyb_method, sw_yield_stain, sw_yield_automated, sw_yield_arc_probe, sw_yield_bac_probe,
       m_hyb_method, m_stain, m_automated, aov_probe, m_bac_probe,
       hyb_method_summary, stain_summary, automated_summary, arc_probe_summary, bac_probe_summary,
       file="data/sw_yield_dfs.RData")
  
  list(p_FISH_method=p_FISH_method, 
       p_stain=p_stain, 
       p_automated=p_automated, 
       p_arc_probe=p_arc_probe, 
       p_bac_probe=p_bac_probe,
       hyb_method_summary=hyb_method_summary,
       stain_summary=stain_summary,
       automated_summary=automated_summary,
       arc_probe_summary=arc_probe_summary,
       bac_probe_summary=bac_probe_summary,
       m_hyb_method=m_hyb_method,
       m_stain=m_stain, 
       m_automated=m_automated, 
       m_probe=aov_probe, 
       m_bac_probe=m_bac_probe)
  
}