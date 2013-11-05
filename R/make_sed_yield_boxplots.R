##' Make boxplots of yield for sediments data, split by various factors
##' 
##' @param corrected_seds Data frame of corrected sediments
##' @param n Number of replicates to be passed on to `aov_perm_test()` 
##' @export


make_sed_yield_boxplots <- function(corrected_seds, n=1000) {
  #browser()
  # All sed yield boxplots will use only data in which permeabilization method containes proteinase K
  #   Omit rows in which Arc.permeabilization or FISH yield is NA
  seds_protK_only <- corrected_seds[!is.na(corrected_seds$Arc.permeabilization) &
                                      !is.na(corrected_seds$FISH.yield), ]
  #   Keep only the samples where Fish.or.cardFish is  
  seds_protK_only <- seds_protK_only[seds_protK_only$Fish.or.cardFish=="FISH" |
                                       (seds_protK_only$Fish.or.cardFish=="CARDFISH" &
                                          seds_protK_only$Arc.permeabilization=="proteinase K"), ]
  
  
  ########
  # Yield by cell stain (AO or DAPI or SYBR)
  ########
  sed_yield_stain <- seds_protK_only[!is.na(seds_protK_only$Cell.stain), ]
  
  cell.stain_p <-aov_perm_test(n, sed_yield_stain, data_var="FISH.yield", cat_var="Cell.stain", 
                                          title="sediments: cell stain",
                                          fn="plots/seds_cell_stain.png")
  
  aov_cell_stain <- aov(FISH.yield ~ Cell.stain, data=sed_yield_stain[sed_yield_stain$Cell.stain!="SYBR Green" & sed_yield_stain$Cell.stain!="SYBR Green II", ])
  print(TukeyHSD(aov_cell_stain))
  cell_stain_text <- data.frame(x=levels(as.factor(as.character(sed_yield_stain[sed_yield_stain$Cell.stain!="SYBR Green" & sed_yield_stain$Cell.stain!="SYBR Green II", ]$Cell.stain))), 
                                y=rep(0.1, 4),
                                label=c("a", "b", "c", "c"))
  
  p_Cell.stain <- single_sw_yield_boxplot(d=sed_yield_stain, xvar="Cell.stain", yvar="FISH.yield", sig_text=significance_labeller(cell.stain_p)) +
    xlab("cell stain") +
    ylab("yield") +
    geom_text(data=cell_stain_text, aes(x=x, y=y, label=label))
  
  ########
  # Compare yield by bacterial probe
  ########
  
  sed_yield_probe <- seds_protK_only[!is.na(seds_protK_only$Bac.Probe) &
                                       seds_protK_only$Bac.Probe!="Eub338 and hAqui1045 and Aqui1197", ]
  
  # Note: no significant differences
  aov_probe <- aov(FISH.yield ~ Bac.Probe, data=sed_yield_probe)
  print(TukeyHSD(aov_probe))
  
  bac.probe_p <- aov_perm_test(n, sed_yield_probe, cat_var="Bac.Probe", data_var="FISH.yield", 
                               title="sediments: Bacterial probe", fn="plots/seds_bac.probe.png")
  p_bac.probe <- single_sw_yield_boxplot(d=sed_yield_probe, xvar="Bac.Probe", yvar="FISH.yield", sig_text=significance_labeller(bac.probe_p)) +
    xlab("Bacterial probe") +
    ylab("yield")
  
  print("**Two data points with Bac.Probe==Eub338 and hAqui1045 and Aqui1197 were excluded, since we can't use those two points for reliable statistics")
  
  
  ####
  # Hybridization method
  ####
  #sed_yield_hyb_method <- seds_protK_only[!is.na(seds_protK_only$Fish.or.cardFish), ]
  sed_yield_hyb_method <- corrected_seds[!is.na(corrected_seds$FISH.yield) &
                                           #!is.na(corrected_seds$Fish.or.cardFish) & 
                                           (corrected_seds$Fish.or.cardFish=="FISH" | 
                                              (corrected_seds$Fish.or.cardFish=="CARDFISH" &
                                              corrected_seds$Arc.permeabilization=="proteinase K")), ]
  sed_yield_hyb_method <- sed_yield_hyb_method[!is.na(sed_yield_hyb_method$Fish.or.cardFish), ]
  
  Fish.or.cardFish_p <- aov_perm_test(n, sed_yield_hyb_method, cat_var="Fish.or.cardFish", data_var="FISH.yield",
                                      title="sediments: Hybridization method",
                                      fn="plots/seds_hyb_method.png")
  
  p_hyb_method <- single_sw_yield_boxplot(d=sed_yield_hyb_method, xvar="Fish.or.cardFish", yvar="FISH.yield",
                                          sig_text=significance_labeller(Fish.or.cardFish_p)) +
    xlab("Hybridization method") +
    ylab("yield")
  

 
  #########
  # Create summary tables
  #########
  
  sed_stain_summary <- ddply(sed_yield_stain, .(Cell.stain), summarise, medianYield=median(FISH.yield), interquartile.range=IQR(FISH.yield))
  sed_bac_probe_summary <- ddply(sed_yield_probe, .(Bac.Probe), summarise, medianYield=median(FISH.yield), interquartile.range=IQR(FISH.yield))
  sed_method_summary <- ddply(sed_yield_hyb_method, .(Fish.or.cardFish), summarise, medianYield=median(FISH.yield), interquartile.range=IQR(FISH.yield))
  
  
  list(p_Cell.stain=p_Cell.stain, p_bac.probe=p_bac.probe, p_hyb_method=p_hyb_method,
       sed_stain_summary=sed_stain_summary,
       sed_bac_probe_summary=sed_bac_probe_summary,
       sed_method_summary=sed_method_summary)
}